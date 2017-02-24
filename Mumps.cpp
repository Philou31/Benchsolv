//! \file functions.cpp
//! \brief Functions for the project
//! \author filou
//! \version 0.1
//! \date 13/11/2016, 18:42
//!

#include "Mumps.h"

Mumps::Mumps(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, int par, int sym, int distr, 
        std::string loc, int loc_option, int format, int comm, MPI_Comm mpi_comm, 
        std::string pb_spec_file, int int_opt_key, int int_opt_value, int nrows,
        int ncols, int nz, int mem_relax, float mem_factor, int erranal):
    Solver(test_id, file_A, n_present_A, file_b, n_present_b, nrows, ncols, nz), 
    _mpi_comm(mpi_comm), _pb_spec_file(pb_spec_file), _distr(distr), 
    _loc(loc), _loc_option(loc_option), _format(format), 
    _opt_key(int_opt_key), _opt_value(int_opt_value)
{
    std::cout << "Mumps initialization:\nWorking host: " << par << 
        "\nSymmetry: " << sym << "\nCommunicator: " << comm << "\nFormat: " << 
        format << "\nDistribution: " << distr << "\nDistributed inputs: " << 
        loc << "\n";
    _id.par = par;
    _id.sym = sym;
    _id.comm_fortran = comm;
    _id.nrhs = 1;
    if (mem_relax != 0)
        _mem_relax = mem_relax;
    if (_mem_factor != 1)
        _mem_factor = mem_factor;
    set_opt(parm::ERR_ANALYSIS, erranal);
    init();
    get_A();
    get_b();
}

Mumps::~Mumps() {
    deallocate_A();
    finalize();
}

////////////////////////////////////////////////////
// MPI communication methods
////////////////////////////////////////////////////  

void Mumps::assemble_A() {
    // Gather all local nz in an array
    int *recvcounts;
    if (is_host())
        recvcounts = new int[_nb_procs];
    MPI_Gather(&_id.nz_loc, 1, MPI_INT, recvcounts, 1, MPI_INT, cst::HOST_ID,
            MPI_COMM_WORLD);
    // The displacements of the gather are the same as the receive counts
    int *displs;
    if (is_host()) {
        displs = new int[_nb_procs];
        int i;
        for(i=0;i<_nb_procs;++i) {
            if (i>0)
                displs[i]+=recvcounts[i-1];
        }
        // Gather the matrix A
        _id.irn=new int[_id.nz];
        _id.jcn=new int[_id.nz];
        _id.a=new double[_id.nz];
    }
    MPI_Gatherv(_id.irn_loc, _id.nz_loc, MPI_INT, _id.irn, recvcounts, displs,
        MPI_INT, cst::HOST_ID, MPI_COMM_WORLD);
    MPI_Gatherv(_id.jcn_loc, _id.nz_loc, MPI_INT, _id.jcn, recvcounts, displs,
        MPI_INT, cst::HOST_ID, MPI_COMM_WORLD);
    MPI_Gatherv(_id.a_loc, _id.nz_loc, MPI_DOUBLE, _id.a, recvcounts, displs,
        MPI_DOUBLE, cst::HOST_ID, MPI_COMM_WORLD);
}

bool Mumps::is_host() {
    return _proc_id == cst::HOST_ID;
}

long long Mumps::total_time(long long *t) {
    long long t_tot = 0;
    MPI_Reduce(t, &t_tot, 1, MPI_LONG_LONG, MPI_SUM, cst::HOST_ID, _mpi_comm);
    return t_tot;
}
    
////////////////////////////////////////////////////
// MATRIX INPUT
////////////////////////////////////////////////////
void Mumps::parse_loc_file() {
    _id.nz_loc=0;
    std::clog << "\nGetting loc file: " << _loc << "\n";
    std::ifstream infile(_loc.c_str());
    std::string strInput;
    int id;
    _loc_beg=std::numeric_limits<int>::max();
    _loc_end=0;
    int beg, end, nz;
    // Read the distribution file line per line
    std::getline(infile, strInput);
    while (infile) {
        std::stringstream stream(strInput);
        // chunk_id
        stream >> id;
        // if the chunk has to be taken by the current process (ex. 128 chunks,
        // 4 processes: process1=>chunks 0 to 31/process2=>chunks 32 to 63,..)
        if (id*_nb_procs>=_proc_id*cst::LOC_CHUNK_NUMBER and 
                id*_nb_procs<(_proc_id+1)*cst::LOC_CHUNK_NUMBER) {
            // the blocks are contiguous between chunks, a process begin with
            // block min(beg_block), ends with block max(end_block)
            stream >> beg;
            if (_loc_beg > beg)
                _loc_beg = beg;
            stream >> end;
            if (_loc_end < end)
                _loc_end = end;
            // Sum number of non-zeros over chunks to get the process local nz
            stream >> nz;
            _id.nz_loc += nz;
        }
        std::getline(infile, strInput);
    }
    std:: clog << _id.nz_loc << "\t" << _loc_end << "\t" << _loc_beg << "\n";
    infile.close();
}

int Mumps::nz_mapping() {
    // If distribution with mapping
    if (_distr == parm::A_DISTR_FACTO_MAPPING) {
        // Get mapping array in processes other than host
        if (!is_host())
            _id.mapping = new int[_id.nz];
        MPI_Bcast(_id.mapping, _id.nz, MPI_INTEGER, cst::HOST_ID, _mpi_comm);
        // Compute local nz from mapping array
        _id.nz_loc = 0;
        for(int i = 0; i < _id.nz; ++i)
            // If the current value corresponds to the proc id in the mapping 
            // array, it must be read in this process
            if (_id.mapping[i] == _proc_id)
                _id.nz_loc += 1;
    }
    return _id.nz_loc;
}

bool Mumps::take_A_value_loc(int m, int n, int i, bool local) {
    bool res = true;
    if (!local) return res;
    // If personal distribution (per row, per col, ...)
    if ((_distr == parm::A_DISTR_ANALYSIS || _distr == parm::A_DISTR_FACTO) && 
            _loc != cst::EMPTY_FILE) {
        // The fact we get the value depends on which distribution is chosen
        switch (_loc_option) {
            case cst::DISTR_ROW_UNEVEN:
                if (_proc_id == _nb_procs-1)
                    res = m>=_proc_id*_id.n/_nb_procs;
                else
                    res = m>=_proc_id*_id.n/_nb_procs && m<(_proc_id+1)*_id.n/_nb_procs;
                break;
            case cst::DISTR_ROW_BLOCK:
                res = m>=_loc_beg && m<=_loc_end;
                break;
            case cst::DISTR_COL_BLOCK:
                res = n>=_loc_beg && n<=_loc_end;
                break;
            case cst::DISTR_ARROW_BLOCK:
                res = n>=_loc_beg && m>=_loc_beg && (m<=_loc_end || n<=_loc_end);
                break;
            default:
                std::cerr << "Unrecognized option";
                break;
        }
    }
    // In the case of mapping distribution, we take the value only if it has the
    // same id as the process
    else if (_distr == parm::A_DISTR_FACTO_MAPPING && 
            _id.mapping[i] != _proc_id) {
        res = false;
    }
    return res;
}

int Mumps::nz_loc(int nz, bool local) {
    if (!local) return nz;
    // If personal distribution, parse the distribution file
    if ((_distr == parm::A_DISTR_ANALYSIS || _distr == parm::A_DISTR_FACTO) &&
            _loc != cst::EMPTY_FILE) {
        parse_loc_file();
        return _id.nz_loc;
    }
    // Else if distribution in mapping handle the mapping array
    else if (_distr == parm::A_DISTR_FACTO_MAPPING)
        return nz_mapping();
    return nz;
}

void Mumps::get_A_loc() {
    if (_distr == parm::A_CENTRALIZED) return;
    // If a personal distribution was provided or we use a mapping, use it to 
    // read the whole matrix file and save only the specific ones
    if (_loc != cst::EMPTY_FILE || _distr == parm::A_DISTR_FACTO_MAPPING) {
        get_MM(_file_A, _id.n, _id.n, _id.nz_loc, &_id.a_loc, 
            {&_id.irn_loc, &_id.jcn_loc}, _n_present_A, false, true);
    }
    // Else a different file is given for each process with the suffix "_id"
    else
        get_MM(_file_A + std::to_string(static_cast<long long>(_proc_id)), _id.n, _id.n, _id.nz_loc, 
            &_id.a_loc, {&_id.irn_loc, &_id.jcn_loc}, _n_present_A, false, true);
    // Gather global nz from local nz
    MPI_Reduce(&_id.nz_loc, &_id.nz, 1, MPI_INT, MPI_SUM, cst::HOST_ID, _mpi_comm);
    if (is_host())
        std::clog << "All fronts gathered, nz: " << _id.nz << "\n";
}

void Mumps::get_A() {
    // If not matrix distributed (structure+values) before analysis, read entire
    // matrix on host
    if (is_host() && _distr != parm::A_DISTR_ANALYSIS) {
        get_MM(_file_A, _id.n, _id.n, _id.nz, &_id.a, {&_id.irn, &_id.jcn},
            _n_present_A, false);
    }
    // If matrix distributed (structure+values) before analysis/facto, read
    // local part in all processes, even host
    if (_distr == parm::A_DISTR_ANALYSIS)
        get_A_loc();
}

void Mumps::get_A_again() {
    if (_distr == parm::A_DISTR_FACTO_MAPPING || _distr == parm::A_DISTR_FACTO) {
        deallocate_A();
        get_A_loc();
    }
}

void Mumps::get_b() {
    // On host read b in file or initialize with all 1 if empty
    if (is_host()) {
        if (_file_b.compare(cst::EMPTY_FILE))
            get_MM(_file_b, _id.n, _id.n, _id.nz, &_id.rhs, {}, _n_present_b, true);
        else alloc_rhs();
    }
}
    
////////////////////////////////////////////////////
// MATRIX OUTPUT
////////////////////////////////////////////////////
void Mumps::display_A_loc(int n) {
    if (is_host()) display_ass(_id.a_loc, n, {_id.irn_loc, _id.jcn_loc});
}

void Mumps::display_A_loc() {
    if (is_host()) display_ass(_id.a_loc, _id.nz_loc, {_id.irn_loc, _id.jcn_loc});
}

void Mumps::display_A(int n) {
    if (is_host()) display_ass(_id.a, n, {_id.irn, _id.jcn});
}

void Mumps::display_A() {
    if (is_host()) display_ass(_id.a, _id.nz, {_id.irn, _id.jcn});
}

void Mumps::display_b(int n) {
    if (is_host()) display_ass(_id.rhs, n, {});    
}

void Mumps::display_b() {
    if (is_host()) display_ass(_id.rhs, _id.n, {});
}
    
////////////////////////////////////////////////////
// RUNNING THE SOLVER
////////////////////////////////////////////////////
void Mumps::mumps(int job) {
    _id.job = job;
    dmumps_c(&_id);
}

void Mumps::set_opt(int key, int value) {
    if (is_host())
        std::cout << "Setting " << key << " to value " << value << "\n";
    _id.ICNTL(key) = value;
}

void Mumps::set_opt(int key, int value, std::string sol_spec_file) {
    set_opt(key, value);

    if (is_host()) {
        std::ofstream myfile;
        myfile.open(sol_spec_file.c_str(), std::ofstream::app);
        myfile << "Setting " << key << " to value " << value << "\n";
        myfile.close();
    }
}

void Mumps::init() {
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen = 0;
    
    // Initialize MPI
    int ierr = MPI_Init (NULL, NULL);
    ierr = MPI_Comm_size(_mpi_comm, &_nb_procs);
    ierr = MPI_Comm_rank(_mpi_comm, &_proc_id);
    ierr = MPI_Get_processor_name(processor_name, &namelen);
    std::clog << "MPI initialization in Mumps, error code: " << ierr << "\n";
    std::clog << "Running on CPU " << sched_getcpu() << " on node " << processor_name << "\n";
    std::clog << "Proc id " << _proc_id << " on a total of " << _nb_procs << " procs\n\n";
    
    // Initialize MUMPS
    mumps(parm::JOB_INIT);
    
    // Display initialized OpenMP
    int nthreads, tid;
    #pragma omp parallel private(tid)
    {
        /* Obtain and print thread id */
       	tid = omp_get_thread_num();
        std::clog << "Initialisation of OpenMP thread = " << tid << " on cpu " << sched_getcpu() << "\n";
        /* Only master thread does this */
        if (tid == 0) {
            nthreads = omp_get_num_threads();
            std::clog << "Number of OpenMP threads = " << nthreads << "\n";
        }  /* All threads join master thread and terminate */
    }
    
    // Distribution and format
    if (is_host()) {
        std::cout << "Setting distribution and format:\n";
        set_opt(parm::A_FORMAT, _format);
        set_opt(parm::A_DISTRIBUTION, _distr);
    }
    
    // Default options
    if (is_host())
        std::cout << "\nDefault options:\n";
    set_opt(parm::SEQPAR_ANALYSIS, parm::ANAL_PAR);
    set_opt(parm::SYMPERM_PAR, parm::SYMPERM_PARAUTO);
    set_opt(parm::NULL_PIVOT, parm::NULL_PIVOT_YES);
    set_opt(parm::SCALING, parm::SCALE_ANALYSIS);
    if (is_host() && parm::MEMORY_PERCENT_INC != 0)
        set_opt(parm::MEMORY_PERCENT_INC, parm::MEMORY_DEFAULT_PERCENT_INC);
    
    //Test Option
    if (_opt_key != cst::EMPTY_INT_OPT_KEY) {
        if (is_host())
            std::cout << "Benchmark option to test:\n";
        set_opt(_opt_key, _opt_value);
    }
    
    //Test ID
    if (is_host())
        std::clog << "TEST ID: " << _test_id << "\n";
}

void Mumps::analyse() {
    mumps(parm::JOB_ANALYSIS);
    if (is_host() && parm::MEMORY_SIZE_FACTOR != 1) {
        std::clog << "LOWER BOUND OF MEMORY SIZE: " << _id.INFOG(parm::MEMORY_LOWER_BOUND) << "\n";
        set_opt(parm::MEMORY_SIZE,
            std::floor(parm::MEMORY_SIZE_FACTOR*_id.INFOG(parm::MEMORY_LOWER_BOUND)));
    }
}

bool Mumps::get_b_before_facto() {
    return _id.ICNTL(parm::FACTO_ELIM) == parm::FELIM_YES;
}

void Mumps::factorize() {
    mumps(parm::JOB_FACTO);
}

void Mumps::solve() {
    mumps(parm::JOB_SOLVE);
}

void Mumps::metrics() {
    // If distributed matrix, assemble on host
    if (_distr != parm::A_CENTRALIZED) {
        assemble_A();
    }
    // Compute metrics on host
    if (is_host()) {
        alloc_solve_residual();
        _metrics.init_metrics(_id.n, _id.n, _id.nz, _id.a, _id.irn, _id.jcn, 
            _r, _id.rhs);
        _metrics.residual_norm(_rnrm);
        _metrics.residual_orth(_onrm);
        _metrics.A_norm(_anrm);
        _metrics.sol_norm(_xnrm);
        _metrics.b_norm(_bnrm);
    }
}

void Mumps::call() {
    analyse();
    factorize();
    solve();
    metrics();
}

void Mumps::finalize() {
    delete[] _id.rhs;
    if (is_host())
        problem_spec_metrics_output();
    mumps(parm::JOB_END);
    int ierr = MPI_Finalize();
    std::cout << "MPI finalization in Mumps, error code: " << ierr << "\n";
}

////////////////////////////////////////////////////
// (DE)ALLOCATION
////////////////////////////////////////////////////
void Mumps::alloc_solve_residual() {
    // Allocate a residual array equal to the rhs in Mumps structure on host
    if (is_host()) {
        _r = new double[_id.n];
        for(int i = 0; i < _id.n; i++) _r[i] = _id.rhs[i];
    }
}

void Mumps::alloc_rhs() {
    // Initialize rhs in Mumps structure to all 1 on host
    if (is_host()) {
        std::clog << "Allocating RHS to vector with all 1 of size " << _id.n << "\n";
        _id.rhs = new double[_id.n];
        for(int i = 0; i < _id.n; i++) _id.rhs[i] = 1.0;
    }
}

void Mumps::deallocate_A() {
    delete[] _id.irn;
    delete[] _id.jcn;
    delete[] _id.a;
    delete[] _id.irn_loc;
    delete[] _id.jcn_loc;
    delete[] _id.a_loc;
}

////////////////////////////////////////////////////
// OUTPUTS
////////////////////////////////////////////////////
void Mumps::problem_spec_metrics_output() {
    if (is_host()) {
        std::ifstream ifstr(_pb_spec_file);
        bool empty_file = ifstr.peek() == std::ifstream::traits_type::eof();
        ifstr.close();

        std::ofstream ofstr;
        ofstr.open(_pb_spec_file.c_str(), std::ofstream::app);
        if (empty_file)
            ofstr << "test_id\tfile_A\tanrm\tm\tn\tnz\tsymetry\tcond1\tcond2\tfile_b\tbnrm\n";
        ofstr << _test_id << "\t" << _file_A << "\t" << _anrm << "\t" << _id.n << "\t" << _id.n << 
            "\t" << _id.nz << "\t" << _id.INFOG(8) << "\t" << _id.RINFOG(10) << 
            "\t" << _id.RINFOG(11) << "\t" << _file_b << "\t" << _bnrm << "\n";
        ofstr.close();
    }
}

void Mumps::output_metrics_init(std::string file) {
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "test_id\tfile_A\tsolver\t#procs\toption\tvalue\tta\ttf\tts\t"
            "ta_tot\ttf_tot\tts_tot\txnrm\trnrm\tonrm\tSR\tSR1\tSR2\tup_rnrm\t"
            "e_elim_flops\te_ass_flops\telim_flops\toff_diag_piv\tdelayed_piv\t"
            "tiny_piv\tnull_piv\titer_ref\te_max_front_size\t#nodes\t"
            "order_largest_front\t#factors_entries\n";
        myfile.close();
    }
}

void Mumps::output_metrics(std::string file, long long ta, 
        long long tf, long long ts, std::string key,
        std::string value) {
    long long ta_tot=total_time(&ta);
    long long tf_tot=total_time(&tf);
    long long ts_tot=total_time(&ts);
    if (is_host()) {
        std::cout << "\n#procs     =  " << _nb_procs << "\n" <<
            "current option        =  " << key << "\n" <<
            "current value         =  " << value << "\n" <<
            "time for analysis     =  " << cst::TIME_RATIO*ta << "\n" <<
            "time for facto        =  " << cst::TIME_RATIO*tf << "\n" <<
            "time for solve        =  " << cst::TIME_RATIO*ts << "\n" <<
            "tot time for analysis =  " << cst::TIME_RATIO*ta_tot << "\n" <<
            "tot time for facto    =  " << cst::TIME_RATIO*tf_tot << "\n" <<
            "tot time for solve    =  " << cst::TIME_RATIO*ts_tot << "\n" <<
            "||A||                 =  " << _anrm << "\n" <<
            "||b||                 =  " << _bnrm << "\n" <<
            "||x||                 =  " << _xnrm << "\n" <<
            "||r||/||A||           =  " << _rnrm << "\n" <<
            "||A^tr||/||r||        =  " << _onrm << "\n" <<
            "scaled residual       = " << _id.RINFOG(6) << "\n" <<
            "componentwise SR1     = " << _id.RINFOG(7) << "\n" <<
            "componentwise SR2     = " << _id.RINFOG(8) << "\n" <<
            "upper bound sol error = " << _id.RINFOG(9) << "\n" <<
            "estimated elim flops  = " << _id.RINFOG(1) << "\n" <<
            "estimated assem flops = " << _id.RINFOG(2) << "\n" <<
            "elimination flops     = " << _id.RINFOG(3) << "\n" <<
            "off diag pivots       = " << _id.INFOG(12) << "\n" <<
            "delayed pivots        = " << _id.INFOG(13) << "\n" <<
            "tiny pivots           = " << _id.INFOG(25) << "\n" <<
            "null pivots           = " << _id.INFOG(28) << "\n" <<
            "iterative refinement  = " << _id.INFOG(15) << "\n" <<
            "est. max front size   = " << _id.INFOG(5) << "\n" <<
            "#nodes in elim tree = " << _id.INFOG(6) << "\n" <<
            "order largest frontal = " << _id.INFOG(11) << "\n" <<
            "#factors entries    = " << _id.INFOG(29) << "\n" <<
            "\n";

        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << _test_id << "\t" << _file_A << "\tmumps\t" << _nb_procs << 
            "\t" << key << "\t" << value << "\t" <<
            cst::TIME_RATIO*ta << "\t" << cst::TIME_RATIO*tf << "\t" << 
            cst::TIME_RATIO*ts << "\t" << cst::TIME_RATIO*ta_tot << "\t" << 
            cst::TIME_RATIO*tf_tot << "\t" << cst::TIME_RATIO*ts_tot << 
            "\t" << _xnrm << "\t" << _rnrm << "\t" << _onrm << 
            "\t" << _id.RINFOG(6) << "\t" << _id.RINFOG(7) << "\t" << 
            _id.RINFOG(8) << "\t" << _id.RINFOG(9) << "\t" << _id.RINFOG(1) << 
            "\t" << _id.RINFOG(2) << "\t" << _id.RINFOG(3) << "\t" << 
            _id.INFOG(12) << "\t" << _id.INFOG(13) << "\t" << _id.INFOG(25) << 
            "\t" << _id.INFOG(28) << "\t" << _id.INFOG(15) << "\t" << 
            _id.INFOG(5) << "\t" << _id.INFOG(6) << "\t" << _id.INFOG(11) <<
            "\t" << _id.INFOG(29) << "\n";
        myfile.close();
    }
}

void Mumps::set_no_output() {
    _id.ICNTL(parm::OUT_ERROR) = -1;
    _id.ICNTL(parm::OUT_DIAGNOSTIC) = -1;
    _id.ICNTL(parm::OUT_GINFO) = -1;
    _id.ICNTL(parm::OUT_LEVEL) = 0;
}
