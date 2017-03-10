//! \file Mumps.h
//! \brief Mumps class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class Mumps inheriting from/implementing methods of class Solver
//!

#include "Mumps.h"
#include <typeinfo>

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
    // Call init and get matrices
    base_construct();
    if (is_host())
        std::clog << "Error Analysis\n";
    set_opt(parm::ERR_ANALYSIS, erranal);
}

Mumps::~Mumps() {
    base_destruct();
    if (_distr != parm::A_CENTRALIZED)
        deallocate_A_loc();
}
//
//////////////////////////////////////////////////////
//// MPI communication methods
//////////////////////////////////////////////////////
void Mumps::assemble_A() {
//    TODO: make it work
//    assemble_MM(_id.nz, &_id.a, {&_id.irn, &_id.jcn}, _id.nz_loc, &_id.a_loc, 
//        {&_id.irn_loc, &_id.jcn_loc}, _mpi_comm);
    _id.irn=NULL;
    _id.jcn=NULL;
    _id.a=NULL;
    // Gather all local nz in an array
    int *recvcounts=NULL;
    if (is_host())
        recvcounts = new int[_nb_procs];
    MPI_Gather(&_id.nz_loc, 1, MPI_INT, recvcounts, 1, MPI_INT, cst::HOST_ID,
            MPI_COMM_WORLD);
    // The displacements of the gather are the same as the receive counts
    int *displs=NULL;
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
    _A_assembled=true;
    if (is_host()) {
        delete[] recvcounts;
        delete[] displs;
    }
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

int Mumps::read_nz_loc(int nz, bool local) {
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

void Mumps::read_A_loc() {
    if (_distr == parm::A_CENTRALIZED) return;
    // If a personal distribution was provided or we use a mapping, use it to 
    // read the whole matrix file and save only the specific ones
    if (_loc != cst::EMPTY_FILE || _distr == parm::A_DISTR_FACTO_MAPPING) {
        read_MM(_file_A, _id.n, _id.n, _id.nz_loc, &_id.a_loc, 
            {&_id.irn_loc, &_id.jcn_loc}, _n_present_A, false, true);
    }
    // Else a different file is given for each process with the suffix "_id"
    else
        read_MM(_file_A + std::to_string(static_cast<long long>(_proc_id)), _id.n, _id.n, _id.nz_loc, 
            &_id.a_loc, {&_id.irn_loc, &_id.jcn_loc}, _n_present_A, false, true);
    // Gather global nz from local nz
    MPI_Reduce(&_id.nz_loc, &_id.nz, 1, MPI_INT, MPI_SUM, cst::HOST_ID, _mpi_comm);
    if (is_host())
        std::clog << "All fronts gathered, nz: " << _id.nz << "\n";
}

void Mumps::read_A() {
    // If not matrix distributed (structure+values) before analysis, read entire
    // matrix on host
    if (_distr != parm::A_DISTR_ANALYSIS) {
        if (is_host())
            read_MM(_file_A, _id.n, _id.n, _id.nz, &_id.a, {&_id.irn, &_id.jcn},
                _n_present_A, false);
            _A_assembled=true;
    }
    // If matrix distributed (structure+values) before analysi, read local part 
    // in all processes, even host
    else read_A_loc();
}

void Mumps::read_A_again() {
    if (_distr == parm::A_DISTR_FACTO_MAPPING || _distr == parm::A_DISTR_FACTO) {
        deallocate_A_loc();
        read_A_loc();
    }
}

void Mumps::read_b_again() {
    deallocate_b();
    read_b();
}

void Mumps::read_b() {
    // On host read b in file or initialize with all 1 if empty
    if (is_host()) {
        if (_file_b.compare(cst::EMPTY_FILE))
            read_MM(_file_b, _id.n, _id.n, _id.nz, &_id.rhs, {}, _n_present_b, true);
        else alloc_rhs();
    }
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
    // Initialize MPI
    init_MPI(_mpi_comm);
    
    // Initialize MUMPS
    mumps(parm::JOB_INIT);
    
    // Display initialized OpenMP
    init_OpenMP();
    
    // Default options
    if (is_host())
        std::cout << "\nDefault options:\n";
    // Force sequential behaviour if only one process
    if (_nb_procs == 1)
        _distr = parm::A_CENTRALIZED;
        _id.par = parm::WORKING_HOST;
    if (_nb_procs * _nthreads == 1) {
        set_opt(parm::SEQPAR_ANALYSIS, parm::ANAL_SEQ);
        set_opt(parm::SYMPERM_SEQ, parm::SYMPERM_SEQAUTO);
    } else {
        set_opt(parm::SEQPAR_ANALYSIS, parm::ANAL_PAR);
        set_opt(parm::SYMPERM_PAR, parm::SYMPERM_PARAUTO);
    }
    set_opt(parm::NULL_PIVOT, parm::NULL_PIVOT_YES);
    set_opt(parm::SCALING, parm::SCALE_AUTO);
    if (is_host() && parm::MEMORY_PERCENT_INC != 0)
        set_opt(parm::MEMORY_PERCENT_INC, parm::MEMORY_DEFAULT_PERCENT_INC);
    
    //Test Option
    if (_opt_key != cst::EMPTY_INT_OPT_KEY) {
        if (is_host())
            std::cout << "Benchmark option to test:\n";
        set_opt(_opt_key, _opt_value);
    }
    
    // Distribution and format
    if (is_host()) 
//    {
        std::cout << "Setting distribution and format:\n";
        set_opt(parm::A_FORMAT, _format);
        set_opt(parm::A_DISTRIBUTION, _distr);
//    }
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

bool Mumps::read_A_before_facto() {
    return _id.ICNTL(parm::A_DISTRIBUTION) == parm::A_DISTR_FACTO_MAPPING
            || _id.ICNTL(parm::A_DISTRIBUTION) == parm::A_DISTR_FACTO;
}

bool Mumps::read_b_before_facto() {
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
    if (_distr == parm::A_DISTR_ANALYSIS && !_A_assembled) {
        if (is_host())
            deallocate_A();
        assemble_A();
    }
    // Compute metrics on host
    if (is_host()) {
        alloc_solve_residual();
        _metrics.init_metrics(_id.n, _id.n, _id.nz, _id.a, _id.irn, _id.jcn, 
            _id.rhs, _r);
        _metrics.compute_metrics(_rnrm, _onrm, _anrm, _xnrm, _bnrm);
        delete[] _r;
    }
}

void Mumps::finalize() {
    if (is_host())
        problem_spec_metrics_output();
    mumps(parm::JOB_END);
    finalize_MPI();
}

////////////////////////////////////////////////////
// (DE)ALLOCATION
////////////////////////////////////////////////////
void Mumps::alloc_solve_residual() {
    alloc_array(_id.n, &_r, true, _id.rhs);
}

void Mumps::alloc_rhs() {
    alloc_array(_id.n, &_id.rhs);
}

void Mumps::deallocate_A() {
    delete[] _id.irn;
    delete[] _id.jcn;
    delete[] _id.a;
}

void Mumps::deallocate_A_loc() {
    delete[] _id.irn_loc;
    delete[] _id.jcn_loc;
    delete[] _id.a_loc;
}

void Mumps::deallocate_b() {
    delete[] _id.rhs;
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
    base_output_metrics_init(file);
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "SR\tSR1\tSR2\tup_rnrm\t"
            "e_elim_flops\te_ass_flops\telim_flops\toff_diag_piv\tdelayed_piv\t"
            "tiny_piv\tnull_piv\titer_ref\te_max_front_size\t#nodes\t"
            "order_largest_front\t#factors_entries\n";
        myfile.close();
    }
}

void Mumps::output_metrics(std::string file, long long ta, long long tf, 
        long long ts, std::string key, std::string value) {
    base_output_metrics(file, _mpi_comm, ta, tf, ts);
    if (is_host()) {
        std::cout << "current option        =  " << key << "\n" <<
            "current value         =  " << value << "\n" <<
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
        myfile << key << "\t" << value << "\t" << _id.RINFOG(6) << "\t" << 
            _id.RINFOG(7) << "\t" << _id.RINFOG(8) << "\t" << _id.RINFOG(9) << 
            "\t" << _id.RINFOG(1) << "\t" << _id.RINFOG(2) << "\t" << 
            _id.RINFOG(3) << "\t" << _id.INFOG(12) << "\t" << _id.INFOG(13) << 
            "\t" << _id.INFOG(25) << "\t" << _id.INFOG(28) << "\t" << 
            _id.INFOG(15) << "\t" << _id.INFOG(5) << "\t" << _id.INFOG(6) << 
            "\t" << _id.INFOG(11) << "\t" << _id.INFOG(29) << "\n";
        myfile.close();
    }
}

void Mumps::set_no_output() {
    _id.ICNTL(parm::OUT_ERROR) = -1;
    _id.ICNTL(parm::OUT_DIAGNOSTIC) = -1;
    _id.ICNTL(parm::OUT_GINFO) = -1;
    _id.ICNTL(parm::OUT_LEVEL) = 0;
}

////////////////////////////////////////////////////
// GETTERS
////////////////////////////////////////////////////
int Mumps::get_m() {
    return _id.n;
}

int Mumps::get_n() {
    return _id.n;
}

int Mumps::get_nz() {
    return _id.nz;
}

int Mumps::get_nz_loc() {
    return _id.nz_loc;
}

double* Mumps::get_a() {
    return _id.a;
}

double* Mumps::get_a_loc() {
    return _id.a_loc;
}

int* Mumps::get_irn() {
    return _id.irn;
}

int* Mumps::get_irn_loc() {
    return _id.irn_loc;
}

int* Mumps::get_jcn() {
    return _id.jcn;
}

int* Mumps::get_jcn_loc() {
    return _id.jcn_loc;
}

double* Mumps::get_rhs() {
    return _id.rhs;
}

double* Mumps::get_x() {
    return _id.rhs;
}

double* Mumps::get_r() {
    return _r;
}
