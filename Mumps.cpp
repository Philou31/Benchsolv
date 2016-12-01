//! \file functions.cpp
//! \brief Functions for the project
//! \author filou
//! \version 0.1
//! \date 13/11/2016, 18:42
//!

#include <streambuf>

#include "Mumps.h"
#include "constants.h"
#include "Benchmark.h"

Mumps::Mumps(std::string test_id, std::string file_A, bool n_present_A, 
    std::string file_b, bool n_present_b, int par, int sym, int distr, 
    bool loc, int format, int comm, MPI_Comm mpi_comm, 
    std::string pb_spec_file, int int_opt_key, int int_opt_value):
    Solver(test_id, file_A, n_present_A, file_b, n_present_b), 
    _mpi_comm(mpi_comm), _pb_spec_file(pb_spec_file), _distr(distr), 
    _loc(loc), _format(format), _opt_key(int_opt_key), _opt_value(int_opt_value)
{
    std::cout << "Mumps initialization:\nWorking host: " << par << 
        "\nSymmetry: " << sym << "\nCommunicator: " << comm << "\nFormat: " << 
        format << "\nDistribution: " << distr << "\nDistributed inputs: " << 
        loc << "\n";
    _id.par = par;
    _id.sym = sym;
    _id.comm_fortran = comm;
    _id.nrhs = 1;
    init();
    get_A();
    get_b();
}

Mumps::~Mumps() {
    finalize_solve();
    finalize();   
}
    
bool Mumps::is_host() {
    return _proc_id == parm::HOST_ID;
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

void Mumps::get_simple() {
    if (is_host())
        Solver::get_simple(_id.n, _id.n, _id.nz, &_id.a, &_id.irn, &_id.jcn, 
            _id.nrhs, _id.n, &_id.rhs);
}

void Mumps::get_A() {
    if (_loc) {
        get_MM(_file_A + std::to_string(_proc_id), _id.n, _id.n, _id.nz_loc, 
            &_id.a_loc, {&_id.irn_loc, &_id.jcn_loc}, _n_present_A, false);
        MPI_Reduce(&_id.nz_loc, &_id.nz, 1, MPI_INT, MPI_SUM, cst::HOST_ID, _mpi_comm);
        if (is_host())
            std::clog << "All fronts gathered, nz: " << _id.nz << "\n";
    } else
        if (is_host())
            get_MM(_file_A, _id.n, _id.n, _id.nz, &_id.a, {&_id.irn, &_id.jcn}, 
                _n_present_A, false);
}

void Mumps::get_b() {
    if (is_host()) {
        if (_file_b.compare(cst::EMPTY_FILE))
            get_MM(_file_b, _id.n, _id.n, _id.nz, &_id.rhs, {}, _n_present_b, true);
        else alloc_rhs();
    }
}

void Mumps::display_A(int n) {
    if (is_host()) display_ass(_id.a, n, {_id.irn, _id.jcn});
}

void Mumps::display_A_loc(int n) {
    if (is_host()) display_ass(_id.a_loc, n, {_id.irn_loc, _id.jcn_loc});
}

void Mumps::display_b(int n) {
    if (is_host()) display_ass(_id.rhs, n, {});    
}

void Mumps::display_A() {
    if (is_host()) display_ass(_id.a, _id.nz, {_id.irn, _id.jcn});
}

void Mumps::display_A_loc() {
    if (is_host()) display_ass(_id.a_loc, _id.nz_loc, {_id.irn_loc, _id.jcn_loc});
}

void Mumps::display_b() {
    if (is_host()) display_ass(_id.rhs, _id.n, {});
}

void Mumps::display_x(int n) {
    if (is_host()) std::cerr << "display_x: Not accessible with Mumps.\n";
}
void Mumps::display_r(int n) {
    if (is_host()) std::cerr << "display_r: Not accessible with Mumps.\n";
}
void Mumps::display_x() {
    if (is_host()) std::cerr << "display_x: Not accessible with Mumps.\n";
}
void Mumps::display_r() {
    if (is_host()) std::cerr << "display_r: Not accessible with Mumps.\n";
}

void Mumps::mumps(int job) {
    _id.job = job;
    dmumps_c(&_id);
}

void Mumps::init() {
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen = 0;
    int ierr = MPI_Init (NULL, NULL);
    ierr = MPI_Comm_size(_mpi_comm, &_nb_procs);
    ierr = MPI_Comm_rank(_mpi_comm, &_proc_id);
    ierr = MPI_Get_processor_name(processor_name, &namelen);
    std::clog << "MPI initialization in Mumps, error code: " << ierr << "\n";
    std::clog << "Running on CPU " << sched_getcpu() << " on node " << processor_name << "\n";
    std::clog << "Proc id " << _proc_id << " on a total of " << _nb_procs << " procs\n\n";
    mumps(parm::JOB_INIT);
    
    if (is_host()) {
        std::cout << "Setting distribution and format:\n";
        set_opt(parm::A_FORMAT, _format);
        set_opt(parm::A_DISTRIBUTION, _distr);
    }
    
    if (is_host())
        std::cout << "\nDefault options:\n";
    set_opt(parm::SEQPAR_ANALYSIS, parm::ANAL_SEQ);
    set_opt(parm::SYMPERM_SEQ, parm::SYMPERM_SCOTCH);
    set_opt(parm::ERR_ANALYSIS, parm::ERRANAL_FULL);
    set_opt(parm::NULL_PIVOT, parm::NULL_PIVOT_YES);
    set_opt(parm::SCALING, parm::SCALE_ANALYSIS);
    
    //Test Parameters
    if (_opt_key != cst::EMPTY_INT_OPT_KEY) {
        if (is_host())
            std::cout << "Benchmark option to test:\n";
        set_opt(_opt_key, _opt_value);
    }
    
    //Test ID
    if (is_host())
        std::clog << "TEST ID: " << _test_id << "\n";
}

long long Mumps::total_time(long long *t) {
    long long t_tot = 0;
    MPI_Reduce(t, &t_tot, 1, MPI_LONG_LONG, MPI_SUM, cst::HOST_ID, _mpi_comm);
    return t_tot;
}

void Mumps::analyse() {
    mumps(parm::JOB_ANALYSIS);
}

void Mumps::factorize() {
    mumps(parm::JOB_FACTO);
}

void Mumps::alloc_solve_residual() {
    if (is_host()) {
        _r = new double[_id.n];
        for(int i = 0; i < _id.n; i++) _r[i] = _id.rhs[i];
    }
}

void Mumps::alloc_rhs() {
    if (is_host()) {
        std::clog << "Allocating RHS to vector with all 1 of size " << _id.n << "\n";
        _id.rhs = new double[_id.n];
        for(int i = 0; i < _id.n; i++) _id.rhs[i] = 1.0;
    }
}

void Mumps::solve() {
    mumps(parm::JOB_SOLVE);
}

void Mumps::assemble_A() {
    int GO = 42;
    MPI_Status status;
    if (is_host()) {
        _id.a = new double[_id.nz];
        _id.irn = new int[_id.nz];
        _id.jcn = new int[_id.nz];
        for (int i=0; i < _id.nz_loc; ++i) {
            _id.a[i] = _id.a_loc[i];
            _id.irn[i] = _id.irn_loc[i];
            _id.jcn[i] = _id.jcn_loc[i];
        }
        int beg = 0;
        int size = 0;
        for (int i=0; i < _nb_procs-1; ++i) {
            MPI_Send(&GO, 1, MPI_INT, _proc_id+i+1, cst::TAG_GO, _mpi_comm);
            MPI_Recv(&size, 1, MPI_INT, _proc_id+i+1, cst::TAG_SIZE,
                _mpi_comm, &status);
            MPI_Recv(_id.a+_id.nz_loc+beg, size, MPI_DOUBLE, _proc_id+i+1,
                cst::TAG_ARRAY, _mpi_comm, &status);
            MPI_Recv(_id.irn+_id.nz_loc+beg, size, MPI_INT, _proc_id+i+1,
                cst::TAG_IRN, _mpi_comm, &status);
            MPI_Recv(_id.jcn+_id.nz_loc+beg, size, MPI_INT, _proc_id+i+1,
                cst::TAG_JCN, _mpi_comm, &status);
            beg += size;
        }
    } else {
        MPI_Recv(&GO, 1, MPI_INT, cst::HOST_ID, cst::TAG_GO, _mpi_comm, 
            &status);
        MPI_Send(&_id.nz_loc, 1, MPI_INT, cst::HOST_ID, cst::TAG_SIZE, 
            _mpi_comm);
        MPI_Send(_id.a_loc, _id.nz_loc, MPI_DOUBLE, cst::HOST_ID, 
            cst::TAG_ARRAY, _mpi_comm);
        MPI_Send(_id.irn_loc, _id.nz_loc, MPI_INT, cst::HOST_ID, 
            cst::TAG_IRN, _mpi_comm);
        MPI_Send(_id.jcn_loc, _id.nz_loc, MPI_INT, cst::HOST_ID, 
            cst::TAG_JCN, _mpi_comm);
    }
    MPI_Barrier(_mpi_comm);
}

void Mumps::metrics() {
    if (_loc) {
        assemble_A();
    }
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

bool Mumps::get_b_before_facto() {
    return _id.ICNTL(parm::FACTO_ELIM) == parm::FELIM_YES;
}

void Mumps::set_no_output() {
    _id.ICNTL(parm::OUT_ERROR) = -1;
    _id.ICNTL(parm::OUT_DIAGNOSTIC) = -1;
    _id.ICNTL(parm::OUT_GINFO) = -1;
    _id.ICNTL(parm::OUT_LEVEL) = 0;
}

void Mumps::finalize_solve() {
    delete[] _id.irn;
    delete[] _id.jcn;
    delete[] _id.a;
    delete[] _id.irn_loc;
    delete[] _id.jcn_loc;
    delete[] _id.a_loc;
    delete[] _id.rhs;
}

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

void Mumps::finalize() {
    if (is_host())
        problem_spec_metrics_output();
    mumps(parm::JOB_END);
    int ierr = MPI_Finalize();
    std::cout << "MPI finalization in Mumps, error code: " << ierr << "\n";
}

void Mumps::output_metrics_init(std::string file) {
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "test_id\tfile_A\tsolver\t#procs\tta\ttf\tts\tta_tot\ttf_tot\tts_tot"
            "\txnrm\trnrm\tonrm\tSR\tSR1\tSR2\tup_rnrm\t"
            "e_elim_flops\te_ass_flops\telim_flops\toff_diag_piv\tdelayed_piv\t"
            "tiny_piv\tnull_piv\titer_ref\te_max_front_size\t#nodes\t"
            "order_largest_front\t#factors_entries\n";
        myfile.close();
    }
}

void Mumps::output_metrics(std::string file, long long ta, 
        long long tf, long long ts, long long ta_tot, 
        long long tf_tot, long long ts_tot) {
    if (is_host()) {
        std::cout << "\n#procs     =  " << _nb_procs << "\n" <<
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
        myfile << _test_id << "\t" << _file_A << "\tmumps\t" << _nb_procs << "\t" << 
        cst::TIME_RATIO*ta << "\t" << cst::TIME_RATIO*tf << "\t" << 
        cst::TIME_RATIO*ts << "\t" << cst::TIME_RATIO*ta_tot << "\t" << 
        cst::TIME_RATIO*tf << "\t" << cst::TIME_RATIO*ts_tot << 
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
