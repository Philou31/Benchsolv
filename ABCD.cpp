////! \file ABCD.cpp
////! \brief Functions for the project
////! \author filou
////! \version 0.1
////! \date 13/11/2016, 18:42
////!
//
//#include "ABCD.h"
//#include "Mumps.h"
//
//ABCD::ABCD(std::string test_id, std::string file_A, bool n_present_A, 
//        std::string file_b, bool n_present_b, int int_opt_key, 
//        int int_opt_value, int nrows, int ncols, int nz, bool sym):
//    Solver(test_id, file_A, n_present_A, file_b, n_present_b, nrows, ncols, nz)
//{
//    std::cout << "ABCD initialization\n";
//    _id.sym = sym;
//    _id.start_index = 0;
//    _id.nrhs = 1;
//    base_construct();
//}
//
//ABCD::~ABCD() {
//    base_destruct();
//}
//    
//////////////////////////////////////////////////////
//// MATRIX INPUT
//////////////////////////////////////////////////////
//void ABCD::read_A() {
//    if (is_host())
//        read_MM(_file_A, _id.m, _id.n, _id.nz, &_id.val, {&_id.irn, &_id.jcn},
//            _n_present_A, false);
//}
//
//void ABCD::read_b() {
//    // On host read b in file or initialize with all 1 if empty
//    if (is_host()) {
//        if (_file_b.compare(cst::EMPTY_FILE))
//            read_MM(_file_b, _id.m, _id.n, _id.nz, &_id.rhs, {}, _n_present_b, true);
//        else alloc_rhs();
//    }
//}
//    
//////////////////////////////////////////////////////
//// RUNNING THE SOLVER
//////////////////////////////////////////////////////
//void ABCD::set_opt(int key, int value) {
//    if (is_host())
//        std::cout << "Setting " << key << " to value " << value << "\n";
//    _id.icntl[key] = value;
//}
//
//void ABCD::set_opt(int key, int value, std::string sol_spec_file) {
//    set_opt(key, value);
//
//    if (is_host()) {
//        std::ofstream myfile;
//        myfile.open(sol_spec_file.c_str(), std::ofstream::app);
//        myfile << "Setting " << key << " to value " << value << "\n";
//        myfile.close();
//    }
//}
//
//void ABCD::init() {
//    init_MPI(_mpi_comm);
//    
//    // Initialize ABCD
//    _id(parab::JOB_INIT);
//    
//    // Display initialized OpenMP
//    #pragma omp parallel
//    {
//        /* Obtain and print thread id */
//       	_tid = omp_get_thread_num();
//        std::clog << "Initialisation of OpenMP thread = " << _tid << " on cpu " << sched_getcpu() << "\n";
//        /* Only master thread does this */
//        _nthreads = omp_get_num_threads();
//        if (_tid == 0) {
//            std::clog << "Number of OpenMP threads = " << _nthreads << "\n";
//        }  /* All threads join master thread and terminate */
//    }
//    
//    //Test Option
//    if (_opt_key != cst::EMPTY_INT_OPT_KEY) {
//        if (is_host())
//            std::cout << "Benchmark option to test:\n";
//        set_opt(_opt_key, _opt_value);
//    }
//    
//    //Test ID
//    if (is_host())
//        std::clog << "TEST ID: " << _test_id << "\n";
//}
//
//void ABCD::analyse() {
//    _id(parab::JOB_ANALYSIS);
//}
//
//void ABCD::factorize() {
//    _id(parab::JOB_FACTO);
//}
//
//void ABCD::solve() {
//    _id(parab::JOB_SOLVE);
//}
//
//void ABCD::metrics() {
//    // Compute metrics on host
//    if (is_host()) {
//        alloc_solve_residual();
//        _metrics.init_metrics(_id.n, _id.n, _id.nz, _id.val, _id.irn, _id.jcn, 
//            _id.sol, _r);
//        _metrics.compute_metrics(_rnrm, _onrm, _anrm, _xnrm, _bnrm);
//        delete[] _r;
//    }
//}
//
//void ABCD::finalize() {
//    _id(parab::JOB_END);
//    int ierr = MPI_Finalize();
//    std::cout << "MPI finalization in Mumps, error code: " << ierr << "\n";
//}
//
//////////////////////////////////////////////////////
//// (DE)ALLOCATION
//////////////////////////////////////////////////////
//void ABCD::alloc_solve_residual() {
//    alloc_array(_id.n, &_r, true, _id.rhs);
//}
//
//void ABCD::alloc_rhs() {
//    alloc_array(_id.n, &_id.rhs);
//}
//
//void ABCD::deallocate_A() {
//    delete[] _id.irn;
//    delete[] _id.jcn;
//    delete[] _id.val;
//}
//
//void ABCD::deallocate_b() {
//    delete[] _id.rhs;
//}
//
//////////////////////////////////////////////////////
//// OUTPUTS
//////////////////////////////////////////////////////
//void ABCD::output_metrics_init(std::string file) {
//    base_output_metrics_init(file);
//    if (is_host()) {
//        std::ofstream myfile;
//        myfile.open(file.c_str(), std::ofstream::app);
//        myfile << "option\tvalue\n";
//        myfile.close();
//    }
//}
//
//void ABCD::output_metrics(std::string file, long long ta, 
//        long long tf, long long ts, std::string key,
//        std::string value) {
//    base_output_metrics(file, _mpi_comm, ta, tf, ts);
//    if (is_host()) {
//        std::cout << "current option        =  " << key << "\n" <<
//            "current value         =  " << value << "\n" <<
//            "\n";
//
//        std::ofstream myfile;
//        myfile.open(file.c_str(), std::ofstream::app);
//        myfile << key << "\t" << value << "\n";
//        myfile.close();
//    }
//}
//
//void ABCD::set_no_output() {
//    set_opt(parab::VERBOSE, parab::VERBOSE_NO);
//}
//
//////////////////////////////////////////////////////
//// GETTERS
//////////////////////////////////////////////////////
//int ABCD::get_m() {
//    return _id.m;
//}
//
//int ABCD::get_n() {
//    return _id.n;
//}
//
//int ABCD::get_nz() {
//    return _id.nz;
//}
//
//int ABCD::get_nz_loc() {
//    return _id.nz;
//}
//
//double* ABCD::get_a() {
//    return _id.val;
//}
//
//double* ABCD::get_a_loc() {
//    return _id.val;
//}
//
//int* ABCD::get_irn() {
//    return _id.irn;
//}
//
//int* ABCD::get_irn_loc() {
//    return _id.irn;
//}
//
//int* ABCD::get_jcn() {
//    return _id.jcn;
//}
//
//int* ABCD::get_jcn_loc() {
//    return _id.jcn;
//}
//
//double* ABCD::get_rhs() {
//    return _id.rhs;
//}
//
//double* ABCD::get_x() {
//    return _id.rhs;
//}
//
//double* ABCD::get_r() {
//    return _r;
//}
//
