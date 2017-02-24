//! \file functions.cpp
//! \brief Functions for the project
//! \author filou
//! \version 0.1
//! \date 13/11/2016, 18:42
//!

#include "QR_Mumps.h"
#include "constants.h"

QR_Mumps::QR_Mumps(std::string test_id, std::string file_A, bool n_present_A, 
    std::string file_b, bool n_present_b, std::string string_opt_key, 
    int int_opt_value, int nrows, int ncols, int nz):
    Solver(test_id, file_A, n_present_A, file_b, n_present_b, nrows, ncols, nz),
    _opt_key(string_opt_key), _opt_value(int_opt_value)
{
    std::cout << "QR_Mumps initialization\n";
    init();
    get_A();
    get_b();
}

QR_Mumps::~QR_Mumps() {
    deallocate_A();
    deallocate_b();
    finalize();   
}

////////////////////////////////////////////////////
// MATRIX INPUT
////////////////////////////////////////////////////    
void QR_Mumps::get_A() {
    get_MM(_file_A, _qrm_mat.m, _qrm_mat.n, _qrm_mat.nz, &_qrm_mat.val, 
        {&_qrm_mat.irn, &_qrm_mat.jcn}, _n_present_A, false);
    // Matrix is transposed only if more columns than lines
    if(_qrm_mat.m < _qrm_mat.n) _transp = 't';
    else _transp = 'n';
}

void QR_Mumps::get_b() {
    if (is_host()) {
        // If the file for b is empty, initialize b with all 1
        if (_file_b.compare(cst::EMPTY_FILE)) {
            get_MM(_file_b, _qrm_mat.m, _qrm_mat.n, _qrm_mat.nz, &_b, {}, _n_present_b,
                true);
        } else alloc_rhs();
        alloc_solve_residual();
    }
}
    
////////////////////////////////////////////////////
// MATRIX OUTPUT
////////////////////////////////////////////////////
void QR_Mumps::display_A(int n) {
    display_ass(_qrm_mat.val, n, {_qrm_mat.irn, _qrm_mat.jcn});
}

void QR_Mumps::display_A() {
    display_ass(_qrm_mat.val, _qrm_mat.nz, {_qrm_mat.irn, _qrm_mat.jcn});
}

void QR_Mumps::display_x(int n) {
    display_ass(_x, n, {});
}

void QR_Mumps::display_x() {
    display_ass(_x, _qrm_mat.m, {});
}

void QR_Mumps::display_b(int n) {
    display_ass(_b, n, {});
}

void QR_Mumps::display_b() {
    display_ass(_b, _qrm_mat.m, {});    
}

void QR_Mumps::display_r(int n) {
    display_ass(_r, n, {});
}

void QR_Mumps::display_r() {
    display_ass(_r, _qrm_mat.m, {});
}
    
////////////////////////////////////////////////////
// RUNNING THE SOLVER
////////////////////////////////////////////////////
void QR_Mumps::set_opt(std::string key, int value) {
    std::cout << "Setting " << key << " to value " << value << "\n";
    // The output units are set with qrm_gseti_c, others with dqrm_pseti_c
    if (!key.compare(parqrm::EUNIT) || !key.compare(parqrm::OUNIT) ||
        !key.compare(parqrm::DUNIT)) {
        qrm_gseti_c(key.c_str(), value);
    } else dqrm_pseti_c(&_qrm_mat, key.c_str(), value);
}

void QR_Mumps::set_opt(std::string key, int value, std::string sol_spec_file) {
    set_opt(key, value);

    std::ofstream myfile;
    myfile.open(sol_spec_file.c_str(), std::ofstream::app);
    myfile << "Setting " << key << " to value " << value << "\n";
    myfile.close();
}

void QR_Mumps::init() {
    qrm_init_c(0);
    dqrm_spmat_init_c(&_qrm_mat);
    
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
    
    //Default parameters
    //Global
    if (is_host())
        std::cout << "\nDefault options:\n";
    set_opt(parqrm::DUNIT.c_str(), 1);
    set_opt(parqrm::EUNIT.c_str(), 2);
    set_opt(parqrm::OUNIT.c_str(), 3);
    //Local
    set_opt(parqrm::ORDERING.c_str(), qrm_scotch_);
    set_opt(parqrm::KEEPH.c_str(), qrm_yes_);
    
    //Set Test Parameters
    if (_opt_key != cst::EMPTY_STRING_OPT_KEY && 
            _opt_value != cst::EMPTY_INT_OPT_VALUE) {
        if (is_host())
            std::cout << "Benchmark option to test:\n";
        set_opt(_opt_key, _opt_value);
    }
    
    //Test ID
    if (is_host())
        std::clog << "TEST ID: " << _test_id << "\n";
}

void QR_Mumps::analyse() {
    dqrm_analyse_c(&_qrm_mat, _transp);
}

bool QR_Mumps::get_b_before_facto() {
    return false;
}

void QR_Mumps::factorize() {
    dqrm_factorize_c(&_qrm_mat, _transp);
}

void QR_Mumps::solve() {
    // Launched differently if transposed matrix
    if(_transp == 'n'){
        dqrm_apply_c(&_qrm_mat, 't', _b, _nrhs);
        dqrm_solve_c(&_qrm_mat, 'n', _b, _x, _nrhs);
    } else if (_transp == 't') {
        dqrm_solve_c(&_qrm_mat, 't', _b, _x, _nrhs);
        dqrm_apply_c(&_qrm_mat, 'n', _x, _nrhs);
    }
}

void QR_Mumps::metrics() {
    dqrm_residual_norm_c(&_qrm_mat, _r, _x, 1, &_rnrm);
    dqrm_residual_orth_c(&_qrm_mat, _r, 1, &_onrm);
    dqrm_vecnrm_c(_x, _qrm_mat.n, 1, '2', &_xnrm);
    dqrm_vecnrm_c(_b, _qrm_mat.m, 1, '2', &_bnrm);
    dqrm_matnrm_c(&_qrm_mat, 'f', &_anrm);
}

void QR_Mumps::call() {
    analyse();
    factorize();
    solve();
    metrics();
}

void QR_Mumps::finalize() {
    qrm_finalize_c();
    delete[] _r;
    delete[] _x;
}
    
////////////////////////////////////////////////////
// (DE)ALLOCATION
////////////////////////////////////////////////////
void QR_Mumps::alloc_rhs() {
    // Initialize rhs on host with all 1
    if (is_host()) {
        _b = new double[_qrm_mat.m];
        int i;
        for(i = 0; i < _qrm_mat.m; i++) _b[i] = (double)1.0;
    }
}

void QR_Mumps::alloc_solve_residual() {
    // Initialized residual array
    _r = new double[_qrm_mat.m];
    for(int i = 0; i < _qrm_mat.m; i++) _r[i] = _b[i];
    // Initialized solution array
    _x = new double[_qrm_mat.n];
    for(int i = 0; i < _qrm_mat.n; i++) _x[i] = (double)0.0;
}

void QR_Mumps::deallocate_A() {
    delete[] _qrm_mat.irn;
    delete[] _qrm_mat.jcn;
    delete[] _qrm_mat.val;
    dqrm_spmat_destroy_c(&_qrm_mat);
}

void QR_Mumps::deallocate_b() {
    delete[] _b;
}

////////////////////////////////////////////////////
// OUTPUTS
////////////////////////////////////////////////////
void QR_Mumps::output_metrics_init(std::string file) {
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "test_id\tfile_A\tsolver\t#procs\toption\tvalue\tta\ttf\tts\t"
            "ta_tot\ttf_tot\tts_tot\txnrm\trnrm\tonrm\tnon0_r\tnon0_h\te_non0_r"
            "\te_non0_h\tfacto_flops\te_mempeak\n";
        myfile.close();
    }
}
    
void QR_Mumps::output_metrics(std::string sol_spec_file, long long ta, 
        long long tf, long long ts, std::string key,
        std::string value) {
    if (is_host()) {
        std::cout << "\ncurrent option    =  " << key << "\n" <<
            "current value     =  " << value << "\n" <<
            "time for analysis =  " << ta << "\n" <<
            "time for facto    =  " << cst::TIME_RATIO*tf << "\n" <<
            "time for solve    =  " << cst::TIME_RATIO*ts << "\n" <<
            "||A||             =  " << _anrm << "\n" <<
            "||b||             =  " << _bnrm << "\n" <<
            "||x||             =  " << _xnrm << "\n" <<
            "||r||/||A||       =  " << _rnrm << "\n" <<
            "||A^tr||/||r||    =  " << _onrm << "\n" <<
            "Nonzeros in R                 : " << _qrm_mat.gstats[qrm_nnz_r_] << "\n" <<
            "Nonzeros in H                 : " << _qrm_mat.gstats[qrm_nnz_h_] << "\n" <<
            "Estimated nonzeros in R       : " << _qrm_mat.gstats[qrm_e_nnz_r_] << "\n" <<
            "Estimated nonzeros in H       : " << _qrm_mat.gstats[qrm_e_nnz_h_] << "\n" <<
            "Total flops at facto          : " << _qrm_mat.gstats[qrm_e_facto_flops_] << "\n" <<
            "Estimated Memory Peak at facto: " << _qrm_mat.gstats[qem_e_facto_mempeak_] << "\n" <<
            "\n";
    
        std::ofstream myfile;
        myfile.open(sol_spec_file.c_str(), std::ofstream::app);
        myfile << _test_id << "\t" << _file_A << "\tqr_mumps\t" << 1 << "\t" <<
            key << "\t" << value << "\t" << 
            cst::TIME_RATIO*ta << "\t" << cst::TIME_RATIO*tf << "\t" << 
            cst::TIME_RATIO*ts << "\t" <<       
            _xnrm << "\t" << _rnrm << "\t" << _onrm << 
            "\t" << _qrm_mat.gstats[qrm_nnz_r_] << "\t" << 
            _qrm_mat.gstats[qrm_nnz_h_] << "\t" << 
            _qrm_mat.gstats[qrm_e_nnz_r_] << "\t" << 
            _qrm_mat.gstats[qrm_e_nnz_h_] << "\t" << 
            _qrm_mat.gstats[qrm_e_facto_flops_] << "\t" << 
            _qrm_mat.gstats[qem_e_facto_mempeak_] << "\n";
        myfile.close();
    }
}

void QR_Mumps::set_no_output() {
    set_opt(parqrm::DUNIT.c_str(), -1);
    set_opt(parqrm::OUNIT.c_str(), -1);
    set_opt(parqrm::EUNIT.c_str(), -1);
}
