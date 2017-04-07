//! \file QR_Mumps.h
//! \brief Mumps class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class QR_Mumps inheriting from/implementing methods of class Solver
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
    base_construct();
}

QR_Mumps::~QR_Mumps() {
    base_destruct();
}

////////////////////////////////////////////////////
// MATRIX INPUT
////////////////////////////////////////////////////    
void QR_Mumps::read_A() {
    read_MM(_file_A, _id.m, _id.n, _id.nz, &_id.val, 
        {&_id.irn, &_id.jcn}, _n_present_A, false);
    // Matrix is transposed only if more columns than lines
    if(_id.m < _id.n) _transp = 't';
    else _transp = 'n';
    _A_assembled=true;
}

void QR_Mumps::read_b() {
    if (is_host()) {
        // If the file for b is empty, initialize b with all 1
        if (_file_b.compare(cst::EMPTY_FILE)) {
            read_MM(_file_b, _id.m, _id.n, _id.nz, &_b, {}, _n_present_b,
                true);
        } else alloc_rhs();
        alloc_solve_residual();
    }
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
    } else dqrm_pseti_c(&_id, key.c_str(), value);
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
    dqrm_spmat_init_c(&_id);
    
    init_OpenMP();
    
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
    dqrm_analyse_c(&_id, _transp);
}

void QR_Mumps::factorize() {
    dqrm_factorize_c(&_id, _transp);
}

void QR_Mumps::solve() {
    // Launched differently if transposed matrix
    if(_transp == 'n'){
        dqrm_apply_c(&_id, 't', _b, _nrhs);
        dqrm_solve_c(&_id, 'n', _b, _x, _nrhs);
    } else if (_transp == 't') {
        dqrm_solve_c(&_id, 't', _b, _x, _nrhs);
        dqrm_apply_c(&_id, 'n', _x, _nrhs);
    }
}

void QR_Mumps::metrics() {
    dqrm_residual_norm_c(&_id, _r, _x, 1, &_rnrm);
    dqrm_residual_orth_c(&_id, _r, 1, &_onrm);
    dqrm_vecnrm_c(_x, _id.n, 1, '2', &_xnrm);
    dqrm_vecnrm_c(_b, _id.m, 1, '2', &_bnrm);
    dqrm_matnrm_c(&_id, 'f', &_anrm);
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
    alloc_array(_id.m, &_b);
}

void QR_Mumps::alloc_solve_residual() {
    alloc_array(_id.m, &_r, true, _b);
    alloc_array(_id.m, &_x, false, NULL, 0.0);
}

void QR_Mumps::deallocate_A() {
    delete[] _id.irn;
    delete[] _id.jcn;
    delete[] _id.val;
    dqrm_spmat_destroy_c(&_id);
}

void QR_Mumps::deallocate_b() {
    delete[] _b;
}

////////////////////////////////////////////////////
// OUTPUTS
////////////////////////////////////////////////////
void QR_Mumps::output_metrics_init(std::string file) {
    base_output_metrics_init(file);
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "\tnon0_r\tnon0_h\te_non0_r\te_non0_h\tfacto_flops\t"
            "e_mempeak\n";
        myfile.close();
    }
}
    
void QR_Mumps::output_metrics(std::string file, long long ta, 
        long long tf, long long ts, std::string key,
        std::string value) {
    base_output_metrics(file, MPI_COMM_WORLD, ta, tf, ts);
    if (is_host()) {
        std::cout << "current option        =  " << key << "\n" <<
            "current value         =  " << value << "\n" <<
            "Nonzeros in R                 : " << _id.gstats[qrm_nnz_r_] << "\n" <<
            "Nonzeros in H                 : " << _id.gstats[qrm_nnz_h_] << "\n" <<
            "Estimated nonzeros in R       : " << _id.gstats[qrm_e_nnz_r_] << "\n" <<
            "Estimated nonzeros in H       : " << _id.gstats[qrm_e_nnz_h_] << "\n" <<
            "Total flops at facto          : " << _id.gstats[qrm_e_facto_flops_] << "\n" <<
            "Estimated Memory Peak at facto: " << _id.gstats[qem_e_facto_mempeak_] << "\n" <<
            "\n";
    
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << key << "\t" << value << "\t" << _id.gstats[qrm_nnz_r_] << 
            "\t" << _id.gstats[qrm_nnz_h_] << "\t" << 
            _id.gstats[qrm_e_nnz_r_] << "\t" << 
            _id.gstats[qrm_e_nnz_h_] << "\t" << 
            _id.gstats[qrm_e_facto_flops_] << "\t" << 
            _id.gstats[qem_e_facto_mempeak_] << "\n";
        myfile.close();
    }
}

void QR_Mumps::set_no_output() {
    set_opt(parqrm::DUNIT.c_str(), -1);
    set_opt(parqrm::OUNIT.c_str(), -1);
    set_opt(parqrm::EUNIT.c_str(), -1);
}

////////////////////////////////////////////////////
// GETTERS
////////////////////////////////////////////////////
int QR_Mumps::get_m() {
    return _id.m;
}

int QR_Mumps::get_n() {
    return _id.n;
}

int QR_Mumps::get_nz() {
    return _id.nz;
}

int QR_Mumps::get_nz_loc() {
    return _id.nz;
}

double* QR_Mumps::get_a() {
    return _id.val;
}

double* QR_Mumps::get_a_loc() {
    return _id.val;
}

int* QR_Mumps::get_irn() {
    return _id.irn;
}

int* QR_Mumps::get_irn_loc() {
    return _id.irn;
}

int* QR_Mumps::get_jcn() {
    return _id.jcn;
}

int* QR_Mumps::get_jcn_loc() {
    return _id.jcn;
}

double* QR_Mumps::get_rhs() {
    return _b;
}

double* QR_Mumps::get_x() {
    return _x;
}

double* QR_Mumps::get_r() {
    return _r;
}
