//! \file Metrics.h
//! \brief Metrics class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class Metrics computing metrics on the problem and the found solution. This
//! class is meant to be fair in the results gotten from the solvers: metrics 
//! are computed in the exact same way.
//!

#include "Metrics.h"

Metrics::Metrics() {
    qrm_init_c(0);
    dqrm_spmat_init_c(&_qrm_mat);
}

void Metrics::init_metrics(int m, int n, int nz, double *values, int *irn,
        int*jcn, double *x, double *r) {
    _qrm_mat.m = m;
    _qrm_mat.n = n;
    _qrm_mat.nz = nz;
    _qrm_mat.val = values;
    _qrm_mat.irn = irn;
    _qrm_mat.jcn = jcn;
    _x = x;
    _r = r;
}
        
Metrics::~Metrics() {    
    qrm_finalize_c();
}

void Metrics::A_norm(double &anrm) {
    dqrm_matnrm_c(&_qrm_mat, 'f', &anrm);
}

void Metrics::b_norm(double &bnrm) {
    dqrm_vecnrm_c(_r, _qrm_mat.m, 1, '2', &bnrm);
}

void Metrics::sol_norm(double &xnrm) {
    dqrm_vecnrm_c(_x, _qrm_mat.n, 1, '2', &xnrm);
}

void Metrics::residual_norm(double &rnrm) {
    dqrm_residual_norm_c(&_qrm_mat, _r, _x, 1, &rnrm);    
}

void Metrics::residual_orth(double &onrm) {
    dqrm_residual_orth_c(&_qrm_mat, _r, 1, &onrm);    
}

void Metrics::compute_metrics(double &rnrm, double &onrm, double &anrm, 
        double &xnrm, double &bnrm) {
    residual_norm(rnrm);
    residual_orth(onrm);
    A_norm(anrm);
    sol_norm(xnrm);
    b_norm(bnrm);
}