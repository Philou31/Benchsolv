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