/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Metrics.h
 * Author: filou
 *
 * Created on 15 novembre 2016, 23:33
 */

#ifndef METRICS_H
#define METRICS_H

#include "QR_Mumps.h"

class Metrics {
private:
    dqrm_spmat_type_c _qrm_mat;
    double *_r, *_x;
public:
    Metrics();
    ~Metrics();
    void init_metrics(int m, int n, int nz, double *values, int *irn,
        int*jcn, double *x, double *r);
    void A_norm(double &anrm);
    void b_norm(double &bnrm);
    void sol_norm(double &xnrm);
    void residual_norm(double &rnrm);
    void residual_orth(double &onrm);
};

#endif /* METRICS_H */

