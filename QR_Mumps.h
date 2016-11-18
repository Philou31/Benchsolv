/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   QR_Mumps.h
 * Author: filou
 *
 * Created on 11 novembre 2016, 14:42
 */

#ifndef QR_MUMPS_H
#define QR_MUMPS_H

#include <fstream>
#include <iostream>
#include <sstream>
#include "Solver.h"
extern "C" {
    #include "dqrm_c.h"
}
#include "Parameters_QR_Mumps.h"

class QR_Mumps: public Solver {
private:
    dqrm_spmat_type_c _qrm_mat;
    char _transp{'n'};
    double *_x, *_r, *_b;
    int _nrhs{1};
public:
    QR_Mumps(std::string file_A, bool n_present_A, std::string file_b, 
        bool n_present_b);
    ~QR_Mumps();
    virtual bool is_host() override;
    void set_opt(std::string key, int value, std::string sol_spec_file);
    void get_simple();
    virtual void get_A() override;
    virtual void get_b() override;
//    virtual void display_ass(double values[], int n, std::vector<int*> indexes);
    virtual void display_A(int n) override;
    virtual void display_x(int n) override;
    virtual void display_b(int n) override;
    virtual void display_r(int n) override;
    virtual void display_A() override;
    virtual void display_x() override;
    virtual void display_b() override;
    virtual void display_r() override;
    virtual void init() override;
    virtual void analyse() override;
    virtual void factorize() override;
    virtual void alloc_solve_residual() override;
    virtual void alloc_rhs() override;
    virtual void solve() override;
    virtual void metrics() override;
    virtual void call() override;
    virtual bool get_b_before_facto() override;
    virtual void set_no_output() override;
    virtual void finalize_solve() override;
    virtual void finalize() override;
    virtual void output_metrics_init(std::string file) override;
    virtual void output_metrics(std::string sol_spec_file) override;
};


#endif /* QR_MUMPS_H */

