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
    std::string _opt_key;
    int _opt_value;
public:
    QR_Mumps(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, std::string string_opt_key, 
        int int_opt_value);
    ~QR_Mumps();
    virtual bool is_host() override;
    void set_opt(std::string key, int value);
    void set_opt(std::string key, int value, std::string sol_spec_file);
    virtual bool take_A_value_loc(int m, int n, int i, bool local) override;
    virtual int nz_loc(int nz, bool local) override;
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
    virtual long long total_time(long long *t) override;
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
    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
        long long tf = 0, long long ts = 0, long long ta_tot = 0, 
        long long tf_tot = 0, long long ts_tot = 0, std::string key = "",
        std::string value = "") override;
};


#endif /* QR_MUMPS_H */

