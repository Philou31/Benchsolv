/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Mumps.h
 * Author: filou
 *
 * Created on 11 novembre 2016, 14:42
 */

#ifndef MUMPS_H
#define MUMPS_H

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define INFO(I) info[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFO(I) rinfo[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#include <iostream>
#include <sstream>
#include "dmumps_c.h"
#include "mpi.h"
#include "Parameters_Mumps.h"
#include "Solver.h"
#include "Metrics.h"

class Mumps : public Solver {
private:
    DMUMPS_STRUC_C _id;
    MPI_Comm _mpi_comm;
    int _nb_procs = 1, _proc_id = 0;
    double *_r;
    Metrics _metrics;
    // Problem specific info
    std::string _pb_spec_file;
    int _distr;
    bool _loc;
    int _format;
    int _opt_key, _opt_value;
public:
    Mumps(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, int par, int sym, int distr, 
        bool loc, int format, int comm, MPI_Comm mpi_comm, 
        std::string pb_spec_file, int int_opt_key, int int_opt_value);
    ~Mumps();
    virtual bool is_host() override;
    void set_opt(int key, int value);
    void set_opt(int key, int value, std::string sol_spec_file);
    void get_simple();
    virtual void get_A() override;
    virtual void get_b() override;
    virtual void display_A(int n) override;
    void display_A_loc(int n);
    virtual void display_x(int n) override;
    virtual void display_b(int n) override;
    virtual void display_r(int n) override;
    virtual void display_A() override;
    void display_A_loc();
    virtual void display_x() override;
    virtual void display_b() override;
    virtual void display_r() override;
    void mumps(int job);
    virtual void init() override;
    virtual long long total_time(long long *t) override;
    virtual void analyse() override;
    virtual void factorize() override;
    virtual void alloc_solve_residual() override;
    virtual void alloc_rhs() override;
    virtual void solve() override;
    void assemble_A();
    virtual void metrics() override;
    virtual void call() override;
    virtual bool get_b_before_facto() override;
    virtual void set_no_output() override;
    virtual void finalize_solve() override;
    void problem_spec_metrics_output();
    virtual void finalize() override;
    virtual void output_metrics_init(std::string file) override;
    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
        long long tf = 0, long long ts = 0, long long ta_tot = 0, 
        long long tf_tot = 0, long long ts_tot = 0, std::string key = "",
        std::string value = "") override;
};

#endif /* MUMPS_H */

