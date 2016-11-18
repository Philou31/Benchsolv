/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Solver.h
 * Author: filou
 *
 * Created on 11 novembre 2016, 14:41
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "Solver.h"
//#include "Metrics.h"

class Solver {
protected:
    std::string _file_A;
    bool _n_present_A;
    std::string _file_b;
    bool _n_present_b;
    double _rnrm{0}, _onrm{0}, _anrm{0}, _bnrm{0}, _xnrm{0};
public:
    Solver() {}
    Solver(std::string file_A, bool n_present_A, std::string file_b, 
        bool n_present_b);
    virtual bool is_host() = 0;
    void get_simple(int &m, int &n, int &nz, double **values, int **irn,
        int **jcn, int &nrhs, int &lrhs, double **rhs);
    void get_MM(std::string file, int &m, int &n, int &nz, 
        double **values, std::vector<int**> indexes, bool n_present, bool rhs);
    virtual void get_A() = 0;
    virtual void get_b() = 0;
    void display_ass(double values[], int n, std::vector<int*> indexes);
    virtual void display_A(int n) = 0;
    virtual void display_x(int n) = 0;
    virtual void display_b(int n) = 0;
    virtual void display_r(int n) = 0;
    virtual void display_A() = 0;
    virtual void display_x() = 0;
    virtual void display_b() = 0;
    virtual void display_r() = 0;
    virtual void init() = 0;
    virtual void analyse() = 0;
    virtual void factorize() = 0;
    virtual void alloc_solve_residual() = 0;
    virtual void alloc_rhs() = 0;
    virtual void solve() = 0;
    virtual void metrics() = 0;
    virtual void call() = 0;
    virtual bool get_b_before_facto() = 0;
    virtual void set_no_output() = 0;
    virtual void finalize_solve() = 0;
    virtual void finalize() = 0;
    virtual void output_metrics_init(std::string file) = 0;
    virtual void output_metrics(std::string sol_spec_file) = 0;
};

#endif /* SOLVER_H */

