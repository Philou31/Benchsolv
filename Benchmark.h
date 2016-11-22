/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Benchmark.h
 * Author: filou
 *
 * Created on 11 novembre 2016, 14:42
 */

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <fstream>
#include <iostream>
#include <chrono>
#include <queue>
#include "Solver.h"
#include "Mumps.h"
#include "QR_Mumps.h"

template <class S, typename K, typename V>
class Benchmark {
private:
    S *_solver;
    std::string _b_file, _o_file, _a_file, _f_file, _s_file, _sol_spec_file;
    std::deque<K> _opts;
    std::deque<std::deque<V>> _opts_vals;
    int _current_vals_size = 0;
    K _current_key;
    V _current_value;
    int _nb_tests = 0;
    long long _ta = 0, _tf = 0, _ts = 0;
    long long _ta_tot = 0, _tf_tot = 0, _ts_tot = 0;
public:
    Benchmark() {}
    Benchmark(S *solver, std::string b_file, std::string o_file, 
        std::string a_file, std::string f_file, std::string s_file, 
        std::string sol_spec_file);
    ~Benchmark();
    bool parse_options(std::string file);
    bool iterate_options();
    void display_options();
    bool iterate_solver(std::string file, bool a=true, bool f=true, 
        bool s=true, bool o=true);
    void init();
    void get_Ab();
    void analysis();
    void factorize();
    void solve();
    void output_metrics();
    void call(bool a=true, bool f=true, bool s=true, bool o=true);
    void benchmark(std::string multiple_bench);
    void single_benchmark();
    void multiple_benchmark();
};

#include "Benchmark.cpp"

#endif /* BENCHMARK_H */

