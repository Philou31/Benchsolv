/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Benchmark_Impl.h
 * Author: filou
 *
 * Created on 15 novembre 2016, 10:33
 */

#ifndef BENCHMARK_IMPL_H
#define BENCHMARK_IMPL_H

#include "Benchmark.h"
#include "constants.h"


template <class S, typename K, typename V>
Benchmark<S,K,V>::Benchmark(S *solver, std::string multiple_bench,
        std::string b_file, std::string o_file, std::string a_file,
        std::string f_file, std::string s_file, std::string sol_spec_file,
        bool output_metrics): 
    _solver(solver), _multiple_bench(multiple_bench), _b_file(b_file), 
    _o_file(o_file), _a_file(a_file), _f_file(f_file), _s_file(s_file), 
    _sol_spec_file(sol_spec_file), _output_metrics(output_metrics)
{
    if (_solver->is_host()) {
        std::cout << "\nOutputs:\n";
        // Single benchmark
        if (b_file != cst::EMPTY_FILE && _multiple_bench == cst::SINGLE)
            std::cout << "benchmark option file: " << b_file << "\n";
        // Multiple benchmark
        if (_multiple_bench == cst::MULTIPLE) {
            if (o_file != cst::EMPTY_FILE)
                std::cout << "output options file: " << o_file << "\n";
            if (a_file != cst::EMPTY_FILE)
                std::cout << "analysis options file:" << a_file << "\n";
            if (f_file != cst::EMPTY_FILE)
                std::cout << "factorization options file: " << f_file << "\n";
            if (s_file != cst::EMPTY_FILE)
                std::cout << "solve options file: " << s_file << "\n";
        }
        // Solution specific metrics
        if (sol_spec_file != cst::EMPTY_FILE)
            std::cout << "solution specific output file: " << sol_spec_file << "\n\n";
    }
}

template <class S, typename K, typename V>
Benchmark<S,K,V>::~Benchmark()
{}

template <class S, typename K, typename V>
std::string Benchmark<S,K,V>::to_string(int i) {
    return std::to_string(static_cast<long long>(i));
}
template <class S, typename K, typename V>
std::string Benchmark<S,K,V>::to_string(std::string s){
    return s;
}
  
    
////////////////////////////////////////////////////
// OPTIONS HANDLING
////////////////////////////////////////////////////  
template <class S, typename K, typename V>
bool Benchmark<S,K,V>::parse_options(std::string file) {
    // If no file, no options to parse
    if (file == cst::EMPTY_FILE) return false;

    bool parsed = false;
    std::ifstream infile(file.c_str());
    std::cout << "\nGetting options in file: " << file << "\n\n";

    // If we couldn't open the output file stream for reading
    if (!infile) {
        // Print an error and exit
        std::cerr << "Uh oh, " << file << " could not be opened for reading!\n";
        exit(1);
    }

    std::string line;
    bool is_key = true;
    char delimiter = cst::OPTS_DELIMITER;
    K key;
    V value;
    // Read line per line the options file
    while (infile) {
        std::getline(infile, line);
        std::istringstream stream(line);
        std::string strToken;
        // In a line there if first the key then all values
        while(std::getline(stream, strToken, delimiter)) {
            std::stringstream stream(strToken);
            // First the key
            if (is_key) {
                stream >> key;
                // Add key and initialize the corresponding array in _opts_val
                _opts.push_back(key);
                _opts_vals.push_back({});
                parsed = true;
            } else {
                // Add value
                stream >> value;
                _opts_vals.back().push_back(value);
            }
            is_key = false;
        }
        is_key = true;
    }
    // Display options if any
    if (parsed) display_options();
    return parsed;
}
    
template <class S, typename K, typename V>
bool Benchmark<S,K,V>::iterate_options() {
    int tmp_current_size = 0;
    // If keys left
    if (_opts.size() > 0) {
        // If vals left
        if (_opts_vals.front().size() >= 1) { 
            // save number of vals for current key
            if (_current_vals_size == 0)
                _current_vals_size = _opts_vals.front().size();
            // save current key/value
            K key = _opts.front();
            V value = _opts_vals.front().front();
            _current_key = key;
            _current_value = value;
            // set current option key/value
            _solver->set_opt(key, value, _sol_spec_file);
            // remove value and, if last value, remove key also
            _opts_vals.front().pop_front();
            if (_opts_vals.front().size() == 0) {
                _opts_vals.pop_front();
                _opts.pop_front();
                tmp_current_size = _current_vals_size;
                _current_vals_size = 0;
                // If last key and last value, launch test
                if (tmp_current_size == 1 && _opts.size() == 0) return true;
                // If than 1 value for key, launch test
                if (tmp_current_size > 1) return true;
                // Else (last value of not last key), no test, recurrent call
                return iterate_options();
            }
            // If not last value, launch test
            return true;
        }
        // No value left, remove array ov values
        else _opts_vals.pop_front();
    }
    // No key, no test
    return false;
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::display_options() {
    std::deque<K> opts(_opts);
    std::deque<std::deque<V>> opts_vals(_opts_vals);
    if (_solver->is_host()) {
        int size = 0;
        std::clog << "\nOptions to test:\n";
        for(unsigned int it=0; it < _opts.size(); ++it) {
            std::clog << "Key: " << opts.front() << ", Values:";
            opts.pop_front();
            size = opts_vals.front().size();
            for(int itt=0; itt < size; ++itt) {
                std::clog << " " << opts_vals.front().front();
                opts_vals.front().pop_front();
            }
            std::clog << "\n";
            opts_vals.pop_front();
        }
        std::clog << "\n";
    }
}
    
template <class S, typename K, typename V>
bool Benchmark<S,K,V>::iterate_solver(std::string file, bool a, 
        bool f, bool s, bool o) {
    bool iterated = false;
    parse_options(file);
    // Loop while there are options to test
    while (iterate_options()) {
        // display option on host
        if (_solver->is_host())
            std::clog << "\nNEW TEST WITH KEY: " << _current_key << ", VALUE: " 
                << _current_value << "\n";
        // If modified, read b again:
        //      - Before factorisation (if icntl32=1: forward elim in facto)
        //      - Else, before solve
        // TODOMAYBE: For compatibility with other solvers, do the following
        // inside Mumps ?
        if (!_got_b) {
            bool before_facto=_solver->get_b_before_facto();
            if (f && before_facto) _solver->get_b_again();
            if (s && !before_facto) _solver->get_b_again();
        }
        // call the solver on phases in arguments
        call(a, f, s, o);
        // Result is true if at least a call
        iterated = true;
        ++_nb_tests;
    }
    return iterated;
}
    
    
////////////////////////////////////////////////////
// RUNNING THE SOLVER
////////////////////////////////////////////////////
template <class S, typename K, typename V>
void Benchmark<S,K,V>::phase(void (Solver::*function)(), const std::string name, 
        long long int &t) {
    // Analysis
    if (_solver->is_host())
	    std::cout << "\nStarting the " << name << "\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    (_solver->*function)();
    
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    if (_solver->is_host()) {
        std::cout << "Time to do the " << name << " : " << cst::TIME_RATIO*t << "\n";
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::output_metrics() {
    long long t;
    if (_solver->is_host())
        std::cout << "Starting the Metrics\n";
    auto start = std::chrono::high_resolution_clock::now();
    _solver->metrics();
    if (!_multiple_bench.compare(cst::OPTION) || (_b_file == _o_file && 
        _o_file == _a_file && _a_file == _f_file && _f_file == _s_file && 
        _s_file == cst::EMPTY_FILE)) {
        _solver->output_metrics(_sol_spec_file, _ta, _tf, _ts);
    }
    // If not single_bench, output the current option 
    else {
        _solver->output_metrics(_sol_spec_file, _ta, _tf, _ts, 
            to_string(_current_key), to_string(_current_value));
    }
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    if (_solver->is_host())
        std::cout << "Time to do the metrics : " << cst::TIME_RATIO*t << "\n";
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::call(bool a, bool f, bool s, bool o) {
    if (a) {
        phase(&Solver::analyse, cst::ANALYSIS_PHASE, _ta);
//         Is the mapping of A dependent on options ?????
//        _got_A=false;
    }
    if (f) {
        // In the case of distributed input, local parts must be read before
        // factorization (facto or mapping distribution). But it only has to be 
        // read once (at each analysis ?????)
        if (!_got_A) {
            _solver->get_A_again();
            _got_A=true;
        }
        phase(&Solver::factorize, cst::FACTORIZATION_PHASE, _tf);
    }
    long long int t;
    if (s) {
        phase(&Solver::solve, cst::SOLVE_PHASE, t);
        _got_b=false;
    }
    if (o && _output_metrics) output_metrics();
}
    
    
////////////////////////////////////////////////////
// RUN BENCHMARKS
////////////////////////////////////////////////////
template <class S, typename K, typename V>
void Benchmark<S,K,V>::run() {
    if (_output_metrics)
        _solver->output_metrics_init(_sol_spec_file);
    // If simple test (OPTION), call the whole solver once
    //      !compare is actually string equality
    if (!_multiple_bench.compare(cst::OPTION) || (_b_file == _o_file && 
        _o_file == _a_file && _a_file == _f_file && _f_file == _s_file && 
        _s_file == cst::EMPTY_FILE))
        call();
    // Single options file launching tests with all phases
    else if (!_multiple_bench.compare(cst::SINGLE)) single_benchmark();
    // Multiple options file and multiple tests depending on file
    else multiple_benchmark();
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::single_benchmark() {
    bool b = iterate_solver(_b_file);
    // If no benchmark was launched, launch a test, at least
    if(!b) call();
    std::cout << "\nSingle benchmark executed.\n";
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::multiple_benchmark() {
    // output options: launch complete test
    iterate_solver(_o_file);
    // analysis options: idem and if no previous test, launch an analysis
    bool a = iterate_solver(_a_file);
    if (!a) call(true, false, false, false);
    // facto options: start test at facto and if no previous test, launch a facto
    bool f = iterate_solver(_f_file) || a;
    if (!f) call(false, true, false, false);
    // solve options: start test at solve and if no previous test, launch a solve
    bool s = iterate_solver(_s_file) || f;
    if (!s) call(false, false);
    // Display number of tests
    if (_solver->is_host())
        std::cout << "\nNumber of benchmarks executed: " << _nb_tests << "\n";
}

#endif /* BENCHMARK_IMPL_H */

