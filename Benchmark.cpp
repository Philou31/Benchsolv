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
Benchmark<S,K,V>::Benchmark(S *solver, std::string b_file, std::string o_file, 
    std::string a_file, std::string f_file, std::string s_file, 
    std::string sol_spec_file): 
    _solver(solver), _b_file(b_file), _o_file(o_file), _a_file(a_file), _f_file(f_file), 
    _s_file(s_file), _sol_spec_file(sol_spec_file)
{
    std::cout << "\nOutputs:\n";
    if (b_file != cst::EMPTY_FILE)
        std::cout << "benchmark option file: " << b_file << "\n";
    if (o_file != cst::EMPTY_FILE)
        std::cout << "output options file: " << o_file << "\n";
    if (a_file != cst::EMPTY_FILE)
        std::cout << "analysis options file:" << a_file << "\n";
    if (f_file != cst::EMPTY_FILE)
        std::cout << "factorization options file: " << f_file << "\n";
    if (s_file != cst::EMPTY_FILE)
        std::cout << "solve options file: " << s_file << "\n";
    if (sol_spec_file != cst::EMPTY_FILE)
        std::cout << "solution specific output file: " << sol_spec_file << "\n\n";
}

template <class S, typename K, typename V>
Benchmark<S,K,V>::~Benchmark()
{}
    
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
    // Read banner and m, n, nz line
    while (infile) {
        std::getline(infile, line);
        std::istringstream stream(line);
        std::string strToken;
        while(std::getline(stream, strToken, delimiter)) {
            std::stringstream stream(strToken);
            if (is_key) {
                stream >> key;
                _opts.push_back(key);
                _opts_vals.push_back({});
                parsed = true;
            } else {
                stream >> value;
                _opts_vals.back().push_back(value);
            }
            is_key = false;
        }
        is_key = true;
    }
    if (parsed) display_options();
    return parsed;
}
    
template <class S, typename K, typename V>
bool Benchmark<S,K,V>::iterate_options() {
    int tmp_current_size = 0;
    if (_opts.size() > 0) {
        if (_opts_vals.front().size() >= 1) { 
            if (_current_vals_size == 0)
                _current_vals_size = _opts_vals.front().size();
            K key = _opts.front();
            V value = _opts_vals.front().front();
            _current_key = key;
            _current_value = value;
            _solver->set_opt(key, value, _sol_spec_file);
            _opts_vals.front().pop_front();
            if (_opts_vals.front().size() == 0) {
                _opts_vals.pop_front();
                _opts.pop_front();
                tmp_current_size = _current_vals_size;
                _current_vals_size = 0;
                if (tmp_current_size == 1 && _opts.size() == 0) return true;
                if (tmp_current_size > 1) return true;
                return iterate_options();
            }
            return true;
        } else _opts_vals.pop_front();
    }
    return false;
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::display_options() {
    std::deque<K> opts(_opts);
    std::deque<std::deque<V>> opts_vals(_opts_vals);
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
    
template <class S, typename K, typename V>
bool Benchmark<S,K,V>::iterate_solver(std::string file, bool a, 
        bool f, bool s, bool o) {
    bool iterated = false;
    parse_options(file);
    while (iterate_options()) {
        std::clog << "\nNEW TEST WITH KEY: " << _current_key << ", VALUE: " 
                << _current_value << "\n";
        call(a, f, s, o);
        iterated = true;
        ++_nb_tests;
    }
    return iterated;
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::analysis() {
    // Analysis
    if (_solver->is_host())
	    std::cout << "\nStarting the analysis\n";
    auto start = std::chrono::high_resolution_clock::now();
    _solver->analyse();
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    _ta = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    _ta_tot += _ta;
    if (_solver->is_host()) {
        std::cout << "Time to do the analysis            : " << _ta << "\n";
        std::cout << "Total time to do the analysis      : " << _ta_tot << "\n";
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::factorize() {
    if (_solver->is_host())
        std::cout << "Starting the factorization\n";
    auto start = std::chrono::high_resolution_clock::now();
     _solver->factorize();
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    _tf = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    _tf_tot += _tf;
    if (_solver->is_host()) {
        std::cout << "Time to do the facto               : " << _tf << "\n";
        std::cout << "Total time to do the facto         : " << _tf_tot << "\n";
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::solve() {
    if (_solver->is_host())
        printf("Starting the solve\n");
    auto start = std::chrono::high_resolution_clock::now();
    _solver->solve();
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    _ts = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    _ts_tot += _ts;
    if (_solver->is_host()) {
        std::cout << "Time to do the solve               : " << _ts << "\n";
        std::cout << "Total time to do the solve         : " << _ts_tot << "\n";
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::output_metrics() {
    _solver->metrics();
    _solver->output_metrics(_sol_spec_file, _ta, _tf, _ts);
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::call(bool a, bool f, bool s, bool o) {
//    std::cout << "A:\n";
//    _solver->display_A(10);
//    std::cout << "b:\n";
//    _solver->display_b(10);
    try {
        bool got_b = false;
        if (a) analysis();
        if (f) {
            if (!got_b && _solver->get_b_before_facto()) {
                _solver->get_b();
                got_b = true;
            }
            factorize();
        }
        if (s) {
            if (!got_b && !_solver->get_b_before_facto()) {
                _solver->get_b();
                got_b = true;
            }
            solve();
        }
        if (o) output_metrics();
    } catch(...) {
        std::cerr << "Error On This Test.\n";
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::benchmark(std::string multiple_bench) {
    _solver->output_metrics_init(_sol_spec_file);
    if (!multiple_bench.compare(cst::OPTION) || (_b_file == _o_file && 
        _o_file == _a_file && _a_file == _f_file && _f_file == _s_file && 
        _s_file == cst::EMPTY_FILE)) {
        call();
    }
    else if (!multiple_bench.compare(cst::SINGLE)) single_benchmark();
    else multiple_benchmark();
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::single_benchmark() {
    bool b = iterate_solver(_b_file);
    if(b) std::cout << "\nSingle benchmark executed.\n";
    else {
        call();
    }
}
    
template <class S, typename K, typename V>
void Benchmark<S,K,V>::multiple_benchmark() {
    iterate_solver(_o_file);
    bool a = iterate_solver(_a_file);
    // If no previous test, run an analysis
    if (!a) call(true, false, false, false);
    bool f = iterate_solver(_f_file);
    // If no previous test, run a facto
    if (!f) call(false, true, false, false);
    bool s = iterate_solver(_s_file);
    // If no previous test, run a solve
    if (!s) call(false, false);
    std::cout << "\nNumber of benchmarks executed: " << _nb_tests << "\n";
}

#endif /* BENCHMARK_IMPL_H */

