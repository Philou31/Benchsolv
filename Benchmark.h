//! \file Benchmark.h
//! \brief Benchmark class
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Benchmark is the class manipulating the Solvers. We can launch the test via
//! the run method of this class.
//!

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <chrono>
#include <queue>
#include "Solver.h"

template <class S, typename K, typename V>
class Benchmark {
private:
    S *_solver; // Current solver to test (Mumps, QR_Mumps,...)
    // type of bench to launch (option, multiple,...)
    std::string _multiple_bench;
    // Input options to test:
    //      Single benchmark (all options in a row)
    std::string _b_file;
    //      Output, analysis, factorization, solve options (per type separately)
    std::string _o_file, _a_file, _f_file, _s_file;
    std::string _sol_spec_file; // File for the metrics specific to a solution
    std::deque<K> _opts;    // Keys of the options to test
    std::deque<std::deque<V>> _opts_vals;   // Values of the options to test
    int _current_vals_size = 0; // #values to test for the current option
    K _current_key; // Current option key
    V _current_value; // Current option value
    int _nb_tests = 0;  // #Tests executed
    // Execution time for analysis, facto and solve
    long long _ta = 0, _tf = 0, _ts = 0;
    // Right hand side read and not modified. Initialized true: b always read at
    // solver initialization
    bool _got_b=true;
    bool _output_metrics=true;
public:
    //!
    //! \brief Constructor of the Benchmark class
    //!
    //! \param solver: The solver of the benchmark (Mumps, QR_Mumps,...)
    //! \param multiple_bench: type of bench to launch (option, multiple,...)
    //! \param b_file: Single benchmark options file
    //! \param o_file, a_file, f_file, s_file: output, analysis, facto, solve 
    //!     options
    //! \return the Metrics instance
    //!
    Benchmark(S *solver, std::string multiple_bench,
        std::string b_file, std::string o_file, std::string a_file,
        std::string f_file, std::string s_file, std::string sol_spec_file,
        bool output_metrics);
    
    //!
    //! \brief Destructor of the Benchmark class
    //!
    ~Benchmark();
    
    //!
    //! \fn std::string to_string(int i)
    //! \brief Cast integer to string
    //!
    //! \param i: integer
    //! \return the integer as string
    //!
    std::string to_string(int i);
    
    //!
    //! \fn std::string to_string(std::string s)
    //! \brief Cast integer to string
    //!
    //! \param s: string
    //! \return the integer as string
    //!
    std::string to_string(std::string s);
    
    
    ////////////////////////////////////////////////////
    // OPTIONS HANDLING
    ////////////////////////////////////////////////////
    //!
    //! \fn bool parse_options(std::string file)
    //! \brief Parse the options in the file and return true if there is one
    //!
    //! Parse the options in the file, these options should be of the form:
    //!     "key value1 value2 ... valueN"
    //! The parser will read each line and save it in the _opts, _opts_val 
    //! arrays.
    //!
    //! \param file: file containing the options
    //! \return true if there was an option to read in the file
    //!
    bool parse_options(std::string file);
    
    //!
    //! \fn iterate_options()
    //! \brief Iterate over options: set option, remove it and return true.
    //!
    //! Iterate over options in the _opts, _opts_val arrays. For each option:
    //!     - set option in the solver,
    //!     - remove option from the list
    //!     - return true
    //! return false if no option is left
    //!
    //! Remarks:
    //!     - For the last value of the current key, a test won't be launched: 
    //! another call to the method is done.
    //!     - The last value of the last key option is always launched, though.
    //!
    //! \return true if there was an option to iterate on
    //!
    bool iterate_options();
    
    //!
    //! \fn display_options()
    //! \brief Display the options to test
    //!
    //! Display all the options to test with a key per line:
    //!     "Options to test:
    //!     Key: key1, Values: value1 value2 ... valueN
    //!     ..."
    //!
    void display_options();
    
    //!
    //! \fn iterate_solver(std::string file, bool a=true, bool f=true, 
    //! bool s=true, bool o=true)
    //! \brief Launch all the tests on the solver with an options file
    //!
    //! This method will launch the tests for all the options in the file
    //! calling parse_options, iterate_options and call. It returns true if at
    //! least a test was launched.
    //!
    //! \param file: file containing the options to call parse_options
    //! \param a: analysis must be launched
    //! \param f: factorization must be launched
    //! \param s: solve must be launched
    //! \param o: output must be launched
    //! \return true if at least a test was launched
    //!
    bool iterate_solver(std::string file, bool a=true, bool f=true, 
        bool s=true, bool o=true);
    
    
    ////////////////////////////////////////////////////
    // RUNNING THE SOLVER
    ////////////////////////////////////////////////////
    //!
    //! \fn void phase(void (*function)(), std::string name, 
    //!     int &t)
    //! \brief Launch one of the phase while saving/outputting execution time
    //!
    //! From the function to call, the name of the phase, and the corresponding
    //! execution time, you can launch with this function a phase of the solver:
    //! analysis, factorization or solve.
    //!
    //! \param function: function call to the phase of the solver
    //! \param name: name of the phase for the outputting
    //! \param t: the corresponding execution time
    //!
    void phase(void (Solver::*function)(), const std::string name, 
        long long int &t);
    
    //!
    //! \fn void output_metrics()
    //! \brief Call the output_metrics method of the solver and output exec time
    //!
    void output_metrics();
    
    //!
    //! \fn void get_b_again()
    //! \brief get the right hand side again
    //!
    //! This function read the right hand side again. It should only be called
    //! if an option in the solver makes reading the right hand side again 
    //! mandatory at one stage.
    //! For example, in Mumps, the right hand side is modified during solve and
    //! to launch a new test with the same data, it is necessary to read b 
    //! again (before factorization or solve depending on the option icntl32).
    //!
    void get_b_again();
    
    //!
    //! \fn void call(bool a=true, bool f=true, bool s=true, bool o=true)
    //! \brief Launch some/all phases of the solver
    //!
    //! This function launch some/all phases of the solver depending on the
    //! parameters per phase
    //!
    //! \param a: launch analysis ?
    //! \param f: launch factorization ?
    //! \param s: launch solving ?
    //! \param o: launch output_metrics ?
    //!
    void call(bool a=true, bool f=true, bool s=true, bool o=true);
    
    
    ////////////////////////////////////////////////////
    // RUN BENCHMARKS
    ////////////////////////////////////////////////////
    //!
    //! \fn run()
    //! \brief Run the benchmark test(s)
    //!
    //! This function is the main one from benchmark. It allows to launch:
    //!     - a single test with the method call
    //!     - or a single benchmark with the method single_benchmark
    //!     - or a multiple benchmark with the method multiple_benchmark
    //!
    void run();
    
    //!
    //! \fn single_benchmark()
    //! \brief Run a single benchmark with all options from a file at once
    //!
    //! This function launch a single benchmark, setting all options from the
    //! attribute _b_file at once with iterate_solver and launching a test via 
    //! the method call. So only one test is launched, the options file must
    //! contain only lines with couples "Key Value"
    //!
    void single_benchmark();
    
    //!
    //! \fn multiple_benchmark()
    //! \brief Run multiple benchmarks from multiple options file
    //!
    //! This function allows to launch multiple tests using different options
    //! (one test per option) thanks to the iterate_options method. There should
    //! be a file per phase (analysis, factorization, solve): this allows to
    //! launch multiple test with only part of the solver launched. The analysis
    //! options file launch the whole solver, the facto file launch from facto
    //! only (no more analysis), ...
    //!
    void multiple_benchmark();
};

#include "Benchmark.cpp"

#endif /* BENCHMARK_H */

