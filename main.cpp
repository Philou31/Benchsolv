//! \file main.cpp
//! \brief Mumps resolution of Ax=b
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Resolution of Ax=b using MUMPS. Details of the parameters:
//!     - blablabla
//!
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "cmdline.h"
#include "Mumps.h"
#include "QR_Mumps.h"
#include "Benchmark.h"
#include "constants.h"

int main(int argc, char **argv) {
//    long long t;
//    auto start = std::chrono::high_resolution_clock::now();

    cmdline::parser a;
    // add specified type of variable.
    // 1st argument is long name
    // 2nd argument is short name (no short name if '\0' specified)
    // 3rd argument is description
    // 4th argument is mandatory (optional. default is false)
    // 5th argument is default value  (optional. it used when mandatory is false)
    // 6th argument is extra constraint.
    // Here, port number must be 1 to 65535 by cmdline::range().
    //    a.add<int>("port", 'p', "port number", false, 80, cmdline::range(1, 65535));
    //  Boolean flags also can be defined. Call add method without a type parameter.
    //  a.add("gzip", '\0', "gzip when transfer");

    // cmdline::oneof() can restrict options.
    
    ////////////////////////////////////////////////////
    // DEBUGGING ARGUMENTS
    ////////////////////////////////////////////////////
//    //Matrix A
//    a.add<std::string>("Amatrix", 'A', "File containing the A matrix", false,
//            cst::A_MHD_10_FILE);
//    a.add<std::string>("A_n_present", 'a', "Dimensions are present or not in A", 
//            false, cst::TRUE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
//    //Matrix B
//    a.add<std::string>("RHS", 'B', "File containing the Right Hand Side matrix",
//            false, cst::RHS_MHD_10_FILE);
////    a.add<std::string>("RHS", 'B', "File containing the Right Hand Side matrix",
////            false, cst::RHS_MHD_10_FILE);
//    a.add<std::string>("b_n_present", 'b', "Dimensions are present or not in b",
//            false, cst::TRUE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
//    //Mumps
//    a.add<int>("A_distribution", '=', "Distribution of the matrix A",
//            false, parm::A_CENTRALIZED, cmdline::oneof<int>(parm::A_DISTR_ANALYSIS,
//            parm::A_DISTR_FACTO, parm::A_DISTR_FACTO_MAPPING));
//    a.add<std::string>("A_loc", '&', "Are getting the matrix A locally or on master",
//            false, cst::FALSE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
//    a.add<int>("A_format", '_', "Format of the matrix A",
//            false, parm::A_ASSEMBLED_FORMAT, cmdline::oneof<int>(parm::A_ASSEMBLED_FORMAT,
//            parm::A_ELEMENTAL_FORMAT));
//    a.add<int>("A_symmetry", 'y', "Symmetry of the matrix A",
//            false, parm::SYM_UNSYM, cmdline::oneof<int>(parm::SYM_UNSYM,
//            parm::SYM_GENERAL, parm::SYM_DEFPOS));
//    a.add<int>("working_host", 'w', "The host is working or not",
//            false, parm::WORKING_HOST, cmdline::oneof<int>(parm::WORKING_HOST, 
//            parm::NOT_WORKING_HOST));
//    //Benchmark
//    ////Test id
//    a.add<std::string>("test_id", '$', "Id of the test, suffix of the output files",
//            false, "42");
//    ////Solver
//    a.add<std::string>(cst::SOLVER, cst::SOLVER_C, "Type of solver  for the test", false,
//            cst::MUMPS, cmdline::oneof<std::string>(cst::MUMPS, cst::QR_MUMPS, cst::MUMPS_C, 
//            cst::QR_MUMPS_C));
//    ////Test type
//    a.add<std::string>(cst::MULTIPLE_BENCH, cst::MULTIPLE_BENCH_C, "Run multiple benchmarks or just options from analysis file",
//            false, cst::OPTION, 
//            cmdline::oneof<std::string>(cst::MULTIPLE, cst::SINGLE, cst::OPTION));
//    ////Option test
//    a.add<std::string>("string_opt_key", '(', "String key of the option to change (qr_mumps)",
//            false, cst::EMPTY_STRING_OPT_KEY);
//    a.add<int>("int_opt_key", ')', "Integer key of the option to change (mumps)",
//            false, cst::EMPTY_INT_OPT_KEY);
//    a.add<int>("int_opt_value", '-', "Integer value of the option to change (mumps/qr_mumps)",
//            false, cst::EMPTY_INT_OPT_VALUE);
//    ////Single test
//    a.add<std::string>("bench_opt", '!', "File containing the options to test in single benchmark",
//            false, "options/run.params");
//    ////Multiple tests (output, analysis, factorization, solve)
//    a.add<std::string>("output", 't', "File containing the options to test for output",
//            false, "options/mumps/output.opt");
//    a.add<std::string>("analysis", 'z', "File containing the options to test in analysis",
//            false, "options/mumps/analysis.opt");
//    a.add<std::string>("facto", 'i', "File containing the options to test in factorisation",
//            false, "options/mumps/factorization.opt");
//    a.add<std::string>("solve", 'p', "File containing the options to test in solve",
//            false, "options/mumps/solve.opt");
//    //Output
//    ////Output files
//    a.add<std::string>("fortran_output", 'm', "File where all fortran outputs will go",
//            false, "res/fortran");
//    a.add<std::string>("output_file", 'o', "File where all normal outputs will go",
//            false, "res/out");
//    a.add<std::string>("error_file", 'e', "File where all error outputs will go",
//            false, "res/err");
//    ////Metrics files
//    a.add<std::string>("sol_spec_metrics", 'l', "File where all the solution specific metrics will go",
//            false, "res/sol_spec.txt");
//    a.add<std::string>("pb_spec_metrics", 'k', "File where all the problem specific metrics will go",
//            false, "res/pb_spec.txt");
    
    
    ////////////////////////////////////////////////////
    // PRODUCTION ARGUMENTS
    ////////////////////////////////////////////////////
    //Matrix A
    a.add<std::string>("Amatrix", 'A', "File containing the A matrix", 
            true, cst::A_MHD_10_FILE);
    a.add<std::string>("A_n_present", 'a', "Dimensions are present or not in A", 
            false, cst::TRUE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
    //Matrix B
    a.add<std::string>("RHS", 'B', "File containing the Right Hand Side matrix",
            true, cst::EMPTY_FILE);
    a.add<std::string>("b_n_present", 'b', "Dimensions are present or not in b",
            false, cst::TRUE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
    //Mumps
    a.add<int>("A_distribution", '=', "Distribution of the matrix A",
            false, parm::A_CENTRALIZED, cmdline::oneof<int>(parm::A_DISTR_ANALYSIS,
            parm::A_DISTR_FACTO, parm::A_DISTR_FACTO_MAPPING));
    a.add<std::string>("A_loc", '&', "Are getting the matrix A locally or on master",
            false, cst::FALSE, cmdline::oneof<std::string>(cst::TRUE, cst::FALSE));
    a.add<int>("A_format", '_', "Format of the matrix A",
            false, parm::A_ASSEMBLED_FORMAT, cmdline::oneof<int>(parm::A_ASSEMBLED_FORMAT,
            parm::A_ELEMENTAL_FORMAT));
    a.add<int>("A_symmetry", 'y', "Symmetry of the matrix A",
            false, parm::SYM_UNSYM, cmdline::oneof<int>(parm::SYM_UNSYM,
            parm::SYM_GENERAL, parm::SYM_DEFPOS));
    a.add<int>("working_host", 'w', "The host is working or not",
            false, parm::WORKING_HOST, cmdline::oneof<int>(parm::WORKING_HOST, 
            parm::NOT_WORKING_HOST));
    //Benchmark
    ////Test id
    a.add<std::string>("test_id", '$', "Id of the test, suffix of the output files",
            true, "");
    ////Solver
    a.add<std::string>(cst::SOLVER, cst::SOLVER_C, "Type of solver  for the test", true,
            cst::MUMPS, cmdline::oneof<std::string>(cst::MUMPS, cst::QR_MUMPS, cst::MUMPS_C, 
            cst::QR_MUMPS_C));
    ////Test type
    a.add<std::string>(cst::MULTIPLE_BENCH, cst::MULTIPLE_BENCH_C, "Run multiple benchmarks or just options from analysis file",
            false, cst::OPTION, 
            cmdline::oneof<std::string>(cst::MULTIPLE, cst::SINGLE, cst::OPTION));
    ////Option test
    a.add<std::string>("string_opt_key", '(', "String key of the option to change (qr_mumps)",
            false, cst::EMPTY_STRING_OPT_KEY);
    a.add<int>("int_opt_key", ')', "Integer key of the option to change (mumps)",
            false, cst::EMPTY_INT_OPT_KEY);
    a.add<int>("int_opt_value", '-', "Integer value of the option to change (mumps/qr_mumps)",
            false, cst::EMPTY_INT_OPT_VALUE);
    ////Single test
    a.add<std::string>("bench_opt", '!', "File containing the options to test in single benchmark",
            false, "options/run.params");
    ////Multiple tests (output, analysis, factorization, solve)
    a.add<std::string>("output", 't', "File containing the options to test for output",
            false, cst::EMPTY_FILE);
    a.add<std::string>("analysis", 'z', "File containing the options to test in analysis",
            false, cst::EMPTY_FILE);
    a.add<std::string>("facto", 'i', "File containing the options to test in factorisation",
            false, cst::EMPTY_FILE);
    a.add<std::string>("solve", 'p', "File containing the options to test in solve",
            false, cst::EMPTY_FILE);
    //Output
    ////Output files
    a.add<std::string>("fortran_output", 'm', "File where all fortran outputs will go",
            false, "res/fortran");
    a.add<std::string>("output_file", 'o', "File where all normal outputs will go",
            false, "res/out");
    a.add<std::string>("error_file", 'e', "File where all error outputs will go",
            false, "res/err");
    ////Metrics files
    a.add<std::string>("sol_spec_metrics", 'l', "File where all the solution specific metrics will go",
            false, "res/sol_spec.txt");
    a.add<std::string>("pb_spec_metrics", 'k', "File where all the problem specific metrics will go",
            false, "res/pb_spec.txt");

    // Run parser.
    // It returns only if command line arguments are valid.
    // If arguments are invalid, a parser output error msgs then exit program.
    // If help flag ('--help' or '-?') is specified, a parser output usage 
    // message then exit program.
    a.parse_check(argc, argv);
    
    // boolean flags are referred by calling exist() method.
    //Matrix A
    std::string A_file = a.get<std::string>("Amatrix");
    bool An = !cst::TRUE.compare(a.get<std::string>("A_n_present"));
    //Matrix B
    std::string b_file = a.get<std::string>("RHS");
    bool bn = !cst::TRUE.compare(a.get<std::string>("b_n_present"));
    //Benchmark
    ////Solver
    std::string solver = a.get<std::string>("solver");
    ////Mumps
    int sym = a.get<int>("A_symmetry");
    int distr = a.get<int>("A_distribution");
    bool loc = !cst::TRUE.compare(a.get<std::string>("A_loc"));
    int format = a.get<int>("A_format");
    int par = a.get<int>("working_host");
    ////Test id
    std::string test_id = a.get<std::string>("test_id");
    ////Test type
    std::string multiple_bench = a.get<std::string>(cst::MULTIPLE_BENCH);
    ////Option test
    std::string string_opt_key = a.get<std::string>("string_opt_key");
    int int_opt_key = a.get<int>("int_opt_key");
    int int_opt_value = a.get<int>("int_opt_value");
    ////Single test
    std::string bench_file = a.get<std::string>("bench_opt");
    ////Mulitple tests (output, analysis, factorization, solve)
    std::string out_file = a.get<std::string>("output");
    std::string anal_file = a.get<std::string>("analysis");
    std::string facto_file = a.get<std::string>("facto");
    std::string sol_file = a.get<std::string>("solve");
    //Output
    ////Output files
    std::string suffix = "";
    if (test_id != "")
        suffix = "_" + test_id;
    std::string fortran_output = a.get<std::string>("fortran_output") + suffix;
    std::string output_file = a.get<std::string>("output_file") + suffix;
    std::string error_file = a.get<std::string>("error_file") + suffix;
    ////Metrics files
    std::string sol_spec_file = a.get<std::string>("sol_spec_metrics");
    std::string pb_spec_file = a.get<std::string>("pb_spec_metrics");
    
    FILE *f = freopen(fortran_output.c_str(), "a", stdout);
    std::cout << "Redirecting fortran output to : " << fortran_output << 
            " at adress: " << f << "\n";
    std::ofstream coutstr(output_file, std::ofstream::app);
    std::cout.rdbuf(coutstr.rdbuf());
    std::ofstream cerrstr(error_file, std::ofstream::app);
    std::cerr.rdbuf(cerrstr.rdbuf());
    
//    auto elapsed = std::chrono::high_resolution_clock::now() - start;
//    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
//    std::clog << "Time to do the init               : " << cst::TIME_RATIO*t << "\n";

    if (solver == "mumps" || solver == "m") {
//    start = std::chrono::high_resolution_clock::now();
        Mumps s(test_id, A_file, An, b_file, bn, par, sym, distr, loc, format,
            cst::USE_COMM_WORLD, MPI_COMM_WORLD, pb_spec_file, int_opt_key, 
            int_opt_value);
//    elapsed = std::chrono::high_resolution_clock::now() - start;
//    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
//    std::clog << "Time to init MUMPS                : " << cst::TIME_RATIO*t << "\n";
//
//    start = std::chrono::high_resolution_clock::now();
        Benchmark<Mumps, int, int> b(&s, bench_file, out_file, anal_file, 
            facto_file, sol_file, sol_spec_file);
//    elapsed = std::chrono::high_resolution_clock::now() - start;
//    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
//    std::clog << "Time to init Benchmark            : " << cst::TIME_RATIO*t << "\n";
//
//    start = std::chrono::high_resolution_clock::now();
        b.benchmark(multiple_bench);
//    elapsed = std::chrono::high_resolution_clock::now() - start;
//    t = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
//    std::clog << "Time for the Bench                : " << cst::TIME_RATIO*t << "\n";
    } else if (solver == "qr_mumps" || solver == "qrm") {
        QR_Mumps s(test_id, A_file, An, b_file, bn, string_opt_key, 
            int_opt_value);
        Benchmark<QR_Mumps, std::string, int> b(&s, bench_file, out_file, 
            anal_file, facto_file, sol_file, sol_spec_file);
        b.benchmark(multiple_bench);
    }
}
