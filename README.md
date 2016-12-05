# Benchsolv
**Framework: Automatic tool to launch benchmarks on a solver (Mumps, QR_mumps, ABCD, HPDDM)**

## Requirements:
* Working solvers+dependencies
* C++11

## Running tests
usage: /home/filou/NetBeansProjects/Benchmark_test/dist/Release/GNU-Linux/benchmark_test --Amatrix=string --test_id=string --solver=string [options] ... 
options:
  -A, --Amatrix             File containing the A matrix (string)
  -a, --A_n_present         Dimensions are present or not in A (string [=true])
  -B, --RHS                 File containing the Right Hand Side matrix (string [=])
  -b, --b_n_present         Dimensions are present or not in b (string [=true])
  -=, --A_distribution      Distribution of the matrix A (int [=0])
  -&, --A_loc               Are getting the matrix A locally or on master (string [=false])
  -_, --A_format            Format of the matrix A (int [=0])
  -y, --A_symmetry          Symmetry of the matrix A (int [=0])
  -w, --working_host        The host is working or not (int [=1])
  -$, --test_id             Id of the test, suffix of the output files (string)
  -s, --solver              Type of solver  for the test (string)
  -q, --multiple_bench      Run multiple benchmarks or just options from analysis file (string [=option])
  -(, --string_opt_key      String key of the option to change (qr_mumps) (string [=])
  -), --int_opt_key         Integer key of the option to change (mumps) (int [=-42])
  --, --int_opt_value       Integer value of the option to change (mumps/qr_mumps) (int [=-42])
  -!, --bench_opt           File containing the options to test in single benchmark (string [=options/run.params])
  -t, --output              File containing the options to test for output (string [=])
  -z, --analysis            File containing the options to test in analysis (string [=])
  -i, --facto               File containing the options to test in factorisation (string [=])
  -p, --solve               File containing the options to test in solve (string [=])
  -m, --fortran_output      File where all fortran outputs will go (string [=res/fortran])
  -o, --output_file         File where all normal outputs will go (string [=res/out])
  -e, --error_file          File where all error outputs will go (string [=res/err])
  -l, --sol_spec_metrics    File where all the solution specific metrics will go (string [=res/sol_spec.txt])
  -k, --pb_spec_metrics     File where all the problem specific metrics will go (string [=res/pb_spec.txt])
  -?, --help                print this message

## Test:
	- matrix A
		- A_n_present correspond to the presence or not of the matrices size as header of the matrix files [true]
		- distribution of A (if distributed, the name of the matrix file is only a prefix to which we add "_$proc_id", we specify then that A is loc.) [centralized]
		- format of A [assembled]
		- Symetry of A [unsymmetric]
	- matrix b (if not given, it will be allocated as  all 1.0s)
		- b_n_present: idem as A_n_present
	- An id specific to the test for the output metrics

## Different kinds of tests possible (with default options)
	1) Option:
		- arguments string_opt_key/int_opt_key/int_opt_value for a specific option to test
		- if none given, just a test with default options
	2) Single bench:
		- argument bench_opt: tab-delimited file with the 2 columns: key\tvalue corresponding to all values to be tested simultaneously. If an option is given a value multiple times, only the last one is taken.
	3) Multiple bench: + FIGURE TESTS LAUNCHED
		- arguments output/analysis/facto/solve: tab-delimited files with columns key\tvalue1\tvalue2...
		- If an option is given multiple values on 1 line, each value corresponds to a test launched. If only one value is given on a line, a test is launched only on the next line. Exception for last line, a test is always launched. The last value of the line will be kept for all following tests.
		- The different files all have different actions:
			- Tests in output and analysis files launch complete test (analysis+factorization+solve phases). There are two different files because output is intended only for test of output options.
			- Tests in factorization launch only factorization+solve phases
			- Tests in solve launch only the solve phase.

## output:
	1) fortran_output contains the output for fortran calculus (currently only MUMPS)
	2) output_file: content of cout stream
	3) error_file: content of cout stream
	4) sol_spec_metrics: output metrics for each test in a tab-delimited file
		- !! For each test a line with metrics names is printed
		- !! When the value of an option is changed, a line with "Setting x to y" is printed
	5) pb_spec_metrics: output metrics specific to the problem Ax=b and not the solve

## Modular organisation:
	- Benchmark launch all the tests
	- Solver is an interface identical for all solvers => Any solver can be added to the framework quite easily by just implementing the different functions. Except the functions to get the matrices from the files which are common to all !
	- You can even use the class Metrics to get the same system for most metrics for solvers. It is based on QR_Mumps => A, b, solution, residual and orthogonal residual norms. (Currently not working)
	- A special mention about the Mumps solver: It is the only one outputting the Ax=b specific metrics for the problem
  
