//! \file Solver.h
//! \brief Abstract class for the solvers classes
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! This is an abstract class (some methods implemented) for the solver classes 
//! implemented by Mumps, QR_Mumps,... It implements only common functions:
//!     - Input matrix: read_MM
//!     - Output matrix: display...
//!

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <vector>
#include <omp.h>
#include <limits>
#include <cmath>
#include "constants.h"
#include "mpi.h"
#include "Solver.h"

class Solver {
protected:
    std::string _test_id;   // Id of the current test
    std::string _file_A;    // File of the matrix A
    bool _n_present_A;  // Presence of the size (m,n,nz) in file header
    std::string _file_b;    // File of the matrix b (can be empty)
    bool _n_present_b;  // Presence of the size (n, nrhs) in file header
    // Norm of the residual, orthogonal residual, A, b and the solution
    double _rnrm{0}, _onrm{0}, _anrm{0}, _bnrm{0}, _xnrm{0};
    int _nrows, _ncols, _nz = 0;    // number of rows, columns and non-zero in A
    int _nb_procs = 1, _proc_id = 0; // #MPI and id MPI
    int _nthreads = 1, _tid = 0;  // #OpenMP and thread id
    bool _A_assembled=false;    // true if A is assembled on host
    double *_r; // residual array
public:
    //!
    //! \brief Constructor of the Solver abstract class
    //!
    //! \param test_id: Id of the current test
    //! \param file_A: File of the matrix A
    //! \param n_present_A: Presence of the size (m,n,nz) in file header
    //! \param file_b: File of the matrix b (can be empty)
    //! \param n_present_b: Presence of the size (n, nrhs) in file header
    //! \param nrows: number of rows in A
    //! \param ncols: number of columns in A
    //! \param nz: number of non-zeros in A
    //! \return the Solver instance
    //!
    Solver(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, int nrows, int ncols, int nz);
    
    //!
    //! \fn void base_construct()
    //! \brief Initialize solver and get matrices
    //!
    //! \param init_before: true if solver init to be done before matrices read
    //!
    virtual void base_construct(bool init_before=true);

    //!
    //! \fn void base_destruct()
    //! \brief Destruct solver and free matrices
    //!
    virtual void base_destruct();
    
    ////////////////////////////////////////////////////
    // MPI communication methods
    ////////////////////////////////////////////////////
    //!
    //! \fn init_MPI()
    //! \brief Initialize and display MPI number of proc and id
    //!
    //! \param comm: the MPI communicator
    //!
    virtual void init_MPI(MPI_Comm &comm);
    
    //!
    //! \fn finalize_MPI()
    //! \brief Finalize MPI
    //!
    virtual void finalize_MPI();
    
    //!
    //! \fn init_OpenMP()
    //! \brief Initialize and display OpenMP nthreads and id
    //!
    virtual void init_OpenMP();
//    
//    //!
//    //! \fn assemble_MM()
//    //! \brief Assemble a Matrix Market on the host from the local parts
//    //!
//    //! This function will assemble on the host the local parts of a Matrix
//    //! Market structure to obtain the global matrix arrays (rows, columns, 
//    //! values). Caution: this can be heavy on memory.
//    //!
//    //! \param nz: number of non-zero values
//    //! \param values: pointer to the array of non-zero values
//    //! \param indexes: vector of pointers to the arrays of row and column index
//    //! \param nz_loc: number of local non-zero values
//    //! \param values_loc: pointer to the local array of non-zero values
//    //! \param indexes_loc: vector of pointers to the local arrays of row and column index
//    //! \param comm: the MPI communicator
//    //!
//    void assemble_MM(int nz, double **values, 
//        std::vector<int**> indexes, int nz_loc, double **values_loc, 
//        std::vector<int**> indexes_loc, MPI_Comm comm);
//    
    //!
    //! \fn assemble_A()
    //! \brief Assemble the matrix A on the host from the local parts
    //!
    //! This function will assemble on the host the local parts of the matrix A
    //! to obtain the global matrix arrays (rows, columns, values). Caution:
    //! this can be computationaly heavy.
    //!
    virtual void assemble_A();
    
    //!
    //! \fn bool is_host()
    //! \brief Check if the current process is the host (id=0)
    //!
    //! \return true if the process is the host
    //!
    bool is_host();
    
    //!
    //! \fn long long total_time(long long *t)
    //! \brief Sums and returns the execution time over all processes
    //!
    //! This function will take the time in parameters in a process then add
    //! it to the times in parameter of the method in other processes to finally
    //! return the total time.
    //!
    //! \param t: the local execution time
    //! \param comm: the MPI communicator
    //! \return the summed, total execution time
    //!
    long long total_time(long long *t, MPI_Comm comm);
    
    ////////////////////////////////////////////////////
    // MATRIX INPUT
    ////////////////////////////////////////////////////
    //!
    //! \fn bool take_A_value_loc(int m, int n, int i, bool local)
    //! \brief Check if the matrix value should be read by the current process
    //!
    //! This function makes sense for distributed matrices.
    //! From the row and column coordinates in the complete matrix and the index 
    //! in the list of non-zero values, this function will check if the current 
    //! process should read and save the value or not.
    //!
    //! \param m: row coordinate
    //! \param n: column coordinate
    //! \param i: index of the non-zero value
    //! \param local: true if the input matrix is distributed
    //! \return true if the value is to be read
    //!
    virtual bool take_A_value_loc(int m, int n, int i, bool local);
    
    //!
    //! \fn int read_nz_loc(int nz, bool local)
    //! \brief Returns the local number of non-zero values
    //!
    //! For distributed matrices, this function returns the local number
    //! of non-zero values (depending on the mapping returned by mumps, or the
    //! manual distribution you input). For centralized matrices, it only 
    //! returns the number of non-zero values.
    //!
    //! \param nz: global number of non-zero values
    //! \param local: true if the input matrix is distributed
    //! \return the local number of non-zeros
    //!
    virtual int read_nz_loc(int nz, bool local);

    //!
    //! \fn read_MM(std::string file, int &m, int &n, int &nz, 
    //!    double **values, std::vector<int**> indexes, bool n_present, bool rhs, 
    //!    bool local=false)
    //! \brief Read a matrix in Matrix Market format
    //!
    //! This function reads a matrix from a file in Matrix Market Format. It is
    //! generic to allow reading matrix A and right hand side b.
    //!
    //! \param file
    //! \param m: number of rows
    //! \param n: number of columns
    //! \param nz: number of non-zero values
    //! \param values: pointer to the array of non-zero values
    //! \param indexes: vector of pointers to the arrays of row and column index
    //! \param n_present: presence of the sizes (m,n,nz) in the file header
    //! \param rhs: true if we are reading a right hand side
    //! \param local: true if the input matrix is distributed
    //!
    void read_MM(std::string file, int &m, int &n, int &nz, 
        double **values, std::vector<int**> indexes, bool n_present, bool rhs, 
        bool local=false);
    
    //!
    //! \fn void read_A()
    //! \brief Read the input matrix
    //!
    //! Use the method read_MM(...) with the attribute _file_A in parameter. This
    //! will have a different behaviour depending on the input matrix:
    //!     - Distributed ?
    //!     - Host ?
    //!
    virtual void read_A() = 0;
    
    //!
    //! \fn void read_A_again()
    //! \brief deallocate and get the matrix A again using read_MM
    //!
    //! This function deallocate then reads again the matrix A. It is useful in 
    //! a context where:
    //!     - A is modified but should be used again
    //!     - A is not read at the solver initialization but later
    //!     - A has a special treatment such as distributed matrix,...
    //!
    virtual void read_A_again();
    
    //!
    //! \fn void read_b()
    //! \brief Read the input right hand side
    //!
    //! Use the method read_MM(...) with the attribute _file_b in parameter to
    //! read the right hand side on the host only. If the attribute _file_b is
    //! empty, the right hand side is initialized with all 1 values.
    //!
    virtual void read_b() = 0;
    
    //!
    //! \fn void read_b_again()
    //! \brief deallocate and get the matrix b again using read_b
    //!
    //! This function deallocate then reads again the right hand side b. It is 
    //! useful in a context where b is modified during a step (typically solve).
    //!
    virtual void read_b_again();

    
    ////////////////////////////////////////////////////
    // MATRIX OUTPUT
    ////////////////////////////////////////////////////
    //!
    //! \fn void display_ass(double values[], int n, std::vector<int*> indexes)
    //! \brief Display part of a matrix
    //!
    //! This function display a matrix. It is generic to allow displaying 
    //! matrix A, right hand side b, solution x,... Depending on the number of
    //! lines n to display, you can display whole or part of a matrix. The
    //! displayed lines are of the type:
    //!     - for matrix A: "row column value"
    //!     - for arrays (b, sol,..): "value"
    //!
    //! \param values: pointer to the array of non-zero values
    //! \param n: number of lines/non-zero values to display
    //! \param indexes: vector of pointers to the arrays of row and column index
    //!
    void display_ass(double values[], int n, std::vector<int*> indexes);
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display part of the matrix A
    //!
    //! This function display n lines/non-zero values of the system matrix A
    //! using the display_ass method. Per line: "row column value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_A(int n);
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole matrix A
    //!
    //! This function display all lines/non-zero values of the system matrix A
    //! using the display_ass method. Per line: "row column value"
    //!
    virtual void display_A();
    //!
    //! \fn void display_A_loc(int n)
    //! \brief Display part of the local part of the distributed matrix A
    //!
    //! This function display n lines/non-zero values of the local part of the 
    //! distributed matrix A using the display_ass method.
    //! Per line: "row column value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_A_loc(int n);
    
    //!
    //! \fn void display_A_loc(int n)
    //! \brief Display the local part of the distributed matrix A
    //!
    //! This function display all lines/non-zero values of the local part of the 
    //! distributed matrix A using the display_ass method.
    //! Per line: "row column value"
    //!
    virtual void display_A_loc();
    
    //!
    //! \fn void display_b(int n)
    //! \brief Display part of the right hand side b
    //!
    //! This function display n lines/non-zero values of the right hand side b
    //! using the display_ass method. Per line: "value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_b(int n);
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole right hand side b
    //!
    //! This function display all lines/non-zero values of the right hand side b
    //! using the display_ass method. Per line: "value"
    //!
    virtual void display_b();
    
    //!
    //! \fn void display_x(int n)
    //! \brief Display part of the solution array x
    //!
    //! This function display n lines/non-zero values of the solution array x
    //! using the display_ass method. Per line: "value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_x(int n);
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole solution array x
    //!
    //! This function display all lines/non-zero values of the solution array x
    //! using the display_ass method. Per line: "value"
    //!
    virtual void display_x();
    
    //!
    //! \fn void display_r(int n)
    //! \brief Display part of the residual array r=Ax-b
    //!
    //! This function display n lines/non-zero values of the residual array 
    //! r=Ax-b using the display_ass method. Per line: "value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_r(int n);
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole residual array r=Ax-b
    //!
    //! This function display all lines/non-zero values of the residual array 
    //! r=Ax-b using the display_ass method. Per line: "value"
    //!
    virtual void display_r();

    
    ////////////////////////////////////////////////////
    // RUNNING THE SOLVER
    ////////////////////////////////////////////////////    
    //!
    //! \fn void init()
    //! \brief Initialize the Solver
    //!
    //! Initialize the solver, for example:
    //!     - initialize data structures,
    //!     - initialize MPI processes,
    //!     - set default parameters,
    //!     - ...
    //!
    virtual void init() = 0;
    
    //!
    //! \fn void analyse()
    //! \brief Analysis phase of the solver
    //!
    //! This function will run the analysis of the solver using the data 
    //! structures in attributes.
    //!
    virtual void analyse() = 0;
    
    //!
    //! \fn bool read_A_before_facto()
    //! \brief Returns true if A should be read again before factorization
    //!
    //! This function will check if an option in the solver makes reading the
    //! matrix A again mandatory before a factorisation.
    //! For example, in Mumps, if the matrix is ditributed using a computed
    //! mapping or if the distributed parts are distributed before facto.
    //!
    //! \return true if the matrix A should be read again before facto
    //!
    virtual bool read_A_before_facto();
    
    //!
    //! \fn bool read_b_before_facto()
    //! \brief Returns true if b should be read again before factorization
    //!
    //! This function will check if an option in the solver makes reading the
    //! right hand side again mandatory before a factorisation.
    //! For example, in Mumps, the optional forward elimination of the right 
    //! hand side during factorization modifies b, thus it should be read again 
    //! before running a factorization again.
    //!
    //! \return true if the right hand side should be read again before facto
    //!
    virtual bool read_b_before_facto();
    
    //!
    //! \fn void factorize()
    //! \brief Factorization phase of the solver
    //!
    //! This function will run the factorization of the solver using the data 
    //! structures in attributes.
    //!
    virtual void factorize() = 0;
    
    //!
    //! \fn void solve()
    //! \brief Solve phase of the solver
    //!
    //! This function will run the solving of the solver using the data 
    //! structures in attributes.
    //!
    virtual void solve() = 0;
    
    //!
    //! \fn void metrics()
    //! \brief Computation of Metrics using QR_Mumps structures and functions
    //!
    //! This function will:
    //!     - assemble the matrix A
    //!     - initilize the QR_Mumps structure with data in attributes
    //!     - compute norms (residual, orthogonal residual, matrix A, right hand
    //! side, solution array)
    //!
    virtual void metrics() = 0;
    
    //!
    //! \fn void call()
    //! \brief Run at once analysis, factorization, solve and metrics.
    //!
    virtual void call();
    
    //!
    //! \fn void finalize()
    //! \brief Finalize the solver
    //!
    //! This function will finalize the solving phase. For example, it can:
    //!     - deallocate the right hand side,
    //!     - output the metrics specific to the system (A and b),
    //!     - call the finalization function of the solver,
    //!     - finalize the MPI processes,
    //!     - ...
    //!
    virtual void finalize() = 0;
    
    
    ////////////////////////////////////////////////////
    // (DE)ALLOCATION
    ////////////////////////////////////////////////////
    //!
    //! \fn void alloc_array()
    //! \brief Allocate an array of double
    //!
    //! This method allocate an array of double values. You can either copy
    //! another array of allocate an array of the same value.
    //!
    //! \param n: size of the array to allocate
    //! \param array: array to allocate
    //! \param copy_array: if true, array=copy of to_copy, else array=all value
    //! \param to_copy: array to be copied
    //! \param value: if not copy_array, the allocated array is filled of value
    //!
    void alloc_array(int n, double **array, bool copy_array=false,
        double *to_copy=NULL, double value=1);

    //!
    //! \fn void alloc_solve_residual()
    //! \brief Allocate the residual array
    //!
    virtual void alloc_solve_residual() = 0;
    
    //!
    //! \fn void alloc_rhs()
    //! \brief Allocate the right hand side array and initialize with values 1
    //!
    virtual void alloc_rhs() = 0;
    //!
    //! \fn void deallocate_A()
    //! \brief Deallocate A
    //!
    //! This function will deallocate the rows, columns and values arrays of A
    //!
    virtual void deallocate_A() = 0;
    
    //!
    //! \fn void deallocate_A_loc()
    //! \brief Deallocate the local part of A
    //!
    //! This function will deallocate the rows, columns and values arrays of the
    //! local part of A. It is only appropriate in a disributed context: does
    //! nothing by default.
    //!
    virtual void deallocate_A_loc();
    //!
    //! \fn void deallocate_b()
    //! \brief Deallocate b
    //!
    //! This function will deallocate the values array of b
    //!
    virtual void deallocate_b() = 0;
    
    
    ////////////////////////////////////////////////////
    // OUTPUTS
    ////////////////////////////////////////////////////
    //!
    //! \fn void base_output_metrics_init(std::string file)
    //! \brief Initialize the common solution specific metrics file
    //!
    //! This function will initialize the file for the solution specific metrics
    //! (norm of the residual, computation time,...) outputting the 
    //! tab-delimited metrics names in the file header. This method corresponds
    //! to metrics common to all solvers.
    //!
    //! \param file
    //!
    void base_output_metrics_init(std::string file);
    
    //!
    //! \fn void base_output_metrics(std::string sol_spec_file, long long ta = 0, 
    //!     long long tf = 0, long long ts = 0, long long ta_tot = 0, 
    //!     long long tf_tot = 0, long long ts_tot = 0, std::string key = "",
    //!     std::string value = "")
    //! \brief Output the common solution specific metrics
    //!
    //! This function will output all the metrics specific to the solution found
    //! in a tab-delimited form in the file. The computation times are given
    //! by the Benchmark class calling this function. This method corresponds
    //! to metrics common to all solvers.
    //!
    //! \param sol_spec_file
    //! \param ta: analysis computation time
    //! \param tf: factorization computation time
    //! \param ts: solve computation time
    //! \param key: key of the current option key tested
    //! \param value: value of the current option tested
    //! \param comm: the MPI communicator
    //!
    void base_output_metrics(std::string sol_spec_file, MPI_Comm comm, 
        long long ta = 0, long long tf = 0, long long ts = 0);
    //!
    //! \fn void output_metrics_init(std::string file)
    //! \brief Initialize the solution specific metrics file
    //!
    //! This function will initialize the file for the solution specific metrics
    //! (norm of the residual, computation time,...) outputting the 
    //! tab-delimited metrics names in the file header.
    //!
    //! \param file
    //!
    virtual void output_metrics_init(std::string file) = 0;
    
    //!
    //! \fn void output_metrics(std::string sol_spec_file, long long ta = 0, 
    //!     long long tf = 0, long long ts = 0, long long ta_tot = 0, 
    //!     long long tf_tot = 0, long long ts_tot = 0, std::string key = "",
    //!     std::string value = "")
    //! \brief Output the solution specific metrics
    //!
    //! This function will output all the metrics specific to the solution found
    //! in a tab-delimited form in the file. The computation times are given
    //! by the Benchmark class calling this function.
    //!
    //! \param sol_spec_file
    //! \param ta: analysis computation time
    //! \param tf: factorization computation time
    //! \param ts: solve computation time
    //! \param key: key of the current option key tested
    //! \param value: value of the current option tested
    //!
    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
        long long tf = 0, long long ts = 0, std::string key = "",
        std::string value = "") = 0;
    
    //!
    //! \fn void set_no_output()
    //! \brief Sets options so there is no output during computation.
    //!
    virtual void set_no_output() = 0;
    
//    
//    template <typename T>
//    void Solver::out(std::ostream str, T t) 
//    {
//        str << t << "\n";
//    }
//    template<typename T, typename... Args>
//    void Solver::out(std::ostream str, T t, Args... args) // recursive variadic function
//    {
//        str << t << "\n";
//        out(args...) ;
//    }
//    template<typename... Args>
//    void Solver::to_clog(Args... args) // recursive variadic function
//    {
//        if (is_host())
//            out(std::clog, args...);
//    }
//    template<typename... Args>
//    void Solver::to_cout(Args... args) // recursive variadic function
//    {
//        if (is_host())
//            out(std::cout, args...);
//    }
//    template<typename... Args>
//    void Solver::to_cerr(Args... args) // recursive variadic function
//    {
//        if (is_host())
//            out(std::cerr, args...);
//    }

    ////////////////////////////////////////////////////
    // GETTERS
    ////////////////////////////////////////////////////
    //!
    //! \fn int get_m()
    //! \brief get number of equations
    //!
    virtual int get_m() = 0;
    
    //!
    //! \fn int get_n()
    //! \brief get number of variables
    //!
    virtual int get_n() = 0;
    
    //!
    //! \fn int get_nz()
    //! \brief get number of non-zero values
    //!
    virtual int get_nz() = 0;
    
    //!
    //! \fn int get_nz_loc()
    //! \brief get number of local non-zero values
    //!
    virtual int get_nz_loc() = 0;
    
    //!
    //! \fn double* get_a()
    //! \brief get array of non-zero values
    //!
    virtual double* get_a() = 0;
    
    //!
    //! \fn double* get_a_loc()
    //! \brief get array of local non-zero values
    //!
    virtual double* get_a_loc() = 0;
    
    //!
    //! \fn int* get_irn()
    //! \brief get array of row indices
    //!
    virtual int* get_irn() = 0;
    
    //!
    //! \fn int* get_irn_loc()
    //! \brief get array of local row indices
    //!
    virtual int* get_irn_loc() = 0;
    
    //!
    //! \fn int* get_jcn()
    //! \brief get array of column indices
    //!
    virtual int* get_jcn() = 0;
    
    //!
    //! \fn int* get_jcn_loc()
    //! \brief get array of local column indices
    //!
    virtual int* get_jcn_loc() = 0;
    
    //!
    //! \fn double* get_b()
    //! \brief get array of right hand side
    //!
    virtual double* get_rhs() = 0;
    
    //!
    //! \fn double* get_x()
    //! \brief get array of the solution
    //!
    virtual double* get_x() = 0;
    
    //!
    //! \fn double* get_r()
    //! \brief get array of the residual
    //!
    virtual double* get_r() = 0;
};

#endif /* SOLVER_H */

