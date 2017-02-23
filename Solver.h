//! \file Solver.h
//! \brief Abstract class for the solvers classes
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! This is an abstract class (some methods implemented) for the solver classes 
//! implemented by Mumps, QR_Mumps,... It implements only common functions:
//!     - Input matrix: get_MM
//!     - Output matrix: display_ass
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
    
    
    ////////////////////////////////////////////////////
    // MPI communication methods
    ////////////////////////////////////////////////////
    //!
    //! \fn bool is_host()
    //! \brief Check if the current process is the host (id=0)
    //!
    //! \return true if the process is the host
    //!
    virtual bool is_host();
    
    //!
    //! \fn long long total_time(long long *t)
    //! \brief Sums and returns the execution time over all processes
    //!
    //! This function will take the time in parameters in a process then add
    //! it to the times in parameter of the method in other processes to finally
    //! return the total time.
    //!
    //! \param t: the local execution time
    //! \return the summed, total execution time
    //!
    virtual long long total_time(long long *t);
    
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
    //! \fn int nz_loc(int nz, bool local)
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
    virtual int nz_loc(int nz, bool local);

    //!
    //! \fn get_MM(std::string file, int &m, int &n, int &nz, 
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
    void get_MM(std::string file, int &m, int &n, int &nz, 
        double **values, std::vector<int**> indexes, bool n_present, bool rhs, 
        bool local=false);
    
    //!
    //! \fn void get_A()
    //! \brief Read the input matrix
    //!
    //! Use the method get_MM(...) with the attribute _file_A in parameter. This
    //! will have a different behaviour depending on the input matrix:
    //!     - Distributed ?
    //!     - Host ?
    //!
    virtual void get_A() = 0;
    
    //!
    //! \fn void get_A_again()
    //! \brief deallocate and get the matrix A again using get_MM
    //!
    //! This function deallocate then reads again the matrix A. It is useful in 
    //! a context where:
    //!     - A is modified but should be used again
    //!     - A is not read at the solver initialization but later
    //!     - A has a special treatment such as distributed matrix,...
    //!
    //! \param key
    //! \param value
    //! \param sol_spec_file: File for the metrics specific to the solution
    //!
    virtual void get_A_again();
    
    //!
    //! \fn void get_b()
    //! \brief Read the input right hand side
    //!
    //! Use the method get_MM(...) with the attribute _file_b in parameter to
    //! read the right hand side on the host only. If the attribute _file_b is
    //! empty, the right hand side is initialized with all 1 values.
    //!
    //!
    virtual void get_b() = 0;

    
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
    virtual void display_A(int n) = 0;
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole matrix A
    //!
    //! This function display all lines/non-zero values of the system matrix A
    //! using the display_ass method. Per line: "row column value"
    //!
    virtual void display_A() = 0;
    
    //!
    //! \fn void display_b(int n)
    //! \brief Display part of the right hand side b
    //!
    //! This function display n lines/non-zero values of the right hand side b
    //! using the display_ass method. Per line: "value"
    //!
    //! \param n: number of lines/non-zero values to display
    //!
    virtual void display_b(int n) = 0;
    
    //!
    //! \fn void display_A(int n)
    //! \brief Display the whole right hand side b
    //!
    //! This function display all lines/non-zero values of the right hand side b
    //! using the display_ass method. Per line: "value"
    //!
    virtual void display_b() = 0;
    
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
    //! \fn bool get_b_before_facto()
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
    virtual bool get_b_before_facto() = 0;
    
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
    virtual void call() = 0;
    
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
    //! \fn void deallocate_A()
    //! \brief Deallocate A
    //!
    //! This function will deallocate the rows, columns and values arrays of A
    //!
    virtual void deallocate_A() = 0;
    
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
    
    
    ////////////////////////////////////////////////////
    // OUTPUTS
    ////////////////////////////////////////////////////
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
};

#endif /* SOLVER_H */

