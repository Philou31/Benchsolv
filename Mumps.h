//! \file Mumps.h
//! \brief Mumps class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class Mumps inheriting from/implementing methods of class Solver
//!

#ifndef MUMPS_H
#define MUMPS_H

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define INFO(I) info[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFO(I) rinfo[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

#include "dmumps_c.h"
#include "Solver.h"
#include "Parameters_Mumps.h"
#include "Metrics.h"

class Mumps : public Solver {
private:
    DMUMPS_STRUC_C _id; // MUMPS data structure
    MPI_Comm _mpi_comm; // MPI communicator
    int _nb_procs = 1, _proc_id = 0; // #MPI and id MPI
    double *_r; // residual array
    Metrics _metrics;   // object computing metrics using QR_Mumps structure
    // Problem specific info
    std::string _pb_spec_file;  // File for the metrics specific to a matrix
    int _distr; // Input matrix distribution
    // File indicating the matrix PERSONAL distribution in the format: per line
    //      "chunk_id   beg_block   end_block   chunk_nz"
    //          chunk_id: from 0 to 127
    //          beg_block: number of the first block
    //          end_block: number of the last block
    //          chunk_nz: number of non-zero values in the chunk
    std::string _loc;
    int _loc_beg, _loc_end = 0; // start/end index of local distributed element
    // Type of distribution of the _loc file, can be either 1) per block of 
    // columns, 2) per block of arrowheads. The blocks are constructed so every
    // blocks contain about the same number of non-zero values.
    int _loc_option;    // Type of distribution (per row, per line,...)
    int _format;    // Input matrix format
    int _opt_key, _opt_value;   // Key and value of the option to test
    int *_mapping;  // Mapping array for distribution
    // Percentage relaxation added to the estimated memory
    int _mem_relax=parm::MEMORY_DEFAULT_PERCENT_INC;
    // Estimated memory + relaxation will be mult by this
    float _mem_factor=parm::MEMORY_SIZE_FACTOR;
public:
    //!
    //! \brief Constructor of the Mumps class
    //!
    //! \param test_id, file_A, n_present_A, file_b, n_present_b, nrows, ncol, nz: class Solver
    //! \param par: host working or not (see PAR in Parameters_Mumps.h)
    //! \param sym: matrix symetric, unsymetric or definit positive (see SYM)
    //! \param distr: Input matrix distribution (See Distribution)
    //! \param loc: File indicating the matrix distribution
    //! \param loc_option: Type of distribution (per row, per line,...)
    //! \param format: Input matrix format
    //! \param comm: integer describing the MPI communicator
    //! \param mpi_comm: MPI communicator
    //! \param pb_spec_file: File for the metrics specific to a matrix
    //! \param int_opt_key, int_opt_value: Key and value of the option to test
    //! \param mem_relax: Percentage relaxation added to the estimated memory
    //! \param mem_factor: Estimated memory + relaxation will be mult by this
    //! \return the Mumps instance
    //!
    Mumps(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, int par, int sym, int distr, 
        std::string loc, int loc_option, int format, int comm, MPI_Comm mpi_comm, 
        std::string pb_spec_file, int int_opt_key, int int_opt_value, int nrows,
        int ncols, int nz, int mem_relax, float mem_factor);
    
    //!
    //! \brief Destructor of the Mumps class
    //!
    ~Mumps();
    
    ////////////////////////////////////////////////////
    // MPI communication methods
    ////////////////////////////////////////////////////
    //!
    //! \fn assemble_A()
    //! \brief Assemble the matrix A on the host from the local parts
    //!
    //! This function will assemble on the host the local parts of the matrix A
    //! to obtain the global matrix arrays (rows, columns, values). Caution:
    //! this can be computationaly heavy.
    //!
    void assemble_A();
    
    virtual bool is_host() override;
    virtual long long total_time(long long *t) override;
    
    ////////////////////////////////////////////////////
    // MATRIX INPUT
    ////////////////////////////////////////////////////
    //!
    //! \fn parse_loc_file()
    //! \brief Parse a file of distribution
    //!
    //! This function parses a file indicating the matrix distribution in the 
    //! per line format: "chunk_id   beg_block   end_block   chunk_nz"
    //!      chunk_id: from 0 to 127
    //!      beg_block: number of the first block
    //!      end_block: number of the last block
    //!      chunk_nz: number of non-zero values in the chunk
    //! First find the matrix distribution which can be: 0) per block of rows, 
    //! 1) per block of columns, 2) per block of arrowheads.
    //!
    //! Considering there is a power of 2 MPI processes, a process will take its
    //! fair share of the chunks (with 4 processes, the first will take care of
    //! the 32 first chunks). As a result, each process gets the first row index
    //! and the last row index corresponding to the local part of the process.
    //!
    void parse_loc_file();
    
    //!
    //! \fn nz_mapping()
    //! \brief Get the local number of non-zeros and the mapping array
    //!
    //! Mumps can output a mapping array on the host during analysis if the good 
    //! distribution is chosen: for each non-zero value it specifies on which
    //! MPI process it should be handled.
    //! This function allows all processes to get the mapping array and compute
    //! the local number of non-zero values.
    //!
    int nz_mapping();
    
    //!
    //! \fn void get_A_again()
    //! \brief deallocate then get the local matrix A using get_MM
    //!
    //! This function should only be called with mapping or facto distribution:
    //! the values of the local matrices should be input locally before 
    //! factorisation using get_A_loc
    //!
    //! \param key
    //! \param value
    //! \param sol_spec_file: File for the metrics specific to the solution
    //!
    virtual void get_A_again() override;
    
    //!
    //! \fn void get_A_loc()
    //! \brief read the local part of matrix A in the distributed case
    //!
    //! This function will read the local part of the matrix A depending on the 
    //! chosen distribution
    //! TODO DOCUMENTATION
    //!
    void get_A_loc();
    
    virtual bool take_A_value_loc(int m, int n, int i, bool local) override;
    virtual int nz_loc(int nz, bool local) override;
    virtual void get_A() override;
    virtual void get_b() override;
    
    ////////////////////////////////////////////////////
    // MATRIX OUTPUT
    ////////////////////////////////////////////////////
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
    void display_A_loc(int n);
    
    //!
    //! \fn void display_A_loc(int n)
    //! \brief Display the local part of the distributed matrix A
    //!
    //! This function display all lines/non-zero values of the local part of the 
    //! distributed matrix A using the display_ass method.
    //! Per line: "row column value"
    //!
    void display_A_loc();
    
    virtual void display_A(int n) override;
    virtual void display_A() override;
    virtual void display_b(int n) override;
    virtual void display_b() override;
    
    ////////////////////////////////////////////////////
    // RUNNING THE SOLVER
    ////////////////////////////////////////////////////
    //!
    //! \fn void mumps(int job)
    //! \brief Calls the mumps main function with the specified job option
    //!
    //! \param job: option to call MUMPS (between -2 and 6)
    //!
    void mumps(int job);
    
    //!
    //! \fn void set_opt(int key, int value)
    //! \brief Sets the option key to the specified value
    //!
    //! \param key
    //! \param value
    //!
    void set_opt(int key, int value);
    
    //!
    //! \fn void set_opt(int key, int value, std::string sol_spec_file)
    //! \brief Sets the option key to the specified value, output option in file
    //!
    //! This function uses set_opt(int key, int value) to set the option key
    //! to the value and output "Setting key to value" in the file with metrics
    //! specific to the solution found.
    //!
    //! \param key
    //! \param value
    //! \param sol_spec_file: File for the metrics specific to the solution
    //!
    void set_opt(int key, int value, std::string sol_spec_file);
    
    virtual void init() override;
    virtual void analyse() override;
    virtual bool get_b_before_facto() override;
    virtual void factorize() override;
    virtual void solve() override;
    virtual void metrics() override;
    virtual void call() override;
    virtual void finalize() override;
    
    ////////////////////////////////////////////////////
    // (DE)ALLOCATION
    ////////////////////////////////////////////////////
    virtual void alloc_solve_residual() override;
    virtual void alloc_rhs() override;
    virtual void deallocate_A() override;
    
    ////////////////////////////////////////////////////
    // OUTPUTS
    ////////////////////////////////////////////////////
    
    //!
    //! \fn void problem_spec_metrics_output()
    //! \brief Output the metrics specific to the system (A, b)
    //!
    //! This function will output all the metrics specific to the system (A, b)
    //! in a tab-delimited form in the file. The computation times are given
    //! by the Benchmark class calling this function.
    //!
    void problem_spec_metrics_output();
    
    virtual void output_metrics_init(std::string file) override;
    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
        long long tf = 0, long long ts = 0, std::string key = "",
        std::string value = "") override;
    virtual void set_no_output() override;
};

#endif /* MUMPS_H */

