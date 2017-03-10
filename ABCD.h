////! \file ABCD.h
////! \brief ABCD class for the solver
////! \author filou
////! \version 0.1
////! \date 22/10/2016, 18:42
////!
////! Class ABCD inheriting from/implementing methods of class Solver
////!
//
//#ifndef ABCD_H
//#define ABCD_H
//
//#include <boost/mpi.hpp>
//#include "abcd.h"
//#include "Solver.h"
//#include "Parameters_ABCD.h"
//#include "Metrics.h"
//
//class ABCD : public Solver {
//private:
//    abcd _id; // MUMPS data structure
//    MPI_Comm _mpi_comm; // MPI communicator
//    double *_r; // residual array
//    Metrics _metrics;   // object computing metrics using QR_Mumps structure
//    int _opt_key, _opt_value;   // Key and value of the option to test
//public:
//    //!
//    //! \brief Constructor of the ABCD class
//    //!
//    //! \param test_id, file_A, n_present_A, file_b, n_present_b, nrows, ncol, nz: class Solver
//    //! \param sym: true if the matrix is symmetric
//    //! \return the ABCD instance
//    //!
//    ABCD(std::string test_id, std::string file_A, bool n_present_A, 
//        std::string file_b, bool n_present_b, int int_opt_key, 
//            int int_opt_value, int nrows, int ncols, int nz, bool sym, 
//            MPI_Comm mpi_comm);
//    
//    //!
//    //! \brief Destructor of the ABCD class
//    //!
//    ~ABCD();
//    
//    ////////////////////////////////////////////////////
//    // MATRIX INPUT
//    ////////////////////////////////////////////////////
//    virtual void read_A() override;
//    virtual void read_b() override;
//    
//    ////////////////////////////////////////////////////
//    // RUNNING THE SOLVER
//    ////////////////////////////////////////////////////
//    //!
//    //! \fn void set_opt(int key, int value)
//    //! \brief Sets the option key to the specified value
//    //!
//    //! \param key
//    //! \param value
//    //!
//    void set_opt(int key, int value);
//    
//    //!
//    //! \fn void set_opt(int key, int value, std::string sol_spec_file)
//    //! \brief Sets the option key to the specified value, output option in file
//    //!
//    //! This function uses set_opt(int key, int value) to set the option key
//    //! to the value and output "Setting key to value" in the file with metrics
//    //! specific to the solution found.
//    //!
//    //! \param key
//    //! \param value
//    //! \param sol_spec_file: File for the metrics specific to the solution
//    //!
//    void set_opt(int key, int value, std::string sol_spec_file);
//    
//    virtual void init() override;
//    virtual void analyse() override;
//    virtual void factorize() override;
//    virtual void solve() override;
//    virtual void metrics() override;
//    virtual void finalize() override;
//    
//    ////////////////////////////////////////////////////
//    // (DE)ALLOCATION
//    ////////////////////////////////////////////////////
//    virtual void alloc_solve_residual() override;
//    virtual void alloc_rhs() override;
//    virtual void deallocate_A() override;
//    virtual void deallocate_b() override;
//    
//    ////////////////////////////////////////////////////
//    // OUTPUTS
//    ////////////////////////////////////////////////////
//    virtual void output_metrics_init(std::string file) override;
//    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
//        long long tf = 0, long long ts = 0, std::string key = "",
//        std::string value = "") override;
//    virtual void set_no_output() override;
//
//    ////////////////////////////////////////////////////
//    // GETTERS
//    ////////////////////////////////////////////////////
//    virtual int get_m() override;
//    virtual int get_n() override;
//    virtual int get_nz() override;
//    virtual int get_nz_loc() override;
//    virtual double* get_a() override;
//    virtual double* get_a_loc() override;
//    virtual int* get_irn() override;
//    virtual int* get_irn_loc() override;
//    virtual int* get_jcn() override;
//    virtual int* get_jcn_loc() override;
//    virtual double* get_rhs() override;
//    virtual double* get_x() override;
//    virtual double* get_r() override;
//};
//
//#endif /* ABCD_H */
//
