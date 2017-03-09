//! \file QR_Mumps.h
//! \brief Mumps class for the solver
//! \author filou
//! \version 0.1
//! \date 22/10/2016, 18:42
//!
//! Class QR_Mumps inheriting from/implementing methods of class Solver
//!

#ifndef QR_MUMPS_H
#define QR_MUMPS_H

#include "Solver.h"
extern "C" {
    #include "dqrm_c.h"
}
#include "Parameters_QR_Mumps.h"

class QR_Mumps: public Solver {
private:
    dqrm_spmat_type_c _id; // QR_Mumps data structure
    char _transp{'n'};  // The matrix is to be used transposed or not
    double *_x, *_b;   // solution, residual and right hand side arrays
    int _nrhs{1};   // number of right hand sides
    std::string _opt_key;   // Key of the option to test
    int _opt_value; // Value of the option to test
public:
    //!
    //! \brief Constructor of the QR_Mumps class
    //!
    //! \param test_id, file_A, n_present_A, file_b, n_present_b: class Solver
    //! \param string_opt_key, int_opt_value: Key, value of the option to test
    //! \return the Mumps instance
    //!
    QR_Mumps(std::string test_id, std::string file_A, bool n_present_A, 
        std::string file_b, bool n_present_b, std::string string_opt_key, 
        int int_opt_value, int nrows, int ncols, int nz);
    
    //!
    //! \brief Destructor of the QR_Mumps class
    //!
    ~QR_Mumps();
    
    ////////////////////////////////////////////////////
    // MATRIX INPUT
    ////////////////////////////////////////////////////
    virtual void read_A() override;
    virtual void read_b() override;
    
    ////////////////////////////////////////////////////
    // RUNNING THE SOLVER
    ////////////////////////////////////////////////////
    // SET OPT SHOULD BE THE SAME SCAFFOLD BETWEEN SOLVERS !!
    //!
    //! \fn void set_opt(int key, int value)
    //! \brief Sets the option key to the specified value
    //!
    //! \param key
    //! \param value
    //!
    void set_opt(std::string key, int value);
    
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
    void set_opt(std::string key, int value, std::string sol_spec_file);
    
    virtual void init() override;
    virtual void analyse() override;
    virtual void factorize() override;
    virtual void solve() override;
    virtual void metrics() override;
    virtual void finalize() override;
    
    ////////////////////////////////////////////////////
    // (DE)ALLOCATION
    ////////////////////////////////////////////////////
    virtual void alloc_solve_residual() override;
    virtual void alloc_rhs() override;
    virtual void deallocate_A() override;
    virtual void deallocate_b() override;
    
    ////////////////////////////////////////////////////
    // OUTPUTS
    ////////////////////////////////////////////////////
    virtual void output_metrics_init(std::string file) override;
    virtual void output_metrics(std::string sol_spec_file, long long ta = 0, 
        long long tf = 0, long long ts = 0, std::string key = "",
        std::string value = "") override;
    virtual void set_no_output() override;

    ////////////////////////////////////////////////////
    // GETTERS
    ////////////////////////////////////////////////////
    virtual int get_m() override;
    virtual int get_n() override;
    virtual int get_nz() override;
    virtual int get_nz_loc() override;
    virtual double* get_a() override;
    virtual double* get_a_loc() override;
    virtual int* get_irn() override;
    virtual int* get_irn_loc() override;
    virtual int* get_jcn() override;
    virtual int* get_jcn_loc() override;
    virtual double* get_rhs() override;
    virtual double* get_x() override;
    virtual double* get_r() override;
};


#endif /* QR_MUMPS_H */

