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

#include "Solver.h"
#include "Mumps.h"
#include <typeinfo>

Solver::Solver(std::string test_id, std::string file_A, bool n_present_A, 
    std::string file_b, bool n_present_b, int nrows, int ncols, int nz):
    _test_id(test_id), _file_A(file_A), _n_present_A(n_present_A), 
    _file_b(file_b), _n_present_b(n_present_b), _nrows(nrows), _ncols(ncols),
    _nz(nz)
    {}

void Solver::base_destruct() {
    if (is_host()) {
        deallocate_A();
        deallocate_b();
    }
    finalize();
}

////////////////////////////////////////////////////
// MPI communication methods
////////////////////////////////////////////////////
//void Solver::assemble_MM(int nz, double **values, std::vector<int**> indexes,
//        int nz_loc, double **values_loc, std::vector<int**> indexes_loc,
//        MPI_Comm comm) {
//    // Initialize to NULL for non host processes
//    for(std::vector<int**>::iterator it = indexes.begin(); 
//            it != indexes.end(); ++it) {
//        **it=NULL;
//    }
//    *values=NULL;
//    // Gather all local nz in an array
//    int *recvcounts=NULL;
//    if (is_host())
//        recvcounts = new int[_nb_procs];
//    MPI_Gather(&nz_loc, 1, MPI_INT, recvcounts, 1, MPI_INT, cst::HOST_ID, comm);
//    // The displacements of the gather are the same as the receive counts
//    int *displs=NULL;
//    if (is_host()) {
//        displs = new int[_nb_procs];
//        int i;
//        for(i=0;i<_nb_procs;++i) {
//            if (i>0)
//                displs[i]+=recvcounts[i-1];
//        }
//        // initialize A arrays (line, col, val)
//        for(std::vector<int**>::iterator it = indexes.begin(); 
//                it != indexes.end(); ++it) {
//            **it=new int[nz];
//        }
//        *values=new double[nz];
//    }
//    // Gather the matrix arrays
//    int ind1=0;
//    int ind2=0;
//    for(std::vector<int**>::iterator it = indexes.begin(); 
//            it != indexes.end(); ++it) {
//        for(std::vector<int**>::iterator it_loc = indexes_loc.begin(); 
//                it_loc != indexes_loc.end(); ++it_loc) {
//            if (ind1 == ind2) {
//                std::cout << typeid (**it).name() << "IT\n";
//                std::cout << typeid (*values).name() << "VALUES\n";
//                std::cout << nz << "\t" << nz_loc << "\t" << _nb_procs << "\n";
//                std::clog << "MWAHAHAHAHAHA\n";
//                display_A(10);
//                display_ass(*values_loc, 10, {});
//                MPI_Gatherv(**it_loc, nz_loc, MPI_INT, **it, recvcounts, displs,
//                    MPI_INT, cst::HOST_ID, comm);
//            }
//            ++ind2;
//        }
//        ++ind1;
//    }
//    std::clog << "MWAHAHAHAHAHA2\n";
//    MPI_Gatherv(*values_loc, nz_loc, MPI_DOUBLE, *values, recvcounts, displs,
//        MPI_DOUBLE, cst::HOST_ID, comm);
//    _A_assembled=true;
//    if (is_host()) {
//        delete[] recvcounts;
//        delete[] displs;
//    }
//}

void Solver::assemble_A() {}

bool Solver::is_host() {
    return _proc_id == cst::HOST_ID;
}

long long Solver::total_time(long long *t, MPI_Comm comm) {
    long long t_tot = 0;
    if (_nb_procs > 1)
        MPI_Reduce(t, &t_tot, 1, MPI_LONG_LONG, MPI_SUM, cst::HOST_ID, comm);
    else t_tot=*t;
    return t_tot;
}

////////////////////////////////////////////////////
// MATRIX INPUT
////////////////////////////////////////////////////
void Solver::base_construct() {
    init();
    read_A();
    read_b();
}

void Solver::init_MPI(MPI_Comm &comm) {
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen = 0;
    int ierr = MPI_Init (NULL, NULL);
    std::clog << "MPI initialization in Mumps, error code: " << ierr << "\n";
    ierr = MPI_Comm_size(comm, &_nb_procs);
    std::clog << "Running on CPU " << sched_getcpu() << " on node " << processor_name << "\n";
    ierr = MPI_Comm_rank(comm, &_proc_id);
    ierr = MPI_Get_processor_name(processor_name, &namelen);
    std::clog << "Proc id " << _proc_id << " on a total of " << _nb_procs << " procs\n\n";
}

void Solver::finalize_MPI() {
    int ierr = MPI_Finalize();
    std::cout << "MPI finalization in Mumps, error code: " << ierr << "\n";
}
    
void Solver::init_OpenMP() {
    #pragma omp parallel
    {
        /* Obtain and print thread id */
       	_tid = omp_get_thread_num();
        std::clog << "Initialisation of OpenMP thread = " << _tid << " on cpu " << sched_getcpu() << "\n";
        /* Only master thread does this */
        _nthreads = omp_get_num_threads();
        if (_tid == 0) {
            std::clog << "Number of OpenMP threads = " << _nthreads << "\n";
        }  /* All threads join master thread and terminate */
    }
}
    
bool Solver::take_A_value_loc(int m, int n, int i, bool local) {
    return true;
}

int Solver::read_nz_loc(int nz, bool local) {
    return nz;
}

void Solver::read_MM(std::string file, int &m, int &n, int &nz, double **values, 
        std::vector<int**> indexes, bool n_present, bool rhs, bool local) {
    std::ifstream infile(file.c_str());
    std::cout << "\nGetting matrix in file: " << file << "\n";

    // If we couldn't open the output file stream for reading
    if (!infile) {
        // Print an error and exit
        std::cerr << "Uh oh, " << file << " could not be opened for reading!\n";
        exit(1);
    }

    std::string strInput;
    // Read banner and m, n, nz line
    while (infile) {
        std::getline(infile, strInput);
        std::stringstream stream(strInput);
        // Discard commented lines
        if (strInput[0] == '%') {
            continue;
        }
        // Only read the size if it is present and we are not reading b
        else if (n_present) {
            if (!rhs) {
                stream >> m;
                stream >> n;
                stream >> nz;
                std::clog << "m: " << m << ", n: " << n << ", nz: " << nz << "\n";
            }
            _nrows=m; _ncols=n; _nz=nz;
            n_present = false;
            continue;
        } else
            _nrows=m; _ncols=n; _nz=nz;
            break;
    }

    // Real size of the matrix to read:
    int size = nz;
    //      rhs => size=number of columns
    if (rhs) size = m;
    //      matrix distributed => read local part => size=local_nz
    else size = read_nz_loc(nz, local);
        
    
    std::clog << "Matrix of size: " << size << "\n";
    
    // Initialize data structures in 3 arrays (matrix) or 1 array (rhs))
    for(std::vector<int**>::iterator it = indexes.begin();
            it != indexes.end(); ++it) {
        **it = new int[size];
    }
    *values = new double[size];
    
    int iii(0), jjj(0), kkk(0);
    // While there's still stuff left to read
    unsigned int i;
    int mn[2];
    while (infile) {
        // each line contains: row column value
        std::stringstream stream(strInput);
        // read row and column
        mn[0] = 0;
        mn[1] = 0;
        for(i = 0; i < indexes.size(); ++i) {
            stream >> mn[i];
        }
        // Do we take the values ? Only if rhs or depending on coordinates
        if (rhs || take_A_value_loc(mn[0], mn[1], jjj, local)) {
            i = 0;
            // Save row and column
            for(std::vector<int**>::iterator it = indexes.begin(); 
                    it != indexes.end(); ++it) {
                *((**it)+iii) = mn[i];
                ++i;
            }
            // Save value
            stream >> *(*values+iii);
            ++iii;
            // On the host: output every 10% of the matrix read
            if (size >= 10000)
            if (is_host() && iii%(size/10)==0) {
                ++kkk;
                std::clog << kkk*10 << "% loaded...\n";
            }
        }
        // Get next line
        std::getline(infile, strInput);
        ++jjj;
    }
    // Finalize
    std::clog << "Done\n";
    infile.close();
}

void Solver::read_A_again() {
    std::cout << "\nNo need to read the matrix A again.\n";
    return;
}

void Solver::read_b_again() {
    std::cout << "\nNo need to read the right hand side b again.\n";
    return;
}

////////////////////////////////////////////////////
// MATRIX OUTPUT
////////////////////////////////////////////////////
void Solver::display_ass(double values[], int n, std::vector<int*> indexes) {
    for(int jjj = 0; jjj < n; ++jjj) {
        // Display one line: "Row Column Value"
        for(std::vector<int*>::iterator it = indexes.begin(); 
                it != indexes.end(); ++it) {
            // A pas peur !! it points to element of the vector, *it is the 
            // array (implicit pointer), (*it)+j points to jth element of array
            std::clog << *((*it)+jjj) << " ";
        }
        std::clog << values[jjj] << "\n";
    }
}

void Solver::display_A_loc(int n) {
    if (is_host()) display_ass(get_a_loc(), n, {get_irn_loc(), get_jcn_loc()});
}

void Solver::display_A_loc() {
    if (is_host()) display_ass(get_a_loc(), get_nz_loc(), {get_irn_loc(), 
        get_jcn_loc()});
}

void Solver::display_A(int n) {
    if (is_host()) display_ass(get_a(), n, {get_irn(), get_jcn()});
}

void Solver::display_A() {
    if (is_host()) display_ass(get_a(), get_nz(), {get_irn(), get_jcn()});
}

void Solver::display_b(int n) {
    if (is_host()) display_ass(get_rhs(), n, {});    
}

void Solver::display_b() {
    if (is_host()) display_ass(get_rhs(), get_n(), {});
}

void Solver::display_r(int n) {
    if (is_host()) display_ass(get_r(), n, {});    
}

void Solver::display_r() {
    if (is_host()) display_ass(get_r(), get_n(), {});
}

void Solver::display_x(int n) {
    if (is_host()) display_ass(get_x(), n, {});    
}

void Solver::display_x() {
    if (is_host()) display_ass(get_x(), get_n(), {});
}

////////////////////////////////////////////////////
// RUNNING THE SOLVER
////////////////////////////////////////////////////
bool Solver::read_b_before_facto() {
    return false;
}
void Solver::call() {
    analyse();
    factorize();
    solve();
    metrics();
}
    
////////////////////////////////////////////////////
// (DE)ALLOCATION
////////////////////////////////////////////////////
void Solver::alloc_array(int n, double **array, bool copy_array,
        double *to_copy, double value) {
    if (is_host()) {
        std::clog << "\nAllocation array of size: " << n << "\n";
        *array = new double[n];
        if (copy_array)
            for(int i = 0; i < n; i++) *(*array+i) = to_copy[i];
        else
            for(int i = 0; i < n; i++) *(*array+i) = value;
    }
}

void Solver::deallocate_A_loc() {}

////////////////////////////////////////////////////
// OUTPUTS
////////////////////////////////////////////////////
void Solver::base_output_metrics_init(std::string file) {
    if (is_host()) {
        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << "test_id\tfile_A\tsolver\t#procs\t#threads\toption\tvalue\tta"
            "\ttf\tts\tta_tot\ttf_tot\tts_tot\txnrm\trnrm\tonrm\t";
        myfile.close();
    }
}

void Solver::base_output_metrics(std::string file, MPI_Comm comm, long long ta, 
        long long tf, long long ts) {
    long long ta_tot=total_time(&ta, comm);
    long long tf_tot=total_time(&tf, comm);
    long long ts_tot=total_time(&ts, comm);
    if (is_host()) {
        std::cout << "\n#procs     =  " << _nb_procs << "\n" <<
            "#threads        =  " << _nthreads << "\n" <<
            "time for analysis     =  " << cst::TIME_RATIO*ta << "\n" <<
            "time for facto        =  " << cst::TIME_RATIO*tf << "\n" <<
            "time for solve        =  " << cst::TIME_RATIO*ts << "\n" <<
            "tot time for analysis =  " << cst::TIME_RATIO*ta_tot << "\n" <<
            "tot time for facto    =  " << cst::TIME_RATIO*tf_tot << "\n" <<
            "tot time for solve    =  " << cst::TIME_RATIO*ts_tot << "\n" <<
            "||A||                 =  " << _anrm << "\n" <<
            "||b||                 =  " << _bnrm << "\n" <<
            "||x||                 =  " << _xnrm << "\n" <<
            "||r||/||A||           =  " << _rnrm << "\n" <<
            "||A^tr||/||r||        =  " << _onrm << "\n";

        std::ofstream myfile;
        myfile.open(file.c_str(), std::ofstream::app);
        myfile << _test_id << "\t" << _file_A << "\tmumps\t" << _nb_procs << 
            "\t" << _nthreads << "\t" << cst::TIME_RATIO*ta << "\t" << 
            cst::TIME_RATIO*tf << "\t" << cst::TIME_RATIO*ts << "\t" << 
            cst::TIME_RATIO*ta_tot << "\t" << cst::TIME_RATIO*tf_tot << "\t" << 
            cst::TIME_RATIO*ts_tot << "\t" << _xnrm << "\t" << _rnrm << "\t" << 
            _onrm << "\t";
        myfile.close();
    }
}