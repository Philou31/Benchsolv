//! \file functions.cpp
//! \brief Functions for the project
//! \author filou
//! \version 0.1
//! \date 13/11/2016, 18:42
//!

#include "Solver.h"

Solver::Solver(std::string test_id, std::string file_A, bool n_present_A, 
    std::string file_b, bool n_present_b, int nrows, int ncols, int nz):
    _test_id(test_id), _file_A(file_A), _n_present_A(n_present_A), 
    _file_b(file_b), _n_present_b(n_present_b), _nrows(nrows), _ncols(ncols),
    _nz(nz)
    {}

////////////////////////////////////////////////////
// MPI communication methods
////////////////////////////////////////////////////
bool Solver::is_host() {
    // No MPI in QR_Mumps
    return true;
}

long long Solver::total_time(long long *t) {
    return *t;
}
    
////////////////////////////////////////////////////
// MATRIX INPUT
////////////////////////////////////////////////////
bool Solver::take_A_value_loc(int m, int n, int i, bool local) {
    return true;
}

int Solver::nz_loc(int nz, bool local) {
    return nz;
}

void Solver::get_MM(std::string file, int &m, int &n, int &nz, double **values, 
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
    //      matrix distributed => read local part => size=local_nz
    int size = nz_loc(nz, local);
    //      rhs => size=number of columns
    if (rhs) size = m;
    
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

void Solver::get_A_again() {
    std::cout << "\nNo need to read the matrix A again.\n";
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

void Solver::display_x(int n) {
    std::cerr << "display_x: Not available with this solver.\n";
}
void Solver::display_x() {
    std::cerr << "display_x: Not available with this solver.\n";
}
void Solver::display_r(int n) {
    std::cerr << "display_r: Not available with this solver.\n";
}
void Solver::display_r() {
    std::cerr << "display_r: Not available with this solver.\n";
}