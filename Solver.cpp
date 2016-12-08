//! \file functions.cpp
//! \brief Functions for the project
//! \author filou
//! \version 0.1
//! \date 13/11/2016, 18:42
//!

#include <deque>

#include "Solver.h"

Solver::Solver(std::string test_id, std::string file_A, bool n_present_A, 
    std::string file_b, bool n_present_b):
    _test_id(test_id), _file_A(file_A), _n_present_A(n_present_A), 
    _file_b(file_b), _n_present_b(n_present_b)
    {}

void Solver::get_simple(int &m, int &n, int &nz, double **values, 
        int **irn, int **jcn, int &nrhs, int &lrhs, double **rhs) {
    m = 2;
    n = 2;
    nz = 2;
    *irn = new int[nz] {
        1, 2
    };
    *jcn = new int[nz] {
        1, 2
    };
    *values = new double[nz] {
        1.0, 2.0
    };

    nrhs = 1;
    lrhs = 2;
    *rhs = new double[lrhs] {
        1.0, 4.0
    };
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
        if (strInput[0] == '%') {
            continue;
        } else if (n_present) {
            if (!rhs) {
                stream >> m;
                stream >> n;
                stream >> nz;
                std::clog << "m: " << m << ", n: " << n << ", nz: " << nz << "\n";
            }
            n_present = false;
            continue;
        } else
            break;
    }

    int size = nz_loc(nz, local);
    if (rhs) size = m;
    
    std::clog << "Matrix of size: " << size << "\n";
    
    // Initialize matrix in 3 arrays
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
        mn[0] = 0;
        mn[1] = 0;
        for(i = 0; i < indexes.size(); ++i) {
            stream >> mn[i];
        }
        if (rhs || take_A_value_loc(mn[0], mn[1], jjj, local)) {
            i = 0;
            for(std::vector<int**>::iterator it = indexes.begin(); 
                    it != indexes.end(); ++it) {
                *((**it)+iii) = mn[i];
                ++i;
            }
            stream >> *(*values+iii);
            ++iii;
            if (is_host() && iii%(size/10)==0) {
                ++kkk;
                std::clog << kkk*10 << "% loaded...\n";
            }
        }
        std::getline(infile, strInput);
        ++jjj;
    }
    std::clog << "Done\n";
    infile.close();
}

void Solver::display_ass(double values[], int n, std::vector<int*> indexes) {
    for(int jjj = 0; jjj < n; ++jjj) {
        for(std::vector<int*>::iterator it = indexes.begin(); 
                it != indexes.end(); ++it) {
            // A pas peur !! it points to element of the vector, *it is the 
            // array (implicit pointer), (*it)+j points to jth element of array
            std::clog << *((*it)+jjj) << " ";
        }
        std::clog << values[jjj] << "\n";
    }
}