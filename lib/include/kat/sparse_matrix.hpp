//  ********************************************************************
//  This file is part of KAT - the K-mer Analysis Toolkit.
//
//  KAT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  KAT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with KAT.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;
using boost::lexical_cast;

#include <kat/str_utils.hpp>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;
using std::string;
using std::vector;
using std::map;

namespace kat{
    
template <class T>
class SparseMatrix {
public:
    typedef map<uint32_t, map<uint32_t, T> > mat_t;
    typedef typename mat_t::iterator row_iter;
    typedef map<uint32_t, T> col_t;
    typedef typename col_t::iterator col_iter;
    
    typedef boost::error_info<struct SparseMatrixError,string> SparseMatrixErrorInfo;
    struct SparseMatrixException: virtual boost::exception, virtual std::exception { };

    SparseMatrix() : SparseMatrix(0) {}
    
    SparseMatrix(uint32_t i) {
        m = i;
        n = i;
    }

    SparseMatrix(uint32_t i, uint32_t j) {
        m = i;
        n = j;
    }
    
    /**
     * @brief SparseMatrix Loads a sparse matrix from file.  NOTE, matrix must contain uint64_t!!
     * @param file_path path to the file containing the sparse matrix
     */
    SparseMatrix(const path& file_path) {
        ifstream infile;
        infile.open(file_path.c_str());

        string line("");

        uint32_t i = 0;
        while (!infile.eof()) {
            getline(infile, line);

            // Only do something if the start of the line isn't a #
            if (!line.empty() && line[0] != '#') {
                vector<uint64_t> parts = kat::splitUInt64(line, ' ');

                n = parts.size();

                for (uint32_t j = 0; j < parts.size(); j++) {
                    mat[i][j] = parts[j];
                }

                i++;
            }
        }

        infile.close();

        m = i;
    }

    inline
    const T& operator()(uint32_t i, uint32_t j) const {
        if (i >= m || j >= n) {
            BOOST_THROW_EXCEPTION(SparseMatrixException() << SparseMatrixErrorInfo(string(
                    "Requested coords exceed limits of matrix.  Coords: ") +
                    lexical_cast<string>(i) + "," + lexical_cast<string>(j) + ".  Limits: " +
                    lexical_cast<string>(m) + "," + lexical_cast<string>(n)));
        }
        return mat[i][j];
    }

    inline
    T operator()(uint32_t i, uint32_t j) {
        if (i >= m || j >= n) {
            BOOST_THROW_EXCEPTION(SparseMatrixException() << SparseMatrixErrorInfo(string(
                    "Requested coords exceed limits of matrix.  Coords: ") +
                    lexical_cast<string>(i) + "," + lexical_cast<string>(j) + ".  Limits: " +
                    lexical_cast<string>(m) + "," + lexical_cast<string>(n)));
        }
        return mat[i][j];
    }

    T inc(uint32_t i, uint32_t j, T val) {
        mat[i][j] += val;
        return mat[i][j];
    }

    T get(uint32_t i, uint32_t j) const {
        if (i >= m || j >= n) {
            BOOST_THROW_EXCEPTION(SparseMatrixException() << SparseMatrixErrorInfo(string(
                    "Requested coords exceed limits of matrix.  Coords: ") +
                    lexical_cast<string>(i) + "," + lexical_cast<string>(j) + ".  Limits: " +
                    lexical_cast<string>(m) + "," + lexical_cast<string>(n)));
        }
        
        map<uint32_t, map<uint32_t, uint64_t>>::const_iterator res1 = mat.find(i);
        
        if (res1 == mat.end()) {
            return 0;
        }
        else {
            map<uint32_t, uint64_t>::const_iterator res2 = res1->second.find(j); 
           
            if (res2 == res1->second.end()) {
                return 0;
            }
            else {
                return res2->second;
            }
        }
            
    }

    uint32_t width() const {
        return m;
    }

    uint32_t height() const {
        return n;
    }

    T getMaxVal() const {
        T maxVal = 0;

        for (uint32_t i = 0; i < m; i++) {
            for (uint32_t j = 0; j < n; j++) {
                T val = get(i,j);
                maxVal = maxVal < val ? val : maxVal;
            }
        }

        return maxVal;
    }

    vector<T> operator*(const vector<T>& x) { //Computes y=A*x
        if (this->m != x.size()) {
            BOOST_THROW_EXCEPTION(SparseMatrixException() << SparseMatrixErrorInfo(string(
                    "Incompatible matrix provided for multiplication.  Provided matrix size is ") +
                    lexical_cast<string>(x.size()) + " where as size of matrix managed by this object is " + lexical_cast<string>(m)));
        }

        vector<T> y(this->m);
        T sum;

        row_iter ii;
        col_iter jj;

        for (ii = this->mat.begin(); ii != this->mat.end(); ii++) {
            sum = 0;
            for (jj = (*ii).second.begin(); jj != (*ii).second.end(); jj++) {
                sum += (*jj).second * x[(*jj).first];
            }
            y[(*ii).first] = sum;
        }

        return y;
    }

    void printMat() const {
        row_iter ii;
        col_iter jj;
        for (ii = this->mat.begin(); ii != this->mat.end(); ii++) {
            for (jj = (*ii).second.begin(); jj != (*ii).second.end(); jj++) {
                cout << (*ii).first << ' ';
                cout << (*jj).first << ' ';
                cout << (*jj).second << endl;
            }
        }
        cout << endl;
    }
    
    void getRow(uint32_t row_idx, vector<T>& row) {
        for (uint32_t i = 0; i < this->height(); i++) {
            row.push_back(mat[i][row_idx]);
        }        
    }
    
    void getColumn(uint32_t col_idx, vector<T>& col) {
        for (uint32_t i = 0; i < this->width(); i++) {
            col.push_back(mat[col_idx][i]);
        }        
    }


    T sumColumn(uint32_t col_idx) {
        return sumColumn(col_idx, 0, this->width() - 1);
    }

    T sumColumn(uint32_t col_idx, uint32_t start, uint32_t end) {
        T sum = 0;
        for (uint32_t i = start; i <= end; i++) {
            sum += mat[col_idx][i];
        }

        return sum;
    }

    T sumRow(uint32_t row_idx) {
        return sumRow(row_idx, 0, this->height() - 1);
    }

    T sumRow(uint32_t row_idx, uint32_t start, uint32_t end) {
        T sum = 0;
        for (uint32_t i = start; i <= end; i++) {
            sum += mat[i][row_idx];
        }

        return sum;
    }

    void printMatrix(ostream &out) const {
        printMatrix(out, false);
    }

    void printMatrix(ostream &out, bool transpose) const {
        if (transpose) {
            // Transpose matrix
            for (uint32_t i = 0; i < n; i++) {

                out << get(0, i);

                for (uint32_t j = 0; j < m; j++) {
                    out << " " << get(j, i);
                }

                out << endl;
            }
        } else {
            for (uint32_t i = 0; i < m; i++) {
                out << get(i, 0);

                for (uint32_t j = 1; j < n; j++) {
                    out << " " << get(i, j);
                }

                out << endl;
            }
        }
    }

private:
    mat_t mat;
    uint32_t m;
    uint32_t n;
};

typedef SparseMatrix<uint64_t> SM64;

class ThreadedSparseMatrix {
private:

    uint16_t width;
    uint16_t height;
    uint16_t threads;

    SM64 final_matrix;
    vector<SM64> threaded_matricies;

public:

    ThreadedSparseMatrix() : ThreadedSparseMatrix(0, 0, 0) {};
    
    ThreadedSparseMatrix(uint16_t _width, uint16_t _height, uint16_t _threads) :
    width(_width), height(_height), threads(_threads) {
        final_matrix = SM64(width, height);
        threaded_matricies = vector<SM64>(threads);

        for (int i = 0; i < threads; i++) {
            threaded_matricies[i] = SM64(width, height);
        }
    }

    virtual ~ThreadedSparseMatrix() {
    }

    const SM64& getFinalMatrix() const {
        return final_matrix;
    }

    const SM64& getThreadMatrix(uint16_t index) const {
        return threaded_matricies[index];
    }

    const SM64& mergeThreadedMatricies() {
		// Merge matrix
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < threads; k++) {
                    final_matrix.inc(i, j, threaded_matricies[k].get(i, j));
                }
            }
        }

        return final_matrix;
    }
    
    uint64_t incTM(uint16_t index, size_t i, size_t j, uint64_t val) {
        return threaded_matricies[index].inc(i, j, val);
    }

};
}
