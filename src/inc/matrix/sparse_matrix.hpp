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

#include <str_utils.hpp>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;
using std::string;
using std::vector;
using std::map;

template <class T>
class SparseMatrix
{
public:
    typedef map<size_t, map<size_t , T> > mat_t;
    typedef typename mat_t::iterator row_iter;
    typedef map<size_t, T> col_t;
    typedef typename col_t::iterator col_iter;

    SparseMatrix(size_t i){ m=i; n=i; }
    SparseMatrix(size_t i, size_t j){ m=i; n=j; }

    /**
     * @brief SparseMatrix Loads a sparse matrix from file.  NOTE, matrix must contain uint64_t!!
     * @param file_path path to the file containing the sparse matrix
     */
    SparseMatrix(string file_path)
    {
        ifstream infile;
        infile.open(file_path.c_str());

        string line("");

        uint32_t i = 0;
        while(!infile.eof())
        {
            getline(infile, line);

            // Only do something if the start of the line isn't a #
            if (!line.empty() && line[0] != '#')
            {
                vector<uint64_t> parts = kat::splitUInt64(line, ' ');

                n = parts.size();

                for(uint32_t j = 0; j < parts.size(); j++)
                {
                    mat[i][j] = parts[j];
                }

                i++;
            }
        }

        infile.close();

        m = i;
    }

    inline
    T& operator()(size_t i, size_t j)
    {
        if(i>=m || j>=n) throw;
        return mat[i][j];
    }
    inline
    T operator()(size_t i, size_t j) const
    {
        if(i>=m || j>=n) throw;
        return mat[i][j];
    }

    T inc(size_t i, size_t j, T val)
    {
        mat[i][j] = mat[i][j] + val;
        return mat[i][j];
    }

    T get(size_t i, size_t j)
    {
        if(i>=m || j>=n) throw;
        return mat[i][j];
    }

    size_t width()   {return m;}
    size_t height()  {return n;}

    T getMaxVal()
    {
        T maxVal = 0;

        for(size_t i = 0; i < m; i++)
        {
            for(size_t j = 0; j < n; j++)
            {
                maxVal = maxVal < mat[i][j] ? mat[i][j] : maxVal;
            }
        }

        return maxVal;
    }

    std::vector<T> operator*(const std::vector<T>& x)
    {  //Computes y=A*x
        if(this->m != x.size()) throw;

        std::vector<T> y(this->m);
        T sum;

        row_iter ii;
        col_iter jj;

        for(ii=this->mat.begin(); ii!=this->mat.end(); ii++){
            sum=0;
            for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
                sum += (*jj).second * x[(*jj).first];
            }
            y[(*ii).first]=sum;
        }

        return y;
    }

    void printMat()
    {
        row_iter ii;
        col_iter jj;
        for(ii=this->mat.begin(); ii!=this->mat.end(); ii++){
            for( jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
                cout << (*ii).first << ' ';
                cout << (*jj).first << ' ';
                cout << (*jj).second << endl;
            }
        }
        cout << endl;
    }


    T sumColumn(size_t col_idx)
    {
        return sumColumn(col_idx, 0, this->width() -1);
    }

    T sumColumn(size_t col_idx, size_t start, size_t end)
    {
        T sum = 0;
        for(size_t i = start; i <= end; i++)
        {
            sum += mat[col_idx][i];
        }

        return sum;
    }


    T sumRow(size_t row_idx)
    {
        return sumRow(row_idx, 0, this->height() - 1);
    }

    T sumRow(size_t row_idx, size_t start, size_t end)
    {
        T sum = 0;
        for(size_t i = start; i <= end; i++)
        {
            sum += mat[i][row_idx];
        }

        return sum;
    }

    void printMatrix(ostream &out)
    {
        printMatrix(out, false);
    }

    void printMatrix(ostream &out, bool transpose)
    {
        if (transpose)
        {
            // Transpose matrix
            for(size_t i = 0; i < n; i++) {

                out << mat[0][i];

                for(size_t j = 0; j < m; j++) {
                    out << " " << mat[j][i];
                }

                out << endl;
            }
        }
        else
        {
            for(size_t i = 0; i < m; i++)
            {
                out << mat[i][0];

                for(size_t j = 1; j < n; j++)
                {
                    out << " " << mat[i][j];
                }

                out << endl;
            }
        }
    }
    
protected:
    SparseMatrix(){}

private:
    mat_t mat;
    size_t m;
    size_t n;
};
