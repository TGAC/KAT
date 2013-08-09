//  ********************************************************************
//  This file is part of KAT - the Kmer Analysis Toolkit.
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

#include "sparse_matrix.hpp"

template <class T>
class ThreadedSparseMatrix
{
private:

    uint16_t width;
    uint16_t height;
    uint16_t threads;

    SparseMatrix<T>    *final_matrix;
    SparseMatrix<T>    **threaded_matricies;

public:

    ThreadedSparseMatrix(uint16_t _width, uint16_t _height, uint16_t _threads) :
        width(_width), height(_height), threads(_threads)
    {
        final_matrix = new SparseMatrix<uint64_t>(width, height);
        threaded_matricies = new SparseMatrix<uint64_t>*[threads];

        for(int i = 0; i < threads; i++) {
            threaded_matricies[i] = new SparseMatrix<uint_t>(width, height);
        }
    }

    ~ThreadedSparseMatrix()
    {
        if(final_matrix)
            delete final_matrix;

        final_matrix = NULL;

        if (threaded_matricies)
        {
            // No way of modifying threads after object initialisation so this is safe
            for(int i = 0; i < threads; i++)
            {
                if (threaded_matricies[i])
                {
                    delete threaded_matricies[i];
                    threaded_matricies[i] = NULL;
                }
            }

            delete threaded_matricies;
            threaded_matricies = NULL;
        }
    }


    SparseMatrix<T>* getFinalMatrix()
    {
        return final_matrix;
    }

    SparseMatrix<T>* getThreadMatrix(uint16_t index)
    {
        return threaded_matricies[index];
    }


    SparseMatrix<T>* mergeThreadedMatricies()
    {
        // Merge matrix
        for(int i = 0; i < width; i++)
        {
            for(int j = 0; j < height; j++)
            {
                for(int k = 0; k < threads; k++)
                {
                    final_matrix->inc(i, j, threaded_matricies[k]->get(i, j));
                }
            }
        }

        return final_matrix;
    }

};
