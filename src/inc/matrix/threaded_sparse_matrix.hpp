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
