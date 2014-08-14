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

#include <stdint.h>

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

#include "sparse_matrix.hpp"

typedef shared_ptr<SparseMatrix<uint64_t> > SM64;

class ThreadedSparseMatrix {
private:

    uint16_t width;
    uint16_t height;
    uint16_t threads;

    SM64 final_matrix;
    vector<SM64> threaded_matricies;

public:

    ThreadedSparseMatrix(uint16_t _width, uint16_t _height, uint16_t _threads) :
    width(_width), height(_height), threads(_threads) {
        final_matrix = SM64(new SparseMatrix<uint64_t>(width, height));
        threaded_matricies = vector<SM64>(threads);

        for (int i = 0; i < threads; i++) {
            threaded_matricies.push_back(SM64(new SparseMatrix<uint64_t>(width, height)));
        }
    }

    virtual ~ThreadedSparseMatrix() {
    }

    SM64 getFinalMatrix() {
        return final_matrix;
    }

    SM64 getThreadMatrix(uint16_t index) {
        return threaded_matricies[index];
    }

    SM64 mergeThreadedMatricies() {
        // Merge matrix
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < threads; k++) {
                    final_matrix->inc(i, j, threaded_matricies[k]->get(i, j));
                }
            }
        }

        return final_matrix;
    }

};
