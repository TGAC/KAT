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
#include <memory>
using std::shared_ptr;

#include <kat/sparse_matrix.hpp>

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
