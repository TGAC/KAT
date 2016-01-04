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

#include <assert.h>
#include <vector>
using std::vector;

namespace kat {
    
/**
 * Abstract base class for calculating distance metrics between two arrays of big unsigned integers
 */
class DistanceMetric {
    
public:
    virtual string getName() = 0;
    virtual double calcDistance(const vector<uint64_t>& s1, const vector<uint64_t>& s2) = 0;
    
    virtual ~DistanceMetric() {}
};

class MinkowskiDistance : public DistanceMetric {
    
protected:
    int p; 
    
public:
    MinkowskiDistance() : MinkowskiDistance(1) {}  // Default to manhattan distance
    MinkowskiDistance(int power) {this->p = power;}
    
    virtual string getName() {return "Minkowski";}
    
    double calcDistance(const vector<uint64_t>& s1, const vector<uint64_t>& s2) {
    
        uint64_t sum = 0;

        for(size_t i = 0; i < s1.size(); i++) {
            uint64_t diff = s1[i] < s2[i] ? s2[i] - s1[i] : s1[i] - s2[i];
            sum += std::pow(diff, p);
        }

        return p == 1 ? (double)sum : std::pow((double)sum, 1.0 / (double)p);
    }
};

class ManhattanDistance : public MinkowskiDistance {
public:   
    ManhattanDistance() {this->p = 1;}
    
    string getName() {return "Manhattan";}    
};

class EuclideanDistance : public MinkowskiDistance {
public:   
    EuclideanDistance() {this->p = 2;}
    
    string getName() {return "Euclidean";}
};

class CosineDistance : public DistanceMetric {
    
    string getName() {return "Cosine";}
    
    double calcDistance(const vector<uint64_t>& s1, const vector<uint64_t>& s2) {
        double dot = 0.0, denom_a = 0.0, denom_b = 0.0 ;
        for(size_t i = 0; i < s1.size(); i++) {
            dot += s1[i] * s2[i] ;
            denom_a += std::pow(s1[i], 2);
            denom_b += std::pow(s2[i], 2);
        }
        return 1.0 -(dot / (std::sqrt(denom_a) * std::sqrt(denom_b))) ;
    }
};

class CanberraDistance : public DistanceMetric {
    
    string getName() {return "Canberra";}
    
    double calcDistance(const vector<uint64_t>& s1, const vector<uint64_t>& s2) {
        
        double sum = 0.0;
        for(size_t i = 0; i < s1.size(); i++) {
            double diff = (double)s1[i] - (double)s2[i];
            double sum_i = s1[i] + s2[i];
            if (sum_i > 0) sum += std::abs(diff) / sum_i;
        }
        
        return sum;
    }
};

class JaccardDistance : public DistanceMetric {
    
    string getName() {return "Jaccard";}
    
    double calcDistance(const vector<uint64_t>& s1, const vector<uint64_t>& s2) {
        
        double sum = 0.0;
        for(size_t i = 0; i < s1.size(); i++) {
            sum += std::min(s1[i],s2[i]);
        }
            
        double sum2 = 0.0;
        for(size_t i = 0; i < s1.size(); i++) {
            sum2 += std::max(s1[i],s2[i]);
        }
         
        return 1.0 - (sum / sum2);
    }
};

}

