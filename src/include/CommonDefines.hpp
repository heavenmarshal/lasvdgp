/**
 * @file CommonDefines.h
 * @author Robert Carnell
 * @copyright Copyright (c) 2014, Robert Carnell
 * 
 * @license <a href="http://www.gnu.org/licenses/lgpl.html">GNU Lesser General Public License (LGPL v3)</a>
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __COMMONDEFINES_HPP__
#define	__COMMONDEFINES_HPP__

#include <cstdlib>
#include <cmath>
#include <exception>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cfloat>
#include <climits>
#include <cstdio>
#include "matrix.hpp"
#include "order.hpp"
#include "CRandom.hpp"

#define PRINT_MACRO printf
#define ERROR_MACRO printf

#define PRINT_RESULT 0

/**
 * @namespace lhslib LHS c++ Library namespace
 */
namespace lhslib 
{
    /**
     * Latin hypercube sample algorithm with maximin criterion
     * @param n number of rows / samples in the lha
     * @param k number parameters / columns in the lhs
     * @param dup A factor that determines the number of candidate points used in the search.
     * @param result the result matrix
     * @param oRandom the random number stream
     */
    void maximinLHS(int n, int k, int dup, bclib::matrix<int> & result, 
            bclib::CRandom<double> & oRandom);
    /**
     * type of size type for use with bclib::matrix<T>
     * @note the type of the matrix (i.e. int) is irrelevant for size_type
     */
    typedef bclib::matrix<int>::size_type msize_type;
    /**
     * type of size type for use with std::vector<T>
     * @note the type of the vector (i.e. int) is irrelevant for size_type
     */
    typedef std::vector<int>::size_type vsize_type;
    
    /**
     * Create a random latin hypercube sample
     * @param n number of rows / samples in the lhs
     * @param k number parameters / columns in the lhs
     * @param bPreserveDraw should the order of the draw be preserved if less columns are selected
     * @param result the lhs
     * @param oRandom the random number stream
     */
}

#endif	/* COMMONDEFINES_H */
