/*************************************************************************
eKATS.c 
Copyright (c) IRISA, Pierre-François Marteau
References:
[i] P. F. Marteau and S. Gibet, "On Recursive Edit Distance Kernels With Application to Time Series Classification," 
in IEEE Transactions on Neural Networks and Learning Systems, vol. 26, no. 6, pp. 1121-1133, June 2015. 
doi: 10.1109/TNNLS.2014.2333876
[ii] Pierre-François Marteau, Times series averaging and denoising from a probabilistic perspective on time-elastic kernels, International Journal of Applied Mathematics and Computer Science, De Gruyter, In press

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
**************************************************************************/
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <climits>

#define Point std::vector<double>
#define _KATS 0
#define _EIP 1
#define INF 1e20

namespace ekats2 {
    class eKats2 {
    public:

    double kdtw(std::vector<std::vector<double>> va, std::vector<std::vector<double>> vb, double nu);
    std::vector< std::vector<double> > interpolate(std::vector< std::vector<double> > ts);
    std::vector<std::vector<double>> iTEKA(std::vector<std::vector<double>> C, std::vector<std::vector<std::vector<double>> > dataset, double nu);
    std::vector<std::vector< std::vector<double> >> iTEKA_stdev(std::vector<std::vector<double>> C, std::vector<std::vector<std::vector<double>> > dataset, double nu);
    double SqEuclideanDistance(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y);
    };
};


