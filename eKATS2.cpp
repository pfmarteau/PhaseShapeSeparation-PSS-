/*************************************************************************
eKATS.cpp
Copyright (c) IRISA, Pierre-François Marteau
References:
[i] P. F. Marteau and S. Gibet, "On Recursive Edit Distance Kernels With Application to Time Series Classification," 
in IEEE Transactions on Neural Networks and Learning Systems, vol. 26, no. 6, pp. 1121-1133, June 2015. 
doi: 10.1109/TNNLS.2014.2333876
[ii] Pierre-François Marteau, Times series averaging and denoising from a probabilistic perspective on time-elastic kernels, 
International Journal of Applied Mathematics and Computer Science, De Gruyter, In press

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
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <vector>
#include "LinearInterpolator.hpp"
#include <iostream>
#include <string>
#include <valarray>
#include <time.h>
#include <unistd.h>
#include "eKATS2.h"

using namespace ekats2;
using namespace std;

std::vector< std::vector<double> > reverseArray(std::vector< std::vector<double> > v){
    unsigned int l=v.size();
    unsigned int d=v[0].size();
    std::vector< std::vector<double> > vout(l, std::vector<double>(d));

    for(unsigned int i = 0; i <l; i++){
        for(unsigned int j = 0; j <d; j++){
            vout[i][j]=v[l-i-1][j];
        }
    }
    return vout;
}

double sq_error(std::vector<double> p1, std::vector<double> p2){
    double s=0.0, d;
    for(unsigned int i=0; i<p1.size(); i++){
        d=p1[i]-p2[i];
        s+=d*d;
    }
    return (s);
}

double eKats2::SqEuclideanDistance(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y) {
    size_t d = x.size();
    size_t d1 = y.size();
    if (d1!=d){
            std::cout << "SqEuclideanDistance: error, the two time series have not the same length" << std::endl;
            return 0;
            }
    size_t dim = x[0].size();
    size_t dim1 = y[0].size();
    if (dim1!=dim){
            std::cout << "SqEuclideanDistance: error, the two time series have not the same dimension" << std::endl;
            return 0;
            }
    double S=0.0, s;
    for (unsigned int i = 0; i < d; i++){
       s=0;
       for (unsigned int k = 0; k < dim; k++){
        s+=(x[i][k] - y[i][k])*(x[i][k] - y[i][k]);
    }
    S+=s;
    }
    return S;//s/sqrt(n1*n2);
}

double eKats2::kdtw(std::vector< std::vector<double> > va, std::vector< std::vector<double> > vb, double nu) {
    // I:dim: dimension of multivariate time series va and vb
    // I:la: length of time series va
    // I:va: multivariate time series va
    // I:lb: length of time series vb
    // I:vb: multivariate time series vb
    // I:nu: parameter of the local kernel (exp(-nu*delta(va(i),vb(j)))
    // O:return value, kdtw(va,vb,nu)
    int r = va.size();
    int c = vb.size();
    int dim = va[0].size();
    double dist;
    int i,j,k;
    double factor=1.0/3.0;
    double minprob=1e-300;
    double deg=1.0;
    nu=nu/dim;

    int rc=r;
    if (c>r)
        rc=c;

    double *D3 = (double *)calloc((unsigned int)rc+1, sizeof(double));
    double **D = (double **)calloc((unsigned int)r+1, sizeof(double*));
    double **D0 = (double **)calloc((unsigned int)r+1, sizeof(double*));
    for(i=0; i<=r; i++){
        D[i]=(double *)calloc((unsigned int)c+1, sizeof(double));
        D0[i]=(double *)calloc((unsigned int)c+1, sizeof(double));
    }

    for(i=1; i<=r; i++){
        for (j = 1; j <= c; j++) {
            dist = 0;
            for (k = 0; k < dim; k++)
                dist += (va[i - 1][k] - vb[j - 1][k]) * (va[i - 1][k] - vb[j - 1][k]);
            D[i][j] = factor*(exp(-nu*(pow(dist, deg)) )+minprob);
            if(i==j)
                D0[i][j]=D[i][j];
        }
    }

    for(j=1; j<=rc; j++) {
        dist=0;
        if(j<=r&&j<=c){
            for(k=0; k<dim; k++)
                dist+=(va[j-1][k]-vb[j-1][k])*(va[j-1][k]-vb[j-1][k]);
            D3[j]=factor*(exp(-nu*(pow(dist, deg)))+minprob);
        }
        else
            D3[j]=0;
    }

    D[0][0] = 1;
    D0[0][0] = 1;
    D3[0]=1;
    for (i = 1; i <= r; i++) {
        D[i][0] = (D[i-1][0])*D[i][1];
        D0[i][0] = (D0[i-1][0])*(D3[i]+D3[1])/2.0;
    }
    for (j = 1; j <= c; j++) {
        D[0][j] = (D[0][j-1])*D[1][j];
        D0[0][j] = (D0[0][j-1])*(D3[j]+D3[1])/2.0;
    }
    double dij;
     
    for (i = 1; i <= r; i++) {
        for (j = 1; j <= c; j++) {
            dij = D[i][j];
            D[i][j] = (D[i - 1][j - 1] + D[i - 1][j] + D[i][j - 1])* dij ;
            D0[i][j] = (D3[i] + D3[j])/2.0 *(D0[i - 1][j - 1] +  D0[i-1][j] + D0[i][j-1]);
        }
    }
     
    dist = D[r][c] + D0[r][c] ;

    for(i=0; i<=r; i++) {
        free(D[i]);
        free(D0[i]);
    }
    free(D0);
    free(D);
    free(D3);
    return dist;
}

std::vector< std::vector<double> > kdtw_mat(std::vector< std::vector<double> > va, std::vector< std::vector<double> > vb, double nu) {
    // I:dim: dimension of multivariate time series va and vb
    // I:la: length of time series va
    // I:va: multivariate time series va
    // I:lb: length of time series vb
    // I:vb: multivariate time series vb
    // I:nu: parameter of the local kernel (exp(-nu*delta(va(i),vb(j)))
    // O:return value, alignment matrix kdtw(va,vb,nu)
    int r = va.size();
    int c = vb.size();
    int dim = va[0].size();
    double dist;
    int i,j,k;
    double factor=1.0/3.0;
    double minprob=1e-300;
    double deg=1.0;
    nu=nu/dim;
    int rc=r;
    if (c>r)
        rc=c;

    std::vector<double>  D3=std::vector<double>(rc+1);
    std::vector< std::vector<double> > D(r+1, std::vector<double>(c+1));
    std::vector< std::vector<double> > D0(r+1, std::vector<double>(c+1));

    for(i=1; i<=r; i++){
        for (j = 1; j <= c; j++) {
            dist = 0;
            for (k = 0; k < dim; k++)
                dist += (va[i - 1][k] - vb[j - 1][k]) * (va[i - 1][k] - vb[j - 1][k]);
            D[i][j] = factor*(exp(-nu*(pow(dist, deg)))+minprob);
            if(i==j)
                D0[i][j]=D[i][j];
        }
    }
    for(j=1; j<=rc; j++) {
        dist=0;
        if(j<=r&&j<=c){
            for(k=0; k<dim; k++)
                dist+=(va[j-1][k]-vb[j-1][k])*(va[j-1][k]-vb[j-1][k]);
            D3[j]=factor*(exp(-nu*(pow(dist, deg)))+minprob);
        }
        else
            D3[j]=0;
    }

    D[0][0] = 1;
    D0[0][0] = 1;
    D3[0]=1;
    for (i = 1; i <= r; i++) {
        D[i][0] = D[i-1][0]*D[i][1];
        D0[i][0] = D0[i-1][0]*(D3[i]);
    }
    for (j = 1; j <= c; j++) {
        D[0][j] = D[0][j-1]*D[1][j];
        D0[0][j] = D0[0][j-1]*D3[j];
    }
    double dij;
    for (i = 1; i <= r; i++) {
        for (j = 1; j <= c; j++) {
            dij = D[i][j];
            D[i][j] = (D[i - 1][j - 1] + D[i - 1][j] + D[i][j - 1])* dij ;
            //if (i == j) {
            //    D0[i][j] = D0[i - 1][j - 1] * dij
            //            +  D0[i-1][j]* (D3[i]) + D0[i][j-1]*D3[j];
            //} else
            //    D0[i][j] = D0[i-1][j]* (D3[i]) + D0[i][j-1]*D3[j];
            D0[i][j] = D0[i - 1][j - 1] *(D3[i]*D3[j])/2 + D0[i-1][j]* (D3[i]) + D0[i][j-1]*D3[j];
        }
    }

    dist = D[r][c] + D0[r][c] ;
    return D;
}


std::vector< std::vector<double> > kdtw_fbmat(std::vector< std::vector<double> >  ts1, std::vector< std::vector<double> >  ts2, double nu)
{
    std::vector< std::vector<double> >  matk, matkr;
    matk = (kdtw_mat(ts1, ts2, nu));
    std::vector< std::vector<double> > ts1r = reverseArray(ts1);
    std::vector< std::vector<double> > ts2r = reverseArray(ts2);
    matkr = (kdtw_mat(ts1r, ts2r, nu));
    int l1=matk.size();
    int l2=matk[0].size();
    for(int i=0; i<l1; i++){
        for(int j=0; j<l2; j++){
            matk[i][j]= matk[i][j] * matkr[l1-i-1][l2-j-1];
        }
    }
    return matk;
}

std::vector< std::vector<double> > eKats2::interpolate(std::vector< std::vector<double> > ts){
int l=ts.size(), dim=ts[0].size();
LinearInterpolator LI;
std::vector< std::vector<double> > tsout(l, std::vector<double>(dim-1, 0.0));
std::vector<double> P=std::vector<double>(dim-1); 
for(int i=0; i<l; i++){
  for(int j=1; j<dim; j++)
     P[j-1]=ts[i][j];
  LI.addDataPoint(ts[i][0], P);
}
for(int i=0; i<l; i++)
  for(int j=0; j<dim-1; j++){
    tsout[i][j]=LI.interpolate(i,j);
  }
return tsout;
}

double getInertia(std::vector< std::vector<double> > C, std::vector<std::vector<Point> > tab){
double inertia=0.0;
for (unsigned int i=0; i< tab[0].size(); i++) {
  for (unsigned int j=0; j< C.size(); j++) {
    inertia+=sq_error(C[j],tab[j][i]);
  } 
}
return inertia;
}

Point barycenter(std::vector<Point> tab) {
    if (tab.size() < 1) {
        throw std::runtime_error("empty point in barycenter");
    }
    unsigned int dim = tab[0].size();
    Point sum(dim,0.0);
    for (unsigned int i=0; i< tab.size(); i++) {
        for(unsigned int j=0; j< dim; j++)
            sum[j] += (tab[i][j]);
    }
    for (unsigned int i=0; i< dim; i++)
        sum[i]=sum[i]/tab.size();
    return sum;
}

double barycenterT(std::vector<double> tab) {
    if (tab.size() < 1) {
        throw std::runtime_error("empty point in barycenterT");
    }
    double sum=0.0;
    for (unsigned int i=0; i< tab.size(); i++) {
        sum += (tab[i]);
    }
    sum=sum/tab.size();
    return sum;
}

Point standard_deviation(std::vector<Point> tab, Point B) {
    if (tab.size() < 1) {
        throw std::runtime_error("empty point in standard_deviation");
    }
    unsigned int dim = tab[0].size();
    Point sum(dim,0.0);
    for (unsigned int i=0; i< tab.size(); i++) {
        for(unsigned int j=0; j< dim; j++)
            sum[j] += (B[j]-tab[i][j])*(B[j]-tab[i][j]);
    }
    for (unsigned int i=0; i< dim; i++)
        sum[i]=sqrt(sum[i]/tab.size());
    return sum;
}

double standard_deviationT(std::vector<double> tab, double B) {
    if (tab.size() < 1) {
        throw std::runtime_error("empty point in standard_deviationT");
    }
    double sum=0.0;
    for (unsigned int i=0; i< tab.size(); i++) {
        sum += (tab[i]-B)*(tab[i]-B);
    }
    sum=sqrt(sum/tab.size());
    return sum;
}


std::vector< std::vector<double> > eKats2::iTEKA(std::vector< std::vector<double> > C, std::vector<std::vector<std::vector<double> > > dataset, double nu)
{
    double _dbl_min=1e-300;
    std::vector< std::vector<double> >  ts2;
    std::vector<std::vector<Point> > tupleAssociation;
    std::vector<std::vector<double> > timeStampsAssociation;
    tupleAssociation.resize(C.size());
    timeStampsAssociation.resize(C.size());
    const int l1 = C.size();
    const int dim=C[0].size();

    for(unsigned int n=0; n<dataset.size(); n++){
        ts2=dataset[n];
        std::vector< std::vector<double> > mat=(kdtw_fbmat(C,ts2,nu));
        
        const int l2 = ts2.size();

        if (l1 == 0 || l2 == 0)
            throw std::runtime_error("time series cannot be empty");

        double z;
        std::vector< std::vector<double> > ts22(l1, std::vector<double>(dim, 0.0));
        std::vector<double> normzj(l1, _dbl_min);
   
        for(int i=0; i<l1; i++){
            double t=0.0;
            for(int j=0; j<l2; j++){
                z=mat[i+1][j+1]+_dbl_min;
                t+=z*(double)j;
                for(int k=0; k<dim; k++){
                    ts22[i][k]+=(ts2[j][k])*z;
                }
                normzj[i]+=z;
            }
            for(int k=0; k<dim; k++){
                ts22[i][k]=ts22[i][k]/normzj[i];
            }
            t=t/normzj[i];
            tupleAssociation[i].push_back(ts22[i]);
            timeStampsAssociation[i].push_back(t);
        }
    }

    std::vector< std::vector<double> > ts11(l1, std::vector<double>(dim+1, 0.0));
    for(int i=0; i<l1; i++){
        Point p=barycenter(tupleAssociation[i]);
            double t=barycenterT(timeStampsAssociation[i]);
            for (int j=1; j<dim+1; j++)
        ts11[i][j]=p[j-1];
        ts11[i][0]=t;
        tupleAssociation[i].clear();
        timeStampsAssociation[i].clear();
    }
      std::sort(ts11.begin(), ts11.end());
    tupleAssociation.clear();
    timeStampsAssociation.clear();
    return ts11;
}

std::vector<std::vector< std::vector<double> >> eKats2::iTEKA_stdev(std::vector< std::vector<double> > C, std::vector<std::vector<std::vector<double> > > dataset, double nu)
{
    double _dbl_min=1e-300;
    std::vector< std::vector<double> >  ts2;
    std::vector<std::vector<Point> > tupleAssociation;
    std::vector<std::vector<double> > timeStampsAssociation;
    tupleAssociation.resize(C.size());
    timeStampsAssociation.resize(C.size());
    const int l1 = C.size();
    const int dim=C[0].size();
    std::vector<std::vector< std::vector<double> >> rout={};
    for(unsigned int n=0; n<dataset.size(); n++){
        ts2=dataset[n];
        std::vector< std::vector<double> > mat=(kdtw_fbmat(C,ts2,nu));
        
        const int l2 = ts2.size();

        if (l1 == 0 || l2 == 0)
            throw std::runtime_error("time series cannot be empty");

        double z;
        std::vector< std::vector<double> > ts22(l1, std::vector<double>(dim, 0.0));
        std::vector<double> normzj(l1, _dbl_min);
   
        for(int i=0; i<l1; i++){
            double t=0.0;
            for(int j=0; j<l2; j++){
                z=mat[i+1][j+1]+_dbl_min;
                t+=z*(double)j;
                for(int k=0; k<dim; k++){
                    ts22[i][k]+=(ts2[j][k])*z;
                }
                normzj[i]+=z;
            }
            for(int k=0; k<dim; k++){
                ts22[i][k]=ts22[i][k]/normzj[i];
            }
            t=t/normzj[i];
            tupleAssociation[i].push_back(ts22[i]);
            timeStampsAssociation[i].push_back(t);
        }
    }
    std::vector< std::vector<double> > out(l1, std::vector<double>(2*dim+2, 0.0));
    std::vector< std::vector<double> > errv(l1, std::vector<double>(1,0.0));
    std::vector< std::vector<double> > inertia(1, std::vector<double>(1, 0.0));
    inertia[0][0]=getInertia(C,tupleAssociation);
    for(int i=0; i<l1; i++){
        Point p=barycenter(tupleAssociation[i]);
        errv[i][0]=sq_error(p,C[i]);
        Point s=standard_deviation(tupleAssociation[i],p);
        double t=barycenterT(timeStampsAssociation[i]);
        double s_t=standard_deviationT(timeStampsAssociation[i],t);
            for (int j=1; j<dim+1; j++)
                out[i][j]=p[j-1];
            for (int j=dim+1; j<2*dim+1; j++)
                out[i][j]=s[j- dim - 1];
        out[i][0]=t;
         out[i][2*dim+1]=s_t;
        tupleAssociation[i].clear();
    }
    std::sort(out.begin(), out.end());
    tupleAssociation.clear();
    rout.push_back(out);
    rout.push_back(timeStampsAssociation);
    rout.push_back(errv);
    rout.push_back(inertia);
    return rout;
}






