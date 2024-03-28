#ifndef MYNR_H
#define MYNR_H

#include <QObject>
#include <QVector>

namespace mynr {

/*******************************************************
My Simple Numerical Recipes
  -- MYNR --
     Copyright 08.08.2018 by Sergey Kotov
********************************************************/
// /////////////////////////////////////////////////////

/*******************************************************
Function: corr (Pearson's correlation coefficient)
INPUT:  QVector x, y;
OUTPUT: double r
********************************************************/
double corr(QVector<double> x, QVector<double> y);

/*******************************************************
Function: mean value  of a QVector
INPUT:  QVector v;
OUTPUT: double mean
********************************************************/
void meanv(QVector<double> v, double & mean);

/*******************************************************
Function: mean value and variance of a QVector
INPUT:  v;
OUTPUT: mean, var
********************************************************/
void meanandvar(QVector<double> v, double & mean, double & var);

/*******************************************************
Function: linear regression by least squares method
INPUT:  QVector<double> x, y - training set ;
OUTPUT: double a, b - coefficients of y = a + b*x
********************************************************/
void lr(QVector<double> x, QVector<double> y, double & a, double  & b);

/*******************************************************
Function: linear piecewise interpolation
INPUT:  QVector<double> x, y - data set ;
        double x1 - point to interpolate
OUTPUT: double y1 - interpolated value
********************************************************/
void linterp(QVector<double> x, QVector<double> y, double & x1, double & y1);

/*******************************************************
Function: linear piecewise interpolation for vector input-output
INPUT:  QVector<double> x, y - data set ;
        QVector<double> x1 - vector of points to be interpolated (must be increasing)
OUTPUT: QVector<double> y1 - vector of interpolated value
********************************************************/
void linterpv(QVector<double> x, QVector<double> y, QVector<double> & x1, QVector<double> & y1);

/*******************************************************
Function: periodogram
INPUT:  QVector<double> x - data set ;
OUTPUT: QVector<double> psd, N/2 bins
********************************************************/
QVector<double> periodogram(QVector<double> x);

/*******************************************************
Function: taner filter
INPUT:  QVector<double> xv - data set ; dt - sampling interval; fl, fh - cutting frequencies, c - rolling factor
OUTPUT: std::vector<double> filt - filter coefficients in N/2 bins

Modified from Astrochron R package,
modified from L. Hinnov FORTRAN code TANER.FOR
also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
********************************************************/
std::vector<double> taner(QVector<double> xv, double dt, double fl, double fh, long long c);

/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: filter
INPUT:  QVector<double> xv - data set ; filt - filter coefficients for positive and negative frequencies
OUTPUT: QVector<double> yv - filtered data;
Modified from Astrochron R package,
modified from L. Hinnov FORTRAN code TANER.FOR
also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
********************************************************/
QVector<double> filter(QVector<double> xv, std::vector<double> filt);

/*******************************************************
Function: ssa (singular spectrum analysis)
INPUT:  QVector<double> y_d - data set ; delay - size of delay-time window
OUTPUT: QVector<double> eigval - eigen values
        QVector<QVector<double>> eigvec - matrix of eigen vectors in rows
        QVector<QVector<double>> rcs - matrix of reconstructed components in rows
The algorithm is based on reworked and modified SciLab code from the 'SSA beginner's guide'
(Claessen and Groth followed Ghil et al, 2002) with improved normalisation of RCs at the beginning.
Uses Eigen library for matrix algebra (http://eigen.tuxfamily.org/)
********************************************************/
void ssa(QVector<double> y_d, int delay, QVector<double> & eigval, QVector<QVector<double>> & eigvec, QVector<QVector<double>> & rcs);

/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: envelope
INPUT:  QVector<double> xv - data set ;
OUTPUT: QVector<double> yv - filtered data;
Based in Matlab hilbert function algorithm:
Marple, S. L. “Computing the Discrete-Time Analytic Signal via FFT.”
IEEE® Transactions on Signal Processing. Vol. 47, 1999, pp. 2600–2603.
********************************************************/
QVector<double> envelope(QVector<double> xv);



}


#endif // MYNR_H
