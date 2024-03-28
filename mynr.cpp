#include "mynr.h"
#include <math.h>
#include <QDebug>
#include <QVector>
#include <unsupported/Eigen/FFT> // http://eigen.tuxfamily.org/
#include <Eigen/Dense>

using namespace Eigen;

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
double corr(QVector<double> x, QVector<double> y) {
    double r, mx, my, vx, vy, sum = 0;
    int n;
    n = x.length();
    meanandvar(x, mx, vx);
    meanandvar(y, my, vy);
    for (int i=0; i<n; ++i) sum = sum + x[i]*y[i];
    r = (sum - n*mx*my)/((n-1)*sqrt(vx)*sqrt(vy));
    return r;
}


/*******************************************************
Function: mean value  of a QVector
INPUT:  QVector v;
OUTPUT: double mean
********************************************************/
void meanv(QVector<double> v, double & mean) {
    int n;
    double sum;
    n = v.length(); sum=0.0;
    for (int i=0; i<n; i++) {sum+=v[i];}; mean = sum/n;
}

/*******************************************************
Function: mean value and variance of a QVector
INPUT:  QVector v;
OUTPUT: double mean, var
********************************************************/
void meanandvar(QVector<double> v, double & mean, double & var) {
    int n;
    double sum;
    n = v.length(); sum=0.0;
    meanv(v, mean);
    for (int i=0; i<n; i++) {sum+=(v[i]-mean)*(v[i]-mean);}; var = sum/(n-1);
}

/*******************************************************
Function: linear regression by least squares method
INPUT:  QVector<double> x, y - training set ;
OUTPUT: double a, b - coefficients of y = a + b*x
********************************************************/
void lr(QVector<double> x, QVector<double> y, double & a, double & b) {
    int n;
    double sxx, sxy, mx, my;
    n = x.length();
    meanandvar(x, mx, sxx); sxx =sxx*n;
    meanv(y, my);

    sxy=0.0; for (int i=0; i<n; i++) {sxy+=(x[i]-mx)*(y[i]-my);}
    b = sxy/sxx;
    a = my - b*mx;

}

/*******************************************************
Function: linear piecewise interpolation
INPUT:  QVector<double> x, y - data set ;
        double x1 - point to interpolate
OUTPUT: double y1 - interpolated value
********************************************************/
void linterp(QVector<double> x, QVector<double> y, double & x1, double & y1) {
    int n, i;
    n = x.length(); i = 0;
    if (x1<x[0])            {y1 = y[0] + (y[1] - y[0]) * (x1 - x[0])/(x[1] - x[0]); }
    else  if  (x1>=x[n-1])   {y1 = y[n-2] + (y[n-1] - y[n-2]) * (x1 - x[n-2])/(x[n-1] - x[n-2]);}
    else                    {do { i++ ;} while (x1 > x[i]);
        y1 = y[i-1] + (y[i] - y[i-1]) * (x1 - x[i-1])/(x[i] - x[i-1]);
    }
}

/*******************************************************
Function: linear piecewise interpolation for vector input-output
INPUT:  QVector<double> x, y - data set ;
        QVector<double> x1 - vector of points to be interpolated (must be increasing)
OUTPUT: QVector<double> y1 - vector of interpolated value
********************************************************/
void linterpv(QVector<double> x, QVector<double> y, QVector<double> & x1, QVector<double> & y1) {
    int n, i;
    n = x.length(); i = 0;
    y1.clear();
    for (int j=0; j<x1.length();++j) {
        if (x1[j]<x[0])            {y1.append(y[0] + (y[1] - y[0]) * (x1[j] - x[0])/(x[1] - x[0])); }
        else  if  (x1[j]>=x[n-1])   {y1.append(y[n-2] + (y[n-1] - y[n-2]) * (x1[j] - x[n-2])/(x[n-1] - x[n-2]));}
        else                    {while (x1[j] >= x[i]) {++i;}
            y1.append( y[i-1] + (y[i] - y[i-1]) * (x1[j] - x[i-1])/(x[i] - x[i-1]));
        }
    }
}


/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: periodogram
INPUT:  QVector<double> x - data set ;
OUTPUT: QVector<double> psd, N/2 bins
********************************************************/
QVector<double> periodogram(QVector<double> xv) {
    int n = xv.length();
    QVector<double> psd;

    if (n % 2 != 0) {n++; xv.append(0);} // check if data length is even and improve if not

    Eigen::FFT<double> fft;
    //std::vector<double> x = xv.toStdVector();
    std::vector<double> x = std::vector<double>(xv.begin(), xv.end());
    std::vector<std::complex<double> > xf;
    fft.fwd(xf, x);

    for (int i=0; i<x.size()/2; i++) {psd.append(sqrt((pow(xf[i].real()/n,2)+pow(xf[i].imag()/n,2))));}

    psd.first() = psd.first()/2.0;

    psd.last() = psd.last()/2.0;

    return psd;

}


/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: taner filter
INPUT:  QVector<double> xv - data set ; dt - sampling interval; fl, fh - cutting frequencies, c - rolling factor
OUTPUT: QVector<double> yv - filtered data; filt - filter coefficients in N/2 bins

modified from L. Hinnov Matlab code
also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
********************************************************/
std::vector<double> taner(QVector<double> xv, double dt, double fl, double fh, long long c) {

    const double PI  =3.141592653589793238463;
    long double twopi, w, dw, bw, wl, wc, wh, amp, darg, arg, arg1, arg2, twod, aw;
    int npts, ncut, pad = 0;
    std::vector<double> filt;
    twopi = 2.0*PI;
    //filt.clear();
    npts =xv.length();
    if (npts % 2 != 0) {npts++; xv.append(0); pad = 1;} // check if data length is even and improve if not
    dw = twopi/(npts*dt);
    if (fl<=0.0) fl=-1.0*fh;
    wl = twopi*fl;
    wc = twopi*(fh+fl)/2.0;
    wh = twopi*fh;
    bw = wh - wl;
    amp = 1.0/sqrt(2.0);
    arg1 = 1.0 - (c*log(10.0)) / (20.0*log(amp));
    arg1 = log(arg1);
    arg2 = pow((bw+2.0)/bw,2);
    arg2 = log(arg2);
    twod = 2.0*arg1/arg2;
    ncut = floor(npts/2.0+1);
    //Filter for positive frequencies
    for (int n=0; n<ncut; n++) {
        w = n*dw;
        arg = pow(2.0*std::abs(w-wc)/bw,twod);
        darg =-1.0*arg;
        filt.push_back(amp * exp(darg));
    }
    //Filter for negative frequencies
    ncut = ncut+1;
    for (int n=ncut; n<=(npts); n++) {
        w = (n-npts)*dw;
        aw = std::abs(w);
        arg = pow(2.0*std::abs(aw-wc)/bw,twod);
        darg =-1.0*arg;
        filt.push_back(amp * exp(darg));
    }
    double max = *std::max_element(filt.begin(), filt.end());

    for (int i=0; i<=filt.size(); i++) filt[i] =filt[i]/max;
    return filt;
}

/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: filter
INPUT:  QVector<double> xv - data set ; filt - filter coefficients for positive and negative frequencies
OUTPUT: QVector<double> yv - filtered data;
Modified from Astrochron R package,
modified from L. Hinnov FORTRAN code TANER.FOR
also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
********************************************************/
QVector<double> filter(QVector<double> xv, std::vector<double> filt){
    int npts =xv.length(), pad = 0;
    if (npts % 2 != 0) {npts++; xv.append(0); pad =1;} // check if data length is even and improve if not
    Eigen::FFT<double> fft;
    std::vector<double> x = std::vector<double>(xv.begin(), xv.end());
    std::vector<double> y;
    std::vector<std::complex<double> > xf;

    fft.fwd(xf, x);
    //Apply the filter over all frequencies
    for (int n=0; n<=(npts-1); n++) {
        //qDebug()<< xf[n].real()<<xf[n].imag()<< filt[n];
        xf[n].real(xf[n].real() * filt[n]);
        xf[n].imag(xf[n].imag() * filt[n]);

    }
    //Inverse FFT
    fft.inv(y, xf);
    QVector<double> yv;
    //yv = QVector<double>::fromStdVector(y);
    yv = QVector<double>(y.begin(), y.end());
    if (pad == 1) yv.removeLast();
    return yv;

}

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
void ssa(QVector<double> y_d, int delay, QVector<double> & eigval, QVector<QVector<double>> & eigvec, QVector<QVector<double>> & rcs){
    int n_d = y_d.length();

    VectorXd  y(n_d), eval;
    MatrixXd my = MatrixXd::Zero(n_d, delay), covM , evec, pcomp, rc = MatrixXd::Zero(n_d, delay), z = MatrixXd::Zero(n_d, delay);
    for (int i=0; i<n_d;++i) {y(i) = y_d[i];} // init VectorXd data

    for (int i=0; i<delay; i++) my.block(0,i,n_d-i,1) = y.segment(i, n_d-i); // delay time matrix constructor

    covM = my.transpose()*my/n_d; // covariance matrix calculation

    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(covM);
    eval = es.eigenvalues().reverse();
    evec = es.eigenvectors().rowwise().reverse(); // EOFs
    pcomp = my*evec; // principal components

    // Reconstructed Components RC
    for (int m=0; m<delay; ++m) {
        z = MatrixXd::Zero(n_d, delay);
        for (int m2=0; m2<delay; ++m2) z.block(m2,m2,n_d-m2,1) = pcomp.block(0,m,n_d-m2,1);
        rc.col(m) = z * evec.col(m) / delay;
        //Improvement for normalization at the beginning
        for (int i=0; i<delay; ++i) rc(i,m)=rc(i,m)*delay/(i+1);
        //Improvement for normalization at the end
        //for (int i=n_d-1; i>n_d-delay; --i) rc(i,m)=rc(i,m)*delay/(n_d-i);
    }

    // Output to QVectors
    eigval.clear();
    for (int i=0; i< delay; ++i) eigval.append(eval(i));

    QVector<double> tmp;
    eigvec.clear();rcs.clear();
    for (int i=0; i< delay; ++i) {for (int j=0; j< delay; ++j) tmp.append(evec(j,i));
        eigvec.append(tmp);
        tmp.clear();
    }

    for (int i=0; i< delay; ++i) {for (int j=0; j< n_d; ++j) tmp.append(rc(j,i));
        rcs.append(tmp);
        tmp.clear();
    }
}


/*******************************************************
Uses Eigen library for FFT (http://eigen.tuxfamily.org/)
Function: envelope
INPUT:  QVector<double> xv - data set ;
OUTPUT: QVector<double> yv - filtered data;
Based in Matlab hilbert function algorithm:
Marple, S. L. “Computing the Discrete-Time Analytic Signal via FFT.”
IEEE® Transactions on Signal Processing. Vol. 47, 1999, pp. 2600–2603.
********************************************************/
QVector<double> envelope(QVector<double> xv){
    int npts =xv.length(), pad = 0;
    if (npts % 2 != 0) {npts++; xv.append(0); pad =1;} // check if data length is even and improve if not
    Eigen::FFT<double> fft;
    std::vector<double> x = std::vector<double>(xv.begin(), xv.end());//.toStdVector(),
    std::vector<double> y;
    std::vector<std::complex<double> > xf, anx;

    fft.fwd(xf, x);
    for (int i = 1; i<npts/2; ++i) xf[i]*=2;
    for (int i = npts/2+1; i<npts; ++i) xf[i]=0;
    fft.inv(anx, xf);
    for (int i = 0; i<npts; ++i) y.push_back(std::abs(anx[i]));

    QVector<double> yv;
    yv = QVector<double>(y.begin(), y.end());
    if (pad == 1) yv.removeLast();
    return yv;
}
}

