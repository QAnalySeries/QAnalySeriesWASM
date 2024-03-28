#include <sstream>

#include "spectral.h"
#include "ssaAS.h"


/*
void	MyMath::common_regular_scale::common_regular_scale( const val1Darray<double_vector1*,1>& in_x, bool on_intersection )
{	
		//	first find common limits and a common step (min of the means)
	double	maxi_n = -1e222;
	double	mini_1 = 1e222;
	double	maxi_1 = -1e222;
	double	mini_n = 1e222;
	double	minstep = 1e222;
	for (size_t i=1; i<=n; i++)
	{	double	xn = (*in_x(i))(in_x(i)->length());
		if (xn > maxi_n)	maxi_n = xn;
		if (xn < mini_n)	mini_n = xn;
		double	x1 = (*in_x(i))(1);
		if (x1 < mini_1)	mini_1 = x1;
		if (x1 > maxi_1)	maxi_1 = x1;
		double	stp = (xn-x1)/(in_x(i)->length()-1);
		if (stp < minstep)	minstep = stp;
	}
	empty_intersect = (maxi_n <= mini_1);
	failed = (empty_intersect && on_intersection);
	if (failed)	return;
	
	double	maxi, mini;
	if (on_intersection)	{	maxi = mini_n;	mini = maxi_1;	}
	else					{	maxi = maxi_n;	mini = mini_1;	}
	
	Boolean reg = in_x(1)->check_regular();
	for (size_t i=2; i<=n && reg; i++)		reg = in_x(1)->check_same_regular_step_as( *in_x(i) );
	
	if ( (maxi_n != mini_n) || (maxi_1 != mini_1) || !reg )		//	needs a new scale
	{		//	rounding of common step
		double	step = MyMath::Arondi( minstep, MyMath::round_nearest, -1 );
		double	dstep	= step * 1e-8;
		size_t	n = 1 + MyMath::long_round((maxi - mini)/step);
		double	newmax = mini + (n-1)*step;
		if (newmax - maxi < dstep && newmax > maxi)			newmax += step;
		else if (maxi - newmax < dstep && newmax < maxi)	newmax -= step;
		size_t	nn = 1 + MyMath::long_round((newmax - mini)/step);
		
		resize( nn );
		init_regular( mini, newmax );
		
	
	}
	else			//	just copy the first one
	{	//x = new MyMath::double_regular_array( (*in_x(1))(1), (*in_x(1))(in_x(1)->length()), in_x(1)->length() );
		//return true;
		
		resize( in_x(1)->length() );
		init_regular( (*in_x(1))(1), (*in_x(1))(in_x(1)->length()) );
	}
};*/

std::string& MyMath::BTSpectrum_array::comment()
{	comm = "";
	comm += "BTukey spectrum using a ";
	comm += window->Name();
	comm += " window\r";
	if (lag >= 2)			//	same code in "CSpectralDialog.h"
	{	ostringstream s;
		s << "The bandwidth is " << bw;
		if (spectra)		//	ie, if cross-analysis
			s << "\rNon-zero coherence is higher than " << window->NonZeroCoherence();
		s << "\rThe error estimation on the power spectrum is " << window->LowError()
            << " < ∆Power / Power < " << window->HighError();
		comm += s.str().c_str();
	}
	return	comm;
};

	
void MyMath::DPSS_array::DoComputeVap()
{	double_power2_array	fft_dpss( 2*nndp );
	double	delta_fre = 1 / (2*nndp);
	size_t ndpss = dim1();
	size_t npo 	 = dim2();
	double_matrix1&	dpss = *this;
	
	for (size_t j=1; j<=ndpss; j++)
	{	size_t	npomax = std::min( npo, 2*nndp );
		for (size_t i=1; i<=npomax; i++)		fft_dpss[i] = dpss(j,i);
		for (size_t i=(npo+1); i<=2*nndp; i++)	fft_dpss[i] = 0;
		realft( fft_dpss );
		double	integ_W = sqr( fft_dpss[1] );
		double	freq = 0;
		size_t	i;
		for (i=2; (i<=nndp)&&(freq<=wwidth); i++)
		{	freq = (i-1) * delta_fre;
			integ_W += ( sqr( fft_dpss[2*i-1] ) + sqr( fft_dpss[2*i] ) );
		};
		size_t	ii = i;
		double	integ_res = sqr( fft_dpss[2] );
		for (size_t i=ii; i<=(nndp-1); i++)
			integ_res += ( sqr( fft_dpss[2*i-1] ) + sqr( fft_dpss[2*i] ) );
		double	integ = 0;
		for (size_t i=1; i<=nndp; i++)
			integ += ( sqr( fft_dpss[2*i-1] ) + sqr( fft_dpss[2*i] ) );
		vap[j] 		= integ_W 	/ integ;
		un_vap[j] 	= integ_res / integ;
	};
};


void	MyMath::DPSS_array::DoCompute()
{	size_t ndpss = dim1();
	size_t npo 	 = dim2();
	double_matrix1&	dpss = *this;

	double_vector1	diago(npo);
	double_vector1	offdiago(npo);
	
	double	cos2pW = cos(2 * Pi * wwidth);
	for (size_t i=0; i<=(npo-1); i++)	diago[i+1] = sqr((npo-1)*0.5-i)*cos2pW;
	for (size_t i=1; i<=(npo-1); i++)	offdiago[i+1] = i * (npo-i) * 0.5;

	tqli( diago, offdiago );
	sort( diago );
	for (size_t j=1; j<=ndpss; j++)
	{	dpss(j,1) = 1;
		dpss(j,2) = dpss(1,j) * (cos2pW * sqr((npo-1) * 0.5) - diago[npo-j+1]);
		dpss(j,2) = - dpss(j,2) / ((npo-1) * 0.5);
		for (size_t i=1; i<=((npo-1)/2)-1; i++)
		{	dpss(j,i+2) = dpss(j,i) * i * (npo - i) * 0.5;
			dpss(j,i+2) += (cos2pW * sqr((npo-1)*0.5-i) - diago[npo-j+1]) * dpss(j,i+1);
			dpss(j,i+2) /= -((i + 1) * (npo - 1 - i) * 0.5);
		};
		div_t	r = div( j+1, 2 );
		int		sj = 1 - 2*r.rem;
		for (size_t i=((npo-1)/2)+1; i<=(npo-1); i++)
			dpss(j,i+1) = sj*dpss(j,npo-i);
		double	sumsq = 0;
		for (size_t i=1; i<=npo; i++)	sumsq += sqr(dpss(j,i));
		for (size_t i=1; i<=npo; i++)	dpss(j,i) /= sqrt(sumsq);
	};
};



MyMath::DPSS_array	MyMath::MTMSpectrum_array::dpss;

void	MyMath::MTMSpectrum_array::freqRange( double delta_f, double delta_t, double f_min, double f_max ){
    fRange_nn = 128;
	double	dtdf = delta_f * delta_t;
    do {
        fRange_nn *= 2;
    } while ((fRange_nn < 1/(2*dtdf)) || (fRange_nn < ((npo+1)/2)));
	
	fRange_delta_fre = 1/(delta_t * 2 * fRange_nn);
    fRange_fcount = 1;
	fRange_nf0 = 1;
    for (size_t i=1; i<=fRange_nn+1; i++) {
        double	freq = (i-1) * fRange_delta_fre;
        if ((freq>=f_min) && (freq<=f_max))	{
            fRange_fcount++;
        } else if ( freq < f_min ) {
            fRange_nf0 = i + 1;
        }
	};
	fRange_fcount--;

	x_scale.resize_but_keep( fRange_fcount );
	size_t	k = 1;
	for (size_t i=1; i<=(fRange_nn+1); i++)
	{	double	freq = (i-1) * fRange_delta_fre;
		if ( (freq >= f_min) && (freq <= f_max) )
			x_scale[k++] = freq;
	};

	fftsig.resize2D( dpss.dim1(), 2*fRange_fcount );
	varlog.resize_but_keep( fRange_fcount );
};



void MyMath::MTMSpectrum_array::eigenft( const double_vector1& y )
{
	double_power2_array	data( 2*fRange_nn );
	
    for (size_t j=1; j<=ndpss; j++){
        for (size_t i=1; i<=npo; i++){
            data[i] = dpss(j,i) * y[i];		// Tapering
        }
        for (size_t i=(npo+1); i<=2*fRange_nn; i++)	{
            data[i] = 0;
        }
		
        realft( data );
        for (size_t i=fRange_nf0; i<=(fRange_nf0+fRange_fcount-1); i++)
		{	int		ii = i-fRange_nf0+1;
			if (i == 1)
			{	fftsig(j,2*ii-1) 	= data[1];
				fftsig(j,2*ii) 	= 0;
			}
			else if ( i == fRange_nn+1 )
			{	fftsig(j,2*ii-1) 	= data[2];
				fftsig(j,2*ii) 	= 0;
			}
			else
			{	fftsig(j,2*ii-1) 	= data[2*i-1];
				fftsig(j,2*ii) 	= data[2*i];
			};
		};
	};
};

void MyMath::MTMSpectrum_array::adapt_MTM( double_vector1& power )
{
	double_vector1	d(ndpss);
	double_vector1	S_i(ndpss);
//				*se = 0;
	for (size_t i=1; i<=fRange_fcount; i++)
	{	double	S_part = 0.5 * (sqr(fftsig(1,2 * i - 1)) + sqr(fftsig(1,2 * i)))
					   + 0.5 * (sqr(fftsig(2,2 * i - 1)) + sqr(fftsig(2,2 * i)));
        for (size_t	j=1; j<=(ndpss-1); j++)
		{	for (size_t k=1; k<=ndpss; k++)
				d[k] = dpss.vap[k] / sqr(dpss.vap[k] * S_part + varian * (dpss.un_vap[k]));
			double	denom_i = 0;
			for (size_t k=1; k<=ndpss; k++)		denom_i += d[k];
			double	numer_i = 0;
			for (size_t k=1; k<=ndpss; k++)
			{	double	xx = d[k];
				xx *= (sqr(fftsig(k,2 * i - 1)) + sqr(fftsig(k,2 * i)));
				numer_i += xx;
			};
			S_part = numer_i / denom_i;
		};
		power[i] = S_part;
		
// Jackknifing   
		double	S_all = 0;
		for (size_t j=1; j<=ndpss; j++)
		{	S_i[j] = 0.0;
			for (size_t k=1; k<=ndpss; k++)
				if (k!=j)
					S_i[j] += ( sqr(fftsig(k,2*i-1)) + sqr(fftsig(k,2*i)) );
			S_i[j] = log( S_i[j]/(ndpss-1.0) );
			S_all += S_i[j];
		};
		S_all /= ndpss;
		varlog[i] = 0;
		for (size_t j=1; j<=ndpss; j++)		varlog[i] += sqr(S_i[j] - S_all);
		varlog[i] *= (ndpss-1.0)/ndpss;
	};
};

void MyMath::MTMSpectrum_array::Harmonic( double_vector1*& amp, double_vector1*& conf )
{	
	amp  = new double_vector1( fRange_fcount );
	conf = new double_vector1( fRange_fcount );
	
	double_vector1	U_k0( ndpss );
	for (size_t i=1; i<=ndpss; i++)
	{	U_k0[i] = 0;
		for (size_t j=1; j<=npo; j++)	U_k0[i] += dpss(i,j);
	};

	double	sum_sq_dpwf = 0.0;
	for (size_t i=1; i<=ndpss; i++)		sum_sq_dpwf += sqr(U_k0[i]);
	
	double_vector1	muhat( 2*fRange_fcount );
	muhat.set_all_to(0.0);
	for (size_t i=1; i<=2*fRange_fcount; i++)
	{	for (size_t j=1; j<=ndpss; j++)		muhat[i] += U_k0[j] * fftsig(j,i);
		muhat[i] /= sum_sq_dpwf;
	};
	
	double_vector1	f_test( fRange_fcount );
	f_test.set_all_to(0.0);
	for (size_t i=1; i<=fRange_fcount; i++)
	{	for (size_t j=1; j<=ndpss; j++)
		{	double	diff1 = fftsig(j,2*i-1) - muhat[2*i-1] * U_k0[j];
			double	diff2 = fftsig(j,2*i)   - muhat[2*i] * U_k0[j];
			f_test[i] += sqr(diff1) + sqr(diff2);
		};
		f_test[i] = (ndpss - 1) * (sum_sq_dpwf / f_test[i]) * ( sqr(muhat[2*i-1]) + sqr(muhat[2*i]) );
	};

//	saving...
	double	nu2 = 2 * ndpss - 2;
	double	nu1 = 2;
	for (size_t i=1; i<=fRange_fcount; i++)
	{	(*amp)[i]  = 2.0 * sqrt( sqr(muhat[2 * i - 1]) + sqr(muhat[2 * i]) );
		(*conf)[i] = f_trans( f_test[i], nu1, nu2 );
	};
};


void	MyMath::BTSpectrum_array::computeFFT( const double_vector1& yy, a_double_power2plus1_array& y )
{	MyMath::double_filter_vector autoc( lag );
	MyMath::autocorrel_fft( yy, autoc, (double)(yy.size()) );			
	y.copy( autoc );
	window->SetSize( autoc.size() );
	for (size_t i=1; i<=autoc.size(); i++)	y(i) *= window->WindowValueAt(i);
	y.DoCosFFT();
	for (size_t i=1; i<=y.size(); i++)
		y[i] *= 4*xstep;
};


void	MyMath::BTSpectrum_array::cross_computeFFT(
    const double_vector1& yy,
    const double_vector1& zz,
    const double_vector1& sy,
    const double_vector1& sz,
    a_double_power2plus1_array& coh,
    a_double_power2plus1_array& ph )
{	
	MyMath::double_filter_vector crosscor( lag, lag );
	MyMath::correl_fft( yy, zz, crosscor, (double)(yy.size()), double_vector1::remove_mean );
	window->SetSize( lag+1 );
	double_power2plus1_array		sym_sum( coh.length() );
	double_power2_array				sym_dif( coh.length()-1 );

	sym_sum(1) = 2*crosscor.value(0) * window->WindowValueAt(1);
	sym_dif(1) = 0;
	for (size_t i=1; i<=lag; i++)
	{	double	wi = window->WindowValueAt(i+1);
		sym_sum(i+1) = (crosscor.value(i)+crosscor.value(-i))*wi;
		sym_dif(i+1) = (crosscor.value(i)-crosscor.value(-i))*wi;
	}

	sym_sum.DoCosFFT();					//	=> co-spectrum
	sym_dif.DoSinFFT();					//	=> quadrature-spectrum
	
	for (size_t i=1; i<=coh.length(); i++)
	{	double	a = sym_sum[i];
		double	b = (i == coh.length() ? sym_dif[i-1] : sym_dif[i]);
		ph[i]  = atan2(-b, a);
		coh[i] = sqrt( (a * a + b * b) * sqr(2*xstep) / fabs(sy[i] * sz[i]) );
	};
};
