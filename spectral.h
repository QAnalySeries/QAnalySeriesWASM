
#pragma once


#include "myMath.h"





namespace MyMath {




/////////////	SPECTRAL ANALYSIS	////////////



class square_spectralwindowfunc : public spectralwindowfunc {
public:
	square_spectralwindowfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new square_spectralwindowfunc();	}

	virtual double DoValueAt( double )			const {	return 1;	};
	virtual double DoDerivValueAt( double )		const {	return 0;	};
	virtual double DoPrimitiv( double t )		const {	return t;	};
	virtual	double FullIntegral()				const {	return 2.0;	};
	virtual	double Energy()						const {	return 2.0;	};				//	fullIntegral of the square !! 
	
	virtual	char*	Name()		{	return	(char *)"square";	};
};

class Bartlett_spectralwindowfunc : public spectralwindowfunc {
public:
	Bartlett_spectralwindowfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new Bartlett_spectralwindowfunc();	}

	virtual double DoValueAt( double t )		const {	return 1-t;			};
	virtual double DoDerivValueAt( double )		const {	return -1;			};
	virtual double DoPrimitiv( double t )		const {	return t*(1-t/2);	};
	virtual	double FullIntegral()				const {	return 1.0;			};
	virtual	double Energy()						const {	return 2.0/3.0;		};				//	fullIntegral of the square !! 
	
	virtual	char*	Name()		{	return	(char *)"Bartlett";	};
};

class Welch_spectralwindowfunc : public spectralwindowfunc {
public:
	Welch_spectralwindowfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new Welch_spectralwindowfunc();	}

	virtual double DoValueAt( double t )		const {	return 1-t*t;	};
	virtual double DoDerivValueAt( double t )	const {	return -2*t;	};
	virtual double DoPrimitiv( double t )		const {	return t*(1 - t*t/3);	};
	virtual	double FullIntegral()				const {	return 4.0/3.0;		};
	virtual	double Energy()						const {	return 16.0/15.0;	};				//	fullIntegral of the square !! 
	
	virtual	char*	Name()		{	return	(char *)"Welch";	};
};

class Hann_spectralwindowfunc : public spectralwindowfunc {
public:
	Hann_spectralwindowfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new Hann_spectralwindowfunc();	}

	virtual double DoValueAt( double t )		const {	return 0.5 * (1.0 + cos(Pi*t));		};
	virtual double DoDerivValueAt( double t )	const {	return -0.5*Pi*sin( Pi*t );			};
	virtual double DoPrimitiv( double t )		const {	return 0.5 * ( t + sin(Pi*t)/Pi);	};
	virtual	double FullIntegral()				const {	return 1.0;			};
	virtual	double Energy()						const {	return 3.0/4.0;		};				//	fullIntegral of the square !! 
	
	virtual	char*	Name()		{	return	(char *)"Hann";	};
};
class Tukey_spectralwindowfunc : public Hann_spectralwindowfunc {
	virtual	char*	Name()		{	return	(char *)"Tukey";	};
};

class Parzen_spectralwindowfunc : public spectralwindowfunc {
public:
	Parzen_spectralwindowfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new Parzen_spectralwindowfunc();	}

	virtual double DoValueAt( double t )		const {	return (t <= 0.5 ? 1 - 6*t*t*(1-t) : 2*(1-t)*(1-t)*(1-t));	};
	virtual double DoDerivValueAt( double t )	const {	return (t <= 0.5 ? -12*t + 18*t*t  : -6*(1-t)*(1-t));	};
	virtual double DoPrimitiv( double t )		const {	return (t <= 0.5 ? t - 2*t*t*t + 3*t*t*t*t/2  : -4.0/32 + (1 - (1-t)*(1-t)*(1-t)*(1-t))/2);	};
	virtual	double FullIntegral()				const {	return 3.0/8.0;			};
	virtual	double Energy()						const {	return 151.0/280.0;		};				//	fullIntegral of the square !! 
	
	virtual	char*	Name()		{	return	(char *)"Parzen";	};
};




template <class T, int off>
class powerSpectrum_array	{
public:
	regular_array<T,off>	x_scale;					//	the frequency scale
	size_t					data_length;
	double 					xstep;
	bool 					resample;
		
	powerSpectrum_array( size_t ndata, T stp, T fromFreq, T toFreq, T stepFreq )						//	resampled freq. scale
	: x_scale( fromFreq, toFreq, stepFreq ), data_length(ndata), xstep(stp), spectra(NULL)
	  	{	resample = true;	};
		
	powerSpectrum_array( size_t ndata, T stp, size_t n )									//	for the "natural" freq. scales
	: x_scale( 0, 1/(2*stp), n ), data_length(ndata), xstep(stp), spectra(NULL)
	  	{	resample = false;	};
	    
	void	ComputePower( const double_vector1& data, double_vector1*& power )
	{	power = new double_vector1( x_scale.size() );
		do_compute( data, *power );
	};
	void	ComputeCrossPower( const double_vector1& data1, const double_vector1& data2, const double_vector1& sp1, const double_vector1& sp2, double_vector1*& coh, double_vector1*& ph )			//	assuming individaul power spectra known
	{	coh = new double_vector1( x_scale.size() );
		ph  = new double_vector1( x_scale.size() );
		do_cross_compute( data1, data2, sp1, sp2, *coh, *ph );
	};
	
	virtual char*	name() 		{	return	(char *)"";	};
	virtual string&	comment()	{	return	comm;	};

//		bidouille pour garder les spectres individuels "bruts" pour l'analyse croisee
	size_t	counter;
	val1Darray<double_vector*,1>*	spectra;		//	saved for cross analysis
	
	void SetUpCrossAnalysis( size_t n )
		{	spectra = new val1Darray<double_vector*,1> (n);
			counter = 1;
		}
	void SaveSpectra( double_vector* sp )
		{	if (counter <= spectra->length())	(*spectra)(counter++) = sp;
		}
	virtual ~powerSpectrum_array()
		{	if (spectra && resample)	for (size_t i=1; i<=spectra->length(); i++)	delete (*spectra)(i);
			delete	spectra;
		}
//
	
protected:
	virtual	void do_compute( const double_vector1& data, double_vector1& power ) = 0;
	virtual	void do_cross_compute( const double_vector1&, const double_vector1&, const double_vector1&, const double_vector1&, double_vector1&, double_vector1&)	{};
	
	void	interpol( const val1Darray<T,off>& x, const val1Darray<T,off>& y, double_vector1& power )
	{	MyMath::lin_interpolfunc	lin_interp( x, y );
		MyMath::double_extra_midpoints_array	mPoints( x_scale, false );
		for (size_t j=1; j<=x_scale.length(); j++)
			power[j] = lin_interp.MeanBetween( mPoints[j-1], mPoints[j] );
	};
	void	interpol_phase( const val1Darray<T,off>& x, const val1Darray<T,off>& y, double_vector1& power )
	{	MyMath::lin_interpol_y_cyclicfunc	lin_interp( -Pi, Pi, x, y );
		MyMath::double_extra_midpoints_array	mPoints( x_scale, false );
		for (size_t j=1; j<=x_scale.length(); j++)
			power[j] = lin_interp.MeanBetween( mPoints[j-1], mPoints[j] );
	};
	
	string	comm;
};

typedef	powerSpectrum_array<double,1>				double_powerSpectrum;




class Periodogram_array : public double_powerSpectrum	{
	
private:
	spectralwindowfunc*	window;
	void	computeFFT( const double_vector1& yy, a_double_power2plus1_array& y )
		{	double_power2_array		v( yy, *window, MyMath::double_power2_array::remove_mean );
			v.DoFFT();
			size_t	no2 = v.length()/2;
			size_t	n	= 2*no2;
			y[1] 		= MyMath::sqr( v[1]/n );
			y[no2+1]	= MyMath::sqr( v[2]/n );
			for (size_t i=2; i<=no2; i++)
				y[i]	= 2*(MyMath::sqr( v[2*i-1]/n ) + MyMath::sqr( v[2*i]/n ));
		};
	void	specifyWindow( size_t n )
		{	window->SetFlag( MyMath::spectralwindowfunc::centered, MyMath::spectralwindowfunc::half, n );
		};

protected:		
	virtual	void do_compute( const double_vector1& data, double_vector1& power )
		{	if (resample)
			{	specifyWindow( data.size() );
				size_t	no2p1 = nextPower2(data.size())/2 + 1;
		  		double_regular_array		tempX( 0, 1/(2*xstep), no2p1 );
		  		double_power2plus1_array	*tempV = new double_power2plus1_array( no2p1 );
				computeFFT( data, *tempV );
				interpol( tempX, *tempV, power );
				if (spectra)	SaveSpectra(tempV);
				else			delete tempV;
			}
			else
			{	specifyWindow( data.size() );
		  		//double_power2plus1_array	tempV( (valarray0<double>&)power, valarray0<double>::dontCopy );
				temp_double_power2plus1_array	tempV( power );
				computeFFT( data, tempV );
				if (spectra)	SaveSpectra(&power);
			}
		};
		
public:
	~Periodogram_array()
		{	delete window;	};
	Periodogram_array( size_t ndata, spectralwindowfunc* w, double stp )
	: double_powerSpectrum( ndata, stp, nextPower2(ndata)/2 + 1 ), window(w)
	  	{};
		
	Periodogram_array( size_t ndata, spectralwindowfunc* w, double stp, double fromFreq, double toFreq, double stepFreq )		//	idem with resampling
	: double_powerSpectrum( ndata, stp, fromFreq, toFreq, stepFreq ), window(w)
		{};
	
	virtual char*	name()		{	return (char *)"Periodogram";	};
	virtual string&	comment()	{	comm = "";
									comm += "Periodogram (Windowed FFT) using a ";
									comm += window->Name();
									comm += " window";
									return	comm;
								};
};


class BTSpectrum_array : public double_powerSpectrum	{
private:
	spectralwindowfunc*	window;
	size_t 				lag;
	double				bw;
	double 		overFactor;
	double 		level;
	
	void	specifyWindow()
		{	window->SetFlag( MyMath::spectralwindowfunc::Right, MyMath::spectralwindowfunc::one, lag );
		};
	void	finalCheck( double_vector1& power )
		{	if (false)			//	bidouille finale... pas necessaire ????
			{	double	mmin = power.min();
				if (mmin < 0)
				{	double	summ = power.sum();
					double	coef = summ / (summ - x_scale.length() * 2 * mmin);
					for (size_t j=1; j<=x_scale.length(); j++)
						power[j] = ( power[j] - 2 * mmin) * coef;
			}	}
			else
				for (size_t j=1; j<=x_scale.length(); j++)
					if (power[j] < 0)		power[j] = 0;
		};
	void	setstats()
		{	double	m_over_n = lag/(double)(data_length);
			window->SetBTStats( level, m_over_n );
			bw = window->BandWidth( lag )/xstep;
		};
	void	computeFFT( const double_vector1& yy, a_double_power2plus1_array& y );
	void	cross_computeFFT( const double_vector1& yy, const double_vector1& zz, const double_vector1& sy, const double_vector1& sz, a_double_power2plus1_array& coh, a_double_power2plus1_array& ph );
		
	double	upperCoherence( double coh, double dcoh )
		{	double	x = 0.5 * log( (1+coh)/(1-coh) );			//	x = Argth( coh )
			double	y = exp( 2 * (x + dcoh) );
			double	z = (y - 1)/(y + 1);						//	th( x + dcoh )
			return	z;
		}
	double	lowerCoherence( double coh, double dcoh )
		{	double	x = 0.5 * log( (1+coh)/(1-coh) );			//	x = Argth( coh )
			double	y = exp( 2 * (x - dcoh) );
			double	z = std::max( (y - 1)/(y + 1), 0.0 );				//	th( x - dcoh )
			return	z;
		}
	double	deltaPhase( double coh, double arg )
		{	double	r = arg * sqrt( 1/(sqr(coh)+(1e-20)) - 1 );	//	rayon du cercle
			double	delta;
			if ( r<2 && r>0 )	delta = 2 * asin( 0.5 * r );
			else				delta = Pi;
			return	delta;
		}
		
protected:
	void	do_compute( const double_vector1& yy, double_vector1& power )
		{	specifyWindow();
			if (resample)
			{	size_t	nn = nextPower2( ((size_t)(std::max(overFactor,1.0)*(lag + 1))) - 1 ) + 1;
				double_power2plus1_array*	y = new double_power2plus1_array( nn );
				double_regular_array		x( 0, 1/(2*xstep), nn );
				computeFFT( yy, *y );
				interpol( x, *y, power );
				if (spectra)	SaveSpectra(y);
				else			delete y;
			}
			else
			{	//double_power2plus1_array	tempV( (valarray0<double>&)power, valarray0<double>::dontCopy );
				temp_double_power2plus1_array	tempV( power );
				computeFFT( yy, tempV );
				if (spectra)	SaveSpectra(&power);
			}
			finalCheck( power );
			
		};
	void	do_cross_compute( const double_vector1& yy, const double_vector1& zz, const double_vector1& sy, const double_vector1& sz, double_vector1& coh, double_vector1& ph)
		{	specifyWindow();
			if (resample)
			{	size_t	nn = nextPower2( ((size_t)(std::max(overFactor,1.0)*(lag + 1))) - 1 ) + 1;
				double_power2plus1_array	c( nn );
				double_power2plus1_array	p( nn );
				double_regular_array		x( 0, 1/(2*xstep), nn );
				cross_computeFFT( yy, zz, sy, sz, c, p );
				interpol( x, c, coh );
				interpol_phase( x, p, ph  );
			}
			else
			{	//double_power2plus1_array	c_sp( (valarray0<double>&)coh, valarray0<double>::dontCopy );
				//double_power2plus1_array	p_sp( (valarray0<double>&)ph, valarray0<double>::dontCopy );
				temp_double_power2plus1_array	c_sp( coh );
				temp_double_power2plus1_array	p_sp( ph  );
				cross_computeFFT( yy, zz, sy, sz, c_sp, p_sp );
			}
		};
		
public:	
	~BTSpectrum_array()
		{	delete window;	};
	BTSpectrum_array( size_t ndata, size_t n, spectralwindowfunc* w, double overF, double stp, double lev )
	: double_powerSpectrum( ndata, stp, nextPower2( ((size_t)(std::max(overF,1.0)*(n + 1))) - 1 ) + 1 ), window(w), lag(n), overFactor(overF), level(lev)
		{	setstats();
		};
		
	BTSpectrum_array( size_t ndata, size_t n, spectralwindowfunc* w, double overF, double stp, double fromFreq, double toFreq, double stepFreq, double lev )		//	idem with resampling
	: double_powerSpectrum( ndata, stp, fromFreq, toFreq, stepFreq ), window(w), lag(n), overFactor(overF), level(lev)
		{	setstats();
		};
	
	void	ComputePowerInfSup( const double_vector1& power, double_vector1*& power_inf, double_vector1*& power_sup )
		{	power_inf = new double_vector1( power.size() );
			power_sup = new double_vector1( power.size() );
			for (size_t i=1; i<=power.size(); i++)
			{	(*power_inf)[i] = power[i] * window->LowError();
				(*power_sup)[i] = power[i] * window->HighError();
			};
		}
	void	ComputeCoherencePowerInfSup( const double_vector1& coherence, double_vector1*& coherence_inf, double_vector1*& coherence_sup )
		{	coherence_inf = new double_vector1( coherence.size() );
			coherence_sup = new double_vector1( coherence.size() );
			for (size_t i=1; i<=coherence.size(); i++)
			{	(*coherence_inf)[i] = lowerCoherence( coherence[i], window->CoherenceError() );
				(*coherence_sup)[i] = upperCoherence( coherence[i], window->CoherenceError() );
			};
		}
	void	ComputePhasePowerInfSup( const double_vector1& ph, const double_vector1& coherence, double_vector1*& ph_inf, double_vector1*& ph_sup )
		{	ph_inf = new double_vector1( ph.size() );
			ph_sup = new double_vector1( ph.size() );
			for (size_t i=1; i<=ph.size(); i++)
			{	(*ph_inf)[i] = ph[i] - deltaPhase( coherence[i], window->PhaseErrorArg() );
				(*ph_sup)[i] = ph[i] + deltaPhase( coherence[i], window->PhaseErrorArg() );
			};
		}
		
	virtual char*	name()		{	return (char *)"BTukey";	};
	virtual string&	comment();
};


class MaxESpectrum_array : public double_powerSpectrum	{
protected:
	size_t	nbcoeff;
	virtual	void do_compute( const double_vector1& data, double_vector1& power )
		{	double_vector1	y( data.length() );		y.do_copy( data, abstract_valarray_0<double>::remove_mean );
			double_vector1	cofH( nbcoeff );
			double	pm;
			MyMath::memcof( y, pm, cofH );
				
			for (size_t i=1; i<=x_scale.length(); i++)
			{	double	fdt = x_scale[i] * xstep;
				power[i] = MyMath::evlmem( fdt, cofH, pm );
				if ( (fdt < 1e-9) && (power[i] > 1e8)  )		//	manière peu élégante d'éliminer	
					power[i] = 0;								//  la singularité éventuelle en zéro	
			};
		};
public:	
	MaxESpectrum_array( size_t ndata, size_t nb, double stp, double fromFreq, double toFreq, double stepFreq )
	: double_powerSpectrum( ndata, stp, fromFreq, toFreq, stepFreq ), nbcoeff(nb)
		{};
		
	virtual char*	name()		{	return (char *)"MaxEntropy";	};
	virtual string&	comment()	{	comm = "";
									comm += "MaxEntropy spectrum";
									return	comm;
								};
};

class DPSS_array : public double_matrix1	{
private:
	void	DoCompute();
	void	DoComputeVap();
	
public:
	double_vector1		vap;		//	vap 	= NR_dvector( 1, ncol );
	double_vector1		un_vap;		//	un_vap 	= NR_dvector( 1, ncol );
	double				wwidth;
	size_t				nndp;
	//size_t	ndpss;	= dim1();
	//size_t	npo;	= dim2();
	
	DPSS_array():double_matrix1(0,0),vap(0),un_vap(0) 	{	wwidth = 0;	};
	void	Compute( size_t m, size_t n, size_t _nndp, double _wwidth )
			{	if (m != dim1() || n != dim2() || _wwidth != wwidth)
				{	resize2D( m, n );
					vap.resize_but_keep( m );
					un_vap.resize_but_keep( m );
					nndp = _nndp;	wwidth = _wwidth;
					DoCompute();
					DoComputeVap();
				}
				else if (_nndp != nndp)
				{	DoComputeVap();
				}
			};
};

class MTMSpectrum_array : public double_powerSpectrum	{

private:
	static	DPSS_array	dpss;
	double	wwidth;		//	aussi dans dpss
	size_t	nndp;		//	aussi dans dpss
	size_t	ndpss;		//	aussi dans dpss
	size_t	npo;		//	aussi dans dpss

	size_t	fRange_nn;
	double	fRange_delta_fre;
	size_t	fRange_fcount;
	size_t	fRange_nf0;
	double	varian;
	
	double_matrix1	fftsig;		//	fftsig 	= NR_dmatrix( 1, m, 1, 2 * fRange_fcount );
	double_vector1	varlog;		//	varlog 	= NR_dvector( 1, fcount );
	
	void	freqRange( double delta_f, double delta_t, double f_min, double f_max );
	void	eigenft( const double_vector1& y );
	void	adapt_MTM( double_vector1& power );

	double	f_trans( double x, double nu1, double nu2 )
		{	return (1 - betai(nu2 * 0.5, nu1 * 0.5, nu2 / (nu2 + nu1 * x)));	};

protected:
	virtual	void do_compute( const double_vector1& data, double_vector1& power )
		{	double_vector1 zy( data.length() );	zy.do_copy( data, abstract_valarray_0<double>::remove_mean );
			
			varian = 0;
			for (size_t i=1; i<=zy.length(); i++)	varian += sqr(zy[i]);			//	 += zy[i];		//	au carré ???
			varian /= (zy.length()-1);

			eigenft( zy );
			adapt_MTM( power );
		};

public:	
	MTMSpectrum_array( size_t ndata, size_t m, double wnx, double stp, double fromFreq, double toFreq, double stepFreq )
	: double_powerSpectrum( ndata, stp, fromFreq, toFreq, stepFreq )
		{	nndp 	= nextPower2( std::max(long_round(1 / (2 * xstep * stepFreq)), 255L ) );
			ndpss 	= m;
			npo 	= data_length;
			wwidth	= wnx/npo;
			dpss.Compute( ndpss, npo, nndp, wwidth );			//	compute dpss, if different from previous call
			freqRange( stepFreq, xstep, fromFreq, toFreq );
			
		};

	virtual char*	name()		{	return (char *)"MTM";	};
	virtual string&	comment()	{	comm = "";
									comm += "MTM spectrum";
									return	comm;
								};
	
	void	ComputePowerInfSup( const double_vector1& power, double_vector1*& power_inf, double_vector1*& power_sup )
		{	power_inf = new double_vector1( fRange_fcount );
			power_sup = new double_vector1( fRange_fcount );
			double	xx1 = 1.96;
			for (size_t i=1; i<=fRange_fcount; i++)
			{	(*power_inf)[i] = power[i] * exp( -xx1 * sqrt(varlog[i]) );
				(*power_sup)[i] = power[i] * exp(  xx1 * sqrt(varlog[i]) );
			};
		};
		
	void	Harmonic( double_vector1*& amp, double_vector1*& conf );
};



class SpectralAnalyzer	{			//	takes care of the different methods, with different options
public:
	enum	method	{ Periodogram, BTukey, MaxEntropy, MTM };
	
public:
	method					theMethod;
	size_t					y_length;
	double					xstep;			//	the time step
	
	spectralwindowfunc*		window;
	size_t					iParam;			//	integer parameter
	double					rParam;			//	real parameter
	
	double			fromFreq, toFreq, stepFreq;		//	freq. scale
	bool			resampleFreq;
	
	double			overSampling;			//	for BTukey
	double			statConfLevel;			//	for BTukey
	
	bool			inf_sup_output;			//	outputs with inf_sup
	bool			MTM_Ftest_output;
	bool			MTM_Ampl_output;
	bool			CrossPhase_output;
	bool			DoCrossAnalysis;
	
	void set_inf_sup_Output( bool b )		{	inf_sup_output = b;		};
	void set_MTM_Ftest_Output( bool b )		{	MTM_Ftest_output = b;	};
	void set_MTM_Ampl_Output( bool b )		{	MTM_Ampl_output = b;	};
	void set_CrossPhase_Output( bool b )	{	CrossPhase_output = b;	};
	
	
private:
	MyMath::Periodogram_array*		period;
	MyMath::BTSpectrum_array*		bt;
	MyMath::MaxESpectrum_array*		maxE;
	MyMath::MTMSpectrum_array*		mtm;
	
public:
	void SetUpMethod()
		{	switch ( theMethod )
			{	case Periodogram:
					if (resampleFreq)	period = new MyMath::Periodogram_array( y_length, window, xstep, fromFreq, toFreq, stepFreq );
					else				period = new MyMath::Periodogram_array( y_length, window, xstep );
					spectralMethod = period;
					break;
				case BTukey:				//	Param = Nb of lags
                    if (resampleFreq)
                        bt = new MyMath::BTSpectrum_array( y_length, iParam, window, overSampling, xstep, fromFreq, toFreq, stepFreq, statConfLevel );
                    else
                        bt = new MyMath::BTSpectrum_array( y_length, iParam, window, overSampling, xstep, statConfLevel );
					spectralMethod = bt;
					break;
				case MaxEntropy:			//	Param = Nb of coeffs
					maxE = new MyMath::MaxESpectrum_array( y_length, iParam, xstep, fromFreq, toFreq, stepFreq );
					spectralMethod = maxE;
					break;
				case MTM:					//	Param = Nb of windows
					mtm = new MyMath::MTMSpectrum_array( y_length, iParam, rParam, xstep, fromFreq, toFreq, stepFreq );
					spectralMethod = mtm;
					break;
			};
			DoCrossAnalysis = false;		//	default
		};
	void SetUpCrossAnalysis( size_t n )
		{	DoCrossAnalysis = true;
			if (DoCrossAnalysis && n>0)
				spectralMethod->SetUpCrossAnalysis( n );
		}
		
public:
	double_powerSpectrum*	spectralMethod;		//	includes spectrum and x_scale

	double_vector1*			spectrum;
	double_vector1*			coherence;
	double_vector1*			phase;
	double_vector1*			power_inf;
	double_vector1*			power_sup;
	double_vector1*			coherence_inf;
	double_vector1*			phase_inf;
	double_vector1*			coherence_sup;
	double_vector1*			phase_sup;
	double_vector1*			MTM_ftest;
	double_vector1*			MTM_ampl;

	SpectralAnalyzer( method m, double stp, size_t n )	
		{	theMethod	= m;
			xstep		= stp;
			y_length	= n;
			window 		= NULL;
			iParam 		= 0;
			rParam 		= 0;
			fromFreq	= 0;
			toFreq		= 0;
			stepFreq	= 0;
			resampleFreq 	= false;
			overSampling	= 1;
			statConfLevel	= 1;
			
			inf_sup_output = false;
			
			spectralMethod  = NULL;
			spectrum  = NULL;
			power_inf = NULL;
			power_sup = NULL;
			coherence = NULL;
			phase  	  = NULL;
			coherence_inf = NULL;
			phase_inf  	  = NULL;
			coherence_sup = NULL;
			phase_sup  	  = NULL;
			MTM_ftest = NULL;
			MTM_ampl  = NULL;
		}
		
	~SpectralAnalyzer()
		{	delete	spectralMethod;
		}
		
	void setWindow( spectralwindowfunc* w )
		{	if ( (theMethod == Periodogram) || (theMethod == BTukey) )
				window 		= w;
		}
	void setIParam( size_t p )
		{	if ( (theMethod == BTukey) || (theMethod == MaxEntropy) || (theMethod == MTM) )
				iParam 		= p;
		}
	void setRParam( double r )
        {	if ( theMethod == MTM )
				rParam 		= r;
		}
	void setFScale( double f0, double f1, double df )
		{	fromFreq	= f0;
			toFreq		= f1;
			stepFreq	= df;
		}
	void setResampleFlag( bool r )
		{	resampleFreq = r;
		}
	void setOverSampling( double over )
        {	if ( theMethod == BTukey )
				overSampling = over;
		}
	void setConfLevel( double c )
        {	if ( theMethod == BTukey )
				statConfLevel = c;
		}
	
	void ComputePowerSpectrum( const double_vector& y )
		{	//spectrum = new double_vector1( spectralMethod->x_scale.length() );
			spectralMethod->ComputePower( y, spectrum );
			
			switch ( theMethod )		//	associated optional outputs
			{	case BTukey:
					if (inf_sup_output)		bt->ComputePowerInfSup( *spectrum, power_inf, power_sup );
					break;
				case MTM:
					if (inf_sup_output)		mtm->ComputePowerInfSup( *spectrum, power_inf, power_sup );
					if (MTM_Ftest_output)
					{	mtm->Harmonic( MTM_ampl, MTM_ftest );
						if (!MTM_Ampl_output)	{	delete	MTM_ampl; 	MTM_ampl = NULL;	};
					}
					break;
				case MaxEntropy:
				case Periodogram:
					break;
			};
		}
		
	void ComputeCrossPowerSpectrum( const double_vector& yy, const double_vector& zz, const double_vector& sy, const double_vector& sz )
		{	switch ( theMethod )
			{	case BTukey:
					bt->ComputeCrossPower( yy, zz, sy, sz, coherence, phase );
					if (!CrossPhase_output)	{	delete	phase; 	phase = NULL;	};
					if (inf_sup_output)
					{	bt->ComputeCoherencePowerInfSup( *coherence, coherence_inf, coherence_sup );
						if (CrossPhase_output)
							bt->ComputePhasePowerInfSup( *phase, *coherence, phase_inf, phase_sup );
					}
					break;
				case MTM:
				case MaxEntropy:
				case Periodogram:
					break;
			};
		}
};






}		//	namespace MyMath


