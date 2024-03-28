
#include "myMath.h"
//#include "NRecipes.h"

#include <time.h>



void	MyMath::ThrowError( char* s )
{	
//#ifdef Debug_Throw
//#pragma unused	(s)
//	Throw_(err_AssertFailed);
//#else
	throw	GenericError(s);
//#endif
};
void	MyMath::ThrowErrorIf( bool x )
{	if (x)	{	ThrowError((char *)"throw if" );		};	};

void	MyMath::ThrowErrorIfNot( bool x )
{	if (!x)	{	ThrowError((char *)"throw if not" );	};	};


void	MyMath::hunt( const a_double_vector1& xx, const double& x, unsigned long& jlo )
{
	xx.hunt( x, jlo );
	
	//	d'apres Numerical Recipes 'hunt'
/*	
	unsigned long n = xx.length();
	unsigned long jm,jhi,inc;
	Boolean	ascnd = (n == 1) || (xx[n] > xx[1]);		//	for n==1, ascnd is assumed true !!
	
	if (jlo <= 0 || jlo > n)		{	jlo=0;	jhi=n+1;	}
	else {
		inc=1;
		if (x >= xx[jlo] == ascnd) {
			if (jlo == n) return;
			jhi=(jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				jlo=jhi;
				inc += inc;
				jhi=(jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (jlo == 1) {
				jlo=0;
				return;
			}
			jhi=(jlo)--;
			while (x < xx[jlo] == ascnd) {
				jhi=(jlo);
				inc <<= 1;
				if (inc >= jhi) {
					jlo=0;
					break;
				}
				else jlo=jhi-inc;
			}
		}
	}
	
	while (jhi-jlo != 1) {
		jm=(jhi+jlo) >> 1;
		if (x > xx[jm] == ascnd)
			jlo=jm;
		else
			jhi=jm;
	}
*/
};


double	MyMath::bissectbracket( const a_double_vector1& x, const double& value, unsigned long& klo )
{	
	return x.bissectbracket( value, klo );
	
/*	hunt( x, value, klo );						//	output: 0 <= klo <= n
	
	unsigned long n = x.length();
	if ( klo >= n )			klo = n-1;
	else if ( klo <= 0 )	klo = 1;			//	output: 1 <= klo <= n-1
	
	return	( x[ klo+1 ] - x[ klo ] );			//	always OK, except for n = 1
*/
};


#pragma mark -

double	MyMath::linear_interp( const a_double_vector1& XA, const a_double_vector1& YA, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowErrorIf(H == 0);
	unsigned long	khi = klo+1;
	return	((XA[khi] - X) * YA[klo] + (X - XA[klo]) * YA[khi]) / H;
};

double	MyMath::linear_interp_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowErrorIf (H == 0);
	unsigned long	khi = klo+1;
	return (YA[khi]-YA[klo]) / H;
};

double	MyMath::linear_interp_integ1( const a_double_vector1& x, const a_double_vector1& y, double a, double b, unsigned long& klo )
{	//		ia, ib:		indices min max tels que a < x[ia] <...< x[ib] < b
	//					ia == 0:	x[n] < a < b
	//					ib == 0:	a < b < x[1]
	//					ia > ib:	x[ib] < a < b < x[ia]
	
	double			ya = linear_interp( x, y, a, klo );			//	valeur pour a
	unsigned long	ia = klo+1;
	if ( x[klo] > a )			ia = klo;		//	ya == extrapolation ( klo == 1, a < x[1] )
	else if ( x[klo+1] < a )	ia = 0;			//	rien entre a et b ( ya et yb == extrapolations: klo+1 == n, x[n] < a < b )
	
	double			yb = linear_interp( x, y, b, klo );			//	valeur pour b
	unsigned long	ib = klo;
	if ( x[klo] > b )			ib = 0;			//	rien entre a et b ( ya et yb == extrapolations: klo == 1, a < b < x[1] )
	else if ( x[klo+1] < b )	ib = klo+1;		//	yb == extrapolation ( klo+1 == n, x[n] < b )
	
	double	s;
	if ( (ia==0) || (ib==0) || (ia > ib) )				//	rien entre a et b
		s = (yb + ya) * (b - a);
	else
	{	s = (y[ia] + ya) * (x[ia] - a) + (yb + y[ib]) * (b - x[ib]);
		for (unsigned long i=ia; i<=(ib-1); i++)
			s += (y[i+1] + y[i]) * (x[i+1] - x[i]);
	}
	s *= 0.5;
	return s;
};


#pragma mark -




void	MyMath::splin_interpolfunc::spline_init( const a_double_vector1& x, const a_double_vector1& y, a_double_vector1& y2, double yp1, double ypn )
{
	unsigned long n = x.length();
	double_vector1	u(n-1);
	if (yp1 == unused_spline_yp)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (unsigned long i=2;i<=n-1;i++) {
		double	sig	= (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		double	p	= sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	
	double	qn,un;
	if (ypn == unused_spline_yp)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (unsigned long k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
};

double	MyMath::splin_interpolfunc::spline_interp( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2, double t, unsigned long& klo )
{	double	H = bissectbracket( x, t, klo );
	ThrowErrorIf (H == 0);
	unsigned long	khi = klo+1;
	
	double	A = (x[khi] - t) / H;
	double	B = (t - x[klo]) / H;
	return	( A*y[klo] + B*y[khi] + (A*(sqr(A)-1)*y2[klo] + B * (sqr(B) - 1) * y2[khi]) * sqr(H) / 6 );
}


double	MyMath::splin_interpolfunc::spline_interp_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2,
			double a, double b, unsigned long& klo )
{	
	//	ia, ib	:	indices min max tels que a < x[ia] <...< x[ib] < b
	double	s = 0;
	
	double	H = bissectbracket( x, a, klo );
	ThrowErrorIf (H == 0);
	unsigned long	ia = klo;
	s -= spline_integ_X1_X( a, x, y, y2, y2, klo );		//	Sum de a -> x[ia]
	
	H = bissectbracket( x, b, klo );
	ThrowErrorIf (H == 0);
	unsigned long	ib = klo+1;
	s -= spline_integ_X_X2( b, x, y, y2, y2, klo );		//	Sum de x[ib] -> b
		
	if (ib < ia)
		for (unsigned long i=ib; i<=(ia-1); i++)					//	-Sum de x[ib] -> x[ia]
			s -= spline_integ_X1_X2( x, y, y2, y2, i );
	else
		for (unsigned long i=ia; i<=(ib-1); i++)					//	Sum de x[ia] -> x[ib]
			s += spline_integ_X1_X2( x, y, y2, y2, i );
	
	return s;
};

double	MyMath::splin_interpolfunc::spline_interp_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowErrorIf (H == 0);
	unsigned long	khi = klo+1;
	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	(YA[khi] - YA[klo]) / H + ((1-3*sqr(A)) * Y2A[klo] - (1-3*sqr(B)) * Y2A[khi]) * H / 6;
};

double	MyMath::splin_interpolfunc::spline_interp_deriv2( const a_double_vector1& XA, const a_double_vector1&, const a_double_vector1& Y2A, double X, unsigned long& klo )
{	return	linear_interp( XA, Y2A, X, klo );
};

double	MyMath::splin_interpolfunc::spline_interp_deriv3( const a_double_vector1& XA, const a_double_vector1&, const a_double_vector1& Y2A, double X, unsigned long& klo )
{	return	linear_interp_deriv1( XA, Y2A, X, klo );
};

#pragma mark -

Boolean	MyMath::inc_splin_interpolfunc::inc_spline_check( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& yR, const a_double_vector1& yL, double eps )
{	Boolean test = true;
	for (unsigned long i=1; (i<=x.length()-1)&&(test); i++)
	{	double	dx = x[i+1]-x[i];
		double	dh = (y[i+1]-y[i])/dx;
		double	R = (dh - eps) * 6/dx;
		double	u = R + yR[i] + 2*yL[i+1];
		double	v = R - 2*yR[i] - yL[i+1];
		double	w = (sqr(yL[i+1]) + sqr(yR[i]) + yL[i+1]*yR[i])/(yL[i+1] - yR[i]);
		if ( dx <= 0 )	test = false;
		if ( v < 0 )	test = false;
		if ( u < 0 )	test = false;
		if ( (yR[i]<0) && (yL[i+1]>0) && ( w > R ) )	test = false;
	}
	return test;
};

void	MyMath::inc_splin_interpolfunc::inc_spline_init( const a_double_vector1& x, const a_double_vector1& y, a_double_vector1& y2L, a_double_vector1& y2R,
						double minRate, int t )
{
	int		spType = (t == inc_splin_interpolfunc::fast ? 1 : 4);			//	1 ou bien 4
	
	if ( !inc_spline_check( x, y, y2L, y2R, minRate ) )
	{	unsigned long	n    = x.length();
		unsigned long	nvar = 2*n;
		unsigned long	ncon = n + spType*(n-1);
			
		double_vector0	R(n+1);			//	0 to n
		double_vector1	Lbd(n-1);		//	1 to n-1
		
		double	dxim1 = x[2]-x[1];
		R[0] = R[n] = 0;
		R[1] = ( ((y[2]-y[1])/dxim1) - minRate ) * 6/dxim1;
		Lbd[n-1] = 1;
		
		for (unsigned long i=2; i<=n-1; i++)
		{	double	dx = x[i+1]-x[i];
			R[i] = ( ((y[i+1]-y[i])/dx) - minRate ) * 6/dx;
			Lbd[i-1] = dx/dxim1;
			dxim1 = dx;
		}
		
		double_matrix1	a( ncon+2, nvar+1 );
		ulong_vector1	posv( ncon );
		ulong_vector1	zerov( nvar );
		
		//	variables:	v[1]..v[N] r[1]..r[N]
		
		unsigned long 	Voffset = 1;
		unsigned long 	Roffset = Voffset + n;
			
		a(1,1) = 0;				//	Maximize v[1] + Sum{(1+Lbd[i-1])v[i]; i=2,n-1} + Lbd[n-1]v[n] - 2 Sum{r[i]; i=1,n}
		a(1,Voffset+1) = 1;														//	v[1]
		for (unsigned long i=2; i<=n-1; i++)	a(1,Voffset+i) = 1 + Lbd[i-1];	//	v[i]
		a(1,Voffset+n) = 1;														//	v[n]: Lbd[n-1] = 1.
		for (unsigned long i=1; i<=n; i++)		a(1,Roffset+i) = -2.0;			//	r[i]
		
		unsigned long	offset = 1;
			
		for (unsigned long i=1; i<=n; i++)
		{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
			a(offset+i,1) = 3*(R[i]+R[i-1]);			//	v[i-1] + 2(1+Lbd[i-1])v[i] + Lbd[i]v[i+1] - 3r[i] < 3 (R[i]+R[i-1])
			
			if (i<n)
			{	a(offset+i,Voffset+i)   -= 2.0;			//	v[i]
				a(offset+i,Voffset+i+1) -= Lbd[i];			//	v[i+1]
			}
			if (i>1)
			{	a(offset+i,Voffset+i-1) -= 1.0;			//	v[i-1]
				a(offset+i,Voffset+i)   -= 2*Lbd[i-1];		//	v[i]
			}
			a(offset+i,Roffset+i) -= -3.0;					//	r[i]
		}
		offset += n;
			
		if (spType == 1)
		{	for (unsigned long i=1; i<=n-1; i++)
			{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
				a(offset+i,1) = 3*R[i];					//	u[i] + v[i] < 3 Ri
				a(offset+i,Voffset+i)   -= 1.0;			//	v[i]
				a(offset+i,Voffset+i+1) -= Lbd[i];			//	u[i] = Lbd[i] v[i+1]
			}
			offset += n-1;
		}
		if (spType == 4)
		{	for (unsigned long i=1; i<=n-1; i++)
			{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
				a(offset+i,1) = 3*R[i];					//	v[i] - u[i] < 3 Ri
				a(offset+i,Voffset+i)   -= 1.0;			//	v[i]
				a(offset+i,Voffset+i+1) -= -Lbd[i];		//	u[i] = Lbd[i] v[i+1]
			}
			offset += n-1;
			
			for (unsigned long i=1; i<=n-1; i++)
			{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
				a(offset+i,1) = 3*R[i];					//	u[i] - v[i] < 3 Ri
				a(offset+i,Voffset+i)   -= -1.0;		//	v[i]
				a(offset+i,Voffset+i+1) -= Lbd[i];		//	u[i] = Lbd[i] v[i+1]
			}
			offset += n-1;
			
			for (unsigned long i=1; i<=n-1; i++)
			{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
				a(offset+i,1) = 9*R[i];					//	2 v[i] + u[i] < 9 Ri
				a(offset+i,Voffset+i)   -= 2.0;			//	v[i]
				a(offset+i,Voffset+i+1) -= Lbd[i];		//	u[i] = Lbd[i] v[i+1]
			}
			offset += n-1;
			
			for (unsigned long i=1; i<=n-1; i++)
			{	for (unsigned long j=1; j<=nvar; j++)		a(offset+i,1+j) = 0;		//	mise a zero
				a(offset+i,1) = 9*R[i];					//	v[i] + 2 u[i] < 9 Ri
				a(offset+i,Voffset+i)   -= 1.0;			//	v[i]
				a(offset+i,Voffset+i+1) -= 2*Lbd[i];	//	u[i] = Lbd[i] v[i+1]
			}
			offset += n-1;
		}
			
		int	icase;
		Boolean minimumFound = simplx( a, ncon, nvar, ncon, 0, 0, icase, zerov, posv );
		ThrowErrorIfNot ( minimumFound );
		
		double_vector1	v( nvar );									//	v[1]..v[N] r[1]..r[N]
		for (unsigned long i=1; i<=nvar; i++)		v[i] = 0;
		for (unsigned long i=1; i<=ncon; i++)
			if (posv[i] <= nvar )	v[ posv[i] ] = a(i+1,1);
	
		y2R[n] = 0;
		y2L[1] = 0;
		for (unsigned long i=1; i<=n-1; i++)
		{	y2L[i+1] = (2*Lbd[i]*v[i+1] + v[i])/3 - R[i];		//	= (2*u[i]+v[i])/3 - R[i]	i=1,n-1
			y2R[i]   = R[i] - (Lbd[i]*v[i+1] + 2*v[i])/3;		//	= R[i] - (u[i]+2*v[i])/3	i=1,n-1
		};
	}			//	if ( !inc_spline_check( x, y, y2L, y2R, minRate ) )
	
};


double	MyMath::splin_extra_lin_interpolfunc::spline_interp_extra_lin( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end )
{	// idem above, but extrapolation is linear for natural splines (if y" == 0 at the ends)
	
	unsigned long	khi;
	int				extra;
	double			H = spline_biss_bracket( XA, X, klo, khi, extra );
	
	if ( extra == 1 && ( dont_check_end || Y2A_L[khi] == 0 ))
	{	double	dy = (YA[khi] - YA[klo]) / H + ( Y2A_R[klo] + 2 * Y2A_L[khi]) * H / 6;		//	derivative at the end
		return	YA[khi] + dy * (X - XA[khi]);
	}
	if ( extra == -1 && ( dont_check_end || Y2A_R[klo] == 0 ))
	{	double	dy = (YA[khi] - YA[klo]) / H - (2 * Y2A_R[klo] + Y2A_L[khi]) * H / 6;		//	derivative at the beginning
		return	YA[klo] + dy * (X - XA[klo]);
	}
	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	( A*YA[klo] + B*YA[khi] + (A*(sqr(A)-1)*Y2A_R[klo] + B * (sqr(B) - 1) * Y2A_L[khi]) * sqr(H) / 6 );
}

double	MyMath::splin_extra_lin_interpolfunc::spline_interp_extra_lin_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end )
{	
	unsigned long	khi;
	int				extra;
	double			H = spline_biss_bracket( XA, X, klo, khi, extra );
	
	if ( extra == 1 && ( dont_check_end || Y2A_L[khi] == 0 ))
		return	(YA[khi] - YA[klo]) / H + ( Y2A_R[klo] + 2 * Y2A_L[khi]) * H / 6;		//	derivative at the end

	if ( extra == -1 && ( dont_check_end || Y2A_R[klo] == 0 ))
		return  (YA[khi] - YA[klo]) / H - (2 * Y2A_R[klo] + Y2A_L[khi]) * H / 6;		//	derivative at the beginning

	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	(YA[khi] - YA[klo]) / H + ((1-3*sqr(A)) * Y2A_R[klo] - (1-3*sqr(B)) * Y2A_L[khi]) * H / 6;
};

double	MyMath::splin_extra_lin_interpolfunc::spline_interp_extra_lin_deriv2( const a_double_vector1& XA, const a_double_vector1&, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end )
{	
	unsigned long	khi;
	int				extra;
	double			H = spline_biss_bracket( XA, X, klo, khi, extra );
	
	if ( extra == 1 && ( dont_check_end || Y2A_L[khi] == 0 ))	return	0;		//	2nd	derivative at the end
	if ( extra == -1 && ( dont_check_end || Y2A_R[klo] == 0 ))	return  0;		//	2nd derivative at the beginning
	
	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	A * Y2A_R[klo] + B * Y2A_L[khi];
};

double	MyMath::splin_extra_lin_interpolfunc::spline_interp_extra_lin_deriv3( const a_double_vector1& XA, const a_double_vector1&, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end )
{	
	unsigned long	khi;
	int				extra;
	double			H = spline_biss_bracket( XA, X, klo, khi, extra );
	
	if ( extra == 1 && ( dont_check_end || Y2A_L[khi] == 0 ))	return	0;		//	3rd	derivative at the end
	if ( extra == -1 && ( dont_check_end || Y2A_R[klo] == 0 ))	return  0;		//	3rd derivative at the beginning
	
	return	( - Y2A_R[klo] + Y2A_L[khi] ) / H;
};

double	MyMath::splin_extra_lin_interpolfunc::spline_interp_extra_lin_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2_R, const a_double_vector1& y2_L, double a, double b, unsigned long& klo, bool dont_check_end )
{	
	//	ia, ib	:	indices min max tels que a < x[ia] <...< x[ib] < b
	double			s = 0;
	unsigned long	khi;
	int				extra;
	
	double			H = spline_biss_bracket( x, a, klo, khi, extra );
	unsigned long	ia = klo;
	
			//	Sum de a -> x[ia]
	if ( extra == 1 && ( dont_check_end || y2_L[khi] == 0 ))
	{	double	dy = (y[khi] - y[klo]) / H + ( y2_R[klo] + 2 * y2_L[khi]) * H / 6;		//	derivative at the end
		ia = khi;
		s += (a - x[ia]) * (y[ia] + 0.5*dy * (a - x[ia]));
	}
	else if ( extra == -1 && ( dont_check_end || y2_R[klo] == 0 ))
	{	double	dy = (y[khi] - y[klo]) / H - (2 * y2_R[klo] + y2_L[khi]) * H / 6;		//	derivative at the beginning
		s -= (a - x[ia]) * (y[ia] + 0.5*dy * (a - x[ia]));
	}
	else	s -= spline_integ_X1_X( a, x, y, y2_R, y2_L, klo );
	
	H = spline_biss_bracket( x, b, klo, khi, extra );
	unsigned long	ib = khi;
	
			//	Sum de x[ib] -> b
	if ( extra == 1 && ( dont_check_end || y2_L[khi] == 0 ))
	{	double	dy = (y[khi] - y[klo]) / H + ( y2_R[klo] + 2 * y2_L[khi]) * H / 6;		//	derivative at the end
		s += (a - x[ib]) * (y[ib] + 0.5*dy * (a - x[ib]));
	}
	else if ( extra == -1 && ( dont_check_end || y2_R[klo] == 0 ))
	{	double	dy = (y[khi] - y[klo]) / H - (2 * y2_R[klo] + y2_L[khi]) * H / 6;		//	derivative at the beginning
		ib = klo;
		s -= (a - x[ib]) * (y[ib] + 0.5*dy * (a - x[ib]));
	}
	else	s -= spline_integ_X_X2( b, x, y, y2_R, y2_L, klo );
		
	if (ib < ia)
		for (unsigned long i=ib; i<=(ia-1); i++)					//	-Sum de x[ib] -> x[ia]
			s -= spline_integ_X1_X2( x, y, y2_R, y2_L, i );
	else
		for (unsigned long i=ia; i<=(ib-1); i++)					//	Sum de x[ia] -> x[ib]
			s += spline_integ_X1_X2( x, y, y2_R, y2_L, i );
	
	return s;
};

/*		//	only use the "extra_lin" options

double	MyMath::inc_splin_interpolfunc::spline_interp_integ1( const double_vector1& x, const double_vector1& y, const double_vector1& y2_R, const double_vector1& y2_L, double a, double b, unsigned long& klo )
{	
	//	ia, ib	:	indices min max tels que a < x[ia] <...< x[ib] < b
	double	s = 0;
	
	double	H = bissectbracket( x, a, klo );
	ThrowIf_ (H == 0);
	unsigned long	ia = klo;
	s -= spline_integ_X1_X( a, x, y, y2_R, y2_L, klo );		//	Sum de a -> x[ia]
	
	H = bissectbracket( x, b, klo );
	ThrowIf_ (H == 0);
	unsigned long	ib = klo+1;
	s -= spline_integ_X_X2( b, x, y, y2_R, y2_L, klo );		//	Sum de x[ib] -> b
		
	if (ib < ia)
		for (unsigned long i=ib; i<=(ia-1); i++)					//	-Sum de x[ib] -> x[ia]
			s -= spline_integ_X1_X2( x, y, y2_R, y2_L, i );
	else
		for (unsigned long i=ia; i<=(ib-1); i++)					//	Sum de x[ia] -> x[ib]
			s += spline_integ_X1_X2( x, y, y2_R, y2_L, i );
	
	return s;
};

double	MyMath::inc_splin_interpolfunc::spline_interp( const double_vector1& x, const double_vector1& y, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double t, unsigned long& klo )
{	double	H = bissectbracket( x, t, klo );
	ThrowIf_ (H == 0);
	unsigned long	khi = klo+1;
	double	A = (x[khi] - t) / H;
	double	B = (t - x[klo]) / H;
	return	( A*y[klo] + B*y[khi] + (A *(sqr(A)-1)*Y2A_R[klo] + B*(sqr(B)-1)*Y2A_L[khi])*sqr(H)/6 );
}

double	MyMath::inc_splin_interpolfunc::spline_interp_deriv1( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowIf_ (H == 0);
	unsigned long	khi = klo+1;
	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	(YA[khi] - YA[klo]) / H + ((1-3*sqr(A)) * Y2A_R[klo] - (1-3*sqr(B)) * Y2A_L[khi]) * H / 6;
};

double	MyMath::inc_splin_interpolfunc::spline_interp_deriv2( const double_vector1& XA, const double_vector1&, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowIf_ (H == 0);
	unsigned long	khi = klo+1;
	double	A = (XA[khi] - X) / H;
	double	B = (X - XA[klo]) / H;
	return	A * Y2A_R[klo] + B * Y2A_L[khi];
};

double	MyMath::inc_splin_interpolfunc::spline_interp_deriv3( const double_vector1& XA, const double_vector1&, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double X, unsigned long& klo )
{	double	H = bissectbracket( XA, X, klo );
	ThrowIf_ (H == 0);
	unsigned long	khi = klo+1;
	return	( - Y2A_R[klo] + Y2A_L[khi] ) / H;
};

*/

#pragma mark -

double	MyMath::stair_interp( const a_double_vector1&, const a_double_vector1& YA, const a_double_vector1& midP, double X, unsigned long& klo )
{	hunt( midP, X, klo );			//	returns 0 <= klo <= N-1 (== midP.length())
	return YA[ klo+1 ];
};

double	MyMath::stair_interp_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& midP,
			double a, double b, unsigned long& klo )
{	//		ia, ib	:	indices min max tels que a < x[ia] <...< x[ib] < b
	unsigned long	n =x.length();
	double			s = 0;
	
	double	ya = stair_interp( x, y, midP, a, klo );	//	returns 0 <= klo <= n-1
	unsigned long	ia = klo+1;							//	midP[ia] juste apres a
	if (ia == n)	ia--;								//		ou bien juste avant...
	s += ya * (midP[ia] - a);							//	Sum de a -> midP[ia]
	
	
	double	yb = stair_interp( x, y, midP, b, klo );
	unsigned long	ib = klo;							//	midP[ib] juste avant b
	if (ib == 0)	ib++;								//		ou bien juste apres...
	s += yb * (b - midP[ib]);							//	Sum de midP[ib] -> b
	
	if (ib < ia)										//	si midP[ib] < a < b < midP[ia]
		for (unsigned long i=ib; i<=(ia-1); i++)
			s -= (midP[i+1] - midP[i]) * y[i+1];
	else
		for (unsigned long i=ia; i<=(ib-1); i++)
			s += (midP[i+1] - midP[i]) * y[i+1];		//	Sum de midP[ia] -> midP[ib]
	
	return s;
};


MyMath::common_regular_scale::common_regular_scale( const val1Darray<double_vector1*,1>& in_x, bool on_intersection )
{	
		//	first find common limits and a common step (min of the means)
	double	maxi_n = -1e222;
	double	mini_1 = 1e222;
	double	maxi_1 = -1e222;
	double	mini_n = 1e222;
	double	minstep = 1e222;
	for (size_t i=1; i<=in_x.length(); i++)
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
	for (size_t i=2; i<=in_x.length() && reg; i++)		reg = in_x(1)->check_same_regular_step_as( *in_x(i) );
		
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
		just_copy = false;
	}
	else			//	just copy the first one
	{	//x = new MyMath::double_regular_array( (*in_x(1))(1), (*in_x(1))(in_x(1)->length()), in_x(1)->length() );
		//return true;
		
		resize( in_x(1)->length() );
		init_regular( (*in_x(1))(1), (*in_x(1))(in_x(1)->length()) );
		just_copy = true;
	}
};

double	MyMath::Arondi( double x, roundType t, int n )
{	if (x == 0)		return 0;
	double	mz = magnitude10(x);
	double	z  = x/mz;								//	1 <= z < 10
	if (n > 0)										//	n significant digits
	{	//double	powerten = 0;	//std::pow( (double)10, n-1 );
		double	powerten = std::pow( (double)10, (double)(n-1) );
		double	y;	std::modf( z * powerten, &y );
		double	rd = mz * y/powerten;
		double	ru = mz * (y+1)/powerten;
		switch (t) {
			case roundP_down:	return	rd;		break;
			case roundP_up:		return	ru;		break;
			case round_up:		return	(mz > 0 ? ru : rd);		break;
			case round_down:	return	(mz > 0 ? rd : ru);		break;
			case round_nearest:	return	(std::fabs(x-rd) <= std::fabs(x-ru) ? rd : ru );	break;
		}
	}
	if (n==0)	return mz;		//	just the power of 10
	if (n==-1)					//	just 1, 2, 5
	{	int	n1, n2;
		if (z < 2)		{	n1 = 1;		n2 = 2;		}
		else if (z < 5)	{	n1 = 2;		n2 = 5;		}
		else			{	n1 = 5;		n2 = 10;	}
		double	rd = mz * n1;
		double	ru = mz * n2;
		switch (t) {
			case roundP_down:	return	rd;		break;
			case roundP_up:		return	ru;		break;
			case round_up:		return	(mz > 0 ? ru : rd);		break;
			case round_down:	return	(mz > 0 ? rd : ru);		break;
			case round_nearest:	return	(std::fabs(x-rd) <= std::fabs(x-ru) ? rd : ru );	break;
		}
	}
	return 0;
};

#pragma mark -
#pragma mark //	LINEAR PROGRAMMING: Simplex method

void	MyMath::simp1( const double_matrix1& a, unsigned long mm, const ulong_vector1& ll, unsigned long nll, int iabf, unsigned long& kp, double& bmax )
{	
	kp		= ll[1];
	bmax	= a( mm+1, kp+1 );
	for (unsigned long k=2; k<=nll; k++)
	{	double			test;
		if (iabf == 0)	test = a( mm+1, ll[k]+1 ) - bmax;
		else			test = std::fabs( a( mm+1, ll[k]+1 ))-std::fabs( bmax );
		if (test > 0.0)
		{	bmax	= a( mm+1, ll[k]+1 );
			kp		= ll[k];
		}
	}
};


void	MyMath::simp2( const double_matrix1& a, unsigned long n, const ulong_vector1& l2, unsigned long nl2, unsigned long& ip, unsigned long kp, double& q1 )
{	
	const double EPS = 1.0e-6;

	ip = 0;
	for (unsigned long i=1; i<=nl2; i++)
	{	if ( a( l2[i]+1, kp+1 ) < -EPS )
		{	q1 = -a( l2[i]+1, 1 ) / a( l2[i]+1, kp+1 );
			ip = l2[i];
			for (i=i+1; i<=nl2; i++)
			{	unsigned long	ii = l2[i];
				if ( a( ii+1, kp+1 ) < -EPS )
				{	double	q = -a( ii+1, 1 ) / a( ii+1, kp+1 );
					if (q < q1)
					{	ip = ii;	q1 = q;		}
					else if (q == q1)
					{	double	qp, q0;
						for (unsigned long k=1; k<=n; k++)
						{	qp = -a( ip+1, k+1 ) / a( ip+1, kp+1 );
							q0 = -a( ii+1, k+1 ) / a( ii+1, kp+1 );
							if (q0 != qp) break;
						}
						if (q0 < qp)	ip=ii;
	}	}	}	}	}
}

void	MyMath::simp3( double_matrix1& a, unsigned long i1, unsigned long k1, unsigned long ip, unsigned long kp )
{
	double	piv = 1.0 / a( ip+1, kp+1 );
	
	for (unsigned long ii=1; ii<=i1+1; ii++)
		if (ii-1 != ip) {
			a( ii, kp+1 ) *= piv;
			for (unsigned long kk=1; kk<=k1+1; kk++)
				if (kk-1 != kp)
					a( ii, kk ) -= a( ip+1, kk ) * a( ii, kp+1 );
		}
	for (unsigned long kk=1; kk<=k1+1; kk++)
		if (kk-1 != kp) a( ip+1, kk ) *= -piv;
	a( ip+1, kp+1 ) = piv;
}


Boolean	MyMath::simplx( double_matrix1& a, unsigned long m, unsigned long n, unsigned long m1, unsigned long m2,
					unsigned long m3, int&	icase, ulong_vector1& izrov, ulong_vector1& iposv )
{
	const double EPS = 1.0e-6;

	//int i,ip,ir,is,k,kh,kp,m12,nl1,nl2;
	//int *l1,*l2,*l3;
	//float q1,bmax;

	if (m != (m1+m2+m3)) ThrowError((char *)"Bad input constraint counts in simplx");
	
	ulong_vector1 l1(n+1);
	unsigned long nl1 = n;
	for (unsigned long k=1; k<=n; k++)	l1[k] = izrov[k] = k;
	
	ulong_vector1 l2(m);
	unsigned long nl2 = m;
	for (unsigned long i=1; i<=m; i++)
	{	if ( a( i+1, 1 ) < 0.0)		ThrowError((char *)"Bad input tableau in simplx");
		l2[i]	 = i;
		iposv[i] = n+i;
	}
	
	ulong_vector1 l3(m);
	for (unsigned long i=1; i<=m2; i++)	l3[i]=1;
	
	unsigned long	ir=0;
	if (m2+m3)
	{	ir = 1;
		double	q1;
		for (unsigned long k=1; k<=(n+1); k++)
		{	q1 = 0.0;
			for (unsigned long i=m1+1;i<=m;i++) q1 += a( i+1, k );
			a( m+2, k ) = -q1;
		}
		do {
			double			bmax;
			unsigned long	kp;
			unsigned long	ip;
			simp1( a, m+1, l1, nl1, 0, kp, bmax );
			
			if (bmax <= EPS && a( m+2, 1 ) < -EPS)
			{	icase = -1;
				return false;	//	no solution satifies the constraints
			}
			else if (bmax <= EPS && a( m+2, 1 ) <= EPS)
			{	unsigned long m12 = m1+m2+1;
				if (m12 <= m) {
					for (ip=m12; ip<=m; ip++)
					{	if (iposv[ip] == (ip+n))
						{
							simp1( a, ip, l1, nl1, 1, kp, bmax );
							if (bmax > 0.0)
								goto one;
				}	}	}
				
				ir = 0;
				--m12;
				if (m1+1 <= m12)
					for (unsigned long i=m1+1; i<=m12; i++)
						if (l3[i-m1] == 1)
							for (unsigned long k=1; k<=n+1; k++)
								a( i+1, k ) = -a( i+1, k );
				break;
			}
			simp2( a, n, l2, nl2, ip, kp, q1 );
			if (ip == 0)
			{	icase = -1;
				return false;	//	no solution satifies the constraints
			}
			
	one:	simp3( a, m+1, n, ip, kp );
			unsigned long is;
			if (iposv[ip] >= (n+m1+m2+1))
			{	unsigned long k;
				for (k=1; k<=nl1; k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k; is<=nl1; is++) l1[is]=l1[is+1];
				++(a( m+2, kp+1 ));
				for (unsigned long i=1;i<=m+2;i++) a( i, kp+1 ) = -a( i, kp+1 );
			} else {
				if (iposv[ip] >= (n+m1+1))
				{	unsigned long kh = iposv[ip]-m1-n;
					if (l3[kh])
					{	l3[kh] = 0;
						++(a( m+2, kp+1 ));
						for (unsigned long i=1;i<=m+2;i++)
							a( i, kp+1 ) = -a( i, kp+1 );
					}
				}
			}
			is	= izrov[kp];
			izrov[kp] = iposv[ip];
			iposv[ip] = is;
		} while (ir);
	}
	for (;;)
	{	unsigned long	kp;
		double			bmax;
		simp1( a, 0, l1, nl1, 0, kp, bmax );
		if (bmax <= EPS) {
			icase=0;
			return true;	//	best solution found
		}
		unsigned long	ip;
		double			q1;
		simp2( a, n, l2, nl2, ip, kp, q1 );
		if (ip == 0) {
			icase=1;		//	cost function unbounded
			return false;
		}
		simp3( a, m, n, ip, kp );
		unsigned long	is = izrov[kp];
		izrov[kp] = iposv[ip];
		iposv[ip] = is;
	}
}

#pragma mark -
#pragma mark //	elliptic integrals


double MyMath::elle( double phi, double ak )
{	double	s	= std::sin(phi);
	double	cc	= sqr(std::cos(phi));
	double	q	= (1.0-s*ak)*(1.0+s*ak);
	return	s * ( rf(cc,q,1.0) - sqr(s*ak) * rd(cc,q,1.0)/3.0 );
}

double MyMath::ellf( double phi, double ak )
{	double	s	= sin(phi);
	return	s * rf( sqr(cos(phi)), (1.0-s*ak)*(1.0+s*ak), 1.0 );
}


double MyMath::ellpi( double phi, double en, double ak )
{	double	s	= sin(phi);
	double	ess	= en*s*s;
	double	cc	= sqr(cos(phi));
	double	q	=(1.0-s*ak)*(1.0+s*ak);
	return	s * ( rf(cc,q,1.0) - ess * rj(cc,q,1.0,1.0+ess)/3.0 );
}


double	MyMath::rd( double x, double y, double z )
{	
	const double	errtol	= 0.05;
	const double	tiny	= 1.0e-25;
	const double	big		= 4.5e21;
	const double	c1		= (3.0/14.0);
	const double	c2		= (1.0/6.0);
	const double	c3		= (9.0/22.0);
	const double	c4		= (3.0/26.0);
	const double	c5		= (0.25*c3);
	const double	c6		= (1.5*c4);
	
	if (std::min(x,y) < 0.0 || std::min(x+y,z) < tiny || std::max(std::max(x,y),z) > big)
		ThrowError((char *)"invalid arguments in rd");
		
	double	xt = x;
	double	yt = y;
	double	zt = z;
	double	sum = 0.0;
	double	fac = 1.0;
	double	delx, dely, delz, ave;
	
	do {
		double	sqrtx = sqrt(xt);
		double	sqrty = sqrt(yt);
		double	sqrtz = sqrt(zt);
		double	alamb = sqrtx*(sqrty+sqrtz) + sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac	= 0.25 * fac;
		xt	= 0.25 * (xt+alamb);
		yt	= 0.25 * (yt+alamb);
		zt	= 0.25 * (zt+alamb);
		ave	= 0.2 * (xt+yt+3.0*zt);
		delx	= (ave-xt)/ave;
		dely	= (ave-yt)/ave;
		delz	= (ave-zt)/ave;
	}
	while (std::max(std::max(fabs(delx),fabs(dely)),fabs(delz)) > errtol);
	
	double	ea	= delx*dely;
	double	eb	= delz*delz;
	double	ec	= ea-eb;
	double	ed	= ea-6.0*eb;
	double	ee	= ed+ec+ec;
	
	return 3.0*sum + fac*(1.0+ed*(-c1+c5*ed-c6*delz*ee)
		+ delz*(c2*ee+delz*(-c3*ec+delz*c4*ea)))/(ave*sqrt(ave));
}

double	MyMath::rf(double x, double y, double z)
{	
	const double	errtol	= 0.08;
	const double	tiny	= 1.5e-38;
	const double	big		= 3.0e37;
	const double	third	= (1.0/3.0);
	const double	c1		= (1.0/24.0);
	const double	c2		= 0.1;
	const double	c3		= (3.0/44.0);
	const double	c4		= (1.0/14.0);

	if (std::min(std::min(x,y),z) < 0.0 || std::min(std::min(x+y,x+z),y+z) < tiny || std::max(std::max(x,y),z) > big)
		ThrowError((char *)"invalid arguments in rf");

	double	xt = x;
	double	yt = y;
	double	zt = z;
	double	delx, dely, delz, ave;
	
	do {
		double	sqrtx = sqrt(xt);
		double	sqrty = sqrt(yt);
		double	sqrtz = sqrt(zt);
		double	alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt		= 0.25*(xt+alamb);
		yt		= 0.25*(yt+alamb);
		zt		= 0.25*(zt+alamb);
		ave		= third*(xt+yt+zt);
		delx	= (ave-xt)/ave;
		dely	= (ave-yt)/ave;
		delz	= (ave-zt)/ave;
	}
	while (std::max(std::max(fabs(delx),fabs(dely)),fabs(delz)) > errtol);
	
	double	e2	= delx*dely-delz*delz;
	double	e3	= delx*dely*delz;
	return	(1.0+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave);
}

double	MyMath::rj(double x, double y, double z, double p)
{	
	const double	errtol	= 0.05;
	const double	tiny	= 2.5e-13;
	const double	big		= 9.0e11;
	const double	c1		= (3.0/14.0);
	const double	c2		= (1.0/3.0);
	const double	c3		= (3.0/22.0);
	const double	c4		= (3.0/26.0);
	const double	c5		= (0.75*c3);
	const double	c6		= (1.5*c4);
	const double	c7		= (0.5*c2);
	const double	c8		= (c3+c3);

	if (std::min(std::min(x,y),z) < 0.0 || std::min(std::min(std::min(x+y,x+z),y+z),fabs(p)) < tiny || std::max(std::max(std::max(x,y),z),fabs(p)) > big)
		ThrowError((char *)"invalid arguments in rj");
		
	double	sum = 0.0;
	double	fac = 1.0;
	
	double	xt, zt, yt, pt, a, b, rho, tau;
	
	if (p > 0.0) {
		xt	= x;
		yt	= y;
		zt	= z;
		pt	= p;
	} else {
		xt	= std::min(std::min(x,y),z);
		zt	= std::max(std::max(x,y),z);
		yt	= x+y+z-xt-zt;
		a	= 1.0/(yt-p);
		b	= a*(zt-yt)*(yt-xt);
		pt	= yt+b;
		rho	= xt*zt/yt;
		tau	= p*pt/yt;
	}
	
	double	delx, dely, delz, delp, ave;
	
	do {
		double	sqrtx	= sqrt(xt);
		double	sqrty	= sqrt(yt);
		double	sqrtz	= sqrt(zt);
		double	alamb	= sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		double	alpha	= sqr(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
		double	beta	= pt*sqr(pt+alamb);
		sum += fac * rc(alpha,beta);
		fac = 0.25*fac;
		xt	= 0.25*(xt+alamb);
		yt	= 0.25*(yt+alamb);
		zt	= 0.25*(zt+alamb);
		pt	= 0.25*(pt+alamb);
		ave		= 0.2*(xt+yt+zt+pt+pt);
		delx	= (ave-xt)/ave;
		dely	= (ave-yt)/ave;
		delz	= (ave-zt)/ave;
		delp	= (ave-pt)/ave;
	}
	while ( std::max(std::max(std::max(fabs(delx),fabs(dely)),fabs(delz)),fabs(delp)) > errtol );
		
	double	ea	= delx*(dely+delz)+dely*delz;
	double	eb	= delx*dely*delz;
	double	ec	= delp*delp;
	double	ed	= ea-3.0*ec;
	double	ee	= eb+2.0*delp*(ea-ec);
	double	ans	= 3.0*sum+fac*(1.0+ed*(-c1+c5*ed-c6*ee)+eb*(c7+delp*(-c8+delp*c4))
		+delp*ea*(c2-delp*c3)-c2*delp*ec)/(ave*sqrt(ave));
		
	if (p <= 0.0) ans = a * (b*ans + 3.0*(rc(rho,tau) - rf(xt,yt,zt)));
	
	return ans;
}

double	MyMath::rc( double x, double y )
{	
	const double	errtol	= 0.04;
	const double	tiny	= 1.69e-38;
	const double	sqrtny	= 1.3e-19;
	const double	big		= 3.e37;
	const double	tnbg	= tiny*big;
	const double	comp1	= (2.236/sqrtny);
	const double	comp2	= (tnbg*tnbg/25.0);
	const double	third	= (1.0/3.0);
	const double	c1		= 0.3;
	const double	c2		= (1.0/7.0);
	const double	c3		= 0.375;
	const double	c4		= (9.0/22.0);

	//float alamb,ave,s,w,xt,yt;
	
	if (x < 0.0 || y == 0.0 || (x+fabs(y)) < tiny || (x+fabs(y)) > big || (y<-comp1 && x > 0.0 && x < comp2))
		ThrowError((char *)"invalid arguments in rc");
		
	double	xt, yt, w;
	if (y > 0.0) {
		xt	= x;
		yt	= y;
		w	= 1.0;
	} else {
		xt	= x-y;
		yt	= -y;
		w	= sqrt(x)/sqrt(xt);
	}
	
	double	ave, s;
	do {
		double	alamb = 2.0*sqrt(xt)*sqrt(yt)+yt;
		xt	= 0.25*(xt+alamb);
		yt	= 0.25*(yt+alamb);
		ave	= third*(xt+yt+yt);
		s	= (yt-ave)/ave;
	}
	while (fabs(s) > errtol);
	
	return w*(1.0+s*s*(c1+s*(c2+s*(c3+s*c4))))/sqrt(ave);
}

#pragma mark -
#pragma mark //	linear algebra

void	MyMath::gaussj( double_matrix1& a, double_matrix1& b, size_t nn )
{
	size_t n = ( nn == 0 ? a.dim1() : ( nn > a.dim1() ? a.dim1() : nn ) );		//	= a.dim2() = b.dim1()
	size_t m = b.dim2();														//	nb of equations
	
	ulong_vector1	indxc(n);
	ulong_vector1	indxr(n);
	ulong_vector1	ipiv(n);
	
	for (size_t j=1; j<=n; j++)	ipiv[j] = 0;
	for (size_t i=1; i<=n; i++)
	{	size_t	icol,irow;
		double	big = 0.0;
		for (size_t j=1; j<=n; j++)
			if (ipiv[j] != 1)
				for (size_t k=1; k<=n; k++)
				{	if (ipiv[k] == 0)
					{	if (std::abs(a(j,k)) >= big)
						{	big = std::abs(a(j,k));
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1) ThrowError((char *)"gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol)
		{	for (size_t l=1; l<=n; l++)	swap(a(irow,l),a(icol,l));
			for (size_t l=1; l<=m; l++)	swap(b(irow,l),b(icol,l));
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a(icol,icol) == 0.0)	ThrowError((char *)"gaussj: Singular Matrix-2");
		
		double pivinv = 1.0/a(icol,icol);
		a(icol,icol)=1.0;
		for (size_t l=1; l<=n; l++)	a(icol,l) *= pivinv;
		for (size_t l=1; l<=m; l++)	b(icol,l) *= pivinv;
		for (size_t ll=1; ll<=n; ll++)
		{	if (ll != icol)
			{	double	dum = a(ll,icol);
				a(ll,icol) = 0.0;
				for (size_t l=1;l<=n;l++) a(ll,l) -= a(icol,l)*dum;
				for (size_t l=1;l<=m;l++) b(ll,l) -= b(icol,l)*dum;
		}	}
	}
	for (size_t l=n; l>=1; l--)
	{	if (indxr[l] != indxc[l])
			for (size_t k=1; k<=n; k++)
				swap( a(k,indxr[l]), a(k,indxc[l]) );
	}
}

void	MyMath::lubksb( const double_matrix1& a, const int_vector1& indx, double_vector1& b )
{	size_t	ii = 0;
	size_t	n = indx.length();
	for (size_t i=1; i<=n; i++)
	{	size_t		ip	= indx[i];
		double	sum	= b[ip];
		b[ip] = b[i];
		if (ii)	for (size_t j=ii; j<=i-1; j++) sum -= a(i,j)*b[j];
		else if (sum) ii=i;
		b[i] = sum;
	}
	for (size_t i=n; i>=1; i--)
	{	double	sum	= b[i];
		for (size_t j=i+1; j<=n; j++) sum -= a(i,j)*b[j];
		b[i] = sum/a(i,i);
	}
}

void	MyMath::ludcmp( double_matrix1& a, int_vector1& indx, Boolean& even )
{	
	unsigned long 	n = indx.length();
	const double	tiny = 1e-20;
	double_vector1	vv(n);

	even = true;
	for (unsigned long i=1; i<=n; i++)
	{	double	big = 0.0;
		double temp;
		for (unsigned long j=1; j<=n; j++)
			if ((temp=fabs(a(i,j))) > big)	big = temp;
		if (big == 0.0) ThrowError((char *)"Singular matrix in routine ludcmp");
		vv[i] = 1.0/big;
	}
	
    for (int j=1; j<=n; j++)
	{
		for (int i=1; i<j; i++)
		{	double	sum = a(i,j);
			for (int k=1; k<i; k++) sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
		}
		
		double	dum;
		int		imax;
		double	big = 0.0;
		for (unsigned long i=j; i<=n; i++)
		{	double	sum = a(i,j);
			for (int k=1; k<j; k++)	sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
			if ( (dum=vv[i]*fabs(sum)) >= big)
			{	big = dum;	imax = i;	}
		}
		
		if (j != imax)
		{	for (unsigned long k=1; k<=n; k++)	std::swap( a(imax,k), a(j,k) );
			even = !even;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		
		if (a(j,j) == 0.0)	a(j,j) = tiny;
		if (j != n)
		{	dum = 1.0/a(j,j);
			for (unsigned long i=j+1; i<=n; i++)	a(i,j) *= dum;
		}
	}
};


void	MyMath::tridag( const double_vector1& a, const double_vector1& b, const double_vector1& c,
			const double_vector1& r, double_vector1& u )
{								//	a sub-diag,	b diag,	c super-diag
	size_t	n = r.length();		//	a (1,n-1),	b (1,n),	c (1,n-1),	r (1,n),	u (1,n)
	double_vector1	gam(n-1);
	
	double			bet = b[1];
	if (bet == 0.0)	ThrowError((char *)"Bad diagonal in tridag");
	u[1] = r[1]/bet;
	
	for (size_t j=1; j<=n-1; j++)
	{	gam[j]= c[j]/bet;
		bet		= b[j+1] - a[j]*gam[j];
		if (bet == 0.0)		ThrowError((char *)"Bad diagonal in tridag");
		u[j+1]	= (r[j+1]-a[j]*u[j])/bet;
	}
	for (size_t j=(n-1); j>=1; j--)
		u[j] -= gam[j]*u[j+1];
}

void	MyMath::svbksb( const a_double_matrix1& u, const a_double_vector1& w, const a_double_matrix1& v,
		 	const a_double_vector1& b, a_double_vector1& x )
{	
	size_t	n = w.length();		//	= u.dim2() = v.dim1() = v.dim2() = x.length()
	size_t	m = b.length();		//	= u.dim1()
	
	double_vector1 tmp(n);
	for (size_t j=1; j<=n; j++)
	{	double	s = 0.0;
		if (w[j])
		{	for (size_t i=1; i<=m; i++)		s += u(i,j)*b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (size_t j=1; j<=n; j++)
	{	double	s = 0.0;
		for (size_t jj=1; jj<=n; jj++) 		s += v(j,jj)*tmp[jj];
		x[j] = s;
	}
};

void	MyMath::svdcmp( a_double_matrix1& a, a_double_vector1& w, a_double_matrix1& v )
{
	size_t	n = a.dim2();		//	= v.dim1() = v.dim2()
	size_t	m = a.dim1();
	
	double_vector1 rv1(n);
	double g 	 = 0;
	double scale = 0;
	double anorm = 0;
	
	double c,f,h,s,x,y,z;
	size_t l;
	
	for (size_t i=1; i<=n; i++)
	{	l = i+1;
		rv1[i]   = scale*g;
		g = s = scale = 0.0;
		
		//CheckRunningDialog( (double)i/(double)n/2 );
		if (i <= m)
		{	for (size_t k=i; k<=m; k++)	scale += std::abs(a(k,i));
			if (scale)
			{	for (size_t k=i; k<=m; k++)
				{	a(k,i) /= scale;	s += a(k,i)*a(k,i);	}
				f = a(i,i);
				g = -set_sign(sqrt(s),f);
				h = f*g - s;
				a(i,i) = f - g;
				for (size_t j=l; j<=n; j++)
				{	s = 0.0;
					for (size_t k=i; k<=m; k++)	s += a(k,i)*a(k,j);
					f = s/h;
					for (size_t k=i; k<=m; k++)	a(k,j) += f*a(k,i);
				}
				for (size_t k=i; k<=m; k++)	a(k,i) *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n)
		{	for (size_t k=l; k<=n; k++) scale += std::abs(a(i,k));
			if (scale)
			{	for (size_t k=l; k<=n; k++)
				{	a(i,k) /= scale;	s += a(i,k)*a(i,k);	}
				f = a(i,l);
				g = -set_sign(sqrt(s),f);
				h = f*g - s;
				a(i,l) = f - g;
				for (size_t k=l; k<=n; k++)	rv1[k] = a(i,k)/h;
				for (size_t j=l; j<=m; j++)
				{	s = 0.0;
					for (size_t k=l; k<=n; k++)	s += a(j,k)*a(i,k);
					for (size_t k=l; k<=n; k++) a(j,k) += s*rv1[k];
				}
				for (size_t k=l; k<=n; k++)	a(i,k) *= scale;
			}
		}
		anorm = std::max( anorm, (std::abs(w[i]) + std::abs(rv1[i])) );
	}
	
	for (size_t i=n; i>=1; i--)
	{	if (i < n)
		{	if (g)
			{	for (size_t j=l; j<=n; j++)		v(j,i) = (a(i,j)/a(i,l))/g;
				for (size_t j=l; j<=n; j++)
				{	s = 0.0;
					for (size_t k=l; k<=n; k++) s += a(i,k)*v(k,j);
					for (size_t k=l; k<=n; k++) v(k,j) += s*v(k,i);
				}
			}
			for (size_t j=l; j<=n; j++)		v(i,j) = v(j,i) = 0.0;
		}
		v(i,i) = 1.0;
		g = rv1[i];
		l = i;
	}
	
	for (size_t i=std::min(m,n); i>=1; i--)
	{	l = i+1;
		g = w[i];
		for (size_t j=l; j<=n; j++)		a(i,j) = 0.0;
		if (g)
		{	g = 1.0/g;
			for (size_t j=l; j<=n; j++)
			{	s = 0.0;
				for (size_t k=l; k<=m; k++)		s += a(k,i)*a(k,j);
				f = (s/a(i,i))*g;
				for (size_t k=i; k<=m; k++)		a(k,j) += f*a(k,i);
			}
			for (size_t j=i; j<=m; j++)		a(j,i) *= g;
		}
		else for (size_t j=i; j<=m; j++) 	a(j,i) = 0.0;
		++a(i,i);
	}
	
	for (size_t k=n; k>=1; k--)
	{	//CheckRunningDialog( 0.5 + (double)(n-k)/(double)n/2 );
		for (size_t its=1; its<=30; its++)
		{	bool	flag = 1;
			size_t	nm;
			for (l=k; l>=1; l--)
			{	nm = l-1;
				if ((double)(std::abs(rv1[l])+anorm) == anorm)	{	flag=0;		break;	}
				if ((double)(std::abs(w[nm])+anorm)  == anorm)		break;
				//if ((float)(abs(rv1[l])+anorm) == (float)anorm)	{	flag=0;		break;	}
				//if ((float)(abs(w[nm])+anorm)  == (float)anorm)		break;
			}
			if (flag)
			{	c = 0.0;
				s = 1.0;
				for (size_t i=l; i<=k; i++)
				{	f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((double)(std::abs(f)+anorm) == anorm) 	break;
					//if ((float)(abs(f)+anorm) == (float)anorm) 	break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (size_t j=1; j<=m; j++)
					{	y = a(j,nm);
						z = a(j,i);
						a(j,nm) = y*c + z*s;
						a(j,i)  = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k)
			{	if (z < 0.0)
				{	w[k] = -z;
					for (size_t j=1; j<=n; j++)		v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == 30)
				ThrowError((char *)"no convergence in 30 svdcmp iterations");
				
			x  = w[l];
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+set_sign(g,f)))-h))/x;
			c = s = 1.0;
			
			for (size_t j=l; j<=nm; j++)
			{	size_t i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for (size_t jj=1; jj<=n; jj++)
				{	x = v(jj,j);
					z = v(jj,i);
					v(jj,j) = x*c + z*s;
					v(jj,i) = z*c - x*s;
				}
				z = pythag(f,h);
				w[j] = z;
				if (z)
				{	z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for (size_t jj=1; jj<=m; jj++)
				{	y = a(jj,j);
					z = a(jj,i);
					a(jj,j) = y*c + z*s;
					a(jj,i) = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k]   = x;
		}
	}
};


void	MyMath::tred2( double_matrix1& a, double_vector1& d, double_vector1& e, Boolean eigV )
{
//	int l,k,j,i;
//	float scale,hh,h,g,f;
	
	size_t n = a.dim1();
	if (a.dim2() != n || d.length() != n || e.length() != n)	ThrowError((char *)"bad dimensions in tred2" );

	for (size_t i=n; i>=2; i--)
	{	size_t	l = i-1;
		double	h = 0;
		double	scale = 0;
		if (l > 1)
		{	for (size_t k=1; k<=l; k++)	scale += fabs(a(i,k));
			if (scale == 0.0)			e(i) = a(i,l);
			else
			{	for (size_t k=1; k<=l; k++)
				{	a(i,k) /= scale;
					h += a(i,k)*a(i,k);
				}
				double	f = a(i,l);
				double	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f*g;
				a(i,l) = f-g;
				f = 0.0;
				for (size_t j=1; j<=l; j++)
				{	a(j,i) = a(i,j)/h;
					g = 0.0;
					for (size_t k=1; k<=j; k++)	g += a(j,k)*a(i,k);
					for (size_t k=j+1;k<=l;k++)	g += a(k,j)*a(i,k);
					e[j] = g/h;
					f += e[j]* a(i,j);
				}
				double hh = f/(h+h);
				for (size_t j=1; j<=l; j++)
				{	f = a(i,j);
					e[j] = g = e[j] - hh*f;
					for (size_t k=1; k<=j; k++)
						a(j,k) -= ( f*e[k] + g*a(i,k) );
				}
			}
		}
		else
			e[i] = a(i,l);
		d[i] = h;
	}
	d[1] = 0.0;
	e[1] = 0.0;
	
	//	Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	if (eigV)
		for (size_t i=1; i<=n; i++)
		{	size_t l = i-1;
			if (d[i])
			{	for (size_t j=1; j<=l; j++)
				{	double g = 0;
					for (size_t k=1; k<=l; k++)		g += a(i,k)*a(k,j);
					for (size_t k=1; k<=l; k++)		a(k,j) -= g*a(k,i);
				}
			}
			d[i] = a(i,i);
			a(i,i) = 1.0;
			for (size_t j=1; j<=l; j++)	a(j,i) = a(i,j)=0.0;
		}
	else
		for (size_t i=1; i<=n; i++)
			d[i] = a(i,i);
};

void	MyMath::tqli( double_vector1& d, double_vector1& e, double_matrix1* z )
{
	size_t n = d.length();
	if (e.length() != n)							ThrowError((char *)"bad dimensions in tqli" );
	if ((z) && (z->dim2() != n || z->dim1() != n))	ThrowError((char *)"bad dimensions in tqli" );

	for (size_t i=2; i<=n; i++)		e[i-1] = e[i];
	e[n] = 0.0;
	for (size_t l=1; l<=n; l++)
	{	size_t iter = 0;
		size_t m;
		do
		{	for (m=l; m<=n-1; m++)
			{	float	dd = fabs(d[m]) + fabs(d[m+1]);
				if ( (fabs((float)(e[m]))+dd) == dd)	break;
			}
			if (m != l) {
				if (iter++ == 30)	ThrowError((char *)"Too many iterations in tqli");
				double	g = (d[l+1]-d[l])/(2.0*e[l]);
				double	r = pythag( g, 1.0 );
				g = d[m] - d[l] + e[l]/(g+set_sign(r,g));
				double	s = 1.0;
				double	c = 1.0;
				double	p = 0.0;
				size_t i;
				for (i=m-1; i>=l; i--)
				{	double	f = s*e[i];
					double	b = c*e[i];
					e[i+1] = (r = pythag(f,g));
					if (r == 0.0)
					{	d[i+1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d[i+1]-p;
					r = (d[i]-g)*s+2.0*c*b;
					d[i+1] = g + (p=s*r);
					g = c*r-b;
					if (z)
						for (size_t k=1; k<=n; k++)
						{	f = (*z)(k,i+1);
							(*z)(k,i+1) = s*(*z)(k,i) + c*f;
							(*z)(k,i)	= c*(*z)(k,i) - s*f;
						}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		}
		while (m != l);
	}
};


void	MyMath::eigsrt( double_vector1& d, double_matrix1& v )
{
	size_t n = v.dim1();
	if (v.dim2() != n || d.length() != n)	ThrowError((char *)"bad dimensions in eigsrt" );

	for (size_t i=1; i<n; i++)
	{	size_t	k = i;
		double	p = d[k];
		for (size_t j=i+1; j<=n; j++)
			if (d[j] >= p)	p = d[k=j];
		if (k != i)
		{	d[k] = d[i];
			d[i] = p;
			for (size_t j=1; j<=n; j++)
			{	p = v(j,i);		v(j,i) = v(j,k);	v(j,k) = p;	}
		}
	}
}


#pragma mark -
#pragma mark //	filters, convol and FFT


void	MyMath::savgol( double_vector1& c, int nl, int nr, int ld, int m )
{
	if (nl < 0 || nr < 0 || ld > m || nl+nr < m)
		ThrowError((char *)"bad args in savgol");
	int np = c.length();
	if (np < nl+nr+1)
	{	c.resize_but_keep(nl+nr+1);	np = c.length();	}
		
	int_vector1		indx(m+1);
	double_matrix1	a(m+1,m+1);
	double_vector1	b(m+1);
	
	for (int ipj=0; ipj<= (m<<1); ipj++)
	{	double	sum = (ipj ? 0.0 : 1.0);
		for (int k=1; k<=nr; k++)	sum += pow((double)k,(double)ipj);
		for (int k=1; k<=nl; k++)	sum += pow((double)-k,(double)ipj);
		int	mm = std::min( ipj, 2*m-ipj );
		for (int imj = -mm; imj<=mm; imj+=2)
			a( 1+(ipj+imj)/2, 1+(ipj-imj)/2 ) = sum;
	}
	
	Boolean	bb;
	ludcmp( a, indx, bb );
	
	for (int j=1; j<=m+1; j++)	b[j] = 0.0;
	b[ld+1] = 1.0;
	lubksb( a, indx, b );
	
	for (int kk=1; kk<=np; kk++)	c[kk] = 0.0;
	
	for (int k = -nl; k<=nr; k++)
	{	double	sum = b[1];
		double	fac = 1.0;
		for (int mm=1; mm<=m; mm++) sum += b[mm+1]*(fac *= k);
		c[((np-k) % np)+1] = sum;
	}
}

void	MyMath::four1( double_power2_array& data, Boolean forward )
{
	unsigned long	n = data.length();
	four1_n( data, n, forward );
};

void	MyMath::four1_n( a_double_vector1& data, size_t n, Boolean forward )
{	//	double wtemp,wr,wpr,wpi,wi,theta;	-> should be long_double ?
	//	float tempr,tempi;

//	unsigned long	n = data.length();		//	n = data.length() << 1;		//	= 2 N
	if ( n != nextPower2(n) )	ThrowError((char *)"Pb in four1_n");
	unsigned long	j = 1;
	for (unsigned long i=1; i<n; i+=2)
	{	if (j > i)
		{	std::swap( data[j], data[i] );
			std::swap( data[j+1], data[i+1] );
		}
		unsigned long	m = n >> 1;				//	= n/2 = N
		while (m >= 2 && j > m)
		{	j -= m;		m >>= 1;	}
		j += m;
	}

	unsigned long	mmax = 2;
	int	isign = (forward ? 1 : -1);
	while (n > mmax)
	{	unsigned long	istep = mmax << 1;
		double	theta = isign*(6.28318530717959/mmax);
		double	wtemp = sin(0.5*theta);
		double	wpr	  = -2.0*wtemp*wtemp;
		double	wpi   = sin(theta);
		double	wr    = 1.0;
		double	wi	  = 0.0;
		for (unsigned long m=1; m<mmax; m+=2)
		{	for (unsigned long i=m; i<=n; i+=istep)
			{	unsigned long j = i+mmax;
				double	tempr = wr*data[j] - wi*data[j+1];
				double	tempi = wr*data[j+1] + wi*data[j];
				data[j]   = data[i]-tempr;
				data[j+1] = data[i+1]-tempi;
				data[i]   += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp=wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
	}
};

void	MyMath::realft( double_power2_array& data, Boolean forward )
{
	unsigned long	n = data.length();
	realft_n( data, n, forward );
};

void	MyMath::realft_n( a_double_vector1& data, size_t n, Boolean forward )
{	
	// on output: data[1] = F(0), 			data[2] = F(n/2)
	//			  data[2k+1] = Re[F(k)],	data[2k+2] = Im[F(k)]	(k=1 a n/2 - 1)
	
	if ( n != nextPower2(n) )	ThrowError((char *)"Pb in realft_n");			//	n = puissance de 2, data.size() >= n
	double	theta = 3.141592653589793/(double) (n>>1);
	double	c1 = 0.5;
	double	c2;
	//	double wtemp,wr,wpr,wpi,wi,theta;	-> should be long_double ?
	
	if (forward)
	{	c2 = -0.5;	four1_n( data, n, true );	}
	else
	{	c2 = 0.5;	theta = -theta;		}
	
	double wtemp = sin(0.5*theta);
	double wpr   = -2.0*wtemp*wtemp;
	double wpi   = sin(theta);
	double wr    = 1.0+wpr;
	double wi    = wpi;
	unsigned long np3 = n + 3;
	
	for (unsigned long i=2; i<=(n>>2); i++)
	{	unsigned long i1 = i + i - 1;
		unsigned long i2 = 1 + i1;
		unsigned long i3 = np3 - i2;
		unsigned long i4 = 1 + i3;
		double	h1r = c1 *(data[i1]+data[i3]);
		double	h1i = c1 *(data[i2]-data[i4]);
		double	h2r = -c2*(data[i2]+data[i4]);
		double	h2i = c2 *(data[i1]-data[i3]);
		data[i1] =  h1r + wr*h2r -wi*h2i;
		data[i2] =  h1i + wr*h2i +wi*h2r;
		data[i3] =  h1r - wr*h2r +wi*h2i;
		data[i4] = -h1i + wr*h2i +wi*h2r;
		wtemp = wr;
		wr = wtemp*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
	}
	if (forward)
	{	double	h1r = data[1];
		data[1] += data[2];
		data[2] = h1r - data[2];
	}
	else
	{	double	h1r = data[1];
		data[1] = c1*( data[1] + data[2] );
		data[2] = c1*( h1r - data[2] );
		four1_n( data, n, false );
	};
}


void	MyMath::cosft1( a_double_power2plus1_array& y )
{
	size_t	n  = y.length() - 1;			//	n = power of 2
	size_t	n2 = n+2;
	
	double	wi	= 0.0;
	double	wr	= 1.0;
	double	theta	= 3.141592653589793/ (double)n;
	double	wtemp	= sin(0.5*theta);
	double	wpr 	= -2.0*wtemp*wtemp;
	double	wpi		= sin(theta);
	double	sum		= 0.5*(y[1]-y[n+1]);
	
	y[1] = 0.5*(y[1]+y[n+1]);
	
	for (size_t j=2; j<=(n>>1); j++)
	{	wr = (wtemp=wr)*wpr-wi*wpi+wr;
		wi = wi*wpr+wtemp*wpi+wi;
		double	y1 = 0.5*(y[j]+y[n2-j]);
		double	y2 = (y[j]-y[n2-j]);
		y[j] = y1-wi*y2;
		y[n2-j] = y1+wi*y2;
		sum += wr*y2;
	}
//	double_power2_array y1( y );
//	realft( y1 );
	realft_n( y, n );
	y[n+1] = y[2];
	y[2]   = sum;
	for (size_t j=4; j<=n; j+=2)
	{	sum += y[j];
		y[j] = sum;
	}
}

void	MyMath::cosft2( double_power2_array& y, Boolean forward )
{
	size_t	n  = y.length();			//	n = power of 2
	
	double theta = 0.5 * 3.14159265358979/n;
	double wi	 = 0.0;
	double wr	 = 1.0;
	double wr1   = cos(theta);
	double wi1   = sin(theta);
	double wpr   = -2.0*wi1*wi1;
	double wpi   = sin(2.0*theta);
	
	if (forward)
	{	for (int i=1; i<=n/2; i++)
		{	double y1 = 0.5*(y[i]+y[n-i+1]);
			double y2 = wi1*(y[i]-y[n-i+1]);
			y[i]     = y1 + y2;
			y[n-i+1] = y1 - y2;
			double wtemp = wr1;
			wr1 = wtemp*wpr - wi1*wpi + wr1;
			wi1 = wi1*wpr + wtemp*wpi + wi1;
		}
		realft( y, forward );
		for (int i=3; i<=n; i+=2)
		{	double wtemp = wr;
			wr = wtemp*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
			double y1 = y[i]*wr-y[i+1]*wi;
			double y2 = y[i+1]*wr+y[i]*wi;
			y[i]   = y1;
			y[i+1] = y2;
		}
		double s = 0.5*y[2];
		for (int i=n; i>=2; i-=2)
		{ 	double sum1 = s;
			s += y[i];
			y[i] = sum1;
	}	}
	else							//	backward
	{	double ytemp = y[n];
		for (int i=n; i>=4; i-=2)
			y[i] = y[i-2] - y[i];
		y[2] = 2.0*ytemp;
		for (int i=3; i<=n; i+=2)
		{	double wtemp = wr;
			wr = wtemp*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
			double y1 = y[i]*wr + y[i+1]*wi;
			double y2 = y[i+1]*wr - y[i]*wi;
			y[i]   = y1;
			y[i+1] = y2;
		}
		realft( y, forward );
		for (int i=1; i<=n/2; i++)
		{	double y1 = y[i]+y[n-i+1];
			double y2 = (0.5/wi1)*(y[i]-y[n-i+1]);
			y[i]     = 0.5*(y1+y2);
			y[n-i+1] = 0.5*(y1-y2);
			double wtemp = wr1;
			wr1 = wtemp*wpr - wi1*wpi + wr1;
			wi1 = wi1*wpr + wtemp*wpi + wi1;
		}
	}
}

void	MyMath::sinft( double_power2_array& y )
{	size_t	n  = y.length();
	size_t	n2 = n + 2;

	double	wi	  = 0.0;
	double	wr	  = 1.0;
	double	theta = 3.14159265358979/(double)n;
	double	wtemp = sin(0.5*theta);
	double	wpr   = -2.0*wtemp*wtemp;
	double	wpi   = sin(theta);
	y[1] = 0.0;
	for (size_t j=2; j<=(n>>1)+1; j++)
	{	wr = (wtemp=wr)*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
		double	y1 = wi*(y[j]+y[n2-j]);
		double	y2 = 0.5*(y[j]-y[n2-j]);
		y[j]	= y1+y2;
		y[n2-j]	= y1-y2;
	}
	realft( y );
	y[1] *= 0.5;
	double	sum = y[2] = 0.0;
	for (size_t j=1; j<=n-1; j+=2)
	{	sum    += y[j];
		y[j]	= y[j+1];
		y[j+1] 	= sum;
	}
}

void	MyMath::convlv_fft( const a_double_vector1& data, double_filter_vector& respns, double_vector1& answer, Boolean forward )
{
	double_power2_array			Data_fft( data, respns.zero_pad_length() );		//	double_power2_array::no );
	unsigned long				n = Data_fft.length();
	const double_power2_array&	Filter_fft = respns.getFFTFilter( n );
	
	realft( Data_fft );
	double_power2_array		answer_fft(n);
	
	unsigned long	no2 = n>>1;				//	= n/2
	for (unsigned long i=4; i<=n; i+=2)
	{	if (forward)
		{	answer_fft[i-1]	= (Data_fft[i-1]*Filter_fft[i-1] - Data_fft[i]*Filter_fft[i]) / no2;
			answer_fft[i]	= (Data_fft[i]*Filter_fft[i-1] + Data_fft[i-1]*Filter_fft[i]) / no2;
		}
		else
		{	double	mag2 = sqr(Filter_fft[i-1]) + sqr(Filter_fft[i]);
			if (mag2 == 0.0)	ThrowError((char *)"Deconvolving at response zero in convlv" );
			answer_fft[i-1]	= (Data_fft[i-1]*Filter_fft[i-1] + Data_fft[i]*Filter_fft[i])/mag2/no2;
			answer_fft[i]	= (Data_fft[i]*Filter_fft[i-1] - Data_fft[i-1]*Filter_fft[i])/mag2/no2;
		}
	}

	if (forward)
	{	answer_fft[1] = Data_fft[1] * Filter_fft[1] / no2;
		answer_fft[2] = Data_fft[2] * Filter_fft[2] / no2;
	}
	else
	{	if (Filter_fft[1] == 0.0)	ThrowError((char *)"Deconvolving at response zero in convlv" );
		answer_fft[1] = Data_fft[1] / Filter_fft[1] / no2;
		if (Filter_fft[2] == 0.0)	ThrowError((char *)"Deconvolving at response zero in convlv" );
		answer_fft[2] = Data_fft[2] / Filter_fft[2] / no2;
	}
	
	realft( answer_fft, false );
	if (answer.length() != data.length())	answer.resize_but_keep(data.length());
	std::copy( answer_fft.firstPtr(), answer_fft.firstPtr()+data.length(), answer.firstPtr() );
}

void	MyMath::convlv_savgol( const double_vector& data, double_savgol_filter_matrix& respns, double_vector& answer,
			Boolean useFFT, Boolean remove_mean, double value )
{	int np = respns.length();				//	= nlm+nrm+1
	int nl = respns[1].l_length();			//	nl = nlm	... a priori. sinon, ThrowError
	int nr = respns[1].r_length();			//	nr = nrm
	if ( (np > 1) && (np != nl+nr+1) )	ThrowError((char *)"Bad filter length in MyMath::convlv_savgol" );
	
	double_vector*	px = (double_vector*)&data;		//	x is const anyway
	double_vector	z;
	if (np==1 && remove_mean)	value = data.mean();
	if (np==1 && (value != 0))	{	z = data;	z -= value;		px = &z;	}
	const double_vector&	x = *px;

	if (useFFT)	MyMath::convlv_fft( x, respns[1], answer, true );
	else		MyMath::convlv( x, respns[1], answer, false );
	
	if (np>1)
	{	for (int k=-nl; k<0; k++)
		{	int n 		= 1 + nl + k;
			double&	a 	= answer[n];	a = 0;
			double_filter_vector&	f = respns[1 + (np-k) % np];
			for (int j=-f.l_length(); j<=f.r_length(); j++)	a += f.value(j) * x[n+j];
		}
		for (int k=nr; k>0; k--)
		{	int n 		= x.length() - nr + k;
			double&	a 	= answer[n];	a = 0;
			double_filter_vector&	f = respns[1 + (np-k) % np];
			for (int j=-f.l_length(); j<=f.r_length(); j++)	a += f.value(j) * x[n+j];
		}
	}
	if (np==1 && (value != 0))	answer += value;
};

void	MyMath::convlv( const a_double_vector1& data, double_filter_vector& respns, double_vector& answer, Boolean renormalize_bounds )
{	if (answer.length() != data.length())
		answer.resize_but_keep(data.length());
	double	c_sum = 0;
	if (renormalize_bounds)
		for (int j=-respns.l_length(); j<=respns.r_length(); j++)
			c_sum += respns.value(j);
	if (c_sum != 0)
	{	for (unsigned long i=1; i<=data.length(); i++)
		{	double	p_sum = 0;
			answer[i] = 0;
			for (int j=-respns.l_length(); j<=respns.r_length(); j++)
				if ( i+j >= 1 && i+j <= data.length())
				{	answer[i] += respns.value(j) * data[i+j];	p_sum += respns.value(j);	}
			if (p_sum != c_sum && p_sum != 0)
				answer[i] *= c_sum/p_sum;
	}	}
	else
	{	for (unsigned long i=1; i<=data.length(); i++)
		{	answer[i] = 0;	
			for (int j=-respns.l_length(); j<=respns.r_length(); j++)
				if ( i+j >= 1 && i+j <= data.length())
					answer[i] += respns.value(j) * data[i+j];
	}	}
}	


void	MyMath::correl( const a_double_vector1& data1, const a_double_vector1& data2, double_filter_vector& answer, double norm, Boolean remove_mean )
{	//if (answer.l_length() > data1.length() || answer.r_length() > data2.length())
	//	ThrowError((char *)"Bad filter length in MyMath::correl" );
	unsigned long nl = answer.l_length();		//	<= data1.length()-1
	unsigned long nr = answer.r_length();		//	<= data2.length()-1
//	unsigned long m0 = std::min( data1.length(), data2.length() );
	if (remove_mean)
	{	double mean1 = data1.mean();
		double mean2 = data2.mean();
		
		for (unsigned long i=0; i<=nl; i++)		//	answer - left
		{	double& a = answer.value(-i);
			unsigned long m = std::min( data1.length()-i, data2.length() );
			for (unsigned long j=1; j<=m; j++)	a += (data1[j+i]-mean1) * (data2[j]-mean2);
			if (norm==0)	a /= m;
			else			a /= norm;
		}
		for (unsigned long i=1; i<=nr; i++)		//	answer - right
		{	double& a = answer.value(i);
			unsigned long m = std::min( data1.length(), data2.length()-i );
			for (unsigned long j=1; j<=m; j++)	a += (data1[j]-mean1) * (data2[j+i]-mean2);
			if (norm==0)	a /= m;
			else			a /= norm;
		}
	}
	else
	{	for (unsigned long i=0; i<=nl; i++)		//	answer - left
		{	double& a = answer.value(-i);
			unsigned long m = std::min( data1.length()-i, data2.length() );
			for (unsigned long j=1; j<=m; j++)	a += data1[j+i] * data2[j];
			if (norm==0)	a /= m;
			else			a /= norm;
		}
		for (unsigned long i=1; i<=nr; i++)		//	answer - right
		{	double& a = answer.value(i);
			unsigned long m = std::min( data1.length(), data2.length()-i );
			for (unsigned long j=1; j<=m; j++)	a += data1[j] * data2[j+i];
			if (norm==0)	a /= m;
			else			a /= norm;
		}
	}
}

void	MyMath::autocorrel( const a_double_vector1& data, double_filter_vector& answer, double norm, Boolean remove_mean )
{	if (answer.length() > data.length() || answer.length() == 0)
		answer.resize_but_keep( data.length() );
	correl( data, data, answer, norm, remove_mean );
}

void	MyMath::autocorrel_fft( const a_double_vector1& data, double_filter_vector& answer, double norm, Boolean remove_mean )
{	if (answer.length() > data.length() || answer.length() == 0)
		answer.resize_but_keep( data.length() );
	correl_fft( data, data, answer, norm, remove_mean );
}

void	MyMath::correl_fft( const a_double_vector1& data1, const a_double_vector1& data2, double_filter_vector& answer, double norm, Boolean removeMean )
{	
	unsigned long			n = nextPower2( data1.length() + data2.length() - 1 );
		//	longueur totale de la fonction de correlation
		//		= data1.length() + data2.length() - 1
	double_power2_array		Data1_fft( data1, n - data1.length(), (removeMean ? double_power2_array::remove_mean : double_power2_array::no) );
	double_power2_array		Data2_fft( data2, n - data2.length(), (removeMean ? double_power2_array::remove_mean : double_power2_array::no) );
	
	realft( Data1_fft );
	realft( Data2_fft );
	double_power2_array		answer_fft(n);
	
	unsigned long	no2 = n>>1;				//	= n/2
	for (unsigned long i=4; i<=n; i+=2)
	{	answer_fft[i-1]	= (Data1_fft[i-1]*Data2_fft[i-1] + Data1_fft[i]*Data2_fft[i]) / no2;
		answer_fft[i]	= (Data1_fft[i]*Data2_fft[i-1] - Data1_fft[i-1]*Data2_fft[i]) / no2;
	}
	answer_fft[1] = Data1_fft[1] * Data2_fft[1] / no2;
	answer_fft[2] = Data1_fft[2] * Data2_fft[2] / no2;
	
	realft( answer_fft, false );
	
	for (unsigned long i=0; i<=answer.l_length(); i++)
	{	if (norm==0)	answer.value(-i) = answer_fft[i+1] / std::min( data1.length()-i, data2.length() );
		else			answer.value(-i) = answer_fft[i+1] / norm;
	}
	for (unsigned long i=1; i<=answer.r_length(); i++)
	{	if (norm==0)	answer.value(i) = answer_fft[n-i+1] / std::min( data1.length(), data2.length()-i );
		else			answer.value(i) = answer_fft[n-i+1] / norm;
	}
}



void	MyMath::memcof( const a_double_vector1& data, double& xms, a_double_vector1& d )
{
	size_t	n = data.length();
	size_t	m = d.length();
	
	double	p = 0.0;
	double_vector	wk1(n);
	double_vector	wk2(n);
	double_vector	wkm(m);
	
	for (size_t j=1; j<=n; j++)		p += sqr(data[j]);
	xms = p/n;
	
	wk1[1]	 = data[1];
	wk2[n-1] = data[n];
	for (size_t j=2; j<=n-1; j++)
	{	wk1[j]=data[j];		wk2[j-1]=data[j];	}
	
	for (size_t k=1; k<=m; k++)
	{	double	num		= 0.0;
		double	denom	= 0.0;
		for (size_t j=1; j<=(n-k); j++)
		{	num	  += wk1[j] * wk2[j];
			denom += sqr(wk1[j])+sqr(wk2[j]);
		}
		d[k] = 2.0*num/denom;
		xms *= (1.0-sqr(d[k]));
		for (size_t i=1; i<=(k-1); i++)		d[i] = wkm[i] - d[k]*wkm[k-i];
		if (k == m)
			return;
		for (size_t i=1; i<=k; i++)			wkm[i]=d[i];
		for (size_t j=1; j<=(n-k-1); j++)
		{	wk1[j] -= wkm[k]*wk2[j];
			wk2[j]  = wk2[j+1] - wkm[k]*wk1[j+1];
		}
	}
	ThrowError((char *)"never get here in memcof." );
}

double	MyMath::evlmem( double fdt,  const a_double_vector1& d, double xms )
{	double	sumr = 1.0;
	double	sumi = 0.0;
	double	wr = 1.0;
	double	wi = 0.0;
	
	size_t	m = d.length();

	const double	theta = 6.28318530717959 * fdt;
	const double	wpr = cos(theta);
	const double	wpi = sin(theta);
	
	for (size_t i=1; i<=m; i++)
	{	double	wtemp;
		wr = (wtemp=wr)*wpr - wi*wpi;
		wi = wi*wpr + wtemp*wpi;
		sumr -= d[i]*wr;
		sumi -= d[i]*wi;
	}
	return xms/(sumr*sumr + sumi*sumi);
}



#pragma mark -
#pragma mark //	fitting



void	MyMath::svdfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, a_double_matrix1& u, a_double_matrix1& v, a_double_vector1& w,
			double& chisq, fit_func* ff )
{
	size_t ndata = x.length();		//	= y.length() = sig.length() = u.dim1()
	size_t ma 	 = a.length();		//	= u.dim2() = v.dim1() = v.dim2() = w.length()
	
	double_vector1 b(ndata);
	double_vector1 afunc(ma);
	
	for (size_t i=1; i<=ndata; i++)
	{	ff->fitF( x[i], afunc );
		double tmp = 1.0/sig[i];
		for (size_t j=1; j<=ma; j++)	u(i,j) = afunc[j]*tmp;
		b[i] = y[i]*tmp;
	}

	svdcmp( u, w, v );

	double wmax = 0.0;
	for (size_t j=1; j<=ma; j++)	if (w[j] > wmax)	wmax = w[j];
		
	//double thresh = svdFitTol*wmax;
	double thresh = 1e-30*wmax;
	for (size_t j=1; j<=ma; j++)	if (w[j] < thresh)	w[j] = 0.0;
	
	svbksb( u, w, v, b, a );

	chisq = 0.0;
	for (size_t i=1; i<=ndata; i++)
	{	ff->fitF( x[i], afunc );
		double sum = 0;
		for (size_t j=1; j<=ma; j++)	sum += a[j]*afunc[j];
		chisq += sqr( (y[i]-sum)/sig[i] );
	}
};

#if 0

//	implementation "ammelioree", qui ne fonctionne pas...

void	MyMath::absfit( const double_vector1& x, const double_vector1& y, const double_vector1& sig,
			double_vector1& a, double& absv, fit_func* ff )
{
	size_t ndat	= x.length();		//	= y.length() = sig.length()
	size_t ma	= a.length();		//
	
	size_t ncon	= ndat - ma;
	size_t nvar	= ndat + ma;
	double_matrix1	aa( ncon+2, nvar+1 );
	double_matrix1	fx( ma, ndat );
	
	for (size_t i=1; i<=ndat; i++)
	{	double_vector1	afunc( ma );
		ff->fitF( x[i], afunc );
		for (size_t j=1; j<=ma; j++)	fx(j,i) = afunc(j);
	};
	
	double_matrix1	lbda( ma, ma );
	double_matrix1	lbdb( ma, 1 );

	for (size_t i=1; i<=ma; i++)		//	or other subset?
	{	for (size_t j=1; j<=ma; j++)	lbda(i,j) = fx(j,i);
		lbdb(i,1) = y[i];
	}
	gaussj( lbda, lbdb );

	const size_t 	Roffset = 1;
	const size_t 	Soffset = Roffset + ma;
	aa(1,1) = 0;					//	Maximize -Sum{(r[k]+s[k]), 1,M} - Sum{A[k](r[k]-s[k]), 1,M} - 2 Sum{s[k], M+1,N}
	for (size_t k=1; k<=ma; k++)
	{	double	Ak = 0;
//		for (size_t i=1; i<=ndat; i++)
		for (size_t i=ma+1; i<=ndat; i++)
		for (size_t j=1; j<=ma; j++)	Ak += lbda(j,k)*sig[k] * fx(j,i)/sig[i];
		aa(1,Roffset+k) = -1.0 - Ak;		//	r[k]
		aa(1,Soffset+k) = -1.0 + Ak;		//	s[k]
	}
	for (size_t k=ma+1; k<=ndat; k++)
		aa(1,Soffset+k) = -2.0;			//	s[k]

	size_t	nconp = 0;
	size_t	nconq = 0;
	for (size_t i=ma+1; i<=ndat; i++)		//	constraints: - s[i] - Sum{ Aik (r[k]-s[k]), 1,M}  <= Ci
	{	double	Ci = y[i]/sig[i];
		for (size_t j=1; j<=ma; j++)	Ci -= lbdb(j,1)*fx(j,i)/sig[i];
		if (Ci >= 0)	nconp++;
		else			nconq++;
		size_t	index = (Ci >= 0 ? 1+nconp : 2+ncon-nconq);
		if (Ci >= 0)	aa(index,1) = Ci;
		else			aa(index,1) = -Ci;
		
		for (size_t k=1; k<=ma; k++)
		{	double	Aik = 0;
			for (size_t j=1; j<=ma; j++)	Aik += lbda(j,k)*sig[k] * fx(j,i)/sig[i];
		//	if (Ci >= 0)	{	aa(index,Roffset+k) = -Aik;		aa(index,Soffset+k) = Aik;	}
		//	else			{	aa(index,Roffset+k) = Aik;		aa(index,Soffset+k) = -Aik;	}
			if (Ci >= 0)	{	aa(index,Roffset+k) = Aik;		aa(index,Soffset+k) = -Aik;	}
			else			{	aa(index,Roffset+k) = -Aik;		aa(index,Soffset+k) = Aik;	}
		}
		for (size_t k=ma+1; k<=ndat; k++)
		{	if (k==i)
		//	{	if (Ci >= 0)	aa(index,Soffset+k) = -1.0;
		//		else			aa(index,Soffset+k) = 1.0;
			{	if (Ci >= 0)	aa(index,Soffset+k) = 1.0;
				else			aa(index,Soffset+k) = -1.0;
			}
			else aa(index,Soffset+k) = 0;
		}
	}
	
	int	icase;
	ulong_vector1	posv( ncon );
	ulong_vector1	zerov( nvar );
	Boolean minimumFound = simplx( aa, ncon, nvar, nconp, nconq, 0, icase, zerov, posv );
	ThrowIfNot_( minimumFound );
	
	double_vector1	v( nvar );									//	r[1]..r[M] s[1]..s[N]
	for (size_t i=1; i<=nvar; i++)		v[i] = 0;
	for (size_t i=1; i<=ncon; i++)
		if (posv[i] <= nvar )	v[ posv[i] ] = aa(i+1,1);
	
	for (size_t i=1; i<=ma; i++)
	{	a[i] = 0;
		for (size_t k=1; k<=ma; k++)	a[i] += lbda(i,k)*(y[k]-sig[k]*(v[k]-v[ma+k]));
	}
	absv = aa(1,1);
};

#else

//	implementation basique, qui marche (mais c'est moins performant)

void	MyMath::absfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, double& absv, fit_func* ff )
{
	size_t ndat	= x.length();		//	= y.length() = sig.length()
	size_t ma	= a.length();		//
	
	size_t ncon	= ndat;
	size_t nvar	= 2*ndat + 2*ma;
	double_matrix1	aa( ncon+2, nvar+1 );
	double_matrix1	fx( ma, ndat );
	
	for (size_t i=1; i<=ndat; i++)
	{	double_vector1	afunc( ma );
		ff->fitF( x[i], afunc );
		for (size_t j=1; j<=ma; j++)	fx(j,i) = afunc(j);
	};

	const size_t 	Roffset = 1;					//	start index for r[i],   i=1,ndat
	const size_t 	Soffset = Roffset + ndat;		//	start index for s[i],   i=1,ndat
	const size_t 	Aoffset = Roffset + 2*ndat;		//	start index for c+[k],  k=1,ma
	const size_t 	AMoffset = Aoffset + ma;		//	start index for c-[k],  k=1,ma
	
	aa(1,1) = 0;					//	Maximize -Sum{ (r[k]+s[k]), 1,N }
	for (size_t k=1; k<=ndat; k++)
	{	aa(1,Roffset+k) = -1.0;		//	r[k]
		aa(1,Soffset+k) = -1.0;		//	s[k]
	}
	for (size_t k=1; k<=ma; k++)
	{	aa(1,Aoffset+k) = 0.0;		//	c+[k]
		aa(1,AMoffset+k) = 0.0;		//	c-[k]
	}
	
	for (size_t i=1; i<=ndat; i++)
	{	aa(1+i,1) = y[i]/sig[i];
		for (size_t k=1; k<=ndat; k++)
		{	if (i==k)
			{	aa(1+i,Roffset+k) = -1.0;
				aa(1+i,Soffset+k) = 1.0;
			}
			else
			{	aa(1+i,Roffset+k) = 0.0;
				aa(1+i,Soffset+k) = 0.0;
		}	}
		for (size_t k=1; k<=ma; k++)
		{	aa(1+i,Aoffset+k) = -fx(k,i)/sig[i];
			aa(1+i,AMoffset+k) = fx(k,i)/sig[i];
		}
	}

	int	icase;
	ulong_vector1	posv( ncon );
	ulong_vector1	zerov( nvar );
	Boolean minimumFound = simplx( aa, ncon, nvar, 0, 0, ndat, icase, zerov, posv );
	ThrowErrorIfNot( minimumFound );
	
	double_vector1	v( nvar );									//	r[i] s[i] c+[k] c-[k]
	for (size_t i=1; i<=nvar; i++)		v[i] = 0;
	for (size_t i=1; i<=ncon; i++)
		if (posv[i] <= nvar )	v[ posv[i] ] = aa(i+1,1);

	for (size_t i=1; i<=ma; i++)
	{	a[i] = v[2*ndat+i]-v[2*ndat+ma+i];				//	a[k] = (c+[k]) - (c-[k])
	}
	absv = -aa(1,1);
};

#endif

void	MyMath::lfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, const int_vector1& ia,
			double_matrix1& covar, double& chisq, fit_func* ff )
{
	size_t ndat	= x.length();		//	= y.length() = sig.length()
	size_t ma	= a.length();		//	= ia.length() = covar.dim1() = covar.dim2()

	double_matrix1	beta(ma,1);
	double_vector1	afunc(ma);
	
	size_t	mfit = 0;
	for (size_t j=1; j<=ma; j++)
		if (ia[j])	mfit++;
	if (mfit == 0) ThrowError((char *)"lfit: no parameters to be fitted");
	
	for (size_t j=1; j<=mfit; j++)
	{	for (size_t k=1; k<=mfit; k++) covar(j,k) = 0.0;
		beta(j,1) = 0.0;
	}
	for (size_t i=1; i<=ndat; i++)
	{	ff->fitF( x[i], afunc );
		double ym = y[i];
		if (mfit < ma)
		{	for (size_t j=1; j<=ma; j++)
				if (!ia[j])
					ym -= a[j]*afunc[j];
		}
		double sig2i = 1/sqr(sig[i]);
		size_t j = 0;
		for (size_t l=1; l<=ma; l++)
		{	if (ia[l])
			{	double wt = afunc[l]*sig2i;
				j++;
				size_t k = 0;
				for (size_t m=1; m<=l; m++)
					if (ia[m]) covar(j,++k) += wt*afunc[m];
				beta(j,1) += ym*wt;
			}
		}
	}
	for (size_t j=2; j<=mfit; j++)
		for (size_t k=1; k<j; k++)
			covar(k,j) = covar(j,k);
			
	gaussj( covar, beta, mfit );
	
	size_t j = 0;
	for (size_t l=1; l<=ma; l++)
		if (ia[l]) a[l] = beta(++j,1);
		
	chisq = 0.0;
	for (size_t i=1; i<=ndat; i++)
	{	ff->fitF( x[i], afunc );
		double	sum = 0;
		for (size_t j=1; j<=ma; j++)	sum += a[j]*afunc[j];
		chisq += sqr( (y[i]-sum)/sig[i] );
	}
	
//	NR_covsrt(covar,ma,ia,mfit);
}



void	MyMath::lin_interpol_fit_func::fitF( double xx, double_vector1& f )
{	size_t nb = f.length();
	for (unsigned long i=1; i<=nb; i++)	f[i] = 0;
	double h = bissectbracket( x, xx, klo );
	double x1 = x[ klo ];
	double x2 = x[ klo+1 ];
	f[ klo ]   = (x2 - xx) / h;
	f[ klo+1 ] = (xx - x1) / h;
};

void	MyMath::stair_interpol_fit_func::fitF( double xx, double_vector1& f )
{	size_t nb = f.length();
	for (unsigned long i=1; i<=nb; i++)	f[i] = 0;
	hunt( x_midPoints(), xx, klo );
	f[ klo+1 ] = 1;
};

void	MyMath::polynom_fit_func::fitF( double xx, double_vector1& f )
{	size_t nb = f.length();
	f[1] = 1.0;
	for (unsigned long i=2; i<=nb; i++)	f[i] = f[i-1]*xx;
};

void	MyMath::splin_interpol_fit_func::fitF( double xx, double_vector1& f )
{	size_t nb = f.length();
	for (unsigned long i=1; i<=nb; i++)
	//	f[i] = spline_interp( x, extY[i], basisDersec[i], xx );
	{	double	h = bissectbracket( x, xx, klo );
		ThrowErrorIf (h == 0);
		size_t	khi = klo+1;
		double	a = (x[khi] - xx) / h;
		double	b = (xx - x[klo]) / h;
		f[i] = ( a*(sqr(a)-1) * ((*basisDersecPtr)[i])[klo]
					+ b * (sqr(b) - 1) * ((*basisDersecPtr)[i])[khi] ) * sqr(h) / 6;
		if (i==klo)		f[i] += a;
		if (i==khi)		f[i] += b;
	}
};



void	MyMath::line_fit_func::fit( const double_vector1& x, const double_vector1& y, const double_vector1& sig, Boolean hasSig, double& a, double& b )
{
	unsigned long		ndata = x.length();
	Boolean	has_sig = hasSig && (sig.length() == ndata);
	
	b = 0.0;
	double	ss;
	double	sx = 0;
	double	sy = 0;
	
	if (has_sig)
	{	ss = 0.0;
		for (size_t i=1; i<=ndata; i++)
		{	double	wt = 1.0/sqr(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	}
	else
	{	for (size_t i=1;i<=ndata;i++)
		{	sx += x[i];
			sy += y[i];
		}
		ss = ndata;
	}
	double	sxoss = sx/ss;
	double	st2 = 0;
	
	if (has_sig)
	{	for (size_t i=1; i<=ndata; i++)
		{	double	t = (x[i]-sxoss)/sig[i];
			st2 += t*t;
			b   += t*y[i]/sig[i];
		}
	}
	else
	{	for (size_t i=1; i<=ndata; i++)
		{	double	t = x[i]-sxoss;
			st2 += t*t;
			b   += t*y[i];
		}
	}
	
	b   /= st2;
	a	 = (sy-sx*b)/ss;
	siga = sqrt( (1.0 + sx*sx/(ss*st2)) / ss );
	sigb = sqrt( 1.0/st2 );
	chisq = 0.0;
	
	if (!has_sig)
	{	for (size_t i=1; i<=ndata; i++)		chisq += sqr( y[i] - a - b*x[i] );
		q = 1.0;
		double	sigdat = sqrt( chisq/(ndata-2) );
		siga *= sigdat;
		sigb *= sigdat;
	}
	else
	{	for (size_t i=1; i<=ndata; i++)		chisq += sqr( (y[i] - a - b*x[i])/sig[i] );
		q = gammq( 0.5 * (ndata-2), 0.5 * chisq );
	}
};

//	globals for routine chixy
double_vector1	MyMath::line_fit_func::chixy_xx;
double_vector1	MyMath::line_fit_func::chixy_yy;
double_vector1	MyMath::line_fit_func::chixy_sx;
double_vector1	MyMath::line_fit_func::chixy_sy;
double_vector1	MyMath::line_fit_func::chixy_ww;
double			MyMath::line_fit_func::chixy_aa;
double			MyMath::line_fit_func::chixy_offs;


double	MyMath::line_fit_func::chixy( double b_ang )
{
	const double	Big = 1.0e30;

	double	b = tan( b_ang );
	double	avex = 0;
	double	avey = 0;
	double	sumw = 0;
	size_t	nn = chixy_xx.length();
	
	for (unsigned long j=1; j<=nn; j++)
	{	chixy_ww[j] = sqr( b * chixy_sx[j] ) + sqr( chixy_sy[j] );
		chixy_ww[j] = ( chixy_ww[j] == 0.0 ? Big : 1.0 / chixy_ww[j] );
		sumw += chixy_ww[j];
		avex += chixy_ww[j] * chixy_xx[j];
		avey += chixy_ww[j] * chixy_yy[j];
	}
	if (sumw == 0.0) sumw = Big;
	avex /= sumw;
	avey /= sumw;
	chixy_aa = avey - b*avex;
	
	double	ans = -chixy_offs;
	for (unsigned long j=1; j<=nn; j++)
		ans += chixy_ww[j] * sqr( chixy_yy[j] - chixy_aa - b * chixy_xx[j] );
	return ans;
};




void	MyMath::line_fit_func::fitexy( const double_vector1& x, const double_vector1& y, const double_vector1& sigx, const double_vector1& sigy,
				double& a, double& b )
{
	const double 	Potn 	= 1.571000;
	const double	Big 	= 1.0e30;
	const double	Acc 	= 1.0e-3;
	const double 	Pi		= 3.14159265358979323846264338;
	
	size_t	ndat = x.length();

	chixy_xx = double_vector1( ndat );
	chixy_yy = double_vector1( ndat );
	chixy_sx = double_vector1( ndat );
	chixy_sy = double_vector1( ndat );
	chixy_ww = double_vector1( ndat );
	
	//double	varx, vary, dum1;
	//avevar( x, dum1, varx );
	//avevar( y, dum1, vary );
	double	varx = x.variance();
	double	vary = y.variance();
	double	scale = sqrt( varx/vary );
	
	for (unsigned long j=1; j<=ndat; j++)
	{	chixy_xx[j] = x[j];
		chixy_yy[j] = y[j] * scale;
		chixy_sx[j] = sigx[j];
		chixy_sy[j] = sigy[j] * scale;
		chixy_ww[j] = sqrt( sqr(chixy_sx[j]) + sqr(chixy_sy[j]) );
	}
	fit( chixy_xx, chixy_yy, chixy_ww, true, a, b );
	
	double_vector1 ang(6);
	double_vector1 ch(6);
	ang[1] = 0.0;
	ang[2] = atan(b);
	ang[4] = 0.0;
	ang[5] = ang[2];
	ang[6] = Potn;
	chixy_offs = ang[1];
	
	for (int j=4; j<=6; j++)	ch[j] = chixy( ang[j] );
	mnbrak( ang[1], ang[2], ang[3], ch[1], ch[2], ch[3], chixy );
	
	chisq = brent( ang[1], ang[2], ang[3], chixy, Acc , b );
	chisq = chixy( b );
	a = chixy_aa;
	q = gammq( 0.5*(ndat-2), chisq*0.5 );
	double	r2 = 0;
	for (unsigned long j=1; j<=ndat; j++)
		r2 += chixy_ww[j];
	r2 = 1.0/r2;
	
	double	bmx = Big;
	double	bmn = Big;
	chixy_offs = chisq + 1.0;
	for (int j=1; j<=6; j++)
	{	if (ch[j] > chixy_offs)
		{	double	d1 = fabs( ang[j] - b );
			while (d1 >= Pi) d1 -= Pi;
			double	d2 = Pi - d1;
			if (ang[j] < b)			{	double	swap = d1;	d1 = d2;	d2 = swap;	}
			if (d1 < bmx) bmx = d1;
			if (d2 < bmn) bmn = d2;
		}
	}
	if (bmx < Big)
	{	bmx 		= zbrent( chixy, b, b+bmx, Acc ) - b;
		double	amx = chixy_aa - a;
		bmn 		= zbrent( chixy, b, b-bmn, Acc ) - b;
		double	amn = chixy_aa - a;
		sigb = sqrt( 0.5*(bmx*bmx+bmn*bmn) ) / (scale * sqr(cos(b)));
		siga = sqrt( 0.5*(amx*amx+amn*amn) + r2) / scale;
	}
	else
		sigb = siga = Big;
	a /= scale;
	b = tan(b)/scale;
};



#pragma mark -
#pragma mark //	ode integration



double	MyMath::ode_integrator::odeint(	double_vector& ystart, double x1, double x2 )
{
	const size_t	MAXSTP = 10000000;		//	10 millions
	const double	TINY = 1.0e-1;

	size_t			nvar = ystart.size();
	double_vector	yscal( nvar );
	double_vector	ytemp( nvar );
	double_vector	y( nvar );
	double_vector	dydx( nvar );
	double			x = x1;
	double			h = MyMath::set_sign( h1, x2-x1 );
	nok = nbad = 0;
	
	for (size_t i=1; i<=nvar; i++)	y[i] = ystart[i];

	double	htemp;
	bool	small_for_saving = false;
	
	for (size_t nstp=1; nstp<=MAXSTP; nstp++)
	{	(*derivs)( x, y, dydx );
		for (size_t i=1; i<=nvar; i++)
			yscal[i] = fabs( y[i] ) + fabs( dydx[i]*h ) + TINY;
			
		if (savor && savor->CheckSave(x))
			savor->SaveState( x, y );

		if (savor)
		{	htemp = h;
			h = savor->NextSavingStep( x, h );				//	don't overshoot saving step
			small_for_saving = (fabs(htemp) > fabs(h));
		}

		if (defaultEndInteg( x+h, x1, x2, y ))
			h = x2-x;										//	don't overshoot final step
		
		double	hnext, hdid;
		for (size_t i=1; i<=nvar; i++)	ytemp[i] = y[i];					//	keep last y (last x is x-hdid)				
		
		(*odeStepF)( y, dydx, x, h, eps, yscal, hdid, hnext, derivs );		//	does change x, y and (hdid, hnext)
		
		if (hdid == h)	++nok;
		else			++nbad;
		if ((hdid == h) && small_for_saving)	hnext = htemp;
		
		//bool	true_end = defaultEndInteg( x, x1, x2, y );
		
		if ((*endInteg)( x, x1, x2, y ))					//	if switching, adjust to a better position
		{	x -= hdid;
			for (size_t i=1; i<=nvar; i++)	y[i] = ytemp[i];		//	rewind ...
			
			double	tol = 1e-6;
			solveFunction	solF( *this, x, hdid, y, dydx, yscal );
			double	fa = solF(0);
			double	fb = solF(hdid);
			h = zbrent( solF, 0, hdid, tol );
			(*odeStepF)( y, dydx, x, h, eps, yscal, hdid, hnext, derivs );
			
			for (size_t i=1; i<=nvar; i++)
				ystart[i] = y[i];				//	for the next start
			return x;							//	end of subroutine HERE.
		}
			
		//if ( true_end || (*endInteg)( x, x1, x2, y ) )		//	check for end AND possible switching
		if ( defaultEndInteg( x, x1, x2, y ) )		//	check for end
		{	
			if (savor && savor->CheckSave(x))
				savor->SaveState( x, y );
	
			for (size_t i=1; i<=nvar; i++)
				ystart[i] = y[i];				//	for the next start
			return x;							//	end of subroutine HERE.
		}
		
		if (fabs(hnext) <= hmin)	MyMath::ThrowError((char *)"Step size too small in odeint");
		h = hnext;
	}
	MyMath::ThrowError((char *)"Too many steps in routine odeint");
}



void	MyMath::rkck( const double_vector& y, const double_vector& dydx, const double& x, const double& h,
						double_vector& yout, double_vector& yerr, derivFuncType derivs )
{	static const double	a2=0.2,	 	a3=0.3, 		a4=0.6, 		a5=1.0, 	a6=0.875,
						b21=0.2,	b31=3.0/40.0,	b32=9.0/40.0,
						b41=0.3,	b42 = -0.9,		b43=1.2,
						b51 = -11.0/54.0,	b52=2.5,			b53 = -70.0/27.0,	b54=35.0/27.0,
						b61=1631.0/55296.0,	b62=175.0/512.0,	b63=575.0/13824.0,	b64=44275.0/110592.0,	b65=253.0/4096.0,
						c1=37.0/378.0,		c3=250.0/621.0,		c4=125.0/594.0,		c6=512.0/1771.0,		dc5 = -277.0/14336.0;

	static const double	dc1=c1-2825.0/27648.0,	dc3=c3-18575.0/48384.0,	dc4=c4-13525.0/55296.0,	dc6=c6-0.25;

	size_t n = y.length();		//	= dydx.length() = yout.length() = yerr.length()
	
	double_vector1	ak2(n);
	double_vector1	ak3(n);
	double_vector1	ak4(n);
	double_vector1	ak5(n);
	double_vector1	ak6(n);
	double_vector1	ytemp(n);
	
	for (size_t i=1;i<=n;i++)		ytemp[i] = y[i] + b21*h*dydx[i];
	(*derivs)( x+a2*h, ytemp, ak2);
	for (size_t i=1;i<=n;i++)		ytemp[i] = y[i] + h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)( x+a3*h, ytemp, ak3);
	for (size_t i=1;i<=n;i++)		ytemp[i] = y[i] + h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)( x+a4*h, ytemp, ak4);
	for (size_t i=1;i<=n;i++)		ytemp[i] = y[i] + h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)( x+a5*h, ytemp, ak5);
	for (size_t i=1;i<=n;i++)		ytemp[i] = y[i] + h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)( x+a6*h, ytemp, ak6);
	for (size_t i=1;i<=n;i++)		yout[i]  = y[i] + h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (size_t i=1;i<=n;i++)		yerr[i]	 = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
};

void	MyMath::rkqs( double_vector& y, const double_vector& dydx, double& x, const double& htry, double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs )
{		//	parameters that can be adjusted:
	static const double	safety =  0.9;
	static const double	pgrow  = -0.2;
	static const double	pshrnk = -0.25;
	static const double	errcon =  1.89e-4;

	size_t n = y.length();		//	= dydx.length() = yscal.length()

	double_vector1	ytemp(n);
	double_vector1	yerr(n);

	double	h = htry;
	Boolean	notDone = true;
	while (notDone)
	{	rkck( y, dydx, x, h, ytemp, yerr, derivs );
		double	errmax = 0.0;
		for (size_t i=1; i<=n; i++)	errmax = std::max( errmax, std::abs(yerr[i]/yscal[i]) );
		errmax /= eps;
		
		if (errmax > 1.0)		//	bad step -> shrink time step
		{	h = safety * h * pow( errmax, pshrnk );
			if (h < 0.1*h) h *= 0.1;
			double	xnew = x + h;
			if (xnew == x) ThrowError((char *)"stepsize underflow in rkqs");
			//continue;
		}
		else					//	good step -> grow time step
		{	if (errmax > errcon) 	hnext = safety * h * pow( errmax, pgrow );
			else 					hnext = 5.0 * h;
			x += (hdid=h);
			for (size_t i=1; i<=n; i++)	y[i] = ytemp[i];
			notDone = false;
		}
	}
};

void	MyMath::bs_extrapol::pzextr( size_t iest, double xest, const double_vector1& yest, double_vector1& yz, double_vector1& dy )
{
	size_t 			nv = yest.length();		// = yz.length() = dy.length()
	double_vector1	c(nv);
	
	x[iest] = xest;
	for (size_t j=1; j<=nv; j++)		dy[j] = yz[j] = yest[j];
	if (iest == 1)
		for (size_t j=1; j<=nv; j++)	d(j,1) = yest[j];
	else
	{	for (size_t j=1; j<=nv; j++) 	c[j] = yest[j];
		for (size_t k1=1; k1<iest; k1++)
		{	double	delta = 1.0/(x[iest-k1]-xest);
			double	f1	  = xest*delta;
			double	f2	  = x[iest-k1]*delta;
			for (size_t j=1; j<=nv; j++)
			{	double	q = d(j,k1);
				d(j,k1) = dy[j];
				delta	= c[j]-q;
				dy[j]	= f1*delta;
				c[j]	= f2*delta;
				yz[j]  += dy[j];
			}
		}
		for (size_t j=1; j<=nv; j++)	d(j,iest) = dy[j];
	}
};


void	MyMath::bs_extrapol::rzextr( size_t iest, double xest, const double_vector1& yest, double_vector1& yz, double_vector1& dy )
{
	size_t 			nv = yest.length();		// = yz.length() = dy.length()

	double_vector1	fx(iest);

	x[iest] = xest;
	if (iest == 1)
		for (size_t j=1; j<=nv; j++)
		{	yz[j]  = yest[j];
			d(j,1) = yest[j];
			dy[j]  = yest[j];
		}
	else
	{	for (size_t k=1; k<iest; k++)	fx[k+1] = x[iest-k]/xest;
		for (size_t j=1; j<=nv; j++)
		{	double	v  = d(j,1);
			double	yy = yest[j];
			double	c  = yest[j];
			double	ddy;
			d(j,1) = yest[j];
			for (size_t k=2; k<=iest; k++)
			{	double	b1 = fx[k]*v;
				double  b = b1-c;
				if (b)
				{	b = (c-v)/b;
					ddy = c*b;
					c = b1*b;
				}
				else
					ddy = v;
				if (k != iest) v = d(j,k);
				d(j,k) = ddy;
				yy += ddy;
			}
			dy[j] = ddy;
			yz[j] = yy;
		}
	}
}

void	MyMath::mmid( const double_vector1& y, const double_vector1& dydx, double xs, double htot, size_t nstep, double_vector1& yout, derivFuncType derivs )
{
	size_t 	nvar = y.length();		// = dydx.length() = yout.length()

	double_vector ym(nvar);
	double_vector yn(nvar);
	
	double	h = htot/nstep;
	for (size_t i=1; i<=nvar; i++)
	{	ym[i] = y[i];	yn[i] = y[i]+h*dydx[i];		};
	double	x = xs+h;
	
	(*derivs)( x, yn, yout );
	double	h2 = 2.0*h;
	for (size_t n=2; n<=nstep; n++)
	{	for (size_t i=1; i<=nvar; i++)
		{	double swap = ym[i]+h2*yout[i];
			ym[i] = yn[i];
			yn[i] = swap;
		}
		x += h;
		(*derivs)( x, yn, yout );
	}
	for (size_t i=1; i<=nvar; i++)
		yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]);
	//NR_free_vector(yn,1,nvar);
	//NR_free_vector(ym,1,nvar);
}

void	MyMath::bsstep( double_vector& y, const double_vector& dydx, double& xx, const double& htry, double eps, const double_vector1& yscal,
	double& hdid, double& hnext, derivFuncType derivs )

{
//	int i,iq,k,kk,km;
//	float eps1,errmax,fact,red,scale,work,wrkmin,xest;
	
	static double	epsold = -1.0;
	static int 		first  = 1;
	static double	xnew;
	static int 		kmax, kopt;
	
	size_t 			nv = y.length();		// = dydx.length() = yscal.length()
	
		//	parameters that can be adjusted:
	static const size_t	k_maxx = 8;
//	static const int	i_maxx = k_maxx + 1;
	static int			nseq[k_maxx+2] = {0,2,4,6,8,10,12,14,16,18};
	static const double	safe1	= 0.25;
	static const double	safe2	= 0.7;
	static const double	red_max	= 1.0e-5;
	static const double	red_min	= 0.7;
	static const double	tiny	= 1.0e-30;
	static const double	scalmx	= 0.1;
	
	static double	a[k_maxx+2];
	static double	alf[k_maxx+1][k_maxx+1];

	bs_extrapol	  bs_ext( nv, k_maxx );
	double_vector err(k_maxx);
	
	double_vector yerr(nv);
	double_vector ysav(nv);
	double_vector yseq(nv);
	
	if (eps != epsold)
	{	hnext = xnew = -1.0e29;
		double	eps1 = safe1*eps;
		a[1] = nseq[1] + 1;
		for (size_t k=1; k<=k_maxx; k++)
			a[k+1]=a[k]+nseq[k+1];
		
		for (size_t iq=2; iq<=k_maxx; iq++)
			for (size_t k=1; k<iq; k++)
				alf[k][iq] = pow( eps1, (a[k+1]-a[iq+1])/ ((a[iq+1]-a[1]+1.0)*(2*k+1)));

		epsold = eps;
		size_t kopt = 2;
		for (; kopt<k_maxx; kopt++)
			if ( a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax = kopt;
	}
	
	double	h = htry;
	for (size_t i=1; i<=nv; i++)	ysav[i] = y[i];
	if (xx != xnew || h != hnext)
	{	first = 1;
		kopt = kmax;
	}
	Boolean	reduct	 = false;
	Boolean	exitflag = false;
	
	size_t	km;
	size_t	k;
	double	errmax, red;
	for (;;)
	{	for (k=1; k<=kmax; k++)
		{	xnew = xx + h;
			if (xnew == xx)	ThrowError((char *)"step size underflow in bsstep");
			
			mmid( ysav, dydx, xx, h, nseq[k], yseq, derivs );
			double xest = sqr( h/nseq[k] );
			bs_ext.pzextr( k, xest, yseq, y, yerr );
			
			if (k != 1)
			{	errmax = tiny;
				for (size_t i=1;i<=nv;i++) errmax = std::max( errmax, std::abs(yerr[i]/yscal[i]) );
				errmax /= eps;
				km = k-1;
				err[km] = pow(errmax/safe1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag = true;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red = safe2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red = 1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red = alf[km][kmax-1]*safe2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red = alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		
		red = min( red,red_min );
		red = std::max( red,red_max );
		h *= red;
		reduct = true;
	}
	xx    = xnew;
	hdid  = h;
	first = 0;
	
	double	wrkmin = 1.0e35;
	double	scale  = 1.0;
	for (size_t kk=1; kk<=km; kk++)
	{	double	fact = std::max( err[kk], scalmx );
		double	work = fact*a[kk+1];
		if (work < wrkmin)
		{	scale = fact;
			wrkmin=work;
			kopt = kk+1;
		}
	}
	
	hnext = h/scale;
	if (kopt >= k && kopt != kmax && !reduct)
	{	double	fact = std::max( scale/alf[kopt-1][kopt], scalmx );
		if (a[kopt+1]*fact <= wrkmin)
		{	hnext = h/fact;
			kopt++;
		}
	}
}


void	MyMath::vFunction_fdjac::fdjac( const double_vector& x, const double_vector& fvec, double_matrix& df ) const
{		//	parameters that can be adjusted:
	static const double	eps = 1.0e-6;
	
	size_t	n = x.length();
	double_vector	f1(n);
	double_vector	x1 = x;
	
	double	h;

	for (size_t j=1; j<=n; j++)
	{	double	temp = x1[j];
		h = eps * std::abs(temp);
		if (h == 0.0) h = eps;
		x1[j] = temp + h;
		h = x1[j]-temp;
		vecF( x1, f1 );
		x1[j]=temp;
		for (size_t i=1; i<=n; i++)
			df(i,j) = (f1[i]-fvec[i])/h;
	}
};

void	MyMath::jacobn( const double& x, const double_vector& y, double_vector& dfdx, double_matrix& dfdy, derivFuncType derivs )
{		//	parameters that can be adjusted:
	static const double	eps = 1.0e-6;
	
	size_t	n = y.length();
	double_vector	dydx0(n);
	double_vector	dydx1(n);
	
	(*derivs)( x, y, dydx0 );
	
	vOdeFunction_fdjac z(derivs,x);
	z.fdjac( y, dydx0, dfdy );
	
	//double	h = eps * abs(temp);
	
	(*derivs)( x + eps, y, dydx1 );
	for (unsigned long i=1; i<=n; i++)
		dfdx[i] = (dydx1[i] - dydx0[i])/eps;
};

void	MyMath::stiff ( double_vector& y, const double_vector& dydx, double& x, const double& htry,
			double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs )
{
	size_t 			n = y.length();		// = dydx.length() = yscal.length()
	
	static const double	safety	=  0.9;
	static const double	grow	=  1.5;
	static const double	pgrow	= -0.25;
	static const double	shrnk	=  0.5;
	static const double	pshrnk	=  (-1.0/3.0);
	static const double	errcon	=  0.1296;
	static const size_t	maxtry	=  40;
		//	Shampine
	static const double	gam		=  (1.0/2.0);
	static const double	a21		=  2.0;
	static const double	a31		=  (48.0/25.0);
	static const double	a32		=  (6.0/25.0);
	static const double	c21		= -8.0;
	static const double	c31		=  (372.0/25.0);
	static const double	c32		=  (12.0/5.0);
	static const double	c41		=  (-112.0/125.0);
	static const double	c42		=  (-54.0/125.0);
	static const double	c43		=  (-2.0/5.0);
	static const double	b1		=  (19.0/9.0);
	static const double	b2		=  (1.0/2.0);
	static const double	b3		=  (25.0/108.0);
	static const double	b4		=  (125.0/108.0);
	static const double	e1		=  (17.0/54.0);
	static const double	e2		=  (7.0/36.0);
	static const double	e3		=  0.0;
	static const double	e4		=  (125.0/108.0);
	static const double	c1x		=  (1.0/2.0);
	static const double	c2x		=  (-3.0/2.0);
	static const double	c3x		=  (121.0/50.0);
	static const double	c4x		=  (29.0/250.0);
	static const double	a2x		=  1.0;
	static const double	a3x		=  (3.0/5.0);
/*		//	Kaps-Rentrop
	static const double	gam		=  0.231;
	static const double	a21		=  2.0;
	static const double	a31		=  4.52470820736;
	static const double	a32		=  4.16352878860;
	static const double	c21		= -5.07167533877;
	static const double	c31		=  6.02015272865;
	static const double	c32		=  0.159750684673;
	static const double	c41		= -1.856343618677;
	static const double	c42		= -8.50538085819;
	static const double	c43		= -2.08407513602;
	static const double	b1		=  3.95750374663;
	static const double	b2		=  4.62489238836;
	static const double	b3		=  0.617477263873;
	static const double	b4		=  1.282612945268;
	static const double	e1		= -2.30215540292;
	static const double	e2		= -3.07363448539;
	static const double	e3		=  0.873280801802;
	static const double	e4		=  1.282612945268; //= b4;
	static const double	c1x		=  gam;
	static const double	c2x		= -0.396296677520e-1;
	static const double	c3x		=  0.550778939579;
	static const double	c4x		= -0.553509845700e-1;
	static const double	a2x		=  0.462;
	static const double	a3x		=  0.880208333333;
*/
	int_vector		indx(n);
	double_matrix	a(n,n);
	double_matrix	dfdy(n,n);
	double_vector	dfdx(n);
	double_vector	dywk(n);
	double_vector	err(n);
	double_vector	g1(n);
	double_vector	g2(n);
	double_vector	g3(n);
	double_vector	g4(n);
	double_vector	ysav(n);
	
	double	xsav = x;
	for (size_t i=1; i<=n; i++)
	{	ysav[i]  = y[i];
		dywk[i] = dydx[i];
	}
	jacobn( xsav, ysav, dfdx, dfdy, derivs );
	double	h = htry;
	for (size_t jtry=1; jtry<=maxtry; jtry++)
	{	for (size_t i=1; i<=n; i++)
		{	for (size_t j=1; j<=n; j++)	a(i,j) = -dfdy(i,j);
			a(i,i) += 1.0/(gam*h);
		}
		Boolean d;
		ludcmp( a, indx, d );
		for (size_t i=1; i<=n; i++)		g1[i] = dywk[i]+h*c1x*dfdx[i];
		lubksb( a, indx, g1 );
		for (size_t i=1; i<=n; i++)		y[i]  = ysav[i]+a21*g1[i];
		x = xsav + a2x*h;
		(*derivs)( x, y, dywk );
		for (size_t i=1; i<=n; i++)		g2[i] = dywk[i]+h*c2x*dfdx[i]+c21*g1[i]/h;
		lubksb( a, indx, g2 );
		for (size_t i=1; i<=n; i++)		y[i]  = ysav[i]+a31*g1[i]+a32*g2[i];
		x = xsav+a3x*h;
		(*derivs)( x, y, dywk );
		for (size_t i=1; i<=n; i++)		g3[i] = dywk[i]+h*c3x*dfdx[i]+(c31*g1[i]+c32*g2[i])/h;
		lubksb( a, indx, g3 );
		for (size_t i=1; i<=n; i++)		g4[i] = dywk[i]+h*c4x*dfdx[i]+(c41*g1[i]+c42*g2[i]+c43*g3[i])/h;
		lubksb( a, indx, g4 );
		for (size_t i=1; i<=n; i++)
		{	y[i]	= ysav[i]+b1*g1[i]+b2*g2[i]+b3*g3[i]+b4*g4[i];
			err[i]	= e1*g1[i]+e2*g2[i]+e3*g3[i]+e4*g4[i];
		}
		x = xsav+h;
		if (x == xsav) ThrowError((char *)"stepsize not significant in stiff");
		double	errmax = 0.0;
		for (size_t i=1; i<=n; i++) errmax = std::max( errmax, std::abs(err[i]/yscal[i]) );
		errmax /= eps;
		if (errmax <= 1.0)
		{	hdid  = h;
			hnext = (errmax > errcon ? safety*h*pow(errmax,pgrow) : grow*h);
			return;
		}
		else
		{	hnext = safety*h*pow(errmax,pshrnk);
			if (hnext < shrnk*h) hnext=shrnk*h;
			h = hnext;
		}
	}
	ThrowError((char *)"exceeded MAXTRY in stiff");
}


void	MyMath::simpr( const double_vector& y,  const double_vector& dydx, const double_vector& dfdx, const double_matrix& dfdy,
	const double& xs, const double& htot, size_t nstep, double_vector& yout, derivFuncType derivs )
{
//	int i,j,nn,*indx;
//	float d,h,x,**a,*del,*ytemp;
	
	size_t		n = y.length();		//	 = yout.length() = dydx.length() = dfdx.length() = dfdy.dim1() = dfdy.dim2() 

	int_vector		indx(n);
	double_matrix	a(n,n);
	double_vector	del(n);
	double_vector	ytemp(n);
	
	double	h = htot/nstep;
	for (size_t i=1; i<=n; i++)
	{	for (size_t j=1; j<=n; j++) a(i,j) = -h*dfdy(i,j);
		++a(i,i);
	}
	Boolean	d;
	ludcmp( a, indx, d );
	for (size_t i=1; i<=n; i++)		yout[i] = h*(dydx[i]+h*dfdx[i]);
	lubksb( a, indx, yout );
	for (size_t i=1; i<=n; i++)		ytemp[i] = y[i]+(del[i]=yout[i]);
	double	x = xs+h;
	(*derivs)( x, ytemp, yout );
	for (size_t nn=2; nn<=nstep; nn++)
	{	for (size_t i=1; i<=n; i++)	yout[i] = h*yout[i]-del[i];
		lubksb( a, indx, yout );
		for (size_t i=1; i<=n; i++)	ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		(*derivs)( x, ytemp, yout );
	}
	for (size_t i=1; i<=n; i++)		yout[i] = h*yout[i]-del[i];
	lubksb( a, indx, yout );
	for (size_t i=1; i<=n; i++)		yout[i] += ytemp[i];
};


void	MyMath::stifbs( double_vector& y, const double_vector& dydx, double& xx, const double& htry, double eps, const double_vector1& yscal,
	double& hdid, double& hnext, derivFuncType derivs )
{
	static double	epsold = -1.0;
	static int 		first  =  1;
	static int 		nvold  = -1;
	static double	xnew;
	static int 		kmax, kopt;
	
	size_t 			nv = y.length();		// = dydx.length() = yscal.length()

		//	parameters that can be adjusted:
	static const int	k_maxx = 7;
//	static const int	i_maxx = k_maxx + 1;
	static int			nseq[k_maxx+2] = {0,2,6,10,14,22,34,50,70};
	static const double	safe1	= 0.25;
	static const double	safe2	= 0.7;
	static const double	red_max	= 1.0e-5;
	static const double	red_min	= 0.7;
	static const double	tiny	= 1.0e-30;
	static const double	scalmx	= 0.1;
	
	static double	a[k_maxx+2];
	static double	alf[k_maxx+1][k_maxx+1];

	bs_extrapol	  bs_ext( nv, k_maxx );
	double_vector err(k_maxx);
	
	double_vector yerr(nv);
	double_vector ysav(nv);
	double_vector yseq(nv);
	double_vector dfdx(nv);
	double_matrix dfdy(nv,nv);
	
	if(eps != epsold || nv != nvold)
	{	hnext = xnew = -1.0e29;
		double	eps1  = safe1*eps;
		a[1]  = nseq[1]+1;
		for (size_t k=1; k<=k_maxx; k++)	a[k+1]=a[k]+nseq[k+1];
		for (size_t iq=2; iq<=k_maxx; iq++)
		{	for (size_t k=1; k<iq; k++)
				alf[k][iq] = pow( eps1, ((a[k+1]-a[iq+1] ) / ((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold = eps;
		nvold  = nv;
		a[1]  += nv;
		for (size_t k=1; k<=k_maxx; k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=2; kopt<k_maxx; kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax = kopt;
	}
	
	double	h = htry;
	for (size_t i=1; i<=nv; i++)	ysav[i] = y[i];
	jacobn( xx, y, dfdx, dfdy, derivs );
	if ( xx != xnew || h != hnext )
	{	first = 1;
		kopt  = kmax;
	}
	
	Boolean	reduct	 = false;
	Boolean	exitflag = false;
	size_t	km, k;
	for (;;)
	{	double	red;
		for (k=1; k<=kmax; k++)
		{	xnew = xx + h;
			if (xnew == xx)	ThrowError((char *)"step size underflow in stifbs");
			simpr( ysav, dydx, dfdx, dfdy, xx, h, nseq[k], yseq, derivs );
			
			double xest = sqr( h/nseq[k] );
			bs_ext.pzextr( k, xest, yseq, y, yerr );
			
			double	errmax;
			if (k != 1)
			{	errmax = tiny;
				for (size_t i=1; i<=nv; i++) errmax = std::max( errmax, std::abs(yerr[i]/yscal[i]) );
				errmax /= eps;
				km = k-1;
				err[km] = pow( errmax/safe1, 1.0/(2*km+1) );
			}
			if (k != 1 && (k >= kopt-1 || first))
			{	if (errmax < 1.0)
				{	exitflag = true;
					break;
				}
				if (k == kmax || k == kopt+1)
				{	red = safe2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km])
				{	red = 1.0/err[km];
					break;
				}
				else if (kopt == kmax && alf[km][kmax-1] < err[km])
				{	red = alf[km][kmax-1]*safe2/err[km];
					break;
				}
				else if (alf[km][kopt] < err[km])
				{	red = alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red = min(red,red_min);
		red = std::max(red,red_max);
		h *= red;
		reduct = true;
	}
	xx 	 = xnew;
	hdid = h;
	first = 0;
	double	wrkmin = 1.0e35;
	double	fact, scale;
	for (size_t kk=1; kk<=km; kk++)
	{	fact = std::max( err[kk], scalmx );
		double	work = fact*a[kk+1];
		if (work < wrkmin)
		{	scale  = fact;
			wrkmin = work;
			kopt   = kk+1;
		}
	}
	hnext = h/scale;
	if (kopt >= k && kopt != kmax && !reduct)
	{	fact = std::max( scale/alf[kopt-1][kopt], scalmx );
		if (a[kopt+1]*fact <= wrkmin)
		{	hnext = h/fact;
			kopt++;
		}
	}
	
};



#pragma mark -
#pragma mark //	quadrature


void	MyMath::gauleg( double x1, double x2, double_vector& x, double_vector& w, size_t n )
{	//	en long double precision?
	
	static const double	eps	= 3.0e-11;

	size_t	m  = (n+1)/2;
	double	xm = 0.5*(x2+x1);
	double	xl = 0.5*(x2-x1);
	for (size_t i=1; i<=m; i++)
	{	double	z = cos( 3.141592654 * (i-0.25)/(n+0.5) );
		double	z1, pp;
		do
		{	double	p1 = 1.0;
			double	p2 = 0.0;
			for (size_t j=1;j<=n;j++)
			{	double p3 = p2;
				p2 = p1;
				p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp = n*(z*p1-p2)/(z*z-1.0);
			z1 = z;
			z  = z1-p1/pp;
		}
		while (std::abs(z-z1) > eps);
		
		x[i]	 = xm - xl*z;
		x[n+1-i] = xm + xl*z;
		w[i]	 = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i] = w[i];
	}
};

#pragma mark -
#pragma mark //	en vrac




void	MyMath::balanc( double_matrix& a )
{	
	size_t 	n = a.dim1();		//	 = a.dim2();

	static const double	radix = 2.0;
	static const double	sqrdx = radix*radix;

	size_t	last = 0;
	while (last == 0)
	{	last = 1;
		for (size_t i=1; i<=n; i++)
		{	double	r = 0.0;
			double	c = 0.0;
			for (size_t j=1; j<=n; j++)
				if (j != i)
					{	c += std::abs(a(j,i));	r += std::abs(a(i,j));	}
			if (c && r)
			{	double	g = r/radix;
				double	f = 1.0;
				double	s = c+r;
				while (c<g)
				{	f *= radix;	c *= sqrdx;	}
				g = r*radix;
				while (c>g)
				{	f /= radix;	c /= sqrdx;	}
				if ((c+r)/f < 0.95*s)
				{	last = 0;
					g = 1.0/f;
					for (size_t j=1; j<=n; j++) a(i,j) *= g;
					for (size_t j=1; j<=n; j++) a(j,i) *= f;
	}	}	}	}
};


void	MyMath::elmhes( double_matrix& a )
{
	size_t 	n = a.dim1();		//	 = a.dim2();
	
	for (size_t m=2; m<n; m++)
	{	double	x = 0.0;
		size_t	i = m;
		for (size_t j=m; j<=n; j++)
		{	if ( std::abs( a(j,m-1) ) > std::abs(x) )
			{	x = a(j,m-1);
				i = j;
			}
		}
		if (i != m)
		{	for (size_t j=m-1; j<=n; j++)	swap( a(i,j), a(m,j) );
			for (size_t j=1; j<=n; j++)		swap( a(j,i), a(j,m) );
		}
		if (x)
		{	for (size_t i=m+1; i<=n; i++)
			{	double	y = a(i,m-1);
				if ( y != 0.0)
				{	y /= x;
					a(i,m-1) = y;
					for (size_t j=m; j<=n; j++)	a(i,j) -= y*a(m,j);
					for (size_t j=1; j<=n; j++)	a(j,m) += y*a(j,i);
	}	}	}	}
};



void	MyMath::hqr( double_matrix& a, double_vector& wr, double_vector& wi )
{
///	int nn,m,l,k,j,its,i,mmin;
//	float z,y,x,w,s,r,q,p;

	size_t 	n = a.dim1();		//	 = a.dim2();

	double	anorm = std::abs(a(1,1));
	for (size_t i=2; i<=n; i++)
		for (size_t j=(i-1); j<=n; j++)
			anorm += std::abs(a(i,j));
			
	size_t	nn = n;
	double	t  = 0.0;
	while (nn >= 1)
	{	size_t	its = 0;
		size_t	l;
		double	s;
		do {
			for (l=nn; l>=2; l--)
			{	s = std::abs(a(l-1,l-1)) + std::abs(a(l,l));
				if (s == 0.0) s=anorm;
				if ((float)(fabs(a(l,l-1)) + s) == s) break;
			}
			double	x = a(nn,nn);
			if (l == nn) 	{	wr[nn] = x+t;	wi[nn--] = 0.0;	}
			else
			{	double	y = a(nn-1,nn-1);
				double	w = a(nn,nn-1) * a(nn-1,nn);
				double	p, q, z, r;
				if (l == (nn-1))
				{	p = 0.5*(y-x);
					q = p*p+w;
					z = sqrt(std::abs(q));
					x += t;
					if (q >= 0.0)
					{	z = p + set_sign(z,p);
						wr[nn-1] = wr[nn] = x+z;
						if (z) wr[nn] = x-w/z;
						wi[nn-1] = wi[nn] = 0.0;
					}
					else
					{	wr[nn-1] = wr[nn] = x+p;
						wi[nn-1] = -(wi[nn]=z);
					}
					nn -= 2;
				}
				else
				{	if (its == 30) ThrowError((char *)"Too many iterations in hqr");
					if (its == 10 || its == 20)
					{	t += x;
						for (size_t i=1; i<=nn; i++) a(i,i) -= x;
						s = std::abs(a(nn,nn-1)) + std::abs(a(nn-1,nn-2));
						y = x = 0.75*s;
						w = -0.4375 *s*s;
					}
					++its;
					size_t	m;
					for (m=(nn-2); m>=l; m--)
					{	z = a(m,m);
						r = x-z;
						s = y-z;
						p = (r*s-w)/a(m+1,m)+a(m,m+1);
						q = a(m+1,m+1)-z-r-s;
						r = a(m+2,m+1);
						s = std::abs(p) + std::abs(q) + std::abs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						double	u = std::abs(a(m,m-1)) * ( std::abs(q) + std::abs(r) );
						double	v = std::abs(p) * ( std::abs(a(m-1,m-1)) + std::abs(z) + std::abs(a(m+1,m+1)));
						if ((float)(u+v) == v) break;
					}
					for (size_t i=m+2; i<=nn; i++)
					{	a(i,i-2) = 0.0;
						if (i != (m+2)) a(i,i-3) = 0.0;
					}
					for (size_t k=m; k<=nn-1; k++)
					{	if (k != m)
						{	p = a(k,k-1);
							q = a(k+1,k-1);
							r = 0.0;
							if (k != (nn-1)) r = a(k+2,k-1);
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0)
							{	p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=set_sign(sqrt(p*p+q*q+r*r),p)) != 0.0)
						{
							if (k == m)
							{	if (l != m)		a(k,k-1) = -a(k,k-1);	}
							else
								a(k,k-1) = -s*x;
							p += s;
							x  = p/s;
							y  = q/s;
							z  = r/s;
							q /= p;
							r /= p;
							for (size_t j=k; j<=nn; j++)
							{	p = a(k,j) + q*a(k+1,j);
								if (k != (nn-1))
								{
									p += r*a(k+2,j);
									a(k+2,j) -= p*z;
								}
								a(k+1,j) -= p*y;
								a(k,j) -= p*x;
							}
							
							size_t	mmin = nn<k+3 ? nn : k+3;
							for (size_t i=l;i<=mmin;i++)
							{	p = x*a(i,k) + y*a(i,k+1);
								if (k != (nn-1))
								{
									p += z*a(i,k+2);
									a(i,k+2) -= p*r;
								}
								a(i,k+1) -= p*q;
								a(i,k) -= p;
			}	}	}	}	}
		}
		while (l < nn-1);
	}
};


void	MyMath::eigvalsrt( double_vector& d, double_vector& v )
{	size_t 	n = d.length();	//	= v.length();

	for (size_t i=1; i<n; i++)
	{	size_t	k = i;
		double	p = d[k];
		for (size_t j=i+1; j<=n; j++)
			if (d[j] >= p)	{	k = j;	p = d[k];	}
		if (k != i)
		{	d[k] = d[i];
			d[i] = p;
			p 	 = v[i];
			v[i] = v[k];
			v[k] = p;
	}	}
};




void	MyMath::newt( double_vector& x, Boolean& notFound, const vectorFuncType& vecfunc )
{
//	int i,its,j,*indx;
//	float d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	size_t 	n = x.length();	//	= v.length();

	int_vector		indx(n);
	double_matrix	fjac(n,n);
	double_vector	g(n);
	double_vector	p(n);
	double_vector	xold(n);
	
	static const double	stpmx	= 100.0;
	static const double	tolx	= 1.0e-7;
	static const double	tolf	= 1.0e-4;
	static const double	tolmin	= 1.0e-6;
	static const size_t	maxits	= 200;
	
//	NR_fmin_fvec=NR_vector(1,n);
//	NR_fmin_nn=n;
//	NR_fmin_nrfuncv=vecfunc;
	vFunction_fmin vFmin( vecfunc, n );
	double	f = vFmin(x);
	
	double	test = 0.0;
	for (size_t i=1; i<=n; i++)
		if ( std::abs(vFmin.fx[i]) > test)	test = std::abs( vFmin.fx[i] );
		
	if (test < 0.01*tolf) return;
	double	sum = 0;
	for (size_t i=1; i<=n; i++) sum += sqr(x[i]);
	
	double	stpmax = stpmx * std::max( sqrt(sum), (double)n );
	for (size_t its=1; its<=maxits; its++)
	{
		fdjac( x, vFmin.fx, fjac, vecfunc );
		for (size_t i=1; i<=n; i++)
		{	double	sum = 0;
			for (size_t j=1; j<=n; j++)	sum += fjac(j,i) * vFmin.fx[j];
			g[i] = sum;
		}
		for (size_t i=1; i<=n; i++)	xold[i] = x[i];
		double	fold = f;
		for (size_t i=1; i<=n; i++)	p[i] = -vFmin.fx[i];
		Boolean	d;
		ludcmp( fjac, indx, d );
		lubksb( fjac, indx, p );
		lnsrch( xold, fold, g, p, x, f, stpmax, notFound, vFmin );
	
		double	test = 0.0;
		for (size_t	i=1;i<=n;i++)
			if ( std::abs(vFmin.fx[i]) > test ) test = std::abs(vFmin.fx[i]);
			
		if (test < tolf) {	notFound = false;		return;		}
		
		if (notFound)
		{	test=0.0;
			double	den = std::max( f, 0.5*n );
			for (size_t i=1; i<=n; i++)
			{	double	temp = std::abs(g[i]) * std::max( std::abs(x[i]), 1.0 )/den;
				if (temp > test) test = temp;
			}
			notFound = (test < tolmin ? 1 : 0);
			return;
		}
		test = 0.0;
		for (size_t i=1; i<=n; i++)
		{	double	temp = ( std::abs(x[i]-xold[i]) )/std::max( std::abs(x[i]), 1.0 );
			if (temp > test)	test = temp;
		}
		if (test < tolx) return;
	}
	ThrowError((char *)"MAXITS exceeded in newt");
}


void	MyMath::lnsrch( const double_vector& xold, const double& fold, const double_vector& g, double_vector& p, double_vector& x,
	double& f, double stpmax, Boolean& notFound, vectorFunction& func )
{
	size_t 	n = xold.length();	//	= g.length();
		
	static const double	alf		= 1.0e-4;
	static const double	tolx	= 1.0e-7;

	notFound = false;
	double	ssum = 0.0;
	for (size_t	i=1; i<=n; i++) 	ssum += p[i]*p[i];
	ssum = sqrt(ssum);
	if (ssum > stpmax)
		for (size_t i=1; i<=n; i++)	p[i] *= stpmax/ssum;
		
	double	slope=0.0;
	for (size_t i=1; i<=n; i++)		slope += g[i]*p[i];
	double	test = 0.0;
	for (size_t i=1; i<=n; i++)
	{	double	temp = std::abs(p[i])/std::max( std::abs(xold[i]), 1.0 );
		if (temp > test) test = temp;
	}
	double	alamin = tolx/test;
	double	alam = 1.0;
	for (;;)
	{	double	tmplam, f2, alam2, fold2;
		for (size_t i=1; i<=n; i++)		x[i] = xold[i] + alam*p[i];
		f = func( x );
		if (alam < alamin)
		{	for (size_t i=1;i<=n;i++)	x[i] = xold[i];
			notFound = true;
			return;						//	notFound = true
		}
		else if (f <= fold + alf*alam*slope) return;	//	OK: notFound = false
		else
		{	if (alam == 1.0)	tmplam = -slope/(2.0*(f-fold-slope));
			else
			{	double	rhs1 = f-fold-alam*slope;
				double	rhs2 = f2-fold2-alam2*slope;
				double	a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				double	b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0)	tmplam = -slope/(2.0*b);
				else
				{	double	disc = b*b-3.0*a*slope;
					if (disc<0.0) ThrowError((char *)"Roundoff problem in lnsrch.");
					else tmplam = (-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam > 0.5*alam)	tmplam = 0.5*alam;
			}
		}
		alam2 = alam;
		f2 = f;
		fold2 = fold;
		alam = std::max( tmplam, 0.1*alam );
	}									//	OK: notFound = false
};


void	MyMath::rotate( double_matrix& r, double_matrix& qt, size_t i, double a, double b )
{
	size_t 	n = r.length();	//	= qt.length();

	double	c, s;
	if (a == 0.0)
	{	c = 0.0;
		s = (b >= 0.0 ? 1.0 : -1.0);
	}
	else if ( std::abs(a) > std::abs(b) )
	{	double	fact = b/a;
		c = set_sign( 1.0/sqrt(1.0+(fact*fact)), a );
		s = fact*c;
	}
	else
	{	double	fact = a/b;
		s = set_sign( 1.0/sqrt(1.0+(fact*fact)), b );
		c = fact*s;
	}
	
	for (size_t	j=i; j<=n; j++)
	{	double	y = r(i,j);
		double	w = r(i+1,j);
		r(i,j)	 = c*y - s*w;
		r(i+1,j) = s*y + c*w;
	}
	for (size_t j=1; j<=n; j++)
	{	double	y = qt(i,j);
		double	w = qt(i+1,j);
		qt(i,j)   = c*y - s*w;
		qt(i+1,j) = s*y + c*w;
	}
}


void	MyMath::rsolv( const double_matrix& a, const double_vector& d, double_vector& b )
{
	size_t 	n = b.length();		//	= d.length() = a.dim1() = a.dim2();

	b[n] /= d[n];
	for (size_t i=n-1; i>=1; i--)
	{	double	sum = 0.0;
		for (size_t j=i+1; j<=n; j++)	sum += a(i,j)*b[j];
		b[i] = (b[i]-sum)/d[i];
	}
}

void	MyMath::qrdcmp( double_matrix& a, double_vector& c, double_vector& d, Boolean& singular )
{
	size_t 	n = d.length();		//	= c.length() = a.dim1() = a.dim2();
	double	scale = 0.0;

	singular = false;
	for (size_t k=1; k<n; k++)
	{	for (size_t i=k; i<=n; i++) scale = std::max( scale, std::abs(a(i,k)) );
		if (scale == 0.0)
		{	singular = true;
			c[k] = d[k] = 0.0;
		}
		else
		{	for (size_t i=k; i<=n; i++)	a(i,k) /= scale;
			double	sum = 0.0;
			for (size_t i=k;i<=n;i++) sum += sqr(a(i,k));
			double	sigma = set_sign( sqrt(sum), a(k,k) );
			a(k,k) += sigma;
			c[k] = sigma*a(k,k);
			d[k] = -scale*sigma;
			for (size_t j=k+1; j<=n; j++)
			{	double	sum = 0.0;
				for (size_t i=k; i<=n; i++) sum += a(i,k)*a(i,j);
				double	tau = sum/c[k];
				for (size_t i=k; i<=n; i++) a(i,j) -= tau*a(i,k);
			}
		}
	}
	d[n] = a(n,n);
	if (d[n] == 0.0)	singular = true;
}


void	MyMath::qrupdt( double_matrix& r, double_matrix& qt, double_vector& u, const double_vector& v )
{
	size_t 	n = u.length();		//	= r.dim1() = r.dim2() = qt.dim1() = qt.dim2() = v.length();

	size_t k;
	for (k=n; k>=1; k--)	if (u[k]) break;
	if (k < 1)	k=1;
	
	for (size_t i=k-1; i>=1; i--)
	{	rotate( r, qt, i, u[i], -u[i+1] );
		if (u[i] == 0.0)					u[i] = std::abs(u[i+1]);
		else if (std::abs(u[i]) > std::abs(u[i+1]))	u[i] = std::abs(u[i])*sqrt(1.0+sqr(u[i+1]/u[i]));
		else								u[i] = std::abs(u[i+1])*sqrt(1.0+sqr(u[i]/u[i+1]));
	}
	for (size_t j=1; j<=n; j++)		r(1,j) += u[1]*v[j];
	for (size_t i=1; i<k; i++)		rotate( r, qt, i, r(i,i), -r(i+1,i) );
}



void	MyMath::broydn( double_vector& x, Boolean& notFound, const vectorFuncType& vecfunc )
{
//	int i,its,j,k,restrt,sing,skip;
//	float den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
//	float *g,*p,**qt,**r,*s,*t,*w,*xold;

	static const double	stpmx	= 100.0;
	static const double	eps		= 1.0e-7;
	static const double	tolx	= 1.0e-7;
	static const double	tolf	= 1.0e-4;
	static const double	tolmin	= 1.0e-6;
	static const size_t	maxits	= 200;
	
/*	c=NR_vector(1,n);
	d=NR_vector(1,n);
	fvcold=NR_vector(1,n);
	g=NR_vector(1,n);
	p=NR_vector(1,n);
	qt=NR_matrix(1,n,1,n);
	r=NR_matrix(1,n,1,n);
	s=NR_vector(1,n);
	t=NR_vector(1,n);
	w=NR_vector(1,n);
	xold=NR_vector(1,n);
*/
	size_t 	n = x.length();	//	= v.length();
	
	double_vector	c(n);
	double_vector	d(n);
	double_vector	fvcold(n);
	double_vector	g(n);
	double_vector	p(n);
	double_matrix	qt(n,n);
	double_matrix	r(n,n);
	double_vector	s(n);
	double_vector	t(n);
	double_vector	w(n);
	double_vector	xold(n);
	
	//NR_fmin_fvec=NR_vector(1,n);
	//NR_fmin_nn=n;
	//NR_fmin_nrfuncv=vecfunc;
	//f=NR_fmin(x);
	vFunction_fmin vFmin( vecfunc, n );
	double	f = vFmin(x);
	
	double	test = 0.0;
	for (size_t i=1; i<=n; i++)
		if ( std::abs(vFmin.fx[i]) > test)	test = std::abs(vFmin.fx[i]);
		
	if (test < 0.01*tolf)	return;
	double	ssum = 0.0;
	for (size_t i=1; i<=n; i++)	ssum += sqr(x[i]);
	double	stpmax = stpmx*std::max( sqrt(ssum), (double)n );
	Boolean	restrt = true;
	
	for (size_t its=1; its<=maxits; its++)
	{	if (restrt)
		{	fdjac( x, vFmin.fx, r, vecfunc );
			Boolean	singular;
			qrdcmp( r, c, d, singular );
			if (singular)	ThrowError((char *)"singular Jacobian in broydn");
			for (size_t i=1; i<=n; i++)
			{	for (size_t j=1; j<=n; j++)	qt(i,j) = 0.0;
				qt(i,i) = 1.0;
			}
			for (size_t k=1; k<n; k++)
			{	if (c[k])
				{	for (size_t j=1; j<=n; j++)
					{	double	sum = 0.0;
						for (size_t i=k; i<=n; i++)	sum += r(i,k)*qt(i,j);
						sum /= c[k];
						for (size_t i=k; i<=n; i++)	qt(i,j) -= sum*r(i,k);
					}
				}
			}
			for (size_t i=1; i<=n; i++)
			{	r(i,i) = d[i];
				for (size_t j=1; j<i; j++)	r(i,j) = 0.0;
			}
		}
		else
		{	for (size_t i=1; i<=n; i++)	s[i] = x[i]-xold[i];
			for (size_t i=1; i<=n; i++)
			{	double	sum = 0.0;
				for (size_t j=i; j<=n; j++)	sum += r(i,j)*s[j];
				t[i] = sum;
			}
			
			Boolean	skip = true;
			for (size_t i=1; i<=n; i++)
			{	double	sum = 0.0;
				for (size_t j=1; j<=n; j++)	sum += qt(j,i)*t[j];
				w[i] = vFmin.fx[i] - fvcold[i] - sum;
				if ( std::abs(w[i]) >= eps * ( std::abs(vFmin.fx[i]) + fabs(fvcold[i]) ))
						skip = false;
				else	w[i] = 0.0;
			}
			if (!skip)
			{	for (size_t i=1;i<=n;i++)
				{	double	sum = 0.0;
					for (size_t j=1; j<=n; j++)	sum += qt(i,j)*w[j];
					t[i] = sum;
				}
				double	den = 0.0;
				for (size_t i=1; i<=n; i++)	den += sqr(s[i]);
				for (size_t i=1; i<=n; i++)	s[i] /= den;
				qrupdt( r, qt, t, s );
				for (size_t i=1; i<=n; i++)
				{	if (r(i,i) == 0.0) ThrowError((char *)"r singular in broydn");
					d[i] = r(i,i);
				}
			}
		}
		for (size_t i=1;i<=n;i++)
		{	double	sum = 0.0;
			for (size_t j=1; j<=n; j++)	sum += qt(i,j)*vFmin.fx[j];
			g[i] = sum;
		}
		for (size_t i=n; i>=1; i--)
		{	double	sum = 0.0;
			for (size_t j=1; j<=i; j++)	sum += r(j,i)*g[j];
			g[i] = sum;
		}
		for (size_t i=1; i<=n; i++)
		{	xold[i] = x[i];
			fvcold[i] = vFmin.fx[i];
		}
		double	fold = f;
		for (size_t i=1;i<=n;i++)
		{	double	sum = 0.0;
			for (size_t j=1; j<=n; j++)	sum += qt(i,j)*vFmin.fx[j];
			p[i] = -sum;
		}
		rsolv( r, d, p );
		lnsrch( xold, fold, g, p, x, f, stpmax, notFound, vFmin );
		test = 0.0;
		for (size_t i=1; i<=n; i++)
			if (std::abs(vFmin.fx[i]) > test)	test = std::abs(vFmin.fx[i]);
		if (test < tolf)
		{	notFound = false;
			return;								//	OK : notFound = false
		}
		if (notFound)
		{	if (restrt) return;					//	notFound = true
			else
			{	test = 0.0;
				double	den = std::max( f, 0.5*n );
				for (size_t i=1; i<=n; i++)
				{	double	temp = std::abs(g[i])*std::max(std::abs(x[i]),1.0)/den;
					if (temp > test)	test = temp;
				}
				if (test < tolmin)	return;		//	notFound = true
				else	restrt = true;
			}
		}
		else
		{	restrt = false;
			test = 0.0;
			for (size_t i=1; i<=n; i++)
			{	double	temp = ( std::abs(x[i]-xold[i]))/std::max(std::abs(x[i]),1.0);
				if (temp > test) test = temp;
			}
			if (test < tolx) return;			//	OK : notFound = false
		}
	}
	ThrowError((char *)"MAXITS exceeded in broydn");
	return;
}


void	MyMath::mnewt( double_vector& x, Boolean& notFound, const vFunction_fdjac& jacF )
{	
	static const double	tolx	= 1.0e-7;
	static const double	tolf	= 1.0e-4;
	static const size_t	maxits	= 200;

	size_t 	n = x.length();		//
	
	int_vector			indx(n);
	double_vector		p(n);
	double_vector		fvec(n);
	double_matrix		fjac(n,n);
	
	for (size_t k=1; k<=maxits; k++)
	{	jacF.vecF( x, fvec );				//	= (*vecfunc)( x, fvec )
		jacF.fdjac( x, fvec, fjac );
		double	errf = 0.0;
		for (size_t i=1; i<=n; i++)	errf += std::abs(fvec[i]);
		if (errf <= tolf) 	{	notFound = false;	return;		};
		for (size_t i=1; i<=n; i++)	p[i] = -fvec[i];
		Boolean	d;
		ludcmp( fjac, indx, d );
		lubksb( fjac, indx, p );
		double	errx = 0.0;
		for (size_t i=1; i<=n; i++)
		{	errx += std::abs(p[i]);
			x[i] += p[i];
		}
		if (errx <= tolx) 	{	notFound = false;	return;		};
	}
	notFound = true;
	return;
}


void	MyMath::mnbrak( double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, doubleFunction& func )
{
//	float ulim,u,r,q,fu,dum;

	const double gold = 1.618034;
	const double tiny = 1.0e-20;
	const double gLimit = 100.0;

	fa = func( ax );
	fb = func( bx );
	
	if (fb > fa)
	{	double	dum = ax;	ax = bx;	bx = dum;
				dum = fb;	fb = fa;	fa = dum;
	}
	cx = bx + gold * (bx-ax);
	fc = func( cx );
	
	while (fb > fc)
	{	double	r = (bx-ax) * (fb-fc);
		double	q = (bx-cx) * (fb-fa);
		double	u = bx - ((bx-cx)*q-(bx-ax)*r)/ (2.0 * set_sign( MaxOf( fabs(q-r), tiny ), q-r ));
		double	ulim = bx + gLimit * (cx-bx);
		double	fu;
		if ((bx-u)*(u-cx) > 0.0)
		{	fu = func( u );
			if (fu < fc)
			{	ax = bx;	bx = u;
				fa = fb;	fb = fu;
				return;
			}
			else if (fu > fb)
			{	cx = u;
				fc = fu;
				return;
			}
			u  = cx + gold * (cx-bx);
			fu = func( u );
		}
		else if ( (cx-u) * (u-ulim) > 0.0)
		{	fu = func( u );
			if (fu < fc)
			{	bx = cx;	cx = u;		u  = cx + gold * (cx-bx);
				fb = fc;	fc = fu;	fu = func( u );
		}	}
		else if ( (u-ulim) * (ulim-cx) >= 0.0)
		{	u  = ulim;
			fu = func( u );
		}
		else
		{	u  = cx + gold * (cx-bx);
			fu = func( u );
		}
		ax = bx;	bx = cx;	cx = u;
		fa = fb;	fb = fc;	fc = fu;
	}
}


//#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


double	MyMath::brent( double ax, double bx, double cx, doubleFunction& f, double tol, double& xmin )
{
//	int iter;
//	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
//	float e=0.0;
	const int		it_max = 100;
	const double 	z_eps  = 1.0e-10;
	const double 	c_gold = 0.3819660;

	double	a = (ax < cx ? ax : cx);
	double	b = (ax > cx ? ax : cx);
	double	x = bx;
	double	v = x;
	double	w = x;
	double	fx = f(x);
	double	fv = fx;
	double	fw = fx;
	double	e = 0.0;
	double	u, fu;			//	not initialized ...
	double	d;				//	not initialized ...
	
	for (int iter=1; iter<=it_max; iter++)
	{	double	xm   = 0.5 * (a+b);
		double	tol1 = tol * fabs(x) + z_eps;
		double	tol2 = 2.0 * tol1;
		if ( fabs( x-xm ) <= (tol2 - 0.5*(b-a)) )
		{	xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{	double	r = (x-w) * (fx-fv);
			double	q = (x-v) * (fx-fw);
			double	p = (x-v) * q - (x-w) * r;
			q = 2.0*(q-r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			double	etemp = e;
			e = d;
			if ( fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x) )
				d = c_gold * (e=(x >= xm ? a-x : b-x));
			else
			{	d = p/q;
				u = x+d;
				if (u-a < tol2 || b-u < tol2)
					d = set_sign(tol1,xm-x);
		}	}
		else
		{	e = (x >= xm ? a-x : b-x);
			d = c_gold * e;
		}
		u  = ( fabs(d) >= tol1 ? x+d : x+set_sign(tol1,d) );
		fu = f(u);
		if (fu <= fx)
		{	if (u >= x)		a = x;
			else			b = x;
			v  = w;		w  = x;		x  = u;		//	SHFT(v,w,x,u)
			fv = fw;	fw = fx;	fx = fu;	//	SHFT(fv,fw,fx,fu)
		}
		else
		{	if (u < x) 		a = u;
			else 			b = u;
			if (fu <= fw || w == x)
			{	v  = w;		w  = u;
				fv = fw;	fw = fu;
			}
			else if (fu <= fv || v == x || v == w)
			{	v  = u;
				fv = fu;
			}
		}
	}
	ThrowError((char *)"Too many iterations in brent" );
	xmin = x;
	return fx;
};



void	MyMath::avevar( const double_vector1& data, double& ave, double& var )
{	
/*	
	unsigned long	n = data.length();
	for (unsigned long j=1, ave=0.0; j<=n; j++)
		ave += data[j];
	ave /= n;
	double	ep = 0.0;
	for (unsigned long j=1, var=0.0; j<=n; j++)
	{	double	s = data[j] - ave;
		ep += s;
		var += s*s;
	}
	var = (var - ep*ep/n)/(n-1);
*/
	ave = data.mean();
	var = data.variance( &ave );
}


//#define ITMAX 100
//#define EPS 3.0e-8


double	MyMath::zbrent( doubleFunction& func, double x1, double x2, double tol )
{
	const int		it_max 	= 100;
	const double 	eps  	= 3.0e-8;
	
	double	a = x1;
	double	b = x2;
	double	c = x2;
	double	fa = func(a);
	double	fb = func(b);
	
//	int iter;
	//float min1,min2;
//	float fc,p,q,r,s,tol1,xm;
	double	d, e;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		ThrowError((char *)"Root must be bracketed in zbrent" );
		
	double	fc = fb;
	
	for (int iter=1; iter<=it_max; iter++)
	{	if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		{	c  = a;
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb))
		{	a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		
		double	tol1 = 2.0 * eps * fabs(b) + 0.5*tol;
		double	xm = 0.5*(c-b);
		
		if ( fabs(xm) <= tol1 || fb == 0.0 )
			return b;
			
		if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) )
		{	double	p, q;
			double	s = fb/fa;
			if (a == c)
			{	p = 2.0*xm*s;
				q = 1.0-s;
			}
			else
			{	double	r = fb/fc;
				q = fa/fc;
				p = s * ( 2.0*xm*q*(q-r) - (b-a)*(r-1.0) );
				q = (q-1.0) * (r-1.0) * (s-1.0);
			}
			if (p > 0.0)	q = -q;
			p = fabs(p);
			
			double	min1 = 3.0*xm*q-fabs(tol1*q);
			double	min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2))
			{	e = d;
				d = p/q;
			}
			else
			{	d = xm;
				e = d;
			}
		}
		else
		{	d = xm;
			e = d;
		}
		a  = b;
		fa = fb;
		
		if (fabs(d) > tol1)		b += d;
		else					b += set_sign( tol1, xm );
		
		fb = func(b);
	}
	ThrowError((char *)"Maximum number of iterations exceeded in zbrent" );
	return 0.0;
}



#pragma mark -
#pragma mark //	statistical functions


double	MyMath::gammln( double xx )
{
	static const double cof[6] = { 76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155,
		0.1208650973866179e-2, -0.5395239384953e-5 };
	
	double	x = xx;
	double	y = x;
	double	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	double	ser = 1.000000000190015;
	for (int j=0; j<=5; j++)	ser += cof[j]/++y;
	
	return	( -tmp + log(2.5066282746310005*ser/x) );
}

double	MyMath::gammp( double a, double x )
{	if (x < 0.0 || a <= 0.0) 	ThrowError((char *)"Invalid arguments in routine gammp");
	if (x < (a+1.0))
	{	double	gamser, gln;
		gser( gamser, a, x, gln );
		return gamser;
	}
	else
	{	double	gammcf, gln;
		gcf( gammcf, a, x, gln );
		return (1.0 - gammcf);
	}
}

double	MyMath::gammq( double a, double x )
{	if (x < 0.0 || a <= 0.0) 	ThrowError((char *)"Invalid arguments in routine gammq");
	if (x < (a+1.0))
	{	double	gamser, gln;
		gser( gamser, a, x, gln );
		return	(1.0 - gamser);
	}
	else
	{	double	gammcf, gln;
		gcf( gammcf, a, x, gln );
		return	gammcf;
	}
}




void	MyMath::gser( double& gamser, double a, double x, double& gln )
{
//	int n;
//	float sum,del,ap;

	gln =  gammln( a );
	
	if (x <= 0.0)
	{	if (x < 0.0)	ThrowError((char *)"x less than 0 in routine gser");
		gamser = 0.0;
		return;
	}
	else
	{	double	ap  = a;
		double	sum = 1.0/a;
		double	del = sum;
		
		const int		it_max 	= 400;		//	default was 100, but unsufficient in AnalySeries
		const double 	eps  	= 3.0e-7;
		
		for (int n=1; n<=it_max; n++)
		{	++ap;
			del *= x/ap;
			sum += del;
			if ( fabs(del) < fabs(sum)*eps )
			{	gamser = sum * exp( -x + a*log(x) - gln );
				return;
			}
		}
		ThrowError((char *)"a too large, ITMAX too small in routine gser");
		return;
	}
}



void	MyMath::gcf( double& gammcf, double a, double x, double& gln )
{
	const int		it_max 	= 100;
	const double 	eps  	= 3.0e-7;
	const double 	fpmin  	= 1.0e-30;
	

	gln = gammln( a );
	
	double	b = x+1.0-a;
	double	c = 1.0/fpmin;
	double	d = 1.0/b;
	double	h = d;
	
	int i;
	for (i=1; i<=it_max; i++)
	{	double	an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < fpmin) d = fpmin;
		c = b+an/c;
		if (fabs(c) < fpmin) c = fpmin;
		d = 1.0/d;
		double	del = d*c;
		h *= del;
		if (fabs(del-1.0) < eps) break;
	}
	if (i > it_max)
		ThrowError((char *)"a too large, ITMAX too small in gcf" );
		
	gammcf = exp( -x + a*log(x) - gln ) * h;
}

double	MyMath::erfcc( double x )
{	double	z = fabs(x);
	double	t = 1.0/(1.0+0.5*z);
	double	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return (x >= 0.0 ? ans : 2.0-ans);
}

double	MyMath::erff( double x )
{	return ( x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x) );
}

double	MyMath::erffc( double x )
{	return ( x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x) );
}

double	MyMath::betacf( double a, double b, double x)
{
//	int m,m2;
//	float aa,c,d,del,h,qab,qam,qap;
	const double fpmin = 1.0e-30;
	const double eps   = 3.0e-7;
	const int	 maxit = 100;

	double	qab = a+b;
	double	qap = a+1.0;
	double	qam = a-1.0;
	double	c = 1.0;
	double	d = 1.0-qab*x/qap;
	if (fabs(d) < fpmin) d=fpmin;
	d = 1.0/d;
	double	h = d;
	int	m;
	for (m=1 ;m<=maxit; m++)
	{	int m2 = 2*m;
		double	aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0+aa*d;
		if (fabs(d) < fpmin) d=fpmin;
		c = 1.0+aa/c;
		if (fabs(c) < fpmin) c=fpmin;
		d = 1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0+aa*d;
		if (fabs(d) < fpmin) d=fpmin;
		c = 1.0+aa/c;
		if (fabs(c) < fpmin) c=fpmin;
		d = 1.0/d;
		double	del = d*c;
		h *= del;
		if (fabs(del-1.0) < eps) break;
	}
	if (m > maxit) ThrowError((char *)"a or b too big, or MAXIT too small in betacf");
	return h;
}

double	MyMath::betai( double a, double b, double x )
{	if (x < 0.0 || x > 1.0) ThrowError((char *)"Bad x in routine betai" );
	
	double	bt;
	if (x == 0.0 || x == 1.0)	bt = 0.0;
	else						bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
		
	if (x < (a+1.0)/(a+b+2.0))	return bt*betacf(a,b,x)/a;
	else						return 1.0-bt*betacf(b,a,1.0-x)/b;
}


long				MyMath::randomGenerator::idum;
unsigned long	MyMath::randomGenerator::udum;

void	MyMath::randomGenerator::initialize()
{	//	initialisation
	clock_t	c = std::clock();					//	time units since program startup
	long	c50 = c / (CLOCKS_PER_SEC/50);		//	c in 1/50 secs units
	time_t	t = std::time(NULL);													//	real time (day, month...)
//	long	t1  = t.tm_sec + 60 * (t.tm_min + 60 * (t.tm_hour + 24 * t.tm_yday));	//	nb of secs since January 1st.
	idum = c50 + t;
	udum = idum;
}

double	MyMath::randomGenerator0::random()
{	const long		IA = 16807;
	const long 		IM = 2147483647;
	const double	AM = (1.0/IM);
	const long		IQ = 127773;
	const long		IR = 2836;
	const long		MASK = 123459876;

	idum	^= MASK;
	long k	= idum/IQ;
	idum	= IA*(idum-k*IQ) - IR*k;
	if (idum < 0) idum += IM;
	double ans	= AM*(idum);
	idum 		^= MASK;
	return ans;
}

double	MyMath::randomGenerator1::random()
{	const long		IA = 16807;
	const long 		IM = 2147483647;
	const double	AM = (1.0/IM);
	const long 		IQ = 127773;
	const long 		IR = 2836;
	const long 		NTAB = 32;
	const long 		NDIV = (1+(IM-1)/NTAB);
	const double	EPS = 1.2e-7;
	const double	RNMX = (1.0-EPS);

	static long iy = 0;
	static long iv[NTAB];

	if (idum <= 0 || !iy)
	{	if (-idum < 1)	idum = 1;
		else 			idum = -idum;
		for (int j=NTAB+7; j>=0; j--)
		{	long	k = idum/IQ;
			idum = IA*(idum-k*IQ) - IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	long	k = idum/IQ;
	idum = IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	int	  j = iy/NDIV;
	iy 		= iv[j];
	iv[j] 	= idum;
	double	temp = AM*iy;
	if (temp > RNMX)	return RNMX;
	else				return temp;
}

double	MyMath::randomGenerator2::random()
{	const long		IM1 = 2147483563;
	const long 		IM2 = 2147483399;
	const double	AM = (1.0/IM1);
	const long 		IMM1 = (IM1-1);
	const long 		IA1 = 40014;
	const long 		IA2 = 40692;
	const long 		IQ1 = 53668;
	const long 		IQ2 = 52774;
	const long 		IR1 = 12211;
	const long 		IR2 = 3791;
	const long 		NTAB = 32;
	const long 		NDIV = (1+IMM1/NTAB);
	const double	EPS = 1.2e-7;
	const double	RNMX = (1.0-EPS);

	static long idum2 	= 123456789;
	static long iy		= 0;
	static long iv[NTAB];

	if (idum <= 0)
	{	if (-idum < 1)	idum = 1;
		else			idum = -idum;
		idum2 = idum;
		for (int j= NTAB+7; j>=0; j--)
		{	long k = idum/IQ1;
			idum = IA1*(idum-k*IQ1) - k*IR1;
			if (idum < 0)	idum += IM1;
			if (j < NTAB) 	iv[j] = idum;
		}
		iy = iv[0];
	}
	long k = idum/IQ1;
	idum = IA1*(idum-k*IQ1) - k*IR1;
	if (idum < 0)	idum += IM1;
	k = idum2/IQ2;
	idum2 	= IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	int j	= iy/NDIV;
	iy	  = iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	double	temp = AM*iy;
	if (temp > RNMX) 	return RNMX;
	else 				return temp;
}


double	MyMath::randomGenerator3::random()
{	const long		MBIG = 1000000000;
	const long 		MSEED = 161803398;
	const long 		MZ = 0;
	const double	FAC = (1.0/MBIG);

	static int	inext,inextp;
	static long	ma[56];
	static int 	iff = 0;
	long 		mj;

	if (idum < 0 || iff == 0)
	{	iff = 1;
		mj = MSEED-(idum < 0 ? -idum : idum);
		mj %= MBIG;
		ma[55] = mj;
		long mk = 1;
		for (int i=1; i<=54; i++)
		{	int ii = (21*i) % 55;
			ma[ii] = mk;
			mk = mj-mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for (int k=1; k<=4; k++)
			for (int i=1; i<=55; i++)
			{	ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext	= 0;
		inextp	= 31;
		idum	= 1;
	}
	if (++inext == 56)	inext  = 1;
	if (++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return	(mj*FAC);
}


void	MyMath::randomGenerator4::psdes( unsigned long& lword, unsigned long& irword )
{	const size_t	NITER = 4;
	unsigned long	itmph = 0;
	unsigned long	itmpl = 0;
	static unsigned long c1[NITER] = { 0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L };
	static unsigned long c2[NITER] = { 0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L };

	for (unsigned long i=0; i<NITER; i++)
	{	unsigned long iswap	= irword;
		unsigned long ia	= iswap ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		unsigned long ib	= itmpl*itmpl + ~(itmph*itmph);
		irword	= lword ^ (((ia = (ib >> 16) | ((ib & 0xffff) << 16)) ^ c2[i]) + itmpl*itmph );
		lword	= iswap;
	}
}

bool	MyMath::randomGenerator4::test_psdes()
{	bool	ok = true;
	unsigned long 	lword	= 1;
	unsigned long 	irword	= 1;
	psdes( lword, irword );
	if ( lword != 0x604d1dceL || irword != 0x509c0c23L )	ok = false;
	lword	= 1;
	irword	= 99;
	psdes( lword, irword );
	if ( lword != 0xd97f8571L || irword != 0xa66cb41aL )	ok = false;
	lword	= 99;
	irword	= 1;
	psdes( lword, irword );
	if ( lword != 0x7822309dL || irword != 0x64300984L )	ok = false;
	lword	= 99;
	irword	= 99;
	psdes( lword, irword );
	if ( lword != 0xd7f376f0L || irword != 0x59ba89ebL )	ok = false;
	return ok;
};

double	MyMath::randomGenerator4::random()
{	static long 	idums = 0;
	
#if 0			// defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static const unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else			//	for IEEE 32-bit reals
	static const unsigned long jflone = 0x3f800000;
	static const unsigned long jflmsk = 0x007fffff;
#endif

	if (idum < 0) {
		idums = -idum;
		idum=1;
	}
	unsigned long 	irword	= idum;
	unsigned long 	lword	= idums;
	psdes( lword, irword );
	unsigned long 	itemp	= jflone | (jflmsk & irword);
	++idum;
	
	return (double)	((*(float*)&itemp) - 1.0);
}

bool	MyMath::randomGeneratorQD2::test_qd2()
{	bool	ok = true;
	udum = 0L;
	random();	if (udum != 0x3c6ef35f)	ok = false;
	random();	if (udum != 0x47502932)	ok = false;
	random();	if (udum != 0xd1ccf6e9)	ok = false;
	random();	if (udum != 0xaaf95334)	ok = false;
	random();	if (udum != 0x6252e503)	ok = false;
	random();	if (udum != 0x9f2ec686)	ok = false;
	random();	if (udum != 0x57fe6c2d)	ok = false;
	random();	if (udum != 0xa3d95fa8)	ok = false;
	random();	if (udum != 0x81fdbee7)	ok = false;
	random();	if (udum != 0x94f0af1a)	ok = false;
	random();	if (udum != 0xcbf633b1)	ok = false;
	return ok;
};

double	MyMath::randomGeneratorQD2::random()
{	
#if 0			// defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static const unsigned long jflone = 0x00004080;
	static const unsigned long jflmsk = 0xffff007f;
#else			//	for IEEE 32-bit reals
	static const unsigned long jflone = 0x3f800000;
	static const unsigned long jflmsk = 0x007fffff;
#endif

	udum = 1664525L * udum + 1013904223L;
	unsigned long 	itemp	= jflone | (jflmsk & udum);
	return (double)	((*(float*)&itemp) - 1.0);
};

double	MyMath::randomGenerator::expdev()
{	double	dum;
	do		dum = random();
	while	(dum == 0.0);
	return	-log(dum);
}

double	MyMath::randomGenerator::doubleexpdev()
{	double	dum;
	do		dum = random();
	while	(dum == 0.0);
	double	s = random();
	if	(s>.5)	return	log(dum);
	else		return	-log(dum);
}

double	MyMath::randomGenerator::lorentzdev()
{	double	dum;
	do		dum = random();
	while	(dum == 0.5);
	return	tan(Pi*dum);
}


double	MyMath::randomGenerator::gasdev()
{	static int 		iset = 0;
	static double 	gset;

	if  (iset == 0)
	{	double fac, rsq, v1, v2;
		do {	v1 = 2.0*random()-1.0;
				v2 = 2.0*random()-1.0;
				rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac		= sqrt(-2.0*log(rsq)/rsq);
		gset	= v1*fac;
		iset	= 1;
		return	v2*fac;
	}
	else
	{	iset = 0;
		return gset;
	}
}


void MyMath::sort( double_vector1& arr )
{	const size_t M = 7;
	const size_t NSTACK = 50;

	size_t	n 		= arr.length();
	size_t	ir 		= n;
	size_t	l 		= 1; 
	int		jstack 	= 0;

	int_vector1 istack( NSTACK );
	
	for (;;)
	{	if (ir-l < M)
		{	for (size_t j=l+1; j<=ir; j++)
			{	double	a = arr[j];
				size_t	i;
				for (i=j-1; i>=1; i--)
				{	if (arr[i] <= a) break;
					arr[i+1] = arr[i];
				}
				arr[i+1] = a;
			}
			if (jstack == 0) break;
			ir = istack[jstack--];
			l  = istack[jstack--];
		}
		else
		{	size_t	k = (l+ir) >> 1;
			swap( arr[k], arr[l+1] );
			if (arr[l+1] > arr[ir])		swap(arr[l+1],arr[ir]);
			if (arr[l] > arr[ir])		swap(arr[l],arr[ir]);
			if (arr[l+1] > arr[l])		swap(arr[l+1],arr[l]);
			size_t	i = l+1;
			size_t	j = ir;
			double	a = arr[l];
			for (;;)
			{	do i++; 	while (arr[i] < a);
				do j--; 	while (arr[j] > a);
				if (j < i) 	break;
				swap(arr[i],arr[j]);
			}
			arr[l] = arr[j];
			arr[j] = a;
			jstack += 2;
			if (jstack > NSTACK)	ThrowError((char *)"NSTACK too small in sort.");
			if (ir-i+1 >= j-l)
			{	istack[jstack]	 = ir;
				istack[jstack-1] = i;
				ir				 = j-1;
			}
			else
			{	istack[jstack]	 = j-1;
				istack[jstack-1] = l;
				l				 = i;
			}
		}
	}
}


void MyMath::sort2( double_vector1& arr, double_vector1& brr )
{
	const int M = 7;
	const int NSTACK = 50;

	size_t	n 		= arr.length();			//	= brr.length();
	size_t	ir 		= n;
	size_t	l 		= 1; 
	int		jstack 	= 0;

	int_vector1 istack( NSTACK );
	
	for (;;)
	{	if (ir-l < M)
		{	for (size_t j=l+1; j<=ir; j++)
			{	double	a = arr[j];
				double	b = brr[j];
				size_t	i;
				for (i=j-1; i>=1; i--)
				{	if (arr[i] <= a) break;
					arr[i+1] = arr[i];
					brr[i+1] = brr[i];
				}
				arr[i+1] = a;
				brr[i+1] = b;
			}
			if (jstack == 0)	return;
			ir = istack[jstack--];
			l  = istack[jstack--];
		}
		else
		{	size_t	k = (l+ir) >> 1;
			swap( arr[k], arr[l+1] );
			swap( brr[k], brr[l+1] );
			if (arr[l+1] > arr[ir])
			{	swap(arr[l+1],arr[ir]);
				swap(brr[l+1],brr[ir]);
			}
			if (arr[l] > arr[ir])
			{	swap(arr[l],arr[ir]);
				swap(brr[l],brr[ir]);
			}
			if (arr[l+1] > arr[l])
			{	swap(arr[l+1],arr[l]);
				swap(brr[l+1],brr[l]);
			}
			size_t	i = l+1;
			size_t	j = ir;
			double	a = arr[l];
			double	b = brr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				swap(arr[i],arr[j]);
				swap(brr[i],brr[j]);
			}
			arr[l] = arr[j];
			arr[j] = a;
			brr[l] = brr[j];
			brr[j] = b;
			jstack += 2;
			if (jstack > NSTACK) 	ThrowError((char *)"NSTACK too small in sort2.");
			if (ir-i+1 >= j-l)
			{	istack[jstack]	 = ir;
				istack[jstack-1] = i;
				ir				 = j-1;
			}
			else
			{	istack[jstack]	 = j-1;
				istack[jstack-1] = l;
				l				 = i;
			}
		}
	}
}


double	MyMath::spear( const double_vector1&  data1, const double_vector1& data2,
	/*float *d, float *zd, float *probd, */ double& probrs )
{
	size_t	n = data1.length();
	if ( data2.length() != n )	throw GenericError( "Pb in MyMath::spear" );
	double_vector1 wksp1( data1 );
	double_vector1 wksp2( data2 );
	MyMath::sort2( wksp1, wksp2 );
	double	sf = MyMath::crank( wksp1 );
	MyMath::sort2( wksp2, wksp1 );
	double	sg = MyMath::crank( wksp2 );
	double	d = 0.0;
	for (size_t j=1; j<=n; j++)	d += MyMath::sqr( wksp1[j]-wksp2[j] );
	double	en = n;
	double	en3n = en*en*en - en;
//	double	aved = en3n/6.0 - (sf+sg)/12.0;
	double	fac = (1.0-sf/en3n)*(1.0-sg/en3n);
//	double	vard = ((en-1.0)*en*en * MyMath::sqr(en+1.0)/36.0)*fac;
//	double	zd = (d-aved)/std::sqrt(vard);
//	probd = MyMath::erffc( fabs(zd)/MyMath::Sqrt2 );
	double	rs = (1.0-(6.0/en3n)*(d+(sf+sg)/12.0))/std::sqrt(fac);
	fac = (rs+1.0)*(1.0-rs);
	if (fac > 0.0)
	{	double	df = en-2.0;
		double	t = rs*sqrt( df/fac );
		probrs = MyMath::betai(0.5*df,0.5,df/(df+t*t));
	}
	else	probrs = 0.0;

	return rs;
}

double MyMath::crank( double_vector1& w )
{
	size_t j=1;
	size_t n = w.size();

	double s = 0.0;
	while (j < n) {
		if (w[j+1] != w[j])		{	w[j]=j;		++j;	}
		else
		{	size_t jt;
			for (jt = j+1; (jt<=n && w[jt]==w[j]); jt++)	{};
			double rank = 0.5*(j+jt-1);
			for (size_t ji = j; ji<=(jt-1); ji++)	w[ji] = rank;
			double t = jt-j;
			s += t*t*t-t;
			j=jt;
	}	}
	if (j == n) w[n]=n;
	return s;
}




double	MyMath::LinearCorrelation( const double_vector1& x1, const double_vector1& y1, const double_vector1& x2, const double_vector1& y2 )
{	//double	xa = std::max( x1.first(), x2.first());
	//double	xb = std::min( x1.last(),  x2.last() );
	double_vector	full_time( x1 );
	full_time.insert_vector( x2 );						//	build a common x-scale, by default with a common range
	
	MyMath::lin_interpolfunc interp2( x2, y2 );			//	linear interpolation on the new common scale
	MyMath::lin_interpolfunc interp1( x1, y1 );
	double_vector	full_2( full_time, interp2 );
	double_vector	full_1( full_time, interp1 );

	MyMath::lin_interpolfunc integ2( full_time, full_2 );		//	prepare integration
	MyMath::lin_interpolfunc integ1( full_time, full_1 );		//	prepare integration
	double	mean2 = integ2.MeanBetween( full_time(1), full_time(full_time.length()) );		//	compute and remove the mean (note: it is not evenly sampled)
	double	mean1 = integ1.MeanBetween( full_time(1), full_time(full_time.length()) );
	full_2 -= mean2;
	full_1 -= mean1;
	
	double_vector	full_2_sqr( full_2, MyMath::sqr );
	double_vector	full_1_sqr( full_1, MyMath::sqr );
	double_vector	full_cross( full_time.length() );	full_cross.set_to( full_2, full_1, MyMath::multiplicate );
	MyMath::lin_interpolfunc s2_inter( full_time, full_2_sqr );		//	prepare integration
	MyMath::lin_interpolfunc s1_inter( full_time, full_1_sqr );		//	prepare integration
	MyMath::lin_interpolfunc cr_inter( full_time, full_cross );		//	prepare integration
	double	s2 = s2_inter.IntegrateBetween( full_time(1), full_time(full_time.length()) );
	double	s1 = s1_inter.IntegrateBetween( full_time(1), full_time(full_time.length()) );
	double	cr = cr_inter.IntegrateBetween( full_time(1), full_time(full_time.length()) );
	return	cr / sqrt(s2 * s1);
};

