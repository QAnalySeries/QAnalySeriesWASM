



#include "myMath.h"



class ssa_autocovariance : public double_matrix {
private:
	//	Divers calculs lies a la matrice de covariance 
	//***********************************************************
	// Calcul de la matrice de covariance						
	// avec la methode de V & G 1989							
	//**********************************************************

	void	autocovariance_VG( const double_vector& ordo )
	{	size_t			ndata = ordo.size();
		size_t			emdim = dim1();
		double_vector1	autoc( emdim );
		int		nn = ndata - emdim + 1;
		for (int i=1; i<=emdim; i++)
		{	autoc[i] = 0;
			for (int j=1; j<=(ndata-i+1); j++)
				autoc[i] += ordo[j] * ordo[j+i-1];
			autoc[i] /= (ndata - i);
			for (int j=1; j<=(emdim-i+1); j++)
			{	(*this)(j,j+i-1) = autoc[i];
				(*this)(j+i-1,j) = autoc[i];
			};
		};
	};
		
		
	//***********************************************************
	// Calcul de la matrice de covariance						 
	//* avec la methode de B & K (1986)							 
	//************************************************************
	void	autocovariance_BK( const double_vector& ordo )
	{	size_t			ndata = ordo.size();
		size_t			emdim = dim1();
		int		nn = ndata - emdim + 1;
		for (int i=1; i<=emdim; i++)
		{	for (int j=i; j<=emdim; j++)
			{	(*this)(i,j) = produit( ordo, i-1, ordo, j-1, nn) / nn;
				if ( j > i )	(*this)(j,i) = (*this)(i,j);
			}
		};
	};

	//	Produit scalaire de deux vecteurs x et y de taille n	

	double	produit( const double_vector& x, int ix, const double_vector& y, int iy, int n )
	{	double	sum = 0;
		for (int i=1; i<=n; i++)	sum += x[i+ix] * y[i+iy];
		return sum;
	};

public:
	typedef enum { VGtype, BKtype } covariance_type;
	
	ssa_autocovariance( size_t emdim, covariance_type t, const double_vector& ordo ) : double_matrix( emdim, emdim )
	{	switch (t) {
			case VGtype:	autocovariance_VG( ordo );	break;
			case BKtype:	autocovariance_BK( ordo );	break;
		}
	}
};


class	vector_vector: public val1Darray<double_vector,1> {
public:
	vector_vector( size_t nv, size_t d ) : val1Darray<double_vector,1>( nv )
	{	for (size_t i=1; i<=nv; i++)	(*this)(i) = double_vector(d);
	}
};

class	ssa : public double_vector {

//Boolean	SSAspec( double* ordo, int ndata, int emdim, int* S, double* eigenVal, double** vep,
//			double** pc, int* index, int nsauv, ssa_autocovariance::covariance_type nmode )
private:
	double_matrix	vep;
public:
	double_vector	eigenValues;
	int				StatDim;
	vector_vector	pc;
	
//	double_vector	recons;				-> *this
//		vep		= NR_dmatrix( 1, emdim, 1, npc);
//		pc		= NR_dmatrix( 1, ndata-emdim+1, 1, npc);
		
public:
	ssa	( /*const double_vector& absi,*/ const double_vector& ordo, size_t emdim,
			const int_vector& index, ssa_autocovariance::covariance_type nmode ): 
		double_vector(ordo.size()),
		//	, /*double* recons,*/ double_matrix& vep, val1Darray<double_vector&,1>& pc, int	S ):
		eigenValues( emdim ), vep( emdim, index.size() ), pc( index.size(), ordo.size()-emdim+1 )
	{
		size_t	ndata = ordo.size();
		size_t	npc = index.size();			//	number of PCs
//		int		nn = ndata - emdim + 1;
		
//	Calcul des npc premieres EOF 
		SSAspec( ordo, emdim, index, npc, nmode );
	//	SSAsave( eigenValues, vep, pc, ndata, emdim, index, npc );
		double_vector	Trend( ndata );
		reconst( Trend, emdim, ndata, npc );
		
		//double	aver = avera( ordo, 0, ndata );
		double	aver = ordo.mean();
		saveres( /*absi,*/ ordo, Trend, aver, ndata );
	};
	
	
	
private:
	void	SSAspec( const double_vector& ordo, int emdim, const int_vector& index, size_t nsauv, ssa_autocovariance::covariance_type nmode )
	{
		ssa_autocovariance	autoco( emdim, nmode, ordo );
		double_vector1	offdiag( emdim );
		MyMath::tred2( autoco, eigenValues, offdiag );			//	mise sous forme tridiagonale de la matrice de covariance 
		MyMath::tqli( eigenValues, offdiag, autoco );			//	calcul des valeurs propres et vecteurs propres de la matrice reduite
		MyMath::eigsrt( eigenValues, autoco );					//	rearangement des valeurs propres et vecteurs propres dans un ordre monotone. Ca aide !	
		
					//	Remplissage des EOF	a sauvegarder
		for (int i=1; i<=emdim; i++)
		for (int j=1; j<=nsauv; j++)
			vep(i,j) = autoco(i,index[j]);
		
		pc_calc( ordo, nsauv );
		
		StatDim = S_calc();
		
	};
		
	// Programme rebuild 
	// Reconstruction d'un signal a partir d'un certain nombre de PC 
	// D'apres la methode de Robert Vautard Janvier 1991
	// Sauvegarde de la serie reconstituee	

	void	saveres( /*const double_vector& X,*/ const double_vector& Y, const double_vector& Ytrend, double aver, int N )
	{	for (int i=1; i<=N; i++)
			(*this)[i] = Ytrend[i] + aver;
	};


	// Reconstruction d'une serie a partir des irho PC selectionnees

	void	reconst( double_vector& X, int M, int N, int irho )
	{
		for (int i=1; i<=N; i++)
		{	X[i] = 0;
			int		imin = MyMath::MinOf(i, M);
			int		imax = MyMath::MaxOf(i - N + M, 1);
			for (int k=1; k<=irho; k++)
			{	double	xjunk = 0;
				for (int j=imax; j<=imin; j++)
					xjunk += pc[k][i-j+1] * vep(j,k);
				X[i] += xjunk;
			};
			X[i] /= (imin - imax + 1);
		};
	};
		


	void	pc_calc( const double_vector& ordo, int nsauv )
	//************************************************************************
	//* calcul et sauvegarde des PCs											
	//************************************************************************
	{	int emdim = vep.dim1();
		int	nn = ordo.size() - emdim + 1;
		for (int i=0; i<=(nn-1); i++)
		for (int j=1; j<=nsauv; j++)
		{	pc[j][i+1] = 0;
			for (int k=1; k<=emdim; k++)
				pc[j][i+1] += ordo[k+i] * vep(k,j);
		};
	};
		

	int	S_calc()
	//***************************************
	//	Calcul du nombre de PC semblables 
	//	a du bruit blanc a 95% de confiance	
	//***************************************
	{	int emdim = eigenValues.size();
		int	i = emdim - 1;
		double	eigenValmin, eigenValmax, stdevl;
		do
		{	stdevl = 0;
			double	aver = avera( eigenValues, i-1, emdim-i+1 );
			for (int j=i; j<=emdim; j++)	stdevl += MyMath::sqr(eigenValues[j] - aver);
			stdevl = sqrt( stdevl / (emdim-i) );
			eigenValmin = eigenValues[i] - 1.96 * stdevl;
			eigenValmax = eigenValues[emdim] + 1.96 * stdevl;
			i--;
		}
		while ( (i>1) && (eigenValmin>eigenValmax) && (eigenValues[i-1] < eigenValues[i] + 1.96 * stdevl) );
		return i;
	};

	double	avera( double_vector& x, int i0, int n )
	//************************************************************
	// Soustraction de la moyenne du signal						 
	//************************************************************
	{	double	moy = 0;
		for (int i=1; i<=n; i++)	moy += x(i+i0);
		moy /= n;
		for (int i=1; i<=n; i++)	x(i+i0) -= moy;
		return moy;
	};

};
	
