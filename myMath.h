

#pragma once

#ifndef	__MyMath_H_
#define	__MyMath_H_

#include <string>

#ifdef __MWERKS__
#include <valarray.h>
#else
#include <valarray>
#endif

// forward declaration ....
/*
namespace MyMath {
using namespace std;
class a_double_vector1;
inline double	bissectbracket( const a_double_vector1& x, const double& value, unsigned long& klo );
};
*/

typedef bool Boolean;
typedef unsigned int FourCharCode;

class	GenericError {
public:
	std::string		txt;
	GenericError( const char* s )	: txt(s)	{};
	GenericError( const std::string& s ) : txt(s)	{};
};

 
template <class T> 
class	TFunction {
protected:
	virtual T	func( const T ) const = 0;
public:
	virtual T	operator()( const T x )	{	return func(x);	};
    virtual ~TFunction<T>()	{};
};
template <class T> 
class	TFunction2 {
protected:
	virtual T	func( const T, const T ) const = 0;
public:
	virtual T	operator()( const T x1, const T x2 )	{	return func(x1,x2);	};
    virtual ~TFunction2<T>()	{};
};
template <class T> 
class	TFunction3 {
protected:
	virtual T	func( const T, const T, const T ) const = 0;
public:
	virtual T	operator()( const T x1, const T x2, const T x3 )	{	return func(x1,x2,x3);	};
    virtual ~TFunction3<T>()	{};
};


template <class T>
class	TStdFunction : public TFunction<T> {
protected:
	T (*f)(T);
	virtual double	func( const T x )	const {	return	(*f)(x);	};
public:
	TStdFunction( T (*f0)(T) ):f(f0)	{}; 
};
template <class T>
class	TStdFunction2 : public TFunction2<T> {
protected:
	T (*f)(T,T);
	virtual double	func( const T x1, const T x2 )	const {	return	(*f)(x1,x2);	};
public:
	TStdFunction2<T>( T (*f0)(T,T) ):f(f0)	{}; 
};
template <class T>
class	TStdFunction3 : public TFunction3<T> {
protected:
	T (*f)(T,T,T);
	virtual double	func( const T x1, const T x2, const T x3 )	const {	return	(*f)(x1,x2,x3);	};
public:
	TStdFunction3<T>( T (*f0)(T,T,T) ):f(f0)	{}; 
};


//	some functional types

typedef	float	(*floatFunc)(float);
typedef	float	(*floatFunc2)(float, float);
typedef	float	(*floatFunc3)(float, float, float);

typedef	double	(*doubleFunc)(double);
typedef	double	(*doubleFunc2)(double, double);
typedef	double	(*doubleFunc3)(double, double, double);

typedef	TFunction<double>	doubleFunction;
typedef	TFunction2<double>	doubleFunction2;
typedef	TFunction3<double>	doubleFunction3;
/*
class	doubleStdFunction : public doubleFunction {
protected:
	doubleFunc		f;
	virtual double	func( const double x )	const {	return	(*f)(x);	};
public:
	doubleStdFunction( doubleFunc f0 ):f(f0)	{}; 
};
class	doubleStdFunction2 : public doubleFunction2 {
protected:
	doubleFunc2		f;
	virtual double	func( const double x1, const double x2 )	const {	return	(*f)(x1,x2);	};
public:
	doubleStdFunction2( doubleFunc2 f0 ):f(f0)	{}; 
};
class	doubleStdFunction3 : public doubleFunction3 {
protected:
	doubleFunc3		f;
	virtual double	func( const double x1, const double x2, const double x3 )	const {	return	(*f)(x1,x2,x3);	};
public:
	doubleStdFunction3( doubleFunc3 f0 ):f(f0)	{}; 
};
*/
typedef	TStdFunction<double>	doubleStdFunction;
typedef	TStdFunction2<double>	doubleStdFunction2;
typedef	TStdFunction3<double>	doubleStdFunction3;




template <class T> 
class abstract_valarray_0 {											//	access by pointers only
public:    
	virtual	T*			firstPtr() 	const	= 0;
	virtual	T*			lastPtr()  	const	= 0;
	virtual	size_t		length()	const	= 0;
	
		//	to be supressed
//	enum	copyFlag_enum		{ dontCopy = true	};				//	only one value. Just a dummy argument for some functions
//	typedef	typename	abstract_valarray_0<T>::copyFlag_enum	copyFlag;
	
	
	size_t	size()	const	{	return	length();	};
	
	T	first() 	const		{	return	*firstPtr();	};
	T	last()  	const		{	return	*lastPtr();		};
	
	void	set_to( const T* t )						{	std::copy( t, t+length(), firstPtr() );				};
    void	get_into( T* t )							{	std::copy( firstPtr(), firstPtr()+length(), t );	};
	void	set_to( const abstract_valarray_0<T>& v )	{	set_to( v.firstPtr() );								};	//	does not change size !!!
	
	void	set_all_to( const T& t )		{	std::fill( firstPtr(), firstPtr()+length(), t );	};

    void	set_to( const std::valarray<T>& v, T (*func)(T) )
		{	TStdFunction<T> f(func);		set_to( v, f );				};
    void	set_to( const std::valarray<T>& v1, const std::valarray<T>& v2, T (*func)(T,T) )
		{	TStdFunction2<T> f(func);		set_to( v1, v2, f );		};
    void	set_to( const std::valarray<T>& v1, const std::valarray<T>& v2, const std::valarray<T>& v3, T (*func)(T,T,T) )
		{	TStdFunction3<T> f(func);		set_to( v1, v2, v3, f );	};
		
    void	set_to( const std::valarray<T>& v, TFunction<T>& func )
		{	if (length() == v.size())
			{	const T* first = &((std::valarray<T>&) v)[0];
				const T* last = first + v.size();
				T* result = firstPtr();
 				while (first != last)	*result++ = func( *first++ );
			}
		};
    void	set_to( const std::valarray<T>& v1, const std::valarray<T>& v2, TFunction2<T>& func )
		{	if (length() == v1.size() && length() == v2.size())
			{	const T*	first1 = &((std::valarray<T>&) v1)[0];	const T* last1 = first1 + v1.size();
				const T*	first2 = &((std::valarray<T>&) v2)[0];	const T* last2 = first2 + v2.size();
				T*	result = firstPtr();
 				while (first1 != last1 && first2 != last2)		*result++ = func( *first1++, *first2++ );
 			}
		};
    void	set_to( const std::valarray<T>& v1, const std::valarray<T>& v2, const std::valarray<T>& v3, TFunction3<T>& func )
		{	if (length() == v1.size() && length() == v2.size() && length() == v3.size())
			{	const T*	first1 = &((std::valarray<T>&) v1)[0];	const T* last1 = first1 + v1.size();
				const T*	first2 = &((std::valarray<T>&) v2)[0];	const T* last2 = first2 + v2.size();
				const T*	first3 = &((std::valarray<T>&) v3)[0];	const T* last3 = first3 + v3.size();
				T*	result = firstPtr();
 				while (first1 != last1 && first2 != last2 && first3 != last3)		*result++ = func( *first1++, *first2++, *first3++ );
 			}
		};
		
	typedef	enum { no, remove_mean, remove_pt_trend, remove_stat_trend } TrendFlag;	
		//	{ do nothing, remove the mean, remove trend defined by first/last points (set them to zero), remove statistical trend }
		
	void	copy( const abstract_valarray_0<T>& v )
		{	copy( v, v.length() );		};
	void	copy( const abstract_valarray_0<T>& v, size_t siz )
		{	std::copy( v.firstPtr(), v.firstPtr()+siz, firstPtr() );			};	//	does not change size !!!
	void	copy_minus_mean( const abstract_valarray_0<T>& v )
		{	T m = v.mean();		T* r = firstPtr();
	    	for (T* p = v.firstPtr(); p < v.firstPtr()+v.size(); p++)	
				*r++ = *p - m;
		}
	void	copy_minus_pt_trend( const abstract_valarray_0<T>& v )
		{	T x0 = *(v.firstPtr());		T w = (*(v.lastPtr()) - x0)/(length()-1);
			T* r = firstPtr();			T w0 = x0;
	    	for (T* p = v.firstPtr(); p < v.firstPtr()+v.size(); p++)	{	*r++ = *p - w0;		w0 += w;	}
		}
	void	do_copy( const abstract_valarray_0<T>& v, TrendFlag f )
		{	switch (f)	{
				case no:				copy( v, v.length() );		break;
				case remove_mean:		copy_minus_mean( v );		break;
				case remove_pt_trend:	copy_minus_pt_trend( v );	break;
				case remove_stat_trend:		throw GenericError( "not-implemented. remove_stat_trend" );	break;
			}
		};
	void	copy_backwards( const abstract_valarray_0<T>& v )
		{	copy_backwards( v, v.length() );		};
	void	copy_backwards( const abstract_valarray_0<T>& v, size_t siz, size_t where = 0 )
		{	T* r = v.firstPtr() + siz-1;
	    	for (T* p = firstPtr()+where; p < firstPtr()+where+siz; p++)		*p = *r--;
		}

	T	somme() const		//	sum()	defined in std::valarray
    	{	T* p = lastPtr();	T s = 0;
    		while (p>=firstPtr())		{	s += (*p);		p--;	}
    		return	s;
		}
    T	mean() const
    	{	return	somme()/length();
    	}
    T	sum_sqr() const
    	{	T* p = lastPtr();	T s = 0;
    		while (p>=firstPtr())		{	s += (*p)*(*p);		p--;	}
    		return	s;
		}
    T	sum_sqr_minus_mean( T* mm = NULL ) const		//	if mm given, dont recompute the mean
    	{	T* p = lastPtr();	T va = 0;	T e = 0;	T m = (mm ? *mm : mean());
    		while (p>=firstPtr())		{	T s = (*p)-m;	e +=s;	va += s*s;	p--;	}
    		return	(va - e*e/length());
		}
    T	variance( T* mm = NULL ) const					//	if mm given, dont recompute the mean
		{	if (length() > 1)	return	sum_sqr_minus_mean(mm)/(length() - 1);
			else	return 0;
		}
    T	cross_sum_sqr_minus_mean( const abstract_valarray_0<T>& v, T* mm = NULL, T* mmv = NULL ) const
		{	if (length() != v.length())		throw GenericError( "vectors with different size" );
			T* p = lastPtr();	T* pv = v.lastPtr();	T ss = 0;
			T m = (mm ? *mm : mean());		T mv = (mmv ? *mmv : v.mean());
    		while (p>=firstPtr())		{	ss += ((*p)-m)*((*pv)-mv);	p--;	pv--;	}
			return	ss;
		}
		
    T		min_step() const
    	{	T* p = lastPtr()-1;		T s = std::fabs(*p - *(p+1));	p--;
    		for (; p>=firstPtr(); p--)
    		{	T x = fabs(*p - *(p+1));	if (x<s)	s = x;	}
    		return	s;
		}
    T		max_step() const
    	{	T* p = lastPtr()-1;		T s = fabs(*p - *(p+1));	p--;
    		for (; p>=firstPtr(); p--)
    		{	T x = fabs(*p - *(p+1));	if (x>s)	s = x;	}
    		return	s;
		}
    void	revert()
    	{	T* p0 = firstPtr();	T* p1 = lastPtr();
    		for (; p1>p0; p1--,p0++)	std::swap( *p0, *p1 );
		}
    T		max_Func( T (*func)(T) ) const
        {	T* p = lastPtr();	T mf = func(*p--);	T m;
            while (p>=firstPtr())	if ((m=func(*p--)) > mf)	mf=m;
    		return	mf;
		}
    T		min_Func( T (*func)(T) ) const
        {	T* p = lastPtr();	T mf = func(*p--);	T m;
            while (p>=firstPtr())	if ((m=func(*p--)) < mf)	mf=m;
    		return	mf;
		}
    
	static	T	valarray0_epsilon()	{	return (T)1e-10;	};
	
    Boolean	check_same_as( const std::valarray<T>& v, T eps = valarray0_epsilon() ) const
		{	if (size() != v.size())	return false;
			T* first = firstPtr();	T* last = lastPtr();	const T* p = &((std::valarray<T>&) v)[0];
    		while (first <= last)
    			if ( std::fabs(*p++ - *first++) > eps ) return false;
    		return	true;
		}
    Boolean	check_contains( T x, T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr();
    		for (; p>=firstPtr(); p--)	if ( fabs(*p - x) <= eps) return true;
    		return	false;
		}
    Boolean	check_regular( T eps = valarray0_epsilon() ) const
    	{	if (size() < 3)	return true;
			T* p = lastPtr();	T d = *p--;	d -= *p--;	T deps = std::fabs(d*eps);	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( std::fabs(d + *p - *(p+1)) <= deps);
    		return	ok;
		}
    Boolean	check_same_regular_step_as( const abstract_valarray_0<T>& v, T eps = valarray0_epsilon() ) const
    	{	T* p1 = firstPtr();		T d1 = std::fabs(*p1 - *(p1+1));
    		T* p2 = v.firstPtr();	T d2 = std::fabs(*p2 - *(p2+1));
    		return ( check_regular(eps) && v.check_regular(eps) && std::fabs(d1-d2) <= d1*eps );
    	}
    Boolean	check_strict_increase( T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr()-1;	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( (*p - *(p+1)) < eps);
    		return	ok;
		}
    Boolean	check_increase( T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr()-1;	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( (*p - *(p+1)) <= eps);
    		return	ok;
		}
    Boolean	check_strict_decrease( T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr()-1;	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( (*p - *(p+1)) > eps);
    		return	ok;
		}
    Boolean	check_decrease( T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr()-1;	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( (*p - *(p+1)) >= eps);
    		return	ok;
		}
    Boolean	check_distinct( T eps = valarray0_epsilon() ) const
    	{	T* p = lastPtr()-1;	Boolean ok = true;
    		for (; p>=firstPtr() && ok; p--)	ok = ( fabs(*p - *(p+1)) <= eps );
    		return	ok;
		}

};

template <class T, int offsetX> 
class abstract_valarray_1D : public abstract_valarray_0<T> {			//	access by index 1D
public:
    T    operator[] (size_t pos) const	{	return	*(this->firstPtr() + pos - offsetX);	};
    T&   operator[] (size_t pos)		{	return	*(this->firstPtr() + pos - offsetX);	};
    T    operator() (size_t pos) const	{	return	*(this->firstPtr() + pos - offsetX);	};
    T&   operator() (size_t pos)		{	return	*(this->firstPtr() + pos - offsetX);	};
	
		
	void	hunt( const T& x, size_t& jlo ) const		//	d'apres Numerical Recipes 'hunt'
	{	size_t	n = this->length();
		size_t	jm,jhi,inc;
		Boolean	ascnd = (n == 1) || (this->last() > this->first());	//	(xx[n] > xx[1]);		//	for n==1, ascnd is assumed true !!
		
		if (jlo <= 0 || jlo > n)		{	jlo=0;	jhi=n+1;	}
		else {	inc=1;
				if (x >= (*this)[jlo] == ascnd) {
					if (jlo == n) return;
					jhi=(jlo)+1;
					while (x >= (*this)[jhi] == ascnd) {
						jlo=jhi;		inc += inc;		jhi=(jlo)+inc;
						if (jhi > n)	{		jhi=n+1;	break;		}
					}
			} else {
				if (jlo == 1) {		jlo=0;		return;		}
				jhi=(jlo)--;
				while (x < (*this)[jlo] == ascnd) {
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
			if (x > (*this)[jm] == ascnd)	jlo=jm;
			else							jhi=jm;
		}
	};
	
	T	bissectbracket( const T& value, size_t& klo )  const
	{	
		hunt( value, klo );						//	output: 0 <= klo <= n
		
		size_t n = this->length();
		if ( klo >= n )			klo = n-1;
		else if ( klo <= 0 )	klo = 1;			//	output: 1 <= klo <= n-1
		
		return	( (*this)[ klo+1 ] - (*this)[ klo ] );			//	always OK, except for n = 1
	};

	

};


template <class T> 
class abstract_valarray_2D_C0 : public abstract_valarray_0<T> {			//	just dimensions, pointer access in C order
protected:
    size_t	sizex, sizey;
    //	fastest index = last index (v11, v12, ...)
public:
    abstract_valarray_2D_C0() : sizex(0), sizey(0)	{};
	abstract_valarray_2D_C0( size_t nx, size_t ny ) : sizex(nx), sizey(ny)	{};
    size_t	dim1()	const	{	return	sizex;	};
    size_t	dim2()	const	{	return	sizey;	};
    size_t	pos_of(size_t posx, size_t posy)	const	{	return	posy+posx*sizey;		};
};
template <class T> 
class abstract_valarray_3D_C0 : public abstract_valarray_0<T> {			//	just dimensions, pointer access in C order
protected:
    size_t	sizex, sizey, sizez;
    //	fastest index = last index (v111, v112, ...)
public:
    abstract_valarray_3D_C0() : sizex(0), sizey(0), sizez(0)	{};
	abstract_valarray_3D_C0( size_t nx, size_t ny, size_t nz ) : sizex(nx), sizey(ny), sizez(nz)	{};
    size_t	dim1()	const	{	return	sizex;	};
    size_t	dim2()	const	{	return	sizey;	};
    size_t	dim3()	const	{	return	sizez;	};
    size_t	pos_of(size_t posx, size_t posy, size_t posz)	const	{	return	posz+sizez*(posy+posx*sizey);		};
};




template <class T, int offsetX, int offsetY> 
class abstract_valarray_2D_C : public abstract_valarray_2D_C0<T> {			//	access by index 2D
public:    
    abstract_valarray_2D_C() : abstract_valarray_2D_C0<T>()	{};
	abstract_valarray_2D_C( size_t nx, size_t ny ) : abstract_valarray_2D_C0<T>(nx,ny)	{};
    T    operator() (size_t posx, size_t posy) const	{	return	*(this->firstPtr() + abstract_valarray_2D_C0<T>::pos_of(posx-offsetX, posy-offsetY));	};
    T&   operator() (size_t posx, size_t posy)			{	return	*(this->firstPtr() + abstract_valarray_2D_C0<T>::pos_of(posx-offsetX, posy-offsetY));	};

    void	set_to_2D( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, T (*func)(T,T) )
		{	TStdFunction2<T> f(func);		set_to_2D( x, y, f );		};

    void	set_to_2D( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, TFunction2<T>& func )
		{	if (abstract_valarray_2D_C0<T>::dim1() == x.size() && abstract_valarray_2D_C0<T>::dim2() == y.size())
			{	const T* first_x = (x.firstPtr());	const T* last_x = first_x + x.size();
				T* result = this->firstPtr();
				while (first_x != last_x)
			    {	const T* first_y = (y.firstPtr());	const T* last_y = first_y + y.size();
					while (first_y != last_y)
					{	*result =func( *first_x, *first_y );	++result, ++first_y;	};
					++first_x;
			    }
			};
		};
};


template <class T, int offsetX, int offsetY, int offsetZ> 
class abstract_valarray_3D_C : public abstract_valarray_3D_C0<T> {			//	access by index 2D
public:    
    abstract_valarray_3D_C() : abstract_valarray_3D_C0<T>()	{};
	abstract_valarray_3D_C( size_t nx, size_t ny, size_t nz ) : abstract_valarray_3D_C0<T>(nx,ny,nz)	{};
    T    operator() (size_t posx, size_t posy, size_t posz) const	{	return	*(this->firstPtr() + abstract_valarray_3D_C0<T>::pos_of(posx-offsetX, posy-offsetY, posz-offsetZ));	};
    T&   operator() (size_t posx, size_t posy, size_t posz)			{	return	*(this->firstPtr() + abstract_valarray_3D_C0<T>::pos_of(posx-offsetX, posy-offsetY, posz-offsetZ));	};

    void	set_to_3D( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, T (*func)(T,T) )
		{	TStdFunction3<T> f(func);		set_to_3D( x, y, z, f );		};

    void	set_to_3D( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, TFunction3<T>& func )	
		{	if (abstract_valarray_3D_C0<T>::dim1() == x.size() && abstract_valarray_3D_C0<T>::dim2() == y.size() && abstract_valarray_3D_C0<T>::dim3() == z.size())
			{	const T* first_x = (x.firstPtr());	const T* last_x = first_x + x.size();
				T* result = this->firstPtr();
				while (first_x != last_x)
			    {	const T* first_y = (y.firstPtr());	const T* last_y = first_y + y.size();
					while (first_y != last_y)
					{	const T* first_z = (z.firstPtr());	const T* last_z = first_z + z.size();
						while (first_z != last_z)
						{	*result =func( *first_x, *first_y, *first_z );	++result, ++first_z;	};
						++first_y;
				    }
					++first_x;
			}	};
		};
};



template <class T>
class valarray0_without_ownership : public std::valarray<T> {
private:
	
	T*			_ptr0;			//	real location
	
	void	setPtr0()				{	_ptr0 = ((T*)(&(std::valarray<T>::operator[](0))));		};

protected:
	T*		ptr() 	const						{	return	_ptr0;	};
	void	change_ptr( T* t, size_t n )		{	*this = std::valarray<T>(t,n);	};
	
//public:
//	static	const T*	access_ptr( const std::valarray<T>& v, size_t i )	{	return ((&(v.operator[](i))));		};
    
public:

	virtual	~valarray0_without_ownership()	{};
      
//		Assignment
	valarray0_without_ownership&	operator= (const valarray0_without_ownership<T>& v)
		{	if (std::valarray<T>::size() != v.size())	resize(v.size());		//	bug in GNU 4.0
			(std::valarray<T>::operator=)(v);
			setPtr0();	return *this;
		};
	valarray0_without_ownership&	operator= (const std::valarray<T>& v)
		{	if (std::valarray<T>::size() != v.size())	resize(v.size());		//	bug in GNU 4.0
			(std::valarray<T>::operator=)(v);
			setPtr0();	return *this;
		};
	 
	void  resize (size_t sz, T c = T ())
		{	std::valarray<T>::resize(sz,c);		setPtr0();	};
		
//		Make new vectors
    valarray0_without_ownership ()			: std::valarray<T>()					{	setPtr0();	};
	valarray0_without_ownership (size_t n)	: std::valarray<T>( n )					{	setPtr0();	};
	valarray0_without_ownership (const T& t, size_t n) : std::valarray<T>( t, n )	{	setPtr0();	};	//	fill all with t
   	
   		

//		Make vectors from existing pointer

			//	don't make a copy.
//	valarray0_without_ownership ( T* t, size_t n, abstract_valarray_0<T>::copyFlag ) : std::valarray<T>(t,n)	{	setPtr0();	};
	
			//	make a copy
	valarray0_without_ownership ( const T* t, size_t n ) : std::valarray<T>( t, n )		{	setPtr0();	};


//		Make vectors from existing ones
//			make a copy
	valarray0_without_ownership ( const std::valarray<T>& v )  : std::valarray<T>( v )				{	setPtr0();	};
	valarray0_without_ownership ( const valarray0_without_ownership<T>& v ) : std::valarray<T>( v )	{	setPtr0();	};

};


#define	valarray00	valarray0_without_ownership



template <class T> 
class valarray0 : public valarray00<T> {

public:
	enum	cutDirection	{ fixed_x, fixed_y	};

	void  resize_but_keep (size_t sz, T c = T ())
		{	if (valarray00<T>::size() != sz)
			{	size_t	min_sz = std::min(sz,valarray00<T>::size());
				T* keep_data = new T [min_sz];
				for (size_t i=0; i<min_sz; i++)		keep_data[i] = std::valarray<T>::operator[](i);
				valarray00<T>::resize(sz,c);
				for (size_t i=0; i<min_sz; i++)		std::valarray<T>::operator[](i) = keep_data[i];
				for (size_t i=min_sz; i<sz; i++)	std::valarray<T>::operator[](i) = c;
				delete [] keep_data;
			};
		};
 
//		Make new vectors
	valarray0 ()			: valarray00<T>() 		{};
	valarray0 (size_t n)	: valarray00<T>( n )	{};
	valarray0 (const T& t, size_t n) : valarray00<T>( t, n )	{};		//	fill all with t
   	
   		

//		Make vectors from existing pointer
			//	don't make a copy.
//	valarray0 ( T* t, size_t n, abstract_valarray_0<T>::copyFlag f) : valarray00<T>(t, n, f)			{};
	
			//	make a copy
	valarray0 ( const T* t, size_t n ) : valarray00<T>( t, n )	{};


//		Make vectors from existing ones
//		make a copy
	valarray0 ( const std::valarray<T>& v ) : valarray00<T>( v )		{};
    // valarray0 ( const valarray0<T>& v ) 	: valarray00<T>( v )		{};
		
public:
   		//	For x-scales: fusion of 2 increasing scales
		//		the 2 routines are parallel
		//	(for stricly increasing series only !!!)
    size_t	nb_diff_values( const valarray0<T>& v, Boolean commonRange = true, T eps = abstract_valarray_0<T>::valarray0_epsilon() ) const
    	{	T* p1 = valarray00<T>::ptr();			T* p1m = valarray00<T>::ptr()+valarray00<T>::size()-1;
    		T* p2 = v.ptr();						T* p2m = v.ptr()+v.size()-1;
    		T vmin = std::max(*p1,*p2);		//	T vmax = std::min(*p1m,*p2m);
    		size_t nbdiff = 0;
	    	if (commonRange)
    		{	while (p1 <= p1m && *p1 < vmin)	p1++;
    			while (p2 <= p2m && *p2 < vmin)	p2++;
    		}
    		while (p1 <= p1m && p2 <= p2m)
	    	{	if	(std::fabs(*p1-*p2)<= eps)	{	nbdiff++;	p1++;	p2++;	}
	    		else if	(*p1-*p2 < eps)		{	nbdiff++;	p1++;	}
	    		else if	(*p2-*p1 < eps)		{	nbdiff++;	p2++;	}
	    	}
	    	if (!commonRange)
    		{	while (p1 <= p1m)	{	nbdiff++;	p1++;	}
    			while (p2 <= p2m)	{	nbdiff++;	p2++;	}
    		}
    		return nbdiff;
    	}
	void  insert_vector( const valarray0<T>& v, Boolean commonRange = true, T eps = abstract_valarray_0<T>::valarray0_epsilon() )
		{	size_t	sz = nb_diff_values( v, commonRange, eps );
			if (valarray00<T>::size() != sz)
			{	T* newdata = new T [sz];			T* p	= newdata;
				T* p1 = valarray00<T>::ptr();		T* p1m	= valarray00<T>::ptr()+valarray00<T>::size()-1;
	    		T* p2 = v.ptr();					T* p2m	= v.ptr()+v.size()-1;
    			T vmin = std::max(*p1,*p2);		//	T vmax = std::min(*p1m,*p2m);
		    	if (commonRange)
	    		{	while (p1 <= p1m && *p1 < vmin)	p1++;			
	    			while (p2 <= p2m && *p2 < vmin)	p2++;
	    		}
	    		while (p1 <= p1m && p2 <= p2m)
		    	{	if	(std::fabs(*p1-*p2)<= eps)	{	*p++ = *p1++;	p2++;	}
		    		else if	(*p1-*p2 < eps)		*p++ = *p1++;
		    		else if	(*p2-*p1 < eps)		*p++ = *p2++;
		    	}
	    		if (!commonRange)
	    		{	while (p1 <= p1m)	*p++ = *p1++;
	    			while (p2 <= p2m)	*p++ = *p2++;
	    		}
				this->change_ptr( newdata, sz );
			};
		};
		
};


template <class T, int offsetX> 
class val1Darray : public valarray0<T>, public abstract_valarray_1D<T,offsetX>	{
public:

	val1Darray() : valarray0<T>() 										{};
	val1Darray( size_t n ) : valarray0<T>( n )							{};
	val1Darray( const T& t, size_t n ) : valarray0<T>( t, n )			{};
	val1Darray( const T* t, size_t n ) : valarray0<T>( t, n )			{};
	val1Darray( const std::valarray<T>& v ) : valarray0<T>( v )			{};
    //val1Darray( const val1Darray<T,offsetX>& v ) : valarray0<T>( v )	{};
	val1Darray( const std::valarray<T>& v, T (*func)(T) ) : valarray0<T>( v.size() )		{	this->set_to( v, func );	};
	val1Darray( const std::valarray<T>& v, TFunction<T>& func ) : valarray0<T>( v.size() )	{	this->set_to( v, func );	};

	val1Darray( const abstract_valarray_2D_C0<T>& v, typename valarray0<T>::cutDirection dir,
				size_t pos, size_t off = 0, size_t siz = 0, size_t jumpsiz = 1 )
		: valarray0<T>( ( siz>0 ? siz : (dir == valarray0<T>::fixed_x ? v.dim2() : v.dim1()) - off ) )

		{	size_t dim = (dir == valarray0<T>::fixed_x ? v.dim2() : v.dim1());
			size_t end_index = ( siz>0 ? std::min( dim-1, off+(siz-1)*jumpsiz ) : dim-1 );
			if (dir == valarray0<T>::fixed_x)
			{	T*	first = v.firstPtr()+v.pos_of(pos,off);	T*	last = v.firstPtr()+v.pos_of(pos,end_index);	T*	result = firstPtr();
				while (first <= last)
				{	*result++ = *first;		first += jumpsiz;	}
			}
			else
			{	T*	first = v.firstPtr()+v.pos_of(off,pos);	T*	last = v.firstPtr()+v.pos_of(end_index,pos);	T*	result = firstPtr();
				while (first <= last)
				{	*result++ = *first;		first += v.dim2()*jumpsiz;	}
			};
		};
	
    T    operator[] (size_t pos) const	{	return abstract_valarray_1D<T,offsetX>::operator[](pos);	}
    T&   operator[] (size_t pos)		{	return abstract_valarray_1D<T,offsetX>::operator[](pos);	}
    T    operator() (size_t pos) const	{	return abstract_valarray_1D<T,offsetX>::operator()(pos);	}
    T&   operator() (size_t pos)		{	return abstract_valarray_1D<T,offsetX>::operator()(pos);	}
    
	virtual T*			firstPtr() 	const		{	return	valarray0<T>::ptr();				};
	virtual T*			lastPtr()  	const		{	return	this->firstPtr()+this->size()-1;	};
	virtual size_t		length()	const		{	return	std::valarray<T>::size();			};
	inline size_t		size()		const		{	return	std::valarray<T>::size();			};

};


template <class T> 
class valarray0_box					{		//	only the access, not the content
protected:
	const T*		temp_ptr;
	const size_t	temp_siz;
public:
    valarray0_box() 									: temp_ptr(0), temp_siz(0) 	    				{};
    valarray0_box( const T* t, size_t n )				: temp_ptr(t), temp_siz(n) 						{};
	valarray0_box( const std::valarray<T>& v )			: temp_siz(v.size()), temp_ptr(v.firstPtr())	{};
	valarray0_box( const abstract_valarray_0<T>& v )	: temp_siz(v.size()), temp_ptr(v.firstPtr())	{};
};


//namespace MyMath { double	bissectbracket( const a_double_vector1& x, const double& value, unsigned long& klo );	}


template <class T, int offsetX> 
class val1Darray_box : public valarray0_box<T>, public abstract_valarray_1D<T,offsetX>	{		//	only the access, not the content
public:
	static size_t	make_first_index_sup( const val1Darray<T,offsetX>& v, T v1 )		//	assuming an increasing scale, get i_start minimum such as v1 <= v[i_start]
		{	size_t klo = 1;
			v.bissectbracket( v1, klo );
			if (v[klo] < v1)	klo++;
			return klo;
		}
	static size_t	make_last_index_inf( const val1Darray<T,offsetX>& v, T v2 )			//	assuming an increasing scale, get i_end max such as v2 >= v[i_end]
		{	size_t klo = v.length()-1;
			v.bissectbracket( v2, klo );	klo++;
			if (v[klo] > v2)	klo--;
			return klo;
		}
public:
	val1Darray_box()							: valarray0_box<T>()		{};
	val1Darray_box( const T* t, size_t n )		: valarray0_box<T>(t,n)	{};
	val1Darray_box( const std::valarray<T>& v )	: valarray0_box<T>(v)		{};
	
	val1Darray_box( const val1Darray<T,offsetX>& v, size_t i_start, size_t i_end )
			: valarray0_box<T>( v.firstPtr()+i_start-1, i_end - i_start + 1 )
		{	if (i_start < 1 || i_start > v.length() || i_end < 1 || i_end > v.length() || i_end < i_start)		
				throw	GenericError("Pb!!");
		};
		
//	val1Darray_box( std::valarray<T>& v, T v_start, T v_end )					//	assuming an increasing scale
//			: valarray0_box( &(v[make_first_index_sup( v, v_start )]), make_last_index_inf( v, v_end ) - make_first_index_sup( v, v_start ) + 1 )	{};
	
    T    operator[] (size_t pos) const	{	return abstract_valarray_1D<T,offsetX>::operator[](pos);	}
    T&   operator[] (size_t pos)		{	return abstract_valarray_1D<T,offsetX>::operator[](pos);	}
    T    operator() (size_t pos) const	{	return abstract_valarray_1D<T,offsetX>::operator()(pos);	}
    T&   operator() (size_t pos)		{	return abstract_valarray_1D<T,offsetX>::operator()(pos);	}
    
	virtual T*			firstPtr() 	const		{	return	(T*)(valarray0_box<T>::temp_ptr);		};
	virtual T*			lastPtr()  	const		{	return	(T*)(this->firstPtr()+this->size()-1);	};
	virtual size_t		length()	const		{	return	valarray0_box<T>::temp_siz;				};
	inline size_t		size()		const		{	return	valarray0_box<T>::temp_siz;				};

};






template <class T> 
class val1DarrayN : public valarray0<T>, public abstract_valarray_1D<T,0>	{
private:
	size_t	offX;
public:
	val1DarrayN() : valarray0<T>(), offX(0) 														{};
	val1DarrayN( size_t n0, size_t n1 ) : valarray0<T>( n1-n0+1 ), offX(n0)							{};
	val1DarrayN( const T& t, size_t n0, size_t n1 ) : valarray0<T>( t, n1-n0+1 ), offX(n0)			{};
	val1DarrayN( const T* t, size_t n0, size_t n1 ) : valarray0<T>( t, n1-n0+1 ), offX(n0)			{};
//	val1DarrayN( T* t, size_t n0, size_t n1, typename valarray0<T>::copyFlag iscopy ) : valarray0<T>( t, n1-n0+1, iscopy ), offX(n0)		{};
//	val1DarrayN( std::valarray<T>& v, size_t n0, typename valarray0<T>::copyFlag iscopy ) : valarray0<T>( v, iscopy ), offX(n0)				{};
	val1DarrayN( const std::valarray<T>& v, size_t n0 ) : valarray0<T>( v ), offX(n0)						{};
	val1DarrayN( const val1DarrayN<T>& v ) : valarray0<T>( v ), offX(v.offX)								{};
	val1DarrayN( const std::valarray<T>& v, size_t n0, T (*func)(T) ) : valarray0<T>( v, func ), offX(n0)	{};

	void resizeN( size_t n0, size_t n1 )
		{	offX = n0;	valarray0<T>::resize_but_keep(n1-n0+1);	};
		
	
    T    operator[] (size_t pos) const	{	return abstract_valarray_1D<T,0>::operator[](pos-offX);	}
    T&   operator[] (size_t pos)		{	return abstract_valarray_1D<T,0>::operator[](pos-offX);	}
    T    operator() (size_t pos) const	{	return abstract_valarray_1D<T,0>::operator()(pos-offX);	}
    T&   operator() (size_t pos)		{	return abstract_valarray_1D<T,0>::operator()(pos-offX);	}
    
	virtual T*			firstPtr() 	const		{	return	valarray0<T>::ptr();				};
	virtual T*			lastPtr()  	const		{	return	this->firstPtr()+this->size()-1;	};
	virtual size_t		length()	const		{	return	std::valarray<T>::size();			};
	inline size_t		size()		const		{	return	std::valarray<T>::size();			};

};

template <class T, int offsetX> 
class zeroMean1Darray : public val1Darray<T,offsetX> {
public:
	zeroMean1Darray( const val1Darray<T,offsetX>& v ) : val1Darray<T,offsetX>( v.size() )
	{	this->do_copy( v, abstract_valarray_0<T>::remove_mean );	};
};

template <class T, int offsetX, int offsetY> 
class val2Darray1 : public valarray0<T>, public abstract_valarray_2D_C<T,offsetX,offsetY>	{
public:
	
	val2Darray1() : valarray0<T>()	{};
	val2Darray1( size_t nx, size_t ny ) : valarray0<T>( nx*ny ), abstract_valarray_2D_C<T,offsetX,offsetY>(nx,ny)	{};
	val2Darray1( const T& t, size_t nx, size_t ny ) : valarray0<T>( t, nx*ny ), abstract_valarray_2D_C<T,offsetX,offsetY>(nx,ny)	{};
	val2Darray1( const T* t, size_t nx, size_t ny ) : valarray0<T>( t, nx, ny ), abstract_valarray_2D_C<T,offsetX,offsetY>(nx,ny)	{};
	val2Darray1( const val2Darray1<T,offsetX,offsetY>& v ) : valarray0<T>( v ), abstract_valarray_2D_C<T,offsetX,offsetY>( v.dim1(), v.dim2() )	{};
	
	val2Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, T (*func)(T,T) )
		: valarray0<T>( x.size() * y.size() ), abstract_valarray_2D_C<T,offsetX,offsetY>( x.size(), y.size() )	{	set_to_2D( x, y, func );	};
	val2Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, TFunction2<T>& func )
		: valarray0<T>( x.size() * y.size() ), abstract_valarray_2D_C<T,offsetX,offsetY>( x.size(), y.size() )	{	this->set_to_2D( x, y, func );	};
		
	val2Darray1&	operator= (const T& t)
   		{	(std::valarray<T>::operator=)(t);	return *this;	};

  	void  resize2D (size_t d1, size_t d2, T c = T ())
		{	valarray0<T>::resize_but_keep( d2 * d1, c );
			abstract_valarray_2D_C0<T>::sizex = d1;	abstract_valarray_2D_C0<T>::sizey = d2;
		};
	
    T    operator() (size_t posx, size_t posy) const	{	return abstract_valarray_2D_C<T,offsetX,offsetY>::operator()(posx,posy);	}
    T&   operator() (size_t posx, size_t posy)			{	return abstract_valarray_2D_C<T,offsetX,offsetY>::operator()(posx,posy);	}
    
	virtual T*			firstPtr() 	const		{	return	valarray0<T>::ptr();				};
	virtual T*			lastPtr()  	const		{	return	this->firstPtr()+this->size()-1;	};
	virtual size_t		length()	const		{	return	std::valarray<T>::size();			};
	inline size_t		size()		const		{	return	std::valarray<T>::size();			};
		
    val1Darray<T,offsetX>    operator[] (size_t posx) const
    	{	return	val1Darray<T,offsetX>( abstract_valarray_2D_C<T,offsetX,offsetY>::access2Ptr(posx,0), abstract_valarray_2D_C0<T>::dim2(), valarray0<T>::dontCopy );	}
    
    void	set_colum( size_t posx, const val1Darray<T,offsetX>& v )
    	{	size_t n_y = std::min( v.length(), abstract_valarray_2D_C0<T>::dim2() );
			T* p  = &((*this)(posx,offsetY));		T* pv = v.firstPtr();
			for (; n_y > 0; --n_y )					{	*p++ = *pv++;	};
    	}
};


template <class T, int offsetX, int offsetY> 
class val2Darray_box : public valarray0_box<T>, public abstract_valarray_2D_C<T,offsetX,offsetY>	{
public:
	val2Darray_box()								: valarray0_box<T>()			{};
	val2Darray_box( T* t, size_t nx, size_t ny )	: valarray0_box<T>(t,nx*ny), abstract_valarray_2D_C<T,offsetX,offsetY>(nx,ny)			{};
	val2Darray_box( abstract_valarray_2D_C0<T>& v )	: valarray0_box<T>(v), abstract_valarray_2D_C<T,offsetX,offsetY>(v.dim1(),v.dim2())	{};
	
    T    operator() (size_t posx, size_t posy) const	{	return abstract_valarray_2D_C<T,offsetX,offsetY>::operator()(posx,posy);	}
    T&   operator() (size_t posx, size_t posy)			{	return abstract_valarray_2D_C<T,offsetX,offsetY>::operator()(posx,posy);	}
    
	virtual T*			firstPtr() 	const		{	return	(T*)(valarray0_box<T>::temp_ptr);		};
	virtual T*			lastPtr()  	const		{	return	(T*)(this->firstPtr()+this->size()-1);	};
	virtual size_t		length()	const		{	return	valarray0_box<T>::temp_siz;				};
	inline size_t		size()		const		{	return	valarray0_box<T>::temp_siz;				};
};










template <class T, int offsetX, int offsetY, int offsetZ>
class val3Darray1 : public valarray0<T>, public abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>	{
public:

//		Make new vectors
    val3Darray1() : valarray0<T>() , abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>()	{};
	val3Darray1 (size_t nx, size_t ny, size_t nz) : valarray0<T>( nx*ny*nz ), abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(nx,ny,nz)	{};
    		
	val3Darray1 (const T& t, size_t nx, size_t ny, size_t nz)				//	fill all with t
			: valarray0<T>( t, nx*ny*nz ), abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(nx,ny,nz)	{};
  
//		Assignment
	val3Darray1&	operator= (const val3Darray1& v)
   		{	(valarray0<T>::operator=)(v);	abstract_valarray_3D_C0<T>::sizex = v.dim1();	abstract_valarray_3D_C0<T>::sizey = v.dim2();	abstract_valarray_3D_C0<T>::sizez = v.dim3();	return *this;	};
  
//		Make vectors from existing pointer
			//	make a copy
	val3Darray1( const T* t, size_t nx, size_t ny, size_t nz )
    	: valarray0<T>( t, nx*ny*nz ), abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(nx,ny,nz)		{};
		
			//	don't make a copy. On output, t = NULL
//	val3Darray1( T* t, size_t nx, size_t ny, size_t nz, abstract_valarray_0<T>::copyFlag iscopy )
//		: valarray0<T>( t, nx*ny*nz, iscopy ), abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(nx,ny,nz)	{};

			//	make a copy
	val3Darray1 ( const val3Darray1& v )
		: valarray0<T>( v ), abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(v.dim1(),v.dim2(),v.dim3())			{};


	val3Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, T (*func)(T,T,T) )
		: valarray0<T>( x.size() * y.size() * z.size() ),abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(x.size(),y.size(),z.size())
		{	set_to_3D( x, y, z, func );	};
	val3Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, TFunction3<T>& func )
		: valarray0<T>( x.size() * y.size() * z.size() ),abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(x.size(),y.size(),z.size())
		{	this->set_to_3D( x, y, z, func );	};

/*
    //	fastest index = last index (v111, v112, ...)
	explicit val3Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, T (*func)(T,T,T) )
		: valarray0<T>( x.size() * y.size() * z.size() ),abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(x.size(),y.size(),z.size())
		{	const T* first_x = (x.firstPtr());	const T* last_x = first_x + x.size();
			T* result = firstPtr();
			while (first_x != last_x)
		    {	const T* first_y = (y.firstPtr());	const T* last_y = first_y + y.size();
				while (first_y != last_y)
				{	const T* first_z = (z.firstPtr());	const T* last_z = first_z + z.size();
					while (first_z != last_z)
					{	*result =func( *first_x, *first_y, *first_z );	++result, ++first_z;	};
					++first_y;
			    }
				++first_x;
		    };
		};
	explicit val3Darray1 ( const abstract_valarray_0<T>& x, const abstract_valarray_0<T>& y, const abstract_valarray_0<T>& z, TFunction3<T>& func )
		: valarray0<T>( x.size() * y.size() * z.size() ),abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>(x.size(),y.size(),z.size())
		{	const T* first_x = (x.firstPtr());	const T* last_x = first_x + x.size();
			T* result = firstPtr();
			while (first_x != last_x)
		    {	const T* first_y = (y.firstPtr());	const T* last_y = first_y + y.size();
				while (first_y != last_y)
				{	const T* first_z = (z.firstPtr());	const T* last_z = first_z + z.size();
					while (first_z != last_z)
					{	*result =func( *first_x, *first_y, *first_z );	++result, ++first_z;	};
					++first_y;
			    }
				++first_x;
		    };
		};
*/

    T    operator() (size_t posx, size_t posy, size_t posz) const	{	return abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>::operator()(posx,posy,posz);	}
    T&   operator() (size_t posx, size_t posy, size_t posz)			{	return abstract_valarray_3D_C<T,offsetX,offsetY,offsetZ>::operator()(posx,posy,posz);	}
    
	virtual T*			firstPtr() 	const		{	return	valarray0<T>::ptr();		};
	virtual T*			lastPtr()  	const		{	return	firstPtr()+size()-1;		};
	virtual size_t		length()	const		{	return	std::valarray<T>::size();	};
	inline size_t		size()		const		{	return	std::valarray<T>::size();	};
	
};



typedef	abstract_valarray_1D<double,0>		a_double_vector0;
typedef	abstract_valarray_1D<double,1>		a_double_vector1;
typedef	abstract_valarray_2D_C<double,0,0>	a_double_matrix0;
typedef	abstract_valarray_2D_C<double,1,1>	a_double_matrix1;


typedef	val1Darray<float,0>			float_vector0;			//	C type vectors, from 0 to n-1
typedef	val1Darray<float,1>			float_vector1;			//	usual vectors, from 1 to n
typedef	val1DarrayN<float>			float_vectorN;
typedef	val2Darray1<float,0,0>		float_matrix0;			//	matrix, from 0 to nx-1, 0 to ny-1
typedef	val2Darray1<float,1,1>		float_matrix1;			//	usual matrix, from 1 to nx, 1 to ny
typedef	val3Darray1<float,1,1,1>	float_3tensor1;			//	usual 3D-tensor, from 1 to nx, 1 to ny, 1 to nz

typedef	val1Darray<double,0>		double_vector0;
typedef	val1Darray<double,1>		double_vector1;
typedef	val1DarrayN<double>			double_vectorN;
typedef	val2Darray1<double,0,0>		double_matrix0;
typedef	val2Darray1<double,1,1>		double_matrix1;
typedef	val2Darray1<double,1,0>		double_matrix10;
typedef	val2Darray1<double,0,1>		double_matrix01;
typedef	val3Darray1<double,1,1,1>	double_3tensor1;

typedef	val1Darray_box<double,0>	temp_double_vector0;
typedef	val1Darray_box<double,1>	temp_double_vector1;
typedef	val2Darray_box<double,1,1>	temp_double_matrix1;

typedef	val1Darray<int,0>			int_vector0;
typedef	val1Darray<int,1>			int_vector1;
typedef	val1DarrayN<int>			int_vectorN;
typedef	val2Darray1<int,0,0>		int_matrix0;
typedef	val2Darray1<int,1,1>		int_matrix1;
typedef	val3Darray1<int,1,1,1>		int_3tensor1;

typedef	val1Darray<long,0>			long_vector0;
typedef	val1Darray<long,1>			long_vector1;
typedef	val1DarrayN<long>			long_vectorN;
typedef	val2Darray1<long,0,0>		long_matrix0;
typedef	val2Darray1<long,1,1>		long_matrix1;
typedef	val3Darray1<long,1,1,1>		long_3tensor1;

typedef	val1Darray<bool,0>			bool_vector0;
typedef	val1Darray<bool,1>			bool_vector1;
typedef	val1DarrayN<bool>			bool_vectorN;
typedef	val2Darray1<bool,0,0>		bool_matrix0;
typedef	val2Darray1<bool,1,1>		bool_matrix1;
typedef	val3Darray1<bool,1,1,1>		bool_3tensor1;

typedef	val1Darray<void*,0>			ptr_vector0;
typedef	val1Darray<void*,1>			ptr_vector1;
typedef	val1DarrayN<void*>			ptr_vectorN;
typedef	val2Darray1<void*,0,0>		ptr_matrix0;
typedef	val2Darray1<void*,1,1>		ptr_matrix1;
typedef	val3Darray1<void*,1,1,1>	ptr_3tensor1;

typedef	val1Darray<size_t,0>		ulong_vector0;
typedef	val1Darray<size_t,1>		ulong_vector1;
typedef	val1DarrayN<size_t>			ulong_vectorN;
typedef	val2Darray1<size_t,0,0>		ulong_matrix0;
typedef	val2Darray1<size_t,1,1>		ulong_matrix1;
typedef	val3Darray1<size_t,1,1,1>	ulong_3tensor1;

typedef	val1Darray<std::string,0>		string_vector0;
typedef	val1Darray<std::string,1>		string_vector1;
typedef	val1DarrayN<std::string>		string_vectorN;
typedef	val2Darray1<std::string,0,0>	string_matrix0;
typedef	val2Darray1<std::string,1,1>	string_matrix1;
typedef	val3Darray1<std::string,1,1,1>	string_3tensor1;

typedef	val1Darray<double_vector1,1>		double_vector_vector;


//	default vectors

typedef	float_vector1		float_vector;			//	default is usual vectors, from 1 to n
typedef	float_matrix1		float_matrix;			//	default is usual matrix, from 1 to nx, 1 to ny
typedef	float_3tensor1		float_3tensor;			//	default is usual matrix, from 1 to nx, 1 to ny

typedef	double_vector1		double_vector;
typedef	double_matrix1		double_matrix;
typedef	double_3tensor1		double_3tensor;

typedef	zeroMean1Darray<double,1>	zeroMean_double_vector;

typedef	int_vector1			int_vector;
typedef	int_matrix1			int_matrix;
typedef int_3tensor1		int_3tensor;

typedef	long_vector1		long_vector;
typedef	long_matrix1		long_matrix;
typedef long_3tensor1		long_3tensor;

typedef	ulong_vector1		ulong_vector;
typedef	ulong_matrix1		ulong_matrix;
typedef ulong_3tensor1		ulong_3tensor;

typedef	bool_vector1		bool_vector;
typedef	bool_matrix1		bool_matrix;
typedef bool_3tensor1		bool_3tensor;

typedef	ptr_vector1			ptr_vector;
typedef	ptr_matrix1			ptr_matrix;
typedef ptr_3tensor1		ptr_3tensor;

typedef	string_vector1		string_vector;
typedef	string_matrix1		string_matrix;
typedef string_3tensor1		string_3tensor;



template <class T> 
class pointer_with_ownership {
private:
	bool	owner;
	T*&		p;
public:
//	pointer_with_ownership()							{	p = NULL;	owner = true;	};
	pointer_with_ownership( T*& pp ):p(pp)				{	p = pp;	owner = true;		};
//	pointer_with_ownership( T*& pp, bool owns ),p(pp)	{	p = pp;	owner = owns;		};
	~pointer_with_ownership()							{	if (owner) delete p;		};
	void	change_ownership( bool owns )				{	owner = owns;	};
};


typedef	pointer_with_ownership<double_vector>	double_vector_ptr;


//	some functional types

typedef	double	(*doubleFuncV)( const double_vector& );
typedef	void	(*vectorFuncV)( const double_vector&, double_vector& );




//	---
#pragma mark -
#pragma mark //	MyMath

namespace MyMath {

using namespace std;

void		ThrowError( char* x );
void		ThrowErrorIf( bool x );
void		ThrowErrorIfNot( bool x );


const double Pi				= 3.14159265358979323846264338;
const double PiSur2			= Pi/2;
const double PiSur4			= Pi/4;
const double PiSurDeux		= Pi/2;
const double DeuxPi			= Pi*2;
const double PiSur180		= Pi/180;
const double PiSur648000	= PiSur180/3600;
const double SqrtPi			= std::sqrt(Pi);
const double Sqrt2			= std::sqrt(2.0);
const double SqrtPiSur2		= SqrtPi/2;


//		UTILITIES
template <class T>
inline T	sqr( T x )					{	return x*x;	};
template <class T>
inline T	multiplicate( T x, T y )	{	return x*y;	};
template <class T>
inline T	sign_of( T x )				{	return (x > 0 ? 1 : (x < 0 ? -1 : 0));		};
template <class T>
inline T	set_sign( T a, T b ) 		{	return (b >= 0.0 ? fabs(a) : -fabs(a));		};

template <class T>
inline T	MinOf( T x, T y )			{	return (x > y ? y : x);		};
template <class T>
inline T	MaxOf( T x, T y )			{	return (x < y ? y : x);		};

template <class T>
inline T	pythag( T a, T b )			{	T absa = abs(a);	T absb = abs(b);
											return (absa > absb ? absa*sqrt(1.0 + sqr(absb/absa))
												: (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb))) );
										};
template <class T>
inline void	swap( T& a, T& b )			{	T temp = a;	a = b;	b = temp;	};
										
inline long		long_floor( double x )	{	return (long)std::floor(x);		};
inline long		long_ceil( double x )	{	return (long)std::ceil(x);		};
inline long		long_round( double x )	{	return std::roundl(x);			};		//	NOTE: not yet in ANSI C

inline unsigned long	ulong_floor( double x )	{	return (unsigned long)std::floor(x);		};
inline unsigned long	ulong_ceil( double x )	{	return (unsigned long)std::ceil(x);		};

inline double	fractional_part( double x, double y )
					{	double	i;	return	( x < 0 ? 1-std::modf( -x/y, &i ): std::modf( x/y, &i ) );	};
inline double	magnitude10( double x )
					//{	return (x == 0 ? 0 : ((x < 0) ? -1 : 1) * std::pow( 10, (int)std::floor(std::log10(std::fabs(x))) ));	};
					{	return (x == 0 ? 0 : ((x < 0) ? -1 : 1) * std::pow( 10, std::floor(std::log10(std::fabs(x))) ));	};


typedef	enum { round_up, round_down, roundP_up, roundP_down, round_nearest } roundType;
//	round_up:		0 < x <= r   or  x <= r < 0		(always x <= r)
//	round_down:		0 < r <= x   or  r <= x < 0		(always x >= r)
//	roundP_up:		0 < x <= r   or  r <= x < 0		(always |x| <= |r|)
//	roundP_down:	0 < r <= x   or  x <= r < 0		(always |x| >= |r|)

double	Arondi( double x, roundType t, int n );
/*{	if (x == 0)		return 0;
	double	mz = magnitude10(x);
	double	z  = x/mz;								//	1 <= z < 10
	if (n > 0)										//	n significant digits
	{	double	powerten = std::pow( (double)10, n-1 );
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
}*/

//inline size_t	nextPower2( size_t n )	{	return	(n==0 ? 0 : (size_t) exp2(ulong_ceil(log2(n))));	};
//HP debug: replace above to avoid Intel/PPC difference of long double treatment (80 vs 64 bits)
inline size_t nextPower2( size_t n )
	{	if (n<1) return 0;
		size_t p = 1;
		while ( p < n )		p *= 2;
		return p;
	};
	



//		BINARY SEARCH

//	hunt: given a vector x of length n (1 to n) of monotonically increasing or decreasing values,
//			a value v, and an initial guess jlo, returns the index jlo such that v is between x[jlo] and x[jlo+1]
//			with jlo = 0 if v 'before' x[1], and jlo = n if v 'after' x[n]
void	hunt( const a_double_vector1& x, const double& v, unsigned long& jlo );

//	bissectbracket:	same as hunt but returns klo between 1 and nb-1, so that x[klo] and x[klo+1] can
//			be used e.g. for inter(extra)polation. Return current step size (x[klo+1] - x[klo])
double	bissectbracket( const a_double_vector1& x, const double& value, unsigned long& klo );



//		SPLINE FUNCTIONS

/*
const double unused_spline_yp = 2e30;

//	spline_init:	given a monotonically increasing vector x and a vector y, computes the second
//			derivative vector y2. yp1 and ypn are the first derivatives at the end points. If unused
//			the second derivative at the end points is set to zero (natural spline).
void	spline_init( const double_vector1& x, const double_vector1& y, double_vector1& y2, double yp1 = unused_spline_yp, double ypn = unused_spline_yp );

//	spline_interp:	given a monotonically increasing vector xa, a vector ya and the second derivative vector y2a, 
//			given an abscissa t and an initial guess for the index klo such that xa[klo] near t, computes the spline
//			interpolation at t.
double	spline_interp( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A, double t, unsigned long& klo );
//	spline_interp_integ1:	...first integral between a and b
double	spline_interp_integ1( const double_vector1& x, const double_vector1& y, const double_vector1& y2,
			double a, double b, unsigned long& klo );
//	spline_interp_deriv1:	...first (second, third) derivative at t
double	spline_interp_deriv1( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A, double t, unsigned long& klo );
double	spline_interp_deriv2( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A, double X, unsigned long& klo );
double	spline_interp_deriv3( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A, double X, unsigned long& klo );


	

//		PSEUDO-SPLINE FUNCTIONS: set to remain monotonously increasing.

Boolean	inc_spline_check( const double_vector1& x, const double_vector1& y, const double_vector1& yR, const double_vector1& yL, double eps );
void	inc_spline_init( const double_vector1& x, const double_vector1& y, double_vector1& y2L, double_vector1& y2R,
						double minRate, int t );
*/


//		LINEAR INTERPOLATING FUNCTIONS

//	linear_interp:	given a monotonically increasing vector xa, a vector ya, 
//			given an abscissa t and an initial guess for the index klo such that xa[klo] near t, computes the linear
//			interpolation at t.
double	linear_interp( const a_double_vector1& XA, const a_double_vector1& YA, double t, unsigned long& klo );
//	linear_interp_deriv1:	...first derivative at t
double	linear_interp_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, double t, unsigned long& klo );
//	linear_interp_deriv1:	...first integral between a and b
double	linear_interp_integ1( const a_double_vector1& XA, const a_double_vector1& YA, double a, double b, unsigned long& klo );


//		STAIR INTERPOLATING FUNCTIONS

double	stair_interp( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& midP, double X, unsigned long& klo );
//	stair_interp_integ1:	...first integral between a and b
double	stair_interp_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& midP, double a, double b, unsigned long& klo );


//		Library procedures

//	LINEAR PROGRAMMING: Simplex method
Boolean	simplx( double_matrix1& a, unsigned long m, unsigned long n, unsigned long m1, unsigned long m2,
					unsigned long m3, int&	icase, ulong_vector1& izrov, ulong_vector1& iposv );
void	simp1( const double_matrix1& a, unsigned long mm, const ulong_vector1& ll, unsigned long nll, int iabf, unsigned long& kp, double& bmax );
void	simp2( const double_matrix1& a, unsigned long n, const ulong_vector1& l2, unsigned long nl2, unsigned long& ip, unsigned long kp, double& q1 );
void	simp3( double_matrix1& a, unsigned long i1, unsigned long k1, unsigned long ip, unsigned long kp );

//	elliptic integrals
double	elle( double phi, double ak );
double	ellf( double phi, double ak );
double	ellpi( double phi, double en, double ak );

double	rc( double x, double y );
double	rd( double x, double y, double z );
double	rf( double x, double y, double z );
double	rj( double x, double y, double z, double p );


//	some statistical functions
double	gammln( double xx );
double	gammp( double a, double x );
double	gammq( double a, double x );
void	gser( double& gamser, double a, double x, double& gln );
void	gcf( double& gammcf, double a, double x, double& gln );
double	erfcc( double x );
double	erff( double x );
double	erffc( double x );
double	betacf( double a, double b, double x);
double	betai( double a, double b, double x );

double	crank( double_vector1& w );
double	spear( const double_vector1&  data1, const double_vector1& data2, /*float *d, float *zd, float *probd, */ double& probrs );

//	random numbers generators

class	randomGenerator {
protected:	static long				idum;
			static unsigned long	udum;
			void	initialize();
public:		virtual double	random() = 0;				//	returns a uniform deviate between 0.0 (exclusive) and 1.0 (exclusive)
			randomGenerator()	{	initialize();	};
			
			double	lorentzdev();
			double	expdev();			//	returns an exponentially distributed (positive) deviate of unit mean
			double	doubleexpdev();		//	... idem symetrized
			double	gasdev();			//	returns a normally distributed deviate with zero mean and unit variance
};

class	randomGenerator0 : public randomGenerator {
public:		virtual double	random();
};
class	randomGenerator1 : public randomGenerator {
public:		virtual double	random();
};
class	randomGenerator2 : public randomGenerator {
public:		virtual double	random();
};
class	randomGenerator3 : public randomGenerator {
public:		virtual double	random();
};
class	randomGenerator4 : public randomGenerator {
public:		virtual double	random();
			randomGenerator4()	{	if (!test_psdes())	ThrowError((char *)"ran4 is not well implemented on this machine" );	};
private:	void	psdes( unsigned long& lword, unsigned long& irword );
			bool	test_psdes();
};
class	randomGeneratorQD2 : public randomGenerator {
public:		virtual double	random();
            randomGeneratorQD2()	{
              if (!test_qd2())	ThrowError((char *)"ranQD2 is not well implemented on this machine" );
              initialize();
            };
private:	bool	test_qd2();
};

class	gaussianWhiteNoise : public val1Darray<double,1>	{
public:		gaussianWhiteNoise( size_t n, randomGenerator& r ) : val1Darray<double,1>(n)					//	with zero mean, unit variance
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = r.gasdev();
 				}
 			gaussianWhiteNoise( size_t n, double s, randomGenerator& r ) : val1Darray<double,1>(n)			//	with zero mean
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = s * r.gasdev();
 				}
 			gaussianWhiteNoise( size_t n, double m, double s, randomGenerator& r ) : val1Darray<double,1>(n)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = m + s * r.gasdev();
 				}
};
class	uniformNoise : public val1Darray<double,1>	{
public:		uniformNoise( size_t n, randomGenerator& r ) : val1Darray<double,1>(n)							//	uniform deviate between 0.0 (exclusive) and 1.0 (exclusive)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = r.random();
 				}
 			uniformNoise( size_t n, double a, double b, randomGenerator& r ) : val1Darray<double,1>(n)		//	uniform deviate between a (exclusive) and b (exclusive)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = a + (b-a) * r.random();
 				}
};
class	exponentialNoise : public val1Darray<double,1>	{
public:		exponentialNoise( size_t n, randomGenerator& r ) : val1Darray<double,1>(n)							//	uniform deviate between 0.0 (exclusive) and 1.0 (exclusive)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = r.expdev();
 				}
};
class	double_exponentialNoise : public val1Darray<double,1>	{
public:		double_exponentialNoise( size_t n, randomGenerator& r ) : val1Darray<double,1>(n)					//	uniform deviate between 0.0 (exclusive) and 1.0 (exclusive)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = r.doubleexpdev();
 				}
};
class	lorentzNoise : public val1Darray<double,1>	{
public:		lorentzNoise( size_t n, randomGenerator& r ) : val1Darray<double,1>(n)							//	uniform deviate between 0.0 (exclusive) and 1.0 (exclusive)
				{	double* p = firstPtr();		const double* last = p + length();
 					while (p != last)			*p++ = r.lorentzdev();
 				}
};


//	quicksort
void	sort( double_vector1& arr );
void	sort2( double_vector1& arr, double_vector1& brr );

//	linear algebra
void	gaussj( double_matrix1& a, double_matrix1& b, size_t nn = 0 );		//	nn = 0 => dimension = a.dim(). Sinon dimension = nn
void	lubksb( const double_matrix1& a, const int_vector1& indx, double_vector1& b );
void	ludcmp( double_matrix1& a, int_vector1& indx, Boolean& even );
void	tridag( const double_vector1& a, const double_vector1& b, const double_vector1& c, const double_vector1& r, double_vector1& u );
void	svbksb( const a_double_matrix1& u, const a_double_vector1& w, const a_double_matrix1& v, const a_double_vector1& b, a_double_vector1& x );
void	svdcmp( a_double_matrix1& a, a_double_vector1& w, a_double_matrix1& v );
void	balanc( double_matrix& a );
void	elmhes( double_matrix& a );
void	hqr( double_matrix& a, double_vector& wr, double_vector& wi );
void	eigvalsrt( double_vector& d, double_vector& v );
void	tred2( double_matrix1& a, double_vector1& d, double_vector1& e, Boolean eigV = true );
void	tqli( double_vector1& d, double_vector1& e, double_matrix1* z );								//	the routine
inline void	tqli( double_vector1& d, double_vector1& e )					{	tqli( d, e, NULL );	};		//	... without eigenvectors
inline void	tqli( double_vector1& d, double_vector1& e, double_matrix1& z )	{	tqli( d, e, &z );	};		//	... with eigenvectors
void	eigsrt( double_vector1& d, double_matrix1& v );

//	minimization


void	avevar( const double_vector1& data, double& ave, double& var );
void	mnbrak( double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, doubleFunction& f );
double	brent( double ax, double bx, double cx, doubleFunction& f, double tol, double& xmin );
double	zbrent( doubleFunction& f, double x1, double x2, double tol );

inline void	mnbrak( double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, doubleFunc f )
			{	return	mnbrak( ax, bx, cx, fa, fb, fc, (doubleFunction&)((const doubleFunction&)(doubleStdFunction(f))) );	};
inline double	brent( double ax, double bx, double cx, doubleFunc f, double tol, double& xmin )
			{	return	brent( ax, bx, cx, (doubleFunction&)((const doubleFunction&)(doubleStdFunction(f))), tol, xmin );	};
inline double	zbrent( doubleFunc f, double x1, double x2, double tol )
			{	return	zbrent( (doubleFunction&)((const doubleFunction&)(doubleStdFunction(f))), x1, x2, tol );	};


//	ODE integration

typedef	void (*derivFuncType)( const double&, const double_vector&, double_vector& );
typedef	void (*odeStepFuncType)( double_vector&, const double_vector&, double&, const double&, double, const double_vector1&, double&, double&, derivFuncType );


void	rkqs  ( double_vector& y, const double_vector& dydx, double& x, const double& htry, double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs );
void	bsstep( double_vector& y, const double_vector& dydx, double& x, const double& htry, double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs );
void	stiff ( double_vector& y, const double_vector& dydx, double& x, const double& htry, double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs );
void	stifbs( double_vector& y, const double_vector& dydx, double& x, const double& htry, double eps, const double_vector1& yscal, double& hdid, double& hnext, derivFuncType derivs );


class	ode_savor	{
private:
	double	asb_dx_sav;
protected:
	double	dx_sav;
	double	last_sav;
public:
	ode_savor( double dx ): dx_sav(dx)
		{	asb_dx_sav = fabs(dx_sav)*.999999;	};				//	numerical saving step > asb_dx_sav
	void	Start( double t )	{	last_sav = t;	};			//	first saving
	void	SaveState( double t, const double_vector& y )
		{	DoSaveState( t, y );
			last_sav = t;
		};
		
	virtual void	DoSaveState( double, const double_vector& ) = 0;
	
	virtual bool	CheckSave( double x )						//	save if time interval exceeds dx_sav
		{	return	(fabs(x-last_sav) >= asb_dx_sav);	};
	virtual double	NextSavingStep( double t, double dt )		//	next integ step should be t + dt
		{	if (CheckSave( t+dt ))
				return	(last_sav + dx_sav) - t;
			return dt;
		};
};

class	std_ode_savor : public ode_savor	{
private:
	int		kmax, kount, nvar;
	bool	vs_time;
public:
	double_vector			xp;
	double_vector_vector	yp;
public:
	std_ode_savor( double dx, double delta_t, size_t nv, bool vs_t = true ): ode_savor(dx), nvar(nv), vs_time(vs_t)
		{	if (dx_sav != 0)								//	with intermediate result saving
			{	kmax = 1 + std::floor( delta_t/dx_sav );
				kount = 0;
				xp = double_vector( kmax );
				if (vs_time)
				{	yp = double_vector_vector( nvar );
                    for (int i=1; i<=nvar; i++)		yp[i] = double_vector( kmax );
				}
				else	yp = double_vector_vector( kmax );
			}	
		};
	virtual void	DoSaveState( double t, const double_vector& y )
		{	if (kmax > 0 && kount < kmax)
			{	xp[++kount] = t;
				if (vs_time)
                {	for (int i=1; i<=nvar; i++)	yp[i][kount] = y[i];
				}
				else
				{	yp[kount] = double_vector( nvar );
                    for (int i=1; i<=nvar; i++)	yp[kount][i] = y[i];
				}
			}
		};
};

class	ode_integrator {
public:
	typedef	bool	(*endIntegType)( double, double, double, const double_vector& );
	typedef	double	(*solveSwitchType)( double, const double_vector& );
	
private:
	int		nok, nbad;				//	number of good steps and bad steps
	double	eps, h1, hmin;			//	small number for odeStepFuncType; first step size; minimal step size

	odeStepFuncType	odeStepF;
	endIntegType	endInteg;		//	function qui signale l'arrêt de l'intégration
	solveSwitchType	solveSwitch;
	derivFuncType	derivs;

	class	solveFunction : public doubleFunction {
	public:
		ode_integrator&		ode;
		double time;
		double step;
		const double_vector& start_v;
		const double_vector& start_dv;
		const double_vector& yscal;
		solveFunction( ode_integrator& o, double t, double h, const double_vector& v, const double_vector& dv, const double_vector& ys ) 
			: ode(o), time(t), step(h), start_v(v), start_dv(dv), yscal(ys)	{};
	protected:
		virtual double	func( const double x )	const
		{	double			hnext, hdid;
			double			h = x;
			double			t = time;
			double_vector	y(start_v);
			(*ode.odeStepF)( y, start_dv, t, h, ode.eps, yscal, hdid, hnext, ode.derivs );		//	does change t, y and (hdid, hnext)
			return	(*ode.solveSwitch)( t, y );
		};
	};
	
protected:
	ode_savor*		savor;
	double	odeint( double_vector& ystart, double t1, double t2 );
	
public:
	static	bool defaultEndInteg( double t, double t1, double t2, const double_vector& )
		{	return ((t-t2)*(t2-t1) >= 0.0);		};
//	static	double defaultEndStep( double t, double h, double t2, const double_vector& )
//		{	return	t2-t;					};

public:
	ode_integrator( derivFuncType derivF, ode_savor* s = NULL )
		:derivs(derivF), savor(s)
		{	setSteps( 10e-6, 0.1, 1e-9 );
			setMethod( MyMath::rkqs );		//	default method
			setEndFunc( defaultEndInteg );
		};
	void	setSteps( double eps_x, double h1_x, double hmin_x )	//	small number for odeStepFuncType; first step size; minimal step size
		{	eps = eps_x;	h1 = h1_x;	hmin = hmin_x;	};
	void	setMethod( odeStepFuncType rkqF )
		{	odeStepF = rkqF;		};
	void	setEndFunc( endIntegType endIntegF = NULL, solveSwitchType	solveSwitchF = NULL )
		{	endInteg = endIntegF;	solveSwitch = solveSwitchF;	};
	virtual void	doIntegrate( double_vector& ystart, double t1, double t2 )
		{	if (savor)	savor->Start( t1 );
			odeint( ystart, t1, t2 );
		};
		
};


typedef	void (*doSwitchFuncType)( const double&, double_vector& );

class	ode_switch_integrator: public ode_integrator {				//	idem ode_integrator, but with discontinuities
private:															//	Here, endInteg means a state switch, not the actual end...
		doSwitchFuncType	doSwitch;
public:
	ode_switch_integrator( derivFuncType derivF, doSwitchFuncType doSwitchF, ode_savor* s = NULL )
		:ode_integrator( derivF, s ), doSwitch(doSwitchF)	{};
	virtual void	doIntegrate( double_vector& ystart, double t1, double t2 )
		{	double	t = t1;
			if (savor)	savor->Start( t );
			bool	firstTime = true;
			while (!defaultEndInteg( t, t1, t2, ystart ))
			{	if (!firstTime)	(*doSwitch)( t, ystart );
				else			firstTime = false;
				t = odeint( ystart, t, t2 );
			};
		};
};


/*		standard setting, with saving trajectory:
	std_ode_savor	savO( dt_sav, t2-t1, nvar );
	ode_integrator	ode( derivF, &savO );
	ode.doIntegrate( ystart, t1, t2 );
	
		...without saving trajectory:
	ode_integrator	ode( derivF );
	ode.doIntegrate( ystart, t1, t2 );
*/


class	vectorFuncType {
public:		virtual void	vecF( const double_vector& x, double_vector& y ) const = 0;
};
class	vStdVectorFuncType : public vectorFuncType	{
private:	vectorFuncV f0;
public:		virtual void	vecF( const double_vector& x, double_vector& y )	const {	(*f0)(x,y);	};
			vStdVectorFuncType( vectorFuncV f ):f0(f)	{};
};

class	bs_extrapol	{
public:
	double_vector1	x;
	double_matrix1	d;
	bs_extrapol( size_t n, size_t kmax ):x(kmax),d(n,kmax) {};
	void	pzextr( size_t iest, double xest, const double_vector1& yest, double_vector1& yz, double_vector1& dy );
	void	rzextr( size_t iest, double xest, const double_vector1& yest, double_vector1& yz, double_vector1& dy );
};
class	vFunction_fdjac : public vectorFuncType	{
public:		//virtual void	vecF( const double_vector& x, double_vector& y ) = 0;
			void	fdjac( const double_vector& x, const double_vector& fvec, double_matrix& df ) const;
};
class	vStdFunction_fdjac : public vFunction_fdjac	{
private:	//vectorFuncV 	f0;
			const vectorFuncType&	f0;
//public:		virtual void	vecF( const double_vector& x, double_vector& y )	{	return (*f0)(x,y);	};
//			vStdFunction_fdjac( vectorFuncV f ):f0(f)	{};
public:		virtual void	vecF( const double_vector& x, double_vector& y )	const {	f0.vecF(x,y);	};
			vStdFunction_fdjac( const vectorFuncType& f ):f0(f)	{};
};
class	vOdeFunction_fdjac : public vFunction_fdjac	{
private:	derivFuncType	f0;
			double		 	time;
public:		virtual void	vecF( const double_vector& x, double_vector& y )		const {	(*f0)(time,x,y);	};
			vOdeFunction_fdjac( derivFuncType df, double t ):f0(df),time(t)	{};
};



inline void	fdjac( const double_vector& x, const double_vector& fvec, double_matrix& df, vectorFuncV vf )
			{	vStdVectorFuncType	svf(vf);	vStdFunction_fdjac z(svf);	z.fdjac( x, fvec, df );	};
inline void	fdjac( const double_vector& x, const double_vector& fvec, double_matrix& df, const vectorFuncType& vf )
			{	vStdFunction_fdjac z(vf);	z.fdjac( x, fvec, df );	};
void	jacobn( const double& x, const double_vector& y, double_vector& dfdx, double_matrix& dfdy, derivFuncType derivs );

void	newt( double_vector& x, Boolean& check, const vectorFuncType& vecfunc );
void	broydn( double_vector& x, Boolean& check, const vectorFuncType& vecfunc );
void	mnewt( double_vector& x, Boolean& notFound, const vFunction_fdjac& funcJ );

inline void	newt( double_vector& x, Boolean& notFound, vectorFuncV vecfunc )
			{	newt( x, notFound, (const vectorFuncType&)(vStdVectorFuncType(vecfunc)) );	};
inline void	broydn( double_vector& x, Boolean& notFound, vectorFuncV vecfunc )
			{	broydn( x, notFound, (const vectorFuncType&)(vStdVectorFuncType(vecfunc)) );	};
inline void	mnewt( double_vector& x, Boolean& notFound, vectorFuncType& vecfunc )
			{	mnewt( x, notFound, (const vFunction_fdjac&)(vStdFunction_fdjac(vecfunc)) );	};
inline void	mnewt( double_vector& x, Boolean& notFound, vectorFuncV vecfunc )
			{	vStdVectorFuncType	svf(vecfunc);	mnewt( x, notFound, (const vFunction_fdjac&)(vStdFunction_fdjac(svf)) );	};

class	vectorFunction {
public:		virtual double	operator()( const double_vector& x ) = 0;
    virtual ~vectorFunction()	{};
};
class	vectorStdFunction : public vectorFunction	{
private:	doubleFuncV 	f0;
public:		vectorStdFunction( doubleFuncV f ):f0(f)	{};
			virtual double	operator()( const double_vector& x )	{	return f0(x);	};
};

class	vFunction_fmin : public vectorFunction	{
private:	const vectorFuncType& f0;
			double	fmin( const double_vector& x )
				{	double	sum = 0;	f0.vecF( x, fx );
					for (size_t i=1; i<=n; i++)	sum += sqr( fx[i] );
					return 0.5*sum;
				}
public:
	double_vector	fx;
	size_t			n;
	vFunction_fmin( const vectorFuncType& f, size_t nn ):f0(f),fx(nn),n(nn)	{};
	virtual double	operator()( const double_vector& x )	{	return fmin(x);	};
};
//class	vStdFunction_fmin : public vFunction_fmin	{
//private:	vStdVectorFuncType sf0;
//public:		vStdFunction_fmin( vectorFuncV& f, size_t nn ):vFunction_fmin(sf0,nn),sf0(f)	{};
//};




void	rkck( const double_vector& y, const double_vector& dydx, const double& x, const double& h, double_vector& yout, double_vector& yerr, derivFuncType derivs );
void	mmid( const double_vector1& y, const double_vector1& dydx, double xs, double htot, size_t nstep, double_vector1& yout, derivFuncType derivs );
void	simpr( const double_vector& y,  const double_vector& dydx, const double_vector& dfdx, const double_matrix& dfdy, const double& xs, const double& htot, size_t nstep, double_vector& yout, derivFuncType derivs );

void	lnsrch( const double_vector& xold, const double& fold, const double_vector& g, double_vector& p, double_vector& x, double& f, double stpmax, Boolean& check, vectorFunction& func );
inline void	lnsrch( const double_vector& xold, const double& fold, const double_vector& g, double_vector& p, double_vector& x, double& f, double stpmax, Boolean& check, doubleFuncV func )
	{	return	lnsrch( xold, fold, g, p, x, f, stpmax, check, (vectorFunction&)((const vectorFunction&)(vectorStdFunction(func))) );	};

void	gauleg( double x1, double x2, double_vector& x, double_vector& w, size_t n );
void	rotate( double_matrix& r, double_matrix& qt, size_t i, double a, double b );
void	rsolv( const double_matrix& a, const double_vector& d, double_vector& b );
void	qrdcmp( double_matrix& a, double_vector& c, double_vector& d, Boolean& singular );
void	qrupdt( double_matrix& r, double_matrix& qt, double_vector& u, const double_vector& v );



//	function codes

const FourCharCode	DerivativeFunctionCode 		= 'Deri';
const FourCharCode	SpectralWindowFunctionCode 	= 'Spec';
const FourCharCode	IdentityFunctionCode 		= 'Iden';
const FourCharCode	PolynomialFunctionCode 		= 'Poly';
const FourCharCode	GaussianSumFunctionCode 	= 'GauS';
const FourCharCode	GaussianFilterFunctionCode 	= 'GauF';
const FourCharCode	LinearInterpFunctionCode 	= 'linI';
const FourCharCode	IncreasingLinearInterpFunctionCode 	= 'il_I';
const FourCharCode	LinearInterpCyclicFunctionCode 		= 'lcyI';
const FourCharCode	SplineInterpFunctionCode 			= 'splI';
const FourCharCode	SplineExtraLinInterpFunctionCode 	= 'exsI';		//	spline interpolation, linear extrapolation
const FourCharCode	IncreasingSplineInterpFunctionCode 	= 'i_sI';
const FourCharCode	StairStartInterpFunctionCode 		= 'sstI';
const FourCharCode	StairEndInterpFunctionCode 			= 'estI';
const FourCharCode	StairInterpFunctionCode 			= 'staI';



class simplefunc : public doubleFunction {
public:
	static	char*	TypeName( FourCharCode c )
		{	switch (c)
			{	case StairInterpFunctionCode:			return (char *)"Stair-mid";			break;
				case StairStartInterpFunctionCode:		return (char *)"Stair-start";		break;
				case StairEndInterpFunctionCode:		return (char *)"Stair-end";			break;
				case LinearInterpFunctionCode:			return (char *)"PiecewiseLinear";	break;
				case SplineInterpFunctionCode:			return (char *)"CubicSpline";		break;
				case PolynomialFunctionCode:			return (char *)"Polynomial";		break;
			}
			return	(char *)"";
		}
protected:
	virtual double	func( const double x ) const {	return	ValueAt( x );	};
public:

	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) = 0;

	virtual double ValueAt( double t ) const = 0;
	virtual double DerivValueAt( double t ) const = 0;
	virtual double IntegrateBetween( double t1, double t2 ) const = 0;
	
	virtual FourCharCode code()	const = 0;
	
	virtual double MeanBetween( double t1, double t2 ) const
				{	return	( t1 == t2 ?  ValueAt( t1 ) : IntegrateBetween( t1, t2 )/(t2 - t1) );	};
};

class simpleDerivfunc : public simplefunc {
private:
	simplefunc& f;
public:
	simpleDerivfunc( simplefunc& f0 ) : f(f0)	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&)
		{	return new simpleDerivfunc(f);	};

	virtual double ValueAt( double t ) 						const {	return	f.DerivValueAt( t );	};
	virtual double DerivValueAt( double )					const {	ThrowError((char *)"not implemented");	return	0;	};
	virtual double IntegrateBetween( double t1, double t2 ) const {	return	f.ValueAt( t2 ) - f.ValueAt( t1 );	};
	
	virtual FourCharCode code()	const	{	return	DerivativeFunctionCode;		};
};


class	fit_func;

const double svdFitTol = 1e-5;
void	svdfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, a_double_matrix1& u, a_double_matrix1& v, a_double_vector1& w,
			double& chisq, fit_func* ff );
void	lfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, const int_vector1& ia,
			double_matrix1& covar, double& chisq, fit_func* ff );
void	absfit( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& sig,
			a_double_vector1& a, double& absv, fit_func* ff );



class spectralwindowfunc : public simplefunc {		//	symetric functions in [-1,1], needs definition only between 0 and 1
public:
	virtual double ValueAt( double t )						const 	{	return (t>1 ? 0 : (t<-1 ? 0 : (t<0 ? DoValueAt(-t) : DoValueAt(t))));				};
	virtual double DerivValueAt( double t )					const 	{	return (t>1 ? 0 : (t<-1 ? 0 : (t<0 ? -DoDerivValueAt(-t) : DoDerivValueAt(t))));	};
	virtual double IntegrateBetween( double t1, double t2 ) const 
		{	if (t2<t1)				return	-IntegrateBetween( t2, t1 );
			if (t2<=-1 || t1>=1 )	return	0;
			double	tt1 = max(-1.0,t1);		double	tt2 = min(1.0,t2);			//	-1 <= tt1 <= tt2 <= 1
			double	p1 = ( tt1 < 0 ? -DoPrimitiv( -tt1 ) : DoPrimitiv( tt1 ) );
			double	p2 = ( tt2 < 0 ? -DoPrimitiv( -tt2 ) : DoPrimitiv( tt2 ) );
			return	p2 - p1;
		};
		
	virtual FourCharCode code()	const	{	return	SpectralWindowFunctionCode;		};
	
	virtual double DoValueAt( double t ) 		const = 0;			//	0 <= t <= 1
	virtual double DoDerivValueAt( double t ) 	const = 0;			//	0 <= t <= 1
	virtual double DoPrimitiv( double t ) 		const = 0;			//	0 <= t <= 1
		
	virtual	double 	FullIntegral()		const {	return IntegrateBetween( -1, 1 );	};
	virtual	double	Energy()			const = 0;				//	fullIntegral of the square !!
	virtual	char*	Name()				= 0;
	
	typedef enum { centered, Right, Left }	positionFlag;
	typedef enum { zero, half, one } 		boundFlag;
		//	centered case:
		//		zero:	v(1) = ValueAt(-1); 	v(n) = ValueAt(1)
		//		one:	v(0) = ValueAt(-1); 	v(n+1) = ValueAt(1) 	(virtuellement, car v[1; n] )
		//		half:	v(1/2) = ValueAt(-1); 	v(n+1/2) = ValueAt(1) 	(virtuellement, car v[1; n] )
		//	Right case:
		//		v(1) = ValueAt(0);
		//	Left case:
		//		v(n) = ValueAt(0);
private:
	positionFlag	pFlag;
	boundFlag		bFlag;
	size_t			n_size;
	double			deno;
	double			offset;
public:
	spectralwindowfunc()	{	SetFlag( centered, half, 10 );	};
	void	SetFlag( positionFlag p, boundFlag b, size_t n )
				{	pFlag = p;	bFlag = b;	n_size = n;
					double	off = (bFlag == zero ? 0 : (bFlag == half ? 1 : 2 ) ) / 2.0;
					deno = (pFlag == centered ? ((double)(n_size-1+2*off))/2.0 : (n_size-1+off) );
					offset = (pFlag == Right ? (-1) : (off - 1 - deno) );
				};
	void	SetFlags( positionFlag p, boundFlag b )
				{	pFlag = p;	bFlag = b;
					double	off = (bFlag == zero ? 0 : (bFlag == half ? 1 : 2 ) ) / 2.0;
					deno = (pFlag == centered ? ((double)(n_size-1+2*off))/2.0 : (n_size-1+off) );
					offset = (pFlag == Right ? (-1) : (off - 1 - deno) );
				};
	void	SetSize( size_t n )
				{	n_size = n;
					double	off = (bFlag == zero ? 0 : (bFlag == half ? 1 : 2 ) ) / 2.0;
					deno = (pFlag == centered ? ((double)(n-1+2*off))/2.0 : (n-1+off) );
					offset = (pFlag == Right ? (-1) : (off - 1 - deno) );
				};
	double	WindowValueAt( size_t i )	{	return	ValueAt( (i+offset)/deno );		}
/*
	void	SetVectorWindow( double_vector& v, boundFlag f = half )
		{	double	off = (f == zero ? 0 : (f == half ? 1 : 2 ) ) / 2.0;
			size_t	n = v.length();
			for (size_t	i = 1; i<=n; i++)	v(i) = ValueAt( 2 * (i-1+off)/((double)(n-1+2*off)) - 1 );
		};
	void	SetHalfVectorWindow( double_vector& v, boundFlag f = half )		//	v(1) = ValueAt(0);
		{	double	off = (f == zero ? 0 : (f == half ? 1 : 2 ) ) / 2.0;
			size_t	n = v.length();
			for (size_t	i = 1; i<=n; i++)	v(i) = ValueAt( (i-1)/((double)(n-1+off)) );
		};
*/
	
private:
	class	solveFunction : public doubleFunction {
	public:	
		double	Nu, Alpha;
		solveFunction( double n, double a ):Nu(n),Alpha(a) {};
	};
	class	solveFunction0 : public solveFunction {
	public:		solveFunction0( double n, double a ):solveFunction(n,a) {};
	protected:	virtual double	func( const double x )		const {	return	gammp( Nu/2, x) - 0.5;			};
	};
	class	solveFunction1 : public solveFunction {
	public:		solveFunction1( double n, double a ):solveFunction(n,a) {};
	protected:	virtual double	func( const double x )		const {	return	gammp( Nu/2, x) - Alpha/2;		};
	};
	class	solveFunction2 : public solveFunction {
	public:		solveFunction2( double n, double a ):solveFunction(n,a) {};
	protected:	virtual double	func( const double x )		const {	return	gammp( Nu/2, x) - 1 + Alpha/2;	};
	};
	class	solveFunction3 : public solveFunction {
	public:		solveFunction3( double n, double a ):solveFunction(n,a) {};
	protected:	virtual double	func( const double x )		const {	return	MyMath::erff(x) - 1 + Alpha;	};
	};
	class	solveFunction4 : public solveFunction {
	public:		solveFunction4( double n, double a ):solveFunction(n,a) {};
	protected:	virtual double	func( const double x )		const {	return	betai( 1, (Nu-2)/2, 2/(2+(Nu-2)*x) ) - 1 + Alpha;	};
	};
	
	double	x1, x2, dcoh, arg;

public:

	double	LowError()			{	return x2;		};
	double	HighError()			{	return x1;		};
	double	CoherenceError()	{	return dcoh;	};
	double	NonZeroCoherence()	{	return tanh(dcoh);	};
	double	PhaseErrorArg()		{	return arg;		};
	
	double	BandWidth( size_t m )		{	return	1/Energy()/(double)m;	};
	
//	void	BTConfid( int m, int n, int window, double level,
//				double* bw, double* x1, double* x2, double* dcoh, double* arg )
	void	SetBTStats( double level, double m_over_n )
	{	//double	isurt;
		const double	tol = 1e-6;
		//double	wi = Energy();
		
		double	Nu  = 2/Energy() / m_over_n;
		double	Alpha = 1 - level;
	
		solveFunction0	f0(Nu,Alpha);
		solveFunction1	f1(Nu,Alpha);
		solveFunction2	f2(Nu,Alpha);
		solveFunction3	f3(Nu,Alpha);
		solveFunction4	f4(Nu,Alpha);
		
		//isurt = Energy() * m_over_n;
		
		
		double	z0 = zbrent( f0, 0, 10000, tol );
		double	z1 = zbrent( f1, 0, z0, tol );
		double	z2 = zbrent( f2, z0, 10000, tol );
		double	z3 = zbrent( f3, 0, 10000, tol ) * sqrt(2.0);
		double	z4 = 0;
		if (Nu > 2)	z4 = zbrent( f4, 0, 10000, tol );
		
		dcoh = z3 * sqrt( Energy() * m_over_n / 2 );
		x1 	 = Nu / (z1 * 2);
		x2   = Nu / (z2 * 2);
		arg  = sqrt( 2/(Nu-2) * z4 );
		
		//	Log( x1 ) < Log( Spectrum ) < Log( x2 )
		//				ie		x1 Spectrum < Spectrum < x2 Spectrum
		//	Arctanh( Coherence ) ± dcoh
		//	tan( Phase ) ± dcoh * sqrt(1/SQR(Coherence) - 1) / SQR( cos(Phase) )
		//	Phase ± arcsin( arg * sqrt(1/SQR(Coherence) - 1) )
	};
};


template <class T, int off, int plus_offset>
class abstract_power2plus_array : public abstract_valarray_1D<T,off> {
public:
	size_t	power2length()		{	return abstract_valarray_0<T>::length() - plus_offset;	};
	void	DoCosFFT()			{	cosft1( *this );	};		//	forward FFT
	void	DoBackwardCosFFT()
		{	cosft1( *this );
			double z = (double)2/(double)power2length();
			*this *= z;
		};
};
template <class T, int off, int plus_offset>
class power2plus_array : public val1Darray<T,off>, public abstract_power2plus_array<T,off,plus_offset>	{
public:
	power2plus_array( size_t n )
		: val1Darray<T,off>( nextPower2(n-plus_offset) + plus_offset )	{};
	power2plus_array( const abstract_valarray_0<T>& v, typename abstract_valarray_0<T>::TrendFlag f = abstract_valarray_0<T>::no )
		: val1Darray<T,off>( nextPower2( v.size() - plus_offset ) + plus_offset )		{	do_copy( v, f );	};
	power2plus_array( const abstract_valarray_0<T>& v, spectralwindowfunc& w, typename abstract_valarray_0<T>::TrendFlag f = abstract_valarray_0<T>::no )
		: val1Darray<T,off>( nextPower2( v.size() - plus_offset ) + plus_offset )
		{	do_copy( v, f );	w.SetSize( v.size() );
	    	for (size_t i=1; i<=v.size(); i++)	(*this)(i) *= w.WindowValueAt(i);
    	};
    	
	virtual T*			firstPtr() 	const		{	return	val1Darray<T,off>::firstPtr();	};
	virtual T*			lastPtr()  	const		{	return	val1Darray<T,off>::lastPtr();	};
	virtual size_t		length()	const		{	return	val1Darray<T,off>::length();	};
	inline size_t		size()		const		{	return	val1Darray<T,off>::size();		};
    T    operator[] (size_t pos)	 const	{	return val1Darray<T,off>::operator[](pos);	}
    T&   operator[] (size_t pos)			{	return val1Darray<T,off>::operator[](pos);	}
    T    operator() (size_t pos)	 const	{	return val1Darray<T,off>::operator()(pos);	}
    T&   operator() (size_t pos)			{	return val1Darray<T,off>::operator()(pos);	}
};
template <class T, int off, int plus_offset>
class power2plus_array_box : public val1Darray_box<T,off>, public abstract_power2plus_array<T,off,plus_offset>	{
public:
	power2plus_array_box( abstract_valarray_0<T>& v )
 		: val1Darray_box<T,off>( v.firstPtr(), ( ( (nextPower2( v.size() - plus_offset ) + plus_offset) > v.size() ) ? (nextPower2( v.size() - plus_offset )/2 + plus_offset) : (nextPower2( v.size() - plus_offset ) + plus_offset) ) )
 		{};
 		
	virtual T*			firstPtr() 	const		{	return	val1Darray_box<T,off>::firstPtr();	};
	virtual T*			lastPtr()  	const		{	return	val1Darray_box<T,off>::lastPtr();	};
	virtual size_t		length()	const		{	return	val1Darray_box<T,off>::length();	};
	inline size_t		size()		const		{	return	val1Darray_box<T,off>::size();		};
    T    operator[] (size_t pos) 	const	{	return val1Darray<T,off>::operator[](pos);	}
    T&   operator[] (size_t pos)			{	return val1Darray<T,off>::operator[](pos);	}
    T    operator() (size_t pos) 	const	{	return val1Darray<T,off>::operator()(pos);	}
    T&   operator() (size_t pos)			{	return val1Darray<T,off>::operator()(pos);	}
};

typedef	abstract_power2plus_array<double,1,1>		a_double_power2plus1_array;		//	length is a power of 2 + 1
typedef	power2plus_array<double,1,1>				double_power2plus1_array;
typedef	power2plus_array_box<double,1,1>			temp_double_power2plus1_array;



template <class T, int off>
class power2_array : public val1Darray<T,off>	{			//	length is a power of 2. usefull for FFT and co.

public:
	power2_array( size_t n ) : val1Darray<T,off>( nextPower2(n) )		{};
	power2_array( power2plus_array<T,off,1>& v ) : val1Darray<T,off>( v.firstPtr(), v.power2length(), abstract_valarray_0<T>::dontCopy )		{};
	//power2_array( const valarray0<T>& v, val1Darray<T,off>::TrendFlag f = val1Darray<T,off>::no )
	power2_array( const abstract_valarray_0<T>& v, typename abstract_valarray_0<T>::TrendFlag f = abstract_valarray_0<T>::no )
	//power2_array( const valarray0<T>& v, TrendFlag f = no )
    	: val1Darray<T,off>( nextPower2(v.size()) )
    	{	this->do_copy( v, f );	};
	//power2_array( const valarray0<T>& v, spectralwindowfunc& w, val1Darray<T,off>::TrendFlag f = val1Darray<T,off>::no )
	power2_array( const abstract_valarray_0<T>& v, spectralwindowfunc& w, typename abstract_valarray_0<T>::TrendFlag f = abstract_valarray_0<T>::no )
	//power2_array( const valarray0<T>& v, spectralwindowfunc& w, TrendFlag f = no )
    	: val1Darray<T,off>( nextPower2(v.size()) )
    	{	this->do_copy( v, f );
    		w.SetFlag( spectralwindowfunc::centered, spectralwindowfunc::half, v.size() );
	    	for (size_t i=1; i<=v.size(); i++)	(*this)(i) *= w.WindowValueAt(i);
    	};
	//power2_array( const valarray0<T>& v, size_t nb_zeros, val1Darray<T,off>::TrendFlag f = val1Darray<T,off>::no )
	power2_array( const abstract_valarray_0<T>& v, size_t nb_zeros, typename abstract_valarray_0<T>::TrendFlag f = abstract_valarray_0<T>::no )
	//power2_array( const valarray0<T>& v, size_t nb_zeros, TrendFlag f = no )	//	make sure there is at least 'nb_zeros' zeros at the end
    	: val1Darray<T,off>( nextPower2(v.size()+nb_zeros) )
    	{	this->do_copy( v, f );	};
    
	void	DoFFT()				{	realft( *this );		};		//	forward FFT
	void	DoBackwardFFT()		{	realft( *this, false );	};		//	backward FFT
	void	Filter( simplefunc* filt, double dx )					//	filter using a real filter function
				{	DoFFT();										//	first fourier transform the data
					size_t	no2 = val1Darray<T,off>::length()/2;
					double	fc = 1/(2*dx);							//	frequency f varies from 0 to fc = 1/(2*dx)
					(*this)[1] *= filt->ValueAt( 0 )/no2;	
					for (size_t i=2; i<=no2; i++)
					{	double f = (i-1)*fc/no2;
						double g = filt->ValueAt( f )/no2;
						(*this)[2*i-1] *= g;
						(*this)[2*i]   *= g;
					}
					(*this)[2] *= filt->ValueAt( fc )/no2;
					DoBackwardFFT();								//	back-transform the data
				};
	void	Filter( simplefunc* R_filt, simplefunc* I_filt, double dx )		//	idem, using a complex filter function
				{	DoFFT();
					size_t	no2 = val1Darray<T,off>::length()/2;
					double	fc = 1/(2*dx);
					(*this)[1] *= R_filt->ValueAt( 0 )/no2;						//	must stay real (output is real)
					for (size_t i=2; i<=no2; i++)
					{	double f = (i-1)*fc/no2;
						(*this)[2*i-1] *= R_filt->ValueAt( f )/no2;
						(*this)[2*i]   *= I_filt->ValueAt( f )/no2;
					}
					(*this)[2] *= R_filt->ValueAt( fc )/no2;					//	must stay real (output is real)
					DoBackwardFFT();
				};
	void	DoSinFFT()			{	sinft( *this );		};			//	forward FFT
	void	DoStaggeredCosFFT()	{	cosft2( *this, true );	};		//	forward FFT
	void	DoBackwardSinFFT()
		{	sinft( *this );
			double z = (double)2/(double)val1Darray<T,off>::length();			//	= power2length() ?
			*this *= z;
		};
	void	DoBackwardStaggeredCosFFT()	{	cosft2( *this, false );	};
};




template <class T, int off>
class filter_vector : public val1Darray<T,off>	{
private:
	size_t	r_len;				//	wrap around order:	c[0], c[-1], c[-2], ..., c[-l_length], c[r_length], ..., c[2], c[1]
								//	l_length = length() - r_length - 1;
	size_t	index0( int n )
		{	int		m = ( n <= 0 ? (-n) % (int)abstract_valarray_0<T>::size() : abstract_valarray_0<T>::size() - (n % abstract_valarray_0<T>::size()) );
			return	m;
		};	//	entre 0 et size()-1
	T*		addr_c( int n )		{	return	this->firstPtr() + index0(n);		};
	val1Darray<T,off>	vn;
	power2_array<T,off>	fft;
	
	void	setsize( size_t sz, val1Darray<T,off>& v )
		{	if (sz != v.size())
			{	v.resize_but_keep( sz );
				std::copy( addr_c(0), addr_c(-l_length()) + 1, v.firstPtr() );
				std::copy( addr_c(r_length()), addr_c(1) + 1, v.lastPtr() - r_length() + 1 );
		}	};
	void	do_fft( size_t sz )
		{	size_t	n = nextPower2(sz);
			if (n != fft.size())	{	setsize( n, fft );	realft( fft );	}
		};
	void	setfilter( size_t n, size_t n0 )
		{	vn.resize_but_keep( n );
			for (size_t i=1;  i<=n0-r_length()-1; 			i++)	vn[i] = 0;
			for (size_t i=n0-r_length(); i<=n0+l_length();	i++)	vn[i] = value(n0-i);
			for (size_t i=n0+l_length()+1; i<=n;			i++)	vn[i] = 0;
		};
		
public:
	void 	setfiltersize( size_t nl, size_t nr = 0 )
		{	valarray0<T>::resize_but_keep(nl+nr+1);	r_len = nr;		};
/*	void	flatten_order()
		{	val1Darray<T,off> temp( addr_c(r_len), r_len );
			std::copy_backward( addr_c(0)-1, addr_c(-l_length()) + 1, addr_c(1) );
			std::copy( temp.firstPtr(), temp.lastPtr() + 1, firstPtr() );
		}
*/	
public:
	filter_vector() : val1Darray<T,off>(), r_len(0), fft(0)	{};
	filter_vector( size_t nl, size_t nr = 0 ) : val1Darray<T,off>(nl+nr+1), r_len(nr), fft(0)	{};
//	size_t	zero_pad_length()		{	return	1 + std::max( r_length(), l_length() );	};
    size_t	zero_pad_length()		{	return	std::max( r_length(), l_length() );	};
	int		r_length()				{	return	r_len;					};
	int		l_length()				{	return	abstract_valarray_0<T>::size() - r_len - 1;		};
	T&		value( int i )			{	return	( *addr_c(i) );			};
	const	val1Darray<T,off>&		getFilter( size_t n, size_t n0 )	{	setfilter(n,n0);	return vn;	};
    const	power2_array<T,off>&	getFFTFilter( size_t n )			{	do_fft(n);			return fft;	};
};

template <class T, int off>
class regular_array : public val1Darray<T,off>	{
private:
	T		stp;
		
protected:	//	for descendents only
    explicit regular_array() : val1Darray<T,off>(0)	 {};
	void	init_regular( const T& xd, const T& xf )
		{	unsigned long n = abstract_valarray_0<T>::size();
			ThrowErrorIf (n <= 1);
			stp = (xf-xd)/(n-1);
			T* p0 = this->firstPtr();
			T* p  = this->firstPtr()+n-1;
			while (n > 1)	*p-- = xd + (--n)*stp;
			*p0 = xd;					//	[1] = xd exactement !!
		};
    		
public:		//	(xd, xf, n, dx) => (xd, xf, n) ou (xd, xf, dx) ou (xd, n, dx);
    explicit regular_array ( T xd, T xf, const unsigned long n )
    		: val1Darray<T,off>( n )								{	init_regular( xd, xf );	};
    explicit regular_array ( T xd, const unsigned long n, T dx )
    		: val1Darray<T,off>( n )								{	init_regular( xd, xd + (n-1)*dx );	};
    explicit regular_array ( T xd, T xf, T dx )
    		: val1Darray<T,off>( 1 + long_round((xf-xd)/dx) )		{	init_regular( xd, xf );	};
    
	T	step()	const {	return stp;	};
};

template <class T, int off>
class midpoints_array : public val1Darray<T,off>	{
private:
	const T*	fromPtr;
	void	make()
		{	unsigned long n = abstract_valarray_0<T>::size();
			ThrowErrorIf (n < 1);
			const T* pf = fromPtr+n;
			T* p = this->firstPtr()+n-1;
			for (; p>=this->firstPtr(); p--)	{	*p = *(pf--);	*p = (*p + *pf)/2;	}
		};
public:
    explicit midpoints_array ( const double_vector1& xx )
    		: val1Darray<T,off>( xx.size()-1 ) 	{	fromPtr = xx.firstPtr();	make();	};
    explicit midpoints_array ( const double* xx, unsigned long n  )				//	n = length of xx
    		: val1Darray<T,off>( n-1 )				{	fromPtr = xx;		make();	};
};

template <class T, int off>
class extra_midpoints_array : public val1Darray<T,off>	{		//	bool extra => extrapolation
private:
	const T*	fromPtr;
	void	make( bool extra )
		{	unsigned long n = abstract_valarray_0<T>::size();
			ThrowErrorIf (n < 1);
			const T* pf = fromPtr+n-2;
			T* p0 = this->firstPtr();
			T* p  = p0+n-1;
			if (extra)	{	*p = -(*--pf);	*p = (*p + 3*(*++pf))/2;	p--;	}
			else		{	*p-- = *pf;	}
			for (; p>p0; p--)	{	*p = *(pf--);	*p = (*p + *pf)/2;	}
			if (extra)	{	*p = -(*++pf);	*p = (*p + 3*(*--pf))/2;	}
			else		{	*p = *pf;	}
		};
public:
    explicit extra_midpoints_array ( const double_vector1& xx, bool extra = true )
    		: val1Darray<T,off>( xx.size()+1 ) 	{	fromPtr = xx.firstPtr();	make(extra);	};
    explicit extra_midpoints_array ( const double* xx, unsigned long n, bool extra = true  )				//	n = length of xx
    		: val1Darray<T,off>( n+1 )				{	fromPtr = xx;		make(extra);	};
};

template <class T, int off>
class savgol_filter : public filter_vector<T,off>	{
public:
	savgol_filter( int nl, int nr, int ld, int m ): filter_vector<T,off>( nl, nr )
		{	savgol( *this, nl, nr, ld, m );		};
};

template <class T, int offx, int offy>
class savgol_filter_matrix : public val1Darray<filter_vector<T,offx>,offy>	{
public:
	savgol_filter_matrix( int nl, int nr, int ld, int m, int nlm=0, int nrm=0 ): val1Darray<filter_vector<T,offx>,offy>( nlm+nrm+1 )
		{	int np = val1Darray<filter_vector<T,offx>,offy>::length();		//	= nlm+nrm+1
			for (int k=-nlm; k<=nrm; k++)
			{	int n = (np-k) % np;
				filter_vector<T,offx>& ck = (*this)[n+1];
				ck.setfiltersize( nl+k, nr-k );
				savgol( ck, nl+k, nr-k, ld, m );
			}
		};
};

typedef	regular_array<float,1>		float_regular_array;			//	usual vectors, from 1 to n
typedef	regular_array<double,1>		double_regular_array;

typedef	midpoints_array<float,1>		float_midpoints_array;
typedef	midpoints_array<double,1>		double_midpoints_array;
typedef	extra_midpoints_array<float,0>	float_extra_midpoints_array;
typedef	extra_midpoints_array<double,0>	double_extra_midpoints_array;

typedef	power2_array<double,1>				double_power2_array;
typedef	filter_vector<double,1>				double_filter_vector;
typedef	savgol_filter<double,1>				double_savgol_filter;
typedef	savgol_filter_matrix<double,1,1>	double_savgol_filter_matrix;


//	Savitzky-Golay filters
void	savgol( double_vector1& c, int nl, int nr, int ld, int m );
//void	savgolmatrix( double_savgol_filter_matrix& c, int nl, int nr, int ld, int m );
void	convlv_savgol( const double_vector& data, double_savgol_filter_matrix& respns, double_vector& answer,
			Boolean useFFT, Boolean remove_mean, double value );



//	FFT	algorithms

void	four1_n( a_double_vector1& data, size_t n, Boolean forward = true );
void	realft_n( a_double_vector1& data, size_t n, Boolean forward = true );
void	four1( double_power2_array& data, Boolean forward = true );
void	realft( double_power2_array& data, Boolean forward = true );
void	convlv_fft( const a_double_vector1& data, double_filter_vector& respns, double_vector1& answer, Boolean forward = true );
void	convlv( const a_double_vector1& data, double_filter_vector& respns, double_vector1& answer, Boolean renormalize_bounds = false );

void	correl( const a_double_vector1& data1, const a_double_vector1& data2, double_filter_vector& answer, double norm=0, Boolean remove_mean = true );
void	autocorrel( const a_double_vector1& data, double_filter_vector& answer, double norm=0, Boolean remove_mean = true );
void	correl_fft( const a_double_vector1& data1, const a_double_vector1& data2, double_filter_vector& answer, double norm=0, Boolean remove_mean = true );
void	autocorrel_fft( const a_double_vector1& data, double_filter_vector& answer, double norm=0, Boolean remove_mean = true );
				//	norm = 0 means "renormalize_bounds" ie "local_normalization" ie "N-i"

void	cosft1( a_double_power2plus1_array& y );
void	cosft2( double_power2_array& y, Boolean forward );
void	sinft( double_power2_array& y );

//	Max Entropy	algorithms
void	memcof( const a_double_vector1& data, double& xms, a_double_vector1& d );
double	evlmem( double fdt,  const a_double_vector1& d, double xms );




class identityfunc : public simplefunc {
public:
	identityfunc()	{};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new identityfunc();	}

	virtual double ValueAt( double t )						const {	return t;	};
	virtual double DerivValueAt( double )					const {	return 1;	};
	virtual double IntegrateBetween( double t1, double t2 )	const {	return sqr(t2-t1)/2;	};
	
	virtual FourCharCode code()	const	{	return	IdentityFunctionCode;		};
};



class polynomfunc : public simplefunc	{
private:
	double_vector0* 	_coef;
protected:
	polynomfunc()	{};			//	empty coef
	
	void		SetCoef( const double_vector0& c )	{	_coef = new double_vector0(c);	};
	void		SetCoef( size_t n )					{	_coef = new double_vector0(n);	};
	double&		Coeff(int i)	{	return (*_coef)(i);				};
public:
//	polynomfunc( const double* xx, unsigned long n ) 	:	coef( (double*)xx, n )	{};
//	polynomfunc( double_vector0& c ) 					:	coef( c )				{};

	polynomfunc( const double* xx, size_t n ) 	{	SetCoef( double_vector0(xx,n) );	};
	polynomfunc( const double_vector0& c ) 				{	SetCoef( c );	};
	
	
	int			degree()			const 	{	return _coef->length()-1;		};
	double		Coefficient(int i)	const	{	return (*_coef)(i);				};
	const double_vector0& Coeffs()	const	{	return *_coef;					};
	
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 		{	return new polynomfunc(*_coef);		};
	
	virtual double ValueAt( double t ) const 
		{	int n = degree();	double p = Coefficient(n);
			while (--n>=0)	p = p*t + Coefficient(n);
			return	p;
		};
	virtual double DerivValueAt( double t ) const 
		{	int n = degree();	double p = n*Coefficient(n);
			while (--n>0)	p = p*t + n*Coefficient(n);
			return	p;
		};
	virtual double IntegrateBetween( double t1, double t2 ) const 
		{	int n = degree();
			double p1 = Coefficient(n)/(n+1);	double p2 = p1;
			while (--n>=0)	{	double c = Coefficient(n)/(n+1);	p1 = p1*t1 + c;		p2 = p2*t2 + c;	}
			return	p2*t2 - p1*t1;
		};

	virtual FourCharCode code()	const	{	return	'Poly';		};
};


class gaussian_sum_func : public simplefunc	{			//	sum of Ai exp( -( (f-fi)/wi )2 )
protected:
	double_vector1 	coef;			//	Ai
	double_vector1 	freq;			//	fi
	double_vector1 	weight;			//	wi
public:
	gaussian_sum_func( double_vector1& c, double_vector1& f, double_vector1& w ) :	coef(c), freq(f), weight(w)		{};
	gaussian_sum_func( double_vector1& f, double_vector1& w ) 	:	coef(1.0,f.length()), /*coef(f.length(),1.0), */freq(f), weight(w)		{};			//	coefs are set to 1
	
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new gaussian_sum_func( coef, freq, weight );	};

	virtual double ValueAt( double t ) const 
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * exp( -sqr( (t-freq[n])/weight[n] ) ) : s += coef[n] );
			return	s;
		};
	virtual double DerivValueAt( double t ) const 
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * 2*(freq[n]-t)/sqr(weight[n]) * exp( -sqr( (t-freq[n])/weight[n] ) ) : s += 0 );
			return	s;
		};
	virtual double IntegrateBetween( double t1, double t2 ) const 
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * weight[n] * SqrtPiSur2 * (erf( (t2-freq[n])/weight[n] ) - erf( (t1-freq[n])/weight[n] )) : s += coef[n] * (t2-t1) );
			return	s;
		};
		
	virtual FourCharCode code()	const	{	return	GaussianSumFunctionCode;		};
};

class gaussian_filter_func : public gaussian_sum_func	{			//	idem 'gaussian_sum_func', with symetry and possible normalizations
protected:
	void	enforce_symetry()					//	add the symetric ones
		{	size_t	n0 = coef.length();
			size_t	n = 2 * n0;
			coef.resize_but_keep( n );	freq.resize_but_keep( n );	weight.resize_but_keep( n );
			for (size_t	i=1; i<=n0; i++)
			{	coef[n0+i] = coef[i];	freq[n0+i] = -freq[i];	weight[n0+i] = weight[i];	}
		};
	void	enforce_max_to_one()
		{	size_t	n = coef.length();
			double	maxX = -1e222;
			for (size_t	i=1; i<=n; i++)
			{	double	x = ValueAt(freq[i]);	if (x>maxX)	maxX = x;	}		//	FAUX: le max. n'est pas 'exactement' a freq(i) des qu'il y a plusieurs gaussiennes.
			for (size_t	i=1; i<=n; i++)
				coef[i] /= maxX;
		};
	void	enforce_all_to_one()				//	recomputes new coef. so that result at freq(j) is one for each j.
		{	size_t	n = coef.length();			//	MAIS: le max. n'est pas 'exactement' a freq(i) des qu'il y a plusieurs gaussiennes.
			double_matrix1	m(n,n);
			double_vector1	v(n);
			int_vector1		indx(n);
			for (size_t	i=1; i<=n; i++)
			{	for (size_t	j=1; j<=n; j++)
					m(i,j) = exp( -sqr( (freq[j]-freq[i])/weight[i] ) );
				v(i) = 1;
			}
			Boolean even;
			ludcmp( m, indx, even );
			lubksb( m, indx, coef );
		};
public:
	gaussian_filter_func( double_vector1& c, double_vector1& f, double_vector1& w, Boolean sym = false ) :	gaussian_sum_func( c, f, w )
		{	if (sym) enforce_symetry();
			enforce_max_to_one();
		};
	gaussian_filter_func( double_vector1& f, double_vector1& w, Boolean sym = false ) :	gaussian_sum_func( f, w )
		{	if (sym) enforce_symetry();
			enforce_all_to_one();
		};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new gaussian_filter_func( coef, freq, weight );	};
	
	virtual FourCharCode code()	const	{	return	GaussianFilterFunctionCode;		};
};


class interpolfunc : public simplefunc	{
protected:
	const temp_double_vector1 	vx;		//	permet d'encapsuler (double*)
	const temp_double_vector1 	vy;
	double_vector1 	rvx;				//	if no input vector(s)....
	double_vector1 	rvy;
protected:
	const a_double_vector1& 	x;		//	access to x
	const a_double_vector1& 	y;
	unsigned long klo;

protected:			//	copy constructors for copySelf()
//	interpolfunc (const interpolfunc& int_f)
//		  : rvx(int_f.rvx), rvy(int_f.rvy),
//			vx( int_f.vx.firstPtr(), int_f.vx.length(), double_vector::dontCopy ),
//			vy( int_f.vy.firstPtr(), int_f.vy.length(), double_vector::dontCopy ),
//			x( (&int_f.x == &int_f.rvx) ? rvx : ((&int_f.x == &int_f.vx) ? vx : int_f.x) ),
//			y( (&int_f.y == &int_f.rvy) ? rvy : ((&int_f.y == &int_f.vy) ? vy : int_f.y) )
//				{	klo = 0;	};
	interpolfunc (const interpolfunc&, const double_vector1& xx, const double_vector1& yy )
		  : vx(), vy(), x( xx ), y( yy )	{	klo = 0;	};
	
public:				
					//	don't make copies of input vectors, if any. Useful for fitting
					
	interpolfunc ( size_t n )
		:	vx(), vy(), rvx(n), rvy(n), 
			x( rvx ), y( rvy )		{	klo = 0;	};
	interpolfunc ( const double_vector1& xx )
		:	vx(), vy(), rvy(xx.length()), 
			x( xx ), y( rvy )		{	klo = 0;	};
			
					//	don't make copies of input vectors	xx, yy
					//	interpolfunc doesn't change the values of xx and yy: const
					
	interpolfunc ( const double* xx, const double* yy, unsigned long n )
		:	vx( (double*)xx, n ),
			vy( (double*)yy, n ),
			x( vx ), y( vy )		{	klo = 0;	};
	interpolfunc ( const double_vector1& xx, const double* yy )
		:	vx(),
			vy( (double*)yy, xx.length() ),
			x( xx ), y( vy )		{	klo = 0;	};
	interpolfunc ( const double* xx, const double_vector1& yy )
		:	vx( (double*&)xx, yy.length() ),
			vy(),
			x( vx ), y( yy )		{	klo = 0;	};
	interpolfunc ( const double_vector1& xx, const double_vector1& yy )
		:	vx(),
			vy(),
			x( xx ), y( yy )		{	klo = 0;	};
			
					//	make copies of input vectors	xx, yy
	interpolfunc ( const double_vector1& xx, const double_vector1& yy, Boolean )
		:	vx(), vy(), rvx(xx), rvy(yy), 
			x( rvx ), y( rvy )		{	klo = 0;	};
			
};

class lin_interpolfunc : public interpolfunc	{
protected:
	lin_interpolfunc(const double_vector1& xx) : interpolfunc(xx)	{};		//	for fitting
public:
	lin_interpolfunc (const double* xx, const double* yy, unsigned long n) 	: interpolfunc (xx, yy, n)	{};
	lin_interpolfunc (const double_vector1& xx, const double* yy) 			: interpolfunc (xx, yy)		{};
	lin_interpolfunc (const double_vector1& xx, const double_vector1& yy) 	: interpolfunc (xx, yy)		{};
	lin_interpolfunc (const double_vector1& xx, const double_vector1& yy, Boolean t) 	: interpolfunc (xx, yy, t)		{};
	
	virtual double ValueAt( double t )			const {	return	MyMath::linear_interp( x, y, t, (unsigned long&)klo );			};
	virtual double DerivValueAt( double t )		const {	return	MyMath::linear_interp_deriv1( x, y, t, (unsigned long&)klo );	};
	virtual double IntegrateBetween( double t1, double t2 ) const 
											{	return	MyMath::linear_interp_integ1( x, y, t1, t2, (unsigned long&)klo );	};
	
	lin_interpolfunc (const lin_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )
		: interpolfunc (int_f,xx,yy)	{};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new lin_interpolfunc(*this,xx,yy);	};
		
	virtual FourCharCode code()	const	{	return	LinearInterpFunctionCode;		};
};


class strictly_increasing_func	{			//	strictly increasing functions have an inverse function
public:
	virtual double BackValueAt( double t ) = 0;
};


class inc_lin_interpolfunc : public lin_interpolfunc, public strictly_increasing_func	{
public:
	inc_lin_interpolfunc (const double* xx, const double* yy, unsigned long n) 			: lin_interpolfunc (xx, yy, n)	{};
	inc_lin_interpolfunc (const double_vector1& xx, const double* yy) 					: lin_interpolfunc (xx, yy)		{};
	inc_lin_interpolfunc (const double_vector1& xx, const double_vector1& yy) 			: lin_interpolfunc (xx, yy)		{};
	inc_lin_interpolfunc (const double_vector1& xx, const double_vector1& yy, Boolean t): lin_interpolfunc (xx, yy, t)	{};
	inc_lin_interpolfunc (const lin_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )	: lin_interpolfunc (int_f,xx,yy)	{};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)	{	return new inc_lin_interpolfunc(*this,xx,yy);	};

	virtual double BackValueAt( double t )			{	return	MyMath::linear_interp( y, x, t, klo );			};
	
	virtual FourCharCode code()	const	{	return	IncreasingLinearInterpFunctionCode;		};
};

class lin_interpol_y_cyclicfunc : public lin_interpolfunc	{
public:
	double	ya, yb;
private:
	double	yc, dy;
	void init()		{	if (ya>yb)	swap(ya,yb);	yc = 0.5*(ya+yb);	dy = yb-ya;	};
public:
	lin_interpol_y_cyclicfunc( double yya, double yyb, const double* xx, const double* yy, unsigned long n )	: lin_interpolfunc (xx, yy, n),ya(yya),yb(yyb)		{	init();	};
	lin_interpol_y_cyclicfunc( double yya, double yyb, const double_vector1& xx, const double* yy ) 			: lin_interpolfunc (xx, yy),ya(yya),yb(yyb)			{	init();	};
	lin_interpol_y_cyclicfunc( double yya, double yyb, const double_vector1& xx, const double_vector1& yy ) 	: lin_interpolfunc (xx, yy),ya(yya),yb(yyb)			{	init();	};
	lin_interpol_y_cyclicfunc( double yya, double yyb, const double_vector1& xx, const double_vector1& yy, Boolean t ) 	: lin_interpolfunc (xx, yy, t),ya(yya),yb(yyb)	{	init();	};
	lin_interpol_y_cyclicfunc( double yya, double yyb, const lin_interpol_y_cyclicfunc& int_f, const double_vector1& xx, const double_vector1& yy )
																												: lin_interpolfunc (int_f,xx,yy),ya(yya),yb(yyb)	{	init();	};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new lin_interpol_y_cyclicfunc(ya,yb,*this,xx,yy);	};
	
	double	ValueModulo( double y ) const
		{	double n;
			double y0 = (y-ya);
			double y1 = (y0 > 0 ? dy * modf( y0/dy, &n ) : dy * (1+modf( y0/dy, &n )));		//	a priori  0 <= y1 < dy
			return ya + y1;
		};
	
	Boolean Continuous( unsigned long k1, double& t0 ) const
		{	//ThrowIf_ (k1 >= x.length() || k1 <= 0);
			ThrowErrorIf(k1 >= x.length() || k1 <= 0);
			unsigned long	k2 = k1+1;
			double	x1 = x(k1);
			double	y1 = y(k1);
			double	z1 = (y1 < yc ? y1+dy : y1-dy );
			double	x2 = x(k2);
			double	y2 = y(k2);
			if (abs(y2-y1) < abs(y2-z1))	return	true;
			else
			{	double	z0 = (y1 < yc ? yb : ya );
				t0 = x1 + (x2-x1)*(z1-z0)/(z1-y2);
				return false;
			};
		};
		
	virtual double ValueAt( double t ) const 
		{	double	H = bissectbracket( x, t, (unsigned long&)klo );
			//ThrowIf_ (H == 0);
			ThrowErrorIf (H == 0);
			unsigned long	khi = klo+1;
			double	t0;
			if (Continuous( klo, t0 ))
				return	((x(khi) - t) * y(klo) + (t - x(klo)) * y(khi)) / H;
			else
			{	double	z1 = (y(klo) < yc ? y(klo)+dy : y(klo)-dy );
				double	z2 = (y(khi) < yc ? y(khi)+dy : y(khi)-dy );
				if (t>t0)	return ((x(khi) - t) * z1 + (t - x(klo)) * y(khi)) / H;
				else		return ((x(khi) - t) * y(klo) + (t - x(klo)) * z2) / H;
			}
		};
	virtual double DerivValueAt( double t ) const 
		{	double	H = bissectbracket( x, t, (unsigned long&)klo );
			ThrowErrorIf (H == 0);
			unsigned long	khi = klo+1;
			double	t0;
			if (Continuous( (unsigned long&)klo, t0 ))
				return	(y(khi)-y(klo)) / H;
			else
			{	double	z1 = (y(klo) < yc ? y(klo)+dy : y(klo)-dy );
				double	z2 = (y(khi) < yc ? y(khi)+dy : y(khi)-dy );
				if (t>t0)	return (y(khi)-z1) / H;
				else		return (y(klo)-z2) / H;
			}
		};
	virtual double IntegrateBetween( double t1, double t2 ) const 				//	without "jumps" from starting point t1
		{	if (t2 < t1)	return (-IntegrateBetween( t2, t1 ));				//	make sure that t1 <= t2
			double	slin = MyMath::linear_interp_integ1( x, y, t1, t2, (unsigned long&)klo );	
			bissectbracket( x, t1, (unsigned long&)klo );
			unsigned long	klo1 = klo;
			bissectbracket( x, t2, (unsigned long&)klo );
			int i = 0;
			for (unsigned long k=klo1; k<=klo; k++)
			{	double	t0;
				double	dx = x(k+1)-x(k);
				if (t1 > x(k))		dx = x(k+1)-t1;
				if (t2 < x(k+1))	dx = t2-x(k);
				if (!Continuous( k, t0 ))		//	found a jump between position k and k+1 (precisely at t0)
				{	if (y(k+1) < yc)	{	slin += (i+0.5) * dx * dy;	i++;	}
					else				{	slin += (i-0.5) * dx * dy;	i--;	}
				}
				else slin += i * dx * dy;
			}
			return slin;
		};
	
	virtual double MeanBetween( double t1, double t2 ) const
		{	double	y = (t1 == t2 ?  ValueAt( t1 ) : IntegrateBetween( t1, t2 )/(t2 - t1));
			return	ValueModulo( y );
		};
				
	virtual FourCharCode code()	const	{	return	LinearInterpCyclicFunctionCode;		};
};

const 	double	unused_spline_yp = 2e30;

class splin_interpolfunc : public interpolfunc {

public:		//	static functions

	static	void	spline_init( const a_double_vector1& x, const a_double_vector1& y, a_double_vector1& y2, double yp1 = unused_spline_yp, double ypn = unused_spline_yp );
//	spline_init:	given a monotonically increasing vector x and a vector y, computes the second
//			derivative vector y2. yp1 and ypn are the first derivatives at the end points. If unused
//			the second derivative at the end points is set to zero (natural spline).

	static	double	spline_interp( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double t, unsigned long& klo );
//	spline_interp:	given a monotonically increasing vector xa, a vector ya and the second derivative vector y2a, 
//			given an abscissa t and an initial guess for the index klo such that xa[klo] near t, computes the spline
//			interpolation at t.

//	spline_interp_integ1:	...first integral between a and b
	static	double	spline_interp_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2,
			double a, double b, unsigned long& klo );
//	spline_interp_deriv1:	...first (second, third) derivative at t
	static	double	spline_interp_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double t, unsigned long& klo );
	static	double	spline_interp_deriv2( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo );
	static	double	spline_interp_deriv3( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo );

//	spline utilities:
//	retourne l'integrale de x1 a x2 de la spline definie par (x1,x2,y1,y2,ys1,ys2)
	inline static	double	spline_integ_X1_X2( const double& x1, const double& x2, const double& y1, const double& y2, const double& ys1, const double& ys2 )
		{	double dx = x2-x1;		return dx * ((y1 + y2)/2 - (ys1 + ys2) * dx*dx/24);		};
	inline static	double	spline_integ_X1_X2( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2_R, const a_double_vector1& y2_L, const unsigned long& i )
		{	return	spline_integ_X1_X2( x[i], x[i+1], y[i], y[i+1], y2_R[i], y2_L[i+1] );		};
//	retourne l'integrale de x1 a x de la spline definie par (x1,x2,y1,y2,ys1,ys2)
	inline static	double	spline_integ_X1_X( const double& x, const double& x1, const double& x2, const double& y1, const double& y2, const double& ys1, const double& ys2 )
		{	double dx = x2-x1;
			double d1 = x-x1;	double d12 = d1*d1;		double d14 = d12*d12;
			double d2 = x2-x;	double d22 = d2*d2;		double d24 = d22*d22;
			return (  y1 * (dx - d22/dx)/2 + y2 * d12/dx/2 + ys1 * ( dx*d22/12 - d24/dx/24 - dx*dx*dx/24 ) + ys2 * ( d14/dx/24 - dx*d12/12 )  );
		};
	inline static	double	spline_integ_X1_X( const double& v, const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2_R, const a_double_vector1& y2_L, const unsigned long& i )
		{	return	spline_integ_X1_X( v, x[i], x[i+1], y[i], y[i+1],  y2_R[i], y2_L[i+1] );	};
//	retourne l'integrale de x a x2 de la spline definie par (x1,x2,y1,y2,ys1,ys2)
	inline static	double	spline_integ_X_X2( const double& x, const double& x1, const double& x2, const double& y1, const double& y2, const double& ys1, const double& ys2 )
		{	return ( spline_integ_X1_X2( x1, x2, y1, y2, ys1, ys2 ) - spline_integ_X1_X( x, x1, x2, y1, y2, ys1, ys2 ) );	};
	inline static	double	spline_integ_X_X2( const double& v, const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& y2_R, const a_double_vector1& y2_L, const unsigned long& i )
		{	return	spline_integ_X_X2( v, x[i], x[i+1], y[i], y[i+1],  y2_R[i], y2_L[i+1] );	};

//		idem bissectbracket, avec controle d'extrapolation (extra = 0 si interpol; = -1 si extrapol gauche; = +1 si extrapol droite)
	inline static	double	spline_biss_bracket( const a_double_vector1& x, const double& t, unsigned long& klo, unsigned long& khi, int& extra )
		{	hunt( x, t, klo );			//	output: 0 <= klo <= n
			unsigned long n = x.length();
			extra = 0;
			if ( klo >= n )			{	klo = n-1;	extra = 1;	}
			else if ( klo <= 0 )	{	klo = 1;	extra = -1;	}
			//	output: 1 <= klo <= n-1
			khi = klo+1;
			double	H = ( x[ khi ] - x[ klo ] );			//	always OK, except for n = 1
			ThrowErrorIf (H == 0);
			return	H;
		};
		
protected:
	double_vector1	dersec;
	splin_interpolfunc(const double_vector1& xx) : interpolfunc(xx), dersec(xx.length())	{};		//	for fitting
	void	spline_initialize()		{	spline_init( x, y, dersec );	};
	
public:
	splin_interpolfunc (const double* xx, const double* yy, unsigned long n) 	: interpolfunc (xx, yy, n), dersec(n)
			{	spline_initialize();	};
	splin_interpolfunc (const double_vector1& xx, const double* yy) 			: interpolfunc (xx, yy), dersec(xx.length())
			{	spline_initialize();	};
	splin_interpolfunc (const double_vector1& xx, const double_vector1& yy) 	: interpolfunc (xx, yy), dersec(xx.length())
			{	spline_initialize();	};			
	splin_interpolfunc (const splin_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy ) : interpolfunc (int_f,xx,yy), dersec(xx.length())
			{	spline_initialize();	};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new splin_interpolfunc(*this,xx,yy);	};
	
	virtual double ValueAt( double t )			const {	return	spline_interp( x, y, dersec, t, (unsigned long&)klo );			};
	virtual double DerivValueAt( double t )		const {	return	spline_interp_deriv1( x, y, dersec, t, (unsigned long&)klo );	};
	virtual double IntegrateBetween( double t1, double t2 )
										const 	{	return	spline_interp_integ1( x, y, dersec, t1, t2, (unsigned long&)klo );	};
										
	virtual FourCharCode code()	const	{	return	SplineInterpFunctionCode;		};
};


class splin_extra_lin_interpolfunc : public splin_interpolfunc {
			// idem "splin_interpolfunc", but extrapolation is linear for natural splines (if y" == 0 at the ends)
public:
		//	implementation is for "inc_splin_interpolfunc", with left and right derivatives		
	static	double	spline_interp_extra_lin( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end = true );
	static	double	spline_interp_extra_lin_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end = true );
	static	double	spline_interp_extra_lin_deriv2( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end = true );
	static	double	spline_interp_extra_lin_deriv3( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double X, unsigned long& klo, bool dont_check_end = true );
	static	double	spline_interp_extra_lin_integ1( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& Y2A_R, const a_double_vector1& Y2A_L, double a, double b, unsigned long& klo, bool dont_check_end = true );

	inline	static	double	spline_interp_extra_lin( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo, bool dont_check_end = false )
		{	return	spline_interp_extra_lin( XA, YA, Y2A, Y2A, X, klo, dont_check_end );	};
	inline	static	double	spline_interp_extra_lin_deriv1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo, bool dont_check_end = false )
		{	return	spline_interp_extra_lin_deriv1( XA, YA, Y2A, Y2A, X, klo, dont_check_end );	};
	inline	static	double	spline_interp_extra_lin_deriv2( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo, bool dont_check_end = false )
		{	return	spline_interp_extra_lin_deriv2( XA, YA, Y2A, Y2A, X, klo, dont_check_end );	};
	inline	static	double	spline_interp_extra_lin_deriv3( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double X, unsigned long& klo, bool dont_check_end = false )
		{	return	spline_interp_extra_lin_deriv3( XA, YA, Y2A, Y2A, X, klo, dont_check_end );	};
	inline	static	double	spline_interp_extra_lin_integ1( const a_double_vector1& XA, const a_double_vector1& YA, const a_double_vector1& Y2A, double a, double b, unsigned long& klo, bool dont_check_end = false )
		{	return	spline_interp_extra_lin_integ1( XA, YA, Y2A, Y2A, a, b, klo, dont_check_end );	};
		
public:
	splin_extra_lin_interpolfunc (const double* xx, const double* yy, unsigned long n) 	: splin_interpolfunc (xx, yy, n)	{};		//	spline_initialize();	};
	splin_extra_lin_interpolfunc (const double_vector1& xx, const double* yy) 			: splin_interpolfunc (xx, yy)		{};		//	spline_initialize();	};
	splin_extra_lin_interpolfunc (const double_vector1& xx, const double_vector1& yy) 	: splin_interpolfunc (xx, yy)		{};		//	spline_initialize();	};			
	splin_extra_lin_interpolfunc (const splin_extra_lin_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )	: splin_interpolfunc (int_f,xx,yy)	{};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)	{	return new splin_extra_lin_interpolfunc(*this,xx,yy);	};
	
	virtual double ValueAt( double t )			const {	return	spline_interp_extra_lin( x, y, dersec, t, (unsigned long&)klo );			};
	virtual double DerivValueAt( double t )		const {	return	spline_interp_extra_lin_deriv1( x, y, dersec, t, (unsigned long&)klo );	};
	virtual double IntegrateBetween( double t1, double t2 )
										const 	{	return	spline_interp_integ1( x, y, dersec, t1, t2, (unsigned long&)klo );	};
										
	virtual FourCharCode code()	const	{	return	SplineExtraLinInterpFunctionCode;		};
};

class inc_splin_interpolfunc : public splin_extra_lin_interpolfunc, public strictly_increasing_func	{
public:
	typedef enum	{ fast, precise }	inc_sp_type;

	static	Boolean	inc_spline_check( const a_double_vector1& x, const a_double_vector1& y, const a_double_vector1& yR, const a_double_vector1& yL, double eps );
	static	void	inc_spline_init( const a_double_vector1& x, const a_double_vector1& y, a_double_vector1& y2L, a_double_vector1& y2R, double minRate, int t );
/*
//	below:	idem as in "splin_interpolfunc", but with left and right second derivatives
	static	double	spline_interp( const double_vector1& x, const double_vector1& y, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double t, unsigned long& klo );
	static	double	spline_interp_deriv1( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double t, unsigned long& klo );
	static	double	spline_interp_deriv2( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double X, unsigned long& klo );
	static	double	spline_interp_deriv3( const double_vector1& XA, const double_vector1& YA, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double X, unsigned long& klo );
	static	double	spline_interp_integ1( const double_vector1& x, const double_vector1& y, const double_vector1& Y2A_R, const double_vector1& Y2A_L, double a, double b, unsigned long& klo );
*/
private:
	double_vector1&	dersecL;
	double_vector1	dersecR;
	double			minimalRate;
	inc_sp_type		sp_type;
	
	void	SetMinimalRate( double minR )
		{	if (minR > 0)	minimalRate = minR;
			else			minimalRate = 10e-6 * (y.lastPtr()-y.firstPtr())/(x.lastPtr()-x.firstPtr());
		}
		
	void	inc_spline_initialize()		{	inc_spline_init( x, y, dersecL, dersecR, minimalRate, sp_type );	};
	
public:				//	don't make copies of input vectors
	
	inc_splin_interpolfunc (const double* xx, const double* yy, unsigned long n, double minR = 0, inc_sp_type t = precise)
		:	splin_extra_lin_interpolfunc (xx, yy, n),
			dersecL( dersec ),						//	left one is not a copy, but a reference to dersec
			dersecR ( dersec ),						//	right one is a copy
			sp_type(t)
			{	SetMinimalRate( minR );	inc_spline_initialize();	};
	inc_splin_interpolfunc (const double_vector1& xx, const double* yy, double minR = 0, inc_sp_type t = precise)
		:	splin_extra_lin_interpolfunc (xx, yy),
			dersecL( dersec ),						//	left one is not a copy, but a reference to dersec
			dersecR ( dersec ),						//	right one is a copy
			sp_type(t)
			{	SetMinimalRate( minR );	inc_spline_initialize();	};
	inc_splin_interpolfunc (const double_vector1& xx, const double_vector1& yy, double minR = 0, inc_sp_type t = precise)
		:	splin_extra_lin_interpolfunc (xx, yy),
			dersecL( dersec ),						//	left one is not a copy, but a reference to dersec
			dersecR ( dersec ),						//	right one is a copy
			sp_type(t)
			{	SetMinimalRate( minR );	inc_spline_initialize();	};
	inc_splin_interpolfunc (const inc_splin_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )
		: splin_extra_lin_interpolfunc (int_f,xx,yy), dersecL(dersec), dersecR(int_f.dersecR), minimalRate(int_f.minimalRate), sp_type(int_f.sp_type)
			{	inc_spline_initialize();	};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new inc_splin_interpolfunc(*this,xx,yy);	};
	
	virtual double ValueAt( double t )			const {	return	spline_interp_extra_lin( x, y, dersecR, dersecL, t, (unsigned long&)klo, true );		};
	virtual double DerivValueAt( double t )		const {	return	spline_interp_extra_lin_deriv1( x, y, dersecR, dersecL, t, (unsigned long&)klo, true );	};
	virtual double IntegrateBetween( double t1, double t2 )
												const {	return	spline_interp_extra_lin_integ1( x, y, dersecR, dersecL, t1, t2, (unsigned long&)klo, true );	};
	
	
	
	class	solveFunction : public doubleFunction {
	private:	
		doubleFunction& f;
		const double x;
	public:	
		solveFunction( doubleFunction& f0, const double x0 ):f(f0), x(x0) {};
	protected:	virtual double	func( const double t )		const {	return	f(t) - x;	};
	};
	
	virtual double BackValueAt( double t )
		{	unsigned long	khi;
			int				extra;
			double			hy = spline_biss_bracket( y, t, klo, khi, extra );
			double			hx = (x[khi] - x[klo]);
			
			if ( extra == 1 )	//	&& dersecL[khi] == 0 )
			{	double	dy = hy / hx + ( dersecR[klo] + 2 * dersecL[khi]) * hx / 6;		//	derivative at the end
				return	x[khi] + (t - y[khi]) / dy;
			}
			if ( extra == -1 )	//	&& dersecR[klo] == 0 )
			{	double	dy = hy / hx - (2 * dersecR[klo] + dersecL[khi]) * hx / 6;		//	derivative at the beginning
				return	x[klo] + (t - y[klo]) / dy;
			}
			solveFunction	sol_f( *this, t );
			return	MyMath::zbrent( sol_f, x[klo], x[khi], 1.0e-6 );
		};
		
	virtual FourCharCode code()	const	{	return	IncreasingSplineInterpFunctionCode;		};
};


class abstract_stair_interpolfunc : public interpolfunc{		//	implementation is: "stair_start_interpolfunc"
protected:
//	virtual const a_double_vector1& x_start() const = 0;
//	virtual const a_double_vector1& y_start() const = 0;
	virtual double value_( size_t klo ) const = 0;
public:
	abstract_stair_interpolfunc (const double* xx, const double* yy, unsigned long n)		: interpolfunc (xx, yy, n)	{};
	abstract_stair_interpolfunc (const double_vector1& xx, const double* yy)				: interpolfunc (xx, yy)		{};
	abstract_stair_interpolfunc (const double_vector1& xx, const double_vector1& yy)		: interpolfunc (xx, yy)		{};
	abstract_stair_interpolfunc (const abstract_stair_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy ) : interpolfunc (int_f,xx,yy)	{};

	virtual double ValueAt( double t ) const
		{		MyMath::hunt( x, t, (unsigned long&)klo );
				//return  ( klo > 0 ? y[klo] : y[1] );
				return	value_( klo );
		};
	virtual double DerivValueAt( double ) const {	return	0;	};
	virtual double IntegrateBetween( double t1, double t2 ) const
		{	//		ia, ib	:	indices min max tels que t1 < x[ia] <...< x[ib] < t2
			unsigned long	n = y.length();
			double			s = 0;
			
			MyMath::hunt( x, t1, (unsigned long&)klo );
			unsigned long	ia = std::min( klo+1, n );							//	x[ia] juste apres t1 ou bien juste avant (si t1>x[n])
			s += value_( klo ) * (x[ia] - t1);									//	Sum de t1 -> x[ia]
			
			MyMath::hunt( x, t2, (unsigned long&)klo );
			unsigned long	ib = std::max( klo, (size_t)1 );							//	x[ib] juste avant t2 ou bien juste apres (si t2<x[1])
			s += value_( klo ) * (t2 - x[ib]);									//	Sum de x[ib] -> t2
	
			if (ib < ia)	for (unsigned long i=ib; i<=(ia-1); i++)
								s -= (x[i+1] - x[i]) * value_(i);
			else			for (unsigned long i=ia; i<=(ib-1); i++)
								s += (x[i+1] - x[i]) * value_(i);				//	Sum de x[ia] -> x[ib]
			return s;
		};
};

class stair_start_interpolfunc : public abstract_stair_interpolfunc {		//	{xi,yi}  f(x) = yi si x(i) <= x < x(i+1);  ou f(x) = y1 si x <= x1
																			//	NOTE : x1 est inutilisé :  f(x) = y1 si x <= x(2)
protected:
//	virtual const a_double_vector1& x_start()	const {	return x;	};
//	virtual const a_double_vector1& y_start()	const {	return y;	};
	
	virtual double value_( size_t klo ) const	{	return  ( klo > 0 ? y[klo] : y[1] );	};
//protected:
//	stair_start_interpolfunc(const double_vector1& xx) : interpolfunc (xx)		{};		//	for fitting
public:
	stair_start_interpolfunc (const double* xx, const double* yy, unsigned long n)		: abstract_stair_interpolfunc (xx, yy, n)	{};
	stair_start_interpolfunc (const double_vector1& xx, const double* yy)				: abstract_stair_interpolfunc (xx, yy)		{};
	stair_start_interpolfunc (const double_vector1& xx, const double_vector1& yy)		: abstract_stair_interpolfunc (xx, yy)		{};
	stair_start_interpolfunc (const stair_start_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )
		: abstract_stair_interpolfunc (int_f,xx,yy)	{};

	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new stair_start_interpolfunc(*this,xx,yy);	};
		
		
/*	virtual double ValueAt( double t ) const
		{		MyMath::hunt( x, t, (unsigned long&)klo );
				return  ( klo > 0 ? y[klo] : y[1] );
		};
	virtual double DerivValueAt( double ) const		{	return	0;	};
	virtual double IntegrateBetween( double t1, double t2 ) const
		{	//		ia, ib	:	indices min max tels que t1 < x[ia] <...< x[ib] < t2
			unsigned long	n = x.length();
			double			s = 0;
			MyMath::hunt( x, t1, (unsigned long&)klo );
			unsigned long	ia = klo+1;							//	x[ia] juste apres a
			if (ia == n)	ia--;								//		ou bien juste avant...
			s += ( klo > 0 ? y[klo] : y[1] ) * (x[ia] - t1);	//	Sum de t1 -> x[ia]
			MyMath::hunt( x, t2, (unsigned long&)klo );
			unsigned long	ib = klo;							//	x[ib] juste avant b
			if (ib == 0)	ib++;								//		ou bien juste apres...
			s += ( klo > 0 ? y[klo] : y[1] ) * (t2 - x[ib]);	//	Sum de x[ib] -> t2
	
			if (ib < ia)	for (unsigned long i=ib; i<=(ia-1); i++)
								s -= (x[i+1] - x[i]) * y[i];
			else			for (unsigned long i=ia; i<=(ib-1); i++)
								s += (x[i+1] - x[i]) * y[i];		//	Sum de x[ia] -> x[ib]
			return s;
		};
*/
	virtual FourCharCode code()	const	{	return	StairStartInterpFunctionCode;		};
};


class stair_end_interpolfunc : public abstract_stair_interpolfunc {		//	{xi,yi}  f(x) = yi si x(i-1) <= x < x(i);  ou f(x) = yn si x > x(n)
																			//	NOTE : xn est inutilisé :  f(x) = yn si x > x(n-1)
protected:
//	virtual const double_vector1& x_start()	const {	return x;	};
//	virtual const double_vector1& y_start()	const {	return y;	};
	virtual double value_( size_t klo ) const	{	return  ( klo < x.size() ? y[klo+1] : y[klo] );	};
public:
	stair_end_interpolfunc (const double* xx, const double* yy, unsigned long n)	: abstract_stair_interpolfunc (xx, yy, n)	{};
	stair_end_interpolfunc (const double_vector1& xx, const double* yy)				: abstract_stair_interpolfunc (xx, yy)		{};
	stair_end_interpolfunc (const double_vector1& xx, const double_vector1& yy)		: abstract_stair_interpolfunc (xx, yy)		{};
	stair_end_interpolfunc (const stair_end_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )
		: abstract_stair_interpolfunc (int_f,xx,yy)	{};

	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new stair_end_interpolfunc(*this,xx,yy);	};
		
//	virtual double ValueAt( double t ) const
//		{		MyMath::hunt( x, t, (unsigned long&)klo );
//				return  ( klo < x.size() ? y[klo+1] : y[klo] );
//		};
		
	virtual FourCharCode code()	const	{	return	StairEndInterpFunctionCode;		};
};




class stair_interpolfunc : public interpolfunc{
private:
	double_midpoints_array	midPoints;
protected:
	stair_interpolfunc(const double_vector1& xx) : interpolfunc (xx), midPoints(xx)		{};		//	for fitting
protected:
	virtual const double_vector1& x_midPoints()	const {	return midPoints;	};
public:
	stair_interpolfunc (const double* xx, const double* yy, unsigned long n)
		: interpolfunc (xx, yy, n), midPoints( xx, n )	{};
	stair_interpolfunc (const double_vector1& xx, const double* yy)					//: interpolfunc (xx, yy), midPoints( xx, true )	{};
		: interpolfunc (xx, yy), midPoints( xx )	{};
	stair_interpolfunc (const double_vector1& xx, const double_vector1& yy)			//: interpolfunc (xx, yy), midPoints( xx, true )	{};
		: interpolfunc (xx, yy), midPoints( xx )	{};
	stair_interpolfunc (const stair_interpolfunc& int_f, const double_vector1& xx, const double_vector1& yy )
		: interpolfunc (int_f,xx,yy), midPoints(int_f.midPoints)	{};
	virtual simplefunc*	copySelf(const double_vector1& xx, const double_vector1& yy)
		{	return new stair_interpolfunc(*this,xx,yy);	};
		
	virtual double ValueAt( double t )			const {	return	MyMath::stair_interp( x, y, midPoints, t, (unsigned long&)klo );	};
	virtual double DerivValueAt( double )		const {	return	0;													};
	virtual double IntegrateBetween( double t1, double t2 )
											const {	return	MyMath::stair_interp_integ1( x, y, midPoints, t1, t2, (unsigned long&)klo );	};

	virtual FourCharCode code()	const	{	return	StairInterpFunctionCode;		};
};


class fit_func {		//	pure virtual, to be added to a simplefunc descendant
public:
	typedef enum { stdFit, absFit } fittingFlag;
//protected:
	double		chisq;		//	or absv if absFit
	fittingFlag	type;
public:	
	fit_func( fittingFlag t )	{	type = t;	chisq = 0;	};
	
		//	if simplefunc depends on parameter set ((a),(b)) with set (a) fixed and set (b) adjustable (in fitting)
		//	fitF( x, b ) find the values of b at x for all basis functions b1=(1,0,0,..), b2=(0,1,0..),....

	virtual void	fitF( double, double_vector1& ) = 0;		//	{};
};

class polynom_fit_func : public polynomfunc, public fit_func	{
private:
	void	Make( fit_func::fittingFlag f, const double_vector1& x, const double_vector1& y, const double_vector1& sigy, int d )
		//	note: coef is not initialized
		{	size_t	n = d+1;	//	coef.length();
			if (f == stdFit)
			{	int ndata = x.length();
				double_matrix1 u(ndata,n);
				double_matrix1 v(n,n);
				double_vector1 w(n);
				double_vector1 c1(n);
				double chisq1;
				svdfit( x, y, sigy, c1, u, v, w, chisq1, this );
				
				double_vector1	c2(n);
				int_vector1 	ic(n);
				for (size_t k=1; k<=n; k++)	ic[k] = 1;
				double chisq2;
				lfit( x, y, sigy, c2, ic, v, chisq2, this );		//	simple linear algorithm is sometimes better (worth trying)
				
				SetCoef( (chisq1 > chisq2 ? c2 : c1) );				//	keep the best one !!
				chisq = (chisq1 > chisq2 ? chisq2 : chisq1);
			}
			else if (f == absFit)
			{	double_vector1 c1(n);
				absfit( x, y, sigy, c1, chisq, this );
				SetCoef( c1 );				//	cannot pass coef to absFit (index starts at 0)
			}
		};
		
public:
	polynom_fit_func( fit_func::fittingFlag f, const double_vector1& x, const double_vector1& y, const double_vector1& sigy, int d ) : fit_func(f)
		{	Make( f, x, y, sigy, d );	};
	polynom_fit_func( fit_func::fittingFlag f, const double_vector1& x, const double_vector1& y, double sy, int d ) : fit_func(f)
		{	double_vector sigy( x.length() );	sigy.set_all_to(sy);
			Make( f, x, y, sigy, d );
		};
	virtual void	fitF( double, double_vector1& );
		//	fitF( x, b ) find the values of b at x for all basis functions b1=(1,0,0,..), b2=(0,1,0..),....
		//	here: (1,x,x2,x3,...xn)
};

class line_fit_func : public polynomfunc, public fit_func	{		//	special case of mean-square straight-line fit, with more options
private:
	static	double_vector1	chixy_xx, chixy_yy, chixy_sx, chixy_sy, chixy_ww;		//	for routine chixy
	static	double			chixy_aa, chixy_offs;									//	for routine chixy
	static	double			chixy( double b_ang );									//	function for fitexy
	
	void	fit( const double_vector1& x, const double_vector1& y, const double_vector1& sig,
				Boolean hasSig, double& a, double& b );	//, double& siga, double& sigb, double& chi2, double& q );		//	standard linear fit
	

	void	fitexy( const double_vector1& x, const double_vector1& y, const double_vector1& sigx, const double_vector1& sigy, double& a, double& b );
public:
	double	siga, sigb, q;
	
	line_fit_func( const double_vector1& x, const double_vector1& y )
		: fit_func( stdFit )
		{	SetCoef(2);
			double_vector1 sig;		//	0-length
			fit( x, y, sig, false, Coeff(0), Coeff(1) );
		}
	line_fit_func( const double_vector1& x, const double_vector1& y, const double_vector1& sig )
		: fit_func( stdFit )
		{	SetCoef(2);
			fit( x, y, sig, true, Coeff(0), Coeff(1) );
		}
	line_fit_func( const double_vector1& x, const double_vector1& y, double sy )
		: fit_func( stdFit )
		{	SetCoef(2);
			double_vector sig( x.length() );	sig.set_all_to(sy);
			fit( x, y, sig, true, Coeff(0), Coeff(1) );
		}
		
/*	line_fit_func( const double_vector1& x, const double_vector1& y, double lambda, const double_vector1& sigy ) 	 // sigx = lambda sigy
		: fit_func( stdFit )
		{	coef.resize_but_keep(2);
		}
*/		
	line_fit_func( const double_vector1& x, const double_vector1& y, const double_vector1& sigx, const double_vector1& sigy )
		: fit_func( stdFit )
		{	SetCoef(2);
			fitexy( x, y, sigx, sigy, Coeff(0), Coeff(1) );
		}
	line_fit_func( const double_vector1& x, const double_vector1& y, double sx, double sy )
		: fit_func( stdFit )
		{	SetCoef(2);
			double_vector sigx( x.length() );	sigx.set_all_to(sx);
			double_vector sigy( x.length() );	sigy.set_all_to(sy);
			fitexy( x, y, sigx, sigy, Coeff(0), Coeff(1) );
		}

	virtual void	fitF( double, double_vector1& )	{};		//	unused in this case
};



class interpol_fit_func : public fit_func	{
public:
	void	Make( const a_double_vector1& x, a_double_vector1& y, fit_func::fittingFlag f, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		{	int n = x.length();
			size_t	i_start = val1Darray_box<double,1>::make_first_index_sup( xn, x[1] );
			size_t	i_end   = val1Darray_box<double,1>::make_last_index_inf( xn, x[n] );
			val1Darray_box<double,1>	xxn( xn, i_start, i_end );
			val1Darray_box<double,1>	yyn( yn, i_start, i_end );
			val1Darray_box<double,1>	ssy( sigy, i_start, i_end );
			int ndata = xxn.length();
			
			if (f == stdFit)
			{	double_matrix1 u(ndata,n);
				double_matrix1 v(n,n);
				double_vector1 w(n);
				svdfit( xxn, yyn, ssy, y, u, v, w, chisq, this );
			}
			else if (f == absFit)
			{	absfit( xxn, yyn, ssy, y, chisq, this );
		}	};
	
	interpol_fit_func( fit_func::fittingFlag t ) : fit_func(t)		{};
};

class lin_interpol_fit_func : public lin_interpolfunc, public interpol_fit_func	{
public:
	lin_interpol_fit_func ( fit_func::fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		: lin_interpolfunc ( xx ), interpol_fit_func(f)
		{	Make( x, rvy, f, xn, yn, sigy );
		};
	lin_interpol_fit_func ( fit_func::fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, double sy = 1.0 )
		: lin_interpolfunc (xx), interpol_fit_func(f)
		{	double_vector sigy( xn.length() );	sigy.set_all_to(sy);
			Make( x, rvy, f, xn, yn, sigy );
		};
	virtual void	fitF( double x, double_vector1& f );
};

class splin_interpol_fit_func : public splin_interpolfunc, public interpol_fit_func	{
private:
	void	MakeBasisF()	// fittingFlag f, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		{		//	construct the set of basis functions to be used by fitF()
			size_t n = x.length();
			//val1Darray<double_vector1,1>	basisDersec(n,n);
			basisDersecPtr = new val1Darray<double_vector1,1>(n,n);
			double_vector1 ty(n);
			for (size_t i=1; i<=n; i++)
			{	for (size_t j=1; j<=n; j++)		ty[j] = 0;
				ty[i] = 1.0;
				spline_init( x, ty, (*basisDersecPtr)[i] );
			}
		};
	void	CleanBasisF()
	{	delete basisDersecPtr;	};
	
/*	void	Make( fittingFlag f, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		{		//	construct the set of basis functions to be used by fitF()
			size_t n = x.length();
			val1Darray<double_vector1,1>	basisDersec(n,n);	basisDersecPtr = &basisDersec;
			double_vector1 ty(n);
			for (size_t i=1; i<=n; i++)
			{	for (size_t j=1; j<=n; j++)		ty[j] = 0;
				ty[i] = 1.0;
				spline_init( x, ty, basisDersec[i] );
			}
				
			if (f == stdFit)
			{	int ndata = xn.length();
				double_matrix1 u(ndata,n);
				double_matrix1 v(n,n);
				double_vector1 w(n);
				svdfit( xn, yn, sigy, rvy, u, v, w, chisq, this );
				spline_initialize();	//	spline_init( x, y, dersec );
			}
			else if (f == absFit)
			{	absfit( xn, yn, sigy, rvy, chisq, this );
				spline_initialize();	//	spline_init( x, y, dersec );
			}
		};*/
public:
	splin_interpol_fit_func ( fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		: splin_interpolfunc (xx), interpol_fit_func(f)
		{	MakeBasisF();
			Make( x, rvy, f, xn, yn, sigy );
			spline_initialize();
			CleanBasisF();
		};
	splin_interpol_fit_func ( fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, double sy = 1.0 )
		: splin_interpolfunc (xx), interpol_fit_func(f)
		{	double_vector sigy( xn.length() );	sigy.set_all_to(sy);
			MakeBasisF();
			Make( x, rvy, f, xn, yn, sigy );
			spline_initialize();
			CleanBasisF();
		};
	virtual void	fitF( double x, double_vector1& f );

private:
	val1Darray<double_vector1,1>*	basisDersecPtr;
};

class stair_interpol_fit_func : public stair_interpolfunc, public interpol_fit_func	{
/*private:
	void	Make( fittingFlag f, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		{	if (f == stdFit)
			{	int ndata = xn.length();
				int n = x.length();
				double_matrix1 u(ndata,n);
				double_matrix1 v(n,n);
				double_vector1 w(n);
				svdfit( xn, yn, sigy, rvy, u, v, w, chisq, this );
			}
			else if (f == absFit)
			{	absfit( xn, yn, sigy, rvy, chisq, this );
			}
		};*/
public:
	stair_interpol_fit_func ( fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, const double_vector1& sigy )
		: stair_interpolfunc (xx), interpol_fit_func(f)
		{	Make( x, rvy, f, xn, yn, sigy );
			/*if (f == stdFit)
			{	int ndata = xn.length();
				int n = x.length();
				double_matrix1 u(ndata,n);
				double_matrix1 v(n,n);
				double_vector1 w(n);
				svdfit( xn, yn, sigy, rvy, u, v, w, chisq, this );
			}
			else if (f == absFit)
			{	absfit( xn, yn, sigy, rvy, chisq, this );
			}*/
		};
	stair_interpol_fit_func ( fittingFlag f, const double_vector1& xx, const double_vector1& xn, const double_vector1& yn, double sy = 1.0 )
		: stair_interpolfunc (xx), interpol_fit_func(f)
		{	double_vector sigy( xn.length() );	sigy.set_all_to(sy);
			Make( x, rvy, f, xn, yn, sigy );
		};
	virtual void	fitF( double x, double_vector1& f );
};



double	LinearCorrelation(  const double_vector1& x1,  const double_vector1& y1,  const double_vector1& x2,  const double_vector1& y2 );


class common_regular_scale : public double_regular_array {
private:
	bool	 empty_intersect;
	bool	 failed;
	bool	 just_copy;
	
public:
	common_regular_scale( const val1Darray<double_vector1*,1>& in_x, bool on_intersection = true )
/*	{
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
	}*/;
	
	bool	 empty_intersection()	{	return empty_intersect;	};
	bool	 has_failed()			{	return failed;			};
	bool	 is_just_a_copy()		{	return just_copy;		};
};

//inline bool		is_NaN( double x )		{	return	__gnu_cxx::__capture_isnan(x);	};
inline bool		is_NaN( double x )		{	return	(x != x);	};
inline double	NaN()					{	return	std::numeric_limits<double>::quiet_NaN();	};

}		//	MyMath




#endif	//	__MyMath_H_
