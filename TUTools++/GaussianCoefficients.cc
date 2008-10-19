/*
 *  $BJ?@.(B14-19$BG/!JFH!K;:6H5;=QAm9g8&5f=j(B $BCx:n8"=jM-(B
 *  
 *  $BAO:n<T!'?"<G=SIW(B
 *
 *  $BK\%W%m%0%i%`$O!JFH!K;:6H5;=QAm9g8&5f=j$N?&0w$G$"$k?"<G=SIW$,AO:n$7!$(B
 *  $B!JFH!K;:6H5;=QAm9g8&5f=j$,Cx:n8"$r=jM-$9$kHkL)>pJs$G$9!%Cx:n8"=jM-(B
 *  $B<T$K$h$k5v2D$J$7$KK\%W%m%0%i%`$r;HMQ!$J#@=!$2~JQ!$Bh;0<T$X3+<($9$k(B
 *  $BEy$N9T0Y$r6X;_$7$^$9!%(B
 *  
 *  $B$3$N%W%m%0%i%`$K$h$C$F@8$8$k$$$+$J$kB;32$KBP$7$F$b!$Cx:n8"=jM-<T$*(B
 *  $B$h$SAO:n<T$O@UG$$rIi$$$^$;$s!#(B
 *
 *  Copyright 2002-2007.
 *  National Institute of Advanced Industrial Science and Technology (AIST)
 *
 *  Creator: Toshio UESHIBA
 *
 *  [AIST Confidential and all rights reserved.]
 *  This program is confidential. Any using, copying, changing or
 *  giving any information concerning with this program to others
 *  without permission by the copyright holder are strictly prohibited.
 *
 *  [No Warranty.]
 *  The copyright holder or the creator are not responsible for any
 *  damages caused by using this program.
 *  
 *  $Id: GaussianCoefficients.cc,v 1.4 2008-09-10 05:10:36 ueshiba Exp $
 */
#include "TU/GaussianConvolver.h"
#include "TU/Minimize.h"

namespace TU
{
/************************************************************************
*  static functions							*
************************************************************************/
static void
coefficients4(float a0, float b0, float omega0, float alpha0,
	      float a1, float b1, float omega1, float alpha1,
	      float c[8])
{
    const float	c0 = cosf(omega0), s0 = sinf(omega0), e0 = expf(-alpha0),
		c1 = cosf(omega1), s1 = sinf(omega1), e1 = expf(-alpha1);
    c[0] =  e0*e1*(e1*(b0*s0 - a0*c0) + e0*(b1*s1 - a1*c1));	// i(n-3)
    c[1] =  a0*e1*e1 + a1*e0*e0
	 -  2.0*e0*e1*((b0*s0 - a0*c0)*c1 + (b1*s1 - a1*c1)*c0);// i(n-2)
    c[2] =  e0*(b0*s0 - (a0 + 2.0*a1)*c0)
	 +  e1*(b1*s1 - (a1 + 2.0*a0)*c1);			// i(n-1)
    c[3] =  a0 + a1;						// i(n)
    c[4] = -e0*e0*e1*e1;					// o(n-4)
    c[5] =  2.0*e0*e1*(e1*c0 + e0*c1);				// o(n-3)
    c[6] = -e0*e0 - e1*e1 - 4.0*e0*e1*c0*c1;			// o(n-2)
    c[7] =  2.0*(e0*c0 + e1*c1);				// o(n-1)
}

/************************************************************************
*  class GaussianCoefficients::Params					*
************************************************************************/
void
GaussianCoefficients::Params::set(double aa, double bb, double tt, double aaa)
{
    a	  = aa;
    b	  = bb;
    theta = tt;
    alpha = aaa;
}

GaussianCoefficients::Params&
GaussianCoefficients::Params::operator -=(const Vector<double>& p)
{
    a	  -= p[0];
    b	  -= p[1];
    theta -= p[2];
    alpha -= p[3];

    return *this;
}

/************************************************************************
*  class GaussianCoefficients::EvenConstraint				*
************************************************************************/
Vector<double>
GaussianCoefficients::EvenConstraint::operator ()(const AT& params) const
{
    Vector<double>	val(1);
    const double	as0 = params[0].alpha/_sigma,
			ts0 = params[0].theta/_sigma,
			as1 = params[1].alpha/_sigma,
			ts1 = params[1].theta/_sigma;
    val[0] = (params[0].a*sinh(as0) + params[0].b* sin(ts0)) *
	     (cosh(as1) - cos(ts1))
	   + (params[1].a*sinh(as1) + params[1].b* sin(ts1)) *
	     (cosh(as0) - cos(ts0));

    return val;
}

Matrix<double>
GaussianCoefficients::EvenConstraint::jacobian(const AT& params) const
{
    const double	c0  = cos(params[0].theta/_sigma),
			s0  = sin(params[0].theta/_sigma),
			ch0 = cosh(params[0].alpha/_sigma),
			sh0 = sinh(params[0].alpha/_sigma),
			c1  = cos(params[1].theta/_sigma),
			s1  = sin(params[1].theta/_sigma),
			ch1 = cosh(params[1].alpha/_sigma),
			sh1 = sinh(params[1].alpha/_sigma);
	
    Matrix<double>	val(1, 8);
    val[0][0] = sh0*(ch1 - c1);
    val[0][1] = s0 *(ch1 - c1);
    val[0][2] = (params[0].b*c0 *(ch1 - c1) +
		 s0 *(params[1].a*sh1 + params[1].b*s1)) / _sigma;
    val[0][3] = (params[0].a*ch0*(ch1 - c1) +
		 sh0*(params[1].a*sh1 + params[1].b*s1)) / _sigma;
    val[0][4] = sh1*(ch0 - c0);
    val[0][5] = s1 *(ch0 - c0);
    val[0][6] = (params[1].b*c1 *(ch0 - c0) +
		 s1 *(params[0].a*sh0 + params[0].b*s0)) / _sigma;
    val[0][7] = (params[1].a*ch1*(ch0 - c0) +
		 sh1*(params[0].a*sh0 + params[0].b*s0)) / _sigma;

    return val;
}

/************************************************************************
*  class GaussianCoefficients::CostFunction				*
************************************************************************/
Vector<double>
GaussianCoefficients::CostFunction::operator ()(const AT& params) const
{
    Vector<double>	val(_ndivisions+1);
    for (int k = 0; k < val.dim(); ++k)
    {
	double	f = 0.0, x = k*_range/_ndivisions;
	for (int i = 0; i < params.dim(); ++i)
	{
	    const Params&	p = params[i];
	    f += (p.a*cos(x*p.theta) + p.b*sin(x*p.theta))*exp(-x*p.alpha);
	}
	val[k] = f - (x*x - 1.0)*exp(-x*x/2.0);
    }
    
    return val;
}
    
Matrix<double>
GaussianCoefficients::CostFunction::jacobian(const AT& params) const
{
    Matrix<double>	val(_ndivisions+1, 4*params.dim());
    for (int k = 0; k < val.nrow(); ++k)
    {
	Vector<double>&	row = val[k];
	double		x = k*_range/_ndivisions;
	
	for (int i = 0; i < params.dim(); ++i)
	{
	    const Params&	p = params[i];
	    const double	c = cos(x*p.theta), s = sin(x*p.theta),
				e = exp(-x*p.alpha);
	    row[4*i]   = c * e;
	    row[4*i+1] = s * e;
	    row[4*i+2] = (-p.a * s + p.b * c) * e;
	    row[4*i+3] = -x * (p.a * c + p.b * s) * e;
	}
    }
    
    return val;
}

void
GaussianCoefficients::CostFunction::update(AT& params,
					   const Vector<double>& dp) const
{
    for (int i = 0; i < params.dim(); ++i)
	params[i] -= dp(4*i, 4);
}

/************************************************************************
*  class GaussianCoefficients						*
************************************************************************/
//! $B$3$N(BGauss$B3K$N=i4|2=$r9T$&(B
/*!
  \param sigma	$B%U%#%k%?%5%$%:$rI=$9@5?t!JBg$-$$$[$I9-$,$j$,Bg$-$$!K(B
  \return	$B$3$N(BGauss$B3K<+?H(B
*/
void
GaussianCoefficients::initialize(float sigma)
{
    coefficients4( 1.80579,   7.42555, 0.676413/sigma, 2.10032/sigma,
		  -0.805838, -1.50785, 1.90174/sigma,  2.15811/sigma, _c0);
    
    coefficients4(-0.628422, -4.68837,  0.666686/sigma, 1.54201/sigma,
		   0.628422,  0.980129, 2.08425/sigma,  1.52152/sigma, _c1);

  /*    coefficients4(-1.27844,   3.24717, 0.752205/sigma, 1.18524/sigma,
	0.278487, -1.54294, 2.21984/sigma,  1.27214/sigma, _c2);*/
    CostFunction::AT	params(CostFunction::D);
    params[0].set(-1.3, 3.6, 0.75, 1.2);
    params[1].set(0.32, -1.7, 2.2, 1.3);
    CostFunction	err(100, 5.0);
    if (sigma < 0.55)
	throw std::runtime_error("GaussianCoefficients::initialize(): sigma must be greater than 0.55");
    EvenConstraint	g(sigma);
    minimizeSquare(err, g, params, 1000, 1.0e-6);
    coefficients4(params[0].a, params[0].b,
		  params[0].theta/sigma, params[0].alpha/sigma,
		  params[1].a, params[1].b,
		  params[1].theta/sigma, params[1].alpha/sigma, _c2);
}

}