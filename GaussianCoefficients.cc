/*!
  \file		GaussianCoefficients.cc
  \author	Toshio UESHIBA
  \brief	クラス TU::GaussianConvolver の実装
*/
#include "TU/GaussianConvolver.h"
#include "TU/Minimize.h"

namespace TU
{
/************************************************************************
*  static functions							*
************************************************************************/
template <class T> static void
coefficients4(double a0, double b0, double omega0, double alpha0,
	      double a1, double b1, double omega1, double alpha1, T c[8])
{
    using namespace	std;
    
    const double	c0 = cos(omega0), s0 = sin(omega0), e0 = exp(-alpha0),
			c1 = cos(omega1), s1 = sin(omega1), e1 = exp(-alpha1);
    c[0] =  e0*e1*(e1*(b0*s0 - a0*c0) + e0*(b1*s1 - a1*c1));	// i(n-3)
    c[1] =  a0*e1*e1 + a1*e0*e0
	 -  2*e0*e1*((b0*s0 - a0*c0)*c1 + (b1*s1 - a1*c1)*c0);	// i(n-2)
    c[2] =  e0*(b0*s0 - (a0 + 2*a1)*c0)
	 +  e1*(b1*s1 - (a1 + 2*a0)*c1);			// i(n-1)
    c[3] =  a0 + a1;						// i(n)
    c[4] = -e0*e0*e1*e1;					// o(n-4)
    c[5] =  2*e0*e1*(e1*c0 + e0*c1);				// o(n-3)
    c[6] = -e0*e0 - e1*e1 - 4*e0*e1*c0*c1;			// o(n-2)
    c[7] =  2*(e0*c0 + e1*c1);					// o(n-1)
}

/************************************************************************
*  class GaussianCoefficients<T>::Params				*
************************************************************************/
template <class T> void
GaussianCoefficients<T>::Params::set(element_type aa, element_type bb,
				     element_type tt, element_type aaa)
{
    a	  = aa;
    b	  = bb;
    theta = tt;
    alpha = aaa;
}

template <class T> typename GaussianCoefficients<T>::Params&
GaussianCoefficients<T>::Params::operator -=(const vector_type& p)
{
    a	  -= p[0];
    b	  -= p[1];
    theta -= p[2];
    alpha -= p[3];

    return *this;
}

/************************************************************************
*  class GaussianCoefficients<T>::EvenConstraint			*
************************************************************************/
template <class T> typename GaussianCoefficients<T>::vector_type
GaussianCoefficients<T>::EvenConstraint::operator ()(const argument_type& params) const
{
    vector_type		val(1u);
    const element_type	as0 = params[0].alpha/_sigma,
			ts0 = params[0].theta/_sigma,
			as1 = params[1].alpha/_sigma,
			ts1 = params[1].theta/_sigma;
    val[0] = (params[0].a*sinh(as0) + params[0].b* sin(ts0)) *
	     (cosh(as1) - cos(ts1))
	   + (params[1].a*sinh(as1) + params[1].b* sin(ts1)) *
	     (cosh(as0) - cos(ts0));

    return val;
}

template <class T> typename GaussianCoefficients<T>::matrix_type
GaussianCoefficients<T>::EvenConstraint::derivative(const argument_type& params) const
{
    using namespace	std;
    
    const element_type	c0  = cos (params[0].theta/_sigma),
			s0  = sin (params[0].theta/_sigma),
			ch0 = cosh(params[0].alpha/_sigma),
			sh0 = sinh(params[0].alpha/_sigma),
			c1  = cos (params[1].theta/_sigma),
			s1  = sin (params[1].theta/_sigma),
			ch1 = cosh(params[1].alpha/_sigma),
			sh1 = sinh(params[1].alpha/_sigma);
	
    matrix_type	val(1, 8);
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
*  class GaussianCoefficients<T>::CostFunction				*
************************************************************************/
template <class T> typename GaussianCoefficients<T>::vector_type
GaussianCoefficients<T>::CostFunction
		       ::operator ()(const argument_type& params) const
{
    vector_type	val(_ndivisions+1);
    for (size_t k = 0; k < val.size(); ++k)
    {
	element_type	f = 0.0, x = k*_range/_ndivisions;
	for (size_t i = 0; i < params.size(); ++i)
	{
	    const auto&	p = params[i];
	    f += (p.a*cos(x*p.theta) + p.b*sin(x*p.theta))*exp(-x*p.alpha);
	}
	val[k] = f - (x*x - 1.0)*exp(-x*x/2.0);
    }
    
    return val;
}
    
template <class T> typename GaussianCoefficients<T>::matrix_type
GaussianCoefficients<T>::CostFunction
		       ::derivative(const argument_type& params) const
{
    matrix_type	val(_ndivisions+1, 4*params.size());
    for (size_t k = 0; k < val.nrow(); ++k)
    {
	auto	row = val[k];
	auto	x = k*_range/_ndivisions;
	
	for (size_t i = 0; i < params.size(); ++i)
	{
	    const auto&		p = params[i];
	    const element_type	c = cos(x*p.theta), s = sin(x*p.theta),
				e = exp(-x*p.alpha);
	    row[4*i]   = c * e;
	    row[4*i+1] = s * e;
	    row[4*i+2] = (-p.a * s + p.b * c) * e;
	    row[4*i+3] = -x * (p.a * c + p.b * s) * e;
	}
    }
    
    return val;
}

template <class T> void
GaussianCoefficients<T>::CostFunction::update(argument_type& params,
					      const vector_type& dp) const
{
    for (size_t i = 0; i < params.size(); ++i)
	params[i] -= dp(4*i, 4);
}

/************************************************************************
*  class GaussianCoefficients<T>					*
************************************************************************/
//! このGauss核の初期化を行う
/*!
  \param sigma	フィルタサイズを表す正数（大きいほど広がりが大きい）
  \return	このGauss核自身
*/
template <class T> void
GaussianCoefficients<T>::initialize(T sigma)
{
    _sigma = sigma;
    
    coefficients4( 1.80579,   7.42555, 0.676413/_sigma, 2.10032/_sigma,
		  -0.805838, -1.50785, 1.90174 /_sigma, 2.15811/_sigma, _c0);
    
    coefficients4(-0.628422, -4.68837,  0.666686/_sigma, 1.54201/_sigma,
		   0.628422,  0.980129, 2.08425 /_sigma, 1.52152/_sigma, _c1);

  /*    coefficients4(-1.27844,   3.24717, 0.752205/_sigma, 1.18524/_sigma,
	0.278487, -1.54294, 2.21984/_sigma,  1.27214/_sigma, _c2);*/
    typename CostFunction::argument_type	params(CostFunction::D);
    params[0].set(-1.3, 3.6, 0.75, 1.2);
    params[1].set(0.32, -1.7, 2.2, 1.3);
    CostFunction	err(100, 5.0);
    if (_sigma < 0.55)
	throw std::runtime_error("GaussianCoefficients::initialize(): _sigma must be greater than 0.55");
    EvenConstraint	g(_sigma);
    minimizeSquare(err, g, params, 1000, 1.0e-6);
    coefficients4(params[0].a, params[0].b,
		  params[0].theta/_sigma, params[0].alpha/_sigma,
		  params[1].a, params[1].b,
		  params[1].theta/_sigma, params[1].alpha/_sigma, _c2);
}

template class GaussianCoefficients<float>;
template class GaussianCoefficients<double>;
    
}
