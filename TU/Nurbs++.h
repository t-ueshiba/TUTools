/*!
  \file		Nurbs++.h
  \author	Toshio UESHIBA
  \brief	非有理/有理B-spline曲線/曲面に関連するクラスの定義と実装
*/
#ifndef TU_NURBSPP_H
#define TU_NURBSPP_H

#include "TU/Vector++.h"

namespace TU
{
/************************************************************************
*  class BSplineKnots<T>						*
************************************************************************/
//! B-spline曲線または曲面のノットを表すクラス
/*!
  \param T	ノットの値の型
*/
template <class T>
class BSplineKnots
{
  public:
    using element_type	= T;
    using knot_array	= Array<element_type>;
    using knot_array2	= Array2<element_type>;
    
  public:
    BSplineKnots(size_t deg, element_type us, element_type ue)		;
    
    auto	degree()		 const	{return _degree;}
    auto	M()			 const	{return size()-1;}
    auto	L()			 const	{return M()-degree()-degree();}
    size_t	findSpan(element_type u) const	;
    size_t	leftmost(size_t k)	 const	;
    size_t	rightmost(size_t k)	 const	;
    size_t	multiplicity(size_t k)	 const	;

    knot_array	basis(element_type u, size_t& I)		 const	;
    knot_array2	derivatives(element_type u, size_t K, size_t& I) const	;

    size_t	insertKnot(element_type u)	;
    size_t	removeKnot(size_t k)		;
    void	elevateDegree()			{++_degree;}

    auto	size()			const	{return _knots.size();}
    auto	data()			const	{return _knots.data();}
    auto&	operator [](size_t i)		{return _knots[i];}
    const auto&	operator [](size_t i)	const	{return _knots[i];}
    
  private:
    knot_array	_knots;
    size_t	_degree;
};

/*
 *  Create a knot sequence of {us, ..., us, ue, ..., ue},
 *			       ^^^^^^^^^^^  ^^^^^^^^^^^
 *			       deg+1 times  deg+1 times
 */
template <class T>
BSplineKnots<T>::BSplineKnots(size_t deg, element_type us, element_type ue)
    :_knots(deg+1+deg+1), _degree(deg)
{
    for (size_t i = 0; i <= degree(); ++i)
	_knots[i] = us;
    for (size_t i = M() - degree(); i <= M(); ++i)
	_knots[i] = ue;
}

/*
 *  Find span such that u_{span} <= u < u_{span+1} for given 'u'.
 */
template <class T> size_t
BSplineKnots<T>::findSpan(element_type u) const
{
    if (u == _knots[M()-degree()])	// special case
	return M()-degree()-1;

  // binary search
    for (size_t low = degree(), high = M()-degree(); low != high; )
    {
	size_t	mid = (low + high) / 2;

	if (u < _knots[mid])
	    high = mid;
	else if (u >= _knots[mid+1])
	    low = mid;
	else
	    return mid;
    }

    throw std::out_of_range("TU::BSplineKnots<T>::findSpan: given parameter is out of range!");
    
    return 0;
}

/*
 *  Return index of the leftmost knot with same value as k-th knot.
 */
template <class T> size_t
BSplineKnots<T>::leftmost(size_t k) const
{
    while (k > 0 && _knots[k-1] == _knots[k])
	--k;
    return k;
}

/*
 *  Return index of the rightmost knot with same value as k-th knot.
 */
template <class T> size_t
BSplineKnots<T>::rightmost(size_t k) const
{
    while (k+1 <= M() && _knots[k+1] == _knots[k])
	++k;
    return k;
}

/*
 *  Return multiplicity of k-th knot.
 */
template <class T> size_t
BSplineKnots<T>::multiplicity(size_t k) const
{
    return rightmost(k) - leftmost(k) + 1;
}

/*
 *  Compute 'I' such that u_{I} <= u < u_{I} and return an array with
 *  values of basis:
 *    array[i] = N_{I-p+i}(u) where 0 <= i <= degree.
 */
template <class T> typename BSplineKnots<T>::knot_array
BSplineKnots<T>::basis(element_type u, size_t& I) const
{
    I = findSpan(u);
    
    knot_array	Npi(degree()+1);
    knot_array	left(degree()), right(degree());
    Npi[0] = 1.0;
    for (size_t i = 0; i < degree(); ++i)
    {
	left[i]	 = u - _knots[I-i];
	right[i] = _knots[I+i+1] - u;
	element_type  saved = 0.0;
	for (size_t j = 0; j <= i; ++j)
	{
	    const element_type	tmp = Npi[j] / (right[j] + left[i-j]);
	    Npi[j] = saved + right[j]*tmp;
	    saved  = left[i-j]*tmp;
	}
	Npi[i+1] = saved;
    }
    return Npi;
}

/*
 *  Compute 'I' such that u_{I} <= u < u_{I} and return an 2D array with
 *  derivative values of basis:
 *    array[k][i] = "k-th derivative of N_{I-p+i}(u)"
 *	where 0 <= k <= K and 0 <= i <= degree.
 */
template <class T> typename BSplineKnots<T>::knot_array2
BSplineKnots<T>::derivatives(element_type u, size_t K, size_t& I) const
{
    I = findSpan(u);
    
    knot_array2	ndu(degree()+1, degree()+1);
    knot_array	left(degree()), right(degree());
    ndu[0][0] = 1.0;
    for (size_t i = 0; i < degree(); ++i)
    {
	left[i]  = u - _knots[I-i];
	right[i] = _knots[I+i+1] - u;
	element_type	saved = 0.0;
	for (size_t j = 0; j <= i; ++j)
	{
	    ndu[j][i+1] = right[j] + left[i-j];		// upper triangle

	    const element_type	tmp = ndu[i][j] / ndu[j][i+1];
	    ndu[i+1][j] = saved + right[j]*tmp;		// lower triangle
	    saved	= left[i-j]*tmp;
	}
	ndu[i+1][i+1] = saved;				// diagonal elements
    }
    
    knot_array2	N(K+1, degree()+1);
    N[0] = ndu[degree()];				// values of basis
    for (size_t i = 0; i <= degree(); ++i)
    {
	knot_array2	a(2, degree()+1);
	int		previous = 0, current = 1;
	a[previous][0] = 1.0;
	for (size_t k = 1; k <= K; ++k)			// k-th derivative
	{
	    N[k][i] = 0.0;
	    for (size_t j = k - std::min(k, i); j <= std::min(k, degree()-i);
		 ++j)
	    {
		a[current][j] = ((j != k ? a[previous][j]   : 0.0) -
				 (j != 0 ? a[previous][j-1] : 0.0))
			      / ndu[i-k+j][degree()-k+1];
		N[k][i] += a[current][j] * ndu[degree()-k][i-k+j];
	    }
	    std::swap(current, previous);
	}
    }

  // Multiply factors    
    size_t	fact = degree();
    for (size_t k = 1; k <= K; ++k)
    {
	for (size_t i = 0; i <= degree(); ++i)
	    N[k][i] *= fact;
	fact *= (degree() - k);
    }
    
    return N;
}

/*
 *  Insert a knot with value 'u' and return its index of location.
 */
template <class T> size_t
BSplineKnots<T>::insertKnot(element_type u)
{
    size_t	l = findSpan(u) + 1;	// insertion point for the new knot
    knot_array	tmp(_knots);
    _knots.resize(size() + 1);
    for (size_t i = 0; i < l; ++i)	// copy unchanged knots
	_knots[i] = tmp[i];
    _knots[l] = u;			// insert a new knot
    for (size_t i = M(); i > l; --i)	// shift unchanged knots
	_knots[i] = tmp[i-1];
    return rightmost(l);		// index of the new knot
}

/*
 *  Remove k-th knot and return its right-most index.
 */
template <class T> size_t
BSplineKnots<T>::removeKnot(size_t k)
{
    k = rightmost(k);			// index of the knot to be removed
    knot_array	tmp(_knots);
    _knots.resize(size() - 1);
    for (size_t i = 0; i < k; ++i)	// copy unchanged knots
	_knots[i] = tmp[i];
    for (size_t i = M(); i >= k; --i)	// shift unchanged knots
	_knots[i] = tmp[i+1];
    return k;				// index of the new knot
}

/************************************************************************
*  class BSplineCurve<C>						*
************************************************************************/
//! 非有理または有理B-spline曲線を表すクラス
/*!
  \param C	制御点座標の型．d次元空間中の非有理曲線であればd次元ベクトル，
		有理曲線であれば(d+1)次元ベクトル．
*/
template <class C>
class BSplineCurve
{
  public:
    using coord_type	= C;
    using coord_array	= Array<coord_type>;
    using element_type	= typename coord_type::element_type;
    using knots_type	= BSplineKnots<element_type>;
    using knot_array	= typename knots_type::knot_array;
    using knot_array2	= typename knots_type::knot_array2;
    
  public:
    explicit BSplineCurve(size_t degree,
			  element_type us=0, element_type ue=1)	;

    static auto	dim()				{ return coord_type::size(); }

    auto	degree()		const	{ return _knots.degree(); }
    auto	M()			const	{ return _knots.M(); }
    auto	L()			const	{ return _knots.L(); }
    auto	N()			const	{ return _c.size()-1; }
    auto	knot(int i)		const	{ return _knots[i]; }
    auto	multiplicity(size_t k)	const	{ return
						     _knots.multiplicity(k); }
    const auto&	knots()			const	{ return _knots; }

    coord_type	operator ()(element_type u)		const	;
    coord_array	derivatives(element_type u, size_t K)	const	;

    size_t	insertKnot(element_type u)	;
    size_t	removeKnot(size_t k)		;
    void	elevateDegree()			;
    auto	data()			const	{ return _c[0].data(); }
    auto&	operator [](size_t i)		{ return _c[i]; }
    const auto&	operator [](size_t i)	const	{ return _c[i]; }
    auto	operator ==(const BSplineCurve& b) const
		{
		    return _c == b._c;
		}
    auto	operator !=(const BSplineCurve& b) const
		{
		    return _c != b._c;
		}
    auto&	save(std::ostream& out)	const	{ return _c.save(out); }
    auto&	restore(std::istream& in)	{ return _c.restore(in); }

  private:
    coord_array	_c;
    knots_type	_knots;
};

template <class C>
BSplineCurve<C>::BSplineCurve(size_t degree, element_type us, element_type ue)
    :_c(degree + 1), _knots(degree, us, ue)
{
}

/*
 *    Evaluate the coodinate of the curve at 'u'.
 */
template <class C> typename BSplineCurve<C>::coord_type
BSplineCurve<C>::operator ()(element_type u) const
{
    size_t	span;
    knot_array	N = _knots.basis(u, span);
    coord_type	c;
    for (size_t i = 0; i <= degree(); ++i)
	c += N[i] * _c[span-degree()+i];
    return c;
}

/*
 *    Evaluate up to K-th derivatives of the curve at 'u':
 *      array[k] = "k-th derivative of the curve at 'u'" where 0 <= k <= K.
 */
template <class C> typename BSplineCurve<C>::coord_array
BSplineCurve<C>::derivatives(element_type u, size_t K) const
{
    size_t	I;
    knot_array2	dN = _knots.derivatives(u, std::min(K,degree()), I);
    coord_array	ders(K+1);
    for (size_t k = 0; k < dN.nrow(); ++k)
	for (size_t i = 0; i <= degree(); ++i)
	    ders[k] += dN[k][i] * _c[I-degree()+i];
    return ders;
}

/*
 *  size_t BSplineCurve<C>::insertKnot(element_type u)
 *
 *    Insert a knot at 'u', recompute control points and return the index
 *    of the new knot.
 */
template <class C> size_t
BSplineCurve<C>::insertKnot(element_type u)
{
    size_t	l = _knots.insertKnot(u);
    coord_array	tmp(_c);
    _c.resize(_c.size() + 1);
    for (size_t i = 0; i < l-degree(); ++i)
	_c[i] = tmp[i];		// copy unchanged control points
    for (size_t i = l-degree(); i < l; ++i)
    {
      //  Note that we have already inserted a new knot at l. So, old
      //  knot(i+degree()) must be accressed as knot(i+degree()+1).

	element_type	alpha = (u - knot(i))
			      / (knot(i+degree()+1) - knot(i));

	_c[i] = (element_type(1) - alpha) * tmp[i-1] + alpha * tmp[i];
    }
    for (size_t i = N(); i >= l; --i)
	_c[i] = tmp[i-1];		// copy unchanged control points

    return l;
}

/*
 *  size_t BSplineCurve<C>::removeKnot(size_t k)
 *
 *    Remove k-th knot, recompute control points and return the index of
 *    the removed knot.
 */
template <class C> size_t
BSplineCurve<C>::removeKnot(size_t k)
{
    size_t		s = multiplicity(k);
    element_type	u = knot(k);
    k = _knots.removeKnot(k);
    coord_array		tmp(_c);
    _c.resize(_c.size() - 1);
    size_t	i, j;
    for (i = 0; i < k - degree(); ++i)
	_c[i] = tmp[i];			// copy unchanged control points
    for (j = N(); j >= k - s; --j)
	_c[j] = tmp[j+1];		// copy unchanged control points
    for (i = k - degree(), j = k - s; i < j - 1; ++i, --j)
    {
	element_type	alpha_i = (u - knot(i))
				/ (knot(i+degree()) - knot(i)),
			alpha_j = (u - knot(j))
				/ (knot(j+degree()) - knot(j));
	_c[i]   = (tmp[i] - (element_type(1) - alpha_i) * _c[i-1]) / alpha_i;
	_c[j-1] = (tmp[j] - alpha_j * _c[j]) / (element_type(1) - alpha_j);
    }
    if (i == j - 1)
    {
	element_type	alpha_i = (u - knot(i))
				/ (knot(i+degree()) - knot(i)),
			alpha_j = (u - knot(j))
				/ (knot(j+degree()) - knot(j));
	_c[i] = ((tmp[i] - (element_type(1) - alpha_i) * _c[i-1]) / alpha_i + 
		 (tmp[j] - alpha_j * _c[j]) / (element_type(1) - alpha_j))
	      / 2.0;
    }
    
    return k;
}

/*
 *  Elevate degree of the curve by one.
 */
template <class C> void
BSplineCurve<C>::elevateDegree()
{
  // Convert to Bezier segments.
    Array<size_t>	mul(L());
    size_t		nsegments = 1;
    for (size_t k = degree() + 1; k < M() - degree(); k += degree())
    {
      // Elevate multiplicity of each internal knot to degree().
	mul[nsegments-1] = multiplicity(k);
	for (size_t n = mul[nsegments-1]; n < degree(); ++n)
	    insertKnot(knot(k));
	++nsegments;
    }
  // Save control points of Bezier segments.    
    coord_array		tmp(_c);

  // Set knots and allocate area for control points.
    for (size_t k = 0; k <= M(); )  // Elevate multiplicity of each knot by one.
	k = _knots.insertKnot(knot(k)) + 1;
    _knots.elevateDegree();
    _c.resize(_c.size() + nsegments);

  // Elevate degree of each Bezier segment.
    for (size_t n = 0; n < nsegments; ++n)
    {
	_c[n*degree()] = tmp[n*(degree()-1)];
	for (size_t i = 1; i < degree(); ++i)
	{
	    element_type	alpha = element_type(i)
				      / element_type(degree());
	    
	    _c[n*degree()+i]
		= alpha * tmp[n*(degree()-1)+i-1]
		+ (element_type(1) - alpha) * tmp[n*(degree()-1)+i];
	}
    }
    _c[nsegments*degree()] = tmp[nsegments*(degree()-1)];

  // Remove redundant internal knots.
    for (size_t k = degree() + 1, n = 0; k < M() - degree(); k += mul[n]+1, ++n)
	for (size_t r = degree(); --r > mul[n]; )
	    removeKnot(k);
}

using BSplineCurve2f		= BSplineCurve<Vector<float, 2> >;
using RationalBSplineCurve2f	= BSplineCurve<Vector<float, 3> >;
using BSplineCurve3f		= BSplineCurve<Vector<float, 3> >;
using RationalBSplineCurve3f	= BSplineCurve<Vector<float, 4> >;
using BSplineCurve2d		= BSplineCurve<Vector<double, 2> >;
using RationalBSplineCurve2d	= BSplineCurve<Vector<double, 3> >;
using BSplineCurve3d		= BSplineCurve<Vector<double, 3> >;
using RationalBSplineCurve3d	= BSplineCurve<Vector<double, 4> >;
    
/************************************************************************
*  class BSplineSurface<C>						*
************************************************************************/
//! 非有理または有理B-spline曲面を表すクラス
/*!
  \param C	制御点の型．d次元空間中の非有理曲面であればd次元ベクトル，
		有理曲面であれば(d+1)次元ベクトル．
*/
template <class C>
class BSplineSurface
{
  public:
    using coord_type	= C;
    using coord_array	= Array<coord_type>;
    using coord_array2	= Array2<coord_type>;
    using element_type	= typename coord_type::element_type;
    using knots_type	= BSplineKnots<element_type>;
    using knot_array	= typename knots_type::knot_array;
    using knot_array2	= typename knots_type::knot_array2;
    
    BSplineSurface(size_t uDegree, size_t vDegree,
		   element_type us=0, element_type ue=1,
		   element_type vs=0, element_type ve=1)	;

    static auto	dim()			{ return coord_type::size(); }

    auto	uDegree()		const	{ return _uKnots.degree(); }
    auto	uM()			const	{ return _uKnots.M(); }
    auto	uL()			const	{ return _uKnots.L(); }
    auto	uN()			const	{ return ncol()-1; }
    auto	vDegree()		const	{ return _vKnots.degree(); }
    auto	vM()			const	{ return _vKnots.M(); }
    auto	vL()			const	{ return _vKnots.L(); }
    auto	vN()			const	{ return nrow()-1; }
    auto	uKnot(int i)		const	{ return _uKnots[i]; }
    auto	vKnot(int j)		const	{ return _vKnots[j]; }
    auto	uMultiplicity(int k)	const	{ return
						     _uKnots.multiplicity(k); }
    auto	vMultiplicity(int l)	const	{ return
						     _vKnots.multiplicity(l); }
    const auto&	uKnots()		const	{ return _uKnots; }
    const auto&	vKnots()		const	{ return _vKnots; }

    coord_type	operator ()(element_type u,
			    element_type v)	const	;
    coord_array2
		derivatives(element_type u, element_type v,
			    size_t D)		const	;

    size_t	uInsertKnot(element_type u)		;
    size_t	vInsertKnot(element_type v)		;
    size_t	uRemoveKnot(size_t k)			;
    size_t	vRemoveKnot(size_t l)			;
    void	uElevateDegree()			;
    void	vElevateDegree()			;

    auto	data()			const	{ return _c[0][0].data(); }
    auto	operator [](size_t i)		{ return _c[i]; }
    auto	operator [](size_t i)	const	{ return _c[i]; }
    auto	ncol()			const	{ return _c.ncol(); }
    auto	nrow()			const	{ return _c.nrow(); }
    auto	operator ==(const BSplineSurface& b) const
		{
		    return _c == b._c;
		}
    auto	operator !=(const BSplineSurface& b) const
		{
		    return _c != b._c;
		}
    auto&	save(std::ostream& out)	const	{ return _c.save(out); }
    auto&	restore(std::istream& in)	{ return _c.restore(in); }
    
    friend std::istream&
		operator >>(std::istream& in, BSplineSurface& b)
		{
		    return in >> b._c;
		}
    friend std::ostream&
		operator <<(std::ostream& out, const BSplineSurface& b)
		{
		    return out << b._c;
		}

  private:
    coord_array2	_c;
    knots_type		_uKnots, _vKnots;
};

template <class C>
BSplineSurface<C>::BSplineSurface(size_t uDeg, size_t vDeg,
				  element_type us, element_type ue,
				  element_type vs, element_type ve)
    :_c(vDeg + 1, uDeg + 1), _uKnots(uDeg, us, ue), _vKnots(vDeg, vs, ve)
{
}

/*
 *    Evaluate the coodinate of the surface at (u, v).
 */
template <class C> C
BSplineSurface<C>::operator ()(element_type u, element_type v) const
{
    size_t	uSpan, vSpan;
    knot_array	Nu = _uKnots.basis(u, uSpan);
    knot_array	Nv = _vKnots.basis(v, vSpan);
    coord_type	c;
    for (size_t j = 0; j <= vDegree(); ++j)
    {
	coord_type	tmp;
	for (size_t i = 0; i <= uDegree(); ++i)
	    tmp += Nu[i] * _c[vSpan-vDegree()+j][uSpan-uDegree()+i];
	c += Nv[j] * tmp;
    }
    return c;
}

/*
 *    Evaluate up derivatives of the surface at (u, v):
 *      array[l][k] = "derivative of order k w.r.t u and order l w.r.t v"
 *        where 0 <= k <= D and 0 <= l <= D.
 */
template <class C> typename BSplineSurface<C>::coord_array2
BSplineSurface<C>::derivatives(element_type u, element_type v, size_t D) const
{
    size_t		I, J;
    knot_array2		udN = _uKnots.derivatives(u, min(D,uDegree()), I),
			vdN = _vKnots.derivatives(v, min(D,vDegree()), J);
    coord_array2	ders(D+1, D+1);
    for (size_t k = 0; k < udN.nrow(); ++k)		// derivatives w.r.t u
    {
	Array<coord_type>	tmp(vDegree()+1);
	for (size_t j = 0; j <= vDegree(); ++j)
	    for (size_t i = 0; i <= uDegree(); ++i)
		tmp[j] += udN[k][i] * _c[J-vDegree()+j][I-uDegree()+i];
	for (size_t l = 0;
	     l < std::min(vdN.nrow(), D-k); ++l)	// derivatives w.r.t v
	    for (size_t j = 0; j <= vDegree(); ++j)
		ders[l][k] += vdN[l][j] * tmp[j];
    }
    return ders;
}

/*
 *  int BSplineSurface<C>::uInsertKnot(element_type u)
 *
 *    Insert a knot in u-direction at 'u', recompute control points and return
 *    the index of the new knot.
 */
template <class C> size_t
BSplineSurface<C>::uInsertKnot(element_type u)
{
    size_t		l = _uKnots.insertKnot(u);
    coord_array2	tmp(_c);
    _c.resize(nrow(), ncol()+1);
    Array<element_type>	alpha(uDegree());
    for (size_t i = l-uDegree(); i < l; ++i)
	alpha[i-l+uDegree()] =
	    (u - uKnot(i)) / (uKnot(i+uDegree()+1) - uKnot(i));
    for (size_t j = 0; j <= vN(); ++j)
    {
	for (size_t i = 0; i < l-uDegree(); ++i)
	    _c[j][i] = tmp[j][i];
	for (size_t i = l-uDegree(); i < l; ++i)
	    _c[j][i] = (element_type(1) - alpha[i-l+uDegree()]) * tmp[j][i-1]
		     + alpha[i-l+uDegree()]			* tmp[j][i];
	for (size_t i = uN(); i >= l; --i)
	    _c[j][i] = tmp[j][i-1];
    }

    return l;
}

/*
 *  size_t BSplineSurface<C>::vInsertKnot(element_type v)
 *
 *    Insert a knot in v-direction at 'v', recompute control points and return
 *    the index of the new knot.
 */
template <class C> size_t
BSplineSurface<C>::vInsertKnot(element_type v)
{
    size_t		l = _vKnots.insertKnot(v);
    coord_array2	tmp(_c);
    _c.resize(nrow()+1, ncol());
    Array<element_type>	alpha(vDegree());
    for (size_t j = l-vDegree(); j < l; ++j)
	alpha[j-l+vDegree()] =
	    (v - vKnot(j)) / (vKnot(j+vDegree()+1) - vKnot(j));
    for (size_t i = 0; i <= uN(); ++i)
    {
	for (size_t j = 0; j < l-vDegree(); ++j)
	    _c[j][i] = tmp[j][i];
	for (size_t j = l-vDegree(); j < l; ++j)
	    _c[j][i] = (element_type(1) - alpha[j-l+vDegree()]) * tmp[j-1][i]
		     + alpha[j-l+vDegree()]			* tmp[j]  [i];
	for (size_t j = vN(); j >= l; --j)
	    _c[j][i] = tmp[j-1][i];
    }

    return l;
}

/*
 *  Remove k-th knot in u-derection, recompute control points and return 
 *  the index of the removed knot.
 */
template <class C> size_t
BSplineSurface<C>::uRemoveKnot(size_t k)
{
    size_t		s = uMultiplicity(k);
    element_type	u = uKnot(k);
    k = _uKnots.removeKnot(k);
    coord_array2	tmp(_c);
    _c.resize(nrow(), ncol()-1);
    for (size_t j = 0; j <= vN(); ++j)
    {
	size_t	is, ie;
	for (is = 0; is < k - uDegree(); ++is)
	    _c[j][is] = tmp[j][is];    // copy unchanged control points
	for (ie = uN(); ie >= k - s; --ie)
	    _c[j][ie] = tmp[j][ie+1];  // copy unchanged control points
	for (is = k - uDegree(), ie = k - s; is < ie - 1; ++is, --ie)
	{
	    element_type	alpha_s = (u - uKnot(is))
					/ (uKnot(is+uDegree()) - uKnot(is)),
				alpha_e = (u - uKnot(ie))
					/ (uKnot(ie+uDegree()) - uKnot(ie));
	    _c[j][is]   = (tmp[j][is] - (element_type(1) - alpha_s)*_c[j][is-1])
			/ alpha_s;
	    _c[j][ie-1] = (tmp[j][ie] - alpha_e * _c[j][ie])
			/ (element_type(1) - alpha_e);
	}
	if (is == ie - 1)
	{
	    element_type	alpha_s = (u - uKnot(is))
					/ (uKnot(is+uDegree()) - uKnot(is)),
				alpha_e = (u - uKnot(ie))
					/ (uKnot(ie+uDegree()) - uKnot(ie));
	    _c[j][is] = ((tmp[j][is] -
			  (element_type(1) - alpha_s)*_c[j][is-1]) / alpha_s +
			 (tmp[j][ie] - alpha_e * _c[j][ie]) /
			 (element_type(1) - alpha_e)) / 2.0;
	}
    }
    
    return k;
}

/*
 *  Remove l-th knot in v-derection, recompute control points and return 
 *  the index of the removed knot.
 */
template <class C> size_t
BSplineSurface<C>::vRemoveKnot(size_t l)
{
    size_t		s = vMultiplicity(l);
    element_type	v = vKnot(l);
    l = _vKnots.removeKnot(l);
    coord_array2	tmp(_c);
    _c.resize(nrow()-1, ncol());
    for (size_t i = 0; i <= uN(); ++i)
    {
	size_t	js, je;
	for (js = 0; js < l - vDegree(); ++js)
	    _c[js][i] = tmp[js][i];    // copy unchanged control points
	for (je = vN(); je >= l - s; --je)
	    _c[je][i] = tmp[je+1][i];  // copy unchanged control points
	for (js = l - vDegree(), je = l - s; js < je - 1; ++js, --je)
	{
	    element_type	alpha_s = (v - vKnot(js))
					/ (vKnot(js+vDegree()) - vKnot(js)),
				alpha_e	= (v - vKnot(je))
					/ (vKnot(je+vDegree()) - vKnot(je));
	    _c[js][i]   = (tmp[js][i] - (element_type(1) - alpha_s)*_c[js-1][i])
			/ alpha_s;
	    _c[je-1][i] = (tmp[je][i] - alpha_e * _c[je][i])
			/ (element_type(1) - alpha_e);
	}
	if (js == je - 1)
	{
	    element_type	alpha_s	= (v - vKnot(js))
					/ (vKnot(js+vDegree()) - vKnot(js)),
				alpha_e = (v - vKnot(je))
					/ (vKnot(je+vDegree()) - vKnot(je));
	    _c[js][i] = ((tmp[js][i] -
			  (element_type(1) - alpha_s)*_c[js-1][i]) / alpha_s +
			 (tmp[je][i] - alpha_e * _c[je][i]) /
			 (element_type(1) - alpha_e)) / 2.0;
	}
    }
    
    return l;
}

/*
 *  Elevate degree of the surface by one in u-direction.
 */
template <class C> void
BSplineSurface<C>::uElevateDegree()
{
  // Convert to Bezier segments.
    Array<size_t>	mul(uL());
    size_t		nsegments = 1;
    for (size_t k = uDegree() + 1; k < uM() - uDegree(); k += uDegree())
    {
      // Elevate multiplicity of each internal knot to uDegree().
	mul[nsegments-1] = uMultiplicity(k);
	for (size_t n = mul[nsegments-1]; n < uDegree(); ++n)
	    uInsertKnot(uKnot(k));
	++nsegments;
    }
    coord_array2	tmp(_c);	// Save Bezier control points.

  // Set knots and allocate area for control points.
    for (size_t k = 0; k <= uM(); )
	k = _uKnots.insertKnot(uKnot(k)) + 1;
    _uKnots.elevateDegree();
    _c.resize(nrow(), ncol() + nsegments);
    
  // Elevate degree of each Bezier segment.
    for (size_t j = 0; j <= vN(); ++j)
    {
	for (size_t n = 0; n < nsegments; ++n)
	{
	    _c[j][n*uDegree()] = tmp[j][n*(uDegree()-1)];
	    for (size_t i = 1; i < uDegree(); ++i)
	    {
		element_type	alpha = element_type(i)
				      / element_type(uDegree());
	    
		_c[j][n*uDegree()+i]
		    = alpha			* tmp[j][n*(uDegree()-1)+i-1]
		    + (element_type(1) - alpha) * tmp[j][n*(uDegree()-1)+i];
	    }
	    _c[j][nsegments*uDegree()] = tmp[j][nsegments*(uDegree()-1)];
	}
    }

  // Remove redundant internal knots.
    for (size_t k = uDegree() + 1, j = 0;
	 k < uM() - uDegree(); k += mul[j]+1, ++j)
	for (size_t r = uDegree(); --r > mul[j]; )
	    uRemoveKnot(k);
}

/*
 *  void BSplineSurface<C>::vElevateDegree()
 *
 *    Elevate degree of the surface by one in v-direction.
 */
template <class C> void
BSplineSurface<C>::vElevateDegree()
{
  // Convert to Bezier segments.
    Array<size_t>	mul(vL());
    size_t		nsegments = 1;
    for (size_t l = vDegree() + 1; l < vM() - vDegree(); l += vDegree())
    {
      // Elevate multiplicity of each internal knot to vDegree().
	mul[nsegments-1] = vMultiplicity(l);
	for (size_t n = mul[nsegments-1]; n < vDegree(); ++n)
	    vInsertKnot(vKnot(l));
	++nsegments;
    }
    coord_array2	tmp(_c);	// Save Bezier control points.

  // Set knots and allocate area for control points.
    for (size_t l = 0; l <= vM(); )
	l = _vKnots.insertKnot(vKnot(l)) + 1;
    _vKnots.elevateDegree();
    _c.resize(nrow() + nsegments, ncol());
    
  // Elevate degree of each Bezier segment.
    for (size_t i = 0; i <= uN(); ++i)
    {
	for (size_t n = 0; n < nsegments; ++n)
	{
	    _c[n*vDegree()][i] = tmp[n*(vDegree()-1)][i];
	    for (size_t j = 1; j < vDegree(); ++j)
	    {
		element_type	alpha = element_type(j)
				      / element_type(vDegree());
	    
		_c[n*vDegree()+j][i]
		    = alpha			* tmp[n*(vDegree()-1)+j-1][i]
		    + (element_type(1) - alpha) * tmp[n*(vDegree()-1)+j][i];
	    }
	    _c[nsegments*vDegree()][i] = tmp[nsegments*(vDegree()-1)][i];
	}
    }

  // Remove redundant internal knots.
    for (size_t l = vDegree() + 1, i = 0;
	 l < vM() - vDegree(); l += mul[i]+1, ++i)
	for (size_t r = vDegree(); --r > mul[i]; )
	    vRemoveKnot(l);
}

using BSplineSurface3f		= BSplineSurface<Vector3f>;
using RationalBSplineSurface3f	= BSplineSurface<Vector4f>;
using BSplineSurface3d		= BSplineSurface<Vector3d>;
using RationalBSplineSurface3d	= BSplineSurface<Vector4d>;

}
#endif	// !TU_NURBSPP_H
