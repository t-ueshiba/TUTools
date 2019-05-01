/*!
 *  \file	Quaternion.h
 *  \author	Toshio UESHIBA
 *  \brief	四元数を表すクラスの定義と実装
 */
#ifndef TU_QUATERNION_H
#define TU_QUATERNION_H

#include <boost/operators.hpp>
#include "TU/Vector++.h"

namespace TU
{
/************************************************************************
*  class Quaternion<T>							*
************************************************************************/
//! 四元数を表すクラス
/*!
    \param T	要素の型
*/
template <class T>
class Quaternion : boost::field_operators<Quaternion<T> >
{
  public:
    using scalar_type	= T;		    //!< スカラー部の型
    using vector_type	= Vector<T, 3>;	    //!< ベクトル部の型
    using vector4_type	= Vector<T, 4>;	    //!< ベクトルとして取り出される型
    using rotation_type	= Matrix<T, 3, 3>;  //!< 回転行列として取り出される型
    
  public:
			Quaternion(const scalar_type& s=scalar_type(),
				   const scalar_type& x=scalar_type(),
				   const scalar_type& y=scalar_type(),
				   const scalar_type& z=scalar_type())
			    :_s(s), _v({x, y, z})			{}
			Quaternion(const scalar_type& s,
				   const vector_type& v)
			    :_s(s), _v(v)				{}
    template <class T_>	Quaternion(const Vector<T_, 4>& a)
			    :_s(a[0]), _v({a[1], a[2], a[3]})		{}
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
			Quaternion(const E_& expr)
			{
			    using E = element_t<E_>;

			    const auto&	R = evaluate(expr);

			    if (size<0>(R) != 3 || size<1>(R) != 3)
				throw std::invalid_argument("TU::Quaternion::Quaternion(): input matrix must be 3x3!");

			    _s = E(0.5) * std::sqrt(trace(R) + E(1));
			    if (_s + E(1) == E(1))	// _s << 1 ?
			    {
				Array<E, 3>	evalues;
				_v = eigen(R, evalues, false)[0];
			    }
			    else
			    {
				const auto	k = E(0.25) / _s;
				const auto	S = R - transpose(R);
				_v = {k*S[2][1], k*S[0][2], k*S[1][0]};
			    }
			}
    
    const scalar_type&	scalar()			const	{ return _s; }
    void		scalar(const scalar_type& s)		{ _s = s; }
    const vector_type&	vector()			const	{ return _v; }
    void		vector(const vector_type& v)		{ _v = v; }

    Quaternion&		normalize()
			{
			    return *this /= std::sqrt(square(*this));
			}
    
    bool		operator ==(const Quaternion& x)
			{
			    return _s == x._s && _v[0] == x._v[0]
					      && _v[1] == x._v[1]
					      && _v[2] == x._v[2];
			}
    
    const Quaternion&	operator +() const
			{
			    return *this;
			}

    Quaternion		operator -() const
			{
			    return {-_s, -_v};
			}

    Quaternion&		operator +=(const Quaternion& x)
			{
			    _s += x._s;
			    _v += x._v;
			    return *this;
			}

    Quaternion&		operator -=(const Quaternion& x)
			{
			    _s -= x._s;
			    _v -= x._v;
			    return *this;
			}

    Quaternion&		operator *=(const scalar_type& c)
			{
			    _s *= c;
			    _v *= c;
			    return *this;
			}

    Quaternion&		operator /=(const scalar_type& c)
			{
			    _s /= c;
			    _v /= c;
			    return *this;
			}

    Quaternion&		operator *=(const Quaternion& x)
			{
			    const auto	s = _s;
			    _s = s * x._s - _v * x._v;
			    _v = evaluate(s * x._v + _v * x._s + (_v ^ x._v));
			    return *this;
			}

    Quaternion&		operator /=(const Quaternion& x)
			{
			    *this *= inverse(x);
			    return *this;
			}

			operator vector4_type() const
			{
			    return {_s, _v[0], _v[1], _v[2]};
			}

    rotation_type	R() const
			{
			    vector_type		q = 2 * _v;
			    rotation_type	Q = q % _v;
			    const auto		c = square(_s) - square(_v);
			    Q[0][0] += c;
			    Q[1][1] += c;
			    Q[2][2] += c;

			    q *= _s;
			    Q[0][1] -= q[2];
			    Q[0][2] += q[1];
			    Q[1][0] += q[2];
			    Q[1][2] -= q[0];
			    Q[2][0] -= q[1];
			    Q[2][1] += q[0];

			    return Q;
			}

    rotation_type	Rt() const
			{
			    return conj(*this).R();
			}
    
    vector_type		rpy() const
			{
			    T	roll, pitch, yaw;
			    rotation_angle(Rt(), roll, pitch, yaw);

			    return {roll, pitch, yaw};
			}

    std::ostream&	put(std::ostream& out) const
			{
			    return _v.put(out) << ' ' << _s;
			}

    std::istream&	get(std::istream& in)
			{
			    return _v.get(in) >> _s;
			}
    
    friend std::istream&
			operator >>(std::istream& in, Quaternion& x)
			{
			    return in >> x._s >> x._v;
			}

  private:
    scalar_type	_s;	//!< scalar part
    vector_type	_v;	//!< vector part
};

template <class T> inline T
square(const Quaternion<T>& x)
{
    return square(x.scalar()) + square(x.vector());
}
    
template <class T> inline Quaternion<T>
conj(const Quaternion<T>& x)
{
    return {x.scalar(), -x.vector()};
}
    
template <class T> inline Quaternion<T>
inverse(const Quaternion<T>& x)
{
    auto	q = conj(x);
    q /= square(q);
    return q;
}
    
template <class T> inline Quaternion<T>
operator *(const Quaternion<T>& x, const T& c)
{
    auto	y = x;
    y *= c;
    return y;
}
    
template <class T> inline Quaternion<T>
operator *(const T& c, const Quaternion<T>& x)
{
    return x * c;
}
    
template <class T> inline Quaternion<T>
operator /(const Quaternion<T>& x, const T& c)
{
    auto	y = x;
    y /= c;
    return y;
}

template <class T> inline std::ostream&
operator <<(std::ostream& out, const Quaternion<T>& x)
{
    return out << x.scalar() << x.vector();
}
    
}
#endif	// !TU_QUATERNION_H
