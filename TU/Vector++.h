/*!
  \file		Vector++.h
  \author	Toshio UESHIBA
  \brief	ベクトルと行列の演算およびそれに関連するクラスの定義と実装
*/
#ifndef TU_VECTORPP_H
#define TU_VECTORPP_H

#include "TU/Array++.h"

namespace TU
{
/************************************************************************
*  generic algorithms for vectors					*
************************************************************************/
//! 3次元ベクトルから3x3反対称行列を生成する．
/*!
  \return	生成された反対称行列，すなわち
  \f[
    \TUskew{u}{} \equiv
    \TUbeginarray{ccc}
      & -u_2 & u_1 \\ u_2 & & -u_0 \\ -u_1 & u_0 &
    \TUendarray
  \f]
  \throw std::invalid_argument	3次元ベクトルでない場合に送出
*/
template <class E, std::enable_if_t<rank<E>() == 1>* = nullptr> inline auto
skew(const E& expr)
{
    using result_t = Array2<element_t<E>, 3, 3>;
    
    assert(size<0>(expr) == 3);

    const auto&	a = evaluate(expr);
    return result_t({{0, -a[2], a[1]}, {a[2], 0, -a[0]}, {-a[1], a[0], 0}});
}

//! 非同次座標を表すベクトルに対し，値1を持つ成分を最後に付加した同次座標ベクトルを返す．
/*!
  \return	同次化されたベクトル
*/
template <class E, std::enable_if_t<rank<E>() == 1>* = nullptr> inline auto
homogeneous(const E& expr)
{
    constexpr size_t	N = size0<E>();
    using element_type	= element_t<E>;
    using result_type	= std::conditional_t<N == 0,
					     Array<element_type, 0>,
					     Array<element_type, N+1> >;

    const auto	n = std::size(expr);
    result_type	r(n + 1);
    slice(r, 0, n) = expr;
    r[n]	   = 1;
    return r;
}

//! 同次座標を表すベクトルに対し，各成分を最後の成分で割った非同次座標ベクトルを返す．
/*!
  \return	非同次化されたベクトル
*/
template <class E, std::enable_if_t<rank<E>() == 1>* = nullptr> inline auto
inhomogeneous(const E& expr)
{
    constexpr size_t	N = size0<E>();
    using element_type	= element_t<E>;
    using result_type	= std::conditional_t<N == 0,
					     Array<element_type, 0>,
					     Array<element_type, N-1> >;

    const auto	n = std::size(expr) - 1;
    return result_type(slice(expr, 0, n) / *(std::cbegin(expr) + n));
}

/************************************************************************
*  generic algorithms for matrices					*
************************************************************************/
//! 全て同一の対角成分値を持つ正方行列を生成する．
/*!
  \param c	対角成分の値
  \return	対角行列，すなわち\f$\TUvec{A}{} \leftarrow \diag(c,\ldots,c)\f$
*/
template <size_t N=0, class T> auto
diag(T c, size_t n=N)
{
    Array2<T, N, N>	a(n, n);
    for (size_t i = 0; i != a.size(); ++i)
	a[i][i] = c;
    return a;
}

//! 正方行列のtraceを返す．
/*!
  \return			trace, すなわち\f$\trace\TUvec{A}{}\f$
  \throw std::invalid_argument	正方行列でない場合に送出
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> auto
trace(const E& expr)
{
    assert(size<0>(expr) == size<1>(expr));

    element_t<E>	val = 0;
    for (size_t i = 0; i < std::size(expr); ++i)
	val += expr[i][i];
    return val;
}

//! 正値対称行列のCholesky分解（上半三角行列）を返す．
/*!
  計算においては，もとの行列の上半部分しか使わない
  \param A	正値対称行列
  \return	\f$\TUvec{A}{} = \TUvec{L}{}\TUtvec{L}{}\f$なる
		\f$\TUtvec{L}{}\f$（上半三角行列）
  \throw std::invalid_argument	正方行列でない場合に送出
  \throw std::runtime_error	正値でない場合に送出
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> auto
cholesky(const E& A)
{
    if (size<0>(A) != size<1>(A))
        throw std::invalid_argument("TU::cholesky(): not square matrix!!");

    Array2<element_t<E>, size0<E>(), size0<E>()>	Lt(A);
    for (size_t i = 0; i < Lt.nrow(); ++i)
    {
	const auto	Li = Lt[i];
	auto		d = Li[i];
	if (d <= 0)
	    throw std::runtime_error("TU::cholesky(): not positive definite matrix!!");
	for (size_t j = 0; j < i; ++j)
	    Li[j] = 0;
	Li[i] = d = std::sqrt(d);
	for (size_t j = i + 1; j < Li.size(); ++j)
	    Li[j] /= d;
	for (size_t j = i + 1; j < Lt.nrow(); ++j)
	    for (size_t k = j; k < Li.size(); ++k)
		Lt[j][k] -= (Li[j] * Li[k]);
    }
    
    return Lt;
}

//! 正方行列の下半三角部分を上半三角部分にコピーして対称化する．
/*!
    \return	この行列
*/
template <class T, size_t N> Array2<T, N, N>&
symmetrize(Array2<T, N, N>& a)
{
    assert(a.nrow() == a.ncol());
    
    for (size_t i = 0; i < a.nrow(); ++i)
	for (size_t j = 0; j < i; ++j)
	    a[j][i] = a[i][j];
    return a;
}

//! この行列の下半三角部分の符号を反転し，上半三角部分にコピーして反対称化する．
/*!
    \return	この行列
*/
template <class T, size_t N> Array2<T, N, N>&
antisymmetrize(Array2<T, N, N>& a)
{
    assert(a.nrow() == a.ncol());
    
    for (size_t i = 0; i < a.nrow(); ++i)
    {
	a[i][i] = 0.0;
	for (size_t j = 0; j < i; ++j)
	    a[j][i] = -a[i][j];
    }
    return a;
}

/************************************************************************
*  class LUDecomposition<T, N>						*
************************************************************************/
//! 正方行列のLU分解を表すクラス
template <class T, size_t N=0>
class LUDecomposition
{
  public:
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    LUDecomposition(const E_& expr)			;

    template <class E_> std::enable_if_t<rank<E_>() == 1, E_&>
			substitute(E_&& expr)	const	;

  //! もとの正方行列の行列式を返す．
  /*!
    \return	もとの正方行列の行列式
  */
    auto		det()			const	{ return _det; }
    auto		size()			const	{ return _A.size(); }
	
  private:
    Array2<T, N, N>	_A;
    Array<size_t>	_indices;
    T			_det;
};

//! 与えられた正方行列のLU分解を生成する．
/*!
  \param m			LU分解する正方行列
  \throw std::invalid_argument	mが正方行列でない場合に送出
*/
template <class T, size_t N>
template <class E_, std::enable_if_t<rank<E_>() == 2>*>
LUDecomposition<T, N>::LUDecomposition(const E_& expr)
    :_A(expr), _indices(size()), _det(1)
{
    if (_A.nrow() != _A.ncol())
        throw std::invalid_argument("TU::LUDecomposition<T>::LUDecomposition: not square matrix!!");

  // Initialize column index for explicit pivotting
    for (size_t j = 0; j < _indices.size(); ++j)
	_indices[j] = j;

  // Find maximum abs. value in each col. for implicit pivotting
    Array<T, N>	scale(size());
    for (size_t j = 0; j < size(); ++j)
    {
	T	max = 0;
	for (const auto a : _A)
	{
	    const T	tmp = std::fabs(a[j]);
	    if (tmp > max)
		max = tmp;
	}
	scale[j] = (max != 0 ? T(1) / max : T(1));
    }

    for (size_t i = 0; i < size(); ++i)
    {
      // Left part (j < i)
	for (size_t j = 0; j < i; ++j)
	{
	    T&	sum = _A[i][j];
	    for (size_t k = 0; k < j; ++k)
		sum -= _A[i][k] * _A[k][j];
	}

      // Diagonal and right part (i <= j)
	size_t	jmax = i;
	T	max  = 0;
	for (size_t j = i; j < size(); ++j)
	{
	    T&	sum = _A[i][j];
	    for (size_t k = 0; k < i; ++k)
		sum -= _A[i][k] * _A[k][j];
	    const T	tmp = std::fabs(sum) * scale[j];
	    if (tmp >= max)
	    {
		max  = tmp;
		jmax = j;
	    }
	}
	if (jmax != i)			// pivotting required ?
	{
	  // Swap i-th and jmax-th column	    
	    for (size_t k = 0; k < size(); ++k)
		std::swap(_A[k][i], _A[k][jmax]);
	    std::swap(_indices[i], _indices[jmax]);	// swap column indices
	    std::swap(scale[i], scale[jmax]);	// swap colum-wise scale factor
	    _det = -_det;
	}

	_det *= _A[i][i];

	if (_A[i][i] == 0)	// singular matrix ?
	    break;

	for (size_t j = i + 1; j < size(); ++j)
	    _A[i][j] /= _A[i][i];
    }
}

//! もとの正方行列を係数行列とした連立1次方程式を解く．
/*!
  \param b			もとの正方行列\f$\TUvec{M}{}\f$と同じ次
				元を持つベクトル．\f$\TUtvec{b}{} =
				\TUtvec{x}{}\TUvec{M}{}\f$の解に変換さ
				れる．
  \throw std::invalid_argument	ベクトルbの次元がもとの正方行列の次元に一致
				しない場合に送出
  \throw std::runtime_error	もとの正方行列が正則でない場合に送出
*/
template <class T, size_t N>
template <class E_> std::enable_if_t<rank<E_>() == 1, E_&>
LUDecomposition<T, N>::substitute(E_&& b) const
{
    if (std::size(b) != size())
	throw std::invalid_argument("TU::LUDecomposition<T>::substitute: Dimension of given vector is not equal to mine!!");

  // tmpはbと異なる実体でなければならないので，
  //     auto tmp = evaluate(b);
  // ではダメ．
    Array<T, N>	tmp = b;
    for (size_t j = 0; j < size(); ++j)
	b[j] = tmp[_indices[j]];

    for (size_t j = 0; j < size(); ++j)		// forward substitution
	for (size_t i = 0; i < j; ++i)
	    b[j] -= b[i] * _A[i][j];
    for (size_t j = size(); j-- > 0; )		// backward substitution
    {
	for (size_t i = size(); --i > j; )
	    b[j] -= b[i] * _A[i][j];
	if (_A[j][j] == 0)			// singular matrix ?
	    throw std::runtime_error("TU::LUDecomposition<T>::substitute: singular matrix !!");
	b[j] /= _A[j][j];
    }

    return b;
}

/************************************************************************
*  generic algorithms using LUDecomposition<T, N>			*
************************************************************************/
//! 行列の行列式を返す．
/*!
  \param A	正方行列
  \return	行列式，すなわち\f$\det\TUvec{A}{}\f$
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> inline auto
det(const E& A)
{
    return LUDecomposition<element_t<E>, size0<E>()>(A).det();
}

//! 行列の小行列式を返す．
/*!
  \param A	正方行列
  \param p	元の行列から取り除く行を指定するindex
  \param q	元の行列から取り除く列を指定するindex
  \return	小行列式，すなわち\f$\det\TUvec{A}{pq}\f$
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> auto
det(const E& A, size_t p, size_t q)
{
    const auto&			A_ = evaluate(A);
    Array2<element_t<E> >	D(A_.nrow()-1, A_.ncol()-1);
    for (size_t i = 0; i < p; ++i)
    {
	for (size_t j = 0; j < q; ++j)
	    D[i][j] = A_[i][j];
	for (size_t j = q; j < D.ncol(); ++j)
	    D[i][j] = A_[i][j+1];
    }
    for (size_t i = p; i < D.nrow(); ++i)
    {
	for (size_t j = 0; j < q; ++j)
	    D[i][j] = A_[i+1][j];
	for (size_t j = q; j < D.ncol(); ++j)
	    D[i][j] = A_[i+1][j+1];
    }
    return det(D);
}

//! 行列の余因子行列を返す．
/*!
  \param A	正方行列
  \return	余因子行列，すなわち
		\f$\TUtilde{A}{} = (\det\TUvec{A}{})\TUinv{A}{}\f$
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> auto
adjoint(const E& A)
{
    constexpr size_t		N = size0<E>();
    Array2<element_t<E>, N, N>	val(size<0>(A), size<1>(A));
    for (size_t i = 0; i < val.nrow(); ++i)
	for (size_t j = 0; j < val.ncol(); ++j)
	    val[i][j] = ((i + j) % 2 ? -det(A, j, i) : det(A, j, i));
    return val;
}

//! 連立1次方程式を解く．
/*!
  \param A	正則な正方行列
  \return	\f$\TUtvec{b}{} = \TUtvec{x}{}\TUvec{A}{}\f$
		の解ベクトル，すなわち \f$\TUtvec{b}{}\TUinv{A}{}\f$
*/
template <class E, class F>
inline std::enable_if_t<rank<E>() == 2 && rank<F>() == 1, F&>
solve(const E& A, F&& b)
{
    return LUDecomposition<element_t<E>, size0<E>()>(A)
	  .substitute(std::forward<F>(b));
}

//! 連立1次方程式を解く．
/*!
  \param A	正則な正方行列
  \param B	行列
  \return	\f$\TUvec{B}{} = \TUvec{X}{}\TUvec{A}{}\f$
		の解を納めた行列，すなわち	\f$\TUvec{B}{}\TUinv{A}{}\f$
*/
template <class E, class F>
inline std::enable_if_t<rank<E>() == 2 && rank<F>() == 2, F&>
solve(const E& A, F&& B)
{
    using	std::begin;
    
    LUDecomposition<element_t<E>, size0<E>()>	lu(A);
    for_each<size0<F>()>([&lu](auto&& b)
			 { lu.substitute(std::forward<decltype(b)>(b)); },
			 std::size(B), begin(B));
    return B;
}
    
//! 行列の逆行列を返す．
/*!
  \param A	正則な正方行列
  \return	逆行列，すなわち\f$\TUinv{A}{}\f$
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> inline auto
inverse(const E& A)
{
    using	element_type = element_t<E>;
    
    constexpr size_t		N = size0<E>();
    Array2<element_type, N, N>	B = diag<N>(element_type(1), std::size(A));
    return solve(A, B);
}
    
/************************************************************************
*  class Householder<T>							*
************************************************************************/
template <class T>	class QRDecomposition;
template <class T>	class TriDiagonal;
template <class T>	class BiDiagonal;

//! Householder変換を表すクラス
template <class T>
class Householder : public Array2<T>
{
  private:
    using super	= Array2<T>;
    
  public:
    using	typename super::element_type;
    
  private:
    Householder(size_t dd, size_t d)
	:super(dd, dd), _d(d), _sigma(super::nrow())		{}
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    Householder(const E_& a, size_t d)				;

    using		super::size;
    
    void		apply_from_left(Array2<T>& a, size_t m)	;
    void		apply_from_right(Array2<T>& a, size_t m);
    void		apply_from_both(Array2<T>& a, size_t m)	;
    void		make_transformation()			;
    const Array<T>&	sigma()				const	{return _sigma;}
    Array<T>&		sigma()					{return _sigma;}
    bool		sigma_is_zero(size_t m, T comp)	const	;

  private:
    const size_t	_d;		// deviation from diagonal element
    Array<T>		_sigma;

    friend class	QRDecomposition<T>;
    friend class	TriDiagonal<T>;
    friend class	BiDiagonal<T>;
};

template <class T>
template <class E_, typename std::enable_if_t<rank<E_>() == 2>*>
Householder<T>::Householder(const E_& a, size_t d)
    :super(a), _d(d), _sigma(size())
{
    if (super::nrow() != super::ncol())
	throw std::invalid_argument("TU::Householder<T>::Householder: Given matrix must be square !!");
}

template <class T> void
Householder<T>::apply_from_left(Array2<T>& a, size_t m)
{
    if (a.nrow() < size())
	throw std::invalid_argument("TU::Householder<T>::apply_from_left: # of rows of given matrix is smaller than my dimension !!");
    
    T	scale = 0.0;
    for (size_t i = m+_d; i < size(); ++i)
	scale += std::fabs(a[i][m]);
	
    if (scale != 0.0)
    {
	T	h = 0.0;
	for (size_t i = m+_d; i < size(); ++i)
	{
	    a[i][m] /= scale;
	    h += a[i][m] * a[i][m];
	}

	const T	s = (a[m+_d][m] > 0.0 ? std::sqrt(h) : -std::sqrt(h));
	h	     += s * a[m+_d][m];			// H = u^2 / 2
	a[m+_d][m]   += s;				// m-th col <== u
	    
	for (size_t j = m+1; j < a.ncol(); ++j)
	{
	    T	p = 0.0;
	    for (size_t i = m+_d; i < size(); ++i)
		p += a[i][m] * a[i][j];
	    p /= h;					// p[j] (p' = u'A / H)
	    for (size_t i = m+_d; i < size(); ++i)
		a[i][j] -= a[i][m] * p;			// A = A - u*p'
	    a[m+_d][j] = -a[m+_d][j];
	}
	    
	for (size_t i = m+_d; i < size(); ++i)
	    (*this)[m][i] = scale * a[i][m];		// copy u
	_sigma[m+_d] = scale * s;
    }
}

template <class T> void
Householder<T>::apply_from_right(Array2<T>& a, size_t m)
{
    if (a.ncol() < size())
	throw std::invalid_argument("Householder<T>::apply_from_right: # of column of given matrix is smaller than my dimension !!");
    
    T	scale = 0.0;
    for (size_t j = m+_d; j < size(); ++j)
	scale += std::fabs(a[m][j]);
	
    if (scale != 0.0)
    {
	T	h = 0.0;
	for (size_t j = m+_d; j < size(); ++j)
	{
	    a[m][j] /= scale;
	    h += a[m][j] * a[m][j];
	}

	const T	s = (a[m][m+_d] > 0.0 ? std::sqrt(h) : -std::sqrt(h));
	h	     += s * a[m][m+_d];			// H = u^2 / 2
	a[m][m+_d]   += s;				// m-th row <== u

	for (size_t i = m+1; i < a.nrow(); ++i)
	{
	    T	p = 0.0;
	    for (size_t j = m+_d; j < size(); ++j)
		p += a[i][j] * a[m][j];
	    p /= h;					// p[i] (p = Au / H)
	    for (size_t j = m+_d; j < size(); ++j)
		a[i][j] -= p * a[m][j];			// A = A - p*u'
	    a[i][m+_d] = -a[i][m+_d];
	}
	    
	for (size_t j = m+_d; j < size(); ++j)
	    (*this)[m][j] = scale * a[m][j];		// copy u
	_sigma[m+_d] = scale * s;
    }
}

template <class T> void
Householder<T>::apply_from_both(Array2<T>& a, size_t m)
{
    auto	u = slice(a[m], m+_d, a.ncol()-m-_d);
    T		scale = 0.0;
    for (size_t j = 0; j < u.size(); ++j)
	scale += std::fabs(u[j]);
	
    if (scale != 0.0)
    {
	u /= scale;

	T	h = u * u;
	const T	s = (u[0] > 0.0 ? std::sqrt(h) : -std::sqrt(h));
	h	+= s * u[0];			// H = u^2 / 2
	u[0]	+= s;				// m-th row <== u

	auto	A = a(m+_d, a.nrow()-m-_d, m+_d, a.ncol()-m-_d);
	auto	p = _sigma(m+_d, size()-m-_d);
	for (size_t i = 0; i < A.size(); ++i)
	    p[i] = (A[i] * u) / h;			// p = Au / H

	const T	k = (u * p) / (h + h);		// K = u*p / 2H
	for (size_t i = 0; i < A.size(); ++i)
	{				// m-th col of 'a' is used as 'q'
	    a[m+_d+i][m] = p[i] - k * u[i];		// q = p - Ku
	    for (size_t j = 0; j <= i; ++j)		// A = A - uq' - qu'
		A[j][i] = (A[i][j] -= (u[i]*a[m+_d+j][m] + a[m+_d+i][m]*u[j]));
	}
	for (size_t j = 1; j < A.size(); ++j)
	    A[j][0] = A[0][j] = -A[0][j];

	for (size_t j = m+_d; j < a.ncol(); ++j)
	    (*this)[m][j] = scale * a[m][j];		// copy u
	_sigma[m+_d] = scale * s;
    }
}

template <class T> void
Householder<T>::make_transformation()
{
    for (size_t m = size(); m-- > 0; )
    {
	for (size_t i = m+1; i < size(); ++i)
	    (*this)[i][m] = 0.0;

	if (_sigma[m] != 0.0)
	{
	    for (size_t i = m+1; i < size(); ++i)
	    {
		T	g = 0.0;
		for (size_t j = m+1; j < size(); ++j)
		    g += (*this)[i][j] * (*this)[m-_d][j];
		g /= (_sigma[m] * (*this)[m-_d][m]);	// g[i] (g = Uu / H)
		for (size_t j = m; j < size(); ++j)
		    (*this)[i][j] -= g * (*this)[m-_d][j];	// U = U - gu'
	    }
	    for (size_t j = m; j < size(); ++j)
		(*this)[m][j] = (*this)[m-_d][j] / _sigma[m];
	    (*this)[m][m] -= 1.0;
	}
	else
	{
	    for (size_t j = m+1; j < size(); ++j)
		(*this)[m][j] = 0.0;
	    (*this)[m][m] = 1.0;
	}
    }
}

template <class T> bool
Householder<T>::sigma_is_zero(size_t m, T comp) const
{
    return (T(std::fabs(_sigma[m])) + comp == comp);
}

/************************************************************************
*  class QRDecomposition<T>						*
************************************************************************/
//! 一般行列のQR分解を表すクラス
/*!
  与えられた行列\f$\TUvec{A}{} \in \TUspace{R}{m\times n}\f$に対して
  \f$\TUvec{A}{} = \TUtvec{R}{}\TUtvec{Q}{}\f$なる下半三角行列
  \f$\TUtvec{R}{} \in \TUspace{R}{m\times n}\f$と回転行列
  \f$\TUtvec{Q}{} \in \TUspace{R}{n\times n}\f$を求める
  （\f$\TUvec{A}{}\f$の各行を\f$\TUtvec{Q}{}\f$の行の線型結合で表現する）．
 */
template <class T>
class QRDecomposition
{
  public:
    using element_type	= T;
    
  public:
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    QRDecomposition(const E_& A)				;

  //! QR分解の下半三角行列を返す．
  /*!
    \return	下半三角行列\f$\TUtvec{R}{}\f$
  */
    const Array2<T>&	Rt()			const	{ return _Rt; }

  //! QR分解の回転行列を返す．
  /*!
    \return	回転行列\f$\TUtvec{Q}{}\f$
  */
    const Array2<T>&	Qt()			const	{ return _Qt; }
    
  private:
    Array2<T>		_Rt;
    Householder<T>	_Qt;			// rotation matrix
};

//! 与えられた一般行列のQR分解を生成する．
/*!
 \param m	QR分解する一般行列
*/
template <class T> template <class E_, std::enable_if_t<rank<E_>() == 2>*>
QRDecomposition<T>::QRDecomposition(const E_& A)
    :_Rt(A), _Qt(_Rt.ncol(), 0)
{
    const size_t	n = std::min(_Rt.nrow(), _Rt.ncol());
    for (size_t j = 0; j < n; ++j)
	_Qt.apply_from_right(_Rt, j);
    _Qt.make_transformation();
    for (size_t i = 0; i < n; ++i)
    {
	auto	r = _Rt[i];
	r[i] = _Qt.sigma()[i];
	for (size_t j = i + 1; j < _Rt.ncol(); ++j)
	    r[j] = 0.0;
    }
}

/************************************************************************
*  class Rotation<T>							*
************************************************************************/
//! 2次元超平面内での回転を表すクラス
/*!
  具体的には
  \f[
    \TUvec{R}{}(p, q; \theta) \equiv
    \begin{array}{r@{}l}
      & \begin{array}{ccccccccccc}
        & & \makebox[4.0em]{} & p & & & \makebox[3.8em]{} & q & & &
      \end{array} \\
      \begin{array}{l}
        \\ \\ \\ \raisebox{1.5ex}{$p$} \\ \\ \\ \\ \raisebox{1.5ex}{$q$} \\ \\ \\
      \end{array} &
      \TUbeginarray{ccccccccccc}
	1 \\
	& \ddots \\
	& & 1 \\
	& & & \cos\theta & & & & -\sin\theta \\
	& & & & 1 \\
	& & & & & \ddots \\
	& & & & & & 1 \\
	& & & \sin\theta & & & & \cos\theta \\
	& & & & & & & & 1\\
	& & & & & & & & & \ddots \\
	& & & & & & & & & & 1
      \TUendarray
    \end{array}
  \f]
  なる回転行列で表される．
*/
template <class T>
class Rotation
{
  public:
    using element_type	= T;	//!< 成分の型
    
  public:
  //! 2次元超平面内での回転を生成する
  /*!
    \param p	p軸を指定するindex
    \param q	q軸を指定するindex
    \param x	回転角を生成する際のx値
    \param y	回転角を生成する際のy値
		\f[
		  \cos\theta = \frac{x}{\sqrt{x^2+y^2}},{\hskip 1em}
		  \sin\theta = \frac{y}{\sqrt{x^2+y^2}}
		\f]
  */
		Rotation(size_t p, size_t q, element_type x, element_type y)
		    :_p(p), _q(q), _l(1), _c(1), _s(0)
		{
		    const element_type	absx = std::abs(x), absy = std::abs(y);
		    _l = (absx > absy ?
			  absx * std::sqrt(1 + (absy*absy)/(absx*absx)) :
			  absy * std::sqrt(1 + (absx*absx)/(absy*absy)));
		    if (_l != 0)
		    {
			_c = x / _l;
			_s = y / _l;
		    }
		}

  //! 2次元超平面内での回転を生成する
  /*!
    \param p		p軸を指定するindex
    \param q		q軸を指定するindex
    \param theta	回転角
  */
		Rotation(size_t p, size_t q, element_type theta)
		    :_p(p), _q(q),
		     _l(1), _c(std::cos(theta)), _s(std::sin(theta))
		{
		}

  //! p軸を返す．
  /*!
    \return	p軸のindex
  */
    size_t	p()				const	{return _p;}

  //! q軸を返す．
  /*!
    \return	q軸のindex
  */
    size_t	q()				const	{return _q;}

  //! 回転角生成ベクトルの長さを返す．
  /*!
    \return	回転角生成ベクトル(x, y)に対して\f$\sqrt{x^2 + y^2}\f$
  */
    T		length()			const	{return _l;}

  //! 回転角のcos値を返す．
  /*!
    \return	回転角のcos値
  */
    T		cos()				const	{return _c;}

  //! 回転角のsin値を返す．
  /*!
    \return	回転角のsin値
  */
    T		sin()				const	{return _s;}
    
  //! この回転行列の転置を与えられた行列の左から掛ける．
  /*!
    \return	与えられた行列，すなわち
		\f$\TUvec{A}{}\leftarrow\TUtvec{R}{}\TUvec{A}{}\f$
  */
    template <class E_>
    std::enable_if_t<rank<E_>() == 2, E_&>
		apply_from_left(E_&& A) const
		{
		    for (size_t j = 0; j < size<1>(A); ++j)
		    {
			const auto	tmp = A[_p][j];
	
			A[_p][j] =  _c*tmp + _s*A[_q][j];
			A[_q][j] = -_s*tmp + _c*A[_q][j];
		    }
		    return A;
		}

  //! この回転行列を与えられた行列の右から掛ける．
  /*!
    \return	与えられた行列，すなわち
		\f$\TUvec{A}{}\leftarrow\TUvec{A}{}\TUvec{R}{}\f$
  */
    template <class E_>
    std::enable_if_t<rank<E_>() == 2, E_&>
		apply_from_right(E_&& A) const
		{
		    for (auto&& a : A)
		    {
			const auto	tmp = a[_p];
	
			a[_p] =  tmp*_c + a[_q]*_s;
			a[_q] = -tmp*_s + a[_q]*_c;
		    }
		    return A;
		}
    
  private:
    const size_t	_p, _q;		// rotation axis
    element_type	_l;		// length of (x, y)
    element_type	_c, _s;		// cos & sin
};

/************************************************************************
*  class TriDiagonal<T>							*
************************************************************************/
//! 対称行列の3重対角化を表すクラス
/*!
  与えられた対称行列\f$\TUvec{A}{} \in \TUspace{R}{d\times d}\f$に対し
  て\f$\TUtvec{U}{}\TUvec{A}{}\TUvec{U}{}\f$が3重対角行列となるような回
  転行列\f$\TUtvec{U}{} \in \TUspace{R}{d\times d}\f$を求める．
 */
template <class T>
class TriDiagonal
{
  public:
    using element_type	= T;	//!< 成分の型
    
  public:
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    TriDiagonal(const E_& a)				;

  //! 3重対角化される対称行列の次元(= 行数 = 列数)を返す．
  /*!
    \return	対称行列の次元
  */
    size_t		size()			const	{return _Ut.nrow();}

  //! 3重対角化を行う回転行列を返す．
  /*!
    \return	回転行列
  */
    const Array2<T>&	Ut()			const	{return _Ut;}

  //! 3重対角行列の対角成分を返す．
  /*!
    \return	対角成分
  */
    const Array<T>&	diagonal()		const	{return _diagonal;}

  //! 3重対角行列の非対角成分を返す．
  /*!
    \return	非対角成分
  */
    const Array<T>&	off_diagonal()		const	{return _Ut.sigma();}

    void		diagonalize(bool abs=true)	;

    template <class E_>
    std::enable_if_t<rank<E_>() == 1, Array2<T> >
			move(E_&& evals)
			{
			    evals = std::move(_diagonal);
			    return std::move(static_cast<Array2<T>&>(_Ut));
			}
    
  private:
    enum		{NITER_MAX = 30};

    bool		off_diagonal_is_zero(size_t n)		const	;
    void		initialize_rotation(size_t m, size_t n,
					    T& x, T& y)		const	;
    
    Householder<T>	_Ut;
    Array<T>		_diagonal;
    Array<T>&		_off_diagonal;
};

//! 与えられた対称行列を3重対角化する．
/*!
  \param a			3重対角化する対称行列
  \throw std::invalid_argument	aが正方行列でない場合に送出
*/
template <class T> template <class E_, std::enable_if_t<rank<E_>() == 2>*>
TriDiagonal<T>::TriDiagonal(const E_& a)
    :_Ut(a, 1), _diagonal(_Ut.nrow()), _off_diagonal(_Ut.sigma())
{
    if (_Ut.nrow() != _Ut.ncol())
        throw std::invalid_argument("TU::TriDiagonal<T>::TriDiagonal: not square matrix!!");

    for (size_t m = 0; m < size(); ++m)
    {
	_Ut.apply_from_both(_Ut, m);
	_diagonal[m] = _Ut[m][m];
    }

    _Ut.make_transformation();
}

//! 3重対角行列を対角化する（固有値，固有ベクトルの計算）．
/*!
  対角成分は固有値となり，\f$\TUtvec{U}{}\f$の各行は固有ベクトルを与える．
  \throw std::runtime_error	指定した繰り返し回数を越えた場合に送出
  \param abs	固有値をその絶対値の大きい順に並べるのであればtrue,
		その値の大きい順に並べるのであればfalse
*/ 
template <class T> void
TriDiagonal<T>::diagonalize(bool abs)
{
    using namespace	std;
    
    for (size_t n = size(); n-- > 0; )
    {
	int	niter = 0;
	
#ifdef TU_DEBUG
	cerr << "******** n = " << n << " ********" << endl;
#endif
	while (!off_diagonal_is_zero(n))
	{					// n > 0 here
	    if (niter++ > NITER_MAX)
		throw runtime_error("TU::TriDiagonal::diagonalize(): Number of iteration exceeded maximum value!!");

	  /* Find first m (< n) whose off-diagonal element is 0 */
	    size_t	m = n;
	    while (!off_diagonal_is_zero(--m))	// 0 <= m < n < size() here
	    {
	    }

	  /* Set x and y which determine initial(i = m+1) plane rotation */
	    T	x, y;
	    initialize_rotation(m, n, x, y);
	  /* Apply rotation P(i-1, i) for each i (i = m+1, n+2, ... , n) */
	    for (size_t i = m; ++i <= n; )
	    {
		Rotation<T>	rot(i-1, i, x, y);
		
		rot.apply_from_left(_Ut);

		if (i > m+1)
		    _off_diagonal[i-1] = rot.length();
		const T w = _diagonal[i] - _diagonal[i-1];
		const T d = rot.sin()*(rot.sin()*w
			       + 2.0*rot.cos()*_off_diagonal[i]);
		_diagonal[i-1]	 += d;
		_diagonal[i]	 -= d;
		_off_diagonal[i] += rot.sin()*(rot.cos()*w
				  - 2.0*rot.sin()*_off_diagonal[i]);
		if (i < n)
		{
		    x = _off_diagonal[i];
		    y = rot.sin()*_off_diagonal[i+1];
		    _off_diagonal[i+1] *= rot.cos();
		}
	    }
#ifdef TU_DEBUG
	    cerr << "  niter = " << niter << ": " << off_diagonal();
#endif	    
	}
    }

  // Sort eigen values and eigen vectors.
    if (abs)
    {
	for (size_t m = 0; m < size(); ++m)
	    for (size_t n = m+1; n < size(); ++n)
		if (std::fabs(_diagonal[n]) >
		    std::fabs(_diagonal[m]))			// abs. values
		{
		    swap(_diagonal[m], _diagonal[n]);
		    for (size_t j = 0; j < size(); ++j)
		    {
			const T	tmp = _Ut[m][j];
			_Ut[m][j] = _Ut[n][j];
			_Ut[n][j] = -tmp;
		    }
		}
    }
    else
    {
	for (size_t m = 0; m < size(); ++m)
	    for (size_t n = m+1; n < size(); ++n)
		if (_diagonal[n] > _diagonal[m])		// raw values
		{
		    swap(_diagonal[m], _diagonal[n]);
		    for (size_t j = 0; j < size(); ++j)
		    {
			const T	tmp = _Ut[m][j];
			_Ut[m][j] = _Ut[n][j];
			_Ut[n][j] = -tmp;
		    }
		}
    }
}

template <class T> bool
TriDiagonal<T>::off_diagonal_is_zero(size_t n) const
{
    return (n == 0 || _Ut.sigma_is_zero(n, std::fabs(_diagonal[n-1]) +
					   std::fabs(_diagonal[n])));
}

template <class T> void
TriDiagonal<T>::initialize_rotation(size_t m, size_t n, T& x, T& y) const
{
    const T	g    = (_diagonal[n] - _diagonal[n-1]) / (2.0*_off_diagonal[n]),
		absg = std::fabs(g),
		gg1  = (absg > 1.0 ?
			absg * std::sqrt(1.0 + (1.0/absg)*(1.0/absg)) :
			std::sqrt(1.0 + absg*absg)),
		t    = (g > 0.0 ? g + gg1 : g - gg1);
    x = _diagonal[m] - _diagonal[n] - _off_diagonal[n]/t;
  //x = _diagonal[m];					// without shifting
    y = _off_diagonal[m+1];
}

/************************************************************************
*  generic algorithms using TriDiagonal<T>				*
************************************************************************/
//! 対称行列の固有値と固有ベクトルを返す．
/*!
  \param A	対称行列 
  \param evals	固有値が返される
  \param abs	固有値を絶対値の大きい順に並べるならtrue, 値の大きい順に
		並べるならfalse
  \return	各行が固有ベクトルから成る回転行列，すなわち
		\f[
		  \TUvec{A}{}\TUvec{U}{} =
		  \TUvec{U}{}\diag(\lambda_0,\ldots,\lambda_{n-1}),
		  {\hskip 1em}\mbox{where}{\hskip 0.5em}
		  \TUtvec{U}{}\TUvec{U}{} = \TUvec{I}{n},~\det\TUvec{U}{} = 1
		\f]
		なる\f$\TUtvec{U}{}\f$
*/
template <class E, class F,
	  std::enable_if_t<rank<E>() == 2 && rank<F>() == 1>* = nullptr> auto
eigen(const E& A, F&& evals, bool abs=true)
{
    TriDiagonal<element_t<E> >	tri(A);
    tri.diagonalize(abs);
    return tri.move(evals);
}

//! 対称行列の一般固有値と一般固有ベクトルを返す．
/*!
  \param A	対称行列 
  \param BB	Aと同一サイズの正値対称行列
  \param evals	一般固有値が返される
  \param abs	一般固有値を絶対値の大きい順に並べるならtrue, 値の大きい順に
		並べるならfalse
  \return	各行が一般固有ベクトルから成る正則行列
		（ただし直交行列ではない），すなわち
		\f[
		  \TUvec{A}{}\TUvec{U}{} =
		  \TUvec{B}{}\TUvec{U}{}\diag(\lambda_0,\ldots,\lambda_{n-1}),
		  {\hskip 1em}\mbox{where}{\hskip 0.5em}
		  \TUtvec{U}{}\TUvec{B}{}\TUvec{U}{} = \TUvec{I}{n}
		\f]
		なる\f$\TUtvec{U}{}\f$
*/
template <class E, class F, class G,
	  std::enable_if_t<rank<E>() == 2 &&
			   rank<F>() == 2 &&
			   rank<G>() == 1>* = nullptr> auto
geigen(const E& A, const F& BB, G&& evals, bool abs=true)
{
    const auto	Ltinv = inverse(BB.cholesky());
    const auto	Linv  = transpose(Ltinv);
    return evaluate(eigen(Linv * A * Ltinv, std::forward<G>(evals), abs) * Linv);
}

/************************************************************************
*  class BiDiagonal<T>							*
************************************************************************/
//! 一般行列の2重対角化を表すクラス
/*!
  与えられた一般行列\f$\TUvec{A}{} \in \TUspace{R}{m\times n}\f$に対し
  て\f$\TUtvec{V}{}\TUvec{A}{}\TUvec{U}{}\f$が2重対角行列となるような2
  つの回転行列\f$\TUtvec{U}{} \in \TUspace{R}{n\times n}\f$,
  \f$\TUtvec{V}{} \in \TUspace{R}{m\times m}\f$を求める．\f$m \le n\f$
  の場合は下半三角な2重対角行列に，\f$m > n\f$の場合は上半三角な2重対角
  行列になる．
 */
template <class T>
class BiDiagonal
{
  public:
    using element_type = T;	//!< 成分の型
    
  public:
    template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr>
    BiDiagonal(const E& a)			;

  //! 2重対角化される行列の行数を返す．
  /*!
    \return	行列の行数
  */
    size_t		nrow()		const	{return _Vt.nrow();}

  //! 2重対角化される行列の列数を返す．
  /*!
    \return	行列の列数
  */
    size_t		ncol()		const	{return _Ut.nrow();}

  //! 2重対角化を行うために右から掛ける回転行列の転置を返す．
  /*!
    \return	右から掛ける回転行列の転置
  */
    const Array2<T>&	Ut()		const	{return _Ut;}

  //! 2重対角化を行うために左から掛ける回転行列を返す．
  /*!
    \return	左から掛ける回転行列
  */
    const Array2<T>&	Vt()		const	{return _Vt;}

  //! 2重対角行列の対角成分を返す．
  /*!
    \return	対角成分
  */
    const Array<T>&	diagonal()	const	{return _Dt.sigma();}

  //! 2重対角行列の非対角成分を返す．
  /*!
    \return	非対角成分
  */
    const Array<T>&	off_diagonal()	const	{return _Et.sigma();}

    void		diagonalize()		;

  private:
    enum		{NITER_MAX = 30};
    
    bool		diagonal_is_zero(size_t n)		const	;
    bool		off_diagonal_is_zero(size_t n)		const	;
    void		initialize_rotation(size_t m, size_t n,
					    T& x, T& y)		const	;

    Householder<T>	_Dt;
    Householder<T>	_Et;
    Array<T>&		_diagonal;
    Array<T>&		_off_diagonal;
    T			_anorm;
    const Array2<T>&	_Ut;
    const Array2<T>&	_Vt;
};

//! 与えられた一般行列を2重対角化する．
/*!
  \param a	2重対角化する一般行列
*/
template <class T> template <class E, std::enable_if_t<rank<E>() == 2>*>
BiDiagonal<T>::BiDiagonal(const E& a)
    :_Dt(std::max(size<0>(a), size<1>(a)), 0),
     _Et(std::min(size<0>(a), size<1>(a)), 1),
     _diagonal(_Dt.sigma()),
     _off_diagonal(_Et.sigma()), _anorm(0),
     _Ut(size<0>(a) < size<1>(a) ? _Dt : _Et),
     _Vt(size<0>(a) < size<1>(a) ? _Et : _Dt)
{
    const auto&	A = evaluate(a);

    if (nrow() < ncol())
	for (size_t i = 0; i < nrow(); ++i)
	    for (size_t j = 0; j < ncol(); ++j)
		_Dt[i][j] = A[i][j];
    else
	for (size_t i = 0; i < nrow(); ++i)
	    for (size_t j = 0; j < ncol(); ++j)
		_Dt[j][i] = A[i][j];

  /* Householder reduction to bi-diagonal (off-diagonal in lower part) form */
    for (size_t m = 0; m < _Et.size(); ++m)
    {
	_Dt.apply_from_right(_Dt, m);
	_Et.apply_from_left(_Dt, m);
    }

    _Dt.make_transformation();	// Accumulate right-hand transformation: V
    _Et.make_transformation();	// Accumulate left-hand transformation: U

    for (size_t m = 0; m < _Et.size(); ++m)
    {
	T	anorm = std::fabs(_diagonal[m]) + std::fabs(_off_diagonal[m]);
	if (anorm > _anorm)
	    _anorm = anorm;
    }
}

//! 2重対角行列を対角化する（特異値分解）．
/*!
  対角成分は特異値となり，\f$\TUtvec{U}{}\f$と\f$\TUtvec{V}{}\f$
  の各行はそれぞれ右特異ベクトルと左特異ベクトルを与える．
  \throw std::runtime_error	指定した繰り返し回数を越えた場合に送出
*/ 
template <class T> void
BiDiagonal<T>::diagonalize()
{
    using namespace	std;
    
    for (size_t n = _Et.size(); n-- > 0; )
    {
	size_t	niter = 0;
	
#ifdef TU_DEBUG
	cerr << "******** n = " << n << " ********" << endl;
#endif
	while (!off_diagonal_is_zero(n))	// n > 0 here
	{
	    if (niter++ > NITER_MAX)
		throw runtime_error("TU::BiDiagonal::diagonalize(): Number of iteration exceeded maximum value");
	    
	  /* Find first m (< n) whose off-diagonal element is 0 */
	    size_t m = n;
	    do
	    {
		if (diagonal_is_zero(m-1))
		{ // If _diagonal[m-1] is zero, make _off_diagonal[m] zero.
		    T	x = _diagonal[m], y = _off_diagonal[m];
		    _off_diagonal[m] = 0.0;
		    for (size_t i = m; i <= n; ++i)
		    {
			Rotation<T>	rotD(m-1, i, x, -y);

			rotD.apply_from_left(_Dt);
			
			_diagonal[i] = -y*rotD.sin()
				     + _diagonal[i]*rotD.cos();
			if (i < n)
			{
			    x = _diagonal[i+1];
			    y = _off_diagonal[i+1]*rotD.sin();
			    _off_diagonal[i+1] *= rotD.cos();
			}
		    }
		    break;	// if _diagonal[n-1] is zero, m == n here.
		}
	    } while (!off_diagonal_is_zero(--m)); // 0 <= m < n < nrow() here.
	    if (m == n)
		break;		// _off_diagonal[n] has been made 0. Retry!

	  /* Set x and y which determine initial(i = m+1) plane rotation */
	    T	x, y;
	    initialize_rotation(m, n, x, y);
#ifdef TU_DEBUG
	    cerr << "--- m = " << m << ", n = " << n << "---"
		 << endl;
	    cerr << "  diagonal:     " << diagonal();
	    cerr << "  off-diagonal: " << off_diagonal();
#endif
	  /* Apply rotation P(i-1, i) for each i (i = m+1, n+2, ... , n) */
	    for (size_t i = m; ++i <= n; )
	    {
	      /* Apply rotation from left */
		Rotation<T>	rotE(i-1, i, x, y);
		
		rotE.apply_from_left(_Et);

		if (i > m+1)
		    _off_diagonal[i-1] = rotE.length();
		T	tmp = _diagonal[i-1];
		_diagonal[i-1]	 =  rotE.cos()*tmp
				 +  rotE.sin()*_off_diagonal[i];
		_off_diagonal[i] = -rotE.sin()*tmp
				 +  rotE.cos()*_off_diagonal[i];
		if (diagonal_is_zero(i))
		    break;		// No more Given's rotation needed.
		y		 =  rotE.sin()*_diagonal[i];
		_diagonal[i]	*=  rotE.cos();

		x = _diagonal[i-1];
		
	      /* Apply rotation from right to recover bi-diagonality */
		Rotation<T>	rotD(i-1, i, x, y);

		rotD.apply_from_left(_Dt);

		_diagonal[i-1] = rotD.length();
		tmp = _off_diagonal[i];
		_off_diagonal[i] =  tmp*rotD.cos() + _diagonal[i]*rotD.sin();
		_diagonal[i]	 = -tmp*rotD.sin() + _diagonal[i]*rotD.cos();
		if (i < n)
		{
		    if (off_diagonal_is_zero(i+1))
			break;		// No more Given's rotation needed.
		    y		        = _off_diagonal[i+1]*rotD.sin();
		    _off_diagonal[i+1] *= rotD.cos();

		    x		        = _off_diagonal[i];
		}
	    }
#ifdef TU_DEBUG
	    cerr << "  niter = " << niter << ": " << off_diagonal();
#endif
	}
    }

    for (size_t m = 0; m < _Et.size(); ++m)  // sort singular values and vectors
	for (size_t n = m+1; n < _Et.size(); ++n)
	    if (std::fabs(_diagonal[n]) > std::fabs(_diagonal[m]))
	    {
		swap(_diagonal[m], _diagonal[n]);
		for (size_t j = 0; j < _Et.size(); ++j)
		{
		    const T	tmp = _Et[m][j];
		    _Et[m][j] = _Et[n][j];
		    _Et[n][j] = -tmp;
		}
		for (size_t j = 0; j < _Dt.size(); ++j)
		{
		    const T	tmp = _Dt[m][j];
		    _Dt[m][j] = _Dt[n][j];
		    _Dt[n][j] = -tmp;
		}
	    }

    const auto l = std::min(_Dt.size() - 1, _Et.size());	// last index
    for (size_t m = 0; m < l; ++m)	// ensure positivity of all singular
	if (_diagonal[m] < 0.0)		// values except for the last one.
	{
	    _diagonal[m] = -_diagonal[m];
	    for (size_t j = 0; j < _Et.size(); ++j)
		_Et[m][j] = -_Et[m][j];

	    if (l < _Et.size())
	    {
		_diagonal[l] = -_diagonal[l];
		for (size_t j = 0; j < _Et.size(); ++j)
		    _Et[l][j] = -_Et[l][j];
	    }
	}
}

template <class T> bool
BiDiagonal<T>::diagonal_is_zero(size_t n) const
{
    return _Dt.sigma_is_zero(n, _anorm);
}

template <class T> bool
BiDiagonal<T>::off_diagonal_is_zero(size_t n) const
{
    return _Et.sigma_is_zero(n, _anorm);
}

template <class T> void
BiDiagonal<T>::initialize_rotation(size_t m, size_t n, T& x, T& y) const
{
    const T	g    = ((_diagonal[n]     + _diagonal[n-1])*
			(_diagonal[n]     - _diagonal[n-1])+
			(_off_diagonal[n] + _off_diagonal[n-1])*
			(_off_diagonal[n] - _off_diagonal[n-1]))
		     / (2.0*_diagonal[n-1]*_off_diagonal[n]),
      // Caution!! You have to ensure that _diagonal[n-1] != 0
      // as well as _off_diagonal[n].
		absg = std::fabs(g),
		gg1  = (absg > 1.0 ?
			absg * std::sqrt(1.0 + (1.0/absg)*(1.0/absg)) :
			std::sqrt(1.0 + absg*absg)),
		t    = (g > 0.0 ? g + gg1 : g - gg1);
    x = ((_diagonal[m] + _diagonal[n])*(_diagonal[m] - _diagonal[n]) -
	 _off_diagonal[n]*(_off_diagonal[n] + _diagonal[n-1]/t)) / _diagonal[m];
  //x = _diagonal[m];				// without shifting
    y = _off_diagonal[m+1];
}

/************************************************************************
*  class SVDecomposition<T>						*
************************************************************************/
//! 一般行列の特異値分解を表すクラス
/*!
  与えられた一般行列\f$\TUvec{A}{} \in \TUspace{R}{m\times n}\f$に対し
  て\f$\TUtvec{V}{}\TUvec{A}{}\TUvec{U}{}\f$が対角行列となるような2つの
  回転行列\f$\TUtvec{U}{} \in \TUspace{R}{n\times n}\f$,
  \f$\TUtvec{V}{} \in \TUspace{R}{m\times m}\f$を求める．
 */
template <class T>
class SVDecomposition : private BiDiagonal<T>
{
  private:
    using super		= BiDiagonal<T>;
    
  public:
    using element_type	= T;			//!< 成分の型

  public:
  //! 与えられた一般行列の特異値分解を求める．
  /*!
    \param a	特異値分解する一般行列
  */
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    SVDecomposition(const E_& a)
	:super(a)				{ super::diagonalize(); }

    using	super::nrow;
    using	super::ncol;
    using	super::Ut;
    using	super::Vt;
    using	super::diagonal;

  //! 特異値を求める．
  /*!
    \param i	絶対値の大きい順に並んだ特異値の1つを指定するindex
    \return	指定されたindexに対応する特異値
  */
    T		operator [](size_t i)	const	{ return diagonal()[i]; }
};

/************************************************************************
*  generic algorithms using SVDecomposition<T>				*
************************************************************************/
//! この行列の疑似逆行列を返す．
/*!
  \param cndnum	最大特異値に対する絶対値の割合がこれに達しない基底は無視
  \return	疑似逆行列，すなわち与えられた行列の特異値分解を
		\f$\TUvec{A}{} = \TUvec{V}{}\diag(\sigma_0,\ldots,\sigma_{n-1})
		\TUtvec{U}{}\f$とすると
		\f[
		  \TUvec{u}{0}\sigma_0^{-1}\TUtvec{v}{0} + \cdots +
		  \TUvec{u}{r}\sigma_{r-1}^{-1}\TUtvec{v}{r-1},
		  {\hskip 1em}\mbox{where}{\hskip 0.5em}
		  \TUabs{\sigma_1} > \epsilon\TUabs{\sigma_0},\ldots,
		  \TUabs{\sigma_{r-1}} > \epsilon\TUabs{\sigma_0}
		\f]
*/
template <class E, std::enable_if_t<rank<E>() == 2>* = nullptr> auto
pseudo_inverse(const E& A, element_t<E> cndnum=1.0e5)
{
    using element_type	= element_t<E>;
    
    SVDecomposition<element_type>	svd(A);
    Array2<element_type>		val(svd.ncol(), svd.nrow());

    for (size_t i = 0; i < svd.diagonal().size(); ++i)
	if (std::fabs(svd[i]) * cndnum > std::fabs(svd[0]))
	    val += (svd.Ut()[i] / svd[i]) % svd.Vt()[i];
    return val;
}

/************************************************************************
*  generic algorithms concerning with rotations				*
************************************************************************/
//! この3次元回転行列から各軸周りの回転角を取り出す．
/*!
  この行列を\f$\TUtvec{R}{}\f$とすると，
  \f[
    \TUvec{R}{} =
    \TUbeginarray{ccc}
      \cos\theta_z & -\sin\theta_z & \\
      \sin\theta_z &  \cos\theta_z & \\
      & & 1
    \TUendarray
    \TUbeginarray{ccc}
       \cos\theta_y & & \sin\theta_y \\
       & 1 & \\
      -\sin\theta_y & & \cos\theta_y
    \TUendarray
    \TUbeginarray{ccc}
      1 & & \\
      & \cos\theta_x & -\sin\theta_x \\
      & \sin\theta_x &  \cos\theta_x
    \TUendarray
  \f]
  なる\f$\theta_x, \theta_y, \theta_z\f$が回転角となる．
 \param theta_x	x軸周りの回転角(\f$ -\pi \le \theta_x \le \pi\f$)を返す．
 \param theta_y	y軸周りの回転角
	(\f$ -\frac{\pi}{2} \le \theta_y \le \frac{\pi}{2}\f$)を返す．
 \param theta_z	z軸周りの回転角(\f$ -\pi \le \theta_z \le \pi\f$)を返す．
 \throw invalid_argument	3次元正方行列でない場合に送出
*/
template <class E> std::enable_if_t<rank<E>() == 2>
rotation_angle(const E& expr, element_t<E>& theta_x,
	       element_t<E>& theta_y, element_t<E>& theta_z)
{
    using	T = element_t<E>;
    
    const auto&	Rt = evaluate(expr);
    
    if (size<0>(Rt) != 3 || size<1>(Rt) != 3)
	throw std::invalid_argument("TU::rot2angle: input matrix must be 3x3!!");

    if (Rt[0][0] == T(0) && Rt[0][1] == T(0))
    {
	theta_x = std::atan2(-Rt[2][1], Rt[1][1]);
	theta_y = (Rt[0][2] < 0 ? M_PI/2 : -M_PI/2);
	theta_z = T(0);
    }
    else
    {
	theta_x =  std::atan2(Rt[1][2], Rt[2][2]);
	theta_y = -std::asin(Rt[0][2]);
	theta_z =  std::atan2(Rt[0][1], Rt[0][0]);
    }
}

//! この3次元回転行列から回転角と回転軸を取り出す．
/*!
  この行列を\f$\TUtvec{R}{}\f$とすると，
  \f[
    \TUtvec{R}{} \equiv \TUvec{I}{3}\cos\theta
    + \TUvec{n}{}\TUtvec{n}{}(1 - \cos\theta)
    - \TUskew{n}{}\sin\theta
  \f]
  なる\f$\theta~(0 \le \theta \le \pi)\f$と\f$\TUvec{n}{}\f$が
  それぞれ回転角と回転軸となる．
  \param expr
  \param c	回転角のcos値，すなわち\f$\cos\theta\f$を返す．
  \param s	回転角のsin値，すなわち\f$\sin\theta (\ge 0)\f$を返す．
  \return	回転軸を表す3次元単位ベクトル，すなわち\f$\TUvec{n}{}\f$
  \throw std::invalid_argument	3x3行列でない場合に送出
*/
template <class E> std::enable_if_t<rank<E>() == 2, Array<element_t<E>, 3> >
rotation_axis(const E& expr, element_t<E>& c, element_t<E>& s)
{
    using	T = element_t<E>;
    
    const auto&	Rt = expr;

    if (size<0>(Rt) != 3 || size<1>(Rt) != 3)
	throw std::invalid_argument("TU::rot2angle: input matrix must be 3x3!!");

  // Compute cosine and sine of rotation angle.
    const T	trace = Rt[0][0] + Rt[1][1] + Rt[2][2];
    c = T(0.5) * (trace - T(1));
    s = T(0.5) * std::sqrt((T(1) + trace)*(T(3) - trace));

  // Compute rotation axis.
    Array<T, 3>	n({Rt[1][2]-Rt[2][1], Rt[2][0]-Rt[0][2], Rt[0][1]-Rt[1][0]});
    return std::move(normalize(n));
}

//! この3次元回転行列から回転角と回転軸を取り出す．
/*!
  この行列を\f$\TUtvec{R}{}\f$とすると，
  \f[
    \TUtvec{R}{} \equiv \TUvec{I}{3}\cos\theta
    + \TUvec{n}{}\TUtvec{n}{}(1 - \cos\theta)
    - \TUskew{n}{}\sin\theta
  \f]
  なる\f$\theta~(0 \le \theta \le \pi)\f$と\f$\TUvec{n}{}\f$が
  それぞれ回転角と回転軸となる．
 \return			回転角と回転軸を表す3次元ベクトル，すなわち
				\f$\theta\TUvec{n}{}\f$
 \throw invalid_argument	3x3行列でない場合に送出
*/
template <class E> std::enable_if_t<rank<E>() == 2, Array<element_t<E>, 3> >
rotation_axis(const E& expr)
{
    using	T = element_t<E>;
    
    const auto&	Rt = expr;

    if (size<0>(Rt) != 3 || size<1>(Rt) != 3)
	throw std::invalid_argument("TU::rot2angle: input matrix must be 3x3!!");

    const T	trace = Rt[0][0] + Rt[1][1] + Rt[2][2],
		s2 = std::sqrt((T(1) + trace)*(T(3) - trace));	// 2*sin
    if (s2 + T(1) == T(1))			// sin << 1 ?
	return Array<T, 3>();			// zero vector
    
    Array<T, 3>	axis({Rt[1][2]-Rt[2][1], Rt[2][0]-Rt[0][2], Rt[0][1]-Rt[1][0]});

    return std::move(axis *= (std::atan2(s2, trace - T(1)) / s2));
}

//! この3次元回転行列から四元数を取り出す．
/*!
  この行列を\f$\TUtvec{R}{}\f$とすると，
  \f[
    \TUtvec{R}{} \equiv \TUvec{I}{3}\cos\theta
    + \TUvec{n}{}\TUtvec{n}{}(1 - \cos\theta)
    - \TUskew{n}{}\sin\theta
  \f]
  なる\f$\theta~(0 \le \theta \le \pi)\f$と\f$\TUvec{n}{}\f$に対して，四元数は
  \f[
    \TUvec{q}{} \equiv
    \TUbeginarray{c}
      \cos\frac{\theta}{2} \\ \TUvec{n}{}\sin\frac{\theta}{2}
    \TUendarray
  \f]
  と定義される．
 \return			四元数を表す4次元単位ベクトル
 \throw invalid_argument	3x3行列でない場合に送出
*/
template <class E> std::enable_if_t<rank<E>() == 2, Array<element_t<E>, 4> >
quaternion(const E& expr)
{
    using	T = element_t<E>;
    
    const auto&	Rt = expr;

    if (size<0>(Rt) != 3 || size<1>(Rt) != 3)
	throw std::invalid_argument("TU::rot2angle: input matrix must be 3x3!!");

    Array<T, 4>	q;
    q[0] = T(0.5) * std::sqrt(trace(Rt) + T(1));
    if (q[0] + T(1) == T(1))	// q[0] << 1 ?
    {
	Array<T, 3>	eval;
	slice(q, 1, 3) = eigen(Rt, eval, false)[0];
    }
    else
    {
	const auto	S = transpose(Rt) - Rt;
	q[1] = T(0.25) * S[2][1] / q[0];
	q[2] = T(0.25) * S[0][2] / q[0];
	q[3] = T(0.25) * S[1][0] / q[0];
    }

    return q;
}

//! 2次元回転行列を生成する．
/*!
  \param theta	回転角
  \return	生成された回転行列，すなわち
		\f[
		  \TUtvec{R}{} \equiv
		  \TUbeginarray{cc}
		  \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta
		  \TUendarray
		\f]
*/
template <class T>
inline std::enable_if_t<std::is_floating_point<T>::value, Array2<T, 2, 2> >
rotation(T theta)
{
    Array2<T, 2, 2>	Qt;
    Qt[0][0] = Qt[1][1] = std::cos(theta);
    Qt[0][1] = std::sin(theta);
    Qt[1][0] = -Qt[0][1];

    return Qt;
}

//! 3次元回転行列を生成する．
/*!
  \param theta_x	x軸まわりの回転角
  \param theta_y	y軸まわりの回転角
  \param theta_z	z軸まわりの回転角
  \return		生成された回転行列，すなわち
			\f[
			\TUtvec{R}{} \equiv \TUtvec{R}{x}(\theta_x)
			\TUtvec{R}{y}(\theta_y)\TUtvec{R}{z}(\theta_z)
			\f]
*/
template <class T>
inline std::enable_if_t<std::is_floating_point<T>::value, Array2<T, 3, 3> >
rotation(T theta_x, T theta_y, T theta_z)
{
    const auto	cx = std::cos(theta_x);
    const auto	sx = std::sin(theta_x);
    const auto	cy = std::cos(theta_y);
    const auto	sy = std::sin(theta_y);
    const auto	cz = std::cos(theta_z);
    const auto	sz = std::sin(theta_z);
    
    Array2<T, 3, 3>	Qt;
    Qt[0][0] =  cy * cz;
    Qt[0][1] =  cy * sz;
    Qt[0][2] = -sy;
    Qt[1][0] =  sx * sy * cz - cx * sz;
    Qt[1][1] =  sx * sy * sz + cx * cz;
    Qt[1][2] =	sx * cy;
    Qt[2][0] =  cx * sy * cz + sx * sz;
    Qt[2][1] =  cx * sy * sz - sx * cz;
    Qt[2][2] =	cx * cy;

    return Qt;
}

//! 3次元回転行列を生成する．
/*!
  \param n	回転軸を表す3次元単位ベクトル
  \param c	回転角のcos値
  \param s	回転角のsin値
  \return	生成された回転行列，すなわち
		\f[
		  \TUtvec{R}{} \equiv \TUvec{I}{3}\cos\theta
		  + \TUvec{n}{}\TUtvec{n}{}(1 - \cos\theta)
		  - \TUskew{n}{}\sin\theta
		\f]
*/
template <class E>
inline std::enable_if_t<rank<E>() == 1, Array2<element_t<E>, 3, 3>>
rotation(const E& n, element_t<E> c, element_t<E> s)
{
    if (size(n) != 3)
	throw std::invalid_argument("TU::Rt: dimension of the argument \'n\' must be 3");
    const auto			nn = evaluate(n);
    Array2<element_t<E>, 3, 3>	Qt = nn % nn;
    Qt *= (1.0 - c);
    Qt[0][0] += c;
    Qt[1][1] += c;
    Qt[2][2] += c;
    Qt[0][1] += nn[2] * s;
    Qt[0][2] -= nn[1] * s;
    Qt[1][0] -= nn[2] * s;
    Qt[1][2] += nn[0] * s;
    Qt[2][0] += nn[1] * s;
    Qt[2][1] -= nn[0] * s;

    return Qt;
}

//! 3次元回転行列を生成する．
/*!
  \param v	回転角と回転軸を表す3次元ベクトルまたは四元数を表す4次元単位ベクトル
  \return	生成された回転行列，すなわち3次元ベクトルの場合は
		\f[
		  \TUtvec{R}{} \equiv \TUvec{I}{3}\cos\theta
		  + \TUvec{n}{}\TUtvec{n}{}(1 - \cos\theta)
		  - \TUskew{n}{}\sin\theta,
		  {\hskip 1em}\mbox{where}{\hskip 0.5em}
		  \theta \equiv \TUnorm{\TUvec{v}{}},~
		  \TUvec{n}{} \equiv \frac{\TUvec{v}{}}{\TUnorm{\TUvec{v}{}}}
		\f]
		4次元単位ベクトルの場合は
		\f[
		  \TUtvec{R}{} \equiv
		  \TUvec{I}{3}(q_0^2 - \TUtvec{q}{}\TUvec{q}{})
		  + 2\TUvec{q}{}\TUtvec{q}{}
		  - 2q_0\TUskew{q}{},
		  {\hskip 1em}\mbox{where}{\hskip 0.5em}
		  q_0 \equiv v_0,~
		  \TUvec{q}{} \equiv [v_1,~v_2,~v_3]^\top
		\f]
*/
template <class E>
inline std::enable_if_t<rank<E>() == 1, Array2<element_t<E>, 3, 3>>
rotation(const E& v)
{
    using	T = element_t<E>;
    
    if (std::size(v) == 4)		// quaternion ?
    {
	const T		q0 = v[0];
	Array<T, 3>	q{{v[1], v[2], v[3]}};
	Array2<T, 3, 3>	Qt = (T(2) * q) % q;
	const T		c = q0*q0 - square(q);
	Qt[0][0] += c;
	Qt[1][1] += c;
	Qt[2][2] += c;
	q *= (T(2) * q0);
	Qt[0][1] += q[2];
	Qt[0][2] -= q[1];
	Qt[1][0] -= q[2];
	Qt[1][2] += q[0];
	Qt[2][0] += q[1];
	Qt[2][1] -= q[0];

	return Qt;
    }

    const T	theta = length(v);
    if (theta + T(1) == T(1))		// theta << 1 ?
	return diag<3>(T(1));
    else
    {
	const T	c = std::cos(theta), s = std::sin(theta);
	return rotation(v / theta, c, s);
    }
}

/************************************************************************
*  type definitions for convenience					*
************************************************************************/
template <class T, size_t N=0, class ALLOC=std::allocator<T> >
using Vector = array<T, ALLOC, N>;				//!< 1次元配列

template <class T, size_t R=0, size_t C=0, class ALLOC=std::allocator<T> >
using Matrix = array<T, ALLOC, R, C>;				//!< 2次元配列

using Vector2s	= Vector<short,  2>;
using Vector2i	= Vector<int,    2>;
using Vector2f	= Vector<float,  2>;
using Vector2d	= Vector<double, 2>;
using Vector3s	= Vector<short,  3>;
using Vector3i	= Vector<int,    3>;
using Vector3f	= Vector<float,  3>;
using Vector3d	= Vector<double, 3>;
using Vector4s	= Vector<short,  4>;
using Vector4i	= Vector<int,    4>;
using Vector4f	= Vector<float,  4>;
using Vector4d	= Vector<double, 4>;
using Matrix22f	= Matrix<float,  2, 2>;
using Matrix22d	= Matrix<double, 2, 2>;
using Matrix23f	= Matrix<float,  2, 3>;
using Matrix23d	= Matrix<double, 2, 3>;
using Matrix33f	= Matrix<float,  3, 3>;
using Matrix33d	= Matrix<double, 3, 3>;
using Matrix34f	= Matrix<float,  3, 4>;
using Matrix34d	= Matrix<double, 3, 4>;
using Matrix44f	= Matrix<float,  4, 4>;
using Matrix44d	= Matrix<double, 4, 4>;
}
#endif	// !TU_VECTORPP_H
