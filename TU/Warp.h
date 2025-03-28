/*!
  \file		Warp.h
  \author	Toshio UESHIBA
  \brief	クラス TU::Warp の定義と実装
*/
#ifndef	TU_WARP_H
#define	TU_WARP_H

#include "TU/simd/Array++.h"
#include "TU/Image++.h"
#include "TU/Camera++.h"
#if defined(USE_TBB)
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#endif

namespace TU
{
/************************************************************************
*  class Warp								*
************************************************************************/
//! 画像を変形するためのクラス
class Warp
{
  private:
#if defined(SIMD)
    template <class P, class T>
    static P	ptr(T* p)
		{
		    return reinterpret_cast<P>(p);
		}
    template <class P, class T, bool ALIGNED>
    static P	ptr(simd::iterator_wrapper<T*, ALIGNED> p)
		{
		    return reinterpret_cast<P>(p.base());
		}
#endif    
    struct FracArray
    {
	FracArray(size_t d=0)
	    :us(d), vs(d), du(d), dv(d), lmost(0)	{}

	size_t		width()			const	{return us.size();}
	void		resize(size_t d)		;
#if defined(SIMD)
	template <class T>
	using allocator	= simd::allocator<T>;
#else
	template <class T>
	using allocator	= std::allocator<T>;
#endif
	Array<short,  0, allocator<short> >	us, vs;
	Array<u_char, 0, allocator<u_char> >	du, dv;
	size_t					lmost;
    };
    
    template <class IN>
    class Interpolate
    {
      public:
	using value_type = typename std::iter_value_t<IN>::value_type;
    
      public:
	Interpolate(IN in)	:_in(in)		{}

	template <class S, class T>
	auto	operator ()(S u, S v, T du, T dv) const
		{
		    const auto	ue = u + S(1);
		    const auto	ve = v + S(1);
		    
		    return interpolate(interpolate(lookup(u,  v),
						   lookup(ue, v), du),
				       interpolate(lookup(u,  ve),
						   lookup(ue, ve), du),
				       dv);
		}
	template <class S, class T>
	auto	operator ()(const std::tuple<S, S, T, T>& arg) const
		{
		    return (*this)(std::get<0>(arg), std::get<1>(arg),
				   std::get<2>(arg), std::get<3>(arg));
		}
    
      private:
	value_type	lookup(int u, int v) const
			{
			    return (*(_in + v))[u];
			}

	template <class S>
	static S	interpolate(S x, S y, S d)
			{
			    return x + (d*(y - x) >> 7);
			}
	template <class E, class T>
	static RGB_<E>	interpolate(const RGB_<E>& x, const RGB_<E>& y, T d)
			{
			    return {interpolate(x.r, y.r, d),
				    interpolate(x.g, y.g, d),
				    interpolate(x.b, y.b, d)};
			}
	template <class T>
	static YUV444	interpolate(const YUV444& x, const YUV444& y, T d)
			{
			    return {interpolate(x.y, y.y, d),
				    interpolate(x.u, y.u, d),
				    interpolate(x.v, y.v, d)};
			}
#if defined(SIMD)
	template <class V=value_type, class T>
	std::enable_if_t<std::is_integral<V>::value, simd::vec<T> >
			lookup(simd::vec<T> u, simd::vec<T> v) const
			{
			    return simd::lookup(std::cbegin(*_in),
						v, u, _in.stride());
			}
	template <class V=value_type>
	std::enable_if_t<!std::is_integral<V>::value, simd::Is32vec>
			lookup(simd::Is32vec u, simd::Is32vec v) const
			{
			    return simd::lookup(
					Warp::ptr<const int32_t*>(
					    std::cbegin(*_in)),
					v, u, _in.stride());
			}
	template <class V=value_type>
	static std::enable_if_t<!std::is_integral<V>::value, simd::Is32vec>
			interpolate(simd::Is32vec x,
				    simd::Is32vec y, simd::Is32vec d)
			{
			    using namespace	simd;

			    const auto	dl = dup<false>(d);
			    const auto	dh = dup<true >(d);

			    return cast<int32_t>(
				cvt<uint8_t>(
				    interpolate(
					cvt<int16_t, false>(cast<uint8_t>(x)),
					cvt<int16_t, false>(cast<uint8_t>(y)),
					cvt<int16_t>(dup<false>(dl),
						     dup<true >(dl))),
				    interpolate(
					cvt<int16_t, true >(cast<uint8_t>(x)),
					cvt<int16_t, true >(cast<uint8_t>(y)),
					cvt<int16_t>(dup<false>(dh),
						     dup<true >(dh)))));
			}
#endif
      private:
	const IN	_in;
    };

#if defined(USE_TBB)
    template <class IN, class OUT>
    class WarpLine
    {
      public:
	WarpLine(const Warp& warp, IN in, OUT out)
	    :_warp(warp), _in(in), _out(out)				{}
	
	void	operator ()(const tbb::blocked_range<size_t>& r) const
		{
		    auto	out = _out + r.begin();
		    for (auto v = r.begin(); v != r.end(); ++v)
		    {
			using	std::begin;
			
			_warp.warpLine(_in, begin(*out), _warp._fracs[v]);
			++out;
		    }
		}

      private:
	const Warp&	_warp;
	const IN	_in;
	const OUT	_out;
    };
#endif

  public:
  //! 画像変形オブジェクトを生成する．
    Warp()	:_fracs(), _width(0)			{}

  //! 出力画像の幅を返す．
  /*!
    return	出力画像の幅
  */
    size_t	width()				const	{return _width;}

  //! 出力画像の高さを返す．
  /*!
    return	出力画像の高さ
  */
    size_t	height()			const	{return _fracs.size();}
    
    size_t	lmost(size_t v)			const	;
    size_t	rmost(size_t v)			const	;

    template <class T>
    void	initialize(const Matrix<T, 3, 3>& Htinv,
			   size_t inWidth,  size_t inHeight,
			   size_t outWidth, size_t outHeight)		;
    template <class I>
    void	initialize(const typename I::matrix33_type& Htinv,
			   const I& intrinsic,
			   size_t inWidth,  size_t inHeight,
			   size_t outWidth, size_t outHeight)		;
    template <class IN, class OUT>
    void	operator ()(IN in, OUT out)			const	;
    Vector2f	operator ()(size_t u, size_t v)			const	;

  private:
    template <class IN, class OUT>
    void	warpLine(IN in, OUT out, const FracArray& frac)	const	;
    
  private:
    Array<FracArray>	_fracs;
    size_t		_width;
};

inline void
Warp::FracArray::resize(size_t d)
{
    us.resize(d);
    vs.resize(d);
    du.resize(d);
    dv.resize(d);
}

//! 出力画像における指定された行の有効左端位置を返す．
/*!
  入力画像が矩形でも出力画像も矩形とは限らないので，出力画像の一部しか
  入力画像の値域(有効領域)とならない．本関数は，出力画像の指定された行
  について，その有効領域の左端となる画素位置を返す．
  \param v	行を指定するindex
  \return	左端位置
*/
inline size_t
Warp::lmost(size_t v) const
{
    return _fracs[v].lmost;
}

//! 出力画像における指定された行の有効右端位置の次を返す．
/*!
  入力画像が矩形でも出力画像も矩形とは限らないので，出力画像の一部しか
  入力画像の値域(有効領域)とならない．本関数は，出力画像の指定された行
  について，その有効領域の右端の右隣となる画素位置を返す．
  \param v	行を指定するindex
  \return	右端位置の次
*/
inline size_t
Warp::rmost(size_t v) const
{
    return _fracs[v].lmost + _fracs[v].width();
}

//! 画像を射影変換するための行列を設定する．
/*!
  入力画像点uは射影変換
  \f[
    \TUbeginarray{c} \TUvec{v}{} \\ 1 \TUendarray \simeq
    \TUvec{H}{} \TUbeginarray{c} \TUvec{u}{} \\ 1 \TUendarray
  \f]
  によって出力画像点vに写される．
  \param Htinv		変形を指定する3x3射影変換行列の逆行列の転置，すなわち
			\f$\TUtinv{H}{}\f$
  \param inWidth	入力画像の幅
  \param inHeight	入力画像の高さ
  \param outWidth	出力画像の幅
  \param outHeight	出力画像の高さ
*/
template <class T> inline void
Warp::initialize(const Matrix<T, 3, 3>& Htinv,
		 size_t inWidth,  size_t inHeight,
		 size_t outWidth, size_t outHeight)
{
    initialize(Htinv, IntrinsicBase<T>(),
	       inWidth, inHeight, outWidth, outHeight);
}

//! 画像の非線形歪みを除去した後に射影変換を行うための行列とカメラ内部パラメータを設定する．
/*!

  canonical座標xから画像座標uへの変換が\f$\TUvec{u}{} = {\cal
  K}(\TUvec{x}{})\f$ と表されるカメラ内部パラメータについて，その線形変
  換部分を表す3x3上半三角行列をKとすると，
  \f[
    \TUbeginarray{c} \TUbar{u}{} \\ 1 \TUendarray =
    \TUvec{K}{}
    \TUbeginarray{c} {\cal K}^{-1}(\TUvec{u}{}) \\ 1 \TUendarray
  \f]
  によって画像の非線形歪みだけを除去できる．本関数は，この歪みを除去した画像点を
  射影変換Hによって出力画像点vに写すように変形パラメータを設定する．すなわち，
  全体の変形は
  \f[
    \TUbeginarray{c} \TUvec{v}{} \\ 1 \TUendarray \simeq
    \TUvec{H}{}\TUvec{K}{}
    \TUbeginarray{c} {\cal K}^{-1}(\TUvec{u}{}) \\ 1 \TUendarray
  \f]
  となる．
  \param Htinv		変形を指定する3x3射影変換行列の逆行列の転置
  \param intrinsic	入力画像に加えれられている放射歪曲を表すカメラ内部パラメータ
  \param inWidth	入力画像の幅
  \param inHeight	入力画像の高さ
  \param outWidth	出力画像の幅
  \param outHeight	出力画像の高さ
*/
template <class I> void
Warp::initialize(const typename I::matrix33_type& Htinv, const I& intrinsic,
		 size_t inWidth,  size_t inHeight,
		 size_t outWidth, size_t outHeight)
{
    using namespace std;

    using intrinsic_type	= I;
    using vector3_type		= Vector<typename intrinsic_type::element_type,
					 3>;
    using matrix33_type		= typename intrinsic_type::matrix33_type;
    
    _fracs.resize(outHeight);
    _width = outWidth;
    
  /* Compute frac for each pixel. */
    const matrix33_type	HKtinv = Htinv * intrinsic.Ktinv();
    vector3_type	leftmost = HKtinv[2];
    for (size_t v = 0; v < height(); ++v)
    {
	auto		x = leftmost;
	FracArray	frac(width());
	size_t		n = 0;
	for (size_t u = 0; u < width(); ++u)
	{
	    const auto	m = intrinsic.u({x[0]/x[2], x[1]/x[2]});
	    if (0.0 <= m[0] && m[0] <= inWidth - 2 &&
		0.0 <= m[1] && m[1] <= inHeight - 2)
	    {
		if (n == 0)
		    frac.lmost = u;
		frac.us[n] = (short)floor(m[0]);
		frac.vs[n] = (short)floor(m[1]);
		frac.du[n] = (u_char)floor((m[0] - floor(m[0])) * 128.0);
		frac.dv[n] = (u_char)floor((m[1] - floor(m[1])) * 128.0);
		++n;
	    }
	    x += HKtinv[0];
	}

	_fracs[v].resize(n);
	_fracs[v].lmost = frac.lmost;

	for (size_t u = 0; u < n; ++u)
	{
	    _fracs[v].us[u] = frac.us[u];
	    _fracs[v].vs[u] = frac.vs[u];
	    _fracs[v].du[u] = frac.du[u];
	    _fracs[v].dv[u] = frac.dv[u];
	}

	leftmost += HKtinv[1];
    }
}

//! 入力画像を変形して出力画像に出力する．
/*!
  \param in	入力画像の最初の行を指す反復子
  \param out	出力画像の最初の行を指す反復子
*/
template <class IN, class OUT> void
Warp::operator ()(IN in, OUT out) const
{
#if defined(USE_TBB)
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _fracs.size(), 1),
		      WarpLine<IN, OUT>(*this, in, out));
#else
    for (const auto& frac : _fracs)
    {
	using	std::begin;
	
	warpLine(in, begin(*out), frac);
	++out;
    }
#endif
}
    
//! 出力画像点を指定してそれにマップされる入力画像点の2次元座標を返す．
/*!
  \param u	出力画像点の横座標
  \param v	出力画像点の縦座標
  \return	出力画像点(u, v)にマップされる入力画像点の2次元座標
*/
inline Vector2f
Warp::operator ()(size_t u, size_t v) const
{
    const auto&	fracs = _fracs[v];
    return {float(fracs.us[u]) + float(fracs.du[u]) / 128.0f,
	    float(fracs.vs[u]) + float(fracs.dv[u]) / 128.0f};
}

template <class IN, class OUT> void
Warp::warpLine(IN in, OUT out, const FracArray& frac) const
{
    Interpolate<IN>	interpolate(in);
    
    auto	u  = frac.us.cbegin();
    auto	ue = frac.us.cend();
    auto	v  = frac.vs.cbegin();
    auto	du = frac.du.cbegin();
    auto	dv = frac.dv.cbegin();
    out += frac.lmost;
#if defined(SIMD)
    using	value_type = typename Interpolate<IN>::value_type;
    using	T = std::conditional_t<(sizeof(value_type) > sizeof(int16_t)),
				       int32_t, int16_t>;
    using	O = std::conditional_t<(sizeof(value_type) > sizeof(int16_t)),
				       int32_t, std::iter_value_t<OUT> >;
	
    const auto	n = simd::vec<u_char>::floor(std::distance(u, ue));

    simd::transform<T>(interpolate,
		       simd::make_accessor(ptr<O*>(out)),
		       simd::make_accessor(u),
		       simd::make_accessor(u + n),
		       simd::make_accessor(v),
		       simd::make_accessor(du),
		       simd::make_accessor(dv));
    simd::empty();

    std::advance(u,   n);
    std::advance(v,   n);
    std::advance(du,  n);
    std::advance(dv,  n);
    std::advance(out, n);
#endif
    while (u != ue)
	*out++ = interpolate(*u++, *v++, *du++, *dv++);
}

}	// namespace TU
#endif	// !TU_WARP_H
