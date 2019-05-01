/*!
  \mainpage	libTUTools++ - 配列，ベクトル，行列，画像等の基本的なデータ型とそれに付随したアルゴリズムを収めたライブラリ
  \anchor	libTUTools
  \author	Toshio UESHIBA
  \copyright	2002-2017年 植芝俊夫

  \section abstract 概要
  libTUTools++は，配列，ベクトル，行列，画像等の基本的なデータ型とそれ
  に付随したアルゴリズムを収めたライブラリである．現在実装されている主
  要なクラスおよび関数はおおまかに以下の分野に分類される．

  <b>多次元配列</b>
  - #TU::array
  - #TU::Array
  - #TU::Array2

  <b>ベクトルと行列および線形計算</b>
  - #TU::Vector
  - #TU::Matrix
  - #TU::LUDecomposition
  - #TU::Householder
  - #TU::QRDecomposition
  - #TU::TriDiagonal
  - #TU::BiDiagonal
  - #TU::SVDecomposition
  - #TU::BlockDiagonalMatrix
  - #TU::SparseMatrix
  - #TU::BandMatrix

  <b>非線形最適化</b>
  - #TU::NullConstraint
  - #TU::ConstNormConstraint
  - #TU::minimizeSquare(const F&, const G&, AT&, size_t, double)
  - #TU::minimizeSquareSparse(const F&, const G&, ATA&, IB, IB, size_t, double)

  <b>RANSAC</b>
  - #TU::ransac(const SAMPLER&, MODEL&, CONFORM&&, T, T)
  - #TU::ransac(IN, IN, MODEL&, CONFORM&&, T, T)

  <b>グラフカット</b>
  - #boost::GraphCuts

  <b>動的計画法</b>
  - #TU::DP

  <b>点，直線，平面等の幾何要素とその変換</b>
  - #TU::Point2
  - #TU::Point3
  - #TU::HyperPlane
  - #TU::Normalize
  - #TU::Projectivity
  - #TU::Affinity
  - #TU::Rigidity
  - #TU::BoundingBox
  
  <b>投影の幾何</b>
  - #TU::IntrinsicBase
  - #TU::IntrinsicWithFocalLength
  - #TU::IntrinsicWithEuclideanImagePlane
  - #TU::Intrinsic
  - #TU::IntrinsicWithDistortion
  - #TU::CanonicalCamera
  - #TU::Camera
  
  <b>画素と画像</b>
  - #TU::RGB
  - #TU::BGR
  - #TU::RGBA
  - #TU::ABGR
  - #TU::ARGB
  - #TU::BGRA
  - #TU::YUV444
  - #TU::YUV422
  - #TU::YUV411
  - #TU::Image
  - #TU::GenericImage
  - #TU::ComplexImage
  - #TU::CCSImage
  - #TU::Movie

  <b>画像処理</b>
  - #TU::Rectify
  - #TU::EdgeDetector
  - #TU::IntegralImage
  - #TU::DiagonalIntegralImage
  
  <b>画像に限らない信号処理</b>
  - #TU::Warp
  - #TU::Filter2
  - #TU::SeparableFilter2
  - #TU::IIRFilter
  - #TU::BidirectionalIIRFilter
  - #TU::BidirectionalIIRFilter2
  - #TU::DericheConvolver
  - #TU::DericheConvolver2
  - #TU::GaussianConvolver
  - #TU::GaussianConvolver2
  - #TU::BoxFilter
  - #TU::BoxFilter2
  - #TU::GuidedFilter
  - #TU::GuidedFilter2
  - #TU::FIRFilter
  - #TU::FIRGaussianConvolver
  - #TU::WeightedMedianFilter
  - #boost::TreeFilter

  <b>特殊データ構造</b>
  - #TU::List
  - #TU::PSTree
  - #TU::NDTree
  
  <b>Bezier曲線とBezier曲面</b>
  - #TU::BezierCurve
  - #TU::BezierSurface
  
  <b>B-spline曲線とB-spline曲面</b>
  - #TU::BSplineKnots
  - #TU::BSplineCurve
  - #TU::BSplineSurface
  
  <b>メッシュ</b>
  - #TU::Mesh

  <b>アルゴリズム</b>
  - #TU::diff(const T&, const T&)
  - #TU::gcd(S, T)
  - #TU::gcd(S, T, ARGS...)
  - #TU::lcm(S, T)
  - #TU::lcm(S, T, ARGS...)
  - #TU::for_each(FUNC, size_t, ITER...)
  - #TU::fill(ITER, size_t, const T&)
  - #TU::copy(IN, size_t, OUT)
  - #TU::inner_product(ITER0, size_t, ITER1, T)
  - #TU::square(ITER, size_t)
  - #TU::square(T)
  - #TU::op3x3(ROW, ROW, OP)
  - #TU::max3x3(COL, COL, COL)
  - #TU::min3x3(COL, COL, COL)
  - #TU::mopOpen(ROW, ROW, size_t)
  - #TU::mopClose(ROW, ROW, size_t)
  - #TU::inclusive_scan(IN, IN, OUT, OP)
  
  <b>関数オブジェクト</b>
  - #TU::identity
  - #TU::plus_assign
  - #TU::minus_assign
  - #TU::multiplies_assign
  - #TU::divides_assign
  - #TU::modulus_assign
  - #TU::bit_and_assign
  - #TU::bit_or_assign
  - #TU::bit_xor_assign

  <b>反復子</b>
  - #TU::make_mbr_iterator()
  - #TU::make_first_iterator()
  - #TU::make_second_iterator()
  - #TU::zip_iterator
  - #TU::assignment_iterator
  - #TU::row2col
  - #TU::column_iterator
  - #TU::ring_iterator
  - #TU::box_filter_iterator
  - #TU::iir_filter_iterator
  
  <b>マニピュレータ</b>
  - #TU::skipl(std::istream&)
  - #TU::IOManip
  - #TU::IManip1
  - #TU::OManip1
  - #TU::IOManip1
  - #TU::IManip2
  - #TU::OManip2

  <b>ストリーム
  - #TU::fdistream
  - #TU::fdostream
  - #TU::fdstream
  
  <b>シリアルインタフェース</b>
  - #TU::Serial
  - #TU::TriggerGenerator
  - #TU::PM16C_04
  - #TU::SHOT602

  <b>SIMD命令</b>
  - #TU::simd::vec
  
  \file		algorithm.h
  \brief	各種アルゴリズムの定義と実装
*/
#ifndef TU_ALGORITHM_H
#define TU_ALGORITHM_H

#include <cmath>		// for std::sqrt()
#include <iterator>		// for std::iterator_traits<ITER>
#include <type_traits>		// for std::common_type<TYPES....>
#include <algorithm>		// for std::copy(), std::copy_n(),...
#include <numeric>		// for std::inner_product()
#include <utility>		// for std::forward()
#include <functional>		// for std::plus<>
#ifdef TU_DEBUG
#  include <iostream>
#endif

namespace std
{
#if __cplusplus < 201700L
namespace detail
{
  template <class IN, class OUT, class SIZE, class GEN> OUT
  sample(IN in, IN ie, input_iterator_tag, OUT out, random_access_iterator_tag,
	 SIZE n, GEN&& gen)
  {
      using distrib_type = std::uniform_int_distribution<SIZE>;
      using param_type   = typename distrib_type::param_type;

      distrib_type	distrib{};
      SIZE		nsampled = 0;

    // 最初のn個をコピー
      for (; in != ie && nsampled != n; ++in, ++nsampled)
	  out[nsampled] = *in;

      for (auto ninputs = nsampled; in != ie; ++in, ++ninputs)
      {
	  const auto	i = distrib(gen, param_type{0, ninputs});
	  if (i < n)
	      out[i] = *in;
      }

      return out + nsampled;
  }

  template<class IN, class OUT, class CAT, class SIZE, class GEN> OUT
  sample(IN in, IN ie, forward_iterator_tag, OUT out, CAT, SIZE n, GEN&& gen)
  {
      using distrib_type = std::uniform_int_distribution<SIZE>;
      using param_type	 = typename distrib_type::param_type;

      distrib_type	distrib{};
      SIZE		nunsampled = std::distance(in, ie);

      if (n > nunsampled)
	  n = nunsampled;
      
      for (; n != 0; ++in)
	  if (distrib(gen, param_type{0, --nunsampled}) < n)
	  {
	      *out++ = *in;
	      --n;
	  }

      return out;
  }
}	// namespace detail
    
/// Take a random sample from a population.
template<class IN, class OUT, class SIZE, class GEN> OUT
sample(IN in, IN ie, OUT out, SIZE n, GEN&& gen)
{
    using in_cat  = typename std::iterator_traits<IN>::iterator_category;
    using out_cat = typename std::iterator_traits<OUT>::iterator_category;

    return detail::sample(in, ie, in_cat{}, out, out_cat{},
			  n, std::forward<GEN>(gen));
}
#endif    
}	// namespace std

//! libTUTools++ のクラスや関数等を収める名前空間
namespace TU
{
#ifdef TU_DEBUG
template <class ITER, size_t SIZE>	class range;
template <class E>			class sizes_holder;

template <class E>
sizes_holder<E>	print_sizes(const E& expr);
template <class E>
std::ostream&	operator <<(std::ostream& out, const sizes_holder<E>& holder);
#endif

/************************************************************************
*  generic algorithms							*
************************************************************************/
//! 2つの引数の差の絶対値を返す．
template <class T> inline constexpr T
diff(const T& a, const T& b)
{
    return (a > b ? a - b : b - a);
}

//! 与えられた二つの整数の最大公約数を求める．
/*!
  \param m	第1の整数
  \param n	第2の整数
  \return	mとnの最大公約数
*/
template <class S, class T> constexpr std::common_type_t<S, T>
gcd(S m, T n)
{
    return (n == 0 ? m : gcd(n, m % n));
}

//! 与えられた三つ以上の整数の最大公約数を求める．
/*!
  \param m	第1の整数
  \param n	第2の整数
  \param args	第3, 第4,...の整数
  \return	m, n, args...の最大公約数
*/
template <class S, class T, class... ARGS>
constexpr std::common_type_t<S, T, ARGS...>
gcd(S m, T n, ARGS... args)
{
    return gcd(gcd(m, n), args...);
}

//! 与えられた二つの整数の最小公倍数を求める．
/*!
  \param m	第1の整数
  \param n	第2の整数
  \return	mとnの最小公倍数
*/
template <class S, class T> constexpr std::common_type_t<S, T>
lcm(S m, T n)
{
    return (m*n == 0 ? 0 : (m / gcd(m, n)) * n);
}

//! 与えられた三つ以上の整数の最小公倍数を求める．
/*!
  \param m	第1の整数
  \param n	第2の整数
  \param args	第3, 第4,...の整数
  \return	m, n, args...の最小公倍数
*/
template <class S, class T, class... ARGS>
constexpr std::common_type_t<S, T, ARGS...>
lcm(S m, T n, ARGS... args)
{
    return lcm(lcm(m, n), args...);
}

/************************************************************************
*  for_each<N>(FUNC func, size_t n, ITER... iter)			*
************************************************************************/
namespace detail
{
  template <class FUNC, class... ITER> inline FUNC
  for_each(std::integral_constant<size_t, 0>,
	   FUNC func, size_t n, ITER... iter)
  {
      if (n)
      {
	  func(*iter...);
	  while (--n)
	      func(*++iter...);
      }
      return func;
  }
  template <class FUNC, class... ITER> inline FUNC
  for_each(std::integral_constant<size_t, 1>, FUNC func, size_t, ITER... iter)
  {
      func(*iter...);
      return func;
  }
  template <size_t N, class FUNC, class... ITER> inline FUNC
  for_each(std::integral_constant<size_t, N>, FUNC func, size_t n, ITER... iter)
  {
      func(*iter...);
      return for_each(std::integral_constant<size_t, N-1>(),
		      func, n, ++iter...);
  }
}	// namespace detail

//! 指定された範囲の各要素に関数を適用する
/*!
  N != 0 の場合，Nで指定した要素数だけ適用し，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param func	適用する関数
  \param n	適用要素数
  \param iter	適用範囲の先頭を指す反復子
*/
template <size_t N, class FUNC, class... ITER> inline FUNC
for_each(FUNC func, size_t n, ITER... iter)
{
    return detail::for_each(std::integral_constant<size_t, N>(),
			    func, n, iter...);
}
    
/************************************************************************
*  fill<N>(ITER iter, size_t n, const T& val)				*
************************************************************************/
//! 指定された範囲を与えられた値で埋める
/*!
  N != 0 の場合，Nで指定した要素数だけ埋め，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param iter	埋める範囲の先頭を指す反復子
  \param n	埋める要素数
  \param val	埋める値
*/
template <size_t N, class ITER, class T> inline void
fill(ITER iter, size_t n, const T& val)
{
    for_each<N>([&val](auto&& dst){ dst = val; }, n, iter);
}
    
/************************************************************************
*  copy<N>(IN in, size_t n, OUT out)					*
************************************************************************/
//! 指定された範囲をコピーする
/*!
  N != 0 の場合，Nで指定した要素数をコピーし，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param in	コピー元の先頭を指す反復子
  \param n	コピーする要素数
  \param out	コピー先の先頭を指す反復子
  \return	コピー先の末尾の次
*/
template <size_t N, class IN, class OUT> inline void
copy(IN in, size_t n, OUT out)
{
#ifdef TU_DEBUG
  //std::cout << "copy<" << N << "> ["
  //	    << print_sizes(range<IN, N>(in, n)) << ']' << std::endl;
#endif
    for_each<N>([](auto&& x, const auto& y){ x = y; }, n, out, in);
}

/************************************************************************
*  inner_product<N>(ITER0 iter0, size_t n, ITER1 iter1, T init)		*
************************************************************************/
template <class X, class Y, class Z>
inline Z	fma(X x, Y y, Z z)	{ return x*y + z; }
    
namespace detail
{
  template <class ITER0, class ITER1, class T> T
  inner_product(ITER0 iter0, size_t n, ITER1 iter1, T init,
		std::integral_constant<size_t, 0>)
  {
      for (; n--; ++iter0, ++iter1)
	  init = fma(*iter0, *iter1, init);
      return init;
  }
  template <class ITER0, class ITER1, class T> inline T
  inner_product(ITER0 iter0, size_t, ITER1 iter1, T init,
		std::integral_constant<size_t, 1>)
  {
      return fma(*iter0, *iter1, init);
  }
  template <class ITER0, class ITER1, class T, size_t N> inline T
  inner_product(ITER0 iter0, size_t n, ITER1 iter1, T init,
		std::integral_constant<size_t, N>)
  {
      return inner_product(std::next(iter0), n, std::next(iter1),
			   fma(*iter0, *iter1, init),
			   std::integral_constant<size_t, N-1>());
  }
}	// namespace detail

//! 指定された範囲の内積の値を返す
/*!
  N != 0 の場合，Nで指定した要素数の範囲の内積を求め，argは無視．
  N = 0 の場合，要素数をnで指定，
  \param iter0	適用範囲の第1変数の先頭を指す反復子
  \param n	要素数
  \param iter1	適用範囲の第2変数の先頭を指す反復子
  \param init	初期値
  \return	内積の値
*/
template <size_t N, class ITER0, class ITER1, class T> inline T
inner_product(ITER0 iter0, size_t n, ITER1 iter1, T init)
{
#ifdef TU_DEBUG
  //std::cout << "inner_product<" << N << "> ["
  //	      << print_sizes(range<ITER0, N>(iter0, n)) << ']' << std::endl;
#endif
    return detail::inner_product(iter0, n, iter1, init,
				 std::integral_constant<size_t, N>());
}
    
/************************************************************************
*  square<N>(ITER iter, size_t n)					*
************************************************************************/
namespace detail
{
  template <class T, std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
  inline auto
  square(T val)
  {
      return val * val;
  }
  template <class ITER> auto
  square(ITER iter, size_t n, std::integral_constant<size_t, 0>)
  {
      using value_type	= decltype(square(*iter));

      value_type	val = 0;
      for (; n--; ++iter)
	  val += square(*iter);
      return val;
  }
  template <class ITER> inline auto
  square(ITER iter, size_t n, std::integral_constant<size_t, 1>)
  {
      return square(*iter);
  }
  template <class ITER, size_t N> inline auto
  square(ITER iter, size_t n, std::integral_constant<size_t, N>)
  {
      const auto	tmp = square(*iter);
      return tmp + square(++iter, n, std::integral_constant<size_t, N-1>());
  }
}	// namespace detail

//! 指定された範囲にある要素の2乗和を返す
/*!
  N != 0 の場合，Nで指定した要素数の範囲の2乗和を求め，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param iter	適用範囲の先頭を指す反復子
  \param arg	適用範囲の末尾の次を指す反復子または要素数
  \return	2乗和の値
*/
template <size_t N, class ITER> inline auto
square(ITER iter, size_t n)
{
    return detail::square(iter, n, std::integral_constant<size_t, N>());
}

//! 与えられた数値の2乘値を返す
/*
  \param val	数値
  \retrun	2乘値
*/ 
template <class T, std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
inline auto
square(T val)
{
    return detail::square(val);
}
    
/************************************************************************
*  conj(T)								*
************************************************************************/
template <class T> inline std::enable_if_t<std::is_arithmetic<T>::value, T>
conj(T val)
{
    return val;
}
    
/************************************************************************
*  inverse(T)								*
************************************************************************/
template <class T> inline std::enable_if_t<std::is_floating_point<T>::value, T>
inverse(T val)
{
    return T(1)/val;
}

/************************************************************************
*  op3x3(ROW row, ROW rowe, OP op)					*
************************************************************************/
//! 2次元データに対して3x3ウィンドウを走査してin-place近傍演算を行う．
/*!
  \param row	最初の行を示す反復子
  \param rowe	最後の行の次を示す反復子
  \param op	3x3ウィンドウを定義域とする演算子
*/
template <class ROW, class OP> void
op3x3(ROW row, ROW rowe, OP op)
{
    --rowe;
    for (const auto row0 = row++; row != rowe; )
    {
	auto	p   = row0->begin();	// 最初の行を一つ前の行のバッファとして使用
	auto	q   = row->begin();
	auto	val = *q;
	auto	re  = (++row)->end();
	--re;
	--re;
	for (auto r = row->begin(); r != re; ++r)
	{
	    auto tmp = op(p, q, r);
	    *p = *q;			// 次行の左上画素 = 左画素をバッファに退避
	    *q = val;			// 左画素における結果を書き込む
	    val = tmp;			// 次ウィンドウの左画素における結果を保存
	    ++p;
	    ++q;
	}
	*p = *q;			// 次行の左上画素 = 左画素をバッファに退避
	*q = val;
	++p;
	++q;
	*p = *q;			// 次行の上画素 = 注目画素をバッファに退避
    }
}

/************************************************************************
*  morphological operations						*
************************************************************************/
//! 3x3ウィンドウ内の最大値を返す．
/*!
  \param p	注目点の左上点を指す反復子
  \param q	注目点の左の点を指す反復子
  \param r	注目点の左下点を指す反復子
  \return	3x3ウィンドウ内の最大値
*/
template <class COL> inline typename std::iterator_traits<COL>::value_type
max3x3(COL p, COL q, COL r)
{
    using namespace	std;
	    
    return max({max({*p, *(p + 1), *(p + 2)}),
		max({*q, *(q + 1), *(q + 2)}),
		max({*r, *(r + 1), *(r + 2)})});
}
    
//! 3x3ウィンドウ内の最小値を返す．
/*!
  \param p	注目点の左上点を指す反復子
  \param q	注目点の左の点を指す反復子
  \param r	注目点の左下点を指す反復子
  \return	3x3ウィンドウ内の最小値
*/
template <class COL> inline typename std::iterator_traits<COL>::value_type
min3x3(COL p, COL q, COL r)
{
    using namespace	std;
	    
    return min({min({*p, *(p + 1), *(p + 2)}),
		min({*q, *(q + 1), *(q + 2)}),
		min({*r, *(r + 1), *(r + 2)})});
}

//! morphological open演算をin-placeで行う．
/*!
  指定された回数だけ収縮(erosion)を行った後，同じ回数だけ膨張(dilation)を行う．
  \param row	最初の行を示す反復子
  \param rowe	最後の行の次を示す反復子
  \param niter	収縮と膨張の回数
*/
template <class ROW> void
mopOpen(ROW row, ROW rowe, size_t niter=1)
{
    using col_iterator
	= typename std::iterator_traits<ROW>::value_type::iterator;

    for (size_t n = 0; n < niter; ++n)
	op3x3(row, rowe, min3x3<col_iterator>);	// 収縮(erosion)
    for (size_t n = 0; n < niter; ++n)
	op3x3(row, rowe, max3x3<col_iterator>);	// 膨張(dilation)
}

//! morphological close演算をin-placeで行う．
/*!
  指定された回数だけ膨張(dilation)を行った後，同じ回数だけ収縮(erosion)を行う．
  \param row	最初の行を示す反復子
  \param rowe	最後の行の次を示す反復子
  \param niter	収縮と膨張の回数
*/
template <class ROW> void
mopClose(ROW row, ROW rowe, size_t niter=1)
{
    using col_iterator
	= typename std::iterator_traits<ROW>::value_type::iterator;
    
    for (size_t n = 0; n < niter; ++n)
	op3x3(row, rowe, max3x3<col_iterator>);	// 膨張(dilation)
    for (size_t n = 0; n < niter; ++n)
	op3x3(row, rowe, min3x3<col_iterator>);	// 収縮(erosion)
}
    
//! 先頭要素がそのまま出力されるscan演算を行う．
/*!
  \param in	入力の先頭を指す反復子
  \param ie	入力の末尾の次を指す反復子
  \param out	出力の先頭を指す反復子
  \param op	結合律を滿たす2変数演算子
  \return	出力の末尾の次を指す反復子
*/
template <class IN, class OUT, class OP=std::plus<> > OUT
inclusive_scan(IN in, IN ie, OUT out, OP op=OP())
{
    if (in != ie)
	for (*out = *in; ++in != ie; )
	{
	    const auto	prev = out;
	    *++out = op(*prev, *in);
	}

    return out;
}
    
}	// namespace TU
#endif	// !TU_ALGORITHM_H
