/*!
  \file		Array++.h
  \author	Toshio UESHIBA
  \brief	SIMD命令を適用できる配列クラスの定義
*/
#if !defined(TU_SIMD_ARRAYPP_H)
#define TU_SIMD_ARRAYPP_H

#include "TU/Array++.h"
#include "TU/simd/simd.h"

namespace TU
{
#if defined(SIMD)
/************************************************************************
*  make_array()								*
************************************************************************/
namespace simd
{
template <size_t SIZE, class ITER, bool ALIGNED> auto
make_array(std::iter_reference_t<ITER> x, iterator_wrapper<ITER, ALIGNED> iter)
    -> array<typename decltype(*make_accessor(iter))::element_type,
	     simd::allocator<
		 typename decltype(*make_accessor(iter))::element_type>,
	     SIZE>;

/************************************************************************
*  algorithms overloaded for simd::iterator_wrapper<ITER, ALIGNED>	*
************************************************************************/
namespace detail
{
  constexpr inline size_t
  nsteps(size_t n, size_t step)
  {
      return (n ? (n - 1)/step + 1 : 0);
  }

  template <class... STRIDE> constexpr inline std::tuple<STRIDE...>
  nsteps(std::tuple<STRIDE...> strides, size_t step)
  {
      return TU::tuple_transform([step](auto stride)
				 { return nsteps(stride, step); },
				 strides);
  }
}	// namespace detail

//! 指定された範囲の各要素に関数を適用する
/*!
  N != 0 の場合，Nで指定した要素数だけ適用し，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param func	適用する関数
  \param n	適用要素数
  \param iter0	第1変数の適用範囲の先頭を指す反復子
  \param iter	第2変数以降の適用範囲の先頭を指す反復子
*/
template <size_t N, class FUNC, class ITER0, bool ALIGNED0,
	  class... ITER, bool... ALIGNED> inline FUNC
for_each(FUNC func, size_t n, iterator_wrapper<ITER0, ALIGNED0> iter0,
	 iterator_wrapper<ITER, ALIGNED>... iter)
{
    auto		map  = make_accessor(make_map_iterator(func, iter0,
							       iter...));
    constexpr auto	STEP = map.step();
    
    TU::for_each<detail::nsteps(N, STEP)>(map, detail::nsteps(n, STEP));

    return func;
}

//! 指定された範囲の内積の値を返す
/*!
  N != 0 の場合，Nで指定した要素数の範囲の内積を求め，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param iter0	第1変数の適用範囲の先頭を指す反復子
  \param n	要素数
  \param iter1	第2変数の適用範囲の先頭を指す反復子
  \param init	初期値
  \return	内積の値
*/
template <size_t N, class ITER0, bool ALIGNED0,
		    class ITER1, bool ALIGNED1, class T> inline T
inner_product(iterator_wrapper<ITER0, ALIGNED0> iter0, size_t n,
	      iterator_wrapper<ITER1, ALIGNED1> iter1, T init)
{
#ifdef TU_DEBUG
    std::cout << "(simd)inner_product<" << N << "> ["
	      << print_sizes(range<iterator_wrapper<ITER0, ALIGNED0>, N>(
				 iter0, n))
	      << "]" << std::endl;
#endif
    constexpr auto	STEP = std::iter_value_t<decltype(make_accessor(iter0))>
				::size;
    
    return hadd(TU::inner_product<detail::nsteps(N, STEP)>(
		    make_accessor(iter0), detail::nsteps(n, STEP),
		    make_accessor(iter1), vec<T>(init)));
}
    
//! 指定された範囲にある要素の2乗和を返す
/*!
  N != 0 の場合，Nで指定した要素数の範囲の2乗和を求め，nは無視．
  N = 0 の場合，要素数をnで指定，
  \param begin	適用範囲の先頭を指す反復子
  \param n	要素数
  \return	2乗和の値
*/
template <size_t N, class ITER, bool ALIGNED> inline auto
square(iterator_wrapper<ITER, ALIGNED> iter, size_t n)
{
    constexpr auto	STEP = std::iter_value_t<decltype(make_accessor(iter))>
				::size;
    
    return hadd(square<detail::nsteps(N, STEP)>(make_accessor(iter),
						detail::nsteps(n, STEP)));
}

}	// namespace simd
	
namespace detail
{
  template <class T, bool ALIGNED>
  struct const_iterator_t<simd::iterator_wrapper<T*, ALIGNED> >
  {
      using type = simd::iterator_wrapper<const T*, ALIGNED>;
  };
}	// namespace detail
    
/************************************************************************
*  traits for Buf<T, ALLOC, SIZE, SIZES...>				*
************************************************************************/
template <class T, class ALLOC>
class BufTraits<simd::vec<T>, ALLOC>
    : public std::allocator_traits<simd::allocator<simd::vec<T> > >
{
  private:
    using super			= std::allocator_traits<
				      simd::allocator<simd::vec<T> > >;

  public:
    using iterator		= simd::store_iterator<T, true>;
    using const_iterator	= simd::load_iterator<T, true>;
    
  protected:
    using			typename super::pointer;

    constexpr static size_t	Alignment = super::allocator_type::Alignment;
    
    static pointer	null()		{ return nullptr; }
    static auto		ptr(pointer p)	{ return p.base(); }
};
#endif	// defined(SIMD)
}	// namespace TU
#endif	// !TU_SIMD_ARRAYPP_H
