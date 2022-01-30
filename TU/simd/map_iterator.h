/*!
  \file		map_iterator.h
  \author	Toshio UESHIBA
  \brief	SIMDベクトル間の型変換関数の定義
*/
#if !defined(TU_SIMD_MAP_ITERATOR_H)
#define	TU_SIMD_MAP_ITERATOR_H

#include "TU/range.h"
#include "TU/simd/cvt.h"
#include "TU/simd/iterator_wrapper.h"

namespace TU
{
namespace simd
{
/************************************************************************
*  class map_iterator<T, MASK, FUNC, ITERS>				*
************************************************************************/
namespace detail
{
  template <class T, class ITERS>
  using map_iterator_argument_element
		= typename std::conditional<
			std::is_void<T>::value,
			iterator_value<std::decay_t<tuple_head<ITERS> > >,
			vec<T> >::type::element_type;

  // ITERS(tuple かもしれない)の各要素を vec<T> に置き換えた型の引数を
  // FUNC に適用したときの結果の型
  template <class T, class FUNC, class ITERS>
  using map_iterator_result
		= decltype(
		    TU::apply(
			std::declval<FUNC>(),
			std::declval<replace_element<ITERS, vec<T> > >()));
}
    
//! 入力を適切に変換してから関数を適用し，結果をconvert downする．
/*!
  戻り値のSIMDベクトルは，vec<S>と入力のうち最下位のSIMDベクトルと同位
  \param T	FUNCの引数となるSIMDベクトルの要素型
  \param MASK
  \param FUNC	適用する関数
 */
template <class T, bool MASK, class FUNC, class ITERS>
class map_iterator
    : public boost::iterator_facade<
		map_iterator<T, MASK, FUNC, ITERS>,
		std::conditional_t<
		    std::is_void<
			detail::map_iterator_result<T, FUNC, ITERS> >::value,
		    vec<T>, detail::map_iterator_result<T, FUNC, ITERS> >, 
		boost::single_pass_traversal_tag,
		std::conditional_t<
		    std::is_void<
			detail::map_iterator_result<T, FUNC, ITERS> >::value,
		    vec<T>, detail::map_iterator_result<T, FUNC, ITERS> > > 
{
  private:
    friend	class boost::iterator_core_access;

    template <class ITER_, class=void>
    struct vsize_impl	// ITER_  が detail::store_proxy の場合にも有効
    {			// にするため，iterator_value<ITER_> は使わない
	constexpr static auto
		max = vsize_impl<typename ITER_::value_type>::max;
    };
    template <class T_>
    struct vsize_impl<vec<T_> >
    {
	constexpr static auto	max = vec<T_>::size;
    };
    template <class DUMMY_>
    struct vsize_impl<std::tuple<>, DUMMY_>
    {
	constexpr static size_t	max = 0;
    };
    template <class VEC_, class... VECS_>
    struct vsize_impl<std::tuple<VEC_, VECS_...> >
    {
	constexpr static auto
		max = std::max(vsize_impl<std::decay_t<VEC_>   >::max,
			       vsize_impl<std::tuple<VECS_...> >::max);
    };
    template <class T_, bool MASK_, class FUNC_, class ITERS_>
    struct vsize_impl<map_iterator<T_, MASK_, FUNC_, ITERS_> >
    {
	constexpr static auto
		max = std::max(vsize_impl<std::decay_t<ITERS_> >::max,
			       vec<T_>::size);
    };

  //! tuple中のベクトルおよび反復子が指すベクトルのうちの最大/最小要素数
    template <class VECS_>
    using vsize = vsize_impl<std::decay_t<VECS_> >;

  //! 入力の型と演算の型がそれぞれ ITERS_, T_ である map_iterator から T_ 型の出力を
    template <class ITERS_, class S_>
    struct result
    {
	constexpr static auto
		size = std::max(vec<S_>::size,
				std::min(vsize<ITERS_>::max, vec<T>::size));
    };

  private:
  /************************ convert up ************************/
  // SIMDベクトルを convert up
    template <bool, size_t N_, class S_>
    static std::enable_if_t<vec<S_>::size == N_, vec<S_> >
		cvtup(vec<S_> x)	// 与えられた vec のサイズが
		{			// 要求値 _N に等しければ
		    return x;		// vec をそのまま返す
		}
    template <bool HI_, size_t N_, class S_,
	      std::enable_if_t<vec<S_>::size != N_>* = nullptr>
    static auto	cvtup(vec<S_> x)	// 与えられた vec のサイズが
		{			// 要求値 _N でなければ直上の vec を返す
		    return cvt<cvt_upper_type<T, S_, MASK>, HI_, MASK>(x);
		}

  // 反復子を convert up
    template <bool, size_t N_, class ITER_,
	      std::enable_if_t<iterator_value<ITER_>::size == N_>* = nullptr>
    static auto	cvtup(ITER_& iter)	// 反復子が指す vec のサイズが
		{			// 要求値 _N に等しければ
		    return *iter++;	// vec を返してインクリメント
		}
    template <bool, size_t N_, class ITER_>
    static std::enable_if_t<iterator_value<ITER_>::size != N_, ITER_&>
		cvtup(ITER_& iter)	// 反復子が指す vec のサイズが
		{			// 要求値 _N でなければ
		    return iter;	// 反復子そのものへの参照を返す
		}

  // map_iterator を convert up
    template <bool HI_, size_t N_,
	      class S_, bool MASK_, class FUNC_, class ITERS_>
    static auto	cvtup(map_iterator<S_, MASK_, FUNC_, ITERS_>& m)
		{
		    return m.template cvtup<HI_, N_>();
		}

  // SIMDベクトル|反復子|map_iterator の tuple を convert up
    template <bool HI_, size_t N_, class... VEC_>
    static auto	cvtup(std::tuple<VEC_...>& t)
		{
		    return tuple_transform(
				[](auto&& x) -> decltype(auto)
				{ return map_iterator::cvtup<HI_, N_>(x); },
				t);
		}

  /************************ convert down ************************/
  // SIMDベクトルを convert down
    template <size_t, class T_=T, class S_>
    static auto	cvtdown(vec<S_> x)		// サイズが不変な場合のみ変換
		{
		    using L = cvt_lower_type<T_, S_, MASK>;
		    using B = std::conditional_t<vec<L>::size == vec<S_>::size,
						 L, S_>;

		    return cvt<B, false, MASK>(x);
		}
    template <class T_, class S_>
    static auto	cvtdown(vec<S_> x, vec<S_> y)
		{
		    return cvt<cvt_lower_type<T_, S_, MASK>, MASK>(x, y);
		}

  // 反復子を convert down
    template <size_t N_, class ITER_,
	      std::enable_if_t<iterator_value<ITER_>::size == N_>* = nullptr>
    static auto	cvtdown(ITER_& iter)		  // 反復子が指す vec のサイズが
		{				  // 要求値 _N に等しければ
		    return cvtdown<N_>(*iter++);  // vec を返してインクリメント
		}
    template <size_t N_, class ITER_,
	      std::enable_if_t<(iterator_value<ITER_>::size < N_)>* = nullptr>
    static auto	cvtdown(ITER_& iter)		// 反復子が指す vec のサイズが
		{				// 要求値 _N よりも小さければ
		    const auto	x = cvtdown<N_/2>(iter);  // 要求値を半分にして
		    const auto	y = cvtdown<N_/2>(iter);  // 2回dereference
    		    return cvtdown<T>(x, y);		  // した結果を融合する
		}

  // map_iteratorを convert down
    template <size_t N_, class S_, bool MASK_, class FUNC_, class ITERS_,
	      std::enable_if_t<result<ITERS_, S_>::size == N_>* = nullptr>
    static auto	cvtdown(map_iterator<S_, MASK_, FUNC_, ITERS_>& m)
		{
		    return m.template eval<T>();
		}
    template <size_t N_, class S_, bool MASK_, class FUNC_, class ITERS_,
	      std::enable_if_t<(result<ITERS_, S_>::size < N_)>* = nullptr>
    static auto	cvtdown(map_iterator<S_, MASK_, FUNC_, ITERS_>& m)
		{
		    const auto	x = cvtdown<N_/2>(m);
		    const auto	y = cvtdown<N_/2>(m);
		    return cvtdown<T>(x, y);
		}

  // detail::store_proxy を convert down
    template <size_t, bool ALIGNED_>
    static const auto&
		cvtdown(const detail::store_proxy<T, ALIGNED_>& x)
		{
		    return x;
		}

  // SIMDベクトル|反復子|map_iterator|detail::store_proxy の tuple を convert down
    template <size_t N_, class... VEC_>
    static auto	cvtdown(std::tuple<VEC_...>& t)
		{
		    return tuple_transform(
				[](auto&& x) -> decltype(auto)
				{ return map_iterator::cvtdown<N_>(x); },
				t);
		}
  
  // ITERS_ の全要素をT型に convert up/down して FUNC を適用
    template <class ITERS_>
    std::enable_if_t<(vsize<ITERS_>::max <= vec<T>::size)>
		exec(ITERS_&& iters) const
		{
		    TU::apply(_func, cvtdown<vec<T>::size>(iters));
		}
    template <class ITERS_>
    std::enable_if_t<(vsize<ITERS_>::max > vec<T>::size)>
		exec(ITERS_&& iters) const
		{
		    constexpr auto	N = vsize<ITERS_>::max;
		    
		    exec(cvtup<false, N/2>(iters));
		    exec(cvtup<true,  N/2>(iters));
		}

  // ITERS_ の全要素をT型に convert up/down して FUNC を適用し，
  // 結果を 最初の ITERS_ のレベル に convert down して返す
    template <class T_, class ITERS_,
	      std::enable_if_t<(vsize<ITERS_>::max <= vec<T>::size)>* = nullptr>
    auto	eval(ITERS_&& iters) const
		{
		    constexpr auto	N = vec<T>::size;

		  // 戻り値は T と同位
		    return cvtdown<N, T_>(TU::apply(_func, cvtdown<N>(iters)));
		}
    template <class T_, class ITERS_,
	      std::enable_if_t<(vsize<ITERS_>::max > vec<T>::size)>* = nullptr>
    auto	eval(ITERS_&& iters) const
		{
		    constexpr auto	N = vsize<ITERS_>::max;

		    const auto	x = eval<T_>(cvtup<false, N/2>(iters));
		    const auto	y = eval<T_>(cvtup<true,  N/2>(iters));
		  // 戻り値は ITERS_ の最低位要素と同位
		    return cvtdown<T_>(x, y);
		}

    auto	dereference() const
		{
		    return const_cast<map_iterator*>(this)->eval<T>();
		}

    void	increment()
		{
		}

    bool	equal(const map_iterator& iter) const
		{
		    return _iters == iter._iters;
		}

  public:
		map_iterator(FUNC&& func, ITERS&& iters)
		    :_iters(std::forward<ITERS>(iters)),
		     _func(std::forward<FUNC>(func))		{}

    constexpr
    static auto	step()			{ return vsize<map_iterator>::max; }
    const auto&	functor()	 const	{ return _func; }
    const auto&	get_iterator_tuple() const
    		{
		    return _iters;
		}
	
    void	operator ()()
		{
		    exec(cvtup<false, step()>(_iters));
		}

    template <class T_>
    auto	eval()
		{
		    return eval<T_>(cvtup<false, step()>(_iters));
		}

    template <bool HI_, size_t N_,
	      std::enable_if_t<vec<T>::size == N_>* = nullptr>
    auto	cvtup()
		{
		    return eval<T>(cvtup<HI_, N_>(_iters));
		}
    template <bool HI_, size_t N_,
	      std::enable_if_t<vec<T>::size != N_>* = nullptr>
    auto	cvtup()
		{
		    using iters_t = decltype(cvtup<HI_, N_>(_iters));
		    
		    return map_iterator<T, MASK, FUNC&, iters_t>(
				_func, cvtup<HI_, N_>(_iters));
		}
    
  private:
    ITERS	_iters;	// 反復子または反復子の tuple
    FUNC	_func;
};

template <class T=void, bool MASK=false, class FUNC,
	  class ITER0, bool ALIGNED0, class... ITER, bool... ALIGNED>
inline auto
make_map_iterator(FUNC&& func,
		  const iterator_wrapper<ITER0, ALIGNED0>& iter0,
		  const iterator_wrapper<ITER,  ALIGNED>&... iter)
{
    using iters_t  = decltype(make_accessor(
				  TU::make_zip_iterator(iter0, iter...)));
    using argelm_t = detail::map_iterator_argument_element<T, iters_t>;

    return wrap_iterator(map_iterator<argelm_t, MASK, FUNC, iters_t>(
			     std::forward<FUNC>(func),
			     make_accessor(
				 TU::make_zip_iterator(iter0, iter...))));
}

template <class T, bool MASK, class FUNC, class ITER, bool ALIGNED,
	  std::enable_if_t<std::is_void<iterator_iterator<ITER> >::value>*
	  = nullptr>
inline auto
make_map_iterator(TU::detail::mapped_tag<T, MASK, FUNC>&& m,
		  const iterator_wrapper<ITER, ALIGNED>& iter)
{
    using S = typename TU::detail::mapped_tag<T, MASK, FUNC>::element_type;
    
    return make_map_iterator<S, MASK>(m.functor(), iter);
}

/*
 *  iterator_stride<ITER> = decltype(stride(std::declval<ITER>())) において
 *    ITER = range_iterator<
 *		simd::iterator_wrapper<simd::map_iterator<...>, true>,
 *		STRIDE, SIZE>
 *  であるとき， simd::map_iterator型である iter.base().base() に対して
 *  ADLを利用して stride() を呼び出すために定義．
 */
template <class T, bool MASK, class FUNC, class ITERS> auto
stride(const map_iterator<T, MASK, FUNC, ITERS>& iter)
    -> decltype(TU::stride(iter.get_iterator_tuple()))
{
    return TU::stride(iter.get_iterator_tuple());
}

}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_MAP_ITERATOR_H
