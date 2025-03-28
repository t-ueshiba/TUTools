/*!
  \file		assignment_iterator.h
  \author	Toshio UESHIBA
  \brief	SIMDベクトル間の型変換関数の定義
*/
#if !defined(TU_SIMD_ASSIGNMENT_ITERATOR_H)
#define	TU_SIMD_ASSIGNMENT_ITERATOR_H

#include "TU/simd/cvt.h"

namespace TU
{
namespace simd
{
/************************************************************************
*  class assignment_iterator<ARG, MASK, FUNC, ITER...>			*
************************************************************************/
namespace detail
{
  template <class ARG, bool MASK, class FUNC, class ITER, class... ITERS>
  class assignment_proxy
  {
    public:
      using iterator	= std::conditional_t<
				sizeof...(ITERS),
				zip_iterator<std::tuple<ITER, ITERS...> >,
				ITER>;
      
    private:
      template <class OP_, class IN_>
      std::enable_if_t<(vsize<IN_>::max <= vec<ARG>::size)>
		exec(OP_ op, IN_&& in)
		{
		    op(*_iter, apply(_func, cvtdown<ARG, MASK>(in)));
		    ++_iter;
		}
      template <class OP_, class IN_>
      std::enable_if_t<(vsize<IN_>::max > vec<ARG>::size)>
		exec(OP_ op, IN_&& in)
		{
		    constexpr auto	N = vsize<IN_>::max;
		    
		    exec(op, cvtup<ARG, false, MASK, N/2>(in));
		    exec(op, cvtup<ARG, true,  MASK, N/2>(in));
		}
    
    public:
		assignment_proxy(iterator& iter, const FUNC& func)
		    :_iter(iter), _func(func)				{}

      template <class IN_>
      auto&	operator =(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y = x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator +=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y += x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator -=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y -= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator *=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y *= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator /=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y /= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator %=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y %= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator &=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y &= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator |=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y |= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}
      template <class IN_>
      auto&	operator ^=(IN_&& in)
		{
		    constexpr auto	N = std::max(vsize<IN_>::max,
						     vec<ARG>::size);
		    
		    exec([](auto&& y, const auto& x){ y ^= x; },
			 cvtup<ARG, false, MASK, N>(in));
		    return *this;
		}

    private:
      iterator&		_iter;
      const FUNC&	_func;
  };
}
    
//! 入力を適切に変換してから関数を適用し，結果をconvert downする．
/*!
  戻り値のSIMDベクトルは，vec<ARG>と入力のうち最下位のSIMDベクトルと同位
  \param FUNC	適用する関数
  \param ARG	FUNCの引数となるSIMDベクトルの要素型
  \param T	FUNCの結果のconvert down先のSIMDベクトルの要素型
  \param MASK
 */
template <class ARG, bool MASK, class FUNC, class... ITER>
class assignment_iterator
    : public boost::iterator_adaptor<
		assignment_iterator<ARG, MASK, FUNC, ITER...>,
		typename detail::assignment_proxy<ARG, MASK, FUNC, ITER...>
			       ::iterator,
		std::iter_value_t<
		    typename detail::assignment_proxy<ARG, MASK, FUNC, ITER...>
				   ::iterator>,
		boost::single_pass_traversal_tag,
		detail::assignment_proxy<ARG, MASK, FUNC, ITER...> >
{
  private:
    using proxy	= detail::assignment_proxy<ARG, MASK, FUNC, ITER...>;
    using super	= boost::iterator_adaptor<
			assignment_iterator,
			typename proxy::iterator,
			std::iter_value_t<typename proxy::iterator>,
			boost::single_pass_traversal_tag,
			proxy>;
    friend	class boost::iterator_core_access;

  public:
    using	typename super::reference;
    
  public:
    assignment_iterator(const FUNC& func, const ITER&... iter)
	:assignment_iterator(
	    std::integral_constant<bool, (sizeof...(ITER) > 1)>(),
	    func, iter...)						{}

    const auto&	functor()	const	{ return _func; }

  private:
    assignment_iterator(std::true_type,  const FUNC& func, const ITER&... iter)
	:super(std::make_tuple(iter...)), _func(func)			{}
    assignment_iterator(std::false_type, const FUNC& func, const ITER&... iter)
	:super(iter...), _func(func)					{}

    reference	dereference() const
		{
		    return {const_cast<assignment_iterator*>(this)
				->base_reference(), _func};
		}
    void	increment()						{}
    
  private:
    FUNC 	_func;	// 代入を可能にするためconstは付けない
};

template <class ARG=void, bool MASK=false, class FUNC, class... ITER>
inline assignment_iterator<ARG, MASK, FUNC, ITER...>
make_assignment_iterator(const FUNC& func, const ITER&... iter)
{
    return {func, iter...};
}

}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_ASSIGNMENT_ITERATOR_H
