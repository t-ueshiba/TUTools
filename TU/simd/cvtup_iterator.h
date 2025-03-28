/*!
  \file		cvtup_iterator.h
  \author	Toshio UESHIBA
  \brief	より大きな成分を持つSIMDベクトルへの型変換を行う反復子の定義
*/
#if !defined(TU_SIMD_CVTUP_ITERATOR_H)
#define TU_SIMD_CVTUP_ITERATOR_H

#include "TU/iterator.h"
#include "TU/simd/cvt.h"

namespace TU
{
namespace simd
{
/************************************************************************
*  class cvtup_iterator<ITER, MASK>					*
************************************************************************/
namespace detail
{
  template <class ITER, bool MASK>
  class cvtup_proxy
  {
    public:
      using value_type	= std::iter_value_t<ITER>;
      using self	= cvtup_proxy;

    private:
      using T		= typename tuple_head<value_type>::element_type;

    private:
      template <class VEC_, class OP_>
      std::enable_if_t<(tuple_head<VEC_>::size == vec<T>::size)>
		exec(const VEC_& x, OP_ op)
		{
		    op(*_iter++, simd::cvtup<T, false, MASK>(x));
		}
      template <class VEC_, class OP_>
      std::enable_if_t<(tuple_head<VEC_>::size > vec<T>::size)>
		exec(const VEC_& x, OP_ op)
		{
		    exec(simd::cvtup<T, false, MASK>(x), op);
		    exec(simd::cvtup<T, true,  MASK>(x), op);
		}

    public:
      cvtup_proxy(const ITER& iter) :_iter(const_cast<ITER&>(iter)) {}
	
      template <class VEC_>
      self&	operator =(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t = s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator +=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t += s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator -=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t -= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator *=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t *= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator /=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t /= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator %=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t %= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator &=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t &= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator |=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t |= s; });
		    return *this;
		}
      template <class VEC_>
      self&	operator ^=(const VEC_& x)
		{
		    exec(x, [](auto&& t, const auto& s){ t ^= s; });
		    return *this;
		}
	
    private:
      ITER&	_iter;
  };
}

//! SIMDベクトルを受け取ってより大きな成分を持つ複数のSIMDベクトルに変換し，それらを指定された反復子を介して書き込む反復子
/*!
  \param ITER	変換されたSIMDベクトルの書き込み先を指す反復子
*/
template <class ITER, bool MASK=false>
class cvtup_iterator
    : public boost::iterator_adaptor<
		 cvtup_iterator<ITER, MASK>,
		 ITER,
		 typename detail::cvtup_proxy<ITER, MASK>::value_type,
		 boost::single_pass_traversal_tag,
		 detail::cvtup_proxy<ITER, MASK> >
{
  private:
    using super	= boost::iterator_adaptor<
			cvtup_iterator,
			ITER,
			typename detail::cvtup_proxy<ITER, MASK>::value_type,
			boost::single_pass_traversal_tag,
			detail::cvtup_proxy<ITER, MASK> >;

  public:
    using typename super::difference_type;
    using typename super::value_type;
    using typename super::reference;

    friend class boost::iterator_core_access;

  public:
    cvtup_iterator(const ITER& iter)	:super(iter)			{}

  private:
    reference		dereference() const
			{
			    return reference(super::base());
			}
    void		advance(difference_type)			{}
    void		increment()					{}
    void		decrement()					{}
    difference_type	distance_to(const cvtup_iterator& iter) const
			{
			    return (iter.base() - super::base())
				 / value_type::size;
			}
};

template <class ITER> cvtup_iterator<ITER, false>
make_cvtup_iterator(ITER iter)
{
    return cvtup_iterator<ITER, false>(iter);
}

template <class ITER>
using cvtup_mask_iterator = cvtup_iterator<ITER, true>;
    
template <class ITER> cvtup_mask_iterator<ITER>
make_cvtup_mask_iterator(ITER iter)
{
    return cvtup_mask_iterator<ITER>(iter);
}

}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_CVTUP_ITERATOR_H
