/*!
  \file		cvtdown_iterator.h
  \author	Toshio UESHIBA
  \brief	より小さな成分を持つSIMDベクトルへの型変換を行う反復子の定義
*/
#if !defined(TU_SIMD_CVTDOWN_ITERATOR_H)
#define TU_SIMD_CVTDOWN_ITERATOR_H

#include "TU/iterator.h"
#include "TU/functional.h"
#include "TU/simd/cvt.h"

namespace TU
{
namespace simd
{
/************************************************************************
*  class cvtdown_iterator<T, ITER, MASK>				*
************************************************************************/
//! SIMDベクトルを出力する反復子を介して複数のSIMDベクトルを読み込み，それをより小さな成分を持つSIMDベクトルに変換する反復子
/*!
  \param T	変換先のSIMDベクトルの成分の型
  \param ITER	SIMDベクトルを出力する反復子
*/
template <class T, class ITER, bool MASK=false>
class cvtdown_iterator
    : public boost::iterator_adaptor<
		 cvtdown_iterator<T, ITER, MASK>,
		 ITER,
		 replace_element<std::iter_value_t<ITER>, vec<T> >,
		 boost::single_pass_traversal_tag,
		 replace_element<std::iter_value_t<ITER>, vec<T> > >
{
  private:
    using src_type	= std::iter_value_t<ITER>;
    using super		= boost::iterator_adaptor<
				cvtdown_iterator,
				ITER,
				replace_element<src_type, vec<T> >,
				boost::single_pass_traversal_tag,
				replace_element<src_type, vec<T> > >;

  // 変換元のベクトルの成分数
    constexpr static size_t	M = tuple_head<src_type>::size;
    
  public:
    using	typename super::difference_type;
    using	typename super::reference;

    friend	class boost::iterator_core_access;

  public:
		cvtdown_iterator(const ITER& iter)	:super(iter)	{}

  private:
    reference	dereference() const
		{
		    return cvtdown<T, MASK>(const_cast<cvtdown_iterator*>(this)
					    ->base_reference());
		}
    void	advance(difference_type)				{}
    void	increment()						{}
    void	decrement()						{}
};
    
template <class T, class ITER> cvtdown_iterator<T, ITER, false>
make_cvtdown_iterator(ITER iter)
{
    return cvtdown_iterator<T, ITER, false>(iter);
}

template <class T, class ITER>
using cvtdown_mask_iterator = cvtdown_iterator<T, ITER, true>;
    
template <class T, class ITER> cvtdown_mask_iterator<T, ITER>
make_cvtdown_mask_iterator(ITER iter)
{
    return cvtdown_mask_iterator<T, ITER>(iter);
}
    
}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_CVTDOWN_ITERATOR_H
