/*!
  \file		shift_iterator.h
  \author	Toshio UESHIBA
  \brief	SIMDベクトルの並びを1成分ずつ右シフトする反復子の定義
*/
#if !defined(TU_SIMD_SHIFT_ITERATOR_H)
#define TU_SIMD_SHIFT_ITERATOR_H

#include "TU/iterator.h"
#include "TU/simd/shift.h"

namespace TU
{
namespace simd
{
/************************************************************************
*  class shift_iterator<ITER>						*
************************************************************************/
template <class ITER>
class shift_iterator
    : public boost::iterator_adaptor<shift_iterator<ITER>,
				     ITER,
				     boost::use_default,
				     boost::forward_traversal_tag,
				     std::iter_value_t<ITER> >
{
  private:
    typedef boost::iterator_adaptor<shift_iterator,
				    ITER,
				    boost::use_default,
				    boost::forward_traversal_tag,
				    std::iter_value_t<ITER> >	super;

  public:
    typedef typename super::value_type		value_type;
    typedef typename super::reference		reference;

    friend class	boost::iterator_core_access;

  public:
		shift_iterator(ITER iter, size_t pos=0)
		    :super(iter), _pos(0), _val(*iter), _next(), _valid(true)
		{
		    while (pos--)
			increment();
		}

  private:
    void	shift() const
		{
		    _val  = shift_r<1>(_next, _val);
		    _next = shift_r<1>(_next);
		}
    void	load_and_shift() const
		{
		    _next  = *super::base();
		    _valid = true;
		    shift();
		}
    reference	dereference() const
		{
		    if (!_valid)		// !_valid なら必ず _pos == 1
			load_and_shift();
		    return _val;
		}
    void	increment()
		{
		    switch (++_pos)
		    {
		      case 1:
			++super::base_reference();
			_valid = false;
			break;
		      case value_type::size:
			_pos = 0;		// default:に落ちる
		      default:
			if (!_valid)		// !_valid なら必ず _pos == 2
			    load_and_shift();
			shift();
			break;
		    }
		}
    bool	equal(const shift_iterator& iter) const
		{
		    return (super::base() == iter.base()) &&
			   (_pos == iter._pos);
		}

  private:
    size_t		_pos;
    mutable value_type	_val, _next;
    mutable bool	_valid;	//!< _nextに入力値が読み込まれていればtrue
};

template <class ITER> shift_iterator<ITER>
make_shift_iterator(ITER iter)
{
    return shift_iterator<ITER>(iter);
}
    
}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_SHIFT_ITERATOR_H
