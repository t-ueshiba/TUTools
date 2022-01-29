/*!
  \file		iterator.h
  \author	Toshio UESHIBA
  \brief	各種反復子の定義と実装
*/
#ifndef TU_ITERATOR_H
#define TU_ITERATOR_H

#include <iterator>
#include <boost/iterator/iterator_adaptor.hpp>
#include "TU/tuple.h"

namespace TU
{
/************************************************************************
*  type aliases								*
************************************************************************/
//! 反復子が指す型
template <class ITER>
using iterator_value	  = typename std::iterator_traits<ITER>::value_type;

//! 反復子が指す型への参照
template <class ITER>
using iterator_reference  = typename std::iterator_traits<ITER>::reference;

//! 反復子が指す型へのポインタ
template <class ITER>
using iterator_pointer	  = typename std::iterator_traits<ITER>::pointer;
    
//! 2つの反復子間の差を表す型
template <class ITER>
using iterator_difference = typename std::iterator_traits<ITER>
					::difference_type;
//! 反復子のカテゴリ
template <class ITER>
using iterator_category	  = typename std::iterator_traits<ITER>
					::iterator_category;

/************************************************************************
*  TU::[begin|end|rbegin|rend](T&&)					*
************************************************************************/
/*
 *  range<ITER, SIZE> 等のproxyオブジェクトの右辺値から
 *  非constな反復子を取り出すために定義
 */ 
template <class T> inline auto
begin(T&& x) -> decltype(std::begin(x))
{
    return std::begin(x);
}
    
template <class T> inline auto
end(T&& x) -> decltype(std::end(x))
{
    return std::end(x);
}
    
template <class T> inline auto
rbegin(T&& x) -> decltype(std::rbegin(x))
{
    return std::rbegin(x);
}
    
template <class T> inline auto
rend(T&& x) -> decltype(std::rend(x))
{
    return std::rend(x);
}

/************************************************************************
*  map_iterator<FUNC, ITER>						*
************************************************************************/
template <class FUNC, class ITER>
class map_iterator
    : public boost::iterator_adaptor<map_iterator<FUNC, ITER>,
	ITER,
	std::decay_t<
	    decltype(TU::apply(std::declval<FUNC>(),
			       std::declval<iterator_reference<ITER> >()))>,
	boost::use_default,
	decltype(TU::apply(std::declval<FUNC>(),
			   std::declval<iterator_reference<ITER> >()))>
{
  private:
    using ref	= decltype(TU::apply(std::declval<FUNC>(),
				     std::declval<iterator_reference<ITER> >()));
    using super	= boost::iterator_adaptor<map_iterator,
					  ITER,
					  std::decay_t<ref>,
					  boost::use_default,
					  ref>;
    friend	class boost::iterator_core_access;

  public:
    using	typename super::reference;
	
  public:
		map_iterator(FUNC&& func, const ITER& iter)
		    :super(iter), _func(std::forward<FUNC>(func))	{}

    const auto&	functor()	const	{ return _func; }
	
  private:
    reference	dereference()	const	{ return TU::apply(_func,
							   *super::base()); }
	
  private:
    FUNC	_func;	//!< 演算子
};

template <class T=void, class FUNC, class... ITER> inline auto
make_map_iterator(FUNC&& func, const ITER&... iter)
{
    using iters_t = decltype(TU::make_zip_iterator(iter...));
    
    return map_iterator<FUNC, iters_t>(std::forward<FUNC>(func),
				       TU::make_zip_iterator(iter...));
}

/************************************************************************
*  make_mbr_iterator<ITER, T>						*
************************************************************************/
//! T型のメンバ変数を持つオブジェクトへの反復子からそのメンバに直接アクセスする反復子を作る．
/*!
  \param iter	ベースとなる反復子
  \param mbr	iterが指すオブジェクトのメンバへのポインタ
*/
template <class ITER, class T> inline auto
make_mbr_iterator(const ITER& iter, T iterator_value<ITER>::* mbr)
{
    return make_map_iterator([mbr](auto&& x){ return x.*mbr; }, iter);
}

//! std::pairへの反復子からその第1要素に直接アクセスする反復子を作る．
/*!
  \param iter	ベースとなる反復子
*/
template <class ITER> inline auto
make_first_iterator(const ITER& iter)
{
    return make_mbr_iterator(iter, &iterator_value<ITER>::first);
}
    
//! std::pairへの反復子からその第2要素に直接アクセスする反復子を作る．
/*!
  \param iter	ベースとなる反復子
*/
template <class ITER> inline auto
make_second_iterator(const ITER& iter)
{
    return make_mbr_iterator(iter, &iterator_value<ITER>::second);
}
    
/************************************************************************
*  class assignment_iterator<FUNC, ITER>				*
************************************************************************/
//! libTUTools++ のクラスや関数の実装の詳細を収める名前空間
namespace detail
{
  template <class FUNC, class ITER>
  class assignment_proxy
  {
    private:
      template <class T_>
      static auto	check_func(ITER iter, const T_& val, FUNC func)
			    -> decltype(func(*iter, val), std::true_type());
      template <class T_>
      static auto	check_func(ITER iter, const T_& val, FUNC func)
			    -> decltype(*iter = func(val), std::false_type());
      template <class T_>
      using is_binary_func	= decltype(check_func(std::declval<ITER>(),
						      std::declval<T_>(),
						      std::declval<FUNC>()));
      
    public:
      assignment_proxy(const ITER& iter, const FUNC& func)
	  :_iter(iter), _func(func)					{}

      template <class T_>
      std::enable_if_t<is_binary_func<T_>::value, const assignment_proxy&>
			operator =(T_&& val) const
			{
			    _func(*_iter, std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      std::enable_if_t<!is_binary_func<T_>::value, const assignment_proxy&>
			operator =(T_&& val) const
			{
			    *_iter  = _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator +=(T_&& val) const
			{
			    *_iter += _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator -=(T_&& val) const
			{
			    *_iter -= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator *=(T_&& val) const
			{
			    *_iter *= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator /=(T_&& val) const
			{
			    *_iter /= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator %=(T_&& val) const
			{
			    *_iter %= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator &=(T_&& val) const
			{
			    *_iter &= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator |=(T_&& val) const
			{
			    *_iter |= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      const auto&	operator ^=(T_&& val) const
			{
			    *_iter ^= _func(std::forward<T_>(val));
			    return *this;
			}

    private:
      const ITER&	_iter;
      const FUNC&	_func;
  };
}

//! operator *()を左辺値として使うときに，この左辺値と右辺値に指定された関数を適用するための反復子
/*!
  \param FUNC	変換を行う関数オブジェクトの型
  \param ITER	変換結果の代入先を指す反復子
*/
template <class FUNC, class ITER>
class assignment_iterator
    : public boost::iterator_adaptor<assignment_iterator<FUNC, ITER>,
				     ITER,
				     iterator_value<ITER>,
				     iterator_category<ITER>,
				     detail::assignment_proxy<FUNC, ITER> >
{
  private:
    using super	= boost::iterator_adaptor<
				assignment_iterator,
				ITER,
				iterator_value<ITER>,
				iterator_category<ITER>,
				detail::assignment_proxy<FUNC, ITER> >;
    friend	class boost::iterator_core_access;

  public:
    using	typename super::reference;
    
  public:
		assignment_iterator(FUNC&& func, const ITER& iter)
		    :super(iter), _func(std::forward<FUNC>(func))	{}

    const auto&	functor()	const	{ return _func; }

  private:
    reference	dereference()	const	{ return {super::base(), _func}; }
    
  private:
    FUNC 	_func;	// 代入を可能にするためconstは付けない
};
    
template <class FUNC, class... ITERS> inline auto
make_assignment_iterator(FUNC&& func, const ITERS&... iters)
{
    using iters_t = decltype(TU::make_zip_iterator(iters...));
    
    return assignment_iterator<FUNC, iters_t>(std::forward<FUNC>(func),
					      TU::make_zip_iterator(iters...));
}

/************************************************************************
*  class row2col<ROW>							*
************************************************************************/
//! 行への参照を与えられると予め指定された列indexに対応する要素への参照を返す関数オブジェクト
/*!
  \param ROW	行を指す反復子
*/ 
template <class ROW>
class row2col
{
  public:
    using argument_type	= iterator_reference<ROW>;
    
  public:
			row2col(size_t col)	:_col(col)		{}
    
    decltype(auto)	operator ()(argument_type row) const
			{
			    return *(begin(row) + _col);
			}
    
  private:
    const size_t	_col;	//!< 列を指定するindex
};

/************************************************************************
*  alias vertical_iterator<ROW>						*
************************************************************************/
template <class ROW>
using vertical_iterator = map_iterator<row2col<ROW>, ROW>;

template <class ROW> inline vertical_iterator<ROW>
make_vertical_iterator(const ROW& row, size_t col)
{
    return {{col}, row};
}

/************************************************************************
*  class ring_iterator<ITER>						*
************************************************************************/
//! 2つの反復子によって指定された範囲を循環バッファとしてアクセスする反復子
/*!
  \param ITER	データ列中の要素を指す反復子の型
*/
template <class ITER>
class ring_iterator : public boost::iterator_adaptor<ring_iterator<ITER>, ITER>
{
  private:
    using super	= boost::iterator_adaptor<ring_iterator, ITER>;
    friend	class boost::iterator_core_access;
    template <class>
    friend	class ring_iterator;
    
  public:
    using	typename super::difference_type;
    
  public:
    ring_iterator()
	:super(), _begin(super::base()), _end(super::base()), _d(0)	{}
    
    ring_iterator(ITER begin, ITER end)
	:super(begin),
	 _begin(begin), _end(end), _d(std::distance(_begin, _end))	{}
    template <class ITER_>
    ring_iterator(const ring_iterator<ITER_>& iter)
	:super(iter.base()),
	 _begin(iter._begin), _end(iter._end), _d(iter._d)		{}

    difference_type	position() const
			{
			    return std::distance(_begin, super::base());
			}
    
  private:
    void		advance(difference_type n)
			{
			    n %= _d;
			    difference_type	i = position() + n;
			    if (i >= _d)
				std::advance(super::base_reference(), n - _d);
			    else if (i < 0)
				std::advance(super::base_reference(), n + _d);
			    else
				std::advance(super::base_reference(), n);
			}
    
    void		increment()
			{
			    if (++super::base_reference() == _end)
				super::base_reference() = _begin;
			}

    void		decrement()
			{
			    if (super::base() == _begin)
				super::base_reference() = _end;
			    --super::base_reference();
			}

    difference_type	distance_to(const ring_iterator& iter) const
			{
			    difference_type	n = iter.base() - super::base();
			    return (n > 0 ? n - _d : n);
			}

    template <class ITER_>
    bool		equal(const ring_iterator<ITER_>& iter) const
			{
			    return super::base() == iter.base();
			}
    
  private:
    ITER		_begin;	// 代入を可能にするためconstは付けない
    ITER		_end;	// 同上
    difference_type	_d;	// 同上
};

template <class ITER> ring_iterator<ITER>
make_ring_iterator(ITER begin, ITER end)	{ return {begin, end}; }

}	// namespace TU
#endif	// !TU_ITERATOR_H
