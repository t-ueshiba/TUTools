/*!
  \file		iterator.h
  \author	Toshio UESHIBA
  \brief	各種反復子の定義と実装
*/
#ifndef TU_ITERATOR_H
#define TU_ITERATOR_H

#include <iterator>
#include <functional>			// for std::function
#include <boost/iterator/iterator_adaptor.hpp>
#include "TU/tuple.h"

namespace std
{
#if __cplusplus < 201402L    
/************************************************************************
*  std::[rbegin|rend|cbegin|cend|crbegin|crend](T)			*
************************************************************************/
template <class T> inline auto
rbegin(const T& x) -> decltype(x.rbegin())
{
    return x.rbegin();
}
    
template <class T> inline auto
rbegin(T& x) -> decltype(x.rbegin())
{
    return x.rbegin();
}
    
template <class T> inline auto
rend(const T& x) -> decltype(x.rend())
{
    return x.rend();
}
    
template <class T> inline auto
rend(T& x) -> decltype(x.rend())
{
    return x.rend();
}
    
template <class T> inline auto
cbegin(const T& x) -> decltype(std::begin(x))
{
    return std::begin(x);
}
    
template <class T> inline auto
cend(const T& x) -> decltype(std::end(x))
{
    return std::end(x);
}

template <class T> inline auto
crbegin(const T& x) -> decltype(std::rbegin(x))
{
    return std::rbegin(x);
}
    
template <class T> inline auto
crend(const T& x) -> decltype(std::rend(x))
{
    return std::rend(x);
}

template <class ITER> inline auto
make_reverse_iterator(ITER iter)
{
    return reverse_iterator<ITER>(iter);
}
#endif
}	// namespace std

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
*  class zip_iterator<ITER_TUPLE>					*
************************************************************************/
namespace detail
{
  struct generic_dereference
  {
    // assignment_iterator<FUNC, ITER...> のように dereference すると
    // その base iterator への参照を内包する proxy を返す反復子もあるので，
    // 引数は const ITER_& 型にする．もしも ITER_ 型にすると，呼出側から
    // コピーされたローカルな反復子 iter への参照を内包する proxy を
    // 返してしまい，dangling reference が生じる．
      template <class ITER_>
      decltype(auto)	operator ()(const ITER_& iter) const
			{
			    return *iter;
			}
  };
}	// namespace detail
    
template <class ITER_TUPLE>
class zip_iterator : public boost::iterator_facade<
			zip_iterator<ITER_TUPLE>,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>())),
			iterator_category<std::tuple_element_t<0, ITER_TUPLE> >,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>()))>
{
  private:
    using super = boost::iterator_facade<
			zip_iterator,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>())),
			iterator_category<std::tuple_element_t<0, ITER_TUPLE> >,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>()))>;
    friend	class boost::iterator_core_access;
    
  public:
    using	typename super::reference;
    using	typename super::difference_type;
    
  public:
		zip_iterator(ITER_TUPLE iter_tuple)
		    :_iter_tuple(iter_tuple)				{}
    template <class ITER_TUPLE_,
	      std::enable_if_t<
		  std::is_convertible<ITER_TUPLE_, ITER_TUPLE>::value>*
	      = nullptr>
		zip_iterator(const zip_iterator<ITER_TUPLE_>& iter)
		    :_iter_tuple(iter.get_iterator_tuple())		{}

    const auto&	get_iterator_tuple()	const	{ return _iter_tuple; }
    
  private:
    reference	dereference() const
		{
		    return tuple_transform(detail::generic_dereference(),
					   _iter_tuple);
		}
    template <class ITER_TUPLE_>
    std::enable_if_t<std::is_convertible<ITER_TUPLE_, ITER_TUPLE>::value, bool>
		equal(const zip_iterator<ITER_TUPLE_>& iter) const
		{
		    return std::get<0>(iter.get_iterator_tuple())
			== std::get<0>(_iter_tuple);
		}
    void	increment()
		{
		    ++_iter_tuple;
		}
    void	decrement()
		{
		    --_iter_tuple;
		}
    void	advance(difference_type n)
		{
		    _iter_tuple += n;
		}
    template <class ITER_TUPLE_>
    std::enable_if_t<std::is_convertible<ITER_TUPLE_, ITER_TUPLE>::value,
		     difference_type>
		distance_to(const zip_iterator<ITER_TUPLE_>& iter) const
		{
		    return std::get<0>(iter.get_iterator_tuple())
			 - std::get<0>(_iter_tuple);
		}

  private:
    ITER_TUPLE	_iter_tuple;
};

template <class ITER_TUPLE>
inline std::enable_if_t<is_tuple<ITER_TUPLE>::value, zip_iterator<ITER_TUPLE> >
make_zip_iterator(ITER_TUPLE iter_tuple)
{
    return {iter_tuple};
}

template <class... ITERS>
inline std::enable_if_t<(sizeof...(ITERS) > 1),
			zip_iterator<std::tuple<ITERS...> > >
make_zip_iterator(ITERS... iters)
{
    return {std::make_tuple(iters...)};
}

/************************************************************************
*  type alias: decayed_iterator_value<ITER>				*
************************************************************************/
namespace detail
{
  template <class ITER>
  struct decayed_iterator_value
  {
      using type = iterator_value<ITER>;
  };
  template <class... ITER>
  struct decayed_iterator_value<zip_iterator<std::tuple<ITER...> > >
  {
      using type = std::tuple<typename decayed_iterator_value<ITER>::type...>;
  };
}	// namespace detail

//! 反復子が指す型を返す．
/*!
  zip_iterator<ITER_TUPLE>::value_type はITER_TUPLE中の各反復子が指す値への
  参照のtupleの型であるが，decayed_iterator_value<zip_iterator<ITER_TUPLE> >
  は，ITER_TUPLE中の各反復子が指す値そのもののtupleの型を返す．
  \param ITER	反復子
*/
template <class ITER>
using decayed_iterator_value = typename detail::decayed_iterator_value<ITER>
					      ::type;

/************************************************************************
*  TU::[begin|end|rbegin|rend](T&&), TU::size(const T&)			*
************************************************************************/
template <class T> inline auto
begin(T&& x) -> decltype(x.begin())
{
    return x.begin();
}

template <class T> inline auto
end(T&& x) -> decltype(x.end())
{
    return x.end();
}
    
template <class T> inline auto
rbegin(T&& x) -> decltype(x.rbegin())
{
    return x.rbegin();
}
    
template <class T> inline auto
rend(T&& x) -> decltype(x.rend())
{
    return x.rend();
}

template <class T> inline auto
size(const T& x) -> decltype(x.size())
{
    return x.size();
}

template <class T, size_t N> constexpr T*
begin(T (&array)[N]) noexcept
{
    return array;
}

template <class T, size_t N> constexpr T*
end(T (&array)[N]) noexcept
{
    return array + N;
}

template <class T, size_t N> std::reverse_iterator<T*>
rbegin(T (&array)[N])
{
    return {array + N};
}

template <class T, size_t N> std::reverse_iterator<T*>
rend(T (&array)[N])
{
    return {array};
}

template <class T, size_t N> inline constexpr size_t
size(const T (&array)[N]) noexcept
{
    return N;
}

template <class TUPLE, std::enable_if_t<is_tuple<TUPLE>::value>* = nullptr>
inline auto
begin(TUPLE&& t)
{
    return TU::make_zip_iterator(tuple_transform(
				     [](auto&& x)
				     { using std::begin; return begin(x); },
				     std::forward<TUPLE>(t)));
}

template <class TUPLE, std::enable_if_t<is_tuple<TUPLE>::value>* = nullptr>
inline auto
end(TUPLE&& t)
{
    return TU::make_zip_iterator(tuple_transform(
				     [](auto&& x)
				     { using std::end; return end(x); },
				     std::forward<TUPLE>(t)));
}

template <class TUPLE, std::enable_if_t<is_tuple<TUPLE>::value>* = nullptr>
inline auto
rbegin(TUPLE&& t)
{
    return std::make_reverse_iterator(end(std::forward<TUPLE>(t)));
}

template <class TUPLE, std::enable_if_t<is_tuple<TUPLE>::value>* = nullptr>
inline auto
rend(TUPLE&& t)
{
    return std::make_reverse_iterator(begin(std::forward<TUPLE>(t)));
}

template <class... T> inline auto
cbegin(const std::tuple<T...>& t)
{
    return begin(t);
}

template <class... T> inline auto
cend(const std::tuple<T...>& t)
{
    return end(t);
}

template <class... T> inline auto
crbegin(const std::tuple<T...>& t)
{
    return rbegin(t);
}

template <class... T> inline auto
crend(const std::tuple<T...>& t)
{
    return rend(t);
}

template <class... T> inline auto
size(const std::tuple<T...>& t)
{
    return TU::size(std::get<0>(t));
}

/************************************************************************
*  map_iterator<FUNC, ITER...>						*
************************************************************************/
template <class FUNC, class... ITER>
class map_iterator
    : public boost::iterator_adaptor<
	map_iterator<FUNC, ITER...>,
	zip_iterator<std::tuple<ITER...> >,
	std::decay_t<std::result_of_t<FUNC(iterator_reference<ITER>...)> >,
	boost::use_default,
	std::result_of_t<FUNC(iterator_reference<ITER>...)> >
{
  private:
    using ref	= std::result_of_t<FUNC(iterator_reference<ITER>...)>;
    using super	= boost::iterator_adaptor<map_iterator,
					  zip_iterator<std::tuple<ITER...> >,
					  std::decay_t<ref>,
					  boost::use_default,
					  ref>;
    friend	class boost::iterator_core_access;

  public:
    using	typename super::reference;
	
  public:
		map_iterator(FUNC&& func,
			     const std::tuple<ITER...>& iter_tuple)
		    :super(iter_tuple), _func(std::forward<FUNC>(func))
		{
		}
    const auto&	functor() const
		{
		    return _func;
		}
	
  private:
    reference	dereference() const
		{
		    return apply(_func, *super::base());
		}
	
  private:
    FUNC	_func;	//!< 演算子
};

template <class FUNC, class... ITER> inline map_iterator<FUNC, ITER...>
make_map_iterator(FUNC&& func, const ITER&... iter)
{
    return {std::forward<FUNC>(func), std::make_tuple(iter...)};
}
    
template <class FUNC, class... ITER> inline map_iterator<FUNC, ITER...>
make_map_iterator(FUNC&& func, const std::tuple<ITER...>& iter_tuple)
{
    return {std::forward<FUNC>(func), iter_tuple};
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
    return make_map_iterator(
	       std::function<std::conditional_t<
				 std::is_same<iterator_pointer<ITER>,
					      iterator_value<ITER>*>::value,
				 T&, const T&>(iterator_reference<ITER>)>(
				     std::mem_fn(mbr)),
	       iter);
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
      std::enable_if_t<is_binary_func<T_>::value, assignment_proxy&>
			operator =(T_&& val)
			{
			    _func(*_iter, std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      std::enable_if_t<!is_binary_func<T_>::value, assignment_proxy&>
			operator =(T_&& val)
			{
			    *_iter  = _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator +=(T_&& val)
			{
			    *_iter += _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator -=(T_&& val)
			{
			    *_iter -= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator *=(T_&& val)
			{
			    *_iter *= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator /=(T_&& val)
			{
			    *_iter /= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator %=(T_&& val)
			{
			    *_iter %= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator &=(T_&& val)
			{
			    *_iter &= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator |=(T_&& val)
			{
			    *_iter |= _func(std::forward<T_>(val));
			    return *this;
			}
      template <class T_>
      assignment_proxy&	operator ^=(T_&& val)
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
    assignment_iterator(const FUNC& func, const ITER& iter)
	:super(iter), _func(func)	{}

    const auto&	functor()	const	{ return _func; }

  private:
    reference	dereference()	const	{ return {super::base(), _func}; }
    
  private:
    FUNC 	_func;	// 代入を可能にするためconstは付けない
};
    
template <class FUNC, class ITER> inline assignment_iterator<FUNC, ITER>
make_assignment_iterator(const FUNC& func, const ITER& iter)
{
    return {func, iter};
}

template <class FUNC, class... ITERS>
inline std::enable_if_t<(sizeof...(ITERS) > 1),
			assignment_iterator<
			    FUNC, zip_iterator<std::tuple<ITERS...> > > >
make_assignment_iterator(const FUNC& func, const ITERS&... iters)
{
    return {func, make_zip_iterator(iters...)};
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
			    using	std::begin;
			    
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