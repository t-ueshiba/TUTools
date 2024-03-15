/*!
  \file		tuple.h
  \author	Toshio UESHIBA
  \brief	std::tupleの用途拡張のためのユティリティ
*/
#ifndef TU_TUPLE_H
#define TU_TUPLE_H

#include <tuple>
#include <utility>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include <iostream>
#include "TU/type_traits.h"	// for TU::any<PRED, T...>

namespace TU
{
/************************************************************************
*  predicate: is_tuple<T>						*
************************************************************************/
//! 与えられた型が std::tuple 又はそれに変換可能であるか判定する
/*!
  \param T	判定対象となる型
*/ 
template <class T>
using is_tuple = is_convertible<T, std::tuple>;

/************************************************************************
*  type alias: tuple_head<T>						*
************************************************************************/
namespace detail
{
  template <class HEAD, class... TAIL>
  HEAD	tuple_head(const std::tuple<HEAD, TAIL...>&)			;
  template <class T>
  T	tuple_head(const T&)						;
}	// namespace detail
    
//! 与えられた型がtupleならばその先頭要素の型を，そうでなければ元の型を返す．
/*!
  \param T	その先頭要素の型を調べるべき型
*/
template <class T>
using tuple_head = decltype(detail::tuple_head(std::declval<T>()));

/************************************************************************
*  type alias: replace_element<S, T>					*
************************************************************************/
namespace detail
{
  template <class S, class T>
  struct replace_element : std::conditional<std::is_void<T>::value, S, T>
  {
  };
  template <class... S, class T>
  struct replace_element<std::tuple<S...>, T>
  {
      using type = std::tuple<typename replace_element<S, T>::type...>;
  };
}	// namespace detail
    
//! 与えられた型がstd::tupleならばその要素の型を，そうでなければ元の型自身を別の型で置き換える．
/*!
  \param S	要素型置換の対象となる型
  \param T	置換後の要素の型．voidならば置換しない．
*/
template <class S, class T>
using replace_element = typename detail::replace_element<S, T>::type;

/************************************************************************
*  make_reference_wrapper(T&&)						*
************************************************************************/
//! 与えられた値の型に応じて実引数を生成する
/*!
  \param x	関数に渡す引数
  \return	xが右辺値参照ならばx自身，定数参照ならばstd::cref(x),
		非定数参照ならばstd::ref(x)
*/
template <class T>
inline std::conditional_t<std::is_lvalue_reference<T>::value,
			  std::reference_wrapper<std::remove_reference_t<T> >,
			  T&&>
make_reference_wrapper(T&& x)
{
    return std::forward<T>(x);
}
    
/************************************************************************
*  tuple_for_each(FUNC, TUPLES&&...)				`	*
************************************************************************/
namespace detail
{
  template <size_t I, class T, std::enable_if_t<!is_tuple<T>::value>* = nullptr>
  inline decltype(auto)
  tuple_get(T&& x)
  {
      return std::forward<T>(x);
  }
  template <size_t I, class T, std::enable_if_t<is_tuple<T>::value>* = nullptr>
  inline decltype(auto)
  tuple_get(T&& x)
  {
      return std::get<I>(std::forward<T>(x));
  }
    
  template <class... T>
  struct first_tuple_size;
  template<>
  struct first_tuple_size<>
  {
      constexpr static size_t	value = 0;
  };
  template <class HEAD, class... TAIL>
  struct first_tuple_size<HEAD, TAIL...>
  {
      template <class T_>
      struct tuple_size
      {
	  constexpr static size_t	value = 0;
      };
      template <class... T_>
      struct tuple_size<std::tuple<T_...> >
      {
	  constexpr static size_t	value = sizeof...(T_);
      };

      using TUPLE = std::decay_t<HEAD>;
      
      constexpr static size_t	value = (tuple_size<TUPLE>::value ?
					 tuple_size<TUPLE>::value :
					 first_tuple_size<TAIL...>::value);
  };
    
  template <class FUNC, class... TUPLES> inline void
  tuple_for_each(std::index_sequence<>, FUNC&&, TUPLES&&...)
  {
  }
  template <size_t I, size_t... IDX, class FUNC, class... TUPLES> inline void
  tuple_for_each(std::index_sequence<I, IDX...>, FUNC&& f, TUPLES&&... x)
  {
      f(tuple_get<I>(std::forward<TUPLES>(x))...);
      tuple_for_each(std::index_sequence<IDX...>(),
		     std::forward<FUNC>(f), std::forward<TUPLES>(x)...);
  }
}	// namespace detail
    
template <class FUNC, class... TUPLES>
inline std::enable_if_t<any<is_tuple, TUPLES...>::value>
tuple_for_each(FUNC&& f, TUPLES&&... x)
{
    detail::tuple_for_each(std::make_index_sequence<
			       detail::first_tuple_size<TUPLES...>::value>(),
			   std::forward<FUNC>(f), std::forward<TUPLES>(x)...);
}

/************************************************************************
*  tuple_transform(FUNC, TUPLES&&...)					*
************************************************************************/
namespace detail
{
  template <class FUNC, class... TUPLES> inline auto
  tuple_transform(std::index_sequence<>, FUNC&&, TUPLES&&...)
  {
      return std::tuple<>();
  }
  template <class FUNC, class... TUPLES, size_t I, size_t... IDX> inline auto
  tuple_transform(std::index_sequence<I, IDX...>, FUNC&& f, TUPLES&&... x)
  {
      auto&&	val = f(tuple_get<I>(std::forward<TUPLES>(x))...);
      return std::tuple_cat(
		std::make_tuple(
		    make_reference_wrapper(std::forward<decltype(val)>(val))),
		tuple_transform(std::index_sequence<IDX...>(),
				std::forward<FUNC>(f),
				std::forward<TUPLES>(x)...));
  }
}	// namespace detail
    
template <class FUNC, class... TUPLES,
	  std::enable_if_t<any<is_tuple, TUPLES...>::value>* = nullptr>
inline auto
tuple_transform(FUNC&& f, TUPLES&&... x)
{
    return detail::tuple_transform(
	       std::make_index_sequence<
		   detail::first_tuple_size<TUPLES...>::value>(),
	       std::forward<FUNC>(f), std::forward<TUPLES>(x)...);
}

/************************************************************************
*  Arithmetic operators							*
************************************************************************/
template <class E, std::enable_if_t<is_tuple<E>::value>* = nullptr> inline auto
operator -(E&& expr)
{
    return tuple_transform([](auto&& x)
			   { return -std::forward<decltype(x)>(x); },
			   std::forward<E>(expr));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator +(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  + std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator -(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  - std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator *(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  * std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator /(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  / std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator %(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  % std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator +=(L&& l, const R& r)
{
    tuple_for_each([](auto&& x, const auto& y){ x += y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator -=(L&& l, const R& r)
{
    tuple_for_each([](auto&& x, const auto& y){ x -= y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator *=(L&& l, const R& r)
{
    tuple_for_each([](auto&& x, const auto& y){ x *= y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator /=(L&& l, const R& r)
{
    tuple_for_each([](auto&& x, const auto& y){ x /= y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator %=(L&& l, const R& r)
{
    tuple_for_each([](auto&& x, const auto& y){ x %= y; }, l, r);
    return l;
}

template <class T> inline std::enable_if_t<is_tuple<T>::value, T&>
operator ++(T&& t)
{
    tuple_for_each([](auto&& x){ ++x; }, t);
    return t;
}

template <class T> inline std::enable_if_t<is_tuple<T>::value, T&>
operator --(T&& t)
{
    tuple_for_each([](auto&& x){ --x; }, t);
    return t;
}

template <class L, class C, class R,
	  std::enable_if_t<any<is_tuple, L, C, R>::value>* = nullptr>
inline auto
fma(L&& l, C&& c, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y, auto&& z)
			   { return fma(std::forward<decltype(x)>(x),
					std::forward<decltype(y)>(y),
					std::forward<decltype(z)>(z)); },
			   std::forward<L>(l),
			   std::forward<C>(c),
			   std::forward<R>(r));
}

/************************************************************************
*  Bit operators							*
************************************************************************/
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator &(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  & std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator |(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  | std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}

template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator ^(L&& l, R&& r)
{
    return tuple_transform([](auto&& x, auto&& y)
			   { return std::forward<decltype(x)>(x)
				  ^ std::forward<decltype(y)>(y); },
			   std::forward<L>(l), std::forward<R>(r));
}
    
template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator &=(L&& l, const R& r)
{
    tuple_for_each([](auto& x, const auto& y){ x &= y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator |=(L&& l, const R& r)
{
    tuple_for_each([](auto& x, const auto& y){ x |= y; }, l, r);
    return l;
}

template <class L, class R> inline std::enable_if_t<is_tuple<L>::value, L&>
operator ^=(L&& l, const R& r)
{
    tuple_for_each([](auto& x, const auto& y){ x ^= y; }, l, r);
    return l;
}

/************************************************************************
*  Logical operators							*
************************************************************************/
template <class... T> inline auto
operator !(const std::tuple<T...>& t)
{
    return tuple_transform([](const auto& x){ return !x; }, t);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator &&(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x && y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator ||(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x || y; }, l, r);
}
    
/************************************************************************
*  Relational operators							*
************************************************************************/
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator ==(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x == y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator !=(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x != y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator <(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x < y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator >(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x > y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator <=(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x <= y; }, l, r);
}
    
template <class L, class R,
	  std::enable_if_t<any<is_tuple, L, R>::value>* = nullptr> inline auto
operator >=(const L& l, const R& r)
{
    return tuple_transform([](const auto& x, const auto& y)
			   { return x >= y; }, l, r);
}

/************************************************************************
*  Selection								*
************************************************************************/
template <class X, class Y> inline auto
select(bool s, X&& x, Y&& y)
{
    return (s ? std::forward<X>(x) : std::forward<Y>(y));
}
    
template <class... S, class X, class Y,
	  std::enable_if_t<any<is_tuple, X, Y>::value>* = nullptr> inline auto
select(const std::tuple<S...>& s, X&& x, Y&& y)
{
    return tuple_transform([](const auto& t, auto&& u, auto&& v)
			   { return select(t,
					   std::forward<decltype(u)>(u),
					   std::forward<decltype(v)>(v)); },
			   s, std::forward<X>(x), std::forward<Y>(y));
}

/************************************************************************
*  class zip_iterator<ITER_TUPLE>					*
************************************************************************/
namespace detail
{
  struct generic_dereference
  {
    // assignment_iterator<FUNC, ITER> のように dereference すると
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
			typename std::iterator_traits<
			    std::tuple_element_t<0, ITER_TUPLE> >
			       ::iterator_category,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>()))>
{
  private:
    using super = boost::iterator_facade<
			zip_iterator,
			decltype(tuple_transform(detail::generic_dereference(),
						 std::declval<ITER_TUPLE>())),
			typename std::iterator_traits<
			    std::tuple_element_t<0, ITER_TUPLE> >
			       ::iterator_category,
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

    const auto&	get_iterator_tuple() const
		{
		    return _iter_tuple;
		}
    
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

template <class... ITERS> inline zip_iterator<std::tuple<ITERS...> >
make_zip_iterator(const std::tuple<ITERS...>& iter_tuple)
{
    return {iter_tuple};
}

template <class ITER> inline ITER
make_zip_iterator(const ITER& iter)
{
    return iter;
}

template <class ITER, class... ITERS> inline auto
make_zip_iterator(const ITER& iter, const ITERS&... iters)
{
    return make_zip_iterator(std::make_tuple(iter, iters...));
}

/************************************************************************
*  TU::[begin|end|rbegin|rend|size](TUPLE&&)				*
************************************************************************/
template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
begin(std::tuple<T...>& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::begin; return begin(x); }, t));
}

template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
end(std::tuple<T...>& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::end; return end(x); }, t));
}

template <class... T> inline auto
rbegin(std::tuple<T...>& t) -> decltype(std::make_reverse_iterator(end(t)))
{
    return std::make_reverse_iterator(end(t));
}

template <class... T> inline auto
rend(std::tuple<T...>& t) -> decltype(std::make_reverse_iterator(begin(t)))
{
    return std::make_reverse_iterator(begin(t));
}

template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
begin(std::tuple<T...>&& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::begin; return begin(x); }, t));
}

template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
end(std::tuple<T...>&& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::end; return end(x); }, t));
}

template <class... T> inline auto
rbegin(std::tuple<T...>&& t) -> decltype(std::make_reverse_iterator(end(t)))
{
    return std::make_reverse_iterator(end(t));
}

template <class... T> inline auto
rend(std::tuple<T...>&& t) -> decltype(std::make_reverse_iterator(begin(t)))
{
    return std::make_reverse_iterator(begin(t));
}

template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
begin(const std::tuple<T...>& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::begin; return begin(x); }, t));
}

template <class... T,
	  std::enable_if_t<all<is_iterable, T...>::value>* = nullptr>
inline auto
end(const std::tuple<T...>& t)
{
    return make_zip_iterator(
		tuple_transform([](auto& x)
				{ using std::end; return end(x); }, t));
}

template <class... T> inline auto
rbegin(const std::tuple<T...>& t)
    -> decltype(std::make_reverse_iterator(end(t)))
{
    return std::make_reverse_iterator(end(t));
}

template <class... T> inline auto
rend(const std::tuple<T...>& t)
    -> decltype(std::make_reverse_iterator(begin(t)))
{
    return std::make_reverse_iterator(begin(t));
}

template <class... T> inline auto
cbegin(const std::tuple<T...>& t) -> decltype(begin(t))
{
    return begin(t);
}

template <class... T> inline auto
cend(const std::tuple<T...>& t) -> decltype(end(t))
{
    return end(t);
}

template <class... T> inline auto
crbegin(const std::tuple<T...>& t) -> decltype(rbegin(t))
{
    return rbegin(t);
}

template <class... T> inline auto
crend(const std::tuple<T...>& t) -> decltype(rend(t))
{
    return rend(t);
}

template <class... T> inline std::size_t
size(const std::tuple<T...>& t)
{
    using	std::size;
    
    return size(std::get<0>(t));
}

/************************************************************************
*  type alias: decayed_iterator_value<ITER>				*
************************************************************************/
namespace detail
{
  template <class ITER>
  struct decayed_iterator_value
  {
      using type = typename std::iterator_traits<ITER>::value_type;
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
*  Applying a multi-input function to a tuple of arguments		*
************************************************************************/
template <class FUNC, class TUPLE,
	  std::enable_if_t<is_tuple<TUPLE>::value>* = nullptr>
inline decltype(auto)
apply(FUNC&& f, TUPLE&& t)
{
    return std::apply(std::forward<FUNC>(f), std::forward<TUPLE>(t));
}

template <class FUNC, class T,
	  std::enable_if_t<!is_tuple<T>::value>* = nullptr>
inline decltype(auto)
apply(FUNC&& f, T&& t)
{
    return f(std::forward<T>(t));
}

}	// namespace TU

namespace std
{
/************************************************************************
*  I/O functions							*
************************************************************************/
template <class... T> inline ostream&
operator <<(ostream& out, const tuple<T...>& t)
{
    out << '(';
    TU::tuple_for_each([&out](const auto& x){ out << ' ' << x; }, t);
    out << ')';

    return out;
}

}	// namespace std
#endif	// !TU_TUPLE_H
