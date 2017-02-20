/*
 *  $Id$
 */
/*!
  \file		range.h
  \brief	2つの反復子によって指定されるレンジの定義と実装
*/
#ifndef __TU_RANGE_H
#define __TU_RANGE_H

#include <iostream>
#include <cassert>
#include <initializer_list>
#include "TU/iterator.h"
#include "TU/functional.h"
#include "TU/algorithm.h"

namespace TU
{
/************************************************************************
*  predicates is_range<E>, is_range1<E>, is_range2<E>			*
************************************************************************/
namespace detail
{
  template <class E> static auto
  is_range(const E& x) -> decltype(std::begin(x), std::true_type())	;
  static std::false_type
  is_range(...)								;

  template <class E> static auto
  is_range2(const E& x) -> decltype(std::begin(*std::begin(x)),
				    std::true_type())			;
  static std::false_type
  is_range2(...)							;
}	// namespace detail

template <class E>
using is_range	= decltype(detail::is_range(std::declval<E>()));

template <class E>
using is_range2	= decltype(detail::is_range2(std::declval<E>()));

template <class E>
using is_range1	= std::integral_constant<bool, (is_range<E>::value &&
						!is_range2<E>::value)>;
    
/************************************************************************
*  type aliases								*
************************************************************************/
namespace detail
{
  template <class E_> static auto
  has_const_iterator(const E_& x) -> decltype(std::begin(x))		;
  static void
  has_const_iterator(...)						;
}	// namespace detail
    
//! 式が持つ定数反復子の型を返す
/*!
  \param E	式の型
  \return	E が定数反復子を持てばその型，持たなければ void
*/
template <class E>
using const_iterator_t	= decltype(detail::has_const_iterator(
				       std::declval<E>()));

namespace detail
{
  template <class T>
  struct identity
  {
      using type = T;
  };

  template <class E>
  struct value_t
  {
      using type = typename std::iterator_traits<const_iterator_t<E> >
			       ::value_type;
  };
      
  template <class E, class=const_iterator_t<E> >
  struct element_t
  {
      using F	 = typename value_t<E>::type;
      using type = typename element_t<F, const_iterator_t<F> >::type;
  };
  template <class E>
  struct element_t<E, void> : identity<E>				{};
}	// namespace detail
    
//! 式が持つ定数反復子が指す型を返す
/*!
  定数反復子を持たない式を与えるとコンパイルエラーとなる.
  \param E	定数反復子を持つ式の型
  \return	E の定数反復子が指す型
*/
template <class E>
using value_t	= typename detail::value_t<E>::type;

//! 式が持つ定数反復子が指す型を再帰的に辿って到達する型を返す
/*!
  \param E	式の型
  \return	E が定数反復子を持てばそれが指す型を再帰的に辿って到達する型，
		持たなければ E 自身
*/
template <class E>
using element_t	= typename detail::element_t<E>::type;

/************************************************************************
*  sizes and strides of multidimensional ranges				*
************************************************************************/
namespace detail
{
  template <class E> static auto
  has_size0(const E& x) -> decltype(E::size0(), std::true_type())	;
  static std::false_type
  has_size0(...)							;

  /*
   *  A ^ b において opnode: product_opnode<bit_xor, L, R> が
   *  生成されるが，これを評価して2次元配列に代入する際に代入先の領域確保のため
   *  size<0>(opnode), size<1>(opnode) が呼び出される．後者はopnodeの反復子が
   *  指す先を評価するが，これは2つの3次元ベクトルのベクトル積を評価することに
   *  等しく，コストが高い処理である．そこで，ベクトル積を評価せずその評価結果の
   *  サイズだけを得るために，以下のオーバーロードを導入する．
   */
  template <class OP, class L, class R>	class product_opnode;
    
  template <class L, class R> constexpr size_t
  size(const product_opnode<bit_xor, L, R>&, std::integral_constant<size_t, 1>)
  {
      return 3;
  }

  template <class E> inline size_t
  size(const E& expr, std::integral_constant<size_t, 0>)
  {
      return std::size(expr);
  }
  template <size_t I, class E> inline size_t
  size(const E& expr, std::integral_constant<size_t, I>)
  {
      return size(*std::begin(expr), std::integral_constant<size_t, I-1>());
  }

  template <class E> inline size_t
  stride(const E& expr, std::integral_constant<size_t, 1>)
  {
      return std::begin(expr).stride();
  }
  template <size_t I, class E> inline size_t
  stride(const E& expr, std::integral_constant<size_t, I>)
  {
      return stride(*std::begin(expr), std::integral_constant<size_t, I-1>());
  }
}	// namespace detail

template <class E>
constexpr inline typename std::enable_if<decltype(
    detail::has_size0(std::declval<E>()))::value, size_t>::type
size0()
{
    return E::size0();
}
template <class E>
constexpr inline typename std::enable_if<!decltype(
    detail::has_size0(std::declval<E>()))::value, size_t>::type
size0()
{
    return 0;
}
    
template <size_t I, class E>
inline typename std::enable_if<is_range<E>::value, size_t>::type
size(const E& expr)
{
    return detail::size(expr, std::integral_constant<size_t, I>());
}

template <size_t I, class E>
inline typename std::enable_if<is_range<E>::value, size_t>::type
stride(const E& expr)
{
    return detail::stride(expr, std::integral_constant<size_t, I>());
}

/************************************************************************
*  manipulator sizes_and_strides					*
************************************************************************/
template <class E>
class sizes_holder
{
  public:
			sizes_holder(const E& expr)	:_expr(expr)	{}

    std::ostream&	operator ()(std::ostream& out) const
			{
			    return print_size(out, _expr);
			}
    
  protected:
    template <class E_>
    static typename std::enable_if<!is_range<E_>::value, std::ostream&>::type
			print_x(std::ostream& out, const E_& expr)
			{
			    return out;
			}
    template <class E_>
    static typename std::enable_if<is_range<E_>::value, std::ostream&>::type
			print_x(std::ostream& out, const E_& expr)
			{
			    return out << 'x';
			}

  private:
    template <class E_>
    static typename std::enable_if<!is_range<E_>::value, std::ostream&>::type
			print_size(std::ostream& out, const E_& expr)
			{
			    return out;
			}
    template <class E_>
    static typename std::enable_if<is_range<E_>::value, std::ostream&>::type
			print_size(std::ostream& out, const E_& expr)
			{
			    const auto&	val = *std::begin(expr);
			    return print_size(print_x(out << std::size(expr),
						      val), val);
			}

  protected:
    const E&	_expr;
};
    
template <class E>
class sizes_and_strides_holder : public sizes_holder<E>
{
  private:
    using super		= sizes_holder<E>;
    
  public:
			sizes_and_strides_holder(const E& expr)
			    :super(expr)				{}

    std::ostream&	operator ()(std::ostream& out) const
			{
			    return print_stride(super::operator ()(out) << ':',
						super::_expr.begin());
			}
    
  private:
    template <class ITER_>
    static typename std::enable_if<!is_range<
				       typename std::iterator_traits<ITER_>
				       ::value_type>::value,
				   std::ostream&>::type
			print_stride(std::ostream& out, const ITER_& iter)
			{
			    return out;
			}
    template <class ITER_>
    static typename std::enable_if<is_range<
				       typename std::iterator_traits<ITER_>
				       ::value_type>::value,
				   std::ostream&>::type
			print_stride(std::ostream& out, const ITER_& iter)
			{
			    return print_stride(super::print_x(
						    out << iter.stride(),
						    *iter->begin()),
						iter->begin());
			}
};

template <class E> sizes_holder<E>
print_sizes(const E& expr)
{
    return sizes_holder<E>(expr);
}

template <class E> std::ostream&
operator <<(std::ostream& out, const sizes_holder<E>& holder)
{
    return holder(out);
}

template <class E> sizes_and_strides_holder<E>
print_sizes_and_strides(const E& expr)
{
    return sizes_and_strides_holder<E>(expr);
}

template <class E> std::ostream&
operator <<(std::ostream& out, const sizes_and_strides_holder<E>& holder)
{
    return holder(out);
}

/************************************************************************
*  class range<ITER, SIZE>						*
************************************************************************/
//! 2つの反復子によって指定される範囲(レンジ)を表すクラス
/*!
  \param ITER	反復子の型
  \param SIZE	レンジに含まれる要素数(0ならば可変長)
*/
template <class ITER, size_t SIZE=0>
class range
{
  public:
    using iterator		= ITER;
    using reverse_iterator	= std::reverse_iterator<iterator>;
    using value_type		= typename std::iterator_traits<iterator>
					      ::value_type;
    using reference		= typename std::iterator_traits<iterator>
					      ::reference;

  public:
		range(iterator begin)		:_begin(begin)	{}
		range(iterator begin, size_t)	:_begin(begin)	{}
    
		range()						= delete;
		range(const range&)				= default;
    range&	operator =(const range& r)
		{
		    copy<SIZE>(r.begin(), size(), begin());
		    return *this;
		}
		range(range&&)					= default;
    range&	operator =(range&&)				= default;
    
    template <class ITER_,
	      typename std::enable_if<
		  std::is_convertible<ITER_, iterator>::value>::type* = nullptr>
		range(const range<ITER_, SIZE>& r)
		    :_begin(r.begin())
		{
		}

    template <class E_>
    typename std::enable_if<is_range<E_>::value, range&>::type
		operator =(const E_& expr)
		{
		    assert(std::size(expr) == size());
		    copy<SIZE>(std::begin(expr), size(), begin());
		    return *this;
		}

		range(std::initializer_list<value_type> args)
		    :_begin(const_cast<iterator>(args.begin()))
    		{
		    assert(args.size() == size());
		}
    range&	operator =(std::initializer_list<value_type> args)
		{
		    assert(args.size() == size());
		    copy<SIZE>(args.begin(), size(), begin());
		    return *this;
		}

    constexpr static
    size_t	size0()			{ return SIZE; }
    constexpr static
    size_t	size()			{ return SIZE; }
    auto	begin()		const	{ return _begin; }
    auto	end()		const	{ return _begin + SIZE; }
    auto	cbegin()	const	{ return begin(); }
    auto	cend()		const	{ return end(); }
    auto	rbegin()	const	{ return reverse_iterator(end()); }
    auto	rend()		const	{ return reverse_iterator(begin()); }
    auto	crbegin()	const	{ return rbegin(); }
    auto	crend()		const	{ return rend(); }
    reference	operator [](size_t i) const
		{
		    assert(i < size());
		    return *(_begin + i);
		}
    
  private:
    const iterator	_begin;
};

//! 可変長レンジ
/*!
  \param ITER	反復子の型
*/
template <class ITER>
class range<ITER, 0>
{
  public:
    using iterator		= ITER;
    using reverse_iterator	= std::reverse_iterator<iterator>;
    using value_type		= typename std::iterator_traits<iterator>
					      ::value_type;
    using reference		= typename std::iterator_traits<iterator>
					      ::reference;

  public:
		range(iterator begin, size_t size)
		    :_begin(begin), _size(size)			{}
    
		range()						= delete;
		range(const range&)				= default;
    range&	operator =(const range& r)
		{
		    assert(r.size() == size());
		    copy<0>(r._begin, size(), begin());
		    return *this;
		}
		range(range&&)					= default;
    range&	operator =(range&&)				= default;
    
    template <class ITER_,
	      typename std::enable_if<
		  std::is_convertible<ITER_, iterator>::value>::type* = nullptr>
		range(const range<ITER_, 0>& r)
		    :_begin(r.begin()), _size(r.size())
		{
		}

    template <class E_>
    typename std::enable_if<is_range<E_>::value, range&>::type
		operator =(const E_& expr)
		{
		    assert(std::size(expr) == size());
		    copy<TU::size0<E_>()>(std::begin(expr), size(), begin());
		    return *this;
		}
		
		range(std::initializer_list<value_type> args)
		    :_begin(const_cast<iterator>(args.begin())),
		     _size(args.size())
    		{
		}
    range&	operator =(std::initializer_list<value_type> args)
		{
		    assert(args.size() == size());
		    copy<0>(args.begin(), size(), begin());
		    return *this;
		}
		
    constexpr static
    size_t	size0()			{ return 0; }
    size_t	size()		const	{ return _size; }
    auto	begin()		const	{ return _begin; }
    auto	end()		const	{ return _begin + _size; }
    auto	cbegin()	const	{ return begin(); }
    auto	cend()		const	{ return end(); }
    auto	rbegin()	const	{ return reverse_iterator(end()); }
    auto	rend()		const	{ return reverse_iterator(begin()); }
    auto	crbegin()	const	{ return rbegin(); }
    auto	crend()		const	{ return rend(); }
    reference	operator [](size_t i) const
		{
		    assert(i < size());
		    return *(_begin + i);
		}

  private:
    const iterator	_begin;
    const size_t	_size;
};

//! 固定長レンジを生成する
/*!
  \param SIZE	レンジ長
  \param iter	レンジの先頭要素を指す反復子
*/
template <size_t SIZE, class ITER> inline range<ITER, SIZE>
make_range(ITER iter)
{
    return {iter};
}
    
//! 可変長レンジを生成する
/*!
  \param iter	レンジの先頭要素を指す反復子
  \param size	レンジの要素数
*/
template <size_t SIZE=0, class ITER> inline range<ITER, SIZE>
make_range(ITER iter, size_t size)
{
    return {iter, size};
}

//! 出力ストリームにレンジの内容を書き出す
/*!
  \param out	出力ストリーム
  \param r	レンジ
  \return	outで指定した出力ストリーム
*/
template <class ITER, size_t SIZE> std::ostream&
operator <<(std::ostream& out, const range<ITER, SIZE>& r)
{
    for (const auto& val : r)
	out << ' ' << val;
    return out << std::endl;
}
    
/************************************************************************
*  class range_iterator<ITER, SIZE, STRIDE>				*
************************************************************************/
//! 実装の詳細を収める名前空間
namespace detail
{
  template <size_t SIZE, size_t STRIDE>
  struct size_and_stride
  {
      constexpr static size_t	size()		{ return SIZE; }
      constexpr static size_t	stride()	{ return STRIDE; }
  };
  template <size_t SIZE>
  struct size_and_stride<SIZE, 0>
  {
      size_and_stride(size_t stride)
	  :_stride(stride)			{}
      constexpr static size_t	size()		{ return SIZE; }
      size_t			stride() const	{ return _stride; }

    private:
      size_t	_stride;
  };
  template <size_t STRIDE>
  struct size_and_stride<0, STRIDE>
  {
      size_and_stride(size_t size)
	  :_size(size)				{}
      size_t			size()	 const	{ return _size; }
      constexpr static size_t	stride()	{ return STRIDE; }

    private:
      size_t	_size;
  };
  template <>
  struct size_and_stride<0, 0>
  {
      size_and_stride(size_t size, size_t stride)
	  :_size(size), _stride(stride)		{}
      size_t			size()	 const	{ return _size; }
      size_t			stride() const	{ return _stride; }

    private:
      size_t	_size;
      size_t	_stride;
  };
}	// namespace detail
    
//! 配列を一定間隔に切り分けたレンジを指す反復子
/*!
  \param ITER	配列の要素を指す反復子の型
  \param SIZE	レンジ長(0ならば可変長)
  \param STRIDE	インクリメントしたときに進める要素数(0ならば可変)
*/
template <class ITER, size_t SIZE=0, size_t STRIDE=0>
class range_iterator
    : public boost::iterator_adaptor<range_iterator<ITER, SIZE, STRIDE>,
				     ITER,
				     range<ITER, SIZE>,
				     boost::use_default,
				     range<ITER, SIZE> >,
      public detail::size_and_stride<SIZE, STRIDE>
{
  private:
    using super	= boost::iterator_adaptor<range_iterator,
					  ITER,
					  range<ITER, SIZE>,
					  boost::use_default,
					  range<ITER, SIZE> >;
    using ss	= detail::size_and_stride<SIZE, STRIDE>;
    
  public:
    using		typename super::reference;
    using		typename super::difference_type;

    friend class	boost::iterator_core_access;
	  
  public:
    template <size_t SIZE_=SIZE, size_t STRIDE_=STRIDE,
	      typename std::enable_if<
		  (SIZE_ != 0) && (STRIDE_ != 0)>::type* = nullptr>
		range_iterator(ITER iter)
		    :super(iter), ss()					{}
    template <size_t SIZE_=SIZE, size_t STRIDE_=STRIDE,
	      typename std::enable_if<
		  (SIZE_ != 0) && (STRIDE_ == 0)>::type* = nullptr>
		range_iterator(ITER iter, size_t stride)
		    :super(iter), ss(stride)				{}
    template <size_t SIZE_=SIZE, size_t STRIDE_=STRIDE,
	      typename std::enable_if<
		  (SIZE_ == 0) && (STRIDE_ != 0)>::type* = nullptr>
		range_iterator(ITER iter, size_t size)
		    :super(iter), ss(size)				{}
    template <size_t SIZE_=SIZE, size_t STRIDE_=STRIDE,
	      typename std::enable_if<
		  (SIZE_ == 0) && (STRIDE_ == 0)>::type* = nullptr>
		range_iterator(ITER iter, size_t size, size_t stride)
		    :super(iter), ss(size, stride)			{}

    template <class ITER_,
	      typename std::enable_if<
		  std::is_convertible<ITER_, ITER>::value>::type* = nullptr>
		range_iterator(
		    const range_iterator<ITER_, SIZE, STRIDE>& iter)
		    :super(iter), ss(iter)				{}

    using	ss::size;
    using	ss::stride;
    
  private:
    reference	dereference() const
		{
		    return {super::base(), size()};
		}
    void	increment()
		{
		    std::advance(super::base_reference(), stride());
		}
    void	decrement()
		{
		    std::advance(super::base_reference(), -stride());
		}
    void	advance(difference_type n)
		{
		    std::advance(super::base_reference(), n*stride());
		}
    difference_type
		distance_to(const range_iterator& iter) const
		{
		    return std::distance(super::base(), iter.base()) / stride();
		}
};

//! 固定長レンジを指し，インクリメント時に固定した要素数だけ進める反復子を生成する
/*!
  \param SIZE	レンジ長
  \param STRIDE	インクリメント時に進める要素数
  \param iter	レンジの先頭要素を指す反復子
*/
template <size_t SIZE, size_t STRIDE, class ITER>
inline range_iterator<ITER, SIZE, STRIDE>
make_range_iterator(ITER iter)
{
    return {iter};
}

//! 固定長レンジを指し，インクリメント時に指定した要素数だけ進める反復子を生成する
/*!
  \param SIZE	レンジ長
  \param iter	レンジの先頭要素を指す反復子
  \param stride	インクリメント時に進める要素数
*/
template <size_t SIZE, class ITER> inline range_iterator<ITER, SIZE>
make_range_iterator(ITER iter, size_t stride)
{
    return {iter, stride};
}
    
//! 指定された長さのレンジを指し，インクリメント時に指定した要素数だけ進める反復子を生成する
/*!
  \param iter	レンジの先頭要素を指す反復子
  \param size	レンジ長
  \param stride	インクリメント時に進める要素数
*/
template <class ITER> inline range_iterator<ITER>
make_range_iterator(ITER iter, size_t size, size_t stride)
{
    return {iter, size, stride};
}
    
/************************************************************************
*  fixed size & fixed stride ranges and associated iterators		*
************************************************************************/
//! 多次元固定長レンジを指し，インクリメント時に固定したブロック数だけ進める反復子を生成する
/*!
  \param SIZE	最上位軸のレンジ長
  \param STRIDE	インクリメント時に進める最上位軸のブロック数
  \param SS	2番目以降の軸の(レンジ長，ストライド)の並び
  \param iter	レンジの先頭要素を指す反復子
*/
template <size_t SIZE, size_t STRIDE, size_t... SS, class ITER,
	  typename std::enable_if<sizeof...(SS) != 0>::type* = nullptr>
inline auto
make_range_iterator(ITER iter)
{
    return make_range_iterator<SIZE, STRIDE>(make_range_iterator<SS...>(iter));
}

template <size_t SIZE, size_t... SS, class ITER,
	  typename std::enable_if<sizeof...(SS) != 0>::type* = nullptr>
inline auto
make_range(ITER iter)
{
    return make_range<SIZE>(make_range_iterator<SS...>(iter));
}

/************************************************************************
*  fixed size & variable stride ranges and associated iterators		*
************************************************************************/
//! 多次元固定長レンジを指し，インクリメント時に指定したブロック数だけ進める反復子を生成する
/*!
  \param SIZE		最上位軸のレンジ長
  \param SIZES		2番目以降の軸のレンジ長の並び
  \param stride		最上位軸のストライド
  \param strides	2番目以降の軸のストライドの並び
  \param iter		レンジの先頭要素を指す反復子
*/
template <size_t SIZE, size_t... SIZES, class ITER, class... STRIDES,
	  typename std::enable_if<
	      sizeof...(SIZES) == sizeof...(STRIDES)>::type* = nullptr>
inline auto
make_range_iterator(ITER iter, size_t stride, STRIDES... strides)
{
    return make_range_iterator<SIZE>(
	       make_range_iterator<SIZES...>(iter, strides...), stride);
}

template <size_t SIZE, size_t... SIZES, class ITER, class... STRIDES,
	  typename std::enable_if<
	      sizeof...(SIZES) == sizeof...(STRIDES)>::type* = nullptr>
inline auto
make_range(ITER iter, STRIDES... strides)
{
    return make_range<SIZE>(make_range_iterator<SIZES...>(iter, strides...));
}

/************************************************************************
*  variable size & variable stride ranges and associated iterators	*
************************************************************************/
//! 多次元可変長レンジを指し，インクリメント時に指定したブロック数だけ進める反復子を生成する
/*!
  \param iter		レンジの先頭要素を指す反復子
  \param size		最上位軸のレンジ長
  \param stride		最上位軸のストライド
  \param ss		2番目以降の軸の(レンジ長, ストライド)の並び
*/
template <class ITER, class... SS> inline auto
make_range_iterator(ITER iter, size_t size, size_t stride, SS... ss)
{
    return make_range_iterator(make_range_iterator(iter, ss...),
			       size, stride);
}

template <class ITER, class... SS> inline auto
make_range(ITER iter, size_t size, SS... ss)
{
    return make_range(make_range_iterator(iter, ss...), size);
}

/************************************************************************
*  ranges with variable but identical size and stride			*
*  and associated iterators						*
************************************************************************/
template <class ITER> inline ITER
make_dense_range_iterator(ITER iter)
{
    return iter;
}
    
template <class ITER, class... SIZES> inline auto
make_dense_range_iterator(ITER iter, size_t size, SIZES... sizes)
{
    return make_range_iterator(make_dense_range_iterator(iter, sizes...),
			       size, size);
}

template <class ITER, class... SIZES> inline auto
make_dense_range(ITER iter, size_t size, SIZES... sizes)
{
    return make_range(make_dense_range_iterator(iter, sizes...), size);
}

/************************************************************************
*  subrange extraction							*
************************************************************************/
template <class ITER> inline ITER
make_subrange_iterator(ITER iter)
{
    return iter;
}

template <class ITER, class... IS> inline auto
make_subrange_iterator(const range_iterator<ITER>& iter,
		       size_t idx, size_t size, IS... is)
{
    return make_range_iterator(make_subrange_iterator(
				   iter->begin() + idx, is...),
			       size, iter.stride());
}

template <class RANGE, class... IS,
	  typename std::enable_if<
	      is_range<typename std::decay<RANGE>::type>::value>::type*
	  = nullptr>
inline auto
make_subrange(RANGE&& r, size_t idx, size_t size, IS... is)
{
    return make_range(make_subrange_iterator(r.begin() + idx, is...), size);
}

template <class RANGE, class... IS,
	  typename std::enable_if<is_range<RANGE>::value>::type* = nullptr>
inline auto
make_subrange(const RANGE& r, size_t idx, size_t size, IS... is)
{
    return make_range(make_subrange_iterator(r.begin() + idx, is...), size);
}

template <size_t SIZE, size_t... SIZES, class ITER, class... INDICES,
	  typename std::enable_if<
	      sizeof...(SIZES) == sizeof...(INDICES)>::type* = nullptr>
inline auto
make_subrange_iterator(const ITER& iter, size_t idx, INDICES... indices)
{
    return make_range_iterator<SIZE>(make_subrange_iterator<SIZES...>(
					 iter->begin() + idx, indices...),
				     iter.stride());
}

template <size_t SIZE, size_t... SIZES, class RANGE, class... INDICES,
	  typename std::enable_if<
	      is_range<typename std::decay<RANGE>::type>::value &&
	      (sizeof...(SIZES) == sizeof...(INDICES))>::type* = nullptr>
inline auto
make_subrange(RANGE&& r, size_t idx, INDICES... indices)
{
    return make_range<SIZE>(make_subrange_iterator<SIZES...>(
				r.begin() + idx, indices...));
}

template <size_t SIZE, size_t... SIZES, class RANGE, class... INDICES,
	  typename std::enable_if<
	      is_range<RANGE>::value &&
	      (sizeof...(SIZES) == sizeof...(INDICES))>::type* = nullptr>
inline auto
make_subrange(const RANGE& r, size_t idx, INDICES... indices)
{
    return make_range<SIZE>(make_subrange_iterator<SIZES...>(
				r.begin() + idx, indices...));
}

/************************************************************************
*  class column_iterator<ROW, NROWS>					*
************************************************************************/
//! 2次元配列の列を指す反復子
/*!
  \param ROW	begin(), end()をサポートするコンテナを指す反復子の型
  \param SIZE	begin()とend()間の距離(0ならば可変長)
*/ 
template <class ROW, size_t NROWS>
class column_iterator
    : public boost::iterator_adaptor<column_iterator<ROW, NROWS>,
				     size_t,
				     range<vertical_iterator<ROW>, NROWS>,
				     std::random_access_iterator_tag,
				     range<vertical_iterator<ROW>, NROWS>,
				     ptrdiff_t>
{
  private:
    using super	= boost::iterator_adaptor<column_iterator,
					  size_t,
					  range<vertical_iterator<ROW>, NROWS>,
					  std::random_access_iterator_tag,
					  range<vertical_iterator<ROW>, NROWS>,
					  ptrdiff_t>;

  public:
    using	typename super::reference;

    friend	class boost::iterator_core_access;

  public:
		column_iterator(ROW row, size_t nrows, size_t col)
		    :super(col), _row(row), _nrows(nrows)
		{
		}

  private:
    reference	dereference() const
		{
		    return {{_row, super::base()}, _nrows};
		}
    
  private:
    ROW		_row;
    size_t	_nrows;
};

template <size_t NROWS=0, class ROW> inline column_iterator<ROW, NROWS>
make_column_iterator(ROW row, size_t nrows, size_t col)
{
    return {row, nrows, col};
}

template <class E> inline auto
column_begin(E& expr)
{
    return make_column_iterator<size0<E>()>(std::begin(expr), std::size(expr),
					    0);
}
    
template <class E> inline auto
column_begin(const E& expr)
{
    return make_column_iterator<size0<E>()>(std::begin(expr), std::size(expr),
					    0);
}
    
template <class E> inline auto
column_cbegin(const E& expr)
{
    return column_begin(expr);
}
    
template <class E> inline auto
column_end(E& expr)
{
    return make_column_iterator<size0<E>()>(std::begin(expr), std::size(expr),
					    size<1>(expr));
}

template <class E> inline auto
column_end(const E& expr)
{
    return make_column_iterator<size0<E>()>(std::begin(expr), std::size(expr),
					    size<1>(expr));
}

template <class E> inline auto
column_cend(const E& expr)
{
    return column_end(expr);
}

/************************************************************************
*  various numeric functions						*
************************************************************************/
namespace detail
{
  /**********************************************************************
  *  type aliases							*
  **********************************************************************/
  template <size_t I, size_t J>
  using max = std::integral_constant<size_t, (I > J ? I : J)>;

  /**********************************************************************
  *  struct opnode							*
  **********************************************************************/
  class opnode;

  template <class E>
  using is_opnode = std::is_convertible<E, opnode>;
    
  //! 演算子のノードを表すクラス
  class opnode
  {
    protected:
      template <class ITER_, size_t SIZE_>
      static std::true_type	check_range(range<ITER_, SIZE_>)	;
      static std::false_type	check_range(...)			;

      template <class E_>
      using is_rangeobj	= decltype(check_range(std::declval<E_>()));

      template <class E_>
      using argument_t	= typename std::conditional<is_opnode<E_>::value ||
						    is_rangeobj<E_>::value,
						    const E_, const E_&>::type;
  };

  /**********************************************************************
  *  class unary_opnode<OP, E>						*
  **********************************************************************/
  //! 配列式に対する単項演算子を表すクラス
  /*!
    \param OP	各成分に適用される単項演算子の型
    \param E	単項演算子の引数となる式の型
  */
  template <class OP, class E>
  class unary_opnode : public opnode
  {
    public:
      using iterator = boost::transform_iterator<OP, const_iterator_t<E> >;

    public:
		unary_opnode(const E& expr, const OP& op)
		    :_expr(expr), _op(op)				{}

      constexpr static size_t
		size0()		{ return TU::size0<E>(); }
      iterator	begin()	const	{ return {std::begin(_expr), _op}; }
      iterator	end()	const	{ return {std::end(_expr),   _op}; }
      size_t	size()	const	{ return std::size(_expr); }

    private:
      argument_t<E>	_expr;
      const OP		_op;
  };
    
  template <class OP, class E> inline unary_opnode<OP, E>
  make_unary_opnode(const E& expr, const OP& op)
  {
      return {expr, op};
  }

  /**********************************************************************
  *  class binary_opnode<OP, L, R>					*
  **********************************************************************/
  //! 配列式に対する2項演算子を表すクラス
  /*!
    \param OP	各成分に適用される2項演算子の型
    \param L	2項演算子の第1引数となる式の型
    \param R	2項演算子の第2引数となる式の型
  */
  template <class OP, class L, class R>
  struct binary_opnode : public opnode
  {
    public:
      using iterator = transform_iterator2<OP, const_iterator_t<L>,
					       const_iterator_t<R> >;

    public:
		binary_opnode(const L& l, const R& r, const OP& op)
		    :_l(l), _r(r), _op(op)
		{
		    assert(std::size(_l) == std::size(_r));
		}

      constexpr static size_t
		size0()
	  	{
		    return max<TU::size0<L>(), TU::size0<R>()>::value;
		}
      iterator	begin()	const	{ return {std::begin(_l), std::begin(_r), _op};}
      iterator	end()	const	{ return {std::end(_l),   std::end(_r),   _op};}
      size_t	size()	const	{ return std::size(_l); }

    private:
      argument_t<L>	_l;
      argument_t<R>	_r;
      const OP		_op;
  };
    
  template <class OP, class L, class R> inline binary_opnode<OP, L, R>
  make_binary_opnode(const L& l, const R& r, const OP& op)
  {
      return {l, r, op};
  }

  template <class ASSIGN, class E> inline E&
  op_assign(E& expr, const element_t<E>& c)
  {
      for (auto& dst : expr)
	  ASSIGN()(dst, c);

      return expr;
  }

  template <class ASSIGN, class L, class R> inline L&
  op_assign(L& l, const R& r)
  {
      assert(std::size(l) == std::size(r));
      auto	src = std::begin(r);
      for (auto& dst : l)
      {
	  ASSIGN()(dst, *src);
	  ++src;
      }

      return l;
  }
}	// namespace detail

//! 与えられた式の各要素の符号を反転する.
/*!
  \param expr	式
  \return	符号反転演算子ノード
*/
template <class E, typename std::enable_if<is_range<E>::value>::type* = nullptr>
inline auto
operator -(const E& expr)
{
    return detail::make_unary_opnode(expr, negate());
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param expr	式
  \param c	乗数
  \return	乗算演算子ノード
*/
template <class E, typename std::enable_if<is_range<E>::value>::type* = nullptr>
inline auto
operator *(const E& expr, const element_t<E>& c)
{
    return detail::make_unary_opnode(expr, std::bind(multiplies(),
						     std::placeholders::_1, c));
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param c	乗数
  \param expr	式
  \return	乗算演算子ノード
*/
template <class E, typename std::enable_if<is_range<E>::value>::type* = nullptr>
inline auto
operator *(const element_t<E>& c, const E& expr)
{
    return detail::make_unary_opnode(expr, std::bind(multiplies(),
						     c, std::placeholders::_1));
}

//! 与えられた式の各要素を定数で割る.
/*!
  \param expr	式
  \param c	除数
  \return	除算演算子ノード
*/
template <class E, typename std::enable_if<is_range<E>::value>::type* = nullptr>
inline auto
operator /(const E& expr, const element_t<E>& c)
{
    return detail::make_unary_opnode(expr, std::bind(divides(),
						     std::placeholders::_1, c));
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param expr	式
  \param c	乗数
  \return	各要素にcが掛けられた結果の式
*/
template <class E>
inline typename std::enable_if<
    is_range<typename std::decay<E>::type>::value, E&>::type
operator *=(E&& expr, const element_t<typename std::decay<E>::type>& c)
{
    return detail::op_assign<multiplies_assign>(expr, c);
}

//! 与えられた式の各要素を定数で割る.
/*!
  \param expr	式
  \param c	除数
  \return	各要素がcで割られた結果の式
*/
template <class E>
inline typename std::enable_if<
    is_range<typename std::decay<E>::type>::value, E&>::type
operator /=(E&& expr, const element_t<typename std::decay<E>::type>& c)
{
    return detail::op_assign<divides_assign>(expr, c);
}

//! 与えられた2つの式の各要素の和をとる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	加算演算子ノード
*/
template <class L, class R,
	  typename std::enable_if<is_range<L>::value &&
				  is_range<R>::value>::type* = nullptr>
inline auto
operator +(const L& l, const R& r)
{
    return detail::make_binary_opnode(l, r, plus());
}

//! 与えられた2つの式の各要素の差をとる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	減算演算子ノード
*/
template <class L, class R,
	  typename std::enable_if<is_range<L>::value &&
				  is_range<R>::value>::type* = nullptr>
inline auto
operator -(const L& l, const R& r)
{
    return detail::make_binary_opnode(l, r, minus());
}

//! 与えられた左辺の式の各要素に右辺の式の各要素を加える.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	各要素が加算された左辺の式
*/
template <class L, class R>
inline typename std::enable_if<is_range<typename std::decay<L>::type>::value &&
			       is_range<R>::value, L&>::type
operator +=(L&& l, const R& r)
{
    return detail::op_assign<plus_assign>(l, r);
}

//! 与えられた左辺の式の各要素から右辺の式の各要素を減じる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	各要素が減じられた左辺の式
*/
template <class L, class R>
inline typename std::enable_if<is_range<typename std::decay<L>::type>::value &&
			       is_range<R>::value, L&>::type
operator -=(L&& l, const R& r)
{
    return detail::op_assign<minus_assign>(l, r);
}

/************************************************************************
*  generic algorithms for ranges					*
************************************************************************/
template <class E, class T>
typename std::enable_if<!is_range<typename std::decay<E>::type>::value>::type
fill(E&& expr, const T& val)
{
    expr = val;
}
template <class E, class T>
typename std::enable_if<is_range<typename std::decay<E>::type>::value>::type
fill(E&& expr, const T& val)
{
    for (auto iter = std::begin(expr); iter != std::end(expr); ++iter)
	fill(*iter, val);
}

//! 与えられた2次元配列式の転置を返す
/*
  \param expr	2次元配列式
  \return	expr を転置した2次元配列式
 */ 
template <class E,
	  typename std::enable_if<is_range2<E>::value>::type* = nullptr>
inline auto
transpose(const E& expr)
{
    constexpr size_t	siz0 = size0<value_t<E> >();
    return make_range<siz0>(column_begin(expr), size<1>(expr));
}

//! 与えられた式の各要素の自乗和を求める.
/*!
  \param x	式
  \return	式の各要素の自乗和
*/
template <class E,
	  typename std::enable_if<is_range<E>::value>::type* = nullptr>
inline auto
square(const E& expr)
{
    return square<size0<E>()>(std::begin(expr), std::size(expr));
}

//! 与えられた式の各要素の自乗和の平方根を求める.
/*!
  \param x	式
  \return	式の各要素の自乗和の平方根
*/
template <class T> inline auto
length(const T& x)
{
    return std::sqrt(square(x));
}
    
//! 与えられた二つの式の各要素の差の自乗和を求める.
/*!
  \param x	第1の式
  \param y	第2の式
  \return	xとyの各要素の差の自乗和
*/
template <class L, class R> inline auto
sqdist(const L& x, const R& y)
{
    return square(x - y);
}
    
//! 与えられた二つの式の各要素の差の自乗和の平方根を求める.
/*!
  \param x	第1の式
  \param y	第2の式
  \return	xとyの各要素の差の自乗和の平方根
*/
template <class L, class R> inline auto
dist(const L& x, const R& y)
{
    return std::sqrt(sqdist(x, y));
}

}	// namespace TU
#endif	// !__TU_RANGE_H
