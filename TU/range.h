/*!
  \file		range.h
  \author	Toshio UESHIBA
  \brief	2つの反復子によって指定されるレンジの定義と実装
*/
#ifndef TU_RANGE_H
#define TU_RANGE_H

#include <cassert>
#include <initializer_list>
#include "TU/iterator.h"
#include "TU/algorithm.h"

namespace TU
{
/************************************************************************
*  type aliases: iterator_t<E>, reverse_iterator_t<E>, value_t<E>	*
*		 element_t<E>, const_iterator_t<ITER>			*
************************************************************************/
namespace detail
{
  //! 式に適用できる反復子の型を返す
  /*!
    E, const E&, E&, E&&型の x に begin() が
    適用できるかチェックし，可能な場合はその型を返す．
  */
  template <class E>
  auto	iterator_t(E&& x) -> decltype(begin(x))				;
  void	iterator_t(...)							;
}	// namespace detail

//! 式に適用できる反復子の型を返す
/*!
  \param E	式またはそれへの参照の型
  \return	E に反復子が適用できればその型，適用できなければ void
*/
template <class E>
using iterator_t = decltype(detail::iterator_t(std::declval<E>()));

//! 反復子が指す式に適用できる反復子の型を返す
/*!
  \param ITER	反復子
  \return	ITER が指す式に反復子が適用できればその型，適用できなければ void
*/
template <class ITER>
using iterator_iterator = iterator_t<iterator_value<ITER> >;

//! 式が持つ逆反復子の型を返す
/*!
  反復子を持たない式を与えるとコンパイルエラーとなる.
  \param E	式またはそれへの参照の型
  \return	E が持つ逆反復子の型
*/
template <class E>
using reverse_iterator_t = std::reverse_iterator<iterator_t<E> >;

//! 式が持つ反復子が指す型を返す
/*!
  反復子を持たない式を与えるとコンパイルエラーとなる.
  \param E	反復子を持つ式またはそれへの参照の型
  \return	E の反復子が指す型
*/
template <class E>
using value_t	= iterator_value<iterator_t<E> >;

namespace detail
{
  template <class T>
  struct identity
  {
      using type = T;
  };

  template <class E, class=TU::iterator_t<E> >
  struct element_t
  {
      using type = typename element_t<value_t<E> >::type;
  };
  template <class E>
  struct element_t<E, void> : identity<E>				{};
}	// namespace detail
    
//! 式が持つ反復子が指す型を再帰的に辿って到達する型を返す
/*!
  \param E	式またはそれへの参照の型
  \return	E が反復子を持てばそれが指す型を再帰的に辿って到達する型，
		持たなければ E 自身
*/
template <class E>
using element_t	= typename detail::element_t<E>::type;

namespace detail
{
  template <class ITER>
  struct const_iterator_t
  {
      using type = ITER;
  };
  template <class T>
  struct const_iterator_t<T*>
  {
      using type = const T*;
  };
  template <class... ITER>
  struct const_iterator_t<std::tuple<ITER...> >
  {
      using type = std::tuple<typename const_iterator_t<ITER>::type...>;
  };
}	// namespace detail

//! 与えられた反復子のconstバージョンを返す
/*!
  \param ITER	反復子の型
  \return	ITER のconstバージョン
*/
template <class ITER>
using const_iterator_t = typename detail::const_iterator_t<ITER>::type;
    
/************************************************************************
*  rank<E>()								*
************************************************************************/
namespace detail
{
  template <class E>
  auto	check_stdbegin(E&& x) -> decltype(std::begin(x),
					  std::true_type())		;
  auto	check_stdbegin(...)   -> std::false_type			;

  template <class E>
  using has_stdbegin = decltype(check_stdbegin(std::declval<E>()));
}	// namespace detail

//! 式の次元数(軸の個数)を返す
/*!
  \param E	式の型
  \return	式の次元数
*/
template <class E>
constexpr std::enable_if_t<!detail::has_stdbegin<E>::value, size_t>
rank()
{
    return 0;
}
template <class E>
constexpr std::enable_if_t<detail::has_stdbegin<E>::value, size_t>
rank()
{
    return 1 + rank<value_t<E> >();
}
    
/************************************************************************
*  size0<E>()								*
************************************************************************/
namespace detail
{
  template <class E>
  auto		  has_size0(const E&)
			-> decltype(E::size0(), std::true_type())	;
  std::false_type has_size0(...)					;
}	// namespace detail

template <class E>
using has_size0 = decltype(detail::has_size0(std::declval<E>()));
    
//! 式の最上位軸の要素数を返す
/*!
  本関数で返されるのは静的なサイズであり，式が可変個の要素を持つ場合，0が返される
  \param E	式またはそれへの参照の型
  \return	最上位軸の静的な要素数(可変サイズの場合0)
 */
template <class E> constexpr std::enable_if_t<!has_size0<E>::value, size_t>
size0()
{
  // Array2<T, R, C> * Array<T, 0, 0> の評価結果の型が Array<T, R, 0>
  // ではなく Array<T, 0, 0> になるために，0ではなく1を返すことが必要．
    return 1;
}
template <class E> constexpr std::enable_if_t<has_size0<E>::value, size_t>
size0()
{
    return std::decay_t<E>::size0();
}
    
/************************************************************************
*  size<I>(const E&)							*
************************************************************************/
namespace detail
{
  /*
   *  A ^ b において演算子ノード product_opnode<bit_xor, L, R> が
   *  生成されるが，これを評価して2次元配列に代入する際に代入先の領域確保のため
   *  size<0>(opnode), size<1>(opnode) が呼び出される．後者はopnodeの反復子が
   *  指す先を評価するが，これは2つの3次元ベクトルのベクトル積を評価することに
   *  等しく，コストが高い処理である．そこで，ベクトル積を評価せずその評価結果の
   *  サイズだけを得るために，以下のオーバーロードを導入する．
   */
  template <class OP, class L, class R>	class product_opnode;
  struct bit_xor;
    
  template <class L, class R> constexpr size_t
  size(const product_opnode<bit_xor, L, R>&, std::integral_constant<size_t, 1>)
  {
      return 3;
  }

  template <class E> inline auto
  size(const E& expr, std::integral_constant<size_t, 0>)
  {
      return std::size(expr);
  }
  template <size_t I, class E> inline auto
  size(const E& expr, std::integral_constant<size_t, I>)
  {
      return size(*std::begin(expr), std::integral_constant<size_t, I-1>());
  }
}	// namespace detail

//! 与えられた式について，指定された軸の要素数を返す
/*!
  \param I	軸を指定するindex (0 <= I < Dimension)
  \param E	式の型
  \return	軸Iの要素数
 */
template <size_t I, class E, std::enable_if_t<rank<E>() != 0>* = nullptr>
inline auto
size(const E& expr)
{
    return detail::size(expr, std::integral_constant<size_t, I>());
}

/************************************************************************
*  manipulator print_sizes() and print_sizes_and_strides()		*
************************************************************************/
template <class E>
class sizes_holder
{
  public:
			sizes_holder(const E& expr)	:_expr(expr)	{}

    std::ostream&	operator ()(std::ostream& out) const
			{
			    return print_size(out, &_expr);
			}
    
  protected:
    template <class ITER_> constexpr
    static size_t	rank()
			{
			    using value_type = iterator_value<ITER_>;
			    return TU::rank<value_type>();
			}
    template <class ITER_>
    static std::enable_if_t<rank<ITER_>() == 0, std::ostream&>
			print_x(std::ostream& out, ITER_)
			{
			    return out;
			}
    template <class ITER_>
    static std::enable_if_t<rank<ITER_>() != 0, std::ostream&>
			print_x(std::ostream& out, ITER_)
			{
			    return out << 'x';
			}

  private:
    template <class ITER_>
    static std::enable_if_t<rank<ITER_>() == 0, std::ostream&>
			print_size(std::ostream& out, ITER_)
			{
			    return out;
			}
    template <class ITER_>
    static std::enable_if_t<rank<ITER_>() != 0, std::ostream&>
			print_size(std::ostream& out, ITER_ iter)
			{
			    using	std::cbegin;
			    using	std::size;
			    
			    const auto&	val = *iter;
			    return print_size(print_x(out << size(val),
						      cbegin(val)),
					      cbegin(val));
			}

  protected:
    const E&	_expr;
};
    
template <class E>
class sizes_and_strides_holder : public sizes_holder<E>
{
  private:
    using super	= sizes_holder<E>;
    
  public:
			sizes_and_strides_holder(const E& expr)
			    :super(expr)				{}

    std::ostream&	operator ()(std::ostream& out) const
			{
			    return print_stride(super::operator ()(out) << ':',
						std::cbegin(super::_expr));
			}
    
  private:
    template <class ITER_>
    static std::enable_if_t<super::template rank<ITER_>() == 0, std::ostream&>
			print_stride(std::ostream& out, ITER_)
			{
			    return out;
			}
    template <class ITER_>
    static std::enable_if_t<super::template rank<ITER_>() != 0, std::ostream&>
			print_stride(std::ostream& out, ITER_ iter)
			{
			    const auto&	val = *iter;
			    return print_stride(super::print_x(
						    out << iter.stride(),
						    std::cbegin(val)),
						std::cbegin(val));
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
//! 先頭要素を指す反復子と要素数によって指定される範囲(レンジ)を表すクラス
/*!
  メンバ変数(_iter, _size)が全てconstなので代入演算子も含め全メンバ関数を
  constにすることも可能だが，非constな反復子を返す begin(), end() は
  敢えて非constメンバ関数にしている．これは，std::cbegin(const T& x) の型が
  decltype(x.begin()) であって decltype(x.cbegin()) でないため，非constな
  rangeに対してstd::cbegin(range)を適用しても非constな反復子が返されて
  しまうことを防ぐためである．
  \param ITER	反復子の型
  \param SIZE	レンジサイズ(0ならば可変長)
*/
template <class ITER, size_t SIZE=0>
class range
{
  public:
    using value_type	 = iterator_value<ITER>;
    using const_iterator = const_iterator_t<ITER>;
    
  public:
		range(ITER begin)		:_begin(begin)	{}
		range(ITER begin, size_t)	:_begin(begin)	{}
    
		range()						= delete;
		range(const range&)				= default;
    range&	operator =(const range& r)
		{
		    copy<SIZE>(r.begin(), SIZE, _begin);
		    return *this;
		}
		range(range&&)					= default;
    range&	operator =(range&&)				= default;
    
    template <class ITER_,
	      std::enable_if_t<std::is_convertible<ITER_, ITER>::value>*
	      = nullptr>
		range(const range<ITER_, SIZE>& r)
		    :_begin(r.begin())
		{
		}

    template <class E_> std::enable_if_t<rank<E_>() != 0, range&>
		operator =(const E_& expr)
		{
		    using	TU::begin;
		    using	std::size;
		    
		    assert(size(expr) == SIZE);
		    copy<SIZE>(begin(expr), SIZE, _begin);
		    return *this;
		}

		range(std::initializer_list<value_type> args)
		    :_begin(const_cast<value_type*>(args.begin()))
    		{
		    assert(args.size() == size());
		}
    range&	operator =(std::initializer_list<value_type> args)
		{
		    assert(args.size() == SIZE);
		  // initializer_list<T> はalignmentされないので，
		  // SIMD命令が使われぬようcopy<N>()は使用しない．
		    std::copy_n(args.begin(), SIZE, _begin);
		    return *this;
		}

    template <class T_> std::enable_if_t<rank<T_>() == 0, range&>
		operator =(const T_& c)
		{
		    fill<SIZE>(_begin, SIZE, c);
		    return *this;
		}

    constexpr static
    size_t	size0()		{ return SIZE; }
    constexpr static
    size_t	size()		{ return SIZE; }
    auto	begin()		{ return _begin; }
    auto	end()		{ return _begin + SIZE; }
    auto	begin()   const	{ return const_iterator(_begin); }
    auto	end()	  const	{ return const_iterator(_begin + SIZE); }
    auto	cbegin()  const	{ return begin(); }
    auto	cend()	  const	{ return end(); }
    auto	rbegin()  	{ return std::make_reverse_iterator(end()); }
    auto	rend()	  	{ return std::make_reverse_iterator(begin()); }
    auto	rbegin()  const { return std::make_reverse_iterator(end()); }
    auto	rend()	  const	{ return std::make_reverse_iterator(begin()); }
    auto	crbegin() const	{ return rbegin(); }
    auto	crend()	  const	{ return rend(); }
    decltype(auto)
		operator [](size_t i) const
		{
		    assert(i < size());
		    return *(_begin + i);
		}
    
  private:
    const ITER	_begin;
};

//! 可変長レンジ
/*!
  \param ITER	反復子の型
*/
template <class ITER>
class range<ITER, 0>
{
  public:
    using value_type	 = iterator_value<ITER>;
    using const_iterator = const_iterator_t<ITER>;
    
  public:
		range(ITER begin, size_t size)
		    :_begin(begin), _size(size)			{}
    
		range()						= delete;
		range(const range&)				= default;
    range&	operator =(const range& r)
		{
		    assert(r.size() == size());
		    copy<0>(r._begin, _size, _begin);
		    return *this;
		}
		range(range&&)					= default;
    range&	operator =(range&&)				= default;
    
    template <class ITER_,
	      std::enable_if_t<std::is_convertible<ITER_, ITER>::value>*
	      = nullptr>
		range(const range<ITER_, 0>& r)
		    :_begin(r.begin()), _size(r.size())
		{
		}

    template <class E_> std::enable_if_t<rank<E_>() != 0, range&>
		operator =(const E_& expr)
		{
		    using	TU::begin;
		    using	std::size;
		    
		    assert(size(expr) == _size);
		    copy<TU::size0<E_>()>(begin(expr), _size, _begin);
		    return *this;
		}
		
		range(std::initializer_list<value_type> args)
		    :_begin(const_cast<value_type*>(args.begin())),
		     _size(args.size())
    		{
		}
    range&	operator =(std::initializer_list<value_type> args)
		{
		    assert(args.size() == _size);
		  // initializer_list<T> はalignmentされないので，
		  // SIMD命令が使われぬようcopy<N>()は使用しない．
		    std::copy_n(args.begin(), _size, _begin);
		    return *this;
		}
		
    template <class T_> std::enable_if_t<rank<T_>() == 0, range&>
		operator =(const T_& c)
		{
		    fill<0>(_begin, _size, c);
		    return *this;
		}

    constexpr static
    size_t	size0()		{ return 0; }
    size_t	size()	  const	{ return _size; }
    auto	begin()		{ return _begin; }
    auto	end()		{ return _begin + _size; }
    auto	begin()	  const	{ return const_iterator(_begin); }
    auto	end()	  const	{ return const_iterator(_begin + _size); }
    auto	cbegin()  const	{ return begin(); }
    auto	cend()	  const	{ return end(); }
    auto	rbegin()  	{ return std::make_reverse_iterator(end()); }
    auto	rend()	  	{ return std::make_reverse_iterator(begin()); }
    auto	rbegin()  const	{ return std::make_reverse_iterator(end()); }
    auto	rend()	  const	{ return std::make_reverse_iterator(begin()); }
    auto	crbegin() const	{ return rbegin(); }
    auto	crend()	  const	{ return rend(); }
    decltype(auto)
		operator [](size_t i) const
		{
		    assert(i < size());
		    return *(_begin + i);
		}

  private:
    const ITER		_begin;
    const size_t	_size;
};

//! 固定長レンジを生成する
/*!
  \param SIZE	レンジサイズ
  \param iter	レンジの先頭要素を指す反復子
  \return	固定長レンジ
*/
template <size_t SIZE, class ITER> inline range<ITER, SIZE>
make_range(ITER iter)
{
    return {iter};
}
    
//! 可変長レンジを生成する
/*!
  \param iter	レンジの先頭要素を指す反復子
  \param size	レンジサイズ
  \return	可変長レンジ
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
  //for (const auto& val : r)
  //	out << ' ' << val;
    for (auto iter = r.begin(); iter != r.end(); ++iter)
	out << ' ' << *iter;
    return out << std::endl;
}
    
/************************************************************************
*  type alias: iterator_stride<ITER>					*
************************************************************************/
template <class T> auto
stride(T*) -> ptrdiff_t							;

template <class ITER> auto
stride(const ITER& iter) -> decltype(stride(iter.base()))		;

template <class... ITER> auto
stride(const std::tuple<ITER...>&)
    -> std::tuple<decltype(stride(std::declval<ITER>()))...>		;

template <class ITER_TUPLE> auto
stride(const zip_iterator<ITER_TUPLE>& iter)
    -> decltype(stride(iter.get_iterator_tuple()))			;

template <class ITER>
using iterator_stride = decltype(stride(std::declval<ITER>()));

/************************************************************************
*  advance_stride(ITER&, STRIDE)					*
************************************************************************/
template <class ITER> inline auto
advance_stride(ITER& iter, ptrdiff_t stride)
    -> void_t<decltype(iter += stride)>
{
    iter += stride;
}

template <class ITER, class... STRIDE> inline auto
advance_stride(ITER& iter, const std::tuple<STRIDE...>& stride)
    -> void_t<decltype(iter.get_iterator_tuple())>
{
    using tuple_t = std::decay_t<decltype(iter.get_iterator_tuple())>;
    
    tuple_for_each([](auto&& x, const auto& y){ advance_stride(x, y); },
		   const_cast<tuple_t&>(iter.get_iterator_tuple()), stride);
}

template <class ITER, class... STRIDE> inline auto
advance_stride(ITER& iter, const std::tuple<STRIDE...>& stride)
    -> void_t<decltype(iter.base())>
{
    using base_t = std::decay_t<decltype(iter.base())>;
    
    advance_stride(const_cast<base_t&>(iter.base()), stride);
}

/************************************************************************
*  class range_iterator<ITER, STRIDE, SIZE>				*
************************************************************************/
namespace detail
{
  template <class DIFF, ptrdiff_t STRIDE, size_t SIZE>
  struct stride_and_size
  {
      using stride_t	= ptrdiff_t;
      
      stride_and_size(stride_t, size_t)		{}
      constexpr static stride_t	stride()	{ return STRIDE; }
      constexpr static size_t	size()		{ return SIZE; }
  };
  template <class DIFF, size_t SIZE>
  struct stride_and_size<DIFF, 0, SIZE>
  {
      using stride_t	= DIFF;
      
      stride_and_size(stride_t stride, size_t)
	  :_stride(stride)			{}
      stride_t			stride() const	{ return _stride; }
      constexpr static size_t	size()		{ return SIZE; }

    private:
      stride_t	_stride;
  };
  template <class DIFF, ptrdiff_t STRIDE>
  struct stride_and_size<DIFF, STRIDE, 0>
  {
      using stride_t	= ptrdiff_t;
      
      stride_and_size(stride_t, size_t size)
	  :_size(size)				{}
      constexpr static stride_t	stride()	{ return STRIDE; }
      size_t			size()	 const	{ return _size; }

    private:
      size_t	_size;
  };
  template <class DIFF>
  struct stride_and_size<DIFF, 0, 0>
  {
      using stride_t	= DIFF;
      
      stride_and_size(stride_t stride, size_t size)
	  :_stride(stride), _size(size)		{}
      stride_t			stride() const	{ return _stride; }
      size_t			size()	 const	{ return _size; }

    private:
      stride_t	_stride;
      size_t	_size;
  };
}	// namespace detail

//! 配列を一定間隔に切り分けたレンジを指す反復子
/*!
  \param ITER	配列の要素を指す反復子の型
  \param STRIDE	インクリメントしたときに進める要素数(0ならば可変)
  \param SIZE	レンジサイズ(0ならば可変長)
*/
template <class ITER, ptrdiff_t STRIDE, size_t SIZE>
class range_iterator
    : public boost::iterator_adaptor<range_iterator<ITER, STRIDE, SIZE>,
				     ITER,
				     range<ITER, SIZE>,
				     boost::use_default,
				     range<ITER, SIZE> >,
      public detail::stride_and_size<iterator_stride<ITER>, STRIDE, SIZE>
{
  private:
    using super	= boost::iterator_adaptor<range_iterator,
					  ITER,
					  range<ITER, SIZE>,
					  boost::use_default,
					  range<ITER, SIZE> >;
    using ss	= detail::stride_and_size<iterator_stride<ITER>,
					  STRIDE, SIZE>;
    
  public:
    using		typename super::reference;
    using		typename super::difference_type;
    using stride_t    =	typename ss::stride_t;
    friend class	boost::iterator_core_access;
	  
  public:
		range_iterator(ITER iter,
			       iterator_stride<ITER> stride=STRIDE,
			       size_t size=SIZE)
		    :super(iter), ss(stride, size)			{}

    template <class ITER_,
	      std::enable_if_t<std::is_convertible<ITER_, ITER>::value>*
	      = nullptr>
		range_iterator(
		    const range_iterator<ITER_, STRIDE, SIZE>& iter)
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
		    advance_stride(super::base_reference(), stride());
		}
    void	decrement()
		{
		    advance_stride(super::base_reference(), -stride());
		}
    void	advance(difference_type n)
		{
		    advance_stride(super::base_reference(), n*stride());
		}
    difference_type
		distance_to(const range_iterator& iter) const
		{
		    return (iter.base() - super::base()) / leftmost(stride());
		}
    static auto	leftmost(ptrdiff_t stride) -> ptrdiff_t
		{
		    return stride;
		}
    template <class STRIDE_>
    static auto	leftmost(const STRIDE_& stride)
		{
		    using	std::get;
		    
		    return leftmost(get<0>(stride));
		}
};
	
/************************************************************************
*  TU::stride(const ITER&)						*
************************************************************************/
//! 反復子が指すレンジが所属する軸のストライドを返す
/*!
  \param iter	レンジを指す反復子
  \return	レンジ軸のストライド
*/
template <class ITER, ptrdiff_t STRIDE, size_t SIZE> inline auto
stride(const range_iterator<ITER, STRIDE, SIZE>& iter)
{
    return iter.stride();
}

template <class ITER> auto
stride(const ITER& iter) -> decltype(stride(iter.base()))
{
    return stride(iter.base());
}

template <class... ITER> inline auto
stride(const std::tuple<ITER...>& iter_tuple)
    -> std::tuple<decltype(stride(std::declval<ITER>()))...>
{
    return tuple_transform([](const auto& iter){ return stride(iter); },
			   iter_tuple);
}

template <class ITER_TUPLE> inline auto
stride(const zip_iterator<ITER_TUPLE>& iter)
    -> decltype(stride(iter.get_iterator_tuple()))
{
    return stride(iter.get_iterator_tuple());
}

template <class ITER0, class ITER1, class... ITERS> inline auto
stride(const ITER0& iter0, const ITER1& iter1, const ITERS&... iters)
{
    return stride(std::make_tuple(iter0, iter1, iters...));
}

/************************************************************************
*  fixed size & fixed stride ranges and associated iterators		*
************************************************************************/
//! 固定長レンジを指し，インクリメント時に固定した要素数だけ進める反復子を生成する
/*!
  \param STRIDE	インクリメント時に進める要素数
  \param SIZE	レンジサイズ
  \param iter	レンジの先頭要素を指す反復子
  \return	レンジ反復子
*/
template <ptrdiff_t STRIDE, size_t SIZE, class ITER>
inline range_iterator<ITER, STRIDE, SIZE>
make_range_iterator(ITER iter)
{
    return {iter};
}

//! 多次元固定長レンジを指し，インクリメント時に固定したブロック数だけ進める反復子を生成する
/*!
  \param STRIDE	インクリメント時に進める最上位軸のブロック数
  \param SIZE	最上位軸のレンジサイズ
  \param SS	2番目以降の軸の{ストライド, レンジサイズ}の並び
  \param iter	レンジの先頭要素を指す反復子
*/
template <size_t STRIDE, size_t SIZE, size_t... SS, class ITER,
	  std::enable_if_t<sizeof...(SS) != 0>* = nullptr>
inline auto
make_range_iterator(ITER iter)
{
    return make_range_iterator<STRIDE, SIZE>(make_range_iterator<SS...>(iter));
}

template <size_t SIZE, size_t... SS, class ITER,
	  std::enable_if_t<sizeof...(SS) != 0>* = nullptr>
inline auto
make_range(ITER iter)
{
    return make_range<SIZE>(make_range_iterator<SS...>(iter));
}

/************************************************************************
*  fixed size & variable stride ranges and associated iterators		*
************************************************************************/
//! 固定長レンジを指し，インクリメント時に指定した要素数だけ進める反復子を生成する
/*!
  \param SIZE	レンジサイズ
  \param iter	レンジの先頭要素を指す反復子
  \param stride	インクリメント時に進める要素数
  \return	レンジ反復子
*/
template <size_t SIZE, class ITER>
inline range_iterator<ITER, 0, SIZE>
make_range_iterator(ITER iter, iterator_stride<ITER> stride)
{
    return {iter, stride};
}
    
//! 多次元固定長レンジを指し，インクリメント時に指定したブロック数だけ進める反復子を生成する
/*!
  \param SIZE		最上位軸のレンジサイズ
  \param SIZES		2番目以降の軸のレンジサイズの並び
  \param stride		最上位軸のストライド
  \param strides	2番目以降の軸のストライドの並び
  \param iter		レンジの先頭要素を指す反復子
  \return		レンジ反復子
*/
template <size_t SIZE, size_t... SIZES, class ITER, class... STRIDES,
	  std::enable_if_t<sizeof...(SIZES) == sizeof...(STRIDES)>* = nullptr>
inline auto
make_range_iterator(ITER iter,
		    iterator_stride<ITER> stride, STRIDES... strides)
{
    return make_range_iterator<SIZE>(
	       make_range_iterator<SIZES...>(iter, strides...), stride);
}

template <size_t SIZE, size_t... SIZES, class ITER, class... STRIDES,
	  std::enable_if_t<sizeof...(SIZES) == sizeof...(STRIDES)>* = nullptr>
inline auto
make_range(ITER iter, STRIDES... strides)
{
    return make_range<SIZE>(make_range_iterator<SIZES...>(iter, strides...));
}

/************************************************************************
*  variable size & variable stride ranges and associated iterators	*
************************************************************************/
//! 指定された長さのレンジを指し，インクリメント時に指定した要素数だけ進める反復子を生成する
/*!
  \param iter	レンジの先頭要素を指す反復子
  \param stride	インクリメント時に進める要素数
  \param size	レンジサイズ
  \return	レンジ反復子
*/
template <class ITER> inline range_iterator<ITER, 0, 0>
make_range_iterator(ITER iter, iterator_stride<ITER> stride, size_t size)
{
    return {iter, stride, size};
}
    
//! 多次元可変長レンジを指し，インクリメント時に指定したブロック数だけ進める反復子を生成する
/*!
  \param iter		レンジの先頭要素を指す反復子
  \param stride		最上位軸のストライド
  \param size		最上位軸のレンジサイズ
  \param ss		2番目以降の軸の{ストライド, レンジサイズ}の並び
  \return		レンジ反復子
*/
template <class ITER, class... SS> inline auto
make_range_iterator(ITER iter,
		    iterator_stride<ITER> stride, size_t size, SS... ss)
{
    return make_range_iterator(make_range_iterator(iter, ss...),
			       stride, size);
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
namespace detail
{
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
}	// namespace detail

template <class ITER, class... SIZES> inline auto
make_dense_range(ITER iter, size_t size, SIZES... sizes)
{
    return make_range(detail::make_dense_range_iterator(iter, sizes...), size);
}

/************************************************************************
*  slice extraction							*
************************************************************************/
namespace detail
{
  template <class ITER> inline ITER
  make_slice_iterator(ITER iter)
  {
      return iter;
  }

  template <class ITER, class... IS> inline auto
  make_slice_iterator(ITER iter, size_t idx, size_t size, IS... is)
  {
      return make_range_iterator(make_slice_iterator(
				     begin(*iter) + idx, is...),
				 stride(iter), size);
  }

  template <size_t SIZE, size_t... SIZES, class ITER, class... INDICES>
  inline auto
  make_slice_iterator(ITER iter, size_t idx, INDICES... indices)
  {
      return make_range_iterator<SIZE>(make_slice_iterator<SIZES...>(
					   begin(*iter) + idx, indices...),
				       stride(iter));
  }
}	// namespace detail
    
template <class RANGE, class... IS> inline auto
slice(RANGE&& r, size_t idx, size_t size, IS... is)
{
    return make_range(detail::make_slice_iterator(begin(r) + idx, is...),
		      size);
}

template <size_t SIZE, size_t... SIZES, class RANGE, class... INDICES>
inline auto
slice(RANGE&& r, size_t idx, INDICES... indices)
{
    return make_range<SIZE>(detail::make_slice_iterator<SIZES...>(
				begin(r) + idx, indices...));
}

/************************************************************************
*  pipeline stuff							*
************************************************************************/
namespace detail
{
  template <class T, bool MASK, class FUNC>
  class mapped_tag
  {
    private:
      template <class FUNC_>
      static typename FUNC_::argument_type::element_type
		element_t(FUNC_)					;
      template <class FUNC_>
      static typename FUNC_::first_argument_type::element_type
		element_t(FUNC_)					;
      static T	element_t(...)						;
      
    public:
      using element_type = decltype(element_t(
					std::declval<std::decay_t<FUNC> >()));

    public:
      mapped_tag(FUNC&& func)	:_func(std::forward<FUNC>(func))	{}

      FUNC&&	functor() 	{ return std::forward<FUNC>(_func); }

    private:
      FUNC	_func;
  };
}	// namespace detail

template <class T=void, bool MASK=false, class FUNC>
inline detail::mapped_tag<T, MASK, FUNC>
mapped(FUNC&& func)
{
    return {std::forward<FUNC>(func)};
}

template <class... ARG> inline std::tuple<const ARG&...>
zip(const ARG&... x)
{
    return std::tie(x...);
}
    
template <class T, bool MASK, class FUNC, class ITER,
	  std::enable_if_t<std::is_void<iterator_iterator<ITER> >::value>*
	  = nullptr>
inline auto
make_map_iterator(detail::mapped_tag<T, MASK, FUNC>&& m, const ITER& iter)
{
    return make_map_iterator(m.functor(), iter);
}

template <class T, bool MASK, class FUNC, class ITER,
	  std::enable_if_t<!std::is_void<iterator_iterator<ITER> >::value>*
	  = nullptr>
inline auto
make_map_iterator(detail::mapped_tag<T, MASK, FUNC>&& m, const ITER& iter)
{
    using	std::size;
    
    return make_range_iterator(make_map_iterator(std::move(m), begin(*iter)),
			       stride(iter), size(*iter));
}
    
template <class ARG, class T, bool MASK, class FUNC> inline auto
operator >>(const ARG& x, detail::mapped_tag<T, MASK, FUNC>&& m)
{
    using	std::size;
    
    return make_range(make_map_iterator(std::move(m), begin(x)), size(x));
}

/************************************************************************
*  class column_iterator<ROW, NROWS>					*
************************************************************************/
//! 2次元配列の列を指す反復子
/*!
  \param ROW	begin(), end()をサポートするコンテナを指す反復子の型
  \param NROWS	begin()とend()間の距離(0ならば可変長)
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
		    return {{{super::base()}, _row}, _nrows};
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
column_begin(E&& expr)
{
    using	std::begin;
    using	std::size;
    
    constexpr auto	N = size0<E>();
    
    return make_column_iterator<N>(begin(expr), size(expr), 0);
}

template <class E> inline auto
column_cbegin(const E& expr)
{
    return column_begin(expr);
}
    
template <class E> inline auto
column_end(E&& expr)
{
    using	std::begin;
    using	std::size;
    
    constexpr auto	N = size0<E>();
    
    return make_column_iterator<N>(begin(expr), size(expr), size<1>(expr));
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
  template <size_t I, size_t... J>
  struct max
  {
      constexpr static size_t	value = (I > max<J...>::value ?
					 I : max<J...>::value);
  };
  template <size_t I>
  struct max<I>
  {
      constexpr static size_t	value = I;
  };

  template <class ITER, size_t SIZE>
  std::true_type	check_range(const range<ITER, SIZE>&)		;
  std::false_type	check_range(...)				;

  template <class E>
  using is_range = decltype(check_range(std::declval<E>()));

  /**********************************************************************
  *  struct opnode<NODE>						*
  **********************************************************************/
  //! 演算子のノードを表すクラス
  template <class NODE> struct opnode
  {
      constexpr static auto
		size0()
		{
		    return NODE::size0();
		}
      auto	begin()	const
		{
		    return static_cast<const NODE*>(this)->begin();
		}
      auto	end() const
		{
		    return static_cast<const NODE*>(this)->end();
		}
      auto	size() const
		{
		    return static_cast<const NODE*>(this)->size();
		}
  };

  template <class NODE> inline std::ostream&
  operator <<(std::ostream& out, const opnode<NODE>& expr)
  {
      for (const auto& x : expr)
	  out << ' ' << x;
      return out << std::endl;
  }

  /**********************************************************************
  *  class generic_opnode<OP, E...>					*
  **********************************************************************/
  //! 配列式に対する演算子を表すクラス
  /*!
    与えられた配列式が一時オブジェクトの場合は，その内容を本クラスオブジェクト
    内に保持した実体にmoveする．そうでない場合は，それへの参照を本クラスオブ
    ジェクト内に保持する．
    \param OP	各成分に適用される演算子の型
    \param E	演算子の引数となる式または式への参照の型
  */
  template <class OP, class E0, class... E>
  class generic_opnode : public opnode<generic_opnode<OP, E0, E...> >
  {
    private:
      using	expr_t = std::conditional_t<sizeof...(E),
					    std::tuple<E0, E...>, E0>;
      
    public:
		generic_opnode(OP&& op, E0&& expr0, E&&... expr)
		    :_op(std::forward<OP>(op)),
		     _expr(std::forward<E0>(expr0), std::forward<E>(expr)...)
		{}

      constexpr static auto
		size0()
		{
		    return max<TU::size0<E0>(), TU::size0<E>()...>::value;
		}
      auto	begin()	const
		{
		    using	TU::begin;

		    return make_map_iterator(_op, begin(_expr));
		}
      auto	end() const
		{
		    using	TU::end;

		    return make_map_iterator(_op, end(_expr));
		}
      auto	size() const
		{
		    using	std::size;
		    
		    return size(_expr);
		}
      decltype(auto)
		operator [](size_t i) const
		{
		    assert(i < size());
		    return *(begin() + i);
		}
      
    private:
      const OP		_op;
      const expr_t	_expr;
  };
    
  template <class OP, class... E> inline generic_opnode<OP, E...>
  make_generic_opnode(OP&& op, E&&... expr)
  {
    // exprの実引数がX&&型(X型の一時オブジェクト)ならば E = X,
    // そうでなければ E = X& または E = const X& となる．
      return {std::forward<OP>(op), std::forward<E>(expr)...};
  }

}	// namespace detail

//! 与えられた式の各要素の符号を反転する.
/*!
  \param expr	式
  \return	符号反転演算子ノード
*/
template <class E, std::enable_if_t<rank<E>() != 0>* = nullptr> inline auto
operator -(E&& expr)
{
    return detail::make_generic_opnode(
		[](auto&& x){ return -std::forward<decltype(x)>(x); },
		std::forward<E>(expr));
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param expr	式
  \param c	乗数
  \return	乗算演算子ノード
*/
template <class E, std::enable_if_t<rank<E>() != 0>* = nullptr> inline auto
operator *(E&& expr, element_t<E> c)
{
    return detail::make_generic_opnode(
		[c](auto&& x){ return std::forward<decltype(x)>(x)*c; },
		std::forward<E>(expr));
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param c	乗数
  \param expr	式
  \return	乗算演算子ノード
*/
template <class E, std::enable_if_t<rank<E>() != 0>* = nullptr> inline auto
operator *(element_t<E> c, E&& expr)
{
    return detail::make_generic_opnode(
		[c](auto&& x){ return c*std::forward<decltype(x)>(x); },
		std::forward<E>(expr));
}

//! 与えられた式の各要素を定数で割る.
/*!
  \param expr	式
  \param c	除数
  \return	除算演算子ノード
*/
template <class E, std::enable_if_t<rank<E>() != 0>* = nullptr> inline auto
operator /(E&& expr, element_t<E> c)
{
    return detail::make_generic_opnode(
		[c](auto&& x){ return std::forward<decltype(x)>(x)/c; },
		std::forward<E>(expr));
}

//! 与えられた式の各要素に定数を掛ける.
/*!
  \param expr	式
  \param c	乗数
  \return	各要素にcが掛けられた結果の式
*/
template <class E> inline std::enable_if_t<rank<E>() != 0, E&>
operator *=(E&& expr, element_t<E> c)
{
    using std::size;
    constexpr size_t	N = size0<E>();
    
    for_each<N>([c](auto&& x){ x *= c; }, size(expr), begin(expr));
    return expr;
}

//! 与えられた式の各要素を定数で割る.
/*!
  \param expr	式
  \param c	除数
  \return	各要素がcで割られた結果の式
*/
template <class E> inline std::enable_if_t<rank<E>() != 0, E&>
operator /=(E&& expr, element_t<E> c)
{
    using std::size;
    constexpr size_t	N = size0<E>();
    
    for_each<N>([c](auto&& x){ x /= c; }, size(expr), begin(expr));
    return expr;
}

//! 与えられた2つの式の各要素の和をとる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	加算演算子ノード
*/
template <class L, class R,
	  std::enable_if_t<rank<L>() != 0 && rank<L>() == rank<R>()>* = nullptr>
inline auto
operator +(L&& l, R&& r)
{
    return detail::make_generic_opnode(
		[](auto&& x, auto&& y){ return std::forward<decltype(x)>(x)
					     + std::forward<decltype(y)>(y); },
		std::forward<L>(l), std::forward<R>(r));
}

//! 与えられた2つの式の各要素の差をとる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	減算演算子ノード
*/
template <class L, class R,
	  std::enable_if_t<rank<L>() != 0 && rank<L>() == rank<R>()>* = nullptr>
inline auto
operator -(L&& l, R&& r)
{
    return detail::make_generic_opnode(
		[](auto&& x, auto&& y){ return std::forward<decltype(x)>(x)
					     - std::forward<decltype(y)>(y); },
		std::forward<L>(l), std::forward<R>(r));
}

//! 与えられた左辺の式の各要素に右辺の式の各要素を加える.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	各要素が加算された左辺の式
*/
template <class L, class R>
inline std::enable_if_t<rank<L>() != 0 && rank<L>() == rank<R>(), L&>
operator +=(L&& l, const R& r)
{
    using std::size;
    constexpr size_t	N = detail::max<size0<L>(), size0<R>()>::value;
    
    for_each<N>([](auto&& x, const auto& y){ x += y; },
		size(l), begin(l), begin(r));
    return l;
}

//! 与えられた左辺の式の各要素から右辺の式の各要素を減じる.
/*!
  \param l	左辺の式
  \param r	右辺の式
  \return	各要素が減じられた左辺の式
*/
template <class L, class R>
inline std::enable_if_t<rank<L>() != 0 && rank<L>() == rank<R>(), L&>
operator -=(L&& l, const R& r)
{
    using std::size;
    constexpr size_t	N = detail::max<size0<L>(), size0<R>()>::value;
    
    for_each<N>([](auto&& x, const auto& y){ x -= y; },
		size(l), begin(l), begin(r));
    return l;
}

/************************************************************************
*  transposition of 2D expressions					*
************************************************************************/
namespace detail
{
  /**********************************************************************
  *  class transpose_opnode<E>						*
  **********************************************************************/
  template <class E>
  class transpose_opnode : public opnode<transpose_opnode<E> >
  {
    public:
		transpose_opnode(E&& expr)
		    :_expr(std::forward<E>(expr))
		{
		}

      constexpr static auto
		size0()
		{
		    return TU::size0<value_t<E> >();
		}
      auto	begin()	const
		{
		    return column_begin(_expr);
		}
      auto	end() const
		{
		    return column_end(_expr);
		}
      auto	size() const
		{
		    return TU::size<1>(_expr);
		}
      decltype(auto)
		operator [](size_t i) const
		{
		    assert(i < size());
		    return *(begin() + i);
		}
      E		transpose()
		{
		    return std::forward<E>(_expr);
		}
      E		transpose() const
		{
		    return _expr;
		}

    private:
      E		_expr;
  };
}
    
//! 2次元配列式が他の2次元配列式を転置したものであるか判定する．
template <class E>
using is_transposed = is_convertible<E, detail::transpose_opnode>;

//! 与えられたスカラもしくは1次元配列式をそのまま返す
/*
  \param expr	スカラまたは1次元配列式
  \return	expr への定数参照
 */
template <class E> inline std::enable_if_t<(rank<E>() < 2), E&&>
transpose(E&& expr)
{
    return std::forward<E>(expr);
}

//! 与えられた2次元配列式の転置を返す
/*
  \param expr	2次元配列式
  \return	expr を転置した2次元配列式
 */ 
template <class E, std::enable_if_t<rank<E>() == 2 &&
				    !is_transposed<E>::value>* = nullptr>
inline auto
transpose(E&& expr)
{
    return detail::transpose_opnode<E>(std::forward<E>(expr));
}

//! 転置された2次元配列式をさらに転置して元に戻す
/*
  \param r	転置された2次元配列式を表すレンジ
  \return	r をさらに転置した2次元配列式
 */ 
template <class E, std::enable_if_t<rank<E>() == 2 &&
				    is_transposed<E>::value>* = nullptr>
inline decltype(auto)
transpose(E&& expr)
{
    return expr.transpose();
}

/************************************************************************
*  generic algorithms for ranges					*
************************************************************************/
//! 与えられた式の各要素の自乗和を求める.
/*!
  \param x	式
  \return	式の各要素の自乗和
*/
template <class E, std::enable_if_t<(rank<E>() != 0)>* = nullptr> inline auto
square(const E& expr)
{
    using	std::size;
    
    return square<size0<E>()>(std::begin(expr), size(expr));
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
square_distance(const L& x, const R& y)
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
distance(const L& x, const R& y)
{
    return std::sqrt(square_distance(x, y));
}

//! 与えられた式のノルムを1に正規化する．
/*!
   \param expr	式 
   \return	正規化された式，すなわち
		\f$
		  \TUvec{e}{}\leftarrow\frac{\TUvec{e}{}}{\TUnorm{\TUvec{e}{}}}
		\f$
*/
template <class E> inline std::enable_if_t<rank<E>() != 0, E&>
normalize(E&& expr)
{
    return expr /= length(expr);
}

template <class E, class T>
inline std::enable_if_t<rank<E>() == 1 && std::is_arithmetic<T>::value,
			element_t<E> >
at(const E& expr, T x)
{
    const auto	x0 = std::floor(x);
    const auto	dx = x - x0;
    const auto	i  = size_t(x0);
    return (dx ? (1 - dx) * expr[i] + dx * expr[i+1] : expr[i]);
}

template <class E, class T>
inline std::enable_if_t<rank<E>() == 2 && std::is_arithmetic<T>::value,
			element_t<E> >
at(const E& expr, T x, T y)
{
    const auto	y0 = std::floor(y);
    const auto	dy = y - y0;
    const auto	i  = size_t(y0);
    const auto	a0 = at(expr[i], x);
    return (dy ? (1 - dy) * a0 + dy * at(expr[i+1], x) : a0);
}

}	// namespace TU
#endif	// !TU_RANGE_H
