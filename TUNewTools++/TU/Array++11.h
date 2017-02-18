/*
 *  $Id$
 */
/*!
  \file		Array++.h
  \brief	多次元配列クラスの定義と実装
*/
#ifndef __TU_ARRAY_H
#define __TU_ARRAY_H

#include <array>
#include <iomanip>		// for std::ws
#include <memory>		// for std::allocator<T>, std::unique_ptr<T>
#include "TU/range++11.h"
#include "TU/utility.h"		// for std::index_sequence<Ints...>
#include "TU/algorithm.h"	// for TU::lcm()

namespace TU
{
/************************************************************************
*  class external_allocator<T>						*
************************************************************************/
template <class T>
class external_allocator
{
  public:
    using value_type		= T;
    using pointer		= T*;
    using const_pointer		= const T*;
    using reference		= T&;
    using const_reference	= const T&;
    using size_type		= size_t;
    using difference_type	= ptrdiff_t;
    
    template <class T_>
    struct rebind	{ using other = external_allocator<T_>; };

  public:
			external_allocator(pointer p, size_type size)
			    :_p(p), _size(size)				{}

    pointer		allocate(size_type n,
				 typename std::allocator<void>
					     ::const_pointer=nullptr) const
			{
			    if (n > _size)
				throw std::runtime_error("TU::external_allocator<T>::allocate(): too large memory requested!");
			    
			    return _p;
			}
    static void		deallocate(pointer, size_type)			{}
    static void		construct(pointer p, const_reference val)	{}
    static void		destroy(pointer p)				{}
    constexpr
    size_type		max_size()		const	{ return _size; }
    static pointer	address(reference r)		{ return &r; }
    static const_pointer
			address(const_reference r)	{ return &r; }

  private:
    const pointer	_p;
    const size_type	_size;
};
    
/************************************************************************
*  class BufTraits<T>							*
************************************************************************/
template <class T, class ALLOC>
struct BufTraits
{
    using allocator_type = ALLOC;
    using pointer	 = typename allocator_type::pointer;
    using const_pointer	 = typename allocator_type::const_pointer;
    using iterator	 = pointer;
    using const_iterator = const_pointer;
    
  protected:
    static pointer	null()
			{
			    return nullptr;
			}
    
    template <class IN_, class OUT_>
    static OUT_		copy(IN_ ib, IN_ ie, OUT_ out)
			{
			    return std::copy(ib, ie, out);
			}

    template <class T_>
    static void		fill(iterator ib, iterator ie, const T_& c)
			{
			    std::fill(ib, ie, c);
			}

    static void		init(iterator ib, iterator ie)
			{
			    std::fill(ib, ie, 0);
			}
};

/************************************************************************
*  class Buf<T, ALLOC, SIZE, SIZES...>					*
************************************************************************/
//! 固定長多次元バッファクラス
/*!
  単独で使用することはなく，#array の内部バッファクラスとして使う．
  \param T	要素の型
  \param ALLOC	アロケータの型
  \param SIZE	最初の軸の要素数
  \param SIZES	2番目以降の各軸の要素数
*/
template <class T, class ALLOC, size_t SIZE, size_t... SIZES>
class Buf : public BufTraits<T, ALLOC>
{
  private:
  // このバッファの総容量をコンパイル時に計算
    template <size_t SIZE_, size_t... SIZES_>
    struct prod
    {
	constexpr static size_t	value = SIZE_ * prod<SIZES_...>::value;
    };
    template <size_t SIZE_>
    struct prod<SIZE_>
    {
	constexpr static size_t value = SIZE_;
    };

    constexpr static size_t	Capacity = prod<SIZE, SIZES...>::value;

  // このバッファの各軸のサイズをコンパイル時に計算
    template <size_t I_, size_t SIZE_, size_t... SIZES_>
    struct nth
    {
	constexpr static size_t	value = nth<I_-1, SIZES_...>::value;
    };
    template <size_t SIZE_, size_t... SIZES_>
    struct nth<0, SIZE_, SIZES_...>
    {
	constexpr static size_t	value = SIZE_;
    };

    template <size_t I_>
    using siz			= nth<I_, SIZE, SIZES...>;
    template <size_t I_>
    using axis			= std::integral_constant<size_t, I_>;
    using super			= BufTraits<T, ALLOC>;

    using buf_type		= std::array<T, Capacity>;
    using base_iterator		= typename buf_type::iterator;
    using const_base_iterator	= typename buf_type::const_iterator;

  public:
    constexpr static size_t	D = 1 + sizeof...(SIZES);

    using sizes_type		= std::array<size_t, D>;
    using value_type		= T;
    using allocator_type	= void;
    using typename super::pointer;
    using typename super::const_pointer;

  private:
    void	init(std::true_type)		{ _a.fill(0); }
    void	init(std::false_type)		{}

    template <size_t I_>
    static bool	check_sizes(const sizes_type& sizes, axis<I_>)
		{
		    return (sizes[I_-1] != size<I_-1>() ? false :
			    check_sizes(sizes, axis<I_-1>()));
		}
    static bool	check_sizes(const sizes_type& sizes, axis<0>)
		{
		    return true;
		}

    template <class ITER_>
    static ITER_
		make_iterator(ITER_ iter)
		{
		    return iter;
		}
    template <size_t SIZE_, size_t... SIZES_, class ITER_>
    static auto	make_iterator(ITER_ iter)
		    -> decltype(make_range_iterator<SIZE_, SIZE_>(
				    Buf::make_iterator<SIZES_...>(iter)))
		{
		    return make_range_iterator<SIZE_, SIZE_>(
			       make_iterator<SIZES_...>(iter));
		}

    using iterator	 = decltype(Buf::make_iterator<SIZES...>(
					std::declval<base_iterator>()));
    using const_iterator = decltype(Buf::make_iterator<SIZES...>(
					std::declval<const_base_iterator>()));
    
  public:
  // 標準コンストラクタ/代入演算子およびデストラクタ
		Buf()
		{
		    init(typename std::is_arithmetic<value_type>::type());
		}
		Buf(const Buf&)			= default;
    Buf&	operator =(const Buf&)		= default;
		Buf(Buf&&)			= default;
    Buf&	operator =(Buf&&)		= default;
    
  // 各軸のサイズと最終軸のストライドを指定したコンストラクタとリサイズ関数
    explicit	Buf(const sizes_type& sizes, size_t=0)
		{
		    if (!check_sizes(sizes, axis<D>()))
			throw std::logic_error("Buf<T, ALLOC, SIZE, SIZES...>::Buf(): mismatched size!");
		}
    void	resize(const sizes_type& sizes, size_t=0)
		{
		    if (!check_sizes(sizes, axis<D>()))
			throw std::logic_error("Buf<T, ALLOC, SIZE, SIZES...>::resize(): mismatched size!");
		}

    template <size_t I_=0>
    constexpr static size_t	size()		{ return siz<I_>::value; }
    template <size_t I_=D-1>
    constexpr static size_t	stride()	{ return siz<I_>::value; }
    constexpr static size_t	nrow()		{ return siz<0>::value; }
    constexpr static size_t	ncol()		{ return siz<1>::value; }

    pointer	data()		{ return _a.data(); }
    const_pointer
		data()	const	{ return _a.data(); }
    iterator	begin()		{ return make_iterator<SIZES...>(_a.begin()); }
    const_iterator
		begin()	const	{ return make_iterator<SIZES...>(_a.begin()); }
    iterator	end()		{ return make_iterator<SIZES...>(_a.end()); }
    const_iterator
		end()	const	{ return make_iterator<SIZES...>(_a.end()); }

    void	fill(const T& c)		{ _a.fill(c); }
    std::istream&
		get(std::istream& in)
		{
		    for (auto& val : _a)
			in >> val;
		    return in;
		}
    
  private:
    alignas(sizeof(T)) buf_type	_a;
};

//! 可変長多次元バッファクラス
/*!
  単独で使用することはなく，#array の内部バッファクラスとして使う．
  \param T	要素の型
  \param ALLOC	アロケータの型
  \param SIZES	ダミー(各軸の要素数は動的に決定される)
*/
template <class T, class ALLOC, size_t... SIZES>
class Buf<T, ALLOC, 0, SIZES...> : public BufTraits<T, ALLOC>
{
  private:
    using super			= BufTraits<T, ALLOC>;
    using base_iterator		= typename super::iterator;
    using const_base_iterator	= typename super::const_iterator;
    template <size_t I_>
    using axis			= std::integral_constant<size_t, I_>;

    template <class ITER_>
    static ITER_
		make_iterator_dummy(ITER_)			;
    template <size_t SIZE_, size_t... SIZES_, class ITER_>
    static auto	make_iterator_dummy(ITER_ iter)
		    -> decltype(make_range_iterator(
				    Buf::make_iterator_dummy<SIZES_...>(iter),
				    0, 0));

  public:
    constexpr static size_t	D = 1 + sizeof...(SIZES);

    using sizes_type		= std::array<size_t, D>;
    using value_type		= T;
    using typename super::allocator_type;
    using typename super::pointer;
    using typename super::const_pointer;
    using iterator
		= decltype(Buf::make_iterator_dummy<SIZES...>(
				       std::declval<base_iterator>()));
    using const_iterator
		= decltype(Buf::make_iterator_dummy<SIZES...>(
			       std::declval<const_base_iterator>()));
    
  public:
  // 標準コンストラクタ/代入演算子およびデストラクタ
		Buf()
		    :_stride(0), _capacity(0), _p(super::null())
		{
		    _sizes.fill(0);
		}
		Buf(const Buf& b)
		    :_sizes(b._sizes), _stride(b._stride),
		     _capacity(b._capacity), _p(alloc(_capacity))
		{
		    super::copy(b.begin(), b.end(), begin());
		}
    Buf&	operator =(const Buf& b)
		{
		    if (this != &b)
		    {
			resize(b._sizes, b._stride);
			super::copy(b.begin(), b.end(), begin());
		    }
		    return *this;
		}
		Buf(Buf&& b)
		    :_sizes(b._sizes), _stride(b._stride),
		     _capacity(b._capacity), _p(b._p)
		{
		  // b の 破壊時に this->_p がdeleteされることを防ぐ．
		    b._p = super::null();
		}
    Buf&	operator =(Buf&& b)
		{
		    _sizes    = b._sizes;
		    _stride   = b._stride;
		    _capacity = b._capacity;
		    _p        = b._p;
		    
		  // b の 破壊時に this->_p がdeleteされることを防ぐ．
		    b._p = super::null();
		}
		~Buf()
		{
		    free(_p, _capacity);
		}

  // 各軸のサイズと最終軸のストライドを指定したコンストラクタとリサイズ関数
    explicit	Buf(const sizes_type& sizes, size_t stride=0)
		    :_sizes(sizes),
		     _stride(stride ? stride : _sizes[D-1]),
		     _capacity(capacity(axis<0>())),
		     _p(alloc(_capacity))
		{
		}
    void	resize(const sizes_type& sizes, size_t stride=0)
		{
		    free(_p, _capacity);
		    _sizes    = sizes;
		    _stride   = (stride ? stride : _sizes[D-1]);
		    _capacity = capacity(axis<0>());
		    _p	      = alloc(_capacity);
		}

		Buf(pointer p, const sizes_type& sizes, size_t stride=0)
		    :_sizes(sizes),
		     _stride(stride ? stride : _sizes[D-1]),
		     _capacity(capacity(axis<0>())),
		     _allocator(p, _capacity),
		     _p(alloc(_capacity))
		{
		}

    const sizes_type&
		sizes()		const	{ return _sizes; }
    template <size_t I_=0>
    size_t	size()		const	{ return _sizes[I_]; }
    template <size_t I_=D-1>
    size_t	stride()	const	{ return stride_impl(axis<I_>()); }
    size_t	nrow()		const	{ return _sizes[0]; }
    size_t	ncol()		const	{ return _sizes[1]; }
    pointer	data()			{ return _p; }
    const_pointer
		data()		const	{ return _p; }
    iterator	begin()
		{
		    return make_iterator<SIZES...>(base_iterator(_p));
		}
    const_iterator
		begin() const
		{
		    return make_iterator<SIZES...>(const_base_iterator(_p));
		}
    iterator	end()
		{
		    return make_iterator<SIZES...>(base_iterator(
						       _p + _capacity));
		}
    const_iterator
		end() const
		{
		    return make_iterator<SIZES...>(const_base_iterator(
						       _p + _capacity));
		}
    void	fill(const T& c)
		{
		    super::fill(_p, _p + _capacity, c);
		}
    std::istream&
		get(std::istream& in)
		{
		    sizes_type	nvalues, sizes;
		    nvalues.fill(0);
		    sizes.fill(0);

		    get(in >> std::ws, nvalues, sizes);

		    return in;
		}
    
  private:
    template <size_t I_>
    size_t	stride_impl(axis<I_>)		const	{ return _sizes[I_]; }
    size_t	stride_impl(axis<D-1>)		const	{ return _stride; }

    size_t	capacity(axis<D-1>) const
		{
		    return _stride;
		}
    template <size_t I_>
    size_t	capacity(axis<I_>) const
		{
		    return _sizes[I_] * capacity(axis<I_+1>());
		}

    pointer	alloc(size_t siz)
		{
		    const auto	p = _allocator.allocate(siz);
		    for (pointer q = p, qe = q + siz; q != qe; ++q)
			_allocator.construct(q, value_type());
		    return p;
		}
    void	free(pointer p, size_t siz)
		{
		    if (p != super::null())
		    {
			for (pointer q = p, qe = q + siz; q != qe; ++q)
			    _allocator.destroy(q);
			_allocator.deallocate(p, siz);
		    }
		}

    template <class ITER_>
    static ITER_
		make_iterator(ITER_ iter)
		{
		    return iter;
		}
    template <size_t SIZE_, size_t... SIZES_, class ITER_>
    auto	make_iterator(ITER_ iter) const
		    -> decltype(make_range_iterator(
				    make_iterator<SIZES_...>(iter),
				    size<D-1-sizeof...(SIZES_)>(),
				    stride<D-1-sizeof...(SIZES_)>()))
		{
		    constexpr size_t	I = D - 1 - sizeof...(SIZES_);
		    
		    return make_range_iterator(
			       make_iterator<SIZES_...>(iter),
			       size<I>(), stride<I>());
		}

    base_iterator
		get(std::istream& in, sizes_type& nvalues, sizes_type& sizes)
		{
		    constexpr size_t	BufSiz = (sizeof(value_type) < 2048 ?
						  2048/sizeof(value_type) : 1);
		    std::unique_ptr<value_type[]>
					tmp(new value_type[BufSiz]);
		    base_iterator	iter;
		    size_t		n = 0;
		    
		    for (size_t d = D - 1; n < BufSiz; )
		    {
			char	c;
			
			while (in.get(c))
			    if (!isspace(c) || c == '\n')
				break;
			
			if (in && c != '\n')	// 現在軸の末尾でなければ...
			{
			    in.putback(c);	// 1文字読み戻して
			    in >> tmp[n++];	// 改めて要素をバッファに読み込む

			    d = D - 1;		// 最下位軸に戻して
			    ++nvalues[d];	// 要素数を1だけ増やす
			}
			else			// 現在軸の末尾に到達したなら...
			{
			    if (nvalues[d] > sizes[d])
				sizes[d] = nvalues[d];	// 現在軸の要素数を記録

			    if (d == 0)		// 最上位軸の末尾ならば...
			    {
				resize(sizes);	// 領域を確保して
				iter = base_iterator(_p + _capacity);
				break;		// その末端をiterにセットして返す
			    }
		
			    nvalues[d] = 0;	// 現在軸を先頭に戻し
			    ++nvalues[--d];	// 直上軸に移動して1つ進める
			}
		    }

		    if (n == BufSiz)		// バッファが一杯ならば...
			iter = get(in, nvalues, sizes);	// 再帰してさらに読み込む

		    while (n--)
			*(--iter) = std::move(tmp[n]);	// バッファの内容を移す

		    return iter;		// 読み込まれた先頭位置を返す
		}

  private:
    sizes_type		_sizes;		//!< 各軸の要素数
    size_t		_stride;	//!< 最終軸のストライド
    size_t		_capacity;	//!< バッファ中に収めらる総要素数
    allocator_type	_allocator;	//!< 要素を確保するアロケータ
    pointer		_p;		//!< 先頭要素へのポインタ
};
    
/************************************************************************
*  class array<T, ALLOC, SIZE, SIZES...>				*
************************************************************************/
//! 多次元配列を表すクラス
/*!
  \param T	要素の型
  \param ALLOC	アロケータの型
  \param SIZE	最初の軸の要素数
  \param SIZES	2番目以降の各軸の要素数
*/
template <class T, class ALLOC, size_t SIZE, size_t... SIZES>
class array : public Buf<T, ALLOC, SIZE, SIZES...>
{
  private:
    using super	= Buf<T, ALLOC, SIZE, SIZES...>;
    
  public:
    using super::D;
    using typename super::sizes_type;
    using typename super::pointer;
    using typename super::const_pointer;
    using typename super::iterator;
    using typename super::const_iterator;
    using element_type		 = T;
    using reverse_iterator	 = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using value_type		 = typename std::iterator_traits<iterator>
					       ::value_type;
    using const_value_type	 = typename std::iterator_traits<const_iterator>
					       ::value_type;
    using reference		 = typename std::iterator_traits<iterator>
					       ::reference;
    using const_reference	 = typename std::iterator_traits<const_iterator>
					       ::reference;
    
  public:
		array()				= default;
		array(const array&)		= default;
    array&	operator =(const array&)	= default;
		array(array&&)			= default;
    array&	operator =(array&&)		= default;
    
    template <class... SIZES_,
	      typename std::enable_if<sizeof...(SIZES_) == D>::type* = nullptr>
    explicit	array(SIZES_... sizes)
		    :super({to_size(sizes)...})
		{
		}
    template <class... SIZES_>
    typename std::enable_if<sizeof...(SIZES_) == D>::type
		resize(SIZES_... sizes)
		{
		    super::resize({to_size(sizes)...});
		}
    
    template <class... SIZES_,
	      typename std::enable_if<sizeof...(SIZES_) == D>::type* = nullptr>
    explicit	array(size_t unit, SIZES_... sizes)
		    :super({to_size(sizes)...}, to_stride(unit, sizes...))
		{
		}
    template <class... SIZES_>
    typename std::enable_if<sizeof...(SIZES_) == D>::type
		resize(size_t unit, SIZES_... sizes)
		{
		    super::resize({to_size(sizes)...},
				  to_stride(unit, sizes...));
		}

    template <class E_,
	      typename std::enable_if<is_range<E_>::value>::type* = nullptr>
		array(const E_& expr)
		    :super(sizes(expr, std::make_index_sequence<D>()))
		{
		    std::copy(std::begin(expr), std::end(expr), begin());
		}
    template <class E_>
    typename std::enable_if<is_range<E_>::value, array&>::type
		operator =(const E_& expr)
		{
		    super::resize(sizes(expr, std::make_index_sequence<D>()));
		    std::copy(std::begin(expr), std::end(expr), begin());

		    return *this;
		}

		array(std::initializer_list<const_value_type> args)
		    :super(sizes(args))
		{
		    std::copy(args.begin(), args.end(), begin());
		}
    array&	operator =(std::initializer_list<const_value_type> args)
		{
		    super::resize(sizes(args));
		    std::copy(args.begin(), args.end(), begin());

		    return *this;
		}

    template <class... SIZES_,
	      typename std::enable_if<sizeof...(SIZES_) == D>::type* = nullptr>
    explicit	array(pointer p, SIZES_... sizes)
		    :super(p, {to_size(sizes)...})
		{
		}

    template <class... SIZES_,
	      typename std::enable_if<sizeof...(SIZES_) == D>::type* = nullptr>
    explicit	array(pointer p, size_t unit, SIZES_... sizes)
		    :super(p, {to_size(sizes)...}, to_stride(unit, sizes...))
		{
		}

    template <class ALLOC_>
    void	write(array<T, ALLOC_, SIZE, SIZES...>& a) const
		{
		    a.resize(sizes(), a.stride());
		    super::copy(begin(), end(), a.begin());
		}

    constexpr static size_t
		dimension()		{ return D; }
    using	super::size;
    using	super::stride;
    using	super::nrow;
    using	super::ncol;
    using	super::data;
    using	super::begin;
    using	super::end;

    range<iterator>
		operator ()()		{ return {begin(), end()}; }
    range<const_iterator>
		operator ()()	const	{ return {begin(), end()}; }
	    
    const_iterator
		cbegin()  const	{ return begin(); }
    const_iterator
		cend()	  const	{ return end(); }
    reverse_iterator
		rbegin()	{ return reverse_iterator(end()); }
    const_reverse_iterator
		rbegin()  const	{ return const_reverse_iterator(end()); }
    const_reverse_iterator
		crbegin() const	{ return rbegin(); }
    reverse_iterator
		rend()		{ return reverse_iterator(begin()); }
    const_reverse_iterator
		rend()	  const	{ return const_reverse_iterator(begin()); }
    const_reverse_iterator
		crend()	  const	{ return rend(); }
    reference	operator [](size_t i)
		{
		    assert(i < size());
		    return *(begin() + i);
		}
    const_reference
		operator [](size_t i) const
		{
		    assert(i < size());
		    return *(begin() + i);
		}
    void	fill(const element_type& val)
		{
		    super::fill(val);
		}
    std::istream&
		restore(std::istream& in)
		{
		    restore(in, begin(), end());
		    return in;
		}
    std::ostream&
		save(std::ostream& out) const
		{
		    save(out, begin(), end());
		    return out;
		}
    
  private:
    using	sizes_iterator = typename sizes_type::iterator;
    
    template <class T_>
    static typename std::enable_if<std::is_integral<T_>::value, size_t>::type
		to_size(const T_& arg)
		{
		    return size_t(arg);
		}

    template <class E_, size_t... I_>
    static sizes_type
		sizes(const E_& expr, std::index_sequence<I_...>)
		{
		    return {TU::size<I_>(expr)...};
		}

    template <class T_>
    static typename std::enable_if<!is_range<T_>::value>::type
		set_sizes(sizes_iterator iter, sizes_iterator end, const T_& val)
		{
		    throw std::runtime_error("array<BUF>::set_sizes(): too shallow initializer list!");
		}
    template <class T_>
    static typename std::enable_if<is_range<T_>::value>::type
		set_sizes(sizes_iterator iter, sizes_iterator end, const T_& r)
		{
		    *iter = r.size();
		    if (++iter != end)
			set_sizes(iter, end, *r.begin());
		}
    static sizes_type
		sizes(std::initializer_list<const_value_type> args)
		{
		    sizes_type	sizs;
		    set_sizes(sizs.begin(), sizs.end(), args);

		    return sizs;
		}
    
    static size_t
		to_stride(size_t unit, size_t size)
		{
		    constexpr auto	elmsiz = sizeof(element_type);

		    if (unit == 0)
			unit = 1;
		    const auto	n = lcm(elmsiz, unit)/elmsiz;

		    return n*((size + n - 1)/n);
		}
    template <class... SIZES_>
    static size_t
		to_stride(size_t unit, size_t size, SIZES_... sizes)
		{
		    return to_stride(unit, sizes...);
		}

    static void	restore(std::istream& in, pointer begin, pointer end)
		{
		    in.read(reinterpret_cast<char*>(begin),
			    sizeof(element_type) * std::distance(begin, end));
		}
    template <class ITER_>
    static void	restore(std::istream& in, ITER_ begin, ITER_ end)
		{
		    for (; begin != end; ++begin)
			restore(in, begin->begin(), begin->end());
		}

    static void	save(std::ostream& out, const_pointer begin, const_pointer end)
		{
		    out.write(reinterpret_cast<const char*>(begin),
			      sizeof(element_type) * std::distance(begin, end));
		}
    template <class ITER_>
    static void	save(std::ostream& out, ITER_ begin, ITER_ end)
		{
		    for (; begin != end; ++begin)
			save(out, begin->begin(), begin->end());
		}
};

template <size_t I, class T, class ALLOC, size_t SIZE, size_t... SIZES>
inline size_t
size(const array<T, ALLOC, SIZE, SIZES...>& a)
{
    return a.template size<I>();
}

template <size_t I, class T, class ALLOC, size_t SIZE, size_t... SIZES>
inline size_t
stride(const array<T, ALLOC, SIZE, SIZES...>& a)
{
    return a.template stride<I>();
}

template <class T, class ALLOC, size_t SIZE, size_t... SIZES> std::ostream&
operator <<(std::ostream& out, const array<T, ALLOC, SIZE, SIZES...>& a)
{
    for (const auto& val : a)
	out << ' ' << val;
    return out << std::endl;
}
    
template <class T, class ALLOC, size_t SIZE, size_t... SIZES>
inline std::istream&
operator >>(std::istream& in, array<T, ALLOC, SIZE, SIZES...>& a)
{
    return a.get(in);
}

/************************************************************************
*  type definitions for convenience					*
************************************************************************/
template <class T, size_t N=0, class ALLOC=std::allocator<T> >
using Array = array<T, ALLOC, N>;

template <class T, size_t R=0, size_t C=0, class ALLOC=std::allocator<T> >
using Array2 = array<T, ALLOC, R, C>;

template <class T,
	  size_t Z=0, size_t Y=0, size_t X=0, class ALLOC=std::allocator<T> >
using Array3 = array<T, ALLOC, Z, Y, X>;

}	// namespace TU
#endif	// !__TU_ARRAY_H