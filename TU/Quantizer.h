/*!
  \file		Quantizer.h
  \author	Toshio UESHIBA
  \brief	クラス TU::Quantizer, TU::Quantizer2 の定義と実装
*/
#ifndef TU_QUANTIZER_H
#define TU_QUANTIZER_H

#include <vector>
#include "TU/Vector++.h"

namespace TU
{
/************************************************************************
*  class QuantizerBase<T>						*
************************************************************************/
template <class T>
class QuantizerBase
{
  public:
    using value_type	= T;

  private:
    template <class PITER>
    class BinProps
    {
      public:
	using element_type	= float;
	using vector_type	= Vector<element_type, 3>;
	
      public:
	BinProps(PITER begin, PITER end)				;
	BinProps(PITER begin, PITER end, const u_char T::* c, size_t n,
		 const vector_type& sum, const vector_type& sqsum)
	    :_begin(begin), _end(end), _c(c), _n(n),
	     _sum(sum), _var((n*sqsum - square(sum)) / (n*n))		{}

	PITER		begin()			const	{ return _begin; }
	PITER		end()			const	{ return _end; }
	value_type	mean() const
			{
			    return value_type(u_char(_sum[0]/_n),
					      u_char(_sum[1]/_n),
					      u_char(_sum[2]/_n));
			}
	element_type	maxVariance() const
			{
			    return std::max({_var[0], _var[1], _var[2]});
			}
	BinProps	split()				;
	
      private:
	static element_type	square(element_type x)	{ return x*x; }
	static vector_type	square(const vector_type& x)
				{
				    return {x[0]*x[0], x[1]*x[1], x[2]*x[2]};
				}
	
      private:
	PITER			_begin;
	PITER			_end;
	const u_char T::*	_c;
	size_t			_n;
	vector_type		_sum;
	vector_type		_var;
    };
    
  public:
    const T&	operator [](size_t i)	const	{ return _bins[i]; }
    size_t	size()			const	{ return _bins.size(); }

  protected:
    template <class PAIR>
    void	quantize(std::vector<PAIR>& in_out, size_t nbins,
			 std::true_type)				;
    template <class PAIR>
    void	quantize(std::vector<PAIR>& in_out, size_t nbins,
			 std::false_type)				;
	    
  private:
    std::vector<T>	_bins;
};

template <class T> template <class PAIR> void
QuantizerBase<T>::quantize(std::vector<PAIR>& io, size_t nbins, std::true_type)
{
    std::sort(io.begin(), io.end(),
	      [](const PAIR& x, const PAIR& y){return *x.first < *y.first;});
    
  // 入力データをnbin個のbinに量子化した場合のbin幅の最大値を求める
    T	binWidthMax = 2*(*io.back().first - *io.front().first)/nbins;
    for (T th = 1e-5, l = 0; binWidthMax - l > th; )
    {
	T	binMin = *io.front().first;	// bin中の最小値
	T	m = 0.5*(l + binWidthMax);
	size_t	n = 0;
	bool	ok = true;

	for (const auto& x : io)
	    if (*x.first > binMin + m)	// bin中の最大値を越えるなら...
	    {
		if (++n == nbins)	// 新たにbinを作ると
		{			// bin数の上限を越えるなら...
		    ok = false;		// bin幅が不十分
		    break;
		}
		binMin = *x.first;	// 新たなbin中の最小値
	    }

	if (ok)			// 現在のmの値でbin幅が十分ならば...
	    binWidthMax = m;	// mを新たな上限値に
	else			// そうでなければ...
	    l = m;		// mを新たな下限値に
    }

  // 入力データの各元をbinの代表元にマップ
    _bins.clear();
    auto	binBase = io.cbegin();
    size_t	idx = 0;
    for (auto x = io.cbegin(); x != io.cend(); ++x)
	if (*x->first > *binBase->first + binWidthMax)
	{
	    _bins.push_back(*binBase[(x - binBase) >> 1].first);  // 代表元
			    
	    for (; binBase != x; ++binBase)
		*binBase->second = idx;
	    ++idx;
	}
    _bins.push_back(*binBase[(io.cend() - binBase) >> 1].first);
    for (; binBase != io.cend(); ++binBase)
	*binBase->second = idx;
}

template <class T> template <class PAIR> void
QuantizerBase<T>::quantize(std::vector<PAIR>& io, size_t nbins, std::false_type)
{
    using piterator	= typename std::vector<PAIR>::iterator;
    using props_type	= BinProps<piterator>;

    const auto	comp = [](const props_type& x, const props_type& y)
			{ return x.maxVariance() < y.maxVariance(); };
	
    std::vector<props_type>	props;
    props.push_back(props_type(io.begin(), io.end()));
    std::make_heap(props.begin(), props.end(), comp);

    for (size_t n = 1; n < nbins; ++n)
    {
	std::pop_heap(props.begin(), props.end(), comp);
	auto	prop1 = props.back();
	props.pop_back();

	auto	prop2 = prop1.split();

	props.push_back(prop1);
	std::push_heap(props.begin(), props.end(), comp);
	props.push_back(prop2);
	std::push_heap(props.begin(), props.end(), comp);
    }

  // 入力データの各元をbinの代表元にマップ
    for (size_t idx = 0; idx < props.size(); ++idx)
    {
	_bins.push_back(props[idx].mean());
	for (auto& x: props[idx])
	    *x.second = idx;
    }
}

/*
 *  QuantizerBase<T>::BinProps<PITER>
 */
template <class T> template <class PITER>
QuantizerBase<T>::BinProps<PITER>::BinProps(PITER begin, PITER end)
    :_begin(begin), _end(end), _c(nullptr), _n(0)
{
    for (auto iter = _begin; iter != _end; ++iter)
    {
	vector_type	x;
	x[0] = iter->first->r;
	x[1] = iter->first->g;
	x[2] = iter->first->b;
	
	++_n;
	_sum += x;
	_var += square(x);
    }

    if (_n)
	_var = (_n*_var - square(_sum)) / (_n*_n);
}
    
template <class T>
template <class PITER> typename QuantizerBase<T>::template BinProps<PITER>
QuantizerBase<T>::BinProps<PITER>::split()
{
    using pair_type	= iterator_value<PITER>;

  // 要素への反復子の配列を作り，昇順にソート
    const u_char T::*	c = (_var[0] > _var[1] ?
			     _var[0] > _var[2] ? &T::r : &T::b :
			     _var[1] > _var[2] ? &T::g : &T::b);
    if (c != _c)
	std::sort(_begin, _end,
		  [=](const pair_type& x, const pair_type& y)
		  { return x.first->*c < y.first->*c; });
    
  // 大津の判別基準により最適なしきい値を決定．
    const size_t	i = (_var[0] > _var[1] ?
			     _var[0] > _var[2] ? 0 : 2 :
			     _var[1] > _var[2] ? 1 : 2);
    auto		border = _begin;
    auto		head   = _begin;
    const auto		mean = _sum[i]/_n;
    size_t		n = 0;			// しきい値以下の要素数
    vector_type		sum;			// しきい値以下の累積値
    vector_type		sqsum;			// しきい値以下の自乗累積値
    size_t		n1 = 0;
    vector_type		sum1;
    vector_type		sqsum1;
    element_type	interVarianceMax = 0;	// クラス間分散の最大値
    for (auto iter = _begin; iter != _end; ++iter)
    {
	vector_type	x;
	x[0] = iter->first->r;
	x[1] = iter->first->g;
	x[2] = iter->first->b;
	
	if (iter->first->*c != head->first->*c)
	{
	    const auto	interVariance = square(sum[i] - n*mean) / (n*(_n - n));
	    if (interVariance > interVarianceMax)
	    {
		interVarianceMax = interVariance;
		border		 = iter;
		n1		 = n;
		sum1		 = sum;
		sqsum1		 = sqsum;
	    }
	    head = iter;
	}
	
	++n;
	sum   += x;
	sqsum += square(x);
    }

    const BinProps	props(border, _end,
			      c, n - n1, sum - sum1, sqsum - sqsum1);
    
    _end = border;
    _c	 = c;
    _n	 = n1;
    _sum = sum1;
    _var = (n1*sqsum1 - square(sum1)) / (n1*n1);
    
    return props;
}
    
/*
 *  QuantizerBase<u_char> : specialized
 */
template <>
class QuantizerBase<u_char>
{
  public:
    size_t	operator [](size_t i)	const	{ return i; }
    size_t	size()			const	{ return 256; }
};
    
/************************************************************************
*  class Quantizer<T>							*
************************************************************************/
template <class T>
class Quantizer : public QuantizerBase<T>
{
  private:
    using super	= QuantizerBase<T>;
    
  public:
    template <class ITER>
    std::enable_if_t<std::is_same<iterator_value<ITER>, u_char>::value,
		     range<ITER> >
	operator ()(ITER ib, ITER ie, size_t)
	{
	    return range<ITER>(ib, std::distance(ib, ie));
	}
    template <class ITER>
    std::enable_if_t<!std::is_same<iterator_value<ITER>, u_char>::value,
		     const Array<size_t>&>
	operator ()(ITER ib, ITER ie, size_t nbins)			;
    friend std::ostream&
	operator <<(std::ostream& out, const Quantizer<T>& quantizer)
	{
	    return out << quantizer._indices;
	}
    
  private:
    Array<size_t>	_indices;
};

template <class T> template <class ITER>
std::enable_if_t<!std::is_same<iterator_value<ITER>, u_char>::value,
		 const Array<size_t>&>
Quantizer<T>::operator ()(ITER ib, ITER ie, size_t nbins)
{
    using pair_type = std::pair<ITER, Array<size_t>::iterator>;

    _indices.resize(std::distance(ib, ie));
    
    std::vector<pair_type>	io;
    for (auto idx = _indices.begin(); ib != ie; ++ib, ++idx)
	io.push_back(pair_type(ib, idx));

    super::quantize(io, nbins, std::is_arithmetic<T>());

    return _indices;
}
    
/************************************************************************
*  class Quantizer2<T>							*
************************************************************************/
template <class T>
class Quantizer2 : public QuantizerBase<T>
{
  private:
    using super	= QuantizerBase<T>;
    
  public:
    template <class ROW>
    std::enable_if_t<std::is_same<value_t<iterator_value<ROW> >,
				  u_char>::value, range<ROW> >
	operator ()(ROW ib, ROW ie, size_t nbins)
	{
	    return range<ROW>(ib, std::distance(ib, ie));
	}
    template <class ROW>
    std::enable_if_t<!std::is_same<value_t<iterator_value<ROW> >,
				   u_char>::value, const Array2<size_t>&>
	operator ()(ROW ib, ROW ie, size_t nbins)			;
    friend std::ostream&
	operator <<(std::ostream& out, const Quantizer2<T>& quantizer)
	{
	    return out << quantizer._indices;
	}

  private:
    Array2<size_t>	_indices;
};

template <class T> template <class ROW>
std::enable_if_t<!std::is_same<value_t<iterator_value<ROW> >, u_char>::value,
		 const Array2<size_t>&>
Quantizer2<T>::operator ()(ROW ib, ROW ie, size_t nbins)
{
    using pair_type = std::pair<iterator_t<iterator_reference<ROW> >,
				Array<size_t>::iterator>;
    
    _indices.resize(std::distance(ib, ie),
		    (ib == ie ? 0 : std::distance(ib->begin(), ib->end())));
    
    std::vector<pair_type>	io;
    for (auto row = _indices.begin(); ib != ie; ++ib, ++row)
    {
	auto	idx = row->begin();
	for (auto col = ib->begin(); col != ib->end(); ++col, ++idx)
	    io.push_back(pair_type(col, idx));
    }

    super::quantize(io, nbins, std::is_arithmetic<T>());

    return _indices;
}

}
#endif	// !TU_QUANTIZER_H
