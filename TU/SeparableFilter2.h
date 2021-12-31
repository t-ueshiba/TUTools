/*!
  \file		SeparableFilter2.h
  \author	Toshio UESHIBA
  \brief	水平/垂直方向に分離可能な2次元フィルタを実装するための基底クラスの定義
*/
#ifndef TU_SEPARABLEFILTER2_H
#define TU_SEPARABLEFILTER2_H

#include "TU/Array++.h"
#if defined(USE_TBB)
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#endif

namespace TU
{
/************************************************************************
*  class SeparableFilter2<F>						*
************************************************************************/
//! 水平/垂直方向に分離可能な2次元フィルタを実装するための基底クラス
/*!
  \param F	1次元フィルタの型
*/
template <class F>
class SeparableFilter2
{
  public:
    using	filter_type = F;
#if defined(USE_TBB)
  private:
    template <class IN_, class OUT_>
    class ConvolveRows
    {
      public:
	ConvolveRows(F const& filter, IN_ in, OUT_ out, size_t offset)
	    :_filter(filter), _in(in), _out(out), _offset(offset)	{}

	void	operator ()(const tbb::blocked_range<size_t>& r) const
		{
		    using	std::cbegin;
		    using	std::cend;
		    
		    auto	in = _in;
		    std::advance(in, r.begin());
		    auto	ie = _in;
		    std::advance(ie, r.end());
		    for (auto col = r.begin() + _offset; in != ie; ++in, ++col)
			_filter.convolve(cbegin(*in), cend(*in),
					 make_vertical_iterator(_out, col));
		}

      private:
	const F   	_filter;  // cache等の内部状態を持ち得るので参照は不可
	const IN_ 	_in;
	const OUT_	_out;
	const size_t	_offset;
    };

    template <class IN_, class OUT_>
    class ConvolveV
    {
      public:
	ConvolveV(const F& filterV, IN_ ib, IN_ ie, OUT_ out)
	    :_filterV(filterV), _ib(ib), _ie(ie), _out(out)		{}

	void	operator ()(const tbb::blocked_range<size_t>& r) const
		{
		    using	std::cbegin;
		    using	std::begin;
		    
		    _filterV.convolve(make_range_iterator(
					  cbegin(*_ib) + r.begin(),
					  stride(_ib), r.size()),
				      make_range_iterator(
					  cbegin(*_ie) + r.begin(),
					  stride(_ie), r.size()),
				      make_range_iterator(
					  begin(*_out) + r.begin(),
					  stride(_out), r.size()));
		}

      private:
	const F		_filterV;  // cache等の内部状態を持ち得るので参照は不可
	const IN_	_ib;
	const IN_	_ie;
	const OUT_	_out;
    };

    template <class IN_, class OUT_>
    class ConvolveH
    {
      public:
	ConvolveH(const F& filterH, IN_ in, OUT_ out, bool shift)
	    :_filterH(filterH), _in(in), _out(out), _shift(shift)	{}

	void	operator ()(const tbb::blocked_range<size_t>& r) const
		{
		    using	std::cbegin;
		    using	std::cend;
		    using	std::begin;
		    
		    auto	in = _in;
		    std::advance(in, r.begin());
		    auto	ie = _in;
		    std::advance(ie, r.end());
		    auto	out = _out;
		    std::advance(out, r.begin());
		    for (; in != ie; ++in, ++out)
			_filterH.convolve(cbegin(*in), cend(*in),
					  begin(*out), _shift);
		}

      private:
	const F		_filterH;  // cache等の内部状態を持ち得るので参照は不可
	const IN_	_in;
	const OUT_	_out;
	const bool	_shift;
    };
#endif
  public:
    SeparableFilter2()	:_filterH(), _filterV(), _grainSize(1)		{}
    template <class ARGH_, class ARGV_>
    SeparableFilter2(const ARGH_& argH, const ARGV_& argV)
	:_filterH(argH), _filterV(argV), _grainSize(1)			{}
    
    template <class IN_, class OUT_>
    void	convolve(IN_ ib, IN_ ie, OUT_ out,
			 bool shift=false)	const	;
    size_t	grainSize()			const	{ return _grainSize; }
    void	setGrainSize(size_t gs)			{ _grainSize = gs; }
    const F&	filterH()			const	{ return _filterH; }
    const F&	filterV()			const	{ return _filterV; }
    size_t	outSizeH(size_t inSizeH) const
		{
		    return _filterH.outSize(inSizeH);
		}
    size_t	outSizeV(size_t inSizeV) const
		{
		    return _filterV.outSize(inSizeV);
		}
    size_t	offsetH() const
		{
		    return _filterH.offset();
		}
    size_t	offsetV() const
		{
		    return _filterV.offset();
		}
	
  protected:
    F&		filterH()				{ return _filterH; }
    F&		filterV()				{ return _filterV; }

  private:
    F		_filterH;
    F		_filterV;
    size_t	_grainSize;
};
    
//! 与えられた2次元配列とこのフィルタの畳み込みを行う
/*!
  \param ib	入力2次元データ配列の先頭行を指す反復子
  \param ie	入力2次元データ配列の末尾の次の行を指す反復子
  \param out	出力2次元データ配列の先頭行を指す反復子
  \param shift	trueならば，入力データと対応するよう，出力位置を水平/垂直
		方向にそれぞれ offsetH(), offsetV() だけシフトする
*/
template <class F> template <class IN_, class OUT_> void
SeparableFilter2<F>::convolve(IN_ ib, IN_ ie, OUT_ out, bool shift) const
{
    using buf_type	= Array2<typename filter_type::element_type>;
    
    if (ib == ie)
	return;

    if (shift)
	std::advance(out, offsetV());
    
#if defined(CACHE_FRIENDLY)
    buf_type	buf(_filterV.outSize(std::distance(ib, ie)), size(*ib));

#  if defined(USE_TBB)
    using convolveV	= ConvolveV<IN_, typename buf_type::iterator>;
    using convolveH	= ConvolveH<typename buf_type::const_iterator, OUT_>;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, buf.ncol(), _grainSize),
		      convolveV(_filterV, ib, ie, buf.begin()));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, buf.nrow(), _grainSize),
		      convolveH(_filterH, buf.cbegin(), out, shift));
#  else
    _filterV.convolve(ib, ie, buf.begin());
    for (const auto& row : buf)
    {
	using	std::begin;
	
	_filterH.convolve(row.begin(), row.end(), begin(*out), shift);
	++out;
    }
#  endif
#else
    buf_type	buf(_filterH.outSize(size(*ib)), std::distance(ib, ie));

#  if defined(USE_TBB)
    using convolveH	= ConvolveRows<IN_, typename buf_type::iterator>;
    using convolveV	= ConvolveRows<typename buf_type::const_iterator, OUT_>;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, buf.ncol(), _grainSize),
		      convolveH(_filterH, ib, buf.begin(), 0));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, buf.nrow(), _grainSize),
		      convolveV(_filterV,
				buf.cbegin(), out, (shift ? offsetH() : 0)));
#  else
    size_t	col = 0;
    for (; ib != ie; ++ib)
	_filterH.convolve(std::cbegin(*ib), std::cend(*ib),
			  make_vertical_iterator(buf.begin(), col++));
    col = (shift ? offsetH() : 0);
    for (const auto& row : buf)
	_filterV.convolve(row.cbegin(), row.cend(),
			  make_vertical_iterator(out, col++));
#  endif
#endif
}

}
#endif	// !TU_SEPARABLEFILTER2_H
