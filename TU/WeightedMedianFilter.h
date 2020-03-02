/*!
  \file		WeightedMedianFilter.h
  \author	Toshio UESHIBA
  \brief	クラス TU::WeightedMedianFilter の定義と実装
*/
#ifndef TU_WEIGHTEDMEDIANFILTER_H
#define TU_WEIGHTEDMEDIANFILTER_H

#include <boost/intrusive/list.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "TU/Quantizer.h"
#include "TU/algorithm.h"	// diff(const T&, const T&)
#if defined(USE_TBB)
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#endif
#include "TU/Profiler.h"

#if defined(PROFILE) && !defined(USE_TBB)
#  define ENABLE_PROFILER
#else
#  define ENABLE_PROFILER	void
#endif

namespace TU
{
/************************************************************************
*  class ExpDiff<S, T>							*
************************************************************************/
//! exp(-abs(x - y)/sigma) を返す関数オブジェクト
/*!
  引数が算術型でない場合はRGBカラーとみなし，abs(x - y)の代わりに3次元ベクトル
	(x.r - y.r, x.g - y.g, x.b - y.b)
  の長さを用いる．
  \param S	引数の型
  \param T	結果の型
*/
template <class S, class T>
class ExpDiff
{
  public:
    using argument_type	= S;	//!< 引数の型
    using result_type	= T;	//!< 結果の型

  public:
		ExpDiff(T sigma=1) :_sigma(sigma)	{}

    void	setSigma(result_type sigma)		{ _sigma = sigma; }
    T		sigma()					{ return _sigma; }
    T		operator ()(S x, S y) const
		{
		    return f(x, y, typename std::is_arithmetic<S>::type());
		}
    
  private:
    T		f(S x, S y, std::true_type) const
		{
		    return std::exp(-diff(x, y) / _sigma);
		}
    T		f(S x, S y, std::false_type) const
		{
		    Vector<T, 3>	a{T(x.r) - T(y.r),
					  T(x.g) - T(y.g), T(x.b) - T(y.b)};
		    return std::exp(-length(a) / _sigma);
		}

  private:
    T		_sigma;
};
    
namespace detail
{
/************************************************************************
*  class WeightedMedianFilterBase<W>					*
************************************************************************/
//! 重み付けメディアンフィルタの実装のベースとなるクラス
/*!
  \param W	重み付け関数オブジェクトの型
*/
template <class W>
class WeightedMedianFilterBase
{
  public:
    using weight_type	= typename W::result_type;	//!< 重みの型
    using warray2_type	= Array2<weight_type>;
    using warray_type	= decltype(*std::declval<warray2_type>().cbegin());

  protected:
  //! 入力信号の各量子化レベルに対して，対応するガイド信号の各量子化レベル毎の頻度
    class Histogram : public Array<size_t>
    {
      public:
	void	resize(size_t nbinsG)
		{
		    Array<size_t>::resize(nbinsG);
		    Array<size_t>::operator =(0);
		    _n = 0;
		}
	void	add(size_t idxG)		{ ++(*this)[idxG]; ++_n; }
	void	remove(size_t idxG)		{ --(*this)[idxG]; --_n; }
	auto	npoints()		const	{ return _n; }
	    
      private:
	size_t	_n = 0;		// ヒストグラムに登録されている点の総数
    };

  //! ガイド信号の各量子化レベルに対して，cut point前後の入力信号量子化レベルの頻度の差
    class BalanceCountingBox : public boost::intrusive::list_base_hook<>
    {
      public:
	void	clear()				{ _diff = 0; _n = 0; }
		operator ptrdiff_t()	const	{ return _diff; }
	auto	add(bool low)
		{
		    if (low)
			++_diff;
		    else
			--_diff;
		    return ++_n;
		}
	auto	remove(bool low)
		{
		    if (low)
			--_diff;
		    else
			++_diff;
		    return --_n;
		}
	auto	operator +=(ptrdiff_t m)	{ return _diff += m; }
	auto	operator -=(ptrdiff_t m)	{ return _diff -= m; }
	    
      private:
	ptrdiff_t	_diff = 0;	// cut point前後の点の数の差
	size_t		_n    = 0;	// この box 中の点の数
    };

  //! 重み付けメディアン値を与えるcut pointの探索器
    class MedianTracker
    {
      private:
	using	nonzero_boxes_t	= boost::intrusive::list<BalanceCountingBox>;
	
      public:
		MedianTracker()	:_median(0)				{}
		MedianTracker(size_t nbinsI, size_t nbinsG)
		    :_histograms(nbinsI), _boxes(nbinsG), _median(0)
		{
		    for (auto& hist : _histograms)
			hist.resize(nbinsG);
		}
	
	void	initialize(size_t nbinsI, size_t nbinsG)
		{
		    _histograms.resize(nbinsI);
		    for (auto& hist : _histograms)
			hist.resize(nbinsG);

		    _boxes.resize(nbinsG);
		    for (auto& box : _boxes)
			box.clear();
		    _nonzero_boxes.clear();
		    _median = 0;
		}

	void	add(size_t idxI, size_t idxG)
		{
		    _histograms[idxI].add(idxG);
		    
		    auto&	box = _boxes[idxG];
		    if (box.add(idxI < _median) == 1)
			_nonzero_boxes.push_back(box);
		}

	void	remove(size_t idxI, size_t idxG)
		{
		    _histograms[idxI].remove(idxG);

		    auto&	box = _boxes[idxG];
		    if (box.remove(idxI < _median) == 0)
			_nonzero_boxes.erase(_nonzero_boxes.iterator_to(box));
		}

	template <class IDX_I, class IDX_G>
	IDX_G	add(IDX_I idxI, IDX_I idxIe, IDX_G idxG)
		{
		    for (; idxI != idxIe; ++idxI, ++idxG)
			add(*idxI, *idxG);
		    return idxG;
		}

	template <class IDX_I, class IDX_G>
	IDX_G	remove(IDX_I idxI, IDX_I idxIe, IDX_G idxG)
		{
		    for (; idxI != idxIe; ++idxI, ++idxG)
			remove(*idxI, *idxG);
		    return idxG;
		}

	auto	median(const warray_type& weights)
		{
		  // 現在の balance 値を計算する．
		    weight_type	balance = 0;
		    for (const auto& box : _nonzero_boxes)
			balance += box * weights[idx(box)];

		    if (balance >= 0)	// balance >= 0 ならば...
		    {			// balance < 0 となるまで
			do		// cut point を左にシフト
			{
			  // 空でないヒストグラムに遭遇するまで
			  // cut point を左にシフト
			    while (_histograms[--_median].npoints() == 0)
				;

			  // 最右の空でないヒストグラム
			    const auto&	hist = _histograms[_median];

			  // 新たな cut point における balance を再計算
			    balance = 0;
			    for (auto& box : _nonzero_boxes)
			    {
				const auto	idxG = idx(box);
				box	-= 2*hist[idxG];
				balance += box * weights[idxG];
			    }
			} while (balance >= 0);

			return _median;
		    }
		    else		// balance < 0 ならば...
		    {			// balance >= 0 となるまで
			do		// cut point を右にシフト
			{
			  // 空でないヒストグラムに遭遇するまで
			  // cut point を右にシフト
			    while (_histograms[_median].npoints() == 0)
				++_median;

			  // 最左の空でないヒストグラム
			    const auto&	hist = _histograms[_median];
			    ++_median;
			    
			  // 新たな cut point における balance を再計算
			    balance = 0;
			    for (auto& box : _nonzero_boxes)
			    {
				const auto	idxG = idx(box);
				box	+= 2*hist[idxG];
				balance += box * weights[idxG];
			    }
			} while (balance < 0);

			return _median - 1;
		    }
		}
	
      private:
	auto	idx(const BalanceCountingBox& box) const
		{
		    return &box - _boxes.data();
		}
	
      private:
	Array<Histogram>		_histograms;
	Array<BalanceCountingBox>	_boxes;
	nonzero_boxes_t			_nonzero_boxes;
	size_t				_median;
    };

  public:
    WeightedMedianFilterBase(const W& wfunc, size_t winSize,
			     size_t nbinsI, size_t nbinsG)		;

    size_t	winSize()		const	{ return _winSize; }
    size_t	outSize(size_t inSize)	const	{ return inSize + 1 - _winSize; }
    size_t	offset()		const	{ return _winSize/2; }
    size_t	nbinsI()		const	{ return _nbinsI; }
    size_t	nbinsG()		const	{ return _nbinsG; }
    void	setWinSize(size_t w)		{ _winSize = w; }
    void	setNBinsI(size_t nbins)		{ _nbinsI = nbins; }
    void	setNBinsG(size_t nbins)		{ _nbinsG = nbins; }
    void	refreshWeights()		{ _initialized = false; }
    
  protected:
    template <class T>
    void	setWeights(const QuantizerBase<T>& quantizer)		;
    auto	weights(size_t i)	const	{ return _weights[i]; }
    
  private:
    const W&		_wfunc;
    size_t		_winSize;
    size_t		_nbinsI;
    size_t		_nbinsG;
    bool		_initialized;
    Array2<weight_type>	_weights;	// weight function
};

template <class W> inline
WeightedMedianFilterBase<W>::WeightedMedianFilterBase(const W& wfunc,
						      size_t winSize,
						      size_t nbinsI,
						      size_t nbinsG)
    :_wfunc(wfunc), _winSize(winSize),
     _nbinsI(nbinsI), _nbinsG(nbinsG), _initialized(false)
{
}

template <class W> template <class T> void
WeightedMedianFilterBase<W>::setWeights(const QuantizerBase<T>& quantizer)
{
    if (_initialized)
	return;

    _weights.resize(quantizer.size(), quantizer.size());
    for (size_t i = 0; i < _weights.nrow(); ++i)
	for (size_t j = i; j < _weights.ncol(); ++j)
	    _weights[i][j] = _weights[j][i]
			   = _wfunc(quantizer[i], quantizer[j]);

    _initialized = std::is_unsigned<T>::value;
}

}
/************************************************************************
*  class WeightedMedianFilter<T, W>					*
************************************************************************/
//! 1次元重み付けメディアンフィルタを表すクラス
/*!
  \param T	出力信号の要素型
  \param W	重み付け関数オブジェクトの型
*/
template <class T, class W>
class WeightedMedianFilter : public detail::WeightedMedianFilterBase<W>
{
  private:
    using value_type	= T;
    using guide_type	= typename W::argument_type;
    using super		= detail::WeightedMedianFilterBase<W>;
    using		typename super::MedianTracker;

  public:
    WeightedMedianFilter(const W& wfunc=W(), size_t winSize=3,
			 size_t nbinsI=256, size_t nbinsG=256)
	:super(wfunc, winSize, nbinsI, nbinsG)				{}

    using	super::winSize;
    using	super::outSize;
    using	super::offset;
    using	super::nbinsI;
    using	super::nbinsG;
    
    template <class IN, class GUIDE, class OUT>
    void	convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
			 OUT out, bool shift=false)			;

  private:
    Quantizer<value_type>	_quantizerI;
    Quantizer<guide_type>	_quantizerG;
    MedianTracker		_tracker;	// 2D histogram
};
    
template <class T, class W>
template <class IN, class GUIDE, class OUT> void
WeightedMedianFilter<T, W>::convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
				     OUT out, bool shift)
{
    if (std::distance(ib, ie) < winSize())
	return;

    const auto&	indicesI = _quantizerI(ib, ie, nbinsI());  // 入力を量子化
    const auto&	indicesG = _quantizerG(gb, ge, nbinsG());  // ガイドを量子化

    super::setWeights(_quantizerG);	// 重みの2次元lookup tableをセット

    auto	headI = indicesI.begin();
    auto	headG = indicesG.begin();
    auto	tailI = headI;
    auto	idxG  = headG;
    std::advance(tailI, winSize()-1);
    std::advance(idxG,  offset());

  // ウィンドウ初期位置におけるヒストグラムをセット
    _tracker.initialize(_quantizerI.size(), _quantizerG.size());
    auto	tailG = _tracker.add(headI, tailI, headG);

    if (shift)
	std::advance(out, offset());
    
  // median点を探索し，その値を出力
    for (; tailI != indicesI.end(); ++tailI)
    {
	_tracker.add(*tailI, *tailG);		// tail点をヒストグラムに追加
	*out = _quantizerI[_tracker.median(super::weights(*idxG))];
	_tracker.remove(*headI, *headG);	// head点をヒストグラムから除去

	++tailG;
	++headI;
	++headG;
	++idxG;
	++out;
    }
}

/************************************************************************
*  class WeightedMedianFilter2<T, W>					*
************************************************************************/
//! 2次元重み付けメディアンフィルタを表すクラス
/*!
  \param T	出力信号の要素型
  \param W	重み付け関数オブジェクトの型
*/
template <class T, class W>
class WeightedMedianFilter2 : public detail::WeightedMedianFilterBase<W>,
			      public Profiler<ENABLE_PROFILER>
{
  private:
    using value_type	= T;
    using guide_type	= typename W::argument_type;
    using pf_type	= Profiler<ENABLE_PROFILER>;
    using super		= detail::WeightedMedianFilterBase<W>;
    using		typename super::MedianTracker;
#if defined(USE_TBB)
    template <class ROW_I, class ROW_G, class ROW_O>
    class Filter
    {
      public:
	Filter(const WeightedMedianFilter2<T, W>& wmf,
	       ROW_I rowI, ROW_G rowG, ROW_O rowO, bool shift)
	    :_wmf(wmf),
	     _rowI(rowI), _rowG(rowG), _rowO(rowO), _shift(shift)	{}
	    
	void	operator ()(const tbb::blocked_range<size_t>& r) const
		{
		    _wmf.filter(_rowI + r.begin(), _rowI + r.end(),
				_rowG + r.begin(), _rowO + r.begin(), _shift);
		}

      private:
	const WeightedMedianFilter2<T, W>&	_wmf;
	ROW_I					_rowI;
	ROW_G					_rowG;
	ROW_O					_rowO;
	const bool				_shift;
    };

    template <class ROW_I, class ROW_G, class ROW_O>
    Filter<ROW_I, ROW_G, ROW_O>
		makeFilter(ROW_I rowI, ROW_G rowG, ROW_O rowO, bool shift) const
		{
		    return Filter<ROW_I, ROW_G, ROW_O>(*this,
						       rowI, rowG, rowO, shift);
		}
#endif
  // std::reverse_iterator<ITER> はITERが指すオブジェクトへの参照を返すため
  // counting_iterator に対しては使えないので，独自に定義
    template <class ITER>
    class reverse_iterator
	: public boost::iterator_adaptor<reverse_iterator<ITER>,
					 ITER,
					 boost::use_default,
					 boost::use_default,
					 iterator_value<ITER> >
    {
      private:
	using super	= boost::iterator_adaptor<reverse_iterator,
						  ITER,
						  boost::use_default,
						  boost::use_default,
						  iterator_value<ITER> >;

      public:
	using		typename super::reference;
	using		typename super::difference_type;
	friend class	boost::iterator_core_access;

      public:
	reverse_iterator(const ITER& iter)	:super(iter)	{}

      private:
	reference	dereference() const
			{
			    auto	tmp = super::base();
			    return *--tmp;
			}
	void		advance(difference_type n)
			{
			    std::advance(super::base_reference(), -n);
			}
	void		increment()
			{
			    --super::base_reference();
			}
	void		decrement()
			{
			    ++super::base_reference();
			}
	difference_type	distance_to(const reverse_iterator& iter) const
			{
			    return std::distance(iter.base(), super::base());
			}
    };

  public:
    WeightedMedianFilter2(const W& wfunc=W(), size_t winSize=3,
			  size_t nbinsI=256, size_t nbinsG=256)
	:super(wfunc, winSize, nbinsI, nbinsG), pf_type(4),
	 _grainSize(100)						{}

    using	super::winSize;
    using	super::outSize;
    using	super::offset;
    using	super::nbinsI;
    using	super::nbinsG;

    auto	winSizeV()			const	{return winSize();}
    auto	winSizeH()			const	{return winSize();}
    auto	outSizeV(size_t nrow)		const	{return outSize(nrow);}
    auto	outSizeH(size_t ncol)		const	{return outSize(ncol);}
    auto	offsetV()			const	{return offset();}
    auto	offsetH()			const	{return offset();}
	
    template <class IN, class GUIDE, class OUT>
    void	convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
			 OUT out, bool shift=false)	;
    size_t	grainSize()			const	{return _grainSize;}
    void	setGrainSize(size_t gs)			{_grainSize = gs;}
    
  private:
    template <class ROW_I, class ROW_G, class ROW_O>
    void	filter(ROW_I rowI, ROW_I rowIe,
		       ROW_G rowG, ROW_O out, bool shift)	const	;
    template <class ROW_I, class ROW_G, class COL_C, class COL_G, class COL_O>
    void	filterRow(MedianTracker& tracker,
			  ROW_I rowI, ROW_G rowG, COL_C c,
			  COL_G colG, COL_O colO)		const	;
    
  private:
    size_t			_grainSize;
    Quantizer2<value_type>	_quantizerI;
    Quantizer2<guide_type>	_quantizerG;
};

template <class T, class W>
template <class IN, class GUIDE, class OUT> void
WeightedMedianFilter2<T, W>::convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
				      OUT out, bool shift)
{
    if (std::distance(ib, ie) < winSize() || ib->size() < winSize())
	return;

    pf_type::start(0);
    const auto&	indicesI = _quantizerI(ib, ie, nbinsI());  // 入力を量子化
    const auto&	indicesG = _quantizerG(gb, ge, nbinsG());  // ガイドを量子化

    pf_type::start(1);
    super::setWeights(_quantizerG);	// 重みの2次元lookup tableをセット

#if defined(USE_TBB)
    tbb::parallel_for(tbb::blocked_range<size_t>(
			  0, indicesI.size() + 1 - winSize(), _grainSize),
		      makeFilter(indicesI.begin(), indicesG.begin(), out, shift));
#else
    filter(indicesI.begin(), indicesI.end() + 1 - winSize(),
	   indicesG.begin(), out, shift);
#endif
    pf_type::nextFrame();
}

template <class T, class W>
template <class ROW_I, class ROW_G, class ROW_O> void
WeightedMedianFilter2<T, W>::filter(ROW_I rowI, ROW_I rowIe,
				    ROW_G rowG, ROW_O rowO, bool shift) const
{
    using col_iterator	= boost::counting_iterator<size_t>;
    using rcol_iterator	= reverse_iterator<col_iterator>;

    pf_type::start(2);
    auto	endI = rowI;
    std::advance(endI, winSize() - 1);	// ウィンドウの最下行
    auto	midG = rowG;

  // ウィンドウ初期位置におけるヒストグラムをセット
    MedianTracker	tracker(_quantizerI.size(), _quantizerG.size());
    for (auto row = rowI; row != endI; ++row, ++midG)
	tracker.add(row->begin(), row->begin() + winSize() - 1, midG->begin());
    
    pf_type::start(3);
    const auto	mid   = offset();
    const auto	rmid  = winSize() - offset() - 1;
    const auto	midO  = (shift ? mid  : 0);
    const auto	rmidO = (shift ? rmid : 0);
    midG = rowG;
    std::advance(midG, mid);		// ウィンドウの中央行
    if (shift)
	std::advance(rowO, mid);

  // 左から右／右から左に交互に走査してmedian点を探索
    for (bool reverse = false; rowI != rowIe; ++rowI)
    {
	if (!reverse)
	{
	    filterRow(tracker, rowI, rowG, col_iterator(0),
		      midG->begin() + mid, rowO->begin() + midO);
	    reverse = true;
	}
	else
	{
	    filterRow(tracker, rowI, rowG,
		      rcol_iterator(col_iterator(rowI->size())),
		      midG->rbegin() + rmid, rowO->rbegin() + rmidO);
	    reverse = false;
	}

	++rowG;
	++midG;
	++rowO;
    }
}

template <class T, class W>
template <class ROW_I, class ROW_G, class COL_C, class COL_G, class COL_O> void
WeightedMedianFilter2<T, W>::filterRow(MedianTracker& tracker,
				       ROW_I rowI, ROW_G rowG,
				       COL_C head,
				       COL_G colG, COL_O colO) const
{
    auto	endI = rowI;
    std::advance(endI, winSize() - 1);		// ウィンドウ最下行
    auto	endG = rowG;
    std::advance(endG, winSize() - 1);		// ウィンドウ最下行
    auto	tail = head;			// ウィンドウ最左／最右列
    auto	end  = head + winSize() - 1;	// ウィンドウ最右／最左列

  // ウィンドウ最下行の点をヒストグラムに追加
    for (; tail != end; ++tail)
	tracker.add(*(endI->begin() + *tail), *(endG->begin() + *tail));
    
    ++endI;					// 最下行の次
    end = head + rowI->size();			// 列の右端／左端
    
    for (; tail != end; ++head, ++tail)
    {
      // tail点をヒストグラムに追加
	tracker.add(make_vertical_iterator(rowI, *tail),
		    make_vertical_iterator(endI, *tail),
		    make_vertical_iterator(rowG, *tail));

      // median点を検出してその値を出力
	*colO = _quantizerI[tracker.median(super::weights(*colG))];

      // head点をヒストグラムから除去
	tracker.remove(make_vertical_iterator(rowI, *head),
		       make_vertical_iterator(endI, *head),
		       make_vertical_iterator(rowG, *head));

	++colG;
	++colO;
    }

  // ウィンドウ最上行の点をヒストグラムから除去
    for (; head != end; ++head)
	tracker.remove(*(rowI->begin() + *head), *(rowG->begin() + *head));
}
    
}
#endif
