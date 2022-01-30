/*!
  \file		GuidedFilter.h
  \author	Toshio UESHIBA
  \brief	guided filterに関するクラスの定義と実装
*/
#ifndef TU_GUIDEDFILTER_H
#define TU_GUIDEDFILTER_H

#include "TU/BoxFilter.h"

namespace TU
{
/************************************************************************
*  class GuidedFilter<T>						*
************************************************************************/
//! 1次元guided filterを表すクラス
template <class T>
class GuidedFilter : public BoxFilter<T>
{
  public:
    using element_type	= T;
    
    struct init_params
    {
	template <class IN_, class GUIDE_>
	auto	operator ()(IN_&& p, GUIDE_&& g) const
		{
		    return std::make_tuple(p, g, std::forward<IN_>(p)*g, g*g);
		}
	template <class IN_>
	auto	operator ()(IN_&& p) const
		{
		    return std::make_tuple(p, p*p);
		}
    };

    class init_coeffs
    {
      public:
	init_coeffs(size_t n, T e)	:_n(n), _n_sq_e(n*e*e)		{}

	template <class VAL_>
	auto	operator ()(const std::tuple<VAL_, T, VAL_, T>& params) const
		{
		    using	std::get;
		    
		    VAL_	a = (_n*get<2>(params) -
				     get<0>(params)*get<1>(params))
				  / (_n*(get<3>(params) + _n_sq_e) -
				     get<1>(params)*get<1>(params));
		    VAL_	b = (get<0>(params) - a*get<1>(params))/_n;
		    return std::make_tuple(std::move(a), std::move(b));
		}
	auto	operator ()(const std::tuple<T, T>& params) const
		{
		    using	std::get;
		    
		    const auto	var = _n*get<1>(params)
				    - get<0>(params)*get<0>(params);
		    const auto	a   = var/(var + _n*_n_sq_e);
		    return std::make_tuple(
				a, (get<0>(params) - a*get<0>(params))/_n);
		}
	
      private:
	const size_t	_n;
	const T		_n_sq_e;
    };

    class trans_guides
    {
      public:
	trans_guides(size_t n)	:_n(n)					{}

	template <class GUIDE_, class OUT_, class COEFF_>
	void	operator ()(std::tuple<GUIDE_, OUT_>&& t,
			    const std::tuple<COEFF_, COEFF_>& coeffs) const
		{
		    using	std::get;

		    get<1>(t) = (get<0>(coeffs)*get<0>(t) + get<1>(coeffs))/_n;
		}

      private:
	const size_t	_n;
    };
    
  private:
    using super	= BoxFilter<T>;

  public:
    GuidedFilter(size_t w, T e) :super(w), _e(e)			{}

    using	super::winSize;
    
    auto	epsilon()		const	{ return _e; }
    auto&	setEpsilon(T e)			{ _e = e; return *this; }
    
    template <class IN, class GUIDE, class OUT>
    void	convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
			 OUT out, bool shift=false)		const	;
    template <class IN, class OUT>
    void	convolve(IN ib, IN ie,
			 OUT out, bool shift=false)		const	;

    auto	outSize(size_t inSize) const
		{
		    return inSize + 2 - 2*winSize();
		}
    auto	offset() const
		{
		    return winSize() - 1;
		}
    
  private:
    T		_e;
};

//! 1次元入力データ列と1次元ガイドデータ列にguided filterを適用する
/*!
  \param ib	1次元入力データ列の先頭を示す反復子
  \param ie	1次元入力データ列の末尾の次を示す反復子
  \param gb	1次元ガイドデータ列の先頭を示す反復子
  \param ge	1次元ガイドデータ列の末尾の次を示す反復子
  \param out	guided filterを適用したデータの出力先を示す反復子
  \param shift	true ならば，入力データと対応するよう，出力位置を
		offset() だけシフトする
*/
template <class T> template <class IN, class GUIDE, class OUT> void
GuidedFilter<T>::convolve(IN ib, IN ie,
			  GUIDE gb, GUIDE ge, OUT out, bool shift) const
{
    using coeff_t	= replace_element<iterator_substance<IN>, T>;
    using coeffs_t	= std::tuple<coeff_t, coeff_t>;
    
  // guided filterの2次元係数ベクトルを計算する．
    Array<coeffs_t>	c(super::outSize(std::distance(ib, ie)));
    super::convolve(make_map_iterator(init_params(), ib, gb),
		    make_map_iterator(init_params(), ie, ge),
		    make_assignment_iterator(init_coeffs(winSize(), _e),
					     c.begin()));
    
  // 係数ベクトルの平均値を求め，それによってガイドデータ列を線型変換する．
    std::advance(gb, offset());
    if (shift)
	std::advance(out, offset());
    super::convolve(c.cbegin(), c.cend(),
		    make_assignment_iterator(trans_guides(winSize()), gb, out));
}

//! 1次元入力データ列にguided filterを適用する
/*!
  ガイドデータ列は与えられた1次元入力データ列に同一とする．
  \param ib	1次元入力データ列の先頭を示す反復子
  \param ie	1次元入力データ列の末尾の次を示す反復子
  \param out	guided filterを適用したデータの出力先を示す反復子
  \param shift	true ならば，入力データと対応するよう，出力位置を
		offset() だけシフトする
*/
template <class T> template <class IN, class OUT> void
GuidedFilter<T>::convolve(IN ib, IN ie, OUT out, bool shift) const
{
    using coeff_t	= replace_element<iterator_substance<IN>, T>;
    using coeffs_t	= std::tuple<coeff_t, coeff_t>;
    
  // guided filterの2次元係数ベクトルを計算する．
    Array<coeffs_t>	c(super::outSize(std::distance(ib, ie)));
    super::convolve(make_map_iterator(init_params(), ib),
		    make_map_iterator(init_params(), ie),
		    make_assignment_iterator(c.begin(),
					     init_coeffs(winSize(), _e)));

  // 係数ベクトルの平均値を求め，それによって入力データ列を線型変換する．
    std::advance(ib, offset());
    if (shift)
	std::advance(out, offset());
    super::convolve(c.cbegin(), c.cend(),
		    make_assignment_iterator(trans_guides(winSize()), ib, out));
}

/************************************************************************
*  class GuidedFilter2<T>						*
************************************************************************/
//! 2次元guided filterを表すクラス
template <class T>
class GuidedFilter2 : private BoxFilter2<T>
{
  public:
    using element_type	= T;
    
  private:
    using init_params	= typename GuidedFilter<T>::init_params;
    using init_coeffs	= typename GuidedFilter<T>::init_coeffs;
    using trans_guides	= typename GuidedFilter<T>::trans_guides;
    using super		= BoxFilter2<T>;
    
  public:
    GuidedFilter2(size_t wrow, size_t wcol, T e)
	:super(wrow, wcol), _e(e)					{}

    using	super::winSizeV;
    using	super::winSizeH;
    using	super::grainSize;
    using	super::setWinSizeV;
    using	super::setWinSizeH;
    using	super::setGrainSize;
    
    auto	epsilon()		const	{ return _e; }
    auto&	setEpsilon(T e)			{ _e = e; return *this; }
    
    template <class IN, class GUIDE, class OUT>
    void	convolve(IN ib, IN ie, GUIDE gb, GUIDE ge,
			 OUT out, bool shift=false)		const	;
    template <class IN, class OUT>
    void	convolve(IN ib, IN ie,
			 OUT out, bool shift=false)		const	;
    
    auto	outSizeV(size_t nrow)	const	{return nrow + 2 - 2*winSizeV();}
    auto	outSizeH(size_t ncol)	const	{return ncol + 2 - 2*winSizeH();}
    auto	offsetV()		const	{return winSizeV() - 1;}
    auto	offsetH()		const	{return winSizeH() - 1;}
    
  private:
    T	_e;
};

//! 2次元入力データと2次元ガイドデータにguided filterを適用する
/*!
  \param ib	2次元入力データの先頭の行を示す反復子
  \param ie	2次元入力データの末尾の次の行を示す反復子
  \param gb	2次元ガイドデータの先頭の行を示す反復子
  \param ge	2次元ガイドデータの末尾の次の行を示す反復子
  \param out	guided filterを適用したデータの出力先の先頭行を示す反復子
  \param shift	trueならば，入力データと対応するよう，出力位置を水平/垂直
		方向にそれぞれ offsetH(), offsetV() だけシフトする
*/
template <class T> template <class IN, class GUIDE, class OUT> void
GuidedFilter2<T>::convolve(IN ib, IN ie,
			   GUIDE gb, GUIDE ge, OUT out, bool shift) const
{
    using		std::cbegin;
    using		std::cend;
    using		std::begin;
    using		std::size;
    using coeff_t	= replace_element<
				iterator_substance<
				    iterator_t<
					decayed_iterator_value<IN> > >, T>;
    using coeffs_t	= std::tuple<coeff_t, coeff_t>;

    if (ib == ie)
	return;

    const auto		n = winSizeV() * winSizeH();
    Array2<coeffs_t>	c(super::outSizeV(std::distance(ib, ie)),
			  super::outSizeH(size(*ib)));
    
  // guided filterの2次元係数ベクトルを計算する．
    super::convolve(make_range_iterator(
			make_map_iterator(init_params(),
					  cbegin(*ib),
					  cbegin(*gb)),
			stride(ib, gb), size(*ib)),
		    make_range_iterator(
			make_map_iterator(init_params(),
					  cbegin(*ie),
					  cbegin(*ge)),
			stride(ie, ge), size(*ie)),
		    make_range_iterator(
			make_assignment_iterator(init_coeffs(n, _e),
						 c.begin()->begin()),
			stride(c.begin()), c.ncol()));

  // 係数ベクトルの平均値を求め，それによってガイドデータ列を線型変換する．
    std::advance(gb, offsetV());
    if (shift)
	std::advance(out, offsetV());
    super::convolve(c.cbegin(), c.cend(),
		    make_range_iterator(
			make_assignment_iterator(
			    trans_guides(n),
			    begin(*gb)   + offsetH(),
			    begin(*out)) + (shift ? offsetH() : 0),
			stride(gb, out), size(*out)));
}

//! 2次元入力データにguided filterを適用する
/*!
  ガイドデータは与えられた2次元入力データに同一とする．
  \param ib	2次元入力データの先頭の行を示す反復子
  \param ie	2次元入力データの末尾の次の行を示す反復子
  \param out	guided filterを適用したデータの出力先の先頭行を示す反復子
  \param shift	trueならば，入力データと対応するよう，出力位置を水平/垂直
		方向にそれぞれ offsetH(), offsetV() だけシフトする
*/
template <class T> template <class IN, class OUT> void
GuidedFilter2<T>::convolve(IN ib, IN ie, OUT out, bool shift) const
{
    using		std::cbegin;
    using		std::cend;
    using		std::begin;
    using		std::size;
    using coeff_t	= replace_element<
				iterator_substance<
				    iterator_t<
					decayed_iterator_value<IN> > >, T>;
    using coeffs_t	= std::tuple<coeff_t, coeff_t>;
    
    if (ib == ie)
	return;

    const auto		n = winSizeV() * winSizeH();
    Array2<coeffs_t>	c(super::outSizeV(std::distance(ib, ie)),
			  super::outSizeH(size(*ib)));

  // guided filterの2次元係数ベクトルを計算する．
    super::convolve(make_range_iterator(
			make_map_iterator(init_params(), cbegin(*ib)),
			stride(ib), size(*ib)),
		    make_range_iterator(
			make_map_iterator(init_params(), cbegin(*ie)),
			stride(ie), size(*ie)),
		    make_range_iterator(
			make_assignment_iterator(init_coeffs(n, _e),
						 c.begin()->begin()),
			stride(c.begin()), c.ncol()));

  // 係数ベクトルの平均値を求め，それによって入力データ列を線型変換する．
    std::advance(ib, offsetV());
    if (shift)
	std::advance(out, offsetV());
    super::convolve(c.cbegin(), c.cend(),
		    make_range_iterator(
			make_assignment_iterator(
			    trans_guides(n),
			    begin(*ib)  + offsetH(),
			    begin(*out) + (shift ? offsetH() : 0)),
			stride(ib, out), size(*out)));
}

}
#endif	// !TU_GUIDEDFILTER_H
