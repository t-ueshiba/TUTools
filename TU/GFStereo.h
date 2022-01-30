/*!
  \file		GFStereo.h
  \author	Toshio UESHIBA
  \brief	Guided Filterステレオマッチングクラスの定義と実装
*/
#ifndef TU_GFSTEREO_H
#define TU_GFSTEREO_H

#include "TU/StereoBase.h"
#include "TU/BoxFilter.h"

namespace TU
{
/************************************************************************
*  class GFStereo<SCORE, DISP>						*
************************************************************************/
template <class SCORE, class DISP>
class GFStereo : public StereoBase<GFStereo<SCORE, DISP> >
{
  public:
    using Score		= SCORE;
    using Disparity	= DISP;

  private:
    using super		= StereoBase<GFStereo<Score, Disparity> >;
#if defined(SIMD)
    using ScoreVec	= simd::vec<Score>;
    using DisparityVec	= simd::vec<Disparity>;
#else
    using ScoreVec	= Score;
    using DisparityVec	= Disparity;
#endif
    using ScoreVecTuple	= std::tuple<ScoreVec, ScoreVec>;
    using GuideElement	= Array<Score, 2>;
    
    struct init_params
    {
	init_params(Score g)	:_g(g)					{}

	auto	operator ()(ScoreVec p) const
		{
		    return ScoreVecTuple(p, _g * p);
		}
	
      private:
	const ScoreVec	_g;
    };

    struct init_params2 : public init_params
    {
	init_params2(Score g, Score blend)
	    :init_params(g), _blend(blend)				{}
	
	auto	operator ()(const ScoreVecTuple& pp) const
		{
		    return init_params::operator ()(_blend(pp));
		}
	
      private:
	const Blend<ScoreVec>	_blend;
    };
    
    struct update_params
    {
	update_params(Score gn, Score gp)	:_gn(gn), _gp(gp)	{}

	auto	operator ()(ScoreVec pn, ScoreVec pp) const
		{
		    return ScoreVecTuple(pn - pp, _gn * pn - _gp * pp);
		}
	auto	operator ()(const ScoreVecTuple& p) const
		{
		    return (*this)(std::get<0>(p), std::get<1>(p));
		}

      private:
	const ScoreVec	_gn;
	const ScoreVec	_gp;
    };

    struct update_params2 : public update_params
    {
	update_params2(Score gn, Score gp, Score blend)
	    :update_params(gn, gp), _blend(blend)			{}

	auto	operator ()(const std::tuple<ScoreVec, ScoreVec,
					     ScoreVec, ScoreVec>& p) const
		{
		    return update_params::operator ()(
				_blend(std::get<0>(p), std::get<1>(p)),
				_blend(std::get<2>(p), std::get<3>(p)));
		}

      private:
	const Blend<ScoreVec>	_blend;
    };

    struct init_coeffs
    {
	init_coeffs(Score g_avg, Score g_sqavg, Score e)
	    :_g_avg(g_avg), _g_rvar(1/(g_sqavg - g_avg*g_avg + e*e))	{}

	auto	operator ()(const ScoreVecTuple& params) const
		{
		    const auto	a = (std::get<1>(params) -
				     std::get<0>(params)*_g_avg) * _g_rvar;
		    return ScoreVecTuple(a, std::get<0>(params) - a*_g_avg);
		}

      private:
	const ScoreVec	_g_avg;
	const ScoreVec	_g_rvar;
    };

    struct trans_guides
    {
	trans_guides(Score g) :_g(g)					{}

	auto	operator ()(const ScoreVecTuple& coeffs) const
		{
		    return std::get<0>(coeffs) * _g + std::get<1>(coeffs);
		}
	
      private:
	const ScoreVec	_g;
    };

    using ScoreVecArray		= Array<ScoreVec>;
    using ScoreVecArray2	= Array2<ScoreVec>;
    using ScoreVecArray2Tuple	= std::tuple<ScoreVecArray2, ScoreVecArray2>;
    using ScoreVecArray2TupleArray
				= Array<ScoreVecArray2Tuple>;
    using col_siterator		= iterator_t<ScoreVecArray2Tuple>;
    using const_col_siterator	= iterator_t<const ScoreVecArray2Tuple>;
    using const_reverse_col_siterator
				= std::reverse_iterator<const_col_siterator>;
    using row_siterator		= iterator_t<ScoreVecArray2TupleArray>;
    using const_row_siterator	= iterator_t<const ScoreVecArray2TupleArray>;
    using row_sring		= ring_iterator<row_siterator>;
    using const_row_sring	= ring_iterator<const_row_siterator>;
    using const_row_sbox	= box_filter_iterator<const_row_sring>;
    using const_col_sbox	= box_filter_iterator<const_col_siterator>;
    using const_reverse_col_sbox
			= box_filter_iterator<const_reverse_col_siterator>;

    using GuideArray		= Array<GuideElement>;
    using GuideArray2		= Array2<GuideElement>;
    using col_giterator		= typename GuideArray::iterator;
    using const_col_giterator	= typename GuideArray::const_iterator;
    using const_col_gbox	= box_filter_iterator<const_col_giterator>;

    using DisparityArray	= Array<Disparity>;
    using DisparityArray2	= Array2<Disparity>;
    using reverse_col_diterator	= typename DisparityArray::reverse_iterator;

    using FloatArray		= Array<float>;
    using reverse_col_fiterator	= FloatArray::reverse_iterator;
    
    struct Buffers
    {
	void	initialize(size_t N, size_t D, size_t W)		;
	void	initialize(size_t N, size_t D, size_t W, size_t H)	;
	
	ScoreVecArray2Tuple		Q;	// W x D
	GuideArray			F;	// 1 x W
	ScoreVecArray2TupleArray	A;
	DisparityArray			dminL;	// 1 x (W - N + 1)
	FloatArray			delta;	// 1 x (W - N + 1)
	DisparityArray			dminR;	// 1 x (W + D - 1)
	ScoreVecArray			RminR;	// 1 x D
	DisparityArray2			dminV;	// (W - N + 1) x (H + D - 1)
	ScoreVecArray2			RminV;	// (W - N + 1) x D
    };

  public:
    struct Parameters : public super::Parameters
    {
	Parameters()	:windowSize(11), intensityDiffMax(20),
			 derivativeDiffMax(20), blend(0), epsilon(20)	{}

	std::istream&	get(std::istream& in)
			{
			    super::Parameters::get(in);
			    return in >> windowSize >> intensityDiffMax;
			}
	std::ostream&	put(std::ostream& out) const
			{
			    using namespace	std;

			    super::Parameters::put(out);
			    cerr << "  window size:                        ";
			    out << windowSize << endl;
			    cerr << "  maximum intensity difference:       ";
			    out << intensityDiffMax << endl;
			    cerr << "  maximum derivative difference:      ";
			    out << derivativeDiffMax << endl;
			    cerr << "  blend ratio:                        ";
			    out << blend << endl;
			    cerr << "  epsilon for guided filtering:       ";
			    return out << epsilon << endl;
			}
			    
	size_t	windowSize;		//!< ウィンドウのサイズ
	size_t	intensityDiffMax;	//!< 輝度差の最大値
	size_t	derivativeDiffMax;	//!< 輝度勾配差の最大値
	Score	blend;			//!< 輝度差と輝度勾配差の按分率
	Score	epsilon;		//!< guided filterの正則化パラメータ
    };

  public:
    GFStereo()	:super(*this, 7), _params()				{}
    GFStereo(const Parameters& params)
	:super(*this, 7), _params(params)				{}

    const Parameters&
		getParameters()					const	;
    void	setParameters(const Parameters& params)			;
    size_t	getOverlap()					const	;
    template <class ROW, class ROW_D>
    void	match(ROW rowL, ROW rowLe, ROW rowR, ROW_D rowD)	;
    template <class ROW, class ROW_D>
    void	match(ROW rowL, ROW rowLe, ROW rowLlast,
		      ROW rowR, ROW rowV, ROW_D rowD)			;

  private:
    using	super::start;
    using	super::nextFrame;
    using	super::selectDisparities;
    using	super::pruneDisparities;

    template <class COL, class COL_RV>
    void	initializeFilterParameters(COL colL, COL colLe,
					   COL_RV colRV,
					   col_siterator colQ,
					   col_giterator colF)	  const	;
    template <class COL, class COL_RV>
    void	updateFilterParameters(COL colL, COL colLe, COL_RV colRV,
				       COL colLp, COL_RV colRVp,
				       col_siterator colQ,
				       col_giterator colF)	  const	;
    void	initializeFilterCoefficients(const_col_siterator colQ,
					     const_col_siterator colQe,
					     const_col_giterator colF,
					     col_siterator colA)  const	;
    template <class COL, class DMIN_RV, class RMIN_RV>
    void	computeDisparities(const_reverse_col_siterator colB,
				   const_reverse_col_siterator colBe,
				   COL colG,
				   reverse_col_diterator dminL,
				   reverse_col_fiterator delta,
				   DMIN_RV dminRV, RMIN_RV RminRV) const;

  private:
    Parameters					_params;
    typename super::template Pool<Buffers>	_bufferPool;
};

template <class SCORE, class DISP>
inline const typename GFStereo<SCORE, DISP>::Parameters&
GFStereo<SCORE, DISP>::getParameters() const
{
    return _params;
}
    
template <class SCORE, class DISP> inline void
GFStereo<SCORE, DISP>::setParameters(const Parameters& params)
{
    _params = params;
#if defined(SIMD)
    _params.disparitySearchWidth
	= simd::vec<Disparity>::ceil(_params.disparitySearchWidth);
#endif
    if (_params.disparityMax < _params.disparitySearchWidth)
	_params.disparityMax = _params.disparitySearchWidth;
}

template <class SCORE, class DISP> inline size_t
GFStereo<SCORE, DISP>::getOverlap() const
{
    return 2*_params.windowSize - 2;
}
    
template <class SCORE, class DISP> template <class ROW, class ROW_D> void
GFStereo<SCORE, DISP>::match(ROW rowL, ROW rowLe, ROW rowR, ROW_D rowD)
{
    using	std::size;
    
    start(0);
    const size_t	H = std::distance(rowL, rowLe),
			W = (H != 0 ? size(*rowL) : 0),
			N = _params.windowSize,
			D = _params.disparitySearchWidth;
    if (H < 2*N || W < 2*N)			// 充分な行数／列数があるか確認
	return;
    
    auto* const	buffers = _bufferPool.get();	// 各種作業領域を確保
    buffers->initialize(N, D, W);
    
    auto		rowLp = rowL;
    auto		rowRp = rowR;
    row_sring		rowA(buffers->A.begin(), buffers->A.end());
    const_row_sbox	boxB;
    auto		rowG  = rowL;
    const auto		rowL0 = rowL  + N - 1;
    const auto		rowL1 = rowL0 + N - 1;

    for (; rowL != rowLe; ++rowL)
    {
	using	std::cbegin;
	using	std::cend;
	using	std::begin;
	
	start(1);
      // 各左画素に対して視差[0, D)の右画素のそれぞれとの間の相違度を計算し，
      // フィルタパラメータ(= 縦横両方向に積算された相違度(コスト)の総和
      // および画素毎のコストとガイド画素の積和)を初期化
	if (rowL <= rowL0)
	    initializeFilterParameters(cbegin(*rowL), cend(*rowL),
				       cbegin(*rowR),
				       begin(buffers->Q),
				       buffers->F.begin());
	else
	{
	    updateFilterParameters(cbegin(*rowL), cend(*rowL),
				   cbegin(*rowR),
				   cbegin(*rowLp), cbegin(*rowRp),
				   begin(buffers->Q),
				   buffers->F.begin());
	    ++rowLp;
	    ++rowRp;
	}

	if (rowL >= rowL0)	// 最初のN行に対してコストPが計算済みならば...
	{
	    start(2);
	  // さらにコストを横方向に積算してフィルタパラメータを計算し，
	  // それを用いてフィルタ係数を初期化
	    initializeFilterCoefficients(cbegin(buffers->Q),
					 cend(buffers->Q),
					 buffers->F.begin(),
					 begin(*rowA));
	    ++rowA;

	    if (rowL >= rowL1)		// rowL0からN行分のフィルタ係数が
	    {				// 計算済みならば...
		using	std::crbegin;
		using	std::crend;
		
		start(3);
		if (rowL == rowL1)
		    boxB.initialize(rowA - N, N);
		const auto&	B = *boxB;
		
		start(4);
	      // さらにフィルタ係数を横方向に積算して最終的な係数を求め，
	      // それにguide画像を適用してウィンドウコストを求め，それを
	      // 用いてそれぞれ左/右/上画像を基準とした最適視差を計算
	  	buffers->RminR = std::numeric_limits<Score>::max();
		computeDisparities(crbegin(B), crend(B),
				   crbegin(*rowG) + N - 1,
				   buffers->dminL.rbegin(),
				   buffers->delta.rbegin(),
				   buffers->dminR.end() - D + 1,
				   make_dummy_iterator(&(buffers->RminR)));
		++boxB;
		
		start(5);
	      // 左/右基準視差が一致する場合のみ，それをサブピクセル補間して
	      // 視差として書き出す
		selectDisparities(buffers->dminL.cbegin(),
				  buffers->dminL.cend(),
				  buffers->dminR.cbegin(),
				  buffers->delta.cbegin(),
				  begin(*rowD) + N - 1);
	    }
	    
	    ++rowG;	// guide画像と視差画像は左画像よりもN-1行だけ遅れる
	    ++rowD;	// 同上
	}

	++rowR;
    }

    _bufferPool.put(buffers);
    nextFrame();
}
    
template <class SCORE, class DISP> template <class ROW, class ROW_D> void
GFStereo<SCORE, DISP>::match(ROW rowL, ROW rowLe, ROW rowLlast,
			     ROW rowR, ROW rowV, ROW_D rowD)
{
    start(0);
    const size_t	H = std::distance(rowL, rowLe),
			W = (H != 0 ? rowL->size() : 0),
			N = _params.windowSize,
			D = _params.disparitySearchWidth;
    if (H < 2*N || W < 2*N)			// 充分な行数／列数があるか確認
	return;

    auto* const	buffers = _bufferPool.get();
    buffers->initialize(N, D, W, H);		// 各種作業領域を確保

    auto		v = H;
    size_t		cV = std::distance(rowL, rowLlast);
    auto		rowLp = rowL;
    auto		rowRp = rowR;
    auto		cVp = cV;
    row_sring		rowA(buffers->A.begin(), buffers->A.end());
    const_row_sbox	boxB;
    const auto		rowD0 = rowD + N - 1;
    auto		rowG  = rowL;
    const auto		rowL0 = rowL  + N - 1;
    const auto		rowL1 = rowL0 + N - 1;

    for (; rowL != rowLe; ++rowL)
    {
	using	std::cbegin;
	using	std::cend;
	using	std::begin;
	
	--v;
	--cV;
	
	start(1);
      // 各左画素に対して視差[0, D)の右画素のそれぞれとの間の相違度を計算し，
      // フィルタパラメータ(= 縦横両方向に積算された相違度(コスト)の総和
      // および画素毎のコストとガイド画素の積和)を初期化
	if (rowL <= rowL0)
	    initializeFilterParameters(cbegin(*rowL), cend(*rowL),
				       make_zip_iterator(
					   cbegin(*rowR),
					   make_vertical_iterator(rowV, cV)),
				       begin(buffers->Q),
				       buffers->F.begin());
	else
	{
	    updateFilterParameters(cbegin(*rowL), cend(*rowL),
				   make_zip_iterator(
				       cbegin(*rowR),
				       make_vertical_iterator(rowV, cV)),
				   cbegin(*rowLp),
				   make_zip_iterator(
				       cbegin(*rowRp),
				       make_vertical_iterator(rowV, --cVp)),
				   begin(buffers->Q),
				   buffers->F.begin());
	    ++rowLp;
	    ++rowRp;
	}

	if (rowL >= rowL0)	// 最初のN行に対してコストPが計算済みならば...
	{
	    start(2);
	  // さらにコストを横方向に積算してフィルタパラメータを計算し，
	  // それを用いてフィルタ係数を初期化
	    initializeFilterCoefficients(cbegin(buffers->Q),
					 cend(buffers->Q),
					 buffers->F.cbegin(),
					 begin(*rowA));
	    ++rowA;

	    if (rowL >= rowL1)		// rowL0からN行分のフィルタ係数が
	    {				// 計算済みならば...
		using	std::crbegin;
		using	std::crend;

		start(3);
	      // フィルタ係数を縦方向に積算
		if (rowL == rowL1)
		    boxB.initialize(rowA - N, N);
		const auto&	B = *boxB;
		
		start(4);
	      // さらにフィルタ係数を横方向に積算して最終的な係数を求め，
	      // それにguide画像を適用してウィンドウコストを求め，それを
	      // 用いてそれぞれ左/右/上画像を基準とした最適視差を計算
		buffers->RminR = std::numeric_limits<Score>::max();
		computeDisparities(crbegin(B), crend(B),
				   crbegin(*rowG) + N - 1,
				   buffers->dminL.rbegin(),
				   buffers->delta.rbegin(),
				   make_zip_iterator(
				       buffers->dminR.end() - D + 1,
				       make_vertical_iterator(
					   buffers->dminV.end(), v)),
				   make_zip_iterator(
				       make_dummy_iterator(&(buffers->RminR)),
				       buffers->RminV.rbegin()));
		++boxB;

		start(5);
	      // 左/右基準視差が一致する場合のみ，それをサブピクセル補間して
	      // 視差として書き出す
		selectDisparities(buffers->dminL.cbegin(),
				  buffers->dminL.cend(),
				  buffers->dminR.cbegin(),
				  buffers->delta.cbegin(),
				  begin(*rowD) + N - 1);
	    }

	    ++rowG;	// guide画像と視差画像は左画像よりもN-1行だけ遅れる
	    ++rowD;	// 同上
	}

	++rowR;
    }

    if (_params.doVerticalBackMatch)
    {
	start(6);
	rowD = rowD0;
	for (v = H - 2*(N - 1); v-- != 0; )
	{
	    using	std::begin;
	    
	    pruneDisparities(make_vertical_iterator(buffers->dminV.cbegin(), v),
			     make_vertical_iterator(buffers->dminV.cend(),   v),
			     begin(*rowD) + N - 1);
	    ++rowD;
	}
    }
    
    _bufferPool.put(buffers);
    nextFrame();
}

template <class SCORE, class DISP>
template <class COL, class COL_RV> void
GFStereo<SCORE, DISP>::initializeFilterParameters(COL colL, COL colLe,
						  COL_RV colRV,
						  col_siterator colQ,
						  col_giterator colF) const
{
    using pixel_t	= iterator_value<COL>;
    using diff_t	= Diff<pixel_t>;

    if (_params.blend > 0)
    {
	using blend_t	= Blend<ScoreVec>;
	using ddiff_t	= Diff<std::make_signed_t<pixel_t> >;
#if defined(SIMD)
	using qiterator	= simd::cvtup_iterator<
			      assignment_iterator<
				  init_params2, subiterator<col_siterator> > >;
#else
	using qiterator	= assignment_iterator<init_params2,
					      subiterator<col_siterator> >;
#endif
	while (++colL != colLe - 1)
	{
	    ++colRV;
	    ++colQ;
	    ++colF;
	
	    const auto	pixL = *colL;
	    auto	P = make_zip_iterator(
				make_map_iterator(
				    diff_t(pixL, _params.intensityDiffMax),
				    make_col_accessor(colRV)),
				make_map_iterator(
				    ddiff_t(*(colL + 1) - *(colL - 1),
					    _params.derivativeDiffMax),
				    make_map_iterator(
					Minus(),
					make_col_accessor(colRV) + 1,
					make_col_accessor(colRV) - 1)));
	    for (qiterator Q( make_assignment_iterator(
				  init_params2(pixL, _params.blend),
				  begin(*colQ))),
			   Qe(make_assignment_iterator(
				  init_params2(pixL, _params.blend),
				  end(*colQ)));
		 Q != Qe; ++Q, ++P)
		*Q += *P;

	    (*colF)[0] += pixL;
	    (*colF)[1] += pixL * pixL;
	}
    }
    else
    {
#if defined(SIMD)
	using qiterator	= simd::cvtup_iterator<
			      assignment_iterator<
				  init_params, subiterator<col_siterator> > >;
#else
	using qiterator	= assignment_iterator<init_params,
					      subiterator<col_siterator> >;
#endif
	for (; colL != colLe; ++colL)
	{
	    const auto		pixL = *colL;
	    const diff_t	diff(pixL, _params.intensityDiffMax);
	    auto		in = make_col_accessor(colRV);
	
	    for (qiterator Q( make_assignment_iterator(init_params(pixL),
						       begin(*colQ))),
			   Qe(make_assignment_iterator(init_params(pixL),
						       end(*colQ)));
		 Q != Qe; ++Q, ++in)
		*Q += diff(*in);
	    
	    (*colF)[0] += pixL;
	    (*colF)[1] += pixL * pixL;
	
	    ++colRV;
	    ++colQ;
	    ++colF;
	}
    }
}

template <class SCORE, class DISP> template <class COL, class COL_RV> void
GFStereo<SCORE, DISP>::updateFilterParameters(COL colL, COL colLe, COL_RV colRV,
					      COL colLp, COL_RV colRVp,
					      col_siterator colQ,
					      col_giterator colF) const
{
    using pixel_t	= iterator_value<COL>;
    using diff_t	= Diff<pixel_t>;

    if (_params.blend > 0)
    {
	using ddiff_t	= Diff<std::make_signed_t<pixel_t> >;
#if defined(SIMD)
	using qiterator	= simd::cvtup_iterator<
			      assignment_iterator<
				  update_params2,
				  subiterator<col_siterator> > >;
#else
	using qiterator	= assignment_iterator<update_params2,
					      subiterator<col_siterator> >;
#endif
	while (++colL != colLe - 1)
	{
	    ++colRV;
	    ++colLp;
	    ++colRVp;
	    ++colQ;
	    ++colF;
	
	    const auto	pixLp = *colLp;
	    const auto	pixL  = *colL;
	    auto	P = make_zip_iterator(
				make_map_iterator(
				    diff_t(*colL, _params.intensityDiffMax),
				    make_col_accessor(colRV)),
				make_map_iterator(
				    ddiff_t(*(colL + 1) - *(colL - 1),
					    _params.derivativeDiffMax),
				    make_map_iterator(
					Minus(),
					make_col_accessor(colRV) + 1,
					make_col_accessor(colRV) - 1)),
				make_map_iterator(
				    diff_t(*colLp, _params.intensityDiffMax),
				    make_col_accessor(colRVp)),
				make_map_iterator(
				    ddiff_t(*(colLp + 1) - *(colLp - 1),
					    _params.derivativeDiffMax),
				    make_map_iterator(
					Minus(),
					make_col_accessor(colRVp) + 1,
					make_col_accessor(colRVp) - 1)));
	    for (qiterator Q( make_assignment_iterator(
				  update_params2(pixL, pixLp, _params.blend),
				  begin(*colQ))),
			   Qe(make_assignment_iterator(
				  update_params2(pixL, pixLp, _params.blend),
				  end(*colQ)));
		 Q != Qe; ++Q, ++P)
		*Q += *P;

	    (*colF)[0] += (pixL - pixLp);
	    (*colF)[1] += (pixL * pixL - pixLp * pixLp);
	}
    }
    else
    {
#if defined(SIMD)
	using qiterator	= simd::cvtup_iterator<
			      assignment_iterator<
				  update_params, subiterator<col_siterator> > >;
#else
	using qiterator	= assignment_iterator<update_params,
					      subiterator<col_siterator> >;
#endif
	for (; colL != colLe; ++colL)
	{
	    const auto		pixLp = *colLp;
	    const auto		pixL  = *colL;
	    const diff_t	diff_p(pixLp, _params.intensityDiffMax),
				diff_n(pixL,  _params.intensityDiffMax);
	    auto		in_p = make_col_accessor(colRVp);
	    auto		in_n = make_col_accessor(colRV);
	
	    for (qiterator Q( make_assignment_iterator(
				  update_params(pixL, pixLp),
				  begin(*colQ))),
			   Qe(make_assignment_iterator(
				  update_params(pixL, pixLp),
				  end(*colQ)));
		 Q != Qe; ++Q, ++in_p, ++in_n)
		*Q += std::make_tuple(diff_n(*in_n), diff_p(*in_p));

	    (*colF)[0] += (pixL - pixLp);
	    (*colF)[1] += (pixL * pixL - pixLp * pixLp);

	    ++colRV;
	    ++colLp;
	    ++colRVp;
	    ++colQ;
	    ++colF;
	}
    }
}

template <class SCORE, class DISP> void
GFStereo<SCORE, DISP>::initializeFilterCoefficients(const_col_siterator colQ,
						    const_col_siterator colQe,
						    const_col_giterator colF,
						    col_siterator colA) const
{
    const auto	n = _params.windowSize * _params.windowSize;

  // 縦方向に積算したParamsを横方向に積算し，Coeffを初期化する．
    const_col_gbox	boxG(colF, _params.windowSize);
    for (const_col_sbox boxR(colQ, _params.windowSize), boxRe(colQe);
	 boxR != boxRe; ++boxR)
    {
	std::transform(cbegin(*boxR), cend(*boxR), begin(*colA),
		       init_coeffs((*boxG)[0]/n, (*boxG)[1]/n,
				   _params.epsilon));
	++boxG;
	++colA;
    }
}

template <class SCORE, class DISP>
template <class COL, class DMIN_RV, class RMIN_RV> void
GFStereo<SCORE, DISP>::computeDisparities(const_reverse_col_siterator colB,
					  const_reverse_col_siterator colBe,
					  COL colG,
					  reverse_col_diterator dminL,
					  reverse_col_fiterator delta,
					  DMIN_RV dminRV, RMIN_RV RminRV) const
{
    ScoreVecArray	R(size(*colB));
    
  // 評価値を横方向に積算し，最小値を与える視差を双方向に探索する．
    for (const_reverse_col_sbox boxC(colB, _params.windowSize), boxCe(colBe);
	 boxC != boxCe; ++boxC)
    {
	std::transform(cbegin(*boxC), cend(*boxC),
		       R.begin(), trans_guides(*colG));
	++colG;

#if defined(SIMD)
	using miterator = simd::cvtdown_mask_iterator<
			      simd::mask_type<Disparity>,
			      simd::mask_iterator<
				  typename ScoreVecArray::const_iterator,
				  subiterator<RMIN_RV> > >;
#else
	using miterator = mask_iterator<typename ScoreVecArray::const_iterator,
					subiterator<RMIN_RV> >;
#endif
	Idx<DisparityVec>	index;
	auto			dminRVt = make_col_accessor(--dminRV);
#if defined(SIMD)
	miterator	maskRV(make_mask_iterator(R.cbegin(), begin(*RminRV)));
	for (miterator maskRVe(make_mask_iterator(R.cend(), end(*RminRV)));
	     maskRV != maskRVe; ++maskRV)
#else
	miterator	maskRV(R.cbegin(), begin(*RminRV));
	for (miterator maskRVe(R.cend(),   end(*RminRV));
	     maskRV != maskRVe; ++maskRV)
#endif
	{
	    using dvalue_t = decayed_iterator_value<decltype(dminRVt)>;

	  //*dminRVt = select(*maskRV, index, dvalue_t(*dminRVt));
	    *dminRVt = fast_select(*maskRV, index, dvalue_t(*dminRVt));

	    ++dminRVt;
	    ++index;
	}
#if defined(SIMD)
      	const auto	dL = maskRV.base().dL();	// 左画像から見た視差
	const auto	Rb = R.cbegin().base();
#else
      	const auto	dL = maskRV.dL();		// 左画像から見た視差
	const auto	Rb = R.cbegin();
#endif
	const auto	dsw1 = _params.disparitySearchWidth - 1;
	*dminL = dL;
	*delta = (dL == 0 || dL == dsw1 ? 0 :
		  0.5f * float(Rb[dL-1] - Rb[dL+1]) /
		  float(std::max(Rb[dL-1] - Rb[dL], Rb[dL+1] - Rb[dL]) + 1));
	++delta;
	++dminL;
	++RminRV;
    }
}

/************************************************************************
*  class GFStereo<SCORE, DISP>::Buffers					*
************************************************************************/
template <class SCORE, class DISP> void
GFStereo<SCORE, DISP>::Buffers::initialize(size_t N, size_t D, size_t W)
{
#if defined(SIMD)
    const auto	DD = D / ScoreVec::size;
#else
    const auto	DD = D;
#endif
  //Q.resize(W, DD);
  //Q = ScoreVecTuple(0, 0);	// Q(u, *; d)
    std::get<0>(Q).resize(W, DD);	// Q(u, *; d)
    std::get<1>(Q).resize(W, DD);	// Q(u, *; d)
    std::get<0>(Q) = 0;
    std::get<1>(Q) = 0;

    F.resize(W);
    F = 0;

    A.resize(N + 1);
    for (auto& rowA : A)
    {
      //rowA.resize(W - N + 1, DD);
	std::get<0>(rowA).resize(W - N + 1, DD);
	std::get<1>(rowA).resize(W - N + 1, DD);
    }
    
    dminL.resize(W - 2*N + 2);
    delta.resize(dminL.size());
    dminR.resize(dminL.size() + D - 1);
    RminR.resize(DD);
}

template <class SCORE, class DISP> void
GFStereo<SCORE, DISP>::Buffers::initialize(size_t N, size_t D,
					   size_t W, size_t H)
{
    initialize(N, D, W);

    dminV.resize(dminL.size(), H + D - 1);
    RminV.resize(dminL.size(), RminR.size());
    RminV = std::numeric_limits<SCORE>::max();
}

}
#endif	// !TU_GFSTEREO_H
