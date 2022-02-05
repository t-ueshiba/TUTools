/*
 *  $Id: main.cc 1246 2012-11-30 06:23:09Z ueshiba $
 */
#include <cstdlib>
#include "TU/algorithm.h"
#include "TU/io.h"
#include "TU/Rectify.h"
#include "TU/GuidedFilter.h"
#include "TU/TreeFilter.h"
#include "TU/StereoUtility.h"
#include "TU/Profiler.h"

#define DEFAULT_PARAM_FILE	"stereo"
#define DEFAULT_CONFIG_DIRS	".:/usr/local/etc"
#define DEFAULT_SCALE		1.0
#define DEFAULT_GRAINSIZE	50

namespace TU
{
#if defined(SIMD)
template <class T>	using allocator	= simd::allocator<T>;
#else
template <class T>	using allocator	= std::allocator<T>;
#endif

enum Algorithm	{SAD, GF, WMF, TF};
    
/************************************************************************
*  class MyDiff<S, T>							*
************************************************************************/
template <class S, class T>
struct MyDiff
{
    using argument_type	= S;
    using result_type	= T;
    
    result_type	operator ()(argument_type x, argument_type y) const
		{
		    return std::abs(x - y);
		}
};

/************************************************************************
*  static functions							*
************************************************************************/
template <class S, class T> static void
doJob(const Image<T>& imageL, const Image<T>& imageR,
      const StereoParameters& params, double scale, Algorithm algo,
      size_t ntrials)
{
    using disparity_type	= float;
    
  // 画像を平行化する．
    Rectify			rectify;
    Image<T, allocator<T> >	rectifiedImageL, rectifiedImageR;
    rectify.initialize(imageL, imageR, scale,
		       params.disparitySearchWidth, params.disparityMax);
    rectify(imageL, imageR, rectifiedImageL, rectifiedImageR);

  // フィルタの入出力行を指す反復子を作る．
    const auto	rowL  = std::cbegin(rectifiedImageL);
    const auto	rowR  = std::cbegin(rectifiedImageR);
    const auto	rowLe = std::cend(rectifiedImageL);
    const auto	rowRe = std::cend(rectifiedImageR);
    const auto	rowI  = make_range_iterator(
			    make_diff_iterator<S>(params.disparitySearchWidth,
						  params.intensityDiffMax,
						  std::cbegin(*rowL),
						  std::cbegin(*rowR)),
			    std::make_tuple(stride(rowL), stride(rowR)),
			    std::size(*rowL));
    const auto	rowIe = make_range_iterator(
			    make_diff_iterator<S>(params.disparitySearchWidth,
						  params.intensityDiffMax,
						  std::cbegin(*rowLe),
						  std::cbegin(*rowRe)),
			    std::make_tuple(stride(rowLe), stride(rowRe)),
			    std::size(*rowLe));

    Image<S, allocator<S> >
		disparityMap(rectify.width(0), rectify.height(0));
    const auto	rowD  = make_range_iterator(
			    make_matching_iterator<S>(
				disparityMap.begin()->begin(),
				params.disparitySearchWidth,
				params.disparityMax,
				params.disparityInconsistency,
				params.doHorizontalBackMatch),
			    stride(disparityMap.begin()),
			    disparityMap.width());
#ifdef SCORE_ARRAY3
    Array3<S, 0, 0, 0, allocator<S> >
		scores(disparityMap.height(), disparityMap.width(),
		       params.disparitySearchWidth);
    const auto	rowO  = std::begin(scores);
#else
    const auto	rowO  = rowD;
#endif
  // ステレオマッチングを行う．
    Profiler<>		profiler(2);
    switch (algo)
    {
      case SAD:
      {
	BoxFilter2<S>	filter(params.windowSize, params.windowSize);

	for (size_t i = 0; i < ntrials; ++i)
	{
	    for (size_t j = 0; j < 10; ++j)
	    {
		profiler.start(0);
		rectify(imageL, imageR, rectifiedImageL, rectifiedImageR);
		profiler.start(1);
		filter.convolve(rowI, rowIe, rowO);
	        profiler.nextFrame();
	    }
	    std::cerr << "-------------------------------------------" << std::endl;
	    profiler.print(std::cerr);
	}
      }
	break;

      case GF:
      {
	GuidedFilter2<S>
	    filter(params.windowSize, params.windowSize, params.sigma);
	
	for (size_t i = 0; i < ntrials; ++i)
	{
	    for (size_t j = 0; j < 10; ++j)
	    {
		profiler.start(0);
		rectify(imageL, imageR, rectifiedImageL, rectifiedImageR);
		profiler.start(1);
		filter.convolve(rowI, rowIe, rowL, rowLe, rowO);
		profiler.nextFrame();
	    }
	    std::cerr << "-------------------------------------------" << std::endl;
	    profiler.print(std::cerr);
	}
      }
	break;

      case TF:
      {
	using wfunc_type = MyDiff<T, S>;

	boost::TreeFilter<Array<S>, wfunc_type>	filter(wfunc_type(),
						       params.sigma);

	for (size_t i = 0; i < ntrials; ++i)
	{
	    for (size_t j = 0; j < 10; ++j)
	    {
		profiler.start(0);
		rectify(imageL, imageR, rectifiedImageL, rectifiedImageR);
		profiler.start(1);
		filter.convolve(rowI, rowIe, rowL, rowLe, rowO, true);
		profiler.nextFrame();
	    }
	    std::cerr << "-------------------------------------------" << std::endl;
	    profiler.print(std::cerr);
	    filter.print(std::cerr);
	}
      }
	break;

      default:
	break;
    }
#ifdef SCORE_ARRAY3
    std::copy(scores.cbegin(), scores.cend(), rowD);
#endif
    disparityMap.save(std::cout);
}

}
/************************************************************************
*  global functions							*
************************************************************************/
int
main(int argc, char* argv[])
{
    using namespace	TU;

    using pixel_type	= u_char;
    using score_type	= float;
    
    Algorithm	algo			= SAD;
    std::string	paramFile		= DEFAULT_PARAM_FILE;
    std::string	configDirs		= DEFAULT_CONFIG_DIRS;
    double	scale			= DEFAULT_SCALE;
    size_t	windowSize		= 0;
    size_t	disparitySearchWidth	= 0;
    size_t	disparityMax		= 0;
    size_t	grainSize		= DEFAULT_GRAINSIZE;
    size_t	ntrials			= 5;
    
  // コマンド行の解析．
    extern char*	optarg;
    for (int c; (c = getopt(argc, argv, "GMTp:c:s:w:d:m:g:n:")) != EOF; )
	switch (c)
	{
	  case 'G':
	    algo = GF;
	    break;
	  case 'M':
	    algo = WMF;
	    break;
	  case 'T':
	    algo = TF;
	    break;
	  case 'p':
	    paramFile = optarg;
	    break;
	  case 'c':
	    configDirs = optarg;
	    break;
	  case 's':
	    scale = atof(optarg);
	    break;
	  case 'w':
	    windowSize = atoi(optarg);
	    break;
	  case 'd':
	    disparitySearchWidth = atoi(optarg);
	    break;
	  case 'm':
	    disparityMax = atoi(optarg);
	    break;
	  case 'g':
	    grainSize = atoi(optarg);
	    break;
	  case 'n':
	    ntrials = atoi(optarg);
	    break;
	}
    
  // 本当のお仕事．
    try
    {
	std::ifstream		in;
	openFile(in, paramFile, configDirs, ".params");
	StereoParameters	params;
	params.get(in);
	    
	if (windowSize != 0)
	    params.windowSize = windowSize;
	if (disparityMax != 0)
	    params.disparityMax = disparityMax;
	if (disparitySearchWidth != 0)
	    params.disparitySearchWidth = disparitySearchWidth;
	params.grainSize = grainSize;

      // ステレオマッチングパラメータを表示．
	std::cerr << "--- Stereo matching parameters ---\n";
	params.put(std::cerr);

      // 画像を読み込む．
	Image<pixel_type>	imageL, imageR;
	imageL.restore(std::cin);
	imageR.restore(std::cin);

	doJob<score_type>(imageL, imageR, params, scale, algo, ntrials);
    }
    catch (const std::exception& err)
    {
	std::cerr << err.what() << std::endl;
	return 1;
    }

    return 0;
}
