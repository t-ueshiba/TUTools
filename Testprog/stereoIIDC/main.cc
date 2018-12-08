/*
 *  平成14-19年（独）産業技術総合研究所 著作権所有
 *  
 *  創作者：植芝俊夫
 *
 *  本プログラムは（独）産業技術総合研究所の職員である植芝俊夫が創作し，
 *  （独）産業技術総合研究所が著作権を所有する秘密情報です．創作者によ
 *  る許可なしに本プログラムを使用，複製，改変，第三者へ開示する等の著
 *  作権を侵害する行為を禁止します．
 *  
 *  このプログラムによって生じるいかなる損害に対しても，著作権所有者お
 *  よび創作者は責任を負いません。
 *
 *  Copyright 2002-2007.
 *  National Institute of Advanced Industrial Science and Technology (AIST)
 *
 *  Creator: Toshio UESHIBA
 *
 *  [AIST Confidential and all rights reserved.]
 *  This program is confidential. Any using, copying, changing or
 *  giving any information concerning with this program to others
 *  without permission by the creator are strictly prohibited.
 *
 *  [No Warranty.]
 *  The copyright holders or the creator are not responsible for any
 *  damages in the use of this program.
 *  
 *  $Id: main.cc 1495 2014-02-27 15:07:51Z ueshiba $
 */
/*!
  \mainpage	stereoIIDC
  本プログラムは，IIDC Specification ver. 1.30に準拠したデジタルカメラ
  を用いて2眼または3眼リアルタイムステレオビジョンを実現するソフトウェアである．
  カメラ間の同期がとれていれば，本プロブラムによってリアルタイムに3次元情報を
  復元することができる．
*/ 
#include <cstdlib>
#include "TU/SADStereo.h"
#include "TU/GFStereo.h"
#include "TU/io.h"
#include "TU/v/vIIDC++.h"
#include "MyCmdWindow.h"

#define DEFAULT_PARAM_FILE	"stereo"
#define DEFAULT_SCALE		1.0
#define DEFAULT_GRAINSIZE	100

//! 本プログラムで定義されたクラスおよび関数を収める名前空間
namespace TU
{
/************************************************************************
*  static functions							*
************************************************************************/
static void
usage(const char* s)
{
    using namespace	std;
    
    cerr << "\nDo binocular/trinocular stereo matching.\n"
	 << endl;
    cerr << " Usage: " << s << " [options]\n"
	 << endl;
    cerr << " configuration options.\n"
	 << "  -c cameraName:  prefix of camera {conf|calib} file\n"
	 << "                  (default: \""
	 << IIDCCameraArray::DEFAULT_CAMERA_NAME
	 << "\")\n"
	 << "  -B:             IEEE1394b mode (default: off)\n"
	 << endl;
    cerr << " stereo options.\n"
	 << "  -p params:      stereo parameter file (default: \""
	 << DEFAULT_PARAM_FILE
	 << "\")\n"
	 << "  -s scale:       positive scale factor (default: "
	 << DEFAULT_SCALE
	 << ")\n"
	 << endl;
    cerr << " viewing options.\n"
	 << "  -x:             use OpenGL texture mapping\n"
	 << "  -q parallax:    parallax for stereo viewing (default: off)\n"
	 << endl;
    cerr << " other options.\n"
	 << "  -g grainSize:   grain size for parallel processing (default: "
	 << DEFAULT_GRAINSIZE
	 << ")\n"
	 << "  -h:             print this\n"
	 << endl;
}

}
/************************************************************************
*  global functions							*
************************************************************************/
int
main(int argc, char* argv[])
{
    using namespace	std;
    using namespace	TU;
    
#if defined(HUGE_IMAGE)
    typedef SADStereo<int,  u_short>	SADStereoType;
    typedef GFStereo<float, u_short>	GFStereoType;
#else    
    typedef SADStereo<short, u_char>	SADStereoType;
    typedef GFStereo<float,  u_char>	GFStereoType;
#endif

    bool	gfstereo	= false;
    const char*	cameraName	= IIDCCameraArray::DEFAULT_CAMERA_NAME;
    const char*	configDirs	= "/usr/local/etc";
    string	paramFile	= DEFAULT_PARAM_FILE;
    double	scale		= DEFAULT_SCALE;
    bool	textureMapping	= false;
    double	parallax	= -1.0;
    size_t	grainSize	= DEFAULT_GRAINSIZE;
    bool	sync		= false;
    float	shutter		= 0.002;
    float	frameRate	= 50.0;
    
  // コマンド行の解析．
    extern char*	optarg;
    for (int c; (c = getopt(argc, argv, "Gc:p:s:xq:g:Sh")) != -1; )
	switch (c)
	{
	  case 'G':
	    gfstereo = true;
	    break;
	  case 'c':
	    cameraName = optarg;
	    break;
	  case 'p':
	    paramFile = optarg;
	    break;
	  case 's':
	    scale = atof(optarg);
	    break;
	  case 'x':
	    textureMapping = true;
	    break;
	  case 'q':
	    parallax = atof(optarg);
	    break;
	  case 'g':
	    grainSize = atoi(optarg);
	    break;
	  case 'S':
	    sync = true;
	    break;
	  case 'h':
	    usage(argv[0]);
	    return 1;
	}
    
  // 本当のお仕事．
    try
    {
	v::App		vapp(argc, argv);	// GUIの初期化．
#if defined(DISPLAY_3D)
      // OpenGLの設定．
	int		attrs[] = {GLX_RGBA,
				   GLX_RED_SIZE,	8,
				   GLX_GREEN_SIZE,	8,
				   GLX_BLUE_SIZE,	8,
				   GLX_DEPTH_SIZE,	24,
				   GLX_DOUBLEBUFFER,
				   GLX_STEREO,
				   None};
      	if (parallax <= 0.0)
	{
	    const int	nattrs = sizeof(attrs) / sizeof(attrs[0]);
	    attrs[nattrs - 2] = None;
	}
	XVisualInfo*	vinfo = glXChooseVisual(vapp.colormap().display(),
#  if defined(DEMO)
						1,
#  else
						vapp.colormap().vinfo().screen,
#  endif
						attrs);
	if (vinfo == 0)
	    throw runtime_error("No appropriate visual!!");
#endif
      // IIDCカメラのオープン．
	IIDCCameraArray	cameras(cameraName);
	cameras.restore();

      // カメラの同期設定
	if (sync)
	{
	    for (auto& camera : cameras)
		camera.reset()
		    .setActive(IIDCCamera::TRIGGER_MODE, false)
		    .setTriggerMode(IIDCCamera::Trigger_Mode14)
		    .setActive(IIDCCamera::Strobe_1, true)
		    .setDuration(IIDCCamera::Strobe_1, 2000)
		    .setOutputVoltage(false)
		    .setActive(IIDCCamera::AUTO_EXPOSURE, false)
		    .setAuto(IIDCCamera::AUTO_EXPOSURE, false)
		    .setAuto(IIDCCamera::SHUTTER, false)
		    .setAbsControl(IIDCCamera::SHUTTER, true)
		    .setValue(IIDCCamera::SHUTTER, shutter)
		    .setActive(IIDCCamera::FRAME_RATE, true)
		    .setAuto(IIDCCamera::FRAME_RATE, false)
		    .setAbsControl(IIDCCamera::FRAME_RATE, true)
		    .setValue(IIDCCamera::FRAME_RATE, frameRate)
		    .embedTimestamp(true);
	  /*
	    switch (trigger)
	    {
	      case None:
		for (auto& camera : cameras)
		    camera.setActive(IIDCCamera::TRIGGER_MODE, false);
		break;
	
	      case Internal:
		for (auto& camera : cameras)
		    camera.setActive(IIDCCamera::TRIGGER_MODE, true);
		std::begin(cameras)->setActive(IIDCCamera::TRIGGER_MODE, false)
		    .setOutputVoltage(true);
		break;
		
	      case External:
		for (auto& camera : cameras)
		    camera.setActive(IIDCCamera::TRIGGER_MODE, true);
		break;
	    }
	  */
	}
	
      // ステレオマッチングパラメータの読み込み．
	ifstream	in;
	openFile(in, paramFile, configDirs, ".params");

	if (gfstereo)
	{
	    GFStereoType::Parameters	params;
	    params.get(in);
	    params.grainSize		 = grainSize;
	
	  // GUIのwidgetを作成．
	    v::MyCmdWindow<GFStereoType, IIDCCameraArray>
		myWin(vapp,
#if defined(DISPLAY_3D)
		      vinfo, textureMapping, parallax,
#endif
		      cameras, params, scale, sync);

	  // GUIのイベントループ．
	    vapp.run();
	}
	else
	{
	    SADStereoType::Parameters	params;
	    params.get(in);
	    params.grainSize		 = grainSize;
	
	  // GUIのwidgetを作成．
	    v::MyCmdWindow<SADStereoType, IIDCCameraArray>
		myWin(vapp,
#if defined(DISPLAY_3D)
		      vinfo, textureMapping, parallax,
#endif
		      cameras, params, scale, sync);

	  // GUIのイベントループ．
	    vapp.run();
	}
    }
    catch (exception& err)
    {
	cerr << err.what() << endl;
	return 1;
    }

    return 0;
}
