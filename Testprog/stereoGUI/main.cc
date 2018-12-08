/*
 *  $Id: main.cc 1246 2012-11-30 06:23:09Z ueshiba $
 */
#include <unistd.h>
#include <stdlib.h>
#include "TU/io.h"
#include "TU/SADStereo.h"
#include "TU/GFStereo.h"
#include "MyCmdWindow.h"

#define DEFAULT_PARAM_FILE	"stereo"
#define DEFAULT_CONFIG_DIRS	".:/usr/local/etc"
#define DEFAULT_SCALE		1.0
#define DEFAULT_GRAINSIZE	100

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
    cerr << " stereo options.\n"
	 << "  -p params:      stereo parameter file (default: \""
	 << DEFAULT_PARAM_FILE << "\")\n"
	 << "  -d configDirs:  list of directories for camera {conf|calib} file\n"
	 << "                  (default: \"" << DEFAULT_CONFIG_DIRS << "\")\n"
	 << "  -s scale:       positive scale factor (default: "
	 << DEFAULT_SCALE << ")\n"
	 << endl;
    cerr << " viewing options.\n"
	 << "  -x:             use OpenGL texture mapping\n"
	 << "  -q parallax:    parallax for stereo viewing (default: off)\n"
	 << endl;
    cerr << " other options.\n"
	 << "  -g grainSize:   grain size for parallel processin (default: "
	 << DEFAULT_GRAINSIZE << ")\n"
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

    typedef SADStereo<float, u_char>	SADStereoType;
    typedef GFStereo<float,  u_char>	GFStereoType;
    
    bool	gfstereo		= false;
    bool	doHorizontalBackMatch	= true;
    bool	doVerticalBackMatch	= true;
    string	paramFile		= DEFAULT_PARAM_FILE;
    string	configDirs		= DEFAULT_CONFIG_DIRS;
    double	scale			= DEFAULT_SCALE;
    bool	textureMapping		= false;
    double	parallax		= -1.0;
    u_int	grainSize		= DEFAULT_GRAINSIZE;
    
  // コマンド行の解析．
    extern char*	optarg;
    for (int c; (c = getopt(argc, argv, "GHVp:d:s:xq:g:h")) != EOF; )
	switch (c)
	{
	  case 'G':
	    gfstereo = true;
	    break;
	  case 'H':
	    doHorizontalBackMatch = false;
	    break;
	  case 'V':
	    doVerticalBackMatch = false;
	    break;
	  case 'p':
	    paramFile = optarg;
	    break;
	  case 'd':
	    configDirs = optarg;
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
	  case 'h':
	    usage(argv[0]);
	    return 1;
	}
    
  // 本当のお仕事．
    try
    {
	v::App		vapp(argc, argv);	// GUIの初期化．

      // OpenGLの設定．
	int		attrs[] = {GLX_RGBA,
				   GLX_RED_SIZE,	1,
				   GLX_GREEN_SIZE,	1,
				   GLX_BLUE_SIZE,	1,
				   GLX_DEPTH_SIZE,	1,
				   GLX_DOUBLEBUFFER,
				   GLX_STEREO,
				   None};
      	if (parallax <= 0.0)
	{
	    const int	nattrs = sizeof(attrs) / sizeof(attrs[0]);
	    attrs[nattrs - 2] = None;
	}
	XVisualInfo*	vinfo = glXChooseVisual(vapp.colormap().display(),
						vapp.colormap().vinfo().screen,
						attrs);
	if (vinfo == 0)
	    throw runtime_error("No appropriate visual!!");

      // ステレオマッチングパラメータの読み込み．
	ifstream	in;
	openFile(in, paramFile, configDirs, ".params");
	
	if (gfstereo)
	{
	    GFStereoType::Parameters	params;
	    params.get(in);
	    params.doHorizontalBackMatch = doHorizontalBackMatch;
	    params.doVerticalBackMatch	 = doVerticalBackMatch;
	    params.grainSize		 = grainSize;

	    params.put(cerr);
	    
	  // GUIのwidgetを作成．
	    v::MyCmdWindow<GFStereoType, u_char, float>
		myWin(vapp, vinfo, textureMapping, parallax, params, scale);

	  // GUIのイベントループ．
	    vapp.run();
	}
	else
	{
	    SADStereoType::Parameters	params;
	    params.get(in);
	    params.doHorizontalBackMatch = doHorizontalBackMatch;
	    params.doVerticalBackMatch	 = doVerticalBackMatch;
	    params.grainSize		 = grainSize;
	
	  // GUIのwidgetを作成．
	    v::MyCmdWindow<SADStereoType, u_char, float>
		myWin(vapp, vinfo, textureMapping, parallax, params, scale);

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
