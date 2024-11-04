/*
 *  $Id$
 */
#include <cstdlib>
#include "TU/Warp.h"
#include "TU/ICIA.h"

namespace TU
{
/************************************************************************
*  static functions							*
************************************************************************/
template <class T> Image<T>
warp(const Image<T>& src, double du, double dv, double theta)
{
    Matrix33d	Htinv;
    Htinv[0][0] = Htinv[1][1] = cos(theta);
    Htinv[1][0] = sin(theta);
    Htinv[0][1] = -Htinv[1][0];
    Htinv[2][0] = -du;
    Htinv[2][1] = -dv;
    Htinv[2][2] = 1.0;

    Image<T>	dst(src.width(), src.height());
    Warp	warp;
    warp.initialize(Htinv,
		    src.width(), src.height(), dst.width(), dst.height());
    warp(src.cbegin(), dst.begin());

    return dst;
}

template <class Map, class T> void
registerImages(const Image<T>& src, const Image<T>& dst,
	       size_t u0, size_t v0, size_t w, size_t h,
	       typename Map::element_type thresh, bool newton)
{
    using namespace	std;
    using Parameters	= typename ICIA<Map>::Parameters;

    Parameters	params;
    params.newton	   = newton;
    params.intensityThresh = thresh;
    params.niter_max	   = 200;

  // à íuçáÇÌÇπÇé¿çsÅD
    ICIA<Map>	registration(params);
    registration.initialize(src);
    Map		map;
    auto	err = registration(src, dst, map, u0, v0, w, h);
    cerr << "RMS-err = " << sqrt(err) << endl;
    cerr << map;

    registration.print(cerr);
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

    typedef float	T;

    const double	DegToRad = M_PI / 180.0;
    enum Algorithm	{PROJECTIVE, AFFINE, RIGID};
    Algorithm		algorithm = PROJECTIVE;
    double		du = 3.0, dv = -2.0, theta = DegToRad * 3.0;
    T			thresh = 15.0;
    bool		newton = false;
    size_t		u0 = 0, v0 = 0;
    size_t		w = 0, h = 0;
    extern char		*optarg;
    for (int c; (c = getopt(argc, argv, "ARu:v:t:nU:V:W:H:T:")) != -1; )
	switch (c)
	{
	  case 'A':
	    algorithm = AFFINE;
	    break;
	  case 'R':
	    algorithm = RIGID;
	    break;
	  case 'u':
	    du = atof(optarg);
	    break;
	  case 'v':
	    dv = atof(optarg);
	    break;
	  case 't':
	    theta = DegToRad * atof(optarg);
	    break;
	  case 'n':
	    newton = true;
	    break;
	  case 'U':
	    u0 = atoi(optarg);
	    break;
	  case 'V':
	    v0 = atoi(optarg);
	    break;
	  case 'W':
	    w = atoi(optarg);
	    break;
	  case 'H':
	    h = atoi(optarg);
	    break;
	  case 'T':
	    thresh = atof(optarg);
	    break;
	}

    try
    {
	cerr << "Restoring image...";
	Image<u_char>	src;
	src.restore(cin);
	cerr << "done." << endl;

	if (u0 >= src.width())
	    u0 = 0;
	if (v0 >= src.height())
	    v0 = 0;

 	if (w == 0)
	    w = src.width();
	if (h == 0)
	    h = src.height();
	if (u0 + w > src.width())
	    w = src.width() - u0;
	if (v0 + h > src.height())
	    h = src.height() - v0;

	cerr << "Warping image...";
	Image<u_char>	dst = warp(src, du, dv, theta);
	cerr << "done." << endl;

	switch (algorithm)
	{
	  case RIGID:
	    registerImages<Rigidity2<T> >(src, dst, u0, v0, w, h,
					 thresh, newton);
	    break;
	  case AFFINE:
	    registerImages<Affinity22<T> >(src, dst, u0, v0, w, h,
					  thresh, newton);
	    break;
	  default:
	    registerImages<Homography<T> >(src, dst, u0, v0, w, h,
					   thresh, newton);
	    break;
	}
    }
    catch (exception& err)
    {
	cerr << err.what() << endl;
	return 1;
    }

    return 0;
}
