/*
 *  $Id: MatchImage.cc,v 1.1.1.1 2010-10-13 01:35:46 ueshiba Exp $
 */
#include "MatchImage.h"

namespace TU
{
/************************************************************************
*  class MatchImage							*
************************************************************************/
MatchImage&
MatchImage::copy(const Image<u_char>& image, size_t u0, size_t v0)
{
    for (size_t v = 0; v < image.height(); ++v)
	for (size_t u = 0; u < image.width(); ++u)
	    (*this)[v + v0][u + u0] = image[v][u];

    return *this;
}
    
MatchImage&
MatchImage::drawLine(const Point2i& p, const Point2i& q, bool green)
{
    using	std::abs;
    using	std::swap;

    if (p == q)		// Line of zero length.
	return *this;

    Point2i	pp(p), qq(q);
    RGB		rgb((green ? 0 : 255), (green ? 255 : 0), 0);

  // Is line more horizontal than vertical?
    if (abs(pp[0] - qq[0]) > abs(pp[1] - qq[1]))
    {
      // Put points in increasing order by column.
	if (pp[0] > qq[0])
	    swap(pp, qq);
	int	du = qq[0] - pp[0], dv = qq[1] - pp[1];
	for (int u = pp[0]; u <= qq[0]; ++u)
	    (*this)[pp[1] + (u - pp[0]) * dv / du][u] = rgb;
    }
    else
    {
	if (pp[1] > qq[1])
	    swap(pp, qq);
	int	du = qq[0] - pp[0], dv = qq[1] - pp[1];
	for (int v = pp[1]; v <= qq[1]; ++v)
	    (*this)[v][pp[0] + (v - pp[1]) * du / dv] = rgb;
    }

    return *this;
}
    
}
