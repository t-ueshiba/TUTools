/*!
  \file		Image++.h
  \author	Toshio UESHIBA
  \brief	画素と画像に関連するクラスの定義と実装
*/
#ifndef	TU_IMAGEPP_H
#define	TU_IMAGEPP_H

#include <cstring>		// for memcpy()
#include <boost/operators.hpp>
#include <sys/types.h>		// for u_int, u_char
#include "TU/pair.h"
#include "TU/Manip.h"
#include "TU/Camera++.h"

using s_char	= signed char;	//!< 符号付き8bit整数

namespace TU
{
namespace detail
{
/************************************************************************
*  class detail::ColorConverter						*
************************************************************************/
class ColorConverter
{
  private:
    constexpr static float	_yr = 0.299f;		// ITU-R BT.601, PAL
    constexpr static float	_yb = 0.114f;		// ITU-R BT.601, PAL
    constexpr static float	_yg = 1.0f - _yr - _yb;	// ITU-R BT.601, PAL
    constexpr static float	_ku = 0.4921f;		// ITU-R BT.601, PAL
    constexpr static float	_kv = 0.877314f;	// ITU-R BT.601, PAL
    
  public:
		ColorConverter()					;

    int		r(int y, int v) const
		{
		    return limit(y + _r[v]);
		}
    int		g(int y, int u, int v) const
		{
		    return limit(y - scaleDown(_gu[u] + _gv[v]));
		}
    int		b(int y, int u) const
		{
		    return limit(y + _b[u]);
		}
    template <class T>
    static T	y(int r, int g, int b)
		{
		    return round<T>(_yr*r + _yg*g + _yb*b);
		}
    int		u(int b, int y) const
		{
		    return _u[255 + b - y];
		}
    int		v(int r, int y) const
		{
		    return _v[255 + r - y];
		}
    
  private:
    static int	limit(int val)
		{
		    return (val < 0 ? 0 : val > 255 ? 255 : val);
		}
    static int	scaleUp(float val)
		{
		    return round<int>(val * (1 << 10));
		}
    static int	scaleDown(int val)
		{
		    return val >> 10;
		}
    template <class T>
    static std::enable_if_t<std::is_integral<T>::value, T>
		round(float val)
		{
		    return T(::round(val));
		}
    template <class T>
    static std::enable_if_t<std::is_floating_point<T>::value, T>
		round(float val)
		{
		    return T(val);
		}
    
  private:
    int		_u[255 + 1 + 255];
    int		_v[255 + 1 + 255];

    int		_r[256];
    int		_gu[256];
    int		_gv[256];
    int		_b[256];
};

extern const ColorConverter	colorConverter;
    
/************************************************************************
*  struct detail::[RGB|BGR|RGBA|ARGB|ABGR|BGRA]				*
************************************************************************/
struct RGB
{
    using element_type = u_char;

    constexpr static size_t	size = 3;
    
    RGB(element_type rr, element_type gg, element_type bb,
	element_type aa=255)	:r(rr), g(gg), b(bb)			{}
    
    element_type			r, g, b;
    constexpr static element_type	a = 255;
};

struct BGR
{
    using element_type = u_char;

    constexpr static size_t	size = 3;
    
    BGR(element_type rr, element_type gg, element_type bb,
	element_type aa=255)	:b(bb), g(gg), r(rr)			{}
    
    element_type			b, g, r;
    constexpr static element_type	a = 255;
};

struct RGBA
{
    using element_type = u_char;

    constexpr static size_t	size = 4;
    
    RGBA(element_type rr, element_type gg, element_type bb,
	 element_type aa=255)	:r(rr), g(gg), b(bb), a(aa)		{}
    
    element_type r, g, b, a;
};

struct ABGR
{
    using element_type = u_char;

    constexpr static size_t	size = 4;
    
    ABGR(element_type rr, element_type gg, element_type bb,
	 element_type aa=255)	:a(aa), b(bb), g(gg), r(rr)		{}
    
    element_type a, b, g, r;
};

struct ARGB
{
    using element_type = u_char;

    constexpr static size_t	size = 4;

    ARGB(element_type rr, element_type gg, element_type bb,
	 element_type aa=255)	:a(aa), r(rr), g(gg), b(bb)		{}
    
    element_type a, r, g, b;
};

struct BGRA
{
    using element_type = u_char;

    constexpr static size_t	size = 4;
    
    BGRA(element_type rr, element_type gg, element_type bb,
	 element_type aa=255)	:b(bb), g(gg), r(rr), a(aa)		{}
    
    element_type b, g, r, a;
};
}	// namespace detail

/************************************************************************
*  struct RGB_<E>							*
************************************************************************/
struct YUV444;
    
template <class E>
struct RGB_ : public E, boost::additive<RGB_<E>,
			boost::multiplicative<RGB_<E>, float,
			boost::equality_comparable<RGB_<E> > > >
{
    using	typename E::element_type;
    
    RGB_()				:E(0, 0, 0)			{}
    RGB_(element_type rr, element_type gg, element_type bb,
	 element_type aa=255)		:E(rr, gg, bb, aa)		{}
    template <class E_>
    RGB_(const RGB_<E_>& p)		:E(p.r, p.g, p.b, p.a)		{}
    template <class T_,
	      std::enable_if_t<std::is_convertible<T_, element_type>::value>*
	      = nullptr>
    RGB_(const T_& p)
	:E(element_type(p), element_type(p), element_type(p))		{}
    RGB_(const YUV444& p)						;
    
    using	E::r;
    using	E::g;
    using	E::b;
    using	E::a;

    template <class T_,
	      std::enable_if_t<std::is_arithmetic<T_>::value>* = nullptr>
		operator T_() const
		{
		    return detail::colorConverter.y<T_>(r, g, b);
		}
    RGB_&	operator +=(const RGB_& p)
		{
		    r += p.r; g += p.g; b += p.b;
		    return *this;
		}
    RGB_&	operator -=(const RGB_& p)
		{
		    r -= p.r; g -= p.g; b -= p.b;
		    return *this;
		}
    RGB_&	operator *=(float c)
		{
		    r *= c; g *= c; b *= c;
		    return *this;
		}
    RGB_&	operator /=(float c)
		{
		    r /= c; g /= c; b /= c;
		    return *this;
		}
    bool	operator ==(const RGB_& p) const
		{
		    return (r == p.r && g == p.g && b == p.b && a == p.a);
		}
};

template <class E> inline std::istream&
operator >>(std::istream& in, RGB_<E>& p)
{
    u_int	r, g, b;
    in >> r >> g >> b;
    p.r = r;
    p.g = g;
    p.b = b;

    return in;
}

template <class E> inline std::ostream&
operator <<(std::ostream& out, const RGB_<E>& p)
{
    return out << u_int(p.r) << ' ' << u_int(p.g) << ' ' << u_int(p.b);
}
    
/************************************************************************
*  type aliases RGB, BGR, RGBA, ABGR, ARGB, BGRA			*
************************************************************************/
//! Red, Green, Blue（各8bit）の順で並んだカラー画素
using RGB = RGB_<detail::RGB>;
//! Blue, Green, Red（各8bit）の順で並んだカラー画素
using BGR = RGB_<detail::BGR>;
//! Red, Green, Blue, Alpha（各8bit）の順で並んだカラー画素
using RGBA = RGB_<detail::RGBA>;
//! Alpha, Blue, Green, Red（各8bit）の順で並んだカラー画素
using ABGR = RGB_<detail::ABGR>;
//! Alpha, Red, Green, Blue（各8bit）の順で並んだカラー画素
using ARGB = RGB_<detail::ARGB>;
//! Blue, Green, Red, Alpha（各8bit）の順で並んだカラー画素
using BGRA = RGB_<detail::BGRA>;

/************************************************************************
*  struct YUV444, YUV422, YUYV422, YUV411				*
************************************************************************/
//! U, Y, V（各8bit）の順で並んだカラー画素
struct YUV444
{
    using element_type = u_char;
    constexpr static size_t	size = 3;

    YUV444(element_type yy=0, element_type uu=128, element_type vv=128)
		    :u(uu), y(yy), v(vv)				{}
    template <class E>
    YUV444(const RGB_<E>& p)
		    :y(detail::colorConverter.y<element_type>(p.r, p.g, p.b))
		{
		    u = detail::colorConverter.u(p.b, y);
		    v = detail::colorConverter.v(p.r, y);
		}
    template <class T,
	      std::enable_if_t<std::is_convertible<T, element_type>::value>*
	      = nullptr>
    YUV444(const T& p)	:u(128), y(p), v(128)				{}
    
    template <class T,
	      std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
		operator T()			const	{return T(y);}
    bool	operator ==(const YUV444& yuv)	const	{return (u == yuv.u &&
								 y == yuv.y &&
								 v == yuv.v);}
    bool	operator !=(const YUV444& yuv)	const	{return !(*this==yuv);}

    element_type	u, y, v;
};
    
inline std::istream&
operator >>(std::istream& in, YUV444& yuv)
{
    u_int	y, u, v;
    in >> y >> u >> v;
    yuv.y = y;
    yuv.u = u;
    yuv.v = v;
    
    return in;
}

inline std::ostream&
operator <<(std::ostream& out, const YUV444& yuv)
{
    return out << u_int(yuv.y) << ' ' << u_int(yuv.u) << ' ' << u_int(yuv.v);
}

template <class E> inline
RGB_<E>::RGB_(const YUV444& p)
    :E(detail::colorConverter.r(p.y, p.v),
       detail::colorConverter.g(p.y, p.u, p.v),
       detail::colorConverter.b(p.y, p.u))
{
}

struct YUYV422;

//! [U, Y0], [V, Y1]（各8bit）の順で並んだカラー画素(16bits/pixel)
struct YUV422
{
    using element_type = u_char;
    constexpr static size_t	size = 2;

    YUV422(element_type yy=0, element_type xx=128)	:x(xx), y(yy)	{}
    YUV422(const YUYV422& p)						;

    template <class T,
	      std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
		operator T()			const	{return T(y);}
    bool	operator ==(const YUV422& p)	const	{return (x == p.x &&
								 y == p.y);}
    bool	operator !=(const YUV422& p)	const	{return !(*this == p);}
    
    element_type	x, y;
};

inline std::istream&
operator >>(std::istream& in, YUV422& yuv)
{
    u_int	y, x;
    in >> y >> x;
    yuv.y = y;
    yuv.x = x;
    
    return in;
}

inline std::ostream&
operator <<(std::ostream& out, const YUV422& yuv)
{
    return out << u_int(yuv.y) << ' ' << u_int(yuv.x);
}

//! [Y0, U], [Y1, V]（各8bit）の順で並んだカラー画素(16bits/pixel)
struct YUYV422
{
    using element_type = u_char;
    constexpr static size_t	size = 2;

    YUYV422(element_type yy=0, element_type xx=128)	:y(yy), x(xx)	{}
    YUYV422(const YUV422& p)				:y(p.y), x(p.x)	{}

    template <class T,
	      std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
		operator T()			const	{return T(y);}
    bool	operator ==(const YUYV422& p)	const	{return (y == p.y &&
								 x == p.x);}
    bool	operator !=(const YUYV422& p)	const	{return !(*this == p);}
    
    element_type	y, x;
};

inline std::istream&
operator >>(std::istream& in, YUYV422& yuv)
{
    u_int	y, x;
    in >> y >> x;
    yuv.y = y;
    yuv.x = x;
    
    return in;
}

inline std::ostream&
operator <<(std::ostream& out, const YUYV422& yuv)
{
    return out << u_int(yuv.y) << ' ' << u_int(yuv.x);
}

inline
YUV422::YUV422(const YUYV422& p)  :x(p.x), y(p.y)	{}

//! [U, Y0, Y1], [V, Y2, Y3]（各8bit）の順で並んだカラー画素(12bits/pixel)
struct YUV411
{
    using element_type = u_char;

    constexpr static size_t	size = 3;

    YUV411(element_type yy0=0, element_type yy1=0, element_type xx=128)
	:x(xx), y0(yy0), y1(yy1)					{}
    template <class T,
	      std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
    YUV411(const T& p)	:x(128), y0(p), y1((&p)[1])			{}

    bool	operator ==(const YUV411& p)	const	{return (x  == p.x  &&
								 y0 == p.y0 &&
								 y1 == p.y1);}
    bool	operator !=(const YUV411& p)	const	{return !(*this == p);}
    
    element_type	x, y0, y1;
};

inline std::istream&
operator >>(std::istream& in, YUV411& yuv)
{
    u_int	y0, y1, x;
    in >> y0 >> y1 >> x;
    yuv.y0 = y0;
    yuv.y1 = y1;
    yuv.x  = x;
    
    return in;
}

inline std::ostream&
operator <<(std::ostream& out, const YUV411& yuv)
{
    return out << u_int(yuv.y0) << ' ' << u_int(yuv.y1) << ' ' << u_int(yuv.x);
}

/************************************************************************
*  class pixel_iterator<ITER>						*
************************************************************************/
namespace detail
{
  template <class ITER, class T=iterator_value<ITER> >
  class pixel_proxy
  {
    public:
      using value_type	= T;

      static constexpr size_t	npixels = 1;
      
    public:
      pixel_proxy(const ITER& iter)	:_iter(const_cast<ITER&>(iter))	{}

      template <class ITER_, class T_>
      pixel_proxy&	operator =(const pixel_proxy<ITER_, T_>& proxy)
			{
			    using N = std::integral_constant<
					  size_t,
					  pixel_proxy<ITER_, T_>::npixels>;

			    assign(proxy.value(N()));
			    return *this;
			}
      
      value_type	value(std::integral_constant<size_t, npixels>) const
			{
			    const auto	val = *_iter;
			    ++_iter;
			    return val;
			}
      template <size_t N>
      TU::pair_tree<T, N>
			value(std::integral_constant<size_t, N>) const
			{
			    using N2 = std::integral_constant<size_t, N/2>;
		    
			    const auto	val0 = value(N2());
			    return std::make_pair(val0, value(N2()));
			}

    private:
      template <class T_>
      void		assign(const T_& val)
			{
			    *_iter = val;
			    ++_iter;
			}
      template <class T_>
      void		assign(const std::pair<T_, T_>& val)
			{
			    assign(val.first);
			    assign(val.second);
			}

    private:
      ITER&	_iter;
  };
    
  template <class ITER>
  class pixel_proxy<ITER, YUV422>
  {
    public:
      using value_type		= TU::pair_tree<YUV444, 2>;
      using element_type	= YUV422::element_type;
      
      static constexpr size_t	npixels = 2;
      
    public:
      pixel_proxy(const ITER& iter)	:_iter(const_cast<ITER&>(iter))	{}

      template <class ITER_, class T_>
      pixel_proxy&	operator =(const pixel_proxy<ITER_, T_>& proxy)
			{
			    constexpr size_t
				Np = pixel_proxy<ITER_, T_>::npixels;
			    using N = std::integral_constant<
					  size_t,
					  (npixels > Np ? npixels : Np)>;
		    
			    assign(proxy.value(N()));
			    return *this;
			}

      value_type	value(std::integral_constant<size_t, npixels>) const
			{
			    const auto	val0 = *_iter;
			    const auto	val1 = *(++_iter);
			    ++_iter;
			    return {{val0.y, val0.x, val1.x},
				    {val1.y, val0.x, val1.x}};
			}
      template <size_t N>
      TU::pair_tree<YUV444, N>
			value(std::integral_constant<size_t, N>) const
			{
			    using N2 = std::integral_constant<size_t, N/2>;
		    
			    const auto	val0 = value(N2());
			    return std::make_pair(val0, value(N2()));
			}

    private:
      template <class T_> std::enable_if_t<!is_pair<T_>::value>
			assign(const std::pair<T_, T_>& val)
			{
			    YUV444	val0(val.first);
			    *_iter     = {val0.y, val0.u};
			    *(++_iter) = {element_type(val.second), val0.v};
			    ++_iter;
			}
      template <class T_> std::enable_if_t<is_pair<T_>::value>
			assign(const std::pair<T_, T_>& val)
			{
			    assign(val.first);
			    assign(val.second);
			}

    private:
      ITER&	_iter;
  };

  template <class ITER>
  class pixel_proxy<ITER, YUYV422>
  {
    public:
      using value_type		= TU::pair_tree<YUV444, 2>;
      using element_type	= YUYV422::element_type;
      
      static constexpr size_t	npixels = 2;
      
    public:
      pixel_proxy(const ITER& iter)	:_iter(const_cast<ITER&>(iter))	{}

      template <class ITER_, class T_>
      pixel_proxy&	operator =(const pixel_proxy<ITER_, T_>& proxy)
			{
			    constexpr size_t
				Np = pixel_proxy<ITER_, T_>::npixels;
			    using N = std::integral_constant<
					  size_t,
					  (npixels > Np ? npixels : Np)>;
		    
			    assign(proxy.value(N()));
			    return *this;
			}

      value_type	value(std::integral_constant<size_t, npixels>) const
			{
			    const auto	val0 = *_iter;
			    const auto	val1 = *(++_iter);
			    ++_iter;
			    return {{val0.y, val0.x, val1.x},
				    {val1.y, val0.x, val1.x}};
			}
      template <size_t N>
      TU::pair_tree<YUV444, N>
			value(std::integral_constant<size_t, N>) const
			{
			    using N2 = std::integral_constant<size_t, N/2>;
		    
			    const auto	val0 = value(N2());
			    return std::make_pair(val0, value(N2()));
			}

    private:
      template <class T_> std::enable_if_t<!is_pair<T_>::value>
			assign(const std::pair<T_, T_>& val)
			{
			    YUV444	val0(val.first);
			    *_iter     = {val0.y, val0.u};
			    *(++_iter) = {element_type(val.second), val0.v};
			    ++_iter;
			}
      template <class T_> std::enable_if_t<is_pair<T_>::value>
			assign(const std::pair<T_, T_>& val)
			{
			    assign(val.first);
			    assign(val.second);
			}

    private:
      ITER&	_iter;
  };

  template <class ITER>
  class pixel_proxy<ITER, YUV411>
  {
    public:
      using value_type		= TU::pair_tree<YUV444, 4>;
      using element_type	= YUV411::element_type;
      
      static constexpr size_t	npixels = 4;
      
    public:
      pixel_proxy(const ITER& iter)	:_iter(const_cast<ITER&>(iter))	{}

      template <class ITER_, class T_>
      pixel_proxy&	operator =(const pixel_proxy<ITER_, T_>& proxy)
			{
			    constexpr size_t
				Np = pixel_proxy<ITER_, T_>::npixels;
			    using N = std::integral_constant<
					  size_t,
					  (npixels > Np ? npixels : Np)>;
		    
			    assign(proxy.value(N()));
			    return *this;
			}

      value_type	value(std::integral_constant<size_t, npixels>) const
			{
			    const auto	val0 = *_iter;
			    const auto	val1 = *(++_iter);
			    ++_iter;
			    return {{{val0.y0, val0.x, val1.x},
				     {val0.y1, val0.x, val1.x}},
				    {{val1.y0, val0.x, val1.x},
				     {val1.y1, val0.x, val1.x}}};
			}

    private:
      template <class T_>
      void		assign(const std::pair<std::pair<T_, T_>,
					       std::pair<T_, T_> >& val)
			{
			    YUV444	val0(val.first.first);
			    *_iter     = {val0.y,
					  element_type(val.first.second),
					  val0.u};
			    *(++_iter) = {element_type(val.second.first),
					  element_type(val.second.second),
					  val0.v};
			    ++_iter;
			}

    private:
      ITER&	_iter;
  };
}

template <class ITER>
class pixel_iterator
    : public boost::iterator_adaptor<pixel_iterator<ITER>,
				     ITER,
				     detail::pixel_proxy<ITER>,
				     std::forward_iterator_tag,
				     detail::pixel_proxy<ITER> >
{
  private:
    using super	= boost::iterator_adaptor<pixel_iterator,
					  ITER,
					  detail::pixel_proxy<ITER>,
					  std::forward_iterator_tag,
					  detail::pixel_proxy<ITER> >;
    
  public:
    using reference	= typename super::reference;

    friend class	boost::iterator_core_access;
    
  public:
		pixel_iterator(ITER iter)	:super(iter)		{}

  private:
    reference	dereference() const
		{
		    return reference(super::base());
		}
    void	increment()
		{
		}
    bool	equal(const pixel_iterator& iter) const
		{
		    return super::base() > iter.base() - reference::npixels;
		}
};

template <class ITER> inline pixel_iterator<ITER>
make_pixel_iterator(ITER iter)
{
    return pixel_iterator<ITER>(iter);
}

/************************************************************************
*  Bayer pattern decoding functions					*
************************************************************************/
namespace detail
{
template <class IN, class OUT, class C> void
bayerDecodeRowXGGY(IN inY, IN inYe, IN inXp, IN inXn, OUT out, C X, C Y)
{
    if (inY == inYe)
	return;
    
    auto	x0 = (*inXp + *inXn) >> 1;
    auto	g0 = *inY;
    auto	y1 = *++inY;
    out->*X = x0;
    out->g  = g0;
    out->*Y = y1;
    ++out;
    
    while (++inY != inYe)
    {
	const auto	g2 = *inY;
	out->g  = (*++inXp + *++inXn + g0 + g2) >> 2;
	const auto	x2 = (*++inXp + *++inXn) >> 1;
	out->*X = (x0 + x2) >> 1;
	out->*Y = y1;
	++out;

	const auto	y3 = *++inY;
	out->*X = x2;
	out->g  = g2;
	out->*Y = (y1 + y3) >> 1;
	++out;

	x0 = x2;
	g0 = g2;
	y1 = y3;
    }

    out->*X = x0;
    out->g  = (*++inXp + *++inXn) >> 1;
    out->*Y = y1;
}
    
template <class IN, class OUT, class C> void
bayerDecodeRowGXYG(IN inY, IN inYe, IN inXp, IN inXn, OUT out, C X, C Y)
{
    if (inY == inYe)
	return;
    
    out->g  = (*inXp + *inXn) >> 1;
    auto	y0 = *inY;
    auto	x1 = (*++inXp + *++inXn) >> 1;
    auto	g1 = *++inY;
    out->*X = x1;
    out->*Y = y0;
    ++out;
    
    while (++inY != inYe)
    {
	const auto	y2 = *inY;
	out->*X = x1;
	out->g  = g1;
	out->*Y = (y0 + y2) >> 1;
	++out;

	const auto	g3 = *++inY;
	out->g  = (*++inXp + *++inXn + g1 + g3) >> 2;
	const auto	x3 = (*++inXp + *++inXn) >> 1;
	out->*X = (x1 + x3) >> 1;
	out->*Y = y2;
	++out;

	y0 = y2;
	x1 = x3;
	g1 = g3;
    }

    out->*X = x1;
    out->g  = g1;
    out->*Y = y0;
}
    
template <class IN, class OUT, class C> void
bayerDecodeBorderRowXGGY(IN inX, IN inXe, IN inY, OUT out, C X, C Y)
{
    if (inX == inXe)
	return;
    
    auto	x0 = *inX;
    auto	g1 = *++inX;
    auto	y1 = *++inY;
    out->*X = x0;
    out->g  = g1;
    out->*Y = y1;
    ++out;
    
    while (++inX != inXe)
    {
	const auto	x2 = *inX;
	out->*X = (x0 + x2) >> 1;
	out->g  = g1;
	out->*Y = y1;
	++out;

	const auto	g3 = *++inX;
	++inY;
	const auto	y3 = *++inY;
	out->*X = x2;
	out->g  = (g1 + g3) >> 1;
	out->*Y = (y1 + y3) >> 1;
	++out;

	x0 = x2;
	g1 = g3;
	y1 = y3;
    }

    out->*X = x0;
    out->g  = g1;
    out->*Y = y1;
}

template <class IN, class OUT, class C> void
bayerDecodeBorderRowGXYG(IN inX, IN inXe, IN inY, OUT out, C X, C Y)
{
    if (inX == inXe)
	return;
    
    auto	g0 = *inX;
    auto	x1 = *++inX;
    auto	y0 = *inY;
    out->*X = x1;
    out->g  = g0;
    out->*Y = y0;
    ++out;
    
    while (++inX != inXe)
    {
	const auto	g2 = *inX;
	++inY;
	const auto	y2 = *++inY;
	out->*X = x1;
	out->g  = (g0 + g2) >> 1;
	out->*Y = (y0 + y2) >> 1;
	++out;

	const auto	x3 = *++inX;
	out->*X = (x1 + x3) >> 1;
	out->g  = g2;
	out->*Y = y2;
	++out;

	g0 = g2;
	x1 = x3;
	y0 = y2;
    }

    out->*X = x1;
    out->g  = g0;
    out->*Y = y0;
}

}	// namespace detail
    
template <class IN, class OUT> OUT
bayerDecodeRGGB(IN in, IN ie, OUT out)
{
    using	COLOR = typename iterator_value<OUT>::value_type;

    if (in == ie)
	return out;

    auto	ic = in;
    ++in;
    detail::bayerDecodeBorderRowXGGY(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*in), begin(*out),
				     &COLOR::r, &COLOR::b);
    ++out;
    auto	ip = ic;
    for (++ic; ++in != ie; ++ip, ++ic)
    {
	detail::bayerDecodeRowXGGY(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::r, &COLOR::b);
	++out;
	++ip;
	++ic;
	++in;
	detail::bayerDecodeRowGXYG(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::b, &COLOR::r);
	++out;
    }
    detail::bayerDecodeBorderRowGXYG(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*ip), begin(*out),
				     &COLOR::b, &COLOR::r);

    return ++out;
}

template <class IN, class OUT> OUT
bayerDecodeBGGR(IN in, IN ie, OUT out)
{
    using	COLOR = typename iterator_value<OUT>::value_type;

    if (in == ie)
	return out;

    auto	ic = in;
    ++in;
    detail::bayerDecodeBorderRowXGGY(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*in), begin(*out),
				     &COLOR::b, &COLOR::r);
    ++out;
    auto	ip = ic;
    for (++ic; ++in != ie; ++ip, ++ic)
    {
	detail::bayerDecodeRowXGGY(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::b, &COLOR::r);
	++out;
	++ip;
	++ic;
	++in;
	detail::bayerDecodeRowGXYG(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::r, &COLOR::b);
	++out;
    }
    detail::bayerDecodeBorderRowGXYG(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*ip), begin(*out),
				     &COLOR::r, &COLOR::b);

    return ++out;
}

template <class IN, class OUT> OUT
bayerDecodeGRBG(IN in, IN ie, OUT out)
{
    using	COLOR = typename iterator_value<OUT>::value_type;

    if (in == ie)
	return out;

    auto	ic = in;
    ++in;
    detail::bayerDecodeBorderRowGXYG(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*in), begin(*out),
				     &COLOR::r, &COLOR::b);
    ++out;
    auto	ip = ic;
    for (++ic; ++in != ie; ++ip, ++ic)
    {
	detail::bayerDecodeRowGXYG(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::r, &COLOR::b);
	++out;
	++ip;
	++ic;
	++in;
	detail::bayerDecodeRowXGGY(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::b, &COLOR::r);
	++out;
    }
    detail::bayerDecodeBorderRowXGGY(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*ip), begin(*out),
				     &COLOR::b, &COLOR::r);

    return ++out;
}

template <class IN, class OUT> OUT
bayerDecodeGBRG(IN in, IN ie, OUT out)
{
    using	COLOR = typename iterator_value<OUT>::value_type;

    if (in == ie)
	return out;

    auto	ic = in;
    ++in;
    detail::bayerDecodeBorderRowGXYG(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*in), begin(*out),
				     &COLOR::b, &COLOR::r);
    ++out;
    auto	ip = ic;
    for (++ic; ++in != ie; ++ip, ++ic)
    {
	detail::bayerDecodeRowGXYG(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::b, &COLOR::r);
	++out;
	++ip;
	++ic;
	++in;
	detail::bayerDecodeRowXGGY(std::cbegin(*ic), std::cend(*ic),
				   std::cbegin(*ip), std::cbegin(*in),
				   begin(*out), &COLOR::r, &COLOR::b);
	++out;
    }
    detail::bayerDecodeBorderRowXGGY(std::cbegin(*ic), std::cend(*ic),
				     std::cbegin(*ip), begin(*out),
				     &COLOR::r, &COLOR::b);

    return ++out;
}

/************************************************************************
*  class ImageFormat							*
************************************************************************/
//! 外部記憶に読み書きする際の画素のタイプ
class ImageFormat
{
  public:
    enum Type
    {
	DEFAULT = 0,			//!< same as internal type
	U_CHAR  = 5,			//!< unsigned mono	 8bit/pixel
	RGB_24  = 6,			//!< RGB		24bit/pixel	
	SHORT,				//!< signed mono	16bit/pixel
	INT,				//!< signed mono	32bit/pixel	
	FLOAT,				//!< float mono		32bit/pixel 
	DOUBLE,				//!< double mono	64bit/pixel
	YUV_444,			//!< YUV444		24bit/pixel
	YUV_422,			//!< YUV422		16bit/pixel
	YUYV_422,			//!< YUYV422		16bit/pixel
	YUV_411,			//!< YUV411		12bit/pixel
	BMP_8,				//!< BMP indexed color   8bit/pixel
	BMP_24,				//!< BMP BGR		24bit/pixel
	BMP_32				//!< BMP BGRA		32bit/pixel
    };

    template <class _T, class _DUMMY=void>
    struct defaultType
    {
	constexpr static auto	value = RGB_24;
    };
    template <class _DUMMY>
    struct defaultType<u_char, _DUMMY>
    {
	constexpr static auto	value = U_CHAR;
    };
    template <class _DUMMY>
    struct defaultType<short, _DUMMY>
    {
	constexpr static auto	value = SHORT;
    };
    template <class _DUMMY>
    struct defaultType<int, _DUMMY>
    {
	constexpr static auto	value = INT;
    };
    template <class _DUMMY>
    struct defaultType<float, _DUMMY>
    {
	constexpr static auto	value = FLOAT;
    };
    template <class _DUMMY>
    struct defaultType<double, _DUMMY>
    {
	constexpr static auto	value = DOUBLE;
    };
    template <class _DUMMY>
    struct defaultType<YUV444, _DUMMY>
    {
	constexpr static auto	value = YUV_444;
    };
    template <class _DUMMY>
    struct defaultType<YUV422, _DUMMY>
    {
	constexpr static auto	value = YUV_422;
    };
    template <class _DUMMY>
    struct defaultType<YUYV422, _DUMMY>
    {
	constexpr static auto	value = YUYV_422;
    };
    template <class _DUMMY>
    struct defaultType<YUV411, _DUMMY>
    {
	constexpr static auto	value = YUV_411;
    };

  public:
    ImageFormat(Type type, bool bottomToTop=false, size_t ncolors=0)
	:_type(type), _bottomToTop(bottomToTop), _ncolors(ncolors)	{}

    Type	type()				const	{ return _type; }
    bool	bottomToTop()			const	{ return _bottomToTop; }
    size_t	ncolors()			const	{ return _ncolors; }
    size_t	depth()				const	;
    size_t	nbytesPerRow(size_t width) const
		{
		    const size_t	nbytes = depth() * width / 8;
		    switch (_type)
		    {
		      case BMP_8:
		      case BMP_24:
		      case BMP_32:
			  return 4 * ((nbytes + 3) / 4);
		      default:
			break;
		    }

		    return nbytes;
		}
    size_t	nbytesForPadding(size_t width) const
		{
		    switch (_type)
		    {
		      case BMP_8:
		      case BMP_24:
		      case BMP_32:
		      {
			  const size_t	nbytes = depth() * width / 8;
			  return 4 * ((nbytes + 3) / 4) - nbytes;
		      }
		      default:
			break;
		    }

		    return 0;
		}

  private:
    Type	_type;		//!< 画素の型
    bool	_bottomToTop;	//!< 行が下から上へ収められているならtrue
    size_t	_ncolors;	//!< カラーマップの色数
};

inline size_t
ImageFormat::depth() const
{
    switch (_type)
    {
      case SHORT:
	return 8*sizeof(short);
      case INT:
	return 8*sizeof(int);
      case FLOAT:
	return 8*sizeof(float);
      case DOUBLE:
	return 8*sizeof(double);
      case RGB_24:
	return 8*sizeof(RGB);
      case YUV_444:
	return 8*sizeof(YUV444);
      case YUV_422:
	return 8*sizeof(YUV422);
      case YUYV_422:
	return 8*sizeof(YUYV422);
      case YUV_411:
	return 12;
      case BMP_24:
	return 8*sizeof(BGR);
      case BMP_32:
	return 8*sizeof(BGRA);
      default:
	break;
    }

    return 8;
}

/************************************************************************
*  class ImageBase<IMAGE>						*
************************************************************************/
//! 画素の2次元配列として定義されたあらゆる画像の基底となるクラス
template <class IMAGE>
class ImageBase
{
  protected:
  //! 画像を生成し投影行列と放射歪曲係数を初期化する．
  /*!
    投影行列は
    \f$\TUbeginarray{cc} \TUvec{I}{3\times 3} & \TUvec{0}{} \TUendarray\f$に，
    2つの放射歪曲係数はいずれも0に初期化される．
  */
    ImageBase()
	:P({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}}), d1(0), d2(0)	{}

  public:
    ImageFormat		restoreHeader(std::istream& in)			;
    ImageFormat::Type	saveHeader(std::ostream& out,
				   ImageFormat::Type type
					=ImageFormat::DEFAULT)	const	;

  //! 画像の幅を返す．
  /*!
    \return	画像の幅
  */
    size_t		width() const
			{
			    return static_cast<const IMAGE*>(this)->ncol();
			}

  //! 画像の高さを返す．
  /*!
    \return	画像の高さ
  */
    size_t		height() const
			{
			    return static_cast<const IMAGE*>(this)->nrow();
			}

    size_t		npixelsToBorder(size_t u, size_t v,
					size_t dir)		const	;
    void		resize(size_t h, size_t w,
			       const ImageFormat& format)
			{
			    return static_cast<IMAGE*>(this)->resize(h, w,
								     format);
			}
    ImageFormat::Type	defaultType() const
			{
			    return static_cast<const IMAGE*>(this)
					->defaultType();
			}
    
  private:
    static bool		isBigEndian()
			{
			    u_short	val = 0x0001;
			    return ((*(u_char*)&val) != 0x01);
			}
    ImageFormat		restorePBMHeader(std::istream& in)		;
    ImageFormat		restoreBMPHeader(std::istream& in)		;
    ImageFormat::Type	savePBMHeader(std::ostream& out,
				      ImageFormat::Type type)
								const	;
    ImageFormat::Type	saveBMPHeader(std::ostream& out,
				      ImageFormat::Type type)
								const	;
    
  public:
    Matrix34d	P;	//!< カメラの3x4投影行列
    double	d1;	//!< 放射歪曲の第1係数
    double	d2;	//!< 放射歪曲の第2係数
};

//! 指定された向きに沿った与えられた点から画像境界までの画素数を返す．
/*!
  \param u	始点の横座標
  \param v	始点の縦座標
  \param dir	8隣接方向
  \return	画像境界までの画素数(始点を含む)
*/
template <class IMAGE> inline size_t
ImageBase<IMAGE>::npixelsToBorder(size_t u, size_t v, size_t dir) const
{
    switch (dir % 8)
    {
      case 0:
	return width() - u;
      case 1:
	return std::min(width() - u, height() - v);
      case 2:
	return height() - v;
      case 3:
	return std::min(u + 1, height() - v);
      case 4:
	return u;
      case 5:
	return std::min(u + 1, v + 1);
      case 6:
	return v + 1;
      default:
	break;
    }

    return std::min(width() - u, v + 1);
}

//! 入力ストリームから画像のヘッダを読み込む．
/*!
  \param in	入力ストリーム
  \return	読み込まれた画像の画素のタイプ
*/
template <class IMAGE> ImageFormat
ImageBase<IMAGE>::restoreHeader(std::istream& in)
{
  // Reset calibration parameters.
    P = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};
    d1 = d2 = 0.0;
    
  // Read the first magic character.
    char	c;
    if (!in.get(c))
	return ImageFormat::DEFAULT;

  // Read image header.
    switch (c)
    {
      case 'P':
	break;
      case 'B':
	return restoreBMPHeader(in);
      default:
	throw std::runtime_error("TU::ImageBase::restoreHeader: neighter PBM nor BMP file!!");
    }

    return restorePBMHeader(in);
}

//! 指定した画素タイプで出力ストリームに画像のヘッダを書き出す．
/*!
  \param out	出力ストリーム
  \param type	画素タイプ．ただし，#DEFAULTを指定した場合は，
		この画像オブジェクトの画素タイプで書き出される．
  \return	実際に書き出す場合の画素タイプ．
*/
template <class IMAGE> ImageFormat::Type
ImageBase<IMAGE>::saveHeader(std::ostream& out, ImageFormat::Type type) const
{
    if (type == ImageFormat::DEFAULT)
	type = defaultType();

    switch (type)
    {
      case ImageFormat::BMP_8:
      case ImageFormat::BMP_24:
      case ImageFormat::BMP_32:
	return saveBMPHeader(out, type);
      default:
	break;
    }

    return savePBMHeader(out, type);
}

template <class IMAGE> ImageFormat
ImageBase<IMAGE>::restorePBMHeader(std::istream& in)
{
  // Read pbm type.
    ImageFormat::Type	type;
    int			c;
    in >> c >> std::ws;		// Read pbm type and skip trailing white spaces.
    switch (c)
    {
      case ImageFormat::U_CHAR:
	type = ImageFormat::U_CHAR;
	break;
      case ImageFormat::RGB_24:
	type = ImageFormat::RGB_24;
	break;
      default:
	throw std::runtime_error("TU::ImageBase::restorePBMHeader: unknown pbm type!!");
    }

  // Process comment lines.
    bool	legacy = false;	// legacy style of dist. param. representation
    for (; (c = in.get()) == '#'; in >> skipl)
    {
	char	key[256], val[256];
	in >> key;
	if (!strcmp(key, "DataType:"))		// pixel data type
	{
	    in >> val;
	    if (!strcmp(val, "Short"))
		type = ImageFormat::SHORT;
	    else if (!strcmp(val, "Int"))
		type = ImageFormat::INT;
	    else if (!strcmp(val, "Float"))
		type = ImageFormat::FLOAT;
	    else if (!strcmp(val, "Double"))
		type = ImageFormat::DOUBLE;
	    else if (!strcmp(val, "YUV444"))
		type = ImageFormat::YUV_444;
	    else if (!strcmp(val, "YUV422"))
		type = ImageFormat::YUV_422;
	    else if (!strcmp(val, "YUYV422"))
		type = ImageFormat::YUYV_422;
	    else if (!strcmp(val, "YUV411"))
		type = ImageFormat::YUV_411;
	}
	else if (!strcmp(key, "Endian:"))	// big- or little-endian
	{
	    in >> val;
	    switch (type)
	    {
	      case ImageFormat::SHORT:
	      case ImageFormat::INT:
	      case ImageFormat::FLOAT:
	      case ImageFormat::DOUBLE:
		if (isBigEndian())
		{
		    if (!strcmp(val, "Little"))
			throw std::runtime_error("TU::ImageBase::restore_epbm: big endian is not supported!!");
		}
		else
		{
		    if (!strcmp(val, "Big"))
			throw std::runtime_error("TU::ImageBase::restore_epbm: little endian is not supported!!");
		}
		break;
	      default:
		break;
	    }
	}
	else if (!strcmp(key, "PinHoleParameterH11:"))
	    in >> P[0][0];
	else if (!strcmp(key, "PinHoleParameterH12:"))
	    in >> P[0][1];
	else if (!strcmp(key, "PinHoleParameterH13:"))
	    in >> P[0][2];
	else if (!strcmp(key, "PinHoleParameterH14:"))
	    in >> P[0][3];
	else if (!strcmp(key, "PinHoleParameterH21:"))
	    in >> P[1][0];
	else if (!strcmp(key, "PinHoleParameterH22:"))
	    in >> P[1][1];
	else if (!strcmp(key, "PinHoleParameterH23:"))
	    in >> P[1][2];
	else if (!strcmp(key, "PinHoleParameterH24:"))
	    in >> P[1][3];
	else if (!strcmp(key, "PinHoleParameterH31:"))
	    in >> P[2][0];
	else if (!strcmp(key, "PinHoleParameterH32:"))
	    in >> P[2][1];
	else if (!strcmp(key, "PinHoleParameterH33:"))
	    in >> P[2][2];
	else if (!strcmp(key, "PinHoleParameterH34:"))
	    in >> P[2][3];
	else if (!strcmp(key, "DistortionParameterD1:"))
	    in >> d1;
	else if (!strcmp(key, "DistortionParameterD2:"))
	    in >> d2;
	else if (!strcmp(key, "DistortionParameterA:"))	// legacy dist. param.
	{
	    in >> d1;
	    legacy = true;
	}
	else if (!strcmp(key, "DistortionParameterB:"))	// legacy dist. param.
	{
	    in >> d2;
	    legacy = true;
	}
    }
    in.putback(c);

    if (legacy)
    {
	Camera<Intrinsic<double> >	camera(P);
	double				k = camera.k();
	d1 *= (k * k);
	d2 *= (k * k * k * k);
    }

    size_t	w, h;
    in >> w >> h;
    resize(h, w, type);				// set width & height
    in >> w >> skipl;				// skip MaxValue

    return type;
}

template <class IMAGE> ImageFormat
ImageBase<IMAGE>::restoreBMPHeader(std::istream& in)
{
  // Read pbm type.
    int	c = in.get();				// Read second magic character.
    if (c != 'M')
	throw std::runtime_error("TU::ImageBase::restoreBMPHeader: not a BMP file!!");

  // Read file header.
    const auto	get16 =	[&in]()
			{
			    const int	c0 = in.get(), c1 = in.get();
			    return c0 + (c1 << 8);
			};
    const auto	get32 =	[&in]()
			{
			    const int	c0 = in.get(), c1 = in.get(),
					c2 = in.get(), c3 = in.get();
			    return c0 + (c1 << 8) + (c2 << 16) + (c3 << 24);
			};
    get32();					// Skip bfSize.
    get16();					// Skip bfReserved1.
    get16();					// Skip bfReserved2.
    get32();					// Skip bfOffBits.

  // Read information header.
    bool	bottomToTop = false;
    size_t	ncolors = 0;
    int		w = 0, h = 0, d = 0;
    switch (c = get32())			// Read bcSize or biSize.
    {
      case 12:	// BMPCoreHeader:
	w = get16();				// Read bcWidth.
	h = get16();				// Read bcHeight.
	if (h > 0)
	    bottomToTop = true;
	else
	{
	    h = -h;
	    bottomToTop = false;
	}
	get16();				// Skip bcPlanes.
	switch (d = get16())			// Read bcBitCount.
	{
	  case 1:
	  case 4:
	  case 8:
	    ncolors = (1 << d);
	    break;
	}
	break;

      case 40:	// BMPInfoHeader:
	w = get32();				// Read biWidth.
	h = get32();				// Read biHeight.
	if (h > 0)
	    bottomToTop = true;
	else
	{
	    h = -h;
	    bottomToTop = false;
	}
	get16();				// Skip biPlanes.
	d = get16();				// Read biBitCount.
	if (get32() != 0)			// Read biCompression.
	    throw std::runtime_error("TUImageBase::restoreBMPHeader: compressed BMP file not supported!!");
	get32();				// Skip biSizeImage.
	get32();				// Skip biXPixPerMeter.
	get32();				// Skip biYPixPerMeter.
	if ((ncolors = get32()) == 0)		// Read biClrUsed.
	    switch (d)
	    {
	      case 1:
		ncolors = 2;
		break;
	      case 4:
		ncolors = 16;
		break;
	      case 8:
		ncolors = 256;
		break;
	    }
	
	get32();				// Read biClrImportant.
	break;

      default:	// Illegal information header size:
	throw std::runtime_error("TU::ImageBase::restoreBMPHeader: information header corrupted!!");
    }

  // Set type of the image.
    ImageFormat::Type	type = ImageFormat::DEFAULT;
    switch (d)
    {
      case 8:
	type = ImageFormat::BMP_8;
	break;
      case 24:
	type = ImageFormat::BMP_24;
	break;
      case 32:
	type = ImageFormat::BMP_32;
	break;
      default:
	throw std::runtime_error("TU::ImageBase::restoreBMPHeader: unsupported depth!!");
    }

    ImageFormat	format(type, bottomToTop, ncolors);
    resize(h, w, format);		// Allocate image area of w*h size.
    
    return format;
}

template <class IMAGE> ImageFormat::Type
ImageBase<IMAGE>::savePBMHeader(std::ostream& out,
				ImageFormat::Type type) const
{
    out << 'P';
    switch (type)
    {
      case ImageFormat::RGB_24:
	out << static_cast<int>(ImageFormat::RGB_24) << std::endl;
	break;
      default:
	out << static_cast<int>(ImageFormat::U_CHAR) << std::endl;
	break;
    }

    ImageFormat	format(type);
    const auto	bitsToBytes = [](size_t n){ return (n + 7)/8; };
    out << "# PixelLength: " << bitsToBytes(format.depth()) << std::endl;
    out << "# DataType: ";
    switch (type)
    {
      default:
	out << "Char" << std::endl;
	break;
      case ImageFormat::RGB_24:
	out << "RGB24" << std::endl;
	break;
      case ImageFormat::SHORT:
	out << "Short" << std::endl;
	break;
      case ImageFormat::INT:
	out << "Int" << std::endl;
	break;
      case ImageFormat::FLOAT:
	out << "Float" << std::endl;
	break;
      case ImageFormat::DOUBLE:
	out << "Double" << std::endl;
	break;
      case ImageFormat::YUV_444:
	out << "YUV444" << std::endl;
	break;
      case ImageFormat::YUV_422:
	out << "YUV422" << std::endl;
	break;
      case ImageFormat::YUYV_422:
	out << "YUYV422" << std::endl;
	break;
      case ImageFormat::YUV_411:
	out << "YUV411" << std::endl;
	break;
    }
    out << "# Sign: ";
    switch (type)
    {
      case ImageFormat::SHORT:
      case ImageFormat::INT:
      case ImageFormat::FLOAT:
      case ImageFormat::DOUBLE:
	out << "Signed" << std::endl;
	break;
      default:
	out << "Unsigned" << std::endl;
	break;
    }
    if (isBigEndian())
	out << "# Endian: Big" << std::endl;
    else
	out << "# Endian: Little" << std::endl;
    out << "# PinHoleParameterH11: " << P[0][0] << std::endl
	<< "# PinHoleParameterH12: " << P[0][1] << std::endl
	<< "# PinHoleParameterH13: " << P[0][2] << std::endl
	<< "# PinHoleParameterH14: " << P[0][3] << std::endl
	<< "# PinHoleParameterH21: " << P[1][0] << std::endl
	<< "# PinHoleParameterH22: " << P[1][1] << std::endl
	<< "# PinHoleParameterH23: " << P[1][2] << std::endl
	<< "# PinHoleParameterH24: " << P[1][3] << std::endl
	<< "# PinHoleParameterH31: " << P[2][0] << std::endl
	<< "# PinHoleParameterH32: " << P[2][1] << std::endl
	<< "# PinHoleParameterH33: " << P[2][2] << std::endl
	<< "# PinHoleParameterH34: " << P[2][3] << std::endl;
    if (d1 != 0.0 || d2 != 0.0)
	out << "# DistortionParameterD1: " << d1 << std::endl
	    << "# DistortionParameterD2: " << d2 << std::endl;
    out << width() << ' ' << height() << '\n'
	<< 255 << std::endl;
    
    return type;
}

template <class IMAGE> ImageFormat::Type
ImageBase<IMAGE>::saveBMPHeader(std::ostream& out,
				ImageFormat::Type type) const
{
  // Write BMP magic characters.
    out << "BM";

  // Write file header.
    const auto	put16 =	[&out](int val)
			{
			    const char	c0 = u_int(val) & 0xff,
					c1 = (u_int(val) >> 8) & 0xff;
			    out.put(c0);
			    out.put(c1);
			};
    const auto	put32 =	[&out](int val)
			{
			    const char	c0 = u_int(val) & 0xff,
					c1 = (u_int(val) >>  8) & 0xff,
					c2 = (u_int(val) >> 16) & 0xff,
					c3 = (u_int(val) >> 24) & 0xff;
			    out.put(c0);
			    out.put(c1);
			    out.put(c2);
			    out.put(c3);
			};
    ImageFormat	format(type);
    put32(14 + 40 + 4*format.ncolors() +
	  format.nbytesPerRow(width())*height());	// Write bfSize.
    put16(0);						// Write bfReserved1.
    put16(0);						// Write bfReserved2.
    put32(14 + 40 + 4*format.ncolors());		// Write bfOffBits.

  // Write information header.
    put32(40);						// Write biSize.
    put32(width());					// Write biWidth.
    put32(height());					// Write biHeight.
    put16(1);						// Write biPlanes.
    put16(format.depth());				// Write biBitCount.
    put32(0);						// Write biCompression.
    put32(format.nbytesPerRow(width())*height());	// Write biSizeImage.
    put32(0);						// Write biXPixPerMeter.
    put32(0);						// Write biYPixPerMeter.
    put32(format.ncolors());				// Write biClrUsed.
    put32(0);						// Write biClrImportant.
    
    return type;
}

/************************************************************************
*  class Image<T, ALLOC>						*
************************************************************************/
//! T型の画素を持つ画像を表すクラス
/*!
  \param T	画素の型
  \param ALLOC	アロケータの型
*/
template <class T, class ALLOC=std::allocator<T> >
class Image : public Array2<T, 0, 0, ALLOC>, public ImageBase<Image<T, ALLOC> >
{
  private:
    using super	= Array2<T, 0, 0, ALLOC>;
    using base	= ImageBase<Image>;
    
  public:
    using	typename super::element_type;
    using	typename super::pointer;

  public:
  //! 幅と高さを指定して画像を生成する．
  /*!
    \param w	画像の幅
    \param h	画像の高さ
  */
    explicit Image(size_t w=0, size_t h=0)
	:super(h, w)						{}

  //! 幅と高さとpaddingを指定して画像を生成する．
  /*!
    \param w	画像の幅
    \param h	画像の高さ
    \param unit	1行あたりのバイト数がこの値の倍数になる
  */
    explicit Image(size_t w, size_t h, size_t unit)
	:super(unit, h, w)					{}

  //! 他の配列と同一要素を持つ画像を作る（コピーコンストラクタの拡張）．
  /*!
    コピーコンストラクタを定義しないと自動的に作られてしまうので，
    このコンストラクタがあってもコピーコンストラクタを別個に定義
    しなければならない．
    \param expr	コピー元の配列
  */
    template <class E_, std::enable_if_t<rank<E_>() == 2>* = nullptr>
    Image(const E_& expr)
	:super(expr)						{}

    Image(pointer p, size_t w, size_t h)
	:super(p, h, w)						{}
    
    Image(pointer p, size_t unit, size_t w, size_t h)
	:super(p, unit, h, w)					{}
    
  //! 他の配列を自分に代入する（標準代入演算子の拡張）．
  /*!
    標準代入演算子を定義しないと自動的に作られてしまうので，この代入演算子が
    あっても標準代入演算子を別個に定義しなければならない．
    \param expr		コピー元の配列
    \return		この配列
  */
    template <class E_> std::enable_if_t<rank<E_>() == 2, Image&>
    operator =(const E_& expr)
	{
	    super::operator =(expr);
	    return *this;
	}

    template <class T_> std::enable_if_t<rank<T_>() == 0, Image&>
    operator =(const T_& c)
	{
	    super::operator =(c);
	    return *this;
	}
    
    using	super::begin;
    using	super::end;
    using	super::rbegin;
    using	super::rend;
    using	super::resize;
    using	base::width;
    using	base::height;
    using	base::restoreHeader;
    using	base::saveHeader;
    
  //! 指定された位置の画素にアクセスする．
  /*!
    \param p	画素の位置
    \return	指定された画素
  */
    template <class S>
    const T&		operator ()(const Array<S, 2>& p) const
			{
			    return (*this)[p[1]][p[0]];
			}

  //! 指定された位置の画素にアクセスする．
  /*!
    \param p	画素の位置
    \return	指定された画素
  */
    template <class S>
    T&			operator ()(const Array<S, 2>& p)
			{
			    return (*this)[p[1]][p[0]];
			}
    
  //! 入力ストリームから画像を読み込む．
  /*!
    \param in	入力ストリーム
    \return	inで指定した入力ストリーム
  */
    std::istream&	restore(std::istream& in)
			{
			    return restoreData(in, restoreHeader(in));
			}

  //! 指定した画素の型で出力ストリームに画像を書き出す．
  /*!
    \param out	出力ストリーム
    \param type	画素の型．ただし，#ImageFormat::DEFAULT を指定した場合は，
		この画像オブジェクトの画素の型で書き出される．
    \return	outで指定した出力ストリーム
  */
    std::ostream&	save(std::ostream& out,
			     ImageFormat::Type type=ImageFormat::DEFAULT) const
			{
			    return saveData(out, saveHeader(out, type));
			}
    std::istream&	restoreData(std::istream& in,
				    const ImageFormat& format)		;
    std::ostream&	saveData(std::ostream& out,
				 ImageFormat::Type type
				     =ImageFormat::DEFAULT)	const	;
    void		resize(size_t h, size_t w, const ImageFormat& format)
			{
			    super::resize(h, w);
			}
    ImageFormat::Type	defaultType() const
			{
			    return ImageFormat::defaultType<T>::value;
			}

  private:
    template <class T_>
    std::istream&	restoreRows(std::istream& in,
				    const ImageFormat& format)		;
    template <class T_, class C_>
    std::istream&	restoreAndLookupRows(std::istream& in,
					     const ImageFormat& format)	;
    template <class T_, class C_>
    std::ostream&	saveRows(std::ostream& out,
				 ImageFormat::Type type)	const	;
};

//! 入力ストリームから画像の画素データを読み込む．
/*!
  \param in	入力ストリーム
  \param format	ストリーム中のデータの画素フォーマット
		(読み込み先の画像の画素フォーマットではない)
  \return	inで指定した入力ストリーム
*/
template <class T, class ALLOC> std::istream&
Image<T, ALLOC>::restoreData(std::istream& in, const ImageFormat& format)
{
    switch (format.type())
    {
      case ImageFormat::DEFAULT:
	break;
      case ImageFormat::U_CHAR:
	return restoreRows<u_char >(in, format);
      case ImageFormat::SHORT:
	return restoreRows<short  >(in, format);
      case ImageFormat::INT:
	return restoreRows<int	  >(in, format);
      case ImageFormat::FLOAT:
	return restoreRows<float  >(in, format);
      case ImageFormat::DOUBLE:
	return restoreRows<double >(in, format);
      case ImageFormat::RGB_24:
	return restoreRows<RGB	  >(in, format);
      case ImageFormat::YUV_444:
	return restoreRows<YUV444 >(in, format);
      case ImageFormat::YUV_422:
	return restoreRows<YUV422 >(in, format);
      case ImageFormat::YUYV_422:
	return restoreRows<YUYV422>(in, format);
      case ImageFormat::YUV_411:
	return restoreRows<YUV411 >(in, format);
      case ImageFormat::BMP_8:
	return restoreAndLookupRows<u_char, BGRA>(in, format);
      case ImageFormat::BMP_24:
	return restoreRows<BGR	  >(in, format);
      case ImageFormat::BMP_32:
	return restoreRows<BGRA	  >(in, format);
      default:
	throw std::runtime_error("Image<T, ALLOC>::restoreData(): unknown pixel type!!");
    }

    return in;
}

//! 指定した画素フォーマットで出力ストリームに画像の画素データを書き出す．
/*!
  \param out	出力ストリーム
  \param type	画素の型．ただし，#PixelType::DEFAULT を指定した場合は，
		この画像オブジェクトのデフォルトの画素型で書き出される．
  \return	outで指定した出力ストリーム
*/
template <class T, class ALLOC> std::ostream&
Image<T, ALLOC>::saveData(std::ostream& out, ImageFormat::Type type) const
{
    if (type == ImageFormat::DEFAULT)
	type = defaultType();

    switch (type)
    {
      case ImageFormat::U_CHAR:
	return saveRows<u_char,	 RGB >(out, type);
      case ImageFormat::SHORT:
	return saveRows<short,	 RGB >(out, type);
      case ImageFormat::INT:
	return saveRows<int,	 RGB >(out, type);
      case ImageFormat::FLOAT:
	return saveRows<float,	 RGB >(out, type);
      case ImageFormat::DOUBLE:
	return saveRows<double,	 RGB >(out, type);
      case ImageFormat::RGB_24:
	return saveRows<RGB,	 RGB >(out, type);
      case ImageFormat::YUV_444:
	return saveRows<YUV444,	 RGB >(out, type);
      case ImageFormat::YUV_422:
	return saveRows<YUV422,	 RGB >(out, type);
      case ImageFormat::YUYV_422:
	return saveRows<YUYV422, RGB >(out, type);
      case ImageFormat::YUV_411:
	return saveRows<YUV411,  RGB >(out, type);
      case ImageFormat::BMP_8:
	return saveRows<u_char,	 BGRA>(out, type);
      case ImageFormat::BMP_24:
	return saveRows<BGR,	 RGB >(out, type);
      case ImageFormat::BMP_32:
	return saveRows<BGRA,	 RGB >(out, type);
      default:
	throw std::runtime_error("Image<T, ALLOC>::saveData(): unknown pixel type!!");
    }

    return out;
}

template <class T, class ALLOC> template <class T_> std::istream&
Image<T, ALLOC>::restoreRows(std::istream& in, const ImageFormat& format)
{
    const auto	npads = format.nbytesForPadding(width());
    Array<T_>	buf(width());
    if (format.bottomToTop())
    {
	for (auto row = rbegin(); row != rend(); ++row)
	{
	    if (!buf.restore(in) || !in.ignore(npads))
		break;

	    std::copy(make_pixel_iterator(buf.cbegin()),
		      make_pixel_iterator(buf.cend()),
		      make_pixel_iterator((*row).begin()));
	}
    }
    else
    {
	for (auto row : *this)
	{
	    if (!buf.restore(in) || !in.ignore(npads))
		break;

	    std::copy(make_pixel_iterator(buf.cbegin()),
		      make_pixel_iterator(buf.cend()),
		      make_pixel_iterator(row.begin()));
	}
    }

    return in;
}

template <class T, class ALLOC> template <class T_, class C_> std::istream&
Image<T, ALLOC>::restoreAndLookupRows(std::istream& in,
				      const ImageFormat& format)
{
    Array<C_>	colormap(format.ncolors());
    colormap.restore(in);
	
    const auto	lookup = [&colormap](T_ i){ return colormap[i]; };
    const auto	npads = format.nbytesForPadding(width());
    Array<T_>	buf(width());
    if (format.bottomToTop())
    {
	for (auto row = rbegin(); row != rend(); ++row)
	{
	    if (!buf.restore(in) || !in.ignore(npads))
		break;

	    std::copy(make_pixel_iterator(make_map_iterator(lookup,
							    buf.cbegin())),
		      make_pixel_iterator(make_map_iterator(lookup,
							    buf.cend())),
		      make_pixel_iterator((*row).begin()));
	}
    }
    else
    {
	for (auto row : *this)
	{
	    if (!buf.restore(in) || !in.ignore(npads))
		break;

	    std::copy(make_pixel_iterator(make_map_iterator(lookup,
							    buf.cbegin())),
		      make_pixel_iterator(make_map_iterator(lookup,
							    buf.cend())),
		      make_pixel_iterator(row.begin()));
	}
    }

    return in;
}

template <class T, class ALLOC> template <class T_, class C_> std::ostream&
Image<T, ALLOC>::saveRows(std::ostream& out, ImageFormat::Type type) const
{
    ImageFormat	format(type);
    Array<C_>	colormap(format.ncolors());
    for (size_t i = 0; i < colormap.size(); ++i)
	colormap[i] = i;
    colormap.save(out);
    
    Array<u_char>	pads(format.nbytesForPadding(width()));
    Array<T_>		buf(width());
    if (format.bottomToTop())
    {
	for (auto row = rbegin(); row != rend(); ++row)
	{
	    using	std::begin;
	    using	std::end;
	    
	    std::copy(make_pixel_iterator(begin(*row)),
		      make_pixel_iterator(end(*row)),
		      make_pixel_iterator(buf.begin()));
	    if (!buf.save(out) || !pads.save(out))
		break;
	}
    }
    else
    {
	for (const auto row : *this)
	{
	    using	std::begin;
	    using	std::end;
	    
	    std::copy(make_pixel_iterator(begin(row)),
		      make_pixel_iterator(end(row)),
		      make_pixel_iterator(buf.begin()));
	    if (!buf.save(out) || !pads.save(out))
		break;
	}
    }

    return out;
}

/************************************************************************
*  class GenericImage							*
************************************************************************/
//! 画素の型を問わない総称画像クラス
/*!
  個々の行や画素にアクセスすることはできない．
*/
class GenericImage : public ImageBase<GenericImage>
{
  private:
    using array2_type	= Array2<char>;

  public:
    using pointer	= array2_type::pointer;
    using const_pointer	= array2_type::const_pointer;
    
  public:
  //! 総称画像を生成する．
    GenericImage() :_a(), _format(ImageFormat::U_CHAR), _colormap(0)	{}

    pointer		data()				{ return _a.data(); }
    const_pointer	data()			const	{ return _a.data(); }
    size_t		nrow()			const	{ return _a.nrow(); }
    size_t		ncol() const
			{
			    return (_a.ncol()*8) / _format.depth();
			}

  //! 現在保持している画像のタイプ情報を返す．
  /*!
    \return	タイプ情報
  */
    const ImageFormat&	format()		const	{ return _format; }

  //! 入力ストリームから画像を読み込む．
  /*!
    \param in	入力ストリーム
    \return	inで指定した入力ストリーム
  */
    std::istream&	restore(std::istream& in)
			{
			    restoreHeader(in);
			    return restoreData(in);
			}

  //! 出力ストリームに画像を書き出す．
  /*!
    \param out	出力ストリーム
    \return	outで指定した出力ストリーム
  */
    std::ostream&	save(std::ostream& out) const
			{
			    saveHeader(out, _format.type());
			    return saveData(out);
			}
    std::istream&	restoreData(std::istream& in)		;
    std::ostream&	saveData(std::ostream& out)	const	;

    ImageFormat::Type	defaultType() const
			{
			    return _format.type();
			}
    void		resize(size_t h, size_t w, const ImageFormat& format)
			{
			    _format = format;
			    w = (_format.depth()*w + 7) / 8;
			    _a.resize(h, w);
			}
			       
  private:
    array2_type		_a;
    ImageFormat		_format;
    Array<BGRA>		_colormap;
};

}
#endif	// !TU_IMAGEPP_H
