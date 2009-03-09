/*
 *  $BJ?@.(B14-19$BG/!JFH!K;:6H5;=QAm9g8&5f=j(B $BCx:n8"=jM-(B
 *  
 *  $BAO:n<T!'?"<G=SIW(B
 *
 *  $BK\%W%m%0%i%`$O!JFH!K;:6H5;=QAm9g8&5f=j$N?&0w$G$"$k?"<G=SIW$,AO:n$7!$(B
 *  $B!JFH!K;:6H5;=QAm9g8&5f=j$,Cx:n8"$r=jM-$9$kHkL)>pJs$G$9!%Cx:n8"=jM-(B
 *  $B<T$K$h$k5v2D$J$7$KK\%W%m%0%i%`$r;HMQ!$J#@=!$2~JQ!$Bh;0<T$X3+<($9$k(B
 *  $BEy$N9T0Y$r6X;_$7$^$9!%(B
 *  
 *  $B$3$N%W%m%0%i%`$K$h$C$F@8$8$k$$$+$J$kB;32$KBP$7$F$b!$Cx:n8"=jM-<T$*(B
 *  $B$h$SAO:n<T$O@UG$$rIi$$$^$;$s!#(B
 *
 *  Copyright 2002-2007.
 *  National Institute of Advanced Industrial Science and Technology (AIST)
 *
 *  Creator: Toshio UESHIBA
 *
 *  [AIST Confidential and all rights reserved.]
 *  This program is confidential. Any using, copying, changing or
 *  giving any information concerning with this program to others
 *  without permission by the copyright holder are strictly prohibited.
 *
 *  [No Warranty.]
 *  The copyright holder or the creator are not responsible for any
 *  damages caused by using this program.
 *  
 *  $Id: Warp.cc,v 1.11 2009-03-09 05:12:32 ueshiba Exp $
 */
#if defined(__INTEL_COMPILER)
#  undef SSE4
#  undef SSSE3
#  undef SSE3
#  undef SSE2

#  define SSE
//#  define SSE2
//#  define SSE3
#endif
#include "TU/Warp.h"

namespace TU
{
#if defined(SSE)
  static inline mmInt16
  mmLinearInterpolate(mmInt16 x, mmInt16 y, mmInt16 d)
  {
      return x + ((d * (y - x)) >> 7);
  }

  template <class T> static inline mmInt16
  mmBilinearInterpolate(const Image<T>& in,
			mmInt16 us, mmInt16 vs, mmInt16 du, mmInt16 dv)
  {
      const mmInt16	ue = us + mmSetAll<mmInt16>(1);
#  if defined(SSE2)
      mmUInt8	uc = mmSet<mmUInt8>(*(int*)&in[mmNth<1>(vs)][mmNth<1>(ue)],
				    *(int*)&in[mmNth<0>(vs)][mmNth<0>(ue)],
				    *(int*)&in[mmNth<1>(vs)][mmNth<1>(us)],
				    *(int*)&in[mmNth<0>(vs)][mmNth<0>(us)]);
#  else
      mmUInt8	uc = mmSet<mmUInt8>(*(int*)&in[mmNth<0>(vs)][mmNth<0>(ue)],
				    *(int*)&in[mmNth<0>(vs)][mmNth<0>(us)]);
#  endif
      mmInt16	ss = mmLinearInterpolate(mmCvt<mmInt16>(uc),
					 mmCvtH<mmInt16>(uc), du);
      vs = vs + mmSetAll<mmInt16>(1);
#  if defined(SSE2)
      uc = mmSet<mmUInt8>(*(int*)&in[mmNth<1>(vs)][mmNth<1>(ue)],
			  *(int*)&in[mmNth<0>(vs)][mmNth<0>(ue)],
			  *(int*)&in[mmNth<1>(vs)][mmNth<1>(us)],
			  *(int*)&in[mmNth<0>(vs)][mmNth<0>(us)]);
#  else
      uc = mmSet<mmUInt8>(*(int*)&in[mmNth<0>(vs)][mmNth<0>(ue)],
			  *(int*)&in[mmNth<0>(vs)][mmNth<0>(us)]);
#  endif
      return mmLinearInterpolate(ss,
				 mmLinearInterpolate(mmCvt<mmInt16>(uc),
						     mmCvtH<mmInt16>(uc), du),
				 dv);
  }
    
  template <> inline mmInt16
  mmBilinearInterpolate(const Image<u_char>& in,
			mmInt16 us, mmInt16 vs, mmInt16 du, mmInt16 dv)
  {
      const mmInt16	ue = us + mmSetAll<mmInt16>(1);
#  if defined(SSE2)
      mmUInt8	uc = mmSet<mmUInt8>(in[mmNth<7>(vs)][mmNth<7>(ue)],
				    in[mmNth<6>(vs)][mmNth<6>(ue)],
				    in[mmNth<5>(vs)][mmNth<5>(ue)],
				    in[mmNth<4>(vs)][mmNth<4>(ue)],
				    in[mmNth<3>(vs)][mmNth<3>(ue)],
				    in[mmNth<2>(vs)][mmNth<2>(ue)],
				    in[mmNth<1>(vs)][mmNth<1>(ue)],
				    in[mmNth<0>(vs)][mmNth<0>(ue)],
				    in[mmNth<7>(vs)][mmNth<7>(us)],
				    in[mmNth<6>(vs)][mmNth<6>(us)],
				    in[mmNth<5>(vs)][mmNth<5>(us)],
				    in[mmNth<4>(vs)][mmNth<4>(us)],
				    in[mmNth<3>(vs)][mmNth<3>(us)],
				    in[mmNth<2>(vs)][mmNth<2>(us)],
				    in[mmNth<1>(vs)][mmNth<1>(us)],
				    in[mmNth<0>(vs)][mmNth<0>(us)]);
#  else
      mmUInt8	uc = mmSet<mmUInt8>(in[mmNth<3>(vs)][mmNth<3>(ue)],
				    in[mmNth<2>(vs)][mmNth<2>(ue)],
				    in[mmNth<1>(vs)][mmNth<1>(ue)],
				    in[mmNth<0>(vs)][mmNth<0>(ue)],
				    in[mmNth<3>(vs)][mmNth<3>(us)],
				    in[mmNth<2>(vs)][mmNth<2>(us)],
				    in[mmNth<1>(vs)][mmNth<1>(us)],
				    in[mmNth<0>(vs)][mmNth<0>(us)]);
#  endif
      mmInt16	ss = mmLinearInterpolate(mmCvt<mmInt16>(uc),
					 mmCvtH<mmInt16>(uc), du);
      vs = vs + mmSetAll<mmInt16>(1);
#  if defined(SSE2)
      uc = mmSet<mmUInt8>(in[mmNth<7>(vs)][mmNth<7>(ue)],
			  in[mmNth<6>(vs)][mmNth<6>(ue)],
			  in[mmNth<5>(vs)][mmNth<5>(ue)],
			  in[mmNth<4>(vs)][mmNth<4>(ue)],
			  in[mmNth<3>(vs)][mmNth<3>(ue)],
			  in[mmNth<2>(vs)][mmNth<2>(ue)],
			  in[mmNth<1>(vs)][mmNth<1>(ue)],
			  in[mmNth<0>(vs)][mmNth<0>(ue)],
			  in[mmNth<7>(vs)][mmNth<7>(us)],
			  in[mmNth<6>(vs)][mmNth<6>(us)],
			  in[mmNth<5>(vs)][mmNth<5>(us)],
			  in[mmNth<4>(vs)][mmNth<4>(us)],
			  in[mmNth<3>(vs)][mmNth<3>(us)],
			  in[mmNth<2>(vs)][mmNth<2>(us)],
			  in[mmNth<1>(vs)][mmNth<1>(us)],
			  in[mmNth<0>(vs)][mmNth<0>(us)]);
#  else
      uc = mmSet<mmUInt8>(in[mmNth<3>(vs)][mmNth<3>(ue)],
			  in[mmNth<2>(vs)][mmNth<2>(ue)],
			  in[mmNth<1>(vs)][mmNth<1>(ue)],
			  in[mmNth<0>(vs)][mmNth<0>(ue)],
			  in[mmNth<3>(vs)][mmNth<3>(us)],
			  in[mmNth<2>(vs)][mmNth<2>(us)],
			  in[mmNth<1>(vs)][mmNth<1>(us)],
			  in[mmNth<0>(vs)][mmNth<0>(us)]);
#  endif
      return mmLinearInterpolate(ss,
				 mmLinearInterpolate(mmCvt<mmInt16>(uc),
						     mmCvtH<mmInt16>(uc), du),
				 dv);
  }
#endif
    
/************************************************************************
*  static functions							*
************************************************************************/
template <class T> static inline T
bilinearInterpolate(const Image<T>& in, int us, int vs, int du, int dv)
{
    T		in00 = in[vs][us],   in01 = in[vs][us+1],
		in10 = in[vs+1][us], in11 = in[vs+1][us+1];
    int		tmp0, tmp1;
    T		out;
    tmp0 = int(in00.r) + ((du * (int(in01.r) - int(in00.r))) >> 7);
    tmp1 = int(in10.r) + ((du * (int(in11.r) - int(in10.r))) >> 7);
    out.r = tmp0 + ((dv * (tmp1 - tmp0)) >> 7);
    tmp0 = int(in00.g) + ((du * (int(in01.g) - int(in00.g))) >> 7);
    tmp1 = int(in10.g) + ((du * (int(in11.g) - int(in10.g))) >> 7);
    out.g = tmp0 + ((dv * (tmp1 - tmp0)) >> 7);
    tmp0 = int(in00.b) + ((du * (int(in01.b) - int(in00.b))) >> 7);
    tmp1 = int(in10.b) + ((du * (int(in11.b) - int(in10.b))) >> 7);
    out.b = tmp0 + ((dv * (tmp1 - tmp0)) >> 7);
    tmp0 = int(in00.a) + ((du * (int(in01.a) - int(in00.a))) >> 7);
    tmp1 = int(in10.a) + ((du * (int(in11.a) - int(in10.a))) >> 7);
    out.a = tmp0 + ((dv * (tmp1 - tmp0)) >> 7);

    return out;
}

template <> inline u_char
bilinearInterpolate(const Image<u_char>& in, int us, int vs, int du, int dv)
{
    int		in00 = in[vs][us],   in01 = in[vs][us+1],
		in10 = in[vs+1][us], in11 = in[vs+1][us+1];
    int		out0 = in00 + ((du * (in01 - in00)) >> 7),
		out1 = in10 + ((du * (in11 - in10)) >> 7);
    
    return out0 + ((dv * (out1 - out0)) >> 7);
}

/************************************************************************
*  class Warp								*
************************************************************************/
//! $B2hA|$r<M1FJQ49$9$k$?$a$N9TNs$r@_Dj$9$k!%(B
/*!
  $BF~NO2hA|E@(Bu$B$O<M1FJQ49(B
  \f[
    \TUbeginarray{c} \TUvec{v}{} \\ 1 \TUendarray \simeq
    \TUvec{H}{} \TUbeginarray{c} \TUvec{u}{} \\ 1 \TUendarray
  \f]
  $B$K$h$C$F=PNO2hA|E@(Bv$B$K<L$5$l$k!%(B
  \param Htinv		$BJQ7A$r;XDj$9$k(B3x3$B<M1FJQ499TNs$N5U9TNs$NE>CV!$$9$J$o$A(B
			\f$\TUtinv{H}{}\f$
  \param inWidth	$BF~NO2hA|$NI}(B
  \param inHeight	$BF~NO2hA|$N9b$5(B
  \param outWidth	$B=PNO2hA|$NI}(B
  \param outWidth	$B=PNO2hA|$N9b$5(B
*/
void
Warp::initialize(const Matrix33d& Htinv,
		 u_int inWidth,  u_int inHeight,
		 u_int outWidth, u_int outHeight)
{
    initialize(Htinv, CameraWithDistortion::Intrinsic(),
	       inWidth, inHeight, outWidth, outHeight);
}

//! $B2hA|$NHs@~7AOD$_$r=|5n$7$?8e$K<M1FJQ49$r9T$&$?$a$N9TNs$H%+%a%iFbIt%Q%i%a!<%?$r@_Dj$9$k!%(B
/*!

  canonical$B:BI8(Bx$B$+$i2hA|:BI8(Bu$B$X$NJQ49$,(B\f$\TUvec{u}{} = {\cal
  K}(\TUvec{x}{})\f$ $B$HI=$5$l$k%+%a%iFbIt%Q%i%a!<%?$K$D$$$F!$$=$N@~7AJQ(B
  $B49ItJ,$rI=$9(B3x3$B>eH>;03Q9TNs$r(BK$B$H$9$k$H!$(B
  \f[
    \TUbeginarray{c} \TUbar{u}{} \\ 1 \TUendarray =
    \TUvec{K}{}
    \TUbeginarray{c} {\cal K}^{-1}(\TUvec{u}{}) \\ 1 \TUendarray
  \f]
  $B$K$h$C$F2hA|$NHs@~7AOD$_$@$1$r=|5n$G$-$k!%K\4X?t$O!$$3$NOD$_$r=|5n$7$?2hA|E@$r(B
  $B<M1FJQ49(BH$B$K$h$C$F=PNO2hA|E@(Bv$B$K<L$9$h$&$KJQ7A%Q%i%a!<%?$r@_Dj$9$k!%$9$J$o$A!$(B
  $BA4BN$NJQ7A$O(B
  \f[
    \TUbeginarray{c} \TUvec{v}{} \\ 1 \TUendarray \simeq
    \TUvec{H}{}\TUvec{K}{}
    \TUbeginarray{c} {\cal K}^{-1}(\TUvec{u}{}) \\ 1 \TUendarray
  \f]
  $B$H$J$k!%(B
  \param Htinv		$BJQ7A$r;XDj$9$k(B3x3$B<M1FJQ499TNs$N5U9TNs$NE>CV(B
  \param Intrinsic	$BF~NO2hA|$K2C$($l$i$l$F$$$kJ|<MOD6J$rI=$9%+%a%iFbIt%Q%i%a!<%?(B
  \param inWidth	$BF~NO2hA|$NI}(B
  \param inHeight	$BF~NO2hA|$N9b$5(B
  \param outWidth	$B=PNO2hA|$NI}(B
  \param outWidth	$B=PNO2hA|$N9b$5(B
*/
void
Warp::initialize(const Matrix33d& Htinv,
		 const CameraBase::Intrinsic& intrinsic,
		 u_int inWidth,  u_int inHeight,
		 u_int outWidth, u_int outHeight)
{
    _fracs.resize(outHeight);
    _width = outWidth;
    
  /* Compute frac for each pixel. */
    const Matrix<double>&	HKtinv = Htinv * intrinsic.Ktinv();
    Vector<double>		leftmost = HKtinv[2];
    for (int v = 0; v < height(); ++v)
    {
	Vector<double>	xc = leftmost;
	FracArray	frac(width());
	int		n = 0;
	for (int u = 0; u < width(); ++u)
	{
	    const Point2d&
		m = intrinsic(Point2d(xc[0]/xc[2], xc[1]/xc[2]));
	    if (0.0 <= m[0] && m[0] <= inWidth - 2 &&
		0.0 <= m[1] && m[1] <= inHeight - 2)
	    {
		if (n == 0)
		    frac.lmost = u;
		frac.us[n] = (short)floor(m[0]);
		frac.vs[n] = (short)floor(m[1]);
		frac.du[n] = (u_char)floor((m[0] - floor(m[0])) * 128.0);
		frac.dv[n] = (u_char)floor((m[1] - floor(m[1])) * 128.0);
		++n;
	    }
	    xc += HKtinv[0];
	}

	_fracs[v].resize(n);
	_fracs[v].lmost = frac.lmost;
	for (int u = 0; u < n; ++u)
	{
	    _fracs[v].us[u] = frac.us[u];
	    _fracs[v].vs[u] = frac.vs[u];
	    _fracs[v].du[u] = frac.du[u];
	    _fracs[v].dv[u] = frac.dv[u];
	}

	leftmost += HKtinv[1];
    }
}

//! $B=PNO2hA|$NHO0O$r;XDj$7$F2hA|$rJQ7A$9$k!%(B
/*!
  \param in	$BF~NO2hA|(B
  \param out	$B=PNO2hA|(B
  \param vs	$BJQ7A7k2L$H$J$kNN0h$N:G=i$N9T$r;XDj$9$k(Bindex
  \param ve	$BJQ7A7k2L$H$J$kNN0h$N:G8e$N9T$N<!$r;XDj$9$k(Bindex$B!%(B0$B$J$i$P=PNO2hA|$N(B
		$B:G8e$N9T$^$GJQ7A7k2L$K$h$C$FKd$a$i$l$k(B
*/
template <class T> void
Warp::operator ()(const Image<T>& in, Image<T>& out, int vs, int ve) const
{
    out.resize(height(), width());
    if (ve == 0)
	ve = out.height();
    
    for (int v = vs; v < ve; ++v)
    {
	const short	*usp  = _fracs[v].us, *vsp  = _fracs[v].vs;
	const u_char	*dup  = _fracs[v].du, *dvp  = _fracs[v].dv;
	T		*outp = out[v] + _fracs[v].lmost;
	T* const	outq  = outp + _fracs[v].width();
#if defined(SSE)
	for (T* const outr = outq - mmUInt8::NElms; outp <= outr; )
	{
	    const u_int	NClrs16 = mmUInt16::NElms/4;
	    mmInt16	uu = mmLoad(usp), vv = mmLoad(vsp);
	    mmUInt8	du = mmLoad(dup), dv = mmLoad(dvp);
	    mmUInt8	du4 = mmQuad0(du), dv4 = mmQuad0(dv);
	    mmStoreU((u_char*)outp,
		     mmCvt<mmUInt8>(
			 mmBilinearInterpolate(
			     in, uu, vv,
			     mmCvt<mmInt16>(du4), mmCvt<mmInt16>(dv4)),
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<NClrs16>(uu),
			     mmShiftElmR<NClrs16>(vv),
			     mmCvtH<mmInt16>(du4), mmCvtH<mmInt16>(dv4))));
	    outp += mmUInt8::NElms/4;

	    du4 = mmQuad1(du);
	    dv4 = mmQuad1(dv);
	    mmStoreU((u_char*)outp,
		     mmCvt<mmUInt8>(
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<2*NClrs16>(uu),
			     mmShiftElmR<2*NClrs16>(vv),
			     mmCvt<mmInt16>(du4), mmCvt<mmInt16>(dv4)),
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<3*NClrs16>(uu),
			     mmShiftElmR<3*NClrs16>(vv),
			     mmCvtH<mmInt16>(du4), mmCvtH<mmInt16>(dv4))));
	    outp += mmUInt8::NElms/4;
	    usp  += mmInt16::NElms;
	    vsp  += mmInt16::NElms;
	    
	    uu  = mmLoad(usp);
	    vv  = mmLoad(vsp);
	    du4 = mmQuad2(du);
	    dv4 = mmQuad2(dv);
	    mmStoreU((u_char*)outp,
		     mmCvt<mmUInt8>(
			 mmBilinearInterpolate(
			     in, uu, vv,
			     mmCvt<mmInt16>(du4), mmCvt<mmInt16>(dv4)),
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<NClrs16>(uu),
			     mmShiftElmR<NClrs16>(vv),
			     mmCvtH<mmInt16>(du4), mmCvtH<mmInt16>(dv4))));
	    outp += mmUInt8::NElms/4;
	    
	    du4 = mmQuad3(du);
	    dv4 = mmQuad3(dv);
	    mmStoreU((u_char*)outp,
		     mmCvt<mmUInt8>(
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<2*NClrs16>(uu),
			     mmShiftElmR<2*NClrs16>(vv),
			     mmCvt<mmInt16>(du4), mmCvt<mmInt16>(dv4)),
			 mmBilinearInterpolate(
			     in,
			     mmShiftElmR<3*NClrs16>(uu),
			     mmShiftElmR<3*NClrs16>(vv),
			     mmCvtH<mmInt16>(du4), mmCvtH<mmInt16>(dv4))));
	    outp += mmUInt8::NElms/4;
	    usp  += mmInt16::NElms;
	    vsp  += mmInt16::NElms;

	    dup  += mmUInt8::NElms;
	    dvp  += mmUInt8::NElms;
	}
#endif
      	while (outp < outq)
	    *outp++ = bilinearInterpolate(in, *usp++, *vsp++, *dup++, *dvp++);
	out[v].setLimits(_fracs[v].lmost, _fracs[v].lmost + _fracs[v].width());
    }
#if defined(SSE)
    mmEmpty();
#endif	
}

//! $B=PNO2hA|$NHO0O$r;XDj$7$F2hA|$rJQ7A$9$k!%(B
/*!
  \param in	$BF~NO2hA|(B
  \param out	$B=PNO2hA|(B
  \param vs	$BJQ7A7k2L$H$J$kNN0h$N:G=i$N9T$r;XDj$9$k(Bindex
  \param ve	$BJQ7A7k2L$H$J$kNN0h$N:G8e$N9T$N<!$r;XDj$9$k(Bindex$B!%(B0$B$J$i$P=PNO2hA|$N(B
		$B:G8e$N9T$^$GJQ7A7k2L$K$h$C$FKd$a$i$l$k(B
*/
template <> void
Warp::operator ()(const Image<u_char>& in, Image<u_char>& out,
		  int vs, int ve) const
{
    out.resize(height(), width());
    if (ve == 0)
	ve = out.height();
    
    for (int v = vs; v < ve; ++v)
    {
	const short	*usp  = _fracs[v].us, *vsp = _fracs[v].vs;
	const u_char	*dup  = _fracs[v].du, *dvp = _fracs[v].dv;
	u_char		*outp = out[v] + _fracs[v].lmost;
	u_char* const	outq  = outp + _fracs[v].width();
#if defined(SSE)
	for (u_char* const outr = outq - mmUInt8::NElms; outp <= outr; )
	{
	    mmUInt8	du = mmLoad(dup), dv = mmLoad(dup);
	    mmInt16	out0 = mmBilinearInterpolate(in,
						     mmLoad(usp), mmLoad(vsp),
						     mmCvt<mmInt16>(du),
						     mmCvt<mmInt16>(dv));
	    usp += mmInt16::NElms;
	    vsp += mmInt16::NElms;
	    mmStoreU(outp, mmCvt<mmUInt8>(out0,
					  mmBilinearInterpolate(in,
					      mmLoad(usp), mmLoad(vsp),
					      mmCvtH<mmInt16>(du),
					      mmCvtH<mmInt16>(dv))));
	    usp  += mmInt16::NElms;
	    vsp  += mmInt16::NElms;
	    dup  += mmUInt8::NElms;
	    dvp  += mmUInt8::NElms;
	    outp += mmUInt8::NElms;
	}
#endif
      	while (outp < outq)
	    *outp++ = bilinearInterpolate(in, *usp++, *vsp++, *dup++, *dvp++);
	out[v].setLimits(_fracs[v].lmost, _fracs[v].lmost + _fracs[v].width());
    }
#if defined(SSE)
    mmEmpty();
#endif	
}

template void	Warp::operator ()(const Image<RGBA>& in,
				  Image<RGBA>& out, int vs, int ve) const;
template void	Warp::operator ()(const Image<ABGR>& in,
				  Image<ABGR>& out, int vs, int ve) const;
}
