/*
 *  $Id: MeasurementMatrix.cc,v 1.4 2002-08-09 00:09:41 ueshiba Exp $
 */
#include "TU/Calib++.h"
#include "TU/Minimize++.h"
#include <iomanip>
#include <stdexcept>

namespace TU
{
/************************************************************************
*  Static functions							*
************************************************************************/
inline u_int	min(u_int a, u_int b)	{return (a < b ? a : b);}
inline int	sign(double x)		{return (x > 0 ? 1 : x < 0 ? -1 : 0);}

static Vector<double>
extpro(const Vector<double>& a, const Vector<double>& b)
{
    const u_int		d = a.dim();
    Vector<double>	v(d * (d+1) / 2);
    int			n = 0;
    for (int i = 0; i < d; ++i)
    {
	v[n++] = a[i]*b[i];
	for (int j = i+1; j < d; ++j)
	    v[n++] = a[i]*b[j] + a[j]*b[i];
    }
    return v;
}

/************************************************************************
*  Global functions							*
************************************************************************/
//! $B4QB,9TNs$r%9%H%j!<%`$K=q$-=P$9!%(B
/*!
  \param out	$B=PNO%9%H%j!<%`!%(B
  \param Wt	$B4QB,9TNs!%(B
  \return	out$B$GM?$($?=PNO%9%H%j!<%`!%(B
*/
std::ostream&
operator <<(std::ostream& out, const MeasurementMatrix& Wt)
{
    using namespace	std;
	
    for (int j = 0; j < Wt.npoints(); ++j)
    {
	for (int i = 0; i < Wt.nframes(); ++i)
	{
	    out << ' ' << Wt[j][3*i  ]/Wt[j][3*i+2]
		<< ' ' << Wt[j][3*i+1]/Wt[j][3*i+2];
	}
	out << endl;
    }
    return out << endl;
}

/************************************************************************
*  class MeasurementMatrix						*
************************************************************************/
//! $B4{B8$N4QB,9TNs$+$i0lIt$N%U%l!<%`$r<h$j=P$7$FItJ,4QB,9TNs$r@8@.$9$k!%(B
/*!
  \param Wt	$B4{B8$N4QB,9TNs(B
  \param index	$B<h$j=P$9%U%l!<%`$N(Bindex$B$r3JG<$7$?G[Ns!%(Bindex[0] = 0,
		index[1] = 1$B$G$"$l$P%U%l!<%`(B0, 1$B$,(BWt$B$+$i<h$j=P$5$l$k!%(B
*/
MeasurementMatrix::MeasurementMatrix(const MeasurementMatrix& Wt,
				     const Array<u_int>&      index)
    :Matrix<double>(Wt.npoints(), min(3*index.dim(), Wt.ncol()))
{
    for (int i = 0; i < nframes(); ++i)
    {
	if (index[i] >= Wt.nframes())
	    throw std::invalid_argument("TU::MeasurementMatrix::MeasurementMatrix: Cannot find specified frame in given measurement matrix!!");
	(*this)(0, 3*i, npoints(), 3) = Wt(0, 3*index[i], Wt.npoints(), 3);
    }
}

//! $BFCD'E@A|$N=E?4$r7W;;$9$k!%(B
/*!
  \return	$B4QB,9TNs$K4^$^$l$kFCD'E@$N=E?4!%3F%U%l!<%`$K$*$1$kFCD'E@$N(B
		$B=E?4$NF1;~:BI8$r%U%l!<%`=g$KJB$Y$?(B3F$B<!85%Y%/%H%k$H$7$FJV$5(B
		$B$l$k!%(B
*/
Vector<double>
MeasurementMatrix::centroid() const
{
    Vector<T>	c(ncol());
    for (int j = 0; j < npoints(); ++j)
	c += (*this)[j];
    return c /= npoints();
}

//! $B;XDj$5$l$?%U%l!<%`$NFCD'E@A|$r<h$j=P$9!%(B
/*!
  \param i	$B%U%l!<%`$N(Bindex$B!%(B
  \return	i$BHVL\$N%U%l!<%`$K4QB,$5$l$?FCD'E@$rI=$9(B npoints() \f$\times\f$
		3$B9TNs!%(B
*/
const Matrix<double>
MeasurementMatrix::frame(u_int i) const
{
    return (*this)(0, 3*i, npoints(), 3);
}

//! $B;XDj$5$l$?%U%l!<%`$NFCD'E@A|$r<h$j=P$9!%(B
/*!
  \param i	$B%U%l!<%`$N(Bindex$B!%(B
  \return	i$BHVL\$N%U%l!<%`$K4QB,$5$l$?FCD'E@$rI=$9(B npoints() \f$\times\f$
		3$B9TNs!%(B
*/
Matrix<double>
MeasurementMatrix::frame(u_int i)
{
    return (*this)(0, 3*i, npoints(), 3);
}

//! 2$BKg$N%U%l!<%`4V$N%"%U%#%s4pAC9TNs$r7W;;$9$k!%(B
/*!
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	frame0, frame1$B4V$N%"%U%#%s4pAC9TNs!%(Bframe0=\f$i\f$,
		frame1=\f$i'\f$$B$N$H$-(B
		\f[
		\forall j~\TUtud{u}{i'j}\TUvec{F}{A}\TUud{u}{ij} \approx 0
		\f]
		$B$rK~$?$9(B\f$\TUvec{F}{A}\in\TUspace{R}{3\times3}\f$
		$B$rJV$9!%$?$@$7!$(B\f$\TUvec{F}{A}\f$$B$N:8>e(B\f$2\times 2\f$
		$BItJ,$O(B0$B!%(B
*/
Matrix<double>
MeasurementMatrix::affineFundamental(u_int frame0, u_int frame1) const
{
    if (frame0 >= nframes() || frame1 >= nframes())
	throw std::invalid_argument("TU::MeasurementMatrix::affineFundamental: Cannot find specified frames!!");
    if (npoints() < 4)
	throw std::invalid_argument("TU::MeasurementMatrix::affineFundamental: At least 4 points needed!!");

    const Vector<T>&	c = centroid();
    Matrix<T>		u(npoints(), 4);
    for (int j = 0; j < npoints(); ++j)
    {
	u[j](0, 2) = (*this)[j](3*frame0, 2) - c(3*frame0, 2);
	u[j](2, 2) = (*this)[j](3*frame1, 2) - c(3*frame1, 2);
    }
    
    Matrix<double>	A(4, 4);
    for (int j = 0; j < npoints(); ++j)
	A += u[j] % u[j];
    Vector<double>	evalue;
    Vector<double>	n = A.eigen(evalue)[3];
#ifdef DEBUG
    std::cerr << "-----------------" << std::endl
	      << "  Eigen values = " << evalue;
#endif
    Matrix<double>	F(3, 3);
    
    F[2][0] = n[0];
    F[2][1] = n[1];
    F[0][2] = n[2];
    F[1][2] = n[3];
    F[2][2] = -(n[0]*c[3*frame0] + n[1]*c[3*frame0+1] +
		n[2]*c[3*frame1] + n[3]*c[3*frame1+1]);

    return F;
}

//! 2$BKg$N%U%l!<%`4V$N4pAC9TNs$r7W;;$9$k!%(B
/*!
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	frame0, frame1$B4V$N4pAC9TNs!%(Bframe0=\f$i\f$, frame1=\f$i'\f$
		$B$N$H$-(B
		\f[
		\forall j~\TUtud{u}{i'j}\TUvec{F}{}\TUud{u}{ij} \approx 0
		\f]
		$B$rK~$?$9(B\f$\TUvec{F}{}\in\TUspace{R}{3\times3}\f$
		$B$rJV$9!%$?$@$7!$(B\f$\det\TUvec{F}{}=0\f$$B!%(B
*/
Matrix<double>
MeasurementMatrix::fundamental(u_int frame0, u_int frame1) const
{
    if (frame0 >= nframes() || frame1 >= nframes())
	throw std::invalid_argument("TU::MeasurementMatrix::fundamental: Cannot find specified frames!!");
    if (npoints() < 8)
	throw std::invalid_argument("TU::MeasurementMatrix::fundamental: At least 8 points needed!!");

#ifdef DEBUG
    std::cerr << "=== BEGIN: TUMeasurementMatrix::fundamental ==="
	      << std::endl;
#endif
  // Normalize the input measurement matrix.
    Normalization	norm0((*this)(0, 3*frame0, npoints(), 3)),
			norm1((*this)(0, 3*frame1, npoints(), 3));
    Matrix<T>		Wt(npoints(), 6);
    Wt(0, 0, npoints(), 3) = (*this)(0, 3*frame0, npoints(), 3) * norm0.Tt();
    Wt(0, 3, npoints(), 3) = (*this)(0, 3*frame1, npoints(), 3) * norm1.Tt();
    
    Matrix<T>	F(3, 3), A(9, 9);
    Vector<T>	weight(npoints());
    T		c = 0.0;
    weight = 1.0 / npoints();
    for (;;)
    {
	Matrix<T>	B(9, 9);
	A = 0.0;
	for (int i = 0; i < npoints(); ++i)
	{
	    Vector<T>		v(9);
	    const Vector<T>	&ldata = Wt[i](0, 3), &rdata = Wt[i](3, 3);
	    
	    v(0, 3) = rdata[0] * ldata;
	    v(3, 3) = rdata[1] * ldata;
	    v(6, 3) = rdata[2] * ldata;
	    A += (weight[i] * v) % v;

	    const Matrix<T>&	tmp = (weight[i] * ldata) % ldata;
	    B(0, 0, 3, 3) += tmp;
	    B(3, 3, 3, 3) += tmp;
	    for (int j = 0; j < 3; ++j)
		for (int k = 0; k < 3; ++k)
		{
		    B[3*j  ][3*k  ] += weight[i] * rdata[j] * rdata[k];
		    B[3*j+1][3*k+1] += weight[i] * rdata[j] * rdata[k];
		}
	}
	A -= c * B;
	Vector<T>	f((T*)F, 9);
	Vector<T>	evalue;
	f = A.eigen(evalue)[8];
#ifdef DEBUG
	std::cerr << " -----------------" << std::endl
		  << "   Weight       = " << weight
		  << "   Eigen values = " << evalue
		  << "   Error        = " << std::sqrt(c) << std::endl;
#endif
	if (fabs(evalue[8]) < 1.0e-12)
	    break;

	c += evalue[8] / (f * B * f);      // Update error.
	
      // Update weights.
      	T	weight_sum = 0.0;
	for (int i = 0; i < npoints(); ++i)
	    weight_sum += (weight[i] =
			   1.0 / ((F(0, 0, 2, 3) * Wt[i](0, 3)).square() +
				  (Wt[i](3, 3) * F(0, 0, 3, 2)).square()));
	weight /= weight_sum;
    }
#ifdef DEBUG
    SVDecomposition<T>	svdF(F);
    std::cerr << " Singular values of uncorrected F = " << svdF.diagonal();
#endif
    const Matrix<T>&	Ainv = A.pinv(1.0e6);
    Vector<T>		f((T*)F, 9);
    
    for (;;)
    {
	const Matrix<T>&	G = F.adj().trns();
	Vector<T>		g((T*)G, 9);
	const T			detF = f * g / 3.0;
#ifdef DEBUG	
	std::cerr << " det(F) = " << detF << std::endl;
#endif
	if (fabs(detF) < 1.0e-12)
	    break;
	
	f -= detF / (g * Ainv * g) * (Ainv * g);
	f.normalize();
    }
#ifdef DEBUG    
    SVDecomposition<T>	svdF2(F);
    std::cerr << " Singular values of corrected F = " << svdF2.diagonal();
    std::cerr << "=== END:   MeasurementMatrix::fundamental ==="
	      << std::endl;
#endif
    return norm1.Tt() * F * norm0.T();
}

//! 2$BKg$N%U%l!<%`4V$N<M1FJQ499TNs$r7W;;$9$k!%(B
/*!
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	frame0, frame1$B4V$N<M1FJQ499TNs!%(Bframe0=\f$i\f$,
		frame1=\f$i'\f$$B$N$H$-(B
		\f[
		\forall j~\TUtud{u}{i'j} \approx \TUvec{H}{}\TUud{u}{ij}
		\f]
		$B$rK~$?$9(B\f$\TUvec{H}{}\in\TUspace{R}{3\times3}\f$$B$rJV$9!%(B
*/
Matrix<double>
MeasurementMatrix::homography(u_int frame0, u_int frame1) const
{
    if (frame0 >= nframes() || frame1 >= nframes())
	throw std::invalid_argument("TU::MeasurementMatrix::homography: Cannot find specified frames!!");
    if (npoints() < 4)
	throw std::invalid_argument("TU::MeasurementMatrix::homography: At least 4 points needed!!");

  // Normalize the input measurement matrix.
    Normalization	norm0((*this)(0, 3*frame0, npoints(), 3)),
			norm1((*this)(0, 3*frame1, npoints(), 3));
    Matrix<T>		Wt(npoints(), 6);
    Wt(0, 0, npoints(), 3) = (*this)(0, 3*frame0, npoints(), 3) * norm0.Tt();
    Wt(0, 3, npoints(), 3) = (*this)(0, 3*frame1, npoints(), 3) * norm1.Tt();

    Matrix<T>	Q(9, 9);
    for (int j = 0; j < npoints(); ++j)
    {
	const Vector<T>	&ldata = Wt[j](0, 3), &rdata = Wt[j](3, 3);
	const Matrix<double>	X = ldata % ldata;
	
	Q(0, 0, 3, 3) += (rdata[1]*rdata[1] + rdata[2]*rdata[2]) * X;
	Q(3, 3, 3, 3) += (rdata[2]*rdata[2] + rdata[0]*rdata[0]) * X;
	Q(6, 6, 3, 3) += (rdata[0]*rdata[0] + rdata[1]*rdata[1]) * X;
	Q(0, 3, 3, 3) = (Q(3, 0, 3, 3) -= rdata[0]*rdata[1] * X);
	Q(0, 6, 3, 3) = (Q(6, 0, 3, 3) -= rdata[2]*rdata[0] * X);
	Q(3, 6, 3, 3) = (Q(6, 3, 3, 3) -= rdata[1]*rdata[2] * X);
    }

    Matrix<double>	H(3, 3);
    Vector<double>	h((double*)H, 9);
    Vector<double>	evalue(9);
    h = Q.eigen(evalue)[8];
#ifdef DEBUG
    std::cerr << "TU::MeasurementMatrix::homography():: eigen values = "
	      << evalue;
#endif
    H = norm1.Tinv() * H * norm0.T();
    double	det = H.det();
    if (det > 0)
	H /= pow(det, 1.0/3.0);
    else
	H /= -pow(-det, 1.0/3.0);
    
    return H;
}

//! 2$BKg$N%U%l!<%`4V$N2sE>9TNs$r7W;;$9$k!%(B
/*!
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	frame0, frame1$B4V$N2sE>9TNs!%(Bframe0=\f$i\f$,
		frame1=\f$i'\f$$B$N$H$-(B
		\f[
		\forall j~\TUtud{u}{i'j} \approx \TUtvec{R}{}\TUud{u}{ij}
		\f]
		$B$rK~$?$9(B\f$\TUtvec{R}{}\in\TUspace{R}{3\times3}\f$$B$rJV$9!%(B
		$B$?$@$7!$(B\f$\TUtvec{R}{}\TUvec{R}{} = \TUvec{I}{3}\f$$B!%(B
*/
Matrix<double>
MeasurementMatrix::rotation(u_int frame0, u_int frame1) const
{
    if (frame0 >= nframes() || frame1 >= nframes())
	throw std::invalid_argument("TU::MeasurementMatrix::rotation: Cannot find specified frames!!");
    if (npoints() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::rotation: At least 3 points needed!!");

    Matrix<double>    M(3, 3);                        // Moment matrix.
    for (int j = 0; j < npoints(); ++j)
    {
	const Vector<T>	&ldata = (*this)[j](3*frame0, 3),
				&rdata = (*this)[j](3*frame1, 3);
	M += ldata.normal() % rdata.normal();
    }
    SVDecomposition<double>   svd(M);
    
    return svd.Ut().trns() * svd.Vt();
}

//! $BJ#?t$N;kE@$+$iJ?LL%Q%?!<%s$r4QB,$7$F%+%a%i%-%c%j%V%l!<%7%g%s(B($BOD$_?dDj$J$7(B)$B$r9T$&!%(B
/*!
  $B%o!<%k%I:BI87O$r;2>HJ?LL$K8GDj$7!$J?LL>e$NFCD'E@$NF1<!:BI8$r(B
  \f$[X_j, Y_j, 0, 1]^\top~(j=0,\ldots,P-1)\f$$B$H$9$k!%$3$l$rFbIt%Q%i%a!<%?(B
  $B$r8GDj$7$F;kE@(Bi$B$+$i4QB,$7$?A|$r(B\f$\TUud{u}{ij}=[u_{ij},v_{ij},1]^\top~
  (i=0,\ldots,F-1)\f$$B$H$9$k!%$3$l$i$N%G!<%?$h$j4QB,9TNs$r(B
  \f[
    \TUtvec{W}{} =
    \TUbeginarray{cccccc}
      X_0     & Y_0     & 1 & \TUtud{u}{00}    & \cdots & \TUtud{u}{0~F-1} \\
      \vdots  & \vdots  & \vdots & \vdots      &	& \vdots	   \\
      X_{P-1} & Y_{P-1} & 1 & \TUtud{u}{P-1~0} & \cdots & \TUtud{u}{P-1~F-1}
    \TUendarray
  \f]
  $B$H9=@.$9$k!%$3$N4QB,9TNs$KBP$7$FK\%"%k%4%j%:%`E,MQ$9$k$H!$%+%a%i$NFbIt(B
  $B%Q%i%a!<%?$*$h$S3F;kE@$G$N%+%a%i$N30It%Q%i%a!<%?$,7W;;$5$l$k!%(B
  \param cameras $B3F;kE@$K$D$$$F!$;2>HJ?LL$K8GDj$5$l$?%o!<%k%I:BI87O$+$i8+$?(B
		 $B%+%a%i$N30It%Q%i%a!<%?$,JV$5$l$k!%(B
  \return	 $B%+%a%i$NFbIt%Q%i%a!<%?(B($BOD$_?dDjCM$r4^$^$J$$(B)$B!%(B
 */
Camera::Intrinsic
MeasurementMatrix::calibrateWithPlanes(Array<CanonicalCamera>& cameras)
    const
{
    Camera::Intrinsic	K(initializeCalibrationWithPlanes(cameras));
    refineCalibrationWithPlanes(K, cameras);

    return K;
}

//! $BJ#?t$N;kE@$+$iJ?LL%Q%?!<%s$r4QB,$7$F%+%a%i%-%c%j%V%l!<%7%g%s(B($BOD$_?dDj$"$j(B)$B$r9T$&!%(B
/*!
  $B%o!<%k%I:BI87O$r;2>HJ?LL$K8GDj$7!$J?LL>e$NFCD'E@$NF1<!:BI8$r(B
  \f$[X_j, Y_j, 0, 1]^\top~(j=0,\ldots,P-1)\f$$B$H$9$k!%$3$l$rFbIt%Q%i%a!<%?(B
  $B$r8GDj$7$F;kE@(Bi$B$+$i4QB,$7$?A|$r(B\f$\TUud{u}{ij}=[u_{ij},v_{ij},1]^\top~
  (i=0,\ldots,F-1)\f$$B$H$9$k!%$3$l$i$N%G!<%?$h$j4QB,9TNs$r(B
  \f[
    \TUtvec{W}{} =
    \TUbeginarray{cccccc}
      X_0     & Y_0     & 1 & \TUtud{u}{00}    & \cdots & \TUtud{u}{0~F-1} \\
      \vdots  & \vdots  & \vdots & \vdots      &	& \vdots	   \\
      X_{P-1} & Y_{P-1} & 1 & \TUtud{u}{P-1~0} & \cdots & \TUtud{u}{P-1~F-1}
    \TUendarray
  \f]
  $B$H9=@.$9$k!%$3$N4QB,9TNs$KBP$7$FK\%"%k%4%j%:%`E,MQ$9$k$H!$%+%a%i$NFbIt(B
  $B%Q%i%a!<%?$*$h$S3F;kE@$G$N%+%a%i$N30It%Q%i%a!<%?$,7W;;$5$l$k!%(B
  \param cameras $B3F;kE@$K$D$$$F!$;2>HJ?LL$K8GDj$5$l$?%o!<%k%I:BI87O$+$i8+$?(B
		 $B%+%a%i$N30It%Q%i%a!<%?$,JV$5$l$k!%(B
  \return	 $B%+%a%i$NFbIt%Q%i%a!<%?(B($BOD$_?dDjCM$r4^$`(B).
 */
CameraWithDistortion::Intrinsic
MeasurementMatrix::calibrateWithPlanesConsideringDistortion
	(Array<CanonicalCamera>& cameras) const
{
    CameraWithDistortion::Intrinsic
	K(initializeCalibrationWithPlanes(cameras));
    refineCalibrationWithPlanes(K, cameras);

    return K;
}

//! $B4QB,9TNs$r%"%U%#%s%+%a%i9TNs$HFCD'E@$N(B3$B<!85%"%U%#%s:BI8$KJ,2r$9$k!%(B
/*!
  \param P	$B$9$Y$F$N;kE@$K$*$1$k%"%U%#%s%+%a%i9TNs$r!$(B
		\f[
		  \TUvec{P}{} =
		  \TUbeginarray{c} \TUvec{P}{0} \\ \vdots \\ \TUvec{P}{F-1}
		  \TUendarray
		\f]
		(\f$\TUvec{P}{i}\in\TUspace{R}{3\times 4}\f$$B$O;kE@(Bi$B$K$*$1$k(B
		$B%"%U%#%s%+%a%i9TNs(B)$B$J$k(B\f$3F\times 4\f$$B9TNs$H$7$FJV$9!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N(B3$B<!85%"%U%#%s6u4V$G$N0LCV$r!$$=$NF1<!:BI8$r(B
		$BBh(Bj$B9T$H$9$k(B\f$P\times 4\f$$B9TNs$H$7$FJV$9!%(B
*/
void
MeasurementMatrix::affineFactorization(Matrix<T>& P, Matrix<T>& Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::affineFactorization: Two or more frames required!!");
    
  // Make affine measurement matrix around the centroid.
    const Vector<T>	c = centroid();
    Matrix<T>		Wa(2*nframes(), npoints());
    for (int i = 0; i < nframes(); ++i)
	for (int j = 0; j < npoints(); ++j)
	{
	    Wa[2*i  ][j] = (*this)[j][3*i  ] - c[3*i  ];
	    Wa[2*i+1][j] = (*this)[j][3*i+1] - c[3*i+1];
	}

  // Factor the affine measurement matrix.
    SVDecomposition<T> svd(Wa);
#ifdef DEBUG
    std::cerr << "  singular vaules of Wa: " << svd.diagonal();
#endif

  // Compute affine motion and shape.
    P.resize(3*nframes(), 4);
    Xt.resize( npoints(), 4);
    for (int i = 0; i < nframes(); ++i)
    {
	for (int n = 0; n < 3; ++n)
	{
	    P[3*i  ][n] = svd.Vt()[n][2*i];
	    P[3*i+1][n] = svd.Vt()[n][2*i+1];
	    P[3*i+2][n] = 0.0;
	}
	P[3*i  ][3] = c[3*i  ];
	P[3*i+1][3] = c[3*i+1];
	P[3*i+2][3] = 1.0;
    }
    for (int j = 0; j < npoints(); ++j)
    {
	for (int n = 0; n < 3; ++n)
	    Xt[j][n] = svd.Ut()[n][j] * svd[n];
	Xt[j][3] = 1.0;
    }
}

//! $B%+%a%i9TNs$HFCD'E@$N(B3$B<!85:BI8$r%"%U%#%s:BI87O$+$i%f!<%/%j%C%I:BI87O$KD>$9!%(B
/*!
  $B0z?t$H$J$k%+%a%i9TNs$HFCD'E@:BI8$O(BaffineFactorization()$B$K$h$C$F7W;;$G$-$k!%(B
  $B$?$@$7!$$=$N$H$-$N4QB,9TNs$K4^$^$l$k2hA|>e$G$NFCD'E@A|$N0LCV$O!$FbIt%Q%i(B
  $B%a!<%?$N1F6A$r%-%c%s%;%k$7$?(Bcanonical$B:BI8$GI=8=$5$l$F$$$J$1$l$P$J$i$J$$!%(B
  \param P	$B$9$Y$F$N;kE@$K$*$1$k%"%U%#%s%+%a%i9TNs$r(B\f$3F\times 4\f$$B9TNs(B
		$B$H$7$FM?$($k!%$3$l$,%f!<%/%j%C%I:BI87O$KJQ49$5$l$FJV$5$l$k!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N%"%U%#%s:BI8$r(B\f$P\times 4\f$$B9TNs$H$7$FM?$($k!%(B
		$B$3$l$,%f!<%/%j%C%I:BI8$KJQ49$5$l$FJV$5$l$k!%(B
*/
void
MeasurementMatrix::affineToMetric(Matrix<T>& P, Matrix<T>& Xt)
{
    const u_int	nframes = P.nrow() / 3;

    if (nframes < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::affineToMetric: Three or more frames required!!");
    
  // Construct constraint matrix C.
    Matrix<T>	C(2*nframes, 6);
    for (int i = 0; i < nframes; ++i)
    {
	const Vector<T>	&p = P[3*i](0, 3), &q = P[3*i+1](0, 3);
	const T			one_uu = 1.0 + P[3*i  ][3]*P[3*i  ][3],
				one_vv = 1.0 + P[3*i+1][3]*P[3*i+1][3],
				uv     = P[3*i][3]*P[3*i+1][3];
	
	C[2*i  ][0] =	 p[0]*p[0]/one_uu - q[0]*q[0]/one_vv;
	C[2*i  ][1] =	 p[1]*p[1]/one_uu - q[1]*q[1]/one_vv;
	C[2*i  ][2] =	 p[2]*p[2]/one_uu - q[2]*q[2]/one_vv;
	C[2*i  ][3] = 2*(p[1]*p[2]/one_uu - q[1]*q[2]/one_vv);
	C[2*i  ][4] = 2*(p[2]*p[0]/one_uu - q[2]*q[0]/one_vv);
	C[2*i  ][5] = 2*(p[0]*p[1]/one_uu - q[0]*q[1]/one_vv);
	C[2*i+1][0] = p[0]*q[0] - uv*(p[0]*p[0]/one_uu + q[0]*q[0]/one_vv)/2;
	C[2*i+1][1] = p[1]*q[1] - uv*(p[1]*p[1]/one_uu + q[1]*q[1]/one_vv)/2;
	C[2*i+1][2] = p[2]*q[2] - uv*(p[2]*p[2]/one_uu + q[2]*q[2]/one_vv)/2;
	C[2*i+1][3] = p[1]*q[2] + p[2]*q[1]
		    - uv*(p[1]*p[2]/one_uu + q[1]*q[2]/one_vv);
	C[2*i+1][4] = p[2]*q[0] + p[0]*q[2]
		    - uv*(p[2]*p[0]/one_uu + q[2]*q[0]/one_vv);
	C[2*i+1][5] = p[0]*q[1] + p[1]*q[0]
		    - uv*(p[0]*p[1]/one_uu + q[0]*q[1]/one_vv);
    }

  // Construct B.
    Vector<T>		evalue;
    const Vector<T>&	b = (C.trns()*C).eigen(evalue)[5];
#ifdef DEBUG
    std::cerr << "  eigen values of CC: " << evalue;
#endif
    Matrix<T>	B(3, 3);
    B[0][0] = b[0];
    B[1][1] = b[1];
    B[2][2] = b[2];
    B[1][2] = B[2][1] = b[3];
    B[2][0] = B[0][2] = b[4];
    B[0][1] = B[1][0] = b[5];
    if (B.det() < 0.0)
	B *= -1.0;
    
  // Factor B into A and At s.t. B = A * At.
    Matrix<T>	At(4, 4);
    At(0, 0, 3, 3) = B.cholesky();
    At[3][3] = 1.0;

  // Transform from affine to Euclidean world.
    P  *= At.trns();
    Xt *= At.inv();

  // Convert Affine projection matrices to Euclidean ones.
    for (int i = 0; i < nframes; ++i)
    {
	const T		xc = P[3*i][3], yc = P[3*i+1][3];
	Matrix<T>	Nt(3, 3);
	Nt[0] = P[3*i  ](0, 3);
	Nt[1] = P[3*i+1](0, 3);
	Nt[2] = -xc * Nt[0] - yc * Nt[1];

	SVDecomposition<T>	svd(Nt);
	P(3*i, 0, 3, 3)  = svd.Vt().trns() * svd.Ut();
	P(3*i, 3, 3, 1) *= ((2.0 + xc*xc + yc*yc) / (svd[0] + svd[1]));
    }

  /* Fix WC to 0th camera
     and set distance between 0th and 1st camera to unity. */
    At(0, 0, 3, 3) = P(0, 0, 3, 3);	// rotation of 0th camera.
    At[3][0] = -P[0][3];
    At[3][1] = -P[1][3];
    At[3][2] = -P[2][3];
    At[3](0, 3) *= P(0, 0, 3, 3);	// t0: translation of 0th camera.
    Vector<T>	t(3);
    t[0] = -P[3][3];
    t[1] = -P[4][3];
    t[2] = -P[5][3];
    (t *= P(3, 0, 3, 3)) -= At[3](0, 3);	// t1 - t0.
    At[3] /= t.length();
    
    P  *= At.trns();
    Xt *= At.inv();
}

//! $B4QB,9TNs$rF);kEj1F%+%a%i9TNs$HFCD'E@$N(B3$B<!85<M1F:BI8$KJ,2r$9$k!%(B
/*!
  \param P	$B$9$Y$F$N;kE@$K$*$1$kF);kEj1F%+%a%i9TNs$r!$(B
		\f[
		  \TUvec{P}{} =
		  \TUbeginarray{c} \TUvec{P}{0} \\ \vdots \\ \TUvec{P}{F-1}
		  \TUendarray
		\f]
		(\f$\TUvec{P}{i}\in\TUspace{R}{3\times 4}\f$$B$O;kE@(Bi$B$K$*$1$k(B
		$BF);kEj1F%+%a%i9TNs(B)$B$J$k(B\f$3F\times 4\f$$B9TNs$H$7$FJV$9!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N(B3$B<!85<M1F6u4V$G$N0LCV$r!$$=$NF1<!:BI8$rBh(Bj$B9T(B
		$B$H$9$k(B\f$P\times 4\f$$B9TNs$H$7$FJV$9!%(B
*/
void
MeasurementMatrix::projectiveFactorization(Matrix<T>& P, Matrix<T>& Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::projectiveFactorization: Two or more frames required!!");
    
  // Normalize the input measurement matrix.
    Array<Normalization>	norm(nframes());
    Matrix<T>			Wt(npoints(), 3*nframes());
    for (int i = 0; i < nframes(); ++i)
    {
	norm[i].initialize((*this)(0, 3*i, npoints(), 3));
	Wt(0, 3*i, npoints(), 3) = (*this)(0, 3*i, npoints(), 3)
				 * norm[i].Tt();
    }

  // Minimize cost function.
    MeasurementMatrix::CostPF	f(Wt);
    Vector<T>	mu((f.nframes() - 1) * (f.npoints() - 1));
    mu = 1.0;
    f.minimize(mu);

  // Decompose measurement matrix into motion(P) and shape(Xt).
    P.resize(3*nframes(), 4);
    Xt.resize( npoints(), 4);
    for (int n = 0; n < 4; ++n)
    {
	for (int i = 0; i < 3*nframes(); ++i)
	    P[i][n] = f.Vt()[n][i];
	for (int j = 0; j < npoints(); ++j)
	    Xt[j][n] = f.s()[n] * f.Ut()[n][j];
    }

  // Unnormalize camera matrices and recover image origins.
    for (int i = 0; i < nframes(); ++i)
	P(3*i, 0, 3, 4) = norm[i].Tinv() * P(3*i, 0, 3, 4);
}

//! $B%+%a%i9TNs$HFCD'E@$N(B3$B<!85:BI8$r<M1F:BI87O$+$i%f!<%/%j%C%I:BI87O$KD>$9!%(B
/*!
  $B0z?t$H$J$k%+%a%i9TNs$HFCD'E@:BI8$O(BprojectiveFactorization()$B$K$h$C$F7W;;$G(B
  $B$-$k!%$?$@$7!$$=$N$H$-$N4QB,9TNs$K4^$^$l$k2hA|>e$G$NFCD'E@A|$N0LCV$O!$FbIt(B
  $B%Q%i%a!<%?$N1F6A$r%-%c%s%;%k$7$?(Bcanonical$B:BI8$GI=8=$5$l$F$$$J$1$l$P$J$i$J$$!%(B
  \param P	$B$9$Y$F$N;kE@$K$*$1$kF);kEj1F%+%a%i9TNs$r(B\f$3F\times 4\f$$B9TNs(B
		$B$H$7$FM?$($k!%$3$l$,%f!<%/%j%C%I:BI87O$KJQ49$5$l$FJV$5$l$k!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N<M1F:BI8$r(B\f$P\times 4\f$$B9TNs$H$7$FM?$($k!%(B
		$B$3$l$,%f!<%/%j%C%I:BI8$KJQ49$5$l$FJV$5$l$k!%(B
*/
void
MeasurementMatrix::projectiveToMetric(Matrix<T>& P, Matrix<T>& Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::projectiveToMetric: Two or more frames required!!");
    
  // Compute the focal point and the inverse projection of 0-th camera.
    SVDecomposition<T>	svd(P(0, 0, 3, 4));
    const Vector<T>&		c0 = svd.Ut()[3];
    const Matrix<T>&		P0tinv = svd.Ut()[0] % svd.Vt()[0] / svd[0]
				       + svd.Ut()[1] % svd.Vt()[1] / svd[1]
				       + svd.Ut()[2] % svd.Vt()[2] / svd[2];

  // Initial guess of the parameter: q.
    Matrix<T>	A(5*(nframes()-1), 4);
    Vector<T>	b(5*(nframes()-1));
    for (int i = 0; i < nframes()-1; ++i)
    {
	const Matrix<T>&	Pi = P(3*(i+1), 0, 3, 4);
	const Matrix<T>&	Nt = Pi * P0tinv;
	const Vector<T>&	e  = Pi * c0;
	
	A[5*i  ](0, 3)	= 2*e[0]*Nt[0] - 2*e[2]*Nt[2];
	A[5*i+1](0, 3)	= 2*e[1]*Nt[1] - 2*e[2]*Nt[2];
	A[5*i+2](0, 3)	= e[1]*Nt[0] + e[0]*Nt[1];
	A[5*i+3](0, 3)	= e[2]*Nt[0] + e[0]*Nt[2];
	A[5*i+4](0, 3)	= e[2]*Nt[1] + e[1]*Nt[2];
	A[5*i  ][3]	= e[0]*e[0] - e[2]*e[2];
	A[5*i+1][3]	= e[1]*e[1] - e[2]*e[2];
	A[5*i+2][3]	= e[0]*e[1];
	A[5*i+3][3]	= e[0]*e[2];
	A[5*i+4][3]	= e[1]*e[2];
	b[5*i  ]	= -Nt[0]*Nt[0] + Nt[2]*Nt[2];
	b[5*i+1]	= -Nt[1]*Nt[1] + Nt[2]*Nt[2];
	b[5*i+2]	= -Nt[0]*Nt[1];
	b[5*i+3]	= -Nt[0]*Nt[2];
	b[5*i+4]	= -Nt[1]*Nt[2];
    }
    Vector<T>	q = b * A;
    q.solve(A.trns() * A);
#ifdef DEBUG	// Check the order of magnitute of last two eigen values.
    std::cerr << "q = " << q;
    std::cerr << "sqrt(fabs(q(0, 3)*q(0, 3) - q[3])) = "
	      << sqrt(fabs(q(0, 3)*q(0, 3) - q[3])) << std::endl;
#endif
    Vector<T>				p = q(0, 3);
    CostPM				f(A, b);
    NullConstraint<T, CostPM::AT>	g;
    tuMinimizeSquare(f, g, p);
    
  // Construct projective transformation H.
    Matrix<T>	Ht(4, 4);
    Ht(0, 0, 3, 4) = P0tinv.trns() + p % c0;
    Ht[3] = c0;
    
  // Transform from projective to Euclidean world.
    P  *= Ht.trns();
    Xt *= Ht.inv();

  // Determine sign of Ht[3] so that points are in front of the cameras.
    int	sig = 0;
    for (int j = 0; j < npoints(); ++j)
	sig += (Xt[j][3] > 0.0 ? 1 : Xt[j][3] < 0.0 ? -1 : 0);
    if (sig < 0)
    {
	P(0, 3, P.nrow(), 1) *= -1.0;
	Xt(0, 3, Xt.nrow(), 1) *= -1.0;
    }
	
  // Normalize so that 4-th component of Xt is 1.
    for (int j = 0; j < npoints(); ++j)
	Xt[j] /= Xt[j][3];

  // Extract translation and rotation components of camera motion.
    Array<CanonicalCamera>	cameras(nframes());
    for (int i = 0; i < nframes(); ++i)
	cameras[i].setProjection(P(3*i, 0, 3, 4));

    bundleAdjustment(cameras, Xt);

  // Set refined camera parameters.
    for (int i = 0; i < nframes(); ++i)
	P(3*i, 0, 3, 4) = cameras[i].P();
}

//! $B8D!9$N%+%a%i$N>GE@5wN%$r?dDj$7$D$D!$%+%a%i9TNs$HFCD'E@$N(B3$B<!85:BI8$r<M1F:BI87O$+$i%f!<%/%j%C%I:BI87O$KD>$9!%(B
/*!
  $B0z?t$H$J$k%+%a%i9TNs$HFCD'E@:BI8$O(BprojectiveFactorization()$B$K$h$C$F7W;;$G(B
  $B$-$k!%$?$@$7!$$=$N$H$-$N4QB,9TNs$K4^$^$l$k2hA|>e$G$NFCD'E@A|$N0LCV$O!$>GE@(B
  $B5wN%0J30$NFbIt%Q%i%a!<%?$N1F6A$r%-%c%s%;%k$7$?:BI8$GI=8=$5$l$F$$$J$1$l$P$J(B
  $B$i$J$$!%>GE@5wN%$O8D!9$N%+%a%iKh$K0[$J$C$F$$$F$h$$!%(B
  \param P	$B$9$Y$F$N;kE@$K$*$1$kF);kEj1F%+%a%i9TNs$r(B\f$3F\times 4\f$$B9TNs(B
		$B$H$7$FM?$($k!%$3$l$,%f!<%/%j%C%I:BI87O$KJQ49$5$l$FJV$5$l$k!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N<M1F:BI8$r(B\f$P\times 4\f$$B9TNs$H$7$FM?$($k!%(B
		$B$3$l$,%f!<%/%j%C%I:BI8$KJQ49$5$l$FJV$5$l$k!%(B
*/
void
MeasurementMatrix::projectiveToMetricWithFocalLengthsEstimation
			(Matrix<T>& P, Matrix<T>& Xt) const
{
    initializeFocalLengthsEstimation(P, Xt);
    
  // Extract translation and rotation components of camera motion.
    Array<CameraWithFocalLength>	cameras(nframes());
    for (int i = 0; i < nframes(); ++i)
	cameras[i].setProjection(P(3*i, 0, 3, 4));
    
    bundleAdjustment(cameras, Xt);

  // Set refined camera parameters.
    for (int i = 0; i < nframes(); ++i)
	P(3*i, 0, 3, 4) = cameras[i].P();
}

//! $B$9$Y$F$N%+%a%i$K6&DL$N>GE@5wN%$r?dDj$7$D$D!$%+%a%i9TNs$HFCD'E@$N(B3$B<!85:BI8$r<M1F:BI87O$+$i%f!<%/%j%C%I:BI87O$KD>$9!%(B
/*!
  $B0z?t$H$J$k%+%a%i9TNs$HFCD'E@:BI8$O(BprojectiveFactorization()$B$K$h$C$F7W;;$G(B
  $B$-$k!%$?$@$7!$$=$N$H$-$N4QB,9TNs$K4^$^$l$k2hA|>e$G$NFCD'E@A|$N0LCV$O!$>GE@(B
  $B5wN%0J30$NFbIt%Q%i%a!<%?$N1F6A$r%-%c%s%;%k$7$?:BI8$GI=8=$5$l$F$$$J$1$l$P$J(B
  $B$i$J$$!%>GE@5wN%$O$9$Y$F$N%+%a%i$G6&DL$G$J$1$l$P$J$i$J$$!%(B
  \param P	$B$9$Y$F$N;kE@$K$*$1$kF);kEj1F%+%a%i9TNs$r(B\f$3F\times 4\f$$B9TNs(B
		$B$H$7$FM?$($k!%$3$l$,%f!<%/%j%C%I:BI87O$KJQ49$5$l$FJV$5$l$k!%(B
  \param Xt	$B$9$Y$F$NFCD'E@$N<M1F:BI8$r(B\f$P\times 4\f$$B9TNs$H$7$FM?$($k!%(B
		$B$3$l$,%f!<%/%j%C%I:BI8$KJQ49$5$l$FJV$5$l$k!%(B
*/
void
MeasurementMatrix::projectiveToMetricWithCommonFocalLengthEstimation
			(Matrix<T>& P, Matrix<T>& Xt) const
{
    initializeFocalLengthsEstimation(P, Xt);
    
  // Extract translation and rotation components of camera motion.
    Array<CanonicalCamera>	cameras(nframes());
    T				k = 0.0;
    for (int i = 0; i < nframes(); ++i)
    {
	Matrix<T>	Nt = P(3*i, 0, 3, 3);
	Nt[0] *= 0.5;
	Nt[1] *= 0.5;
	SVDecomposition<T>	svd(Nt);
	if (svd.diagonal()[2] < 0.0)
	    throw std::domain_error("TU::MeasurementMatrix::projectiveToMetricWithCommonFocalLengthEstimation: cannot extract camera rotation due to nevative singular value!");
	const Matrix<T>&	Rt = svd.Vt().trns() * svd.Ut();
	T			p = Rt[0]*Nt[0] + Rt[1]*Nt[1], q = Rt[2]*Nt[2];
	Vector<T>	t(3);
	t[0] = -P[3*i  ][3] / p;
	t[1] = -P[3*i+1][3] / p;
	t[2] = -P[3*i+2][3] / q;
	t *= Rt;
	cameras[i].setTranslation(t).setRotation(Rt);
	k += p / q;
    }

    CameraWithFocalLength::Intrinsic	K(k / nframes());
    bundleAdjustment(K, cameras, Xt);

  // Set refined camera parameters.
    for (int i = 0; i < nframes(); ++i)
    {
	P(3*i, 0, 3, 4) = cameras[i].P();
	P[3*i  ] *= K.k();
	P[3*i+1] *= K.k();
    }
}

//! $BJ?LL>e$N%Q%?!<%s$rJ#?t$N;kE@$+$i4QB,$7$?2hA|$NA|$rMQ$$$F!$%+%a%i$NFbIt%Q%i%a!<%?(B($BOD$_$J$7(B)$B$H30It%Q%i%a!<%?$r(Brefine$B$9$k!%(B
/*!
  $B4QB,9TNs$N7A<0$K$D$$$F$O(BcalibrateWithPlanes()$B$r;2>H!%(B
  \param K	 $B%+%a%i$NFbIt%Q%i%a!<%?$N=i4|CM(B($BOD$_?dDjCM$r4^$^$J$$(B)$B$rM?$($k!%(B
		 $B$=$N:GE,?dDjCM$,JV$5$l$k!%(B
  \param cameras $B3F;kE@$K$D$$$F!$;2>HJ?LL$K8GDj$5$l$?%o!<%k%I:BI87O$+$i8+$?(B
		 $B%+%a%i$N30It%Q%i%a!<%?$N=i4|CM$rM?$($k!%$=$N:GE,?dDjCM$,(B
		 $BJV$5$l$k!%(B
 */
void
MeasurementMatrix::refineCalibrationWithPlanes
    (Camera::Intrinsic& K, Array<CanonicalCamera>& cameras) const
{
    if (nframes() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::refineCalibrationWithPlanes: Two or more frames required!!");

    CostCP				err(*this, 5);
    NullConstraint<T, CostCP::ATA>	g;

    tuMinimizeSquareSparse(err, g, K, cameras, 30);
}

//! $BJ?LL>e$N%Q%?!<%s$rJ#?t$N;kE@$+$i4QB,$7$?2hA|$NA|$rMQ$$$F!$%+%a%i$NFbIt%Q%i%a!<%?(B($BOD$_$"$j(B)$B$H30It%Q%i%a!<%?$r(Brefine$B$9$k!%(B
/*!
  $B4QB,9TNs$N7A<0$K$D$$$F$O(BcalibrateWithPlanes()$B$r;2>H!%(B
  \param K	 $B%+%a%i$NFbIt%Q%i%a!<%?$N=i4|CM(B($BOD$_?dDjCM$r4^$`(B)$B$rM?$($k!%(B
		 $B$=$N:GE,?dDjCM$,JV$5$l$k!%(B
  \param cameras $B3F;kE@$K$D$$$F!$;2>HJ?LL$K8GDj$5$l$?%o!<%k%I:BI87O$+$i8+$?(B
		 $B%+%a%i$N30It%Q%i%a!<%?$N=i4|CM$rM?$($k!%$=$N:GE,?dDjCM$,(B
		 $BJV$5$l$k!%(B
 */
void
MeasurementMatrix::refineCalibrationWithPlanes
    (CameraWithDistortion::Intrinsic& K,
     Array<CanonicalCamera>& cameras) const
{
    if (nframes() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::refineCalibrationWithPlanes: Two or more frames required!!");

    CostCP				err(*this, 7);
    NullConstraint<T, CostCP::ATA>	g;

    tuMinimizeSquareSparse(err, g, K, cameras, 500);
}

//! $B%+%a%i$N30It%Q%i%a!<%?$HFCD'E@0LCV$N=i4|CM$rHs@~7?:GE,2=$K$h$j(Brefine$B$9$k!%(B
/*!
  \param cameras	$B8D!9$N%+%a%i$N30It%Q%i%a!<%?$N=i4|CM$rM?$($k!%$=$N:G(B
			$BE,?dDjCM$,JV$5$l$k!%(B
  \param Xt		$BFCD'E@$N%f!<%/%j%C%I:BI87O$K$*$1$kF1<!:BI8$N=i4|CM$r(B
			$BM?$($k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
*/
void
MeasurementMatrix::bundleAdjustment(Array<CanonicalCamera>& cameras,
				      Matrix<T>& Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::bundleAdjustment: Two or more frames required!!");

    CostBA		err(*this, 0);
    CostBA::CostCD	g(cameras);
    Matrix<T>		shape(Xt, 0, 0, Xt.nrow(), 3);

    tuMinimizeSquareSparse(err, g, cameras, shape, 30);
}

//! $B%+%a%i$N30It%Q%i%a!<%?!$>GE@5wN%$*$h$SFCD'E@0LCV$N=i4|CM$rHs@~7?:GE,2=$K$h$j(Brefine$B$9$k!%(B
/*!
  \param cameras	$B8D!9$N%+%a%i$N30It%Q%i%a!<%?$H>GE@5wN%$N=i4|CM$rM?$((B
			$B$k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
  \param Xt		$BFCD'E@$N%f!<%/%j%C%I:BI87O$K$*$1$kF1<!:BI8$N=i4|CM$r(B
			$BM?$($k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
*/
void
MeasurementMatrix::bundleAdjustment
	(Array<CameraWithFocalLength>&	cameras,
	 Matrix<T>&				Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::bundleAdjustment (with focal lengths estimation): Two or more frames required!!");

    CostBA		err(*this, 1);
    CostBA::CostCD	g(cameras);
    Matrix<T>		shape(Xt, 0, 0, Xt.nrow(), 3);

    tuMinimizeSquareSparse(err, g, cameras, shape, 30);
}

//! $B%+%a%i$N30It%Q%i%a!<%?!$A4$F$N%+%a%i$K6&DL$J>GE@5wN%$*$h$SFCD'E@0LCV$N=i4|CM$rHs@~7?:GE,2=$K$h$j(Brefine$B$9$k!%(B
/*!
  \param K		$BA4$F$N%+%a%i$K6&DL$J>GE@5wN%$N=i4|CM$rM?$($k!%$=$N:G(B
			$BE,?dDjCM$,JV$5$l$k!%(B
  \param cameras	$B8D!9$N%+%a%i$N30It%Q%i%a!<%?$N=i4|CM$rM?$($k!%$=$N:G(B
			$BE,?dDjCM$,JV$5$l$k!%(B
  \param Xt		$BFCD'E@$N%f!<%/%j%C%I:BI87O$K$*$1$kF1<!:BI8$N=i4|CM$r(B
			$BM?$($k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
*/
void
MeasurementMatrix::bundleAdjustment(CameraWithFocalLength::Intrinsic& K,
				      Array<CanonicalCamera>& cameras,
				      Matrix<T>&		  Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::bundleAdjustment (with common focal length estimation): Two or more frames required!!");

    CostBACFE		err(*this);
    CostBACFE::ATA	params(K, cameras);
    CostBACFE::CostCD	g(params);
    Matrix<T>		shape(Xt, 0, 0, Xt.nrow(), 3);

    tuMinimizeSquareSparse(err, g, params, shape, 30);

    K = params.K;
    for (int i = 0; i < cameras.dim(); ++i)
	cameras[i] = params[i];
}

//! $B%+%a%i$N0LCV$,ITJQ$H$N2>Dj$N$b$H$G!$8D!9$N%+%a%i$N;Q@*$*$h$SFCD'E@0LCV$N=i4|CM$rHs@~7?:GE,2=$K$h$j(Brefine$B$9$k!%(B
/*!
  \param cameras	$B8D!9$N%+%a%i$N;Q@*$N=i4|CM$rM?$($k!%$=$N:GE,?dDjCM$,(B
			$BJV$5$l$k!%(B
  \param Xt		$BFCD'E@$N%f!<%/%j%C%I:BI87O$K$*$1$kF1<!:BI8$N=i4|CM$r(B
			$BM?$($k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
*/
void
MeasurementMatrix::bundleAdjustmentWithFixedCameraCenters
	(Array<CanonicalCamera>& cameras, Matrix<double>& Xt) const
{
    if (nframes() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::bundleAdjustment (with fixed camera centers): Three or more frames required!!");

    CostBA				err(*this, 0, 1);
    NullConstraint<T, CostBA::ATA>	g;
    Matrix<double>			shape(Xt, 0, 0, Xt.nrow(), 3);

    tuMinimizeSquareSparse(err, g, cameras, shape, 30);
}

//! $B%+%a%i$N0LCV$,ITJQ$H$N2>Dj$N$b$H$G!$8D!9$N%+%a%i$N;Q@*!$>GE@5wN%$*$h$SFCD'E@0LCV$N=i4|CM$rHs@~7?:GE,2=$K$h$j(Brefine$B$9$k!%(B
/*!
  \param cameras	$B8D!9$N%+%a%i$N;Q@*$H>GE@5wN%$N=i4|CM$rM?$($k!%$=$N:G(B
			$BE,?dDjCM$,JV$5$l$k!%(B
  \param Xt		$BFCD'E@$N%f!<%/%j%C%I:BI87O$K$*$1$kF1<!:BI8$N=i4|CM$r(B
			$BM?$($k!%$=$N:GE,?dDjCM$,JV$5$l$k!%(B
*/
void
MeasurementMatrix::bundleAdjustmentWithFixedCameraCenters
	(Array<CameraWithFocalLength>&	cameras,
	 Matrix<double>&			Xt) const
{
    if (nframes() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::bundleAdjustment (and focal lengths estimation with fixed camera centers): Three or more frames required!!");

    CostBA				err(*this, 1, 1);
    NullConstraint<T, CostBA::ATAF>	g;
    Matrix<double>			shape(Xt, 0, 0, Xt.nrow(), 3);

    tuMinimizeSquareSparse(err, g, cameras, shape, 30);
}

//! $BJ#?t$N2hA|>e$K$*$1$kFCD'E@$NA|$N0LCV$+$i!$$=$N(B3$B<!856u4V$G$N0LCV$r5a$a$k!%(B
/*!
  $B%+%a%i9TNs(BP$B$O!$$I$N$h$&$J6u4V(B($B<M1F6u4V!$%"%U%#%s6u4V!$%f!<%/%j%C%I6u4V(B)
  $B$K$*$1$kI=8=$G$b$h$$!%JV$5$l$kFCD'E@$N(B3$B<!850LCV$O!$(BP$B$HF1$86u4V$K$*$1$kI=8=(B
  $B$H$J$k!%(B
  \param P		$B%+%a%i9TNs(B($B%U%l!<%`?t$r(BF$B$H$9$k$H(B\f$3F\times 4\f$$B9TNs(B)
  \param inhomogeneous	$BJV$5$l$kF1<!:BI8$NBh(B4$BMWAG$r(B1$B$K$7$?$$>l9g$O(Btrue$B$rM?$((B
			$B$k!%(B
  \return		$BFCD'E@$N(B3$B<!850LCV$rI=$9F1<!:BI8!%(B
*/
Matrix<double>
MeasurementMatrix::reconstruction(const Matrix<T>& P,
				    bool inhomogeneous) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::reconstruction: Two or more frames required!!");

    Matrix<T>	Xt(npoints(), 4);
    for (int j = 0; j < npoints(); ++j)
    {
	Matrix<T>	Q(2*nframes(), 4);
	for (int i = 0; i < nframes(); ++i)
	{
	    Q[2*i  ] = P[3*i  ] - (*this)[j][3*i  ] * P[3*i+2];
	    Q[2*i+1] = P[3*i+1] - (*this)[j][3*i+1] * P[3*i+2];
	}
	Vector<T>	evalue;
	Xt[j] = (Q.trns() * Q).eigen(evalue)[3];
	if (inhomogeneous)
	    Xt[j] /= Xt[j][3];
    }
    
    return Xt;
}

//! $B4pAC9TNs$N@:EY$rI>2A$9$k!%(B
/*!
  $B@:EY$O!$BP1~$9$k(B2$BE@(B\f$\TUvec{u}{}\f$$B$HE@(B\f$\TUdvec{u}{}\f$$B$K$D$$$F(B
  \f$\TUvec{u}{}\f$$B$H$=$l$KBP1~$9$k%(%T%]!<%i@~(B\f$\TUdtvec{u}{}\TUvec{F}{}\f$
  $B$N5wN%$*$h$S(B\f$\TUdvec{u}{}\f$$B$H$=$l$KBP1~$9$k%(%T%]!<%i@~(B
  \f$\TUvec{F}{}\TUvec{u}{}\f$$B$N5wN%$r$9$Y$F$NFCD'E@$KBP$7$F7W;;$7!$$=$N(B
  $BJ?6Q5wN%$H$7$FI>2A$9$k!%(B
  \param F	$BI>2A$NBP>]$H$J$k4pAC9TNs!%(B
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	$BFCD'E@A|$H%(%T%]!<%i@~$H$N5wN%$NJ?6QCM!%(B
 */
double
MeasurementMatrix::assessFundamental(const Matrix<T>& F,
				       u_int frame0, u_int frame1) const
{
    T	d_sum = 0.0;
    for (int j = 0; j < npoints(); ++j)
    {
	const Vector<T>	&ldata = (*this)[j](3*frame0, 3),
				&rdata = (*this)[j](3*frame1, 3);
	Vector<T>		l = F * ldata;
	d_sum += fabs(rdata * l) / l(0, 2).length();
	l = rdata * F;
	d_sum += fabs(l * ldata) / l(0, 2).length();
    }

    return d_sum / (2*npoints());
}

//! $B<M1FJQ499TNs$N@:EY$rI>2A$9$k!%(B
/*!
  $B@:EY$O!$BP1~$9$k(B2$BE@(B\f$\TUvec{u}{}\f$$B$HE@(B\f$\TUdvec{u}{}\f$$B$K$D$$$F(B
  $BJQ49$7$?FCD'E@A|(B\f$\TUvec{H}{}\TUvec{u}{}\f$$B$H<B:]$NFCD'E@A|(B
  \f$\TUdvec{u}{}\f$$B$N5wN%$r$9$Y$F$NFCD'E@$KBP$7$F7W;;$7!$$=$NJ?6Q5wN%$H$7$F(B
  $BI>2A$9$k!%(B
  \param H	$BI>2A$NBP>]$H$J$k<M1FJQ499TNs!%(B
  \param frame0	$B0lJ}$N%U%l!<%`$N(Bindex.
  \param frame1	$B$b$&0lJ}$N%U%l!<%`$N(Bindex.
  \return	$BJQ49$7$?FCD'E@A|$H<B:]$NFCD'E@A|$H$N5wN%$NJ?6QCM!%(B
 */
double
MeasurementMatrix::assessHomography(const Matrix<T>& H,
				    u_int frame0, u_int frame1) const
{
    T	d_sum = 0.0;
    for (int j = 0; j < npoints(); ++j)
    {
	const Vector<T>	&ldata = (*this)[j](3*frame0, 3),
				&rdata = (*this)[j](3*frame1, 3);
	Vector<T>		r = H * ldata;
	r /= r[2];
	d_sum += r.dist(rdata);
    }

    return d_sum / npoints();
}

//! $B%+%a%i9TNs$N@:EY$rI>2A$9$k!%(B
/*!
  $B@:EY$O!$0[$J$k(B2$B$D$N%U%l!<%`4V$G(BassessFundamental()
  $B$K$h$C$FFCD'E@A|$H%(%T%]!<%i@~$N5wN%$NJ?6QCM$r7W;;$7!$$3$l$r$9$Y$F$N(B
  2$B%U%l!<%`$NAH9g$o$;$KBP$7$F9T$C$F$=$NJ?6QCM$K$h$jI>2A$9$k!%(B
  \param P	$BI>2A$NBP>]$H$J$k(B\f$3F\times 4\f$$B%+%a%i9TNs!%(B
  \return	$BFCD'E@A|$H%(%T%]!<%i@~$H$N5wN%$NJ?6QCM!%(B
 */
double
MeasurementMatrix::assessError(const Matrix<T>& P) const
{
    u_int	ncombinations = 0;
    T		d_sum = 0.0;
    for (int i = 0; i < nframes(); ++i)
    {
	const Matrix<T>&	Pi = P(3*i, 0, 3, 4);
	
	for (int j = i+1; j < nframes(); ++j)
	{
	    const Matrix<T>&		Pj = P(3*j, 0, 3, 4);
	    SVDecomposition<T>	svd(Pi);
	    const Vector<T>&		c = svd.Ut()[3];
	    const Matrix<T>&		F = (Pj * c).skew() * Pj * Pi.trns()
					  * (Pi * Pi.trns()).inv();

	    d_sum += assessFundamental(F, i, j);
	    ++ncombinations;
	}
    }
    
    return d_sum / ncombinations;
}

/*
 *  private member functions.
 */
std::istream&
MeasurementMatrix::get(std::istream& in, int j, u_int nfrms)
{
    char	c;

    while (in.get(c))			// Skip white spaces other than '\n'.
	if (!isspace(c) || c == '\n')
	    break;

    if (!in || c == '\n')
    {
	resize(j, 3*nfrms);
	return in;
    }
    in.putback(c);
    Vector<T> val;
    in >> val;
    get(in, j + 1, val.dim() / 2);
    for (int i = 0; i < nframes(); ++i)
    {
	(*this)[j][3*i]   = val[2*i];
	(*this)[j][3*i+1] = val[2*i+1];
	(*this)[j][3*i+2] = 1.0;
    }
    return in;
}    

Matrix<double>
MeasurementMatrix::initializeCalibrationWithPlanes
	(Array<CanonicalCamera>& cameras) const
{
    if (nframes() < 3)
	throw std::invalid_argument("TU::MeasurementMatrix::initializeCalibrationWithPlanes(): not enough frames!!");
	
  // Compute homographies.
    Array<Matrix<T> >	Ht(nframes()-1);
    for (int i = 0; i < Ht.dim(); ++i)
	Ht[i] = homography(0, i+1).trns();
    
  // Compute IAG(Image of Absolute Conic).
    Matrix<T>	A(2*Ht.dim(), 6);
    for (int i = 0; i < Ht.dim(); ++i)
    {
	A[2*i  ] = extpro(Ht[i][0], Ht[i][0]) - extpro(Ht[i][1], Ht[i][1]);
	A[2*i+1] = extpro(Ht[i][0], Ht[i][1]);
    }
    Vector<T>		evalue;
    const Matrix<T>&	evector = (A.trns() * A).eigen(evalue);
#ifdef DEBUG
    std::cerr << "TU::MeasurementMatrix::initializeCalibrationWithPlanes():: eigen values = "
	      << evalue;
#endif
    Vector<T>		a;
    if (Ht.dim() < 3)			// two reference planes.
	a = evector[5][1] * evector[4] - evector[4][1] * evector[5];
    else
	a = evector[5];
    Matrix<T>	omega(3, 3);
    omega[0][0]		      = a[0];
    omega[0][1] = omega[1][0] = a[1];
    omega[0][2] = omega[2][0] = a[2];
    omega[1][1]		      = a[3];
    omega[1][2] = omega[2][1] = a[4];
    omega[2][2]		      = a[5];
    if (omega.det() < 0.0)
	omega *= -1.0;

  // Extract intrinsic parameters from IAG.
    const Matrix<T>	&Kinv = omega.cholesky(), &Ktinv = Kinv.trns();
    cameras.resize(Ht.dim());
    for (int i = 0; i < cameras.dim(); ++i)
    {
	const Matrix<T>&	tmp = Ht[i] * Ktinv;
	SVDecomposition<T>	svd(tmp(0, 0, 2, 3));
	Matrix<T>		R(3, 3);
	R(0, 0, 2, 3) = svd.Vt().trns() * svd.Ut()(0, 0, 2, 3);
	R[2] = R[0] ^ R[1];
	
	const Vector<T>&	t = -2.0 / (svd[0] + svd[1]) * R * tmp[2];

	cameras[i].setRotation(R.trns()).setTranslation(t);
    }

    Matrix<T>	K(Kinv.inv());
    K /= K[2][2];

    return K;
}

void
MeasurementMatrix::initializeFocalLengthsEstimation(Matrix<T>& P,
						      Matrix<T>& Xt) const
{
    if (nframes() < 2)
	throw std::invalid_argument("TU::MeasurementMatrix::initializeFocalLengthsEstimation: Two or more frames required!!");
    
  // Compute B = H * Ht.
    Matrix<T>	Q(4*nframes(), 10);
    for (int i = 0; i < nframes(); ++i)
    {
	Q[4*i  ] = extpro(P[3*i  ], P[3*i  ])
		 - extpro(P[3*i+1], P[3*i+1]);		// aa - bb
	Q[4*i+1] = extpro(P[3*i  ], P[3*i+1]);		// ab
	Q[4*i+2] = extpro(P[3*i  ], P[3*i+2]);		// ac
	Q[4*i+3] = extpro(P[3*i+1], P[3*i+2]);		// bc
    }
    Vector<T>		evalue;
    const Vector<T>	b = (Q.trns()*Q).eigen(evalue)[9];
#ifdef DEBUG	// Check the order of magnitute of last two eigen values.
    std::cerr << "evalue of QQ = " << evalue;
#endif
    Matrix<T>	B(4, 4);
    B[0][0]	      = b[0];
    B[0][1] = B[1][0] = b[1];
    B[0][2] = B[2][0] = b[2];
    B[0][3] = B[3][0] = b[3];
    B[1][1]	      = b[4];
    B[1][2] = B[2][1] = b[5];
    B[1][3] = B[3][1] = b[6];
    B[2][2]	      = b[7];
    B[2][3] = B[3][2] = b[8];
    B[3][3]	      = b[9];

  // Factor B into H and Ht.
    Matrix<T>	Ht = B.eigen(evalue);
    if (evalue[0] < 0.0)
	evalue *= -1.0;
#ifdef DEBUG	// First three eigen values should be positive.
    std::cerr << "evalue of B = " << evalue;
#endif
    for (int i = 0; i < 3; ++i)
	if (evalue[i] >= 0.0)
	    Ht[i] *= sqrt(evalue[i]);
	else
	    throw std::domain_error("TU::MeasurementMatirix::initializeFocalLengtsEstimation: cannot decompose B into H and Ht due to nevative eigen value!");

  // Compute inverse of H.
    Matrix<T>	Hinv(4, 4);
    (Hinv[0] = Ht[0]) /= evalue[0];
    (Hinv[1] = Ht[1]) /= evalue[1];
    (Hinv[2] = Ht[2]) /= evalue[2];
    Hinv[3] = Ht[3];

  // Transform from projective to Euclidean world.
    P  *= Ht.trns();
    Xt *= Hinv.trns();

  // Ensure left 3*3 part of the projection matrix having positive determinant.
    for (int i = 0; i < nframes(); ++i)
	if (P(3*i, 0, 3, 3).det() < 0.0)
	    P(3*i, 0, 3, 4) *= -1.0;
    
  // Normalize so that 4-th component of Xt is 1.
    for (int j = 0; j < npoints(); ++j)
	Xt[j] /= Xt[j][3];

  // Determine sign of Ht[3] so that object points are in front of the cameras.
    int	sig = 0;
    for (int i = 0; i < nframes(); ++i)
	for (int j = 0; j < npoints(); ++j)
	{
#ifdef DEBUG
	    std::cerr << "depth[" << i << "][" << j << "] = " << P[3*i+2]*Xt[j]
		      << std::endl;
#endif
	    sig += sign(P[3*i+2] * Xt[j]);
	}
    if (sig < 0)
    {
	P(0, 3, P.nrow(), 1)   *= -1.0;
	Xt(0, 0, Xt.nrow(), 3) *= -1.0;
    }
}

/************************************************************************
*  class MeasurementMatrix::CostPF					*
************************************************************************/
MeasurementMatrix::CostPF::CostPF(const Matrix<double>& Wt0,
				    u_int niter_max)
    :_niter_max(niter_max), _Wt0(Wt0), _s(), _Ut(), _Vt()
{
}

inline int
MeasurementMatrix::CostPF::frame_index(u_int p) const
{
    return p / (npoints()-1) + 1;
}

inline int
MeasurementMatrix::CostPF::point_index(u_int p) const
{
    return p % (npoints()-1) + 1;
}

/*
 *  Value of the cost function: J = \sum_{n=4}^N \sigma^2
 */
double
MeasurementMatrix::CostPF::operator ()(const Vector<double>& mu)
{
  /* Compute projective measurement matrix. */
    Matrix<double>	Wt(_Wt0);
    for (int p = 0; p < mu.dim(); ++p)
	Wt[point_index(p)](3*frame_index(p), 3) *= mu[p];

  /* Decompose it by using SVD. */
    SVDecomposition<double>	svd(Wt.trns());
    _s  = svd.diagonal()(0, N());
    _Ut = svd.Ut();
    _Vt = svd.Vt();

    return _s(4, N()-4).square();
}

/*
 *  \TUvec{S}{m\dot np} =
 *	\TUtvec{v}{n}\TUpartial{\TUvec{W}{}}{\mu_p}\TUvec{u}{m}
 */
Matrix<double>
MeasurementMatrix::CostPF::S(int m) const
{
    Matrix<double>	Sm(_Vt.nrow(), (nframes()-1)*(npoints()-1));
    
    for (int p = 0; p < Sm.ncol(); ++p)
    {
	int	i = frame_index(p), j = point_index(p);
	for (int n = 0; n < _Vt.nrow(); ++n)
	    Sm[n][p] = _Vt[n](3*i, 3) * _Wt0[j](3*i, 3) * _Ut[m][j];
    }
    return Sm;
}

/*
 *  \TUvec{T}{m\dot np} =
 *	\TUtvec{v}{m}\TUpartial{\TUvec{W}{}}{\mu_p}\TUvec{u}{n}
 */
Matrix<double>
MeasurementMatrix::CostPF::T(int m) const
{
    Matrix<double>	Tm(_Ut.nrow(), (nframes()-1)*(npoints()-1));

    for (int p = 0; p < Tm.ncol(); ++p)
    {
	int	i = frame_index(p), j = point_index(p);
	for (int n = 0; n < _Ut.nrow(); ++n)
	    Tm[n][p] = _Vt[m](3*i, 3) * _Wt0[j](3*i, 3) * _Ut[n][j];
    }
    return Tm;
}

double
MeasurementMatrix::CostPF::minimize(Vector<double>& mu)
{
    double	val = (*this)(mu);
    
    for (int iter = 0; iter < _niter_max; ++iter)
    {
	print(iter, val, mu);
      //      	std::cerr << "    sigma = " << _s;
	
      /* Make a gradient vector and a Hessian matrix. */
	Vector<double>	grad(mu.dim());
	Matrix<double>	H(mu.dim(), mu.dim()), dH(mu.dim(), mu.dim());
	
	for (int p = 0; p < H.nrow(); ++p)
	{
	    H[p][p] = _Wt0[point_index(p)](3*frame_index(p), 3).square();
	    grad[p] = mu[p] * H[p][p];
	}
	
	for (int m = 0; m < 4; ++m)	// for first four sigular values...
	{
	    const Matrix<double>	&Sm = S(m), &Tm = T(m);
	    
	    for (int p = 0; p < H.nrow(); ++p)
	    {
		for (int q = 0; q <= p; ++q)
		{
		    double	ddsigma = 0.0;	// 2nd derivative of S-value
		    for (int n = 0; n < N(); ++n)
			if (n != m)
			    ddsigma += ((_s[m] * (Sm[n][p] * Sm[n][q] +
						  Tm[n][p] * Tm[n][q]) +
					 _s[n] * (Sm[n][p] * Tm[n][q] +
						  Tm[n][p] * Sm[n][q])) /
					((_s[m] + _s[n]) * (_s[m] - _s[n])));
		    H[p][q] -= (Sm[m][p] * Sm[m][q] + _s[m] * ddsigma);

		    ddsigma = 0.0;		// small term of ddsigma
		    for (int n = N(); n < _Vt.nrow(); ++n)
			ddsigma += (Sm[n][p] * Sm[n][q] / _s[m]);
		    for (int n = N(); n < _Ut.nrow(); ++n)
			ddsigma += (Tm[n][p] * Tm[n][q] / _s[m]);
		    dH[p][q] -= (_s[m] * ddsigma);
		}
		grad[p] -= (_s[m] * Sm[m][p]);
	    }
	}
	H.symmetrize();
	dH.symmetrize();

      /* Update projective depths: mu */
	Vector<double>	dmu = grad;
	double			val_next = (*this)(mu - dmu.solve(H+dH));
	if (val_next > val)
	{
	    dmu = grad;
	    val_next = (*this)(mu - dmu.solve(H));
	}
	mu -= dmu;

	if (fabs(val_next - val) < 1.0e-12)
	    return val_next;
	
	val = val_next;
    }

    std::cerr << "TU::MeasurementMatrix::CostPF::minimize: Too many iterations!!"
	      << std::endl;
    return val;
}

void
MeasurementMatrix::CostPF::print(int i, double val,
				   const Vector<double>& mu) const
{
#ifdef DEBUG
    std::cerr << std::setw(3) << i << ": (" << val << ')' << mu;
#else
    std::cerr << i << '\t' << val << std::endl;
#endif
}

/************************************************************************
*  class MeasurementMatrix::CostPM					*
************************************************************************/
Vector<double>
MeasurementMatrix::CostPM::operator ()(const AT& p) const
{
    Vector<T>	q(4);
    q(0, 3) = p;
    q[3] = p * p;

    return _AA * q - _b;
}

Matrix<double>
MeasurementMatrix::CostPM::jacobian(const AT& p) const
{
    Matrix<T>	J;
    J = _AA(0, 0, _AA.nrow(), 3);
    for (int i = 0; i < J.nrow(); ++i)
	J[i] += (2.0 * _AA[i][3] * p);

    return J;
}

/************************************************************************
*  class MeasurementMatrix::CostCP					*
************************************************************************/
Vector<double>
MeasurementMatrix::CostCP::operator ()(const ATA& K,
					 const ATB& camera, int i) const
{
  // Compute errors for all points.
    Vector<T>	val(2 * npoints());
    for (int j = 0; j < npoints(); ++j)
    {
	Vector<T>	X(3);
	X(0, 2) = _Wt[j](0, 2);		// 3D reference points (X[2] = 0.0)
      	const Point2<T>&	u = K(camera.xc(X));
	val[2*j]   = u[0] - _Wt[j][3*(i+1)];
	val[2*j+1] = u[1] - _Wt[j][3*(i+1)+1];
    }

    return val;
}

Matrix<double>
MeasurementMatrix::CostCP::jacobianA(const ATA& K, const ATB& camera) const
{
    Matrix<T>		J(2 * npoints(), adim());
    for (int j = 0; j < npoints(); ++j)
    {
	Vector<T>	X(3);
	X(0, 2) = _Wt[j](0, 2);		// 3D reference points (X[2] = 0.0)
	J(2*j, 0, 2, J.ncol()) = K.jacobianK(camera.xc(X));
    }
    
    return J;
}

Matrix<double>
MeasurementMatrix::CostCP::jacobianB(const ATA& K, const ATB& camera) const
{
    Matrix<T>	KK(2 * npoints(), 6);
    for (int j = 0; j < npoints(); ++j)
    {
	Vector<T>	X(3);
	X(0, 2) = _Wt[j](0, 2);		// 3D reference points (X[2] = 0.0)
	KK(2*j, 0, 2, KK.ncol())
	    = K.jacobianXC(camera.xc(X)) * camera.jacobianPc(X);
    }
    
    return KK;
}

void
MeasurementMatrix::CostCP::updateA(ATA& K, const Vector<T>& dK) const
{
    K.update(dK);
}

void
MeasurementMatrix::CostCP::updateB(ATB& camera,
				     const Vector<T>& dcamera) const
{
    camera.update(dcamera);
}

/************************************************************************
*  class MeasurementMatrix::CostBA					*
************************************************************************/
MeasurementMatrix::CostBA::CostBA(const MeasurementMatrix& Wt,
				    u_int dofIntrinsic, int fixCameraCenter)
    :_Wt(Wt), _fcc(fixCameraCenter), _adim(0), _adims(nframes())
{
    _adims[0] = dofIntrinsic + (_fcc ? 3 : 0);
    _adim += _adims[0];
    for (int i = 1; i < _adims.dim(); ++i)
    {
	_adims[i] = dofIntrinsic + (_fcc ? 3 : 6);
	_adim += _adims[i];
    }
}

Vector<double>
MeasurementMatrix::CostBA::operator ()(const ATA& p,
					 const ATB& x, int j) const
{
    Vector<T>	val(2 * nframes());
    for (int i = 0; i < nframes(); ++i)
    {
	const Point2<T>&	u = p[i](x);
	val[2*i]   = u[0] - _Wt[j][3*i];
	val[2*i+1] = u[1] - _Wt[j][3*i+1];
    }

    return val;
}

BlockMatrix<double>
MeasurementMatrix::CostBA::jacobianA(const ATA& p, const ATB& x) const
{
    BlockMatrix<T>	J(nframes());
    J[0] = (_fcc ? p[0].jacobianFCC(x) : p[0].jacobianK(x));
    for (int i = 1; i < nframes(); ++i)
	J[i] = (_fcc ? p[i].jacobianFCC(x) : p[i].jacobianP(x));
    
    return J;
}

Matrix<double>
MeasurementMatrix::CostBA::jacobianB(const ATA& p, const ATB& x) const
{
    Matrix<T>	K(2 * nframes(), 3);
    for (int i = 0; i < nframes(); ++i)
	K(2*i, 0, 2, K.ncol()) = p[i].jacobianX(x);

    return K;
}

void
MeasurementMatrix::CostBA::updateA(ATA& p, const Vector<T>& dp) const
{
    int	d = 0;
    if (_fcc)
	p[0].updateFCC(dp(d, _adims[0]));
    else
	p[0].updateIntrinsic(dp(d, _adims[0]));
    d += _adims[0];
    for (int i = 1; i < _adims.dim(); ++i)
    {
      // Update camera parameters.
	if (_fcc)
	    p[i].updateFCC(dp(d, _adims[i]));
	else
	    p[i].update(dp(d, _adims[i]));
	d += _adims[i];
    }

    if (!_fcc)
    {
      // Force the distance between 0th and 1st cameras keep constant.
	Vector<T>	t = p[1].t() - p[0].t();
	t.normalize() *= p[0].t().dist(p[1].t());
	p[1].setTranslation(p[0].t() + t);
    }
}

void
MeasurementMatrix::CostBA::updateB(ATB& x, const Vector<T>& dx) const
{
    x -= dx;
}

Vector<double>
MeasurementMatrix::CostBA::operator ()(const ATAF& p,
					 const ATB& x, int j) const
{
    Vector<T>	val(2 * nframes());
    for (int i = 0; i < nframes(); ++i)
    {
	const Point2<T>&	u = p[i](x);
	val[2*i]   = u[0] - _Wt[j][3*i];
	val[2*i+1] = u[1] - _Wt[j][3*i+1];
    }
    
    return val;
}

BlockMatrix<double>
MeasurementMatrix::CostBA::jacobianA(const ATAF& p, const ATB& x) const
{
    BlockMatrix<T>	J(nframes());
    J[0] = (_fcc ? p[0].jacobianFCC(x) : p[0].jacobianK(x));
    for (int i = 1; i < nframes(); ++i)
	J[i] = (_fcc ? p[i].jacobianFCC(x) : p[i].jacobianP(x));
    
    return J;
}

Matrix<double>
MeasurementMatrix::CostBA::jacobianB(const ATAF& p, const ATB& x) const
{
    Matrix<T>	K(2 * nframes(), 3);
    for (int i = 0; i < nframes(); ++i)
	K(2*i, 0, 2, K.ncol()) = p[i].jacobianX(x);

    return K;
}

void
MeasurementMatrix::CostBA::updateA(ATAF& p, const Vector<T>& dp) const
{
    int	d = 0;
    if (_fcc)
	p[0].updateFCC(dp(d, _adims[0]));
    else
	p[0].updateIntrinsic(dp(d, _adims[0]));
    d += _adims[0];
    for (int i = 1; i < _adims.dim(); ++i)
    {
      // Update camera parameters.
	if (_fcc)
	    p[i].updateFCC(dp(d, _adims[i]));
	else
	    p[i].update(dp(d, _adims[i]));
	d += _adims[i];
    }

    if (!_fcc)
    {
      // Force the distance between 0th and 1st cameras keep constant.
	Vector<T>	t = p[1].t() - p[0].t();
	t.normalize() *= p[0].t().dist(p[1].t());
	p[1].setTranslation(p[0].t() + t);
    }
}

/************************************************************************
*  class MeasurementMatrix::CostBACFE					*
************************************************************************/
Vector<double>
MeasurementMatrix::CostBACFE::operator ()
    (const ATA& p, const ATB& x, int j) const
{
    Vector<T>	val(2 * nframes());
    for (int i = 0; i < nframes(); ++i)
    {
	const Point2<T>&	u = p.K(p[i].xc(x));
	val[2*i]   = u[0] - _Wt[j][3*i];
	val[2*i+1] = u[1] - _Wt[j][3*i+1];
    }

    return val;
}

Matrix<double>
MeasurementMatrix::CostBACFE::jacobianA(const ATA& p, const ATB& x) const
{
    Matrix<T>		J(2 * nframes(), adim());
    J(0, 0, 2, 1) = p.K.jacobianK(p[0].xc(x));
    for (int i = 1; i < nframes(); ++i)
    {
	const Point2<T>& xc = p[i].xc(x);
	J(2*i, 0, 2, 1) = p.K.jacobianK(xc);
	J(2*i, 1 + 6*(i-1), 2, 6) = p.K.jacobianXC(xc) * p[i].jacobianPc(x);
    }
    
    return J;
}

Matrix<double>
MeasurementMatrix::CostBACFE::jacobianB(const ATA& p,	const ATB& x) const
{
    Matrix<T>	K(2 * nframes(), 3);
    for (int i = 0; i < nframes(); ++i)
	K(2*i, 0, 2, K.ncol())
	    = p.K.jacobianXC(p[i].xc(x)) * p[i].jacobianXc(x);

    return K;
}

void
MeasurementMatrix::CostBACFE::updateA(ATA& p,	const Vector<T>& dp) const
{
    p.K.update(dp(0, 1));
    for (int i = 1; i < nframes(); ++i)
	p[i].update(dp(1 + 6*(i-1), 6));

  // Force the distance between 0th and 1st cameras keep constant.
    Vector<T>	t = p[1].t() - p[0].t();
    t.normalize() *= (p[0].t().dist(p[1].t()));
    p[1].setTranslation(p[0].t() + t);
}

void
MeasurementMatrix::CostBACFE::updateB(ATB& x, const Vector<T>& dx) const
{
    x -= dx;
}

}
#ifdef __GNUG__
#  include "TU/Array++.cc"

namespace TU
{
template 
Array<CanonicalCamera>::Array<CanonicalCamera>
	(const Array<CanonicalCamera>& )			;
template 
Array<CameraWithFocalLength>::Array<CameraWithFocalLength>
	(const Array<CameraWithFocalLength>& )			;
}

#endif
