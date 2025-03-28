/*!
  \file		FIRFilter.h
  \author	Toshio UESHIBA
  \brief	一般的なfinite impulse response filterを表すクラスの定義と実装
*/
#ifndef	TU_FIRFILTER_H
#define	TU_FIRFILTER_H

#include "TU/SeparableFilter2.h"

namespace TU
{
/************************************************************************
*  class fir_filter_iterator<D, COEFF, ITER>				*
************************************************************************/
//! データ列中の指定された要素に対してfinite impulse response filterを適用した結果を返す反復子
/*!
  \param D	フィルタの階数
  \param COEFF	フィルタのz変換係数
  \param ITER	データ列中の要素を指す定数反復子の型
*/
template <size_t D, class COEFF, class ITER>
class fir_filter_iterator
    : public boost::iterator_adaptor<fir_filter_iterator<D, COEFF, ITER>,
				     ITER,
				     replace_element<iterator_substance<ITER>,
						     std::iter_value_t<COEFF> >,
				     boost::forward_traversal_tag,
				     replace_element<iterator_substance<ITER>,
						     std::iter_value_t<COEFF> > >
{
  private:
    using super	= boost::iterator_adaptor<
			fir_filter_iterator,
			ITER,
			replace_element<iterator_substance<ITER>,
					std::iter_value_t<COEFF> >,
			boost::forward_traversal_tag,
			replace_element<iterator_substance<ITER>,
					std::iter_value_t<COEFF> > >;
    
  public:
    using	typename super::value_type;
    using	typename super::reference;

    friend	class boost::iterator_core_access;

  private:
    using buf_type	= std::array<value_type, D>;
    template <size_t I_>
    using index		= std::integral_constant<size_t, I_>;
    
  public:
		fir_filter_iterator(const ITER& iter, COEFF c)
		    :super(iter), _c(c), _ibuf(), _n(D-1)
		{
		    set_inpro(index<0>());
		    copy<D-1>(super::base(), D-1, _ibuf.begin());
		    std::advance(super::base_reference(), D-1);
		}
		fir_filter_iterator(const ITER& iter)
		    :super(iter), _c(), _ibuf(), _n(0)
		{
		}
		fir_filter_iterator(const fir_filter_iterator& iter)
		    :super(iter), _c(iter._c), _ibuf(iter._ibuf), _n(iter._n)
		{
		    set_inpro(index<0>());
		}
    fir_filter_iterator&
		operator =(const fir_filter_iterator& iter)
		{
		    if (this != &iter)
		    {
			super::operator =(iter);
			_c    = iter._c;
			_ibuf = iter._ibuf;
			_n    = iter._n;
			set_inpro(index<0>());
		    }

		    return *this;
		}
    
  private:
    reference	dereference() const
		{
		    _ibuf[_n] = *super::base();
		    return (this->*_inpro[_n])(index<0>());
		}
    void	increment()
		{
		    ++super::base_reference();
		    if (++_n == D)
			_n = 0;
		}

    void	set_inpro(index<D>)
		{
		}
    template <size_t N_>
    void	set_inpro(index<N_>)
		{
		    _inpro[N_] = &fir_filter_iterator::inpro<N_>;
		    set_inpro(index<N_+1>());
		}
    
    template <size_t N_>
    value_type	inpro(index<D-1>) const
		{
		    constexpr size_t	J = N_%D;
		    return _c[D-1]*_ibuf[J];
		}
    template <size_t N_, size_t I_>
    value_type	inpro(index<I_>) const
		{
		    constexpr size_t	J = (N_ + I_ + 1)%D;
		    return _c[I_]*_ibuf[J] + inpro<N_>(index<I_+1>());
		}
    
  private:
    using	fptr = value_type (fir_filter_iterator::*)(index<0>) const;
    
    std::array<fptr, D>	_inpro;	//!< 内積関数へのポインタテーブル
    const COEFF		_c;	//!< 先頭のフィルタ係数を指す反復子
    mutable buf_type	_ibuf;	//!< 過去D時点の入力データ
    size_t		_n;	//!< 最新の入力データへのindex
};

//! finite impulse response filter反復子を生成する
/*!
  \param iter	コンテナ中の要素を指す定数反復子
  \param c	先頭の入力フィルタ係数を指す反復子
  \return	finite impulse response filter反復子
*/
template <size_t D, class COEFF, class ITER> fir_filter_iterator<D, COEFF, ITER>
make_fir_filter_iterator(const ITER& iter, COEFF c)
{
    return {iter, c};
}

//! finite impulse response filter反復子(終端)を生成する
/*!
  \param iter	コンテナ中の要素を指す定数反復子
  \return	finite impulse response filter反復子(終端)
*/
template <size_t D, class COEFF, class ITER> fir_filter_iterator<D, COEFF, ITER>
make_fir_filter_iterator(ITER iter)
{
    return {iter};
}

/************************************************************************
*  class FIRFilter<D, T>						*
************************************************************************/
//! 片側Infinite Inpulse Response Filterを表すクラス
template <size_t D, class T=float>
class FIRFilter
{
  public:
    using element_type	= T;
    using coeffs_type	= std::array<T, D>;

    FIRFilter&	initialize(const T c[D])				;
    void	limits(T& limit0, T& limit1, T& limit2)		  const	;
    template <class IN, class OUT>
    OUT		convolve(IN ib, IN ie, OUT out, bool shift=false) const	;

    const coeffs_type&		c()		const	{return _c;}

  //! 与えられた長さの入力データ列に対する出力データ列の長さを返す
  /*!
    \param inSize	入力データ列の長さ
    \return		出力データ列の長さ
  */
    constexpr static size_t	outSize(size_t inSize)	{return inSize + 1 - D;}
    constexpr static size_t	offset()		{return D/2;}
    constexpr static size_t	winSize()		{return D;}
	
  private:
    coeffs_type	_c;			//!< フィルタ係数
};

//! フィルタのz変換係数をセットする
/*!
  \param c	z変換係数. z変換関数は
		\f[
		  H(z^{-1}) = c_{D-1} + c_{D-2}z^{-1} + c_{D-3}z^{-2} +
		  \cdots + c_{0}z^{-(D-1)}
		\f]
  \return	このフィルタ自身
*/
template <size_t D, class T> FIRFilter<D, T>&
FIRFilter<D, T>::initialize(const T c[D])
{
    copy<D>(c, D, _c.begin());

    return *this;
}

//! 特定の入力データ列に対してフィルタを適用した場合の極限値を求める
/*!
  \param limit0		一定入力 in(n) = 1 を与えたときの出力極限値を返す．
  \param limit1		傾き一定入力 in(n) = n を与えたときの出力極限値を返す．
  \param limit2		2次入力 in(n) = n^2 を与えたときの出力極限値を返す．
*/
template <size_t D, class T> void
FIRFilter<D, T>::limits(T& limit0, T& limit1, T& limit2) const
{
    T	x0 = 0, x1 = 0, x2 = 0;
    for (size_t i = 0; i < D; ++i)
    {
	x0 +=	      _c[i];
	x1 +=	    i*_c[D-1-i];
	x2 += (i-1)*i*_c[D-1-i];
    }
    limit0 =  x0;
    limit1 = -x1;
    limit2 =  x1 + x2;
}

//! フィルタを適用する
/*!
  \param ib	入力データ列の先頭を指す反復子
  \param ie	入力データ列の末尾の次を指す反復子
  \param out	出力データ列の先頭を指す反復子
  \param shift	true ならば，入力データと対応するよう，出力位置を
		offset() だけシフトする
  \return	出力データ列の末尾の次を指す反復子
*/
template <size_t D, class T> template <class IN, class OUT> OUT
FIRFilter<D, T>::convolve(IN ib, IN ie, OUT out, bool shift) const
{
    using citerator	= typename coeffs_type::const_iterator;

    if (shift)
	std::advance(out, offset());
    
    return std::copy(make_fir_filter_iterator<D>(ib, _c.begin()),
		     make_fir_filter_iterator<D, citerator>(ie), out);
}

/************************************************************************
*  class FIRFilter2<D, T>						*
************************************************************************/
//! 2次元Finite Inpulse Response Filterを表すクラス
template <size_t D, class T=float>
class FIRFilter2 : public SeparableFilter2<FIRFilter<D, T> >
{
  private:
    using fir_type	= FIRFilter<D, T>;
    using super		= SeparableFilter2<fir_type>;

  public:
    using element_type	= typename fir_type::element_type;
    using coeffs_type	= typename fir_type::coeffs_type;

  public:
    FIRFilter2&		initialize(const T cH[], const T cV[])		;

    using		super::filterH;
    using		super::filterV;

    const coeffs_type&		cH()	const	{return filterH().c();}
    const coeffs_type&		cV()	const	{return filterV().c();}
    constexpr static size_t	winSizeH()	{return fir_type::winSize();}
    constexpr static size_t	winSizeV()	{return fir_type::winSize();}
};
    
//! フィルタのz変換係数をセットする
/*!
  \param cH	横方向z変換係数
  \param cV	縦方向z変換係数
  \return	このフィルタ自身
*/
template <size_t D, class T> inline FIRFilter2<D, T>&
FIRFilter2<D, T>::initialize(const T cH[], const T cV[])
{
    filterH().initialize(cH);
    filterV().initialize(cV);

    return *this;
}

}	// namespace TU
#endif	// !TU_FIRFILTER_H
