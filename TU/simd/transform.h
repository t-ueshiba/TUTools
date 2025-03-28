/*!
  \file		transform.h
  \author	Toshio UESHIBA
  \brief	複数のSIMDベクトルを成分数が等しくなるように同時に変換する関数の定義
*/
#if !defined(TU_SIMD_TRANSFORM_H)
#define	TU_SIMD_TRANSFORM_H

#include "TU/iterator.h"
#include "TU/functional.h"
#include "TU/simd/cvt.h"

namespace TU
{
namespace simd
{
namespace detail
{
  //! 複数の入力反復子からのデータを関数に入力し，その結果を出力反復子に書き出す
  /*!
    入力反復子毎に異なるSIMDベクトルの型をvec<T>型に統一して変換関数に入力し，
    その結果を再び変換して出力反復子に書き出す．
    \param T		関数に入力するSIMDベクトルの要素型
    \param ASSIGN	変換結果を出力反復子に書き出す際に用いる代入演算子
    \param MASK		型変換をマスクベクトルとして行うならtrue, そうでなければfalse
    \param ITER_TUPLE	複数の入力反復子をまとめたタプル
    \param OUT		関数が出力するSIMDベクトルの変換結果を書き出す出力反復子
    \param FUNC		変換関数
  */
  template <class T, class ASSIGN, bool MASK,
	    class ITER_TUPLE, class OUT, class FUNC>
  class transformer
  {
    public:
      using head_iterator = tuple_head<ITER_TUPLE>;
      
    private:
      template <class T_, class=void>
      struct vec_element
      {
	  using type = typename tuple_head<T_>::element_type;
      };
      template <class DUMMY_>
      struct vec_element<void, DUMMY_>
      {
	  using type = void;
      };
      
    //using E = typename vec_element<std::iter_value_t<OUT> >::type;
      using E = typename vec_element<typename std::iterator_traits<OUT>::value_type>::type;
      
    //! 出力反復子に書き出すSIMDベクトルの要素型
      using O = std::conditional_t<std::is_void<E>::value, T, E>;
    //! 変換関数に入力するSIMDベクトルの要素型
      using I = std::conditional_t<std::is_void<T>::value, O, T>;

    private:
    // _funcの適用結果を vec<O> にconvert upしてstoreする
      template <class TUPLE_>
      std::enable_if_t<(vec<O>::size == tuple_head<TUPLE_>::size)>
		upResult_store(const TUPLE_& x)
		{
		    ASSIGN()(*_out, simd::cvtup<O, false, MASK>(x));
		    ++_out;
		}
      template <class TUPLE_>
      std::enable_if_t<(vec<O>::size < tuple_head<TUPLE_>::size)>
		upResult_store(const TUPLE_& x)
		{
		    upResult_store(simd::cvtup<O, false, MASK>(x));
		    upResult_store(simd::cvtup<O, true,  MASK>(x));
		}

    // 既に vec<O> と同位に到達している入力をさらに vec<I> にconvert upして
    // _funcを適用し，その結果を vec<O> までconvert downして返す
      template <class TUPLE_,
		std::enable_if_t<(vsize<TUPLE_>::max == vec<I>::size)>*
		= nullptr>
      auto	upArg_downResult(TUPLE_&& x)
		{
		    return simd::cvtdown<O, MASK>(
				_func(simd::cvtdown<I, MASK>(x)));
		}
      template <class TUPLE_,
		std::enable_if_t<(vsize<TUPLE_>::max > vec<I>::size)>*
		= nullptr>
      auto	upArg_downResult(TUPLE_&& x)
		{
		    constexpr auto	N = vsize<TUPLE_>::max;

		    const auto	y = upArg_downResult(
					simd::cvtup<I, false, MASK, N/2>(x));
		    const auto	z = upArg_downResult(
					simd::cvtup<I, true,  MASK, N/2>(x));

		  // 戻り値のベクトルは TUPLE_ と同位
		    return simd::cvtdown<O, MASK>(y, z);
		}

    // vec<I> と同位に達した入力に _func を適用し，
    // その結果をさらに vec<O> にconvert upしてstoreする．
      template <class TUPLE_>
      std::enable_if_t<(vsize<TUPLE_>::max == vec<I>::size &&
			vsize<TUPLE_>::max >  vec<O>::size)>
		exec(TUPLE_&& x)
		{
		    upResult_store(_func(simd::cvtdown<I, MASK>(x)));
		}

    // vec<O> と同位に達した入力をさらに vec<I> にconvert upして _func を
    // 適用し，その結果を vec<O> にconvert downしてstoreする．
      template <class TUPLE_>
      std::enable_if_t<(vsize<TUPLE_>::max == vec<O>::size)>
		exec(TUPLE_&& x)
		{
		    ASSIGN()(*_out, upArg_downResult(x));
		    ++_out;
		}

    // x の中に vec<O> より下位のベクトルがあれば読み込んでconvert upする．
      template <class TUPLE_>
      std::enable_if_t<(vsize<TUPLE_>::max > vec<I>::size &&
			vsize<TUPLE_>::max > vec<O>::size)>
		exec(TUPLE_&& x)
		{
		    constexpr auto	N = vsize<TUPLE_>::max;
		    
		    exec(simd::cvtup<I, false, MASK, N/2>(x));
		    exec(simd::cvtup<I, true,  MASK, N/2>(x));
		}

    public:
		transformer(ITER_TUPLE t, head_iterator end, OUT out, FUNC func)
		    :_t(t), _end(end), _out(out), _func(func)		{}

      OUT	operator ()()
		{
		    while (std::get<0>(_t) != _end)
		    {
		      // 入力ベクトルと出力ベクトルの中で最下位のベクトルの要素数
			constexpr auto
			    N = (vsize<ITER_TUPLE>::max > vec<O>::size ?
				 vsize<ITER_TUPLE>::max : vec<O>::size);

			exec(simd::cvtup<I, false, MASK, N>(_t));
		    }
		    return _out;
		}

    private:
      ITER_TUPLE		_t;
      const head_iterator	_end;
      OUT			_out;
      FUNC			_func;
  };
}	// namespace detail

template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class FUNC, class OUT, class IN, class... IN_EXTRA> inline OUT
transform(FUNC func, OUT out, IN in, IN end, IN_EXTRA... in_extra)
{
    detail::transformer<T, ASSIGN, MASK,
			std::tuple<IN, IN_EXTRA...>, OUT, FUNC>
	trns(std::make_tuple(in, in_extra...), end, out, func);
    return trns();
}
    
template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class ITER_TUPLE, class OUT, class FUNC> inline OUT
transform(ITER_TUPLE ib, ITER_TUPLE ie, OUT out, FUNC func)
{
    detail::transformer<T, ASSIGN, MASK, ITER_TUPLE, OUT, FUNC>
	trns(ib, std::get<0>(ie), out, func);
    return trns();
}
    
template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class ITER_TUPLE, class OUT, class FUNC> inline OUT
transform(zip_iterator<ITER_TUPLE> ib,
	  zip_iterator<ITER_TUPLE> ie, OUT out, FUNC func)
{
    return transform<T, ASSIGN, MASK>(ib.get_iterator_tuple(),
				      ie.get_iterator_tuple(), out, func);
}
    
template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class OUT, class IN, class... IN_EXTRA> inline OUT
copy(OUT out, IN in, IN end, IN_EXTRA... in_extra)
{
    detail::transformer<T, ASSIGN, MASK,
			std::tuple<IN, IN_EXTRA...>, OUT, identity>
	trns(std::make_tuple(in, in_extra...), end, out, identity());
    return trns();
}
    
template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class ITER_TUPLE, class OUT> inline OUT
copy(ITER_TUPLE ib, ITER_TUPLE ie, OUT out)
{
    detail::transformer<T, ASSIGN, MASK, ITER_TUPLE, OUT, identity>
	trns(ib, std::get<0>(ie), out, identity());
    return trns();
}
    
template <class T=void, class ASSIGN=assign, bool MASK=false,
	  class ITER_TUPLE, class OUT> inline OUT
copy(const zip_iterator<ITER_TUPLE>& ib,
     const zip_iterator<ITER_TUPLE>& ie, OUT out)
{
    return copy<T, ASSIGN, MASK>(ib.get_iterator_tuple(),
				 ie.get_iterator_tuple(), out);
}
    
}	// namespace simd
}	// namespace TU
#endif	// !TU_SIMD_TRANSFORM_H
