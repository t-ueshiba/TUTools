/*!
  \file		DP.h
  \author	Toshio UESHIBA
  \brief	動的計画法による最適化器の定義と実装
*/
#ifndef TU_DP_H
#define TU_DP_H

#include "TU/Array++.h"

namespace TU
{
/************************************************************************
*  class DP<DOM, T>							*
************************************************************************/
//! 動的計画法によって評価関数の最小化を行うクラス
/*!
  \param DOM	評価関数の各ステージにおける引数の集合を指す反復子の型
  \param T	評価関数の値の型
*/
template <class DOM, class T>
class DP
{
  public:
    using argument_iterator = typename std::iter_value_t<DOM>::const_iterator;
    using argument_type	    = std::iter_value_t<argument_iterator>;
    using value_type	    = T;

  private:
    struct Node
    {
	Node()	:arg(), val(0), prev(0)			{}

	bool		operator <(const Node& node) const
			{
			    return val < node.val;
			}
	std::ostream&	put(std::ostream& out) const
			{
			    return out << *arg << ':' << val;
			}
	
	argument_iterator	arg;
	value_type		val;
	const Node*		prev;
    };

    typedef Array<Node>					node_array;
    typedef Array<node_array>				stage_array;
    typedef typename node_array::iterator		node_iterator;
    typedef typename node_array::const_iterator		const_node_iterator;
    typedef typename stage_array::iterator		stage_iterator;
    typedef typename stage_array::const_iterator	const_stage_iterator;
    
  public:
    explicit DP(value_type lambda=1)					;

    value_type		lambda()				const	;
    DP&			setLambda(value_type lambda)			;
    DP&			initialize(DOM args, DOM args_end)		;
    template <class FITER, class OUT>
    value_type		operator ()(FITER f, OUT out)			;
    std::ostream&	put(std::ostream& out)			const	;
    
  private:
    value_type		_lambda;	//!< 平滑化項の係数
    stage_array		_stages;
};

//! 動的計画法による最適化器を生成する．
/*!
  \param lambda	平滑化項の係数
*/
template <class DOM, class T> inline
DP<DOM, T>::DP(value_type lambda)
    :_lambda(lambda), _stages()
{
}

//! 平滑化項の係数を返す．
/*!
  \return	平滑化項の係数
*/ 
template <class DOM, class T> inline typename DP<DOM, T>::value_type
DP<DOM, T>::lambda() const
{
    return _lambda;
}

//! 平滑化項の係数をセットする．
/*!
  \param lambda	平滑化項の係数
  \return	この最適化器
*/ 
template <class DOM, class T> inline DP<DOM, T>&
DP<DOM, T>::setLambda(value_type lambda)
{
    _lambda = lambda;
    return *this;
}

//! 最適化する関数の定義域をセットする．
/*!
  \param args	最初のステージの引数が取り得る値の集合を指す反復子
  \param args	末尾の次のステージの引数が取り得る値の集合を指す反復子
  \return	この最適化器
*/ 
template <class DOM, class T> DP<DOM, T>&
DP<DOM, T>::initialize(DOM args, DOM args_end)
{
    _stages.resize(std::distance(args, args_end));
    
    for (stage_iterator stage = _stages.begin(); stage != _stages.end();
	 ++stage)
    {
      // 第iステージの変数 x_i がとり得る値の数だけノードを確保する．
	stage->resize(std::distance(args->begin(), args->end()));
	if (stage->size() == 0)
	    throw std::invalid_argument("DP<DOM, T>::initialize(): domain of each stage cannot be empty");
	
      // x_i がとり得るそれぞれの値をノードにセットする．
	argument_iterator	arg = args->begin();
	for (node_iterator node = stage->begin(); node != stage->end(); ++node)
	{
	    node->arg = arg;	// x_i がとり得る値を指す反復子をセット
	    ++arg;
	}
	++args;			// x_{i+1} がとり得る値の集合に移る
    }

    return *this;
}

//! 与えられた関数列を最適化する引数列を求める．
/*!
  \param f	最初のステージを最適化する関数を指す反復子
  \param out	最適値を与える引数列を出力先を指す反復子．末尾のステージの引数から
		先に出力されるので，逆反復子を与えること．
  \return	最適化された最小値
*/ 
template <class DOM, class T> template <class FITER, class OUT>
typename DP<DOM, T>::value_type
DP<DOM, T>::operator ()(FITER f, OUT out)
{
    if (_stages.size() == 0)
	return 0;

  // 最初のステージの変数x_iがとり得る値それぞれについてデータ項の値を求める．
    stage_iterator	stage = _stages.begin();
    for (node_iterator node = stage->begin(); node != stage->end(); ++node)
	node->val = (*f)(*node->arg);

    stage_iterator	stage_p = stage;	// 前ステージ
    while (++stage != _stages.end())		// 2番目以降のステージについて...
    {
	++f;				// 評価関数を次のステージに進める．
	
      // 現ステージの変数x_iがとり得る値それぞれについて...
	for (node_iterator node = stage->begin(); node != stage->end(); ++node)
	{
	  // 前ステージの変数x_{i-1}をどの値にすれば現ステージまでの
	  // 評価関数が最小化されるか調べる．
	    node_iterator	node_p = stage_p->begin();

	    node->val  = node_p->val + _lambda * (*f)(*node_p->arg, *node->arg);
	    node->prev = node_p;

	    while (++node_p != stage_p->end())
	    {
		value_type	val = node_p->val
				    + _lambda * (*f)(*node_p->arg, *node->arg);
		
		if (val < node->val)
		{
		    node->val  = val;
		    node->prev = node_p;
		}
	    }

	    node->val += (*f)(*node->arg);
	}
	
	++stage_p;
    }

    --stage;					// 最終ステージに戻る
    const_node_iterator	node = std::min_element(stage->begin(), stage->end());
    value_type		minval = node->val;
    do
    {
	*out = *node->arg;
	++out;
    } while ((node = node->prev) != 0);

    return minval;
}

//! 各ステージが取り得る各引数に対して，その値とその値に対する関数の最小値を出力する．
/*!
  \param out	出力ストリーム
  \return	outで指定した出力ストリーム
*/ 
template <class DOM, class T> std::ostream&
DP<DOM, T>::put(std::ostream& out) const
{
    for (const_stage_iterator stage  = _stages.begin();
			      stage != _stages.end(); ++stage)
    {
	for (const_node_iterator node  = stage->begin();
				 node != stage->end(); ++node)
	{
	    out << ' ';
	    node->put(out);
	}
	out << std::endl;
    }

    return out << std::endl;
}

}
#endif	// !TU_DP_H
