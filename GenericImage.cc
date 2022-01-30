/*!
  \file		GenericImage.cc
  \author	Toshio UESHIBA
  \brief	任意の画素型をとれる画像クラスの実装
*/
#include "TU/Image++.h"

namespace TU
{
/************************************************************************
*  class GenericImage							*
************************************************************************/
//! 入力ストリームから画像の画素データを読み込む．
/*!
  \param in	入力ストリーム
  \return	inで指定した入力ストリーム
*/
std::istream&
GenericImage::restoreData(std::istream& in)
{
    using	std::begin;
    using	std::size;
    
    _colormap.resize(_format.ncolors());
    _colormap.restore(in);

    const auto	npads = _format.nbytesForPadding(width());
    if (_format.bottomToTop())
    {
	for (auto row = _a.rbegin(); row != _a.rend(); ++row)
	    if (!in.read(begin(*row), size(*row)) ||
		!in.ignore(npads))
		break;
    }
    else
    {
	for (auto row : _a)
	    if (!in.read(begin(row), size(row)) ||
		!in.ignore(npads))
		break;
    }

    return in;
}

//! 出力ストリームに画像の画素データを書き出す．
/*!
  \param out	出力ストリーム
  \return	outで指定した出力ストリーム
*/
std::ostream&
GenericImage::saveData(std::ostream& out) const
{
    using	std::cbegin;
    using	std::size;
    
    if (_colormap.size() > 0)
    {
	_colormap.save(out);
	for (size_t i = _colormap.size(); i < 256; ++i)
	    out.put(0).put(0).put(0).put(0);
    }
    
    Array<u_char>	pads(_format.nbytesForPadding(width()));
    if (_format.bottomToTop())
    {
	for (auto row = _a.rbegin(); row != _a.rend(); ++row)
	    if (!out.write(cbegin(*row), size(*row)) ||
		!pads.save(out))
		break;
    }
    else
    {
	for (const auto row : _a)
	    if (!out.write(cbegin(row), size(row)) ||
		!pads.save(out))
		break;
    }
    
    return out;
}

}
