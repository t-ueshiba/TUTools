/*!
 *  \file	main.cc
 *  \author	Toshio UESHIBA
 *  \brief	Test program for class DualNumber<T>
 */
#include <boost/math/quaternion.hpp>
#include "TU/Quaternion.h"
#include "TU/DualNumber.h"

namespace boost
{
namespace math
{
template <class T> inline boost::math::quaternion<T>
inverse(const boost::math::quaternion<T>& x)
{
    return conj(x) / norm(x);
}

}
}

namespace TU
{
template <class T> void
doJob()
{
    for (;;)
    {
      //using	q_type = Quaternion<DualNumber<T> >;
	using	q_type = DualNumber<Quaternion<T> >;
      //using	q_type = boost::math::quaternion<DualNumber<T> >;
      //using	q_type = DualNumber<boost::math::quaternion<T> >;
	
	std::cerr << "q >> ";
	char	c;
	if (!(std::cin >> c))
	    break;

	std::cin.putback(c);
	q_type	q;
	std::cin >> q;

	std::cout << "q              = " << q << std::endl;
	std::cout << "inverse(q)     = " << inverse(q) << std::endl;
	std::cout << "q * inverse(q) = " << q * inverse(q) << std::endl;
	std::cout << "inverse(q) * q = " << inverse(q) * q << std::endl;
    }
}
    
}

int
main(int argc, char* argv[])
{
    try
    {
	using value_type = double;
	
	TU::doJob<value_type>();
    }
    catch (const std::exception& err)
    {
	std::cerr << err.what() << std::endl;
	return 1;
    }

    return 0;
}
