/*!
 *  \file	main.cc
 *  \author	Toshio UESHIBA
 *  \brief	Test program for class Quaternion<T>
 */
#include "TU/Quaternion.h"
#include <boost/math/quaternion.hpp>

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
	const auto	degree = 180 / M_PI;
	
	std::cerr << "q >> ";
	char	c;
	if (!(std::cin >> c))
	    break;

	std::cin.putback(c);
	Quaternion<T>	q;
	q.get(std::cin);

	std::cout << "q              = " << q;
	std::cout << "inverse(q)     = " << inverse(q);
	std::cout << "q * inverse(q) = " << q * inverse(q);
	std::cout << "inverse(q) * q = " << inverse(q) * q;

	q.normalize();
	std::cout << "\n--- R ---\n"  << q.R();
	std::cout << "--- Rt ---\n"   << q.Rt();
	std::cout << "--- Rt*R ---\n" << q.Rt() * q.R();

	std::cout << "q           = " << q;
	std::cout << "q from R    = " << Quaternion<T>(q.R());

	const auto	rpy = degree * q.rpy();
	std::cout << "rpy         = " << rpy;

	const Quaternion<T>	p(1, 2, 3, 4), r(5, 6, 7, 8);
	std::cout << 0.4*p + 0.6*r;
    }
}

template <class T> void
doJob1()
{
    for (;;)
    {
	std::cerr << "q >> ";
	char	c;
	if (!(std::cin >> c))
	    break;

	std::cin.putback(c);
	boost::math::quaternion<T>	q;
	std::cin >> q;

	std::cout << "q              = " << q << std::endl;
	std::cout << "inverse(q)     = " << inverse(q) << std::endl;
	std::cout << "q * inverse(q) = " << q * inverse(q) << std::endl;
	std::cout << "inverse(q) * q = " << inverse(q) * q << std::endl;
	std::cout << "conj(q)        = " << conj(q) << std::endl;
    }
}
    
}

int
main(int argc, char* argv[])
{
    try
    {
	using value_type = double;
	
      //TU::doJob<value_type>();
	TU::doJob1<value_type>();
    }
    catch (const std::exception& err)
    {
	std::cerr << err.what() << std::endl;
	return 1;
    }

    return 0;
}
