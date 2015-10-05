/*
 *  $Id$
 */
#if !defined(__TU_SIMD_INTEL_ZERO_H)
#define __TU_SIMD_INTEL_ZERO_H

namespace TU
{
namespace simd
{
#define SIMD_ZERO(type)							\
    SIMD_SPECIALIZED_FUNC(vec<type> zero<type>(),			\
			  setzero, (), void, type, SIMD_BASE)

SIMD_ZERO(int8_t)
SIMD_ZERO(int16_t)
SIMD_ZERO(int32_t)
SIMD_ZERO(int64_t)
SIMD_ZERO(u_int8_t)
SIMD_ZERO(u_int16_t)
SIMD_ZERO(u_int32_t)
SIMD_ZERO(u_int64_t)
    
#if defined(SSE)
  SIMD_ZERO(float)
#endif
#if defined(SSE2)
  SIMD_ZERO(double)
#endif

#undef SIMD_ZERO

}	// namespace simd
}	// namespace TU
#endif	// !__TU_SIMD_INTEL_ZERO_H
