#
#  $Id$
#
#################################
#  User customizable macros	#
#################################
PROGRAM		= $(shell basename $(PWD))
#LIBRARY		= lib$(shell basename $(PWD))

VPATH		=

IDLS		=
MOCHDRS		=

INCDIRS		= -I../../..
CPPFLAGS	= -DNDEBUG -DSCORE_ARRAY3 #-DTU_DEBUG #-DUSE_TBB
CFLAGS		= -g -O3
NVCCFLAGS	= -g
SANITIZE	= -fsanitize=address -fsanitize-address-use-after-scope
ifeq ($(shell arch), armv7l)
  CPPFLAGS     += -DNEON
else ifeq ($(shell arch), aarch64)
  CPPFLAGS     += -DNEON
else
  CPPFLAGS     += -DSSE4
  CFLAGS       += -msse4
#  CPPFLAGS     += -DAVX2
#  CFLAGS       += -mavx2
endif
CCFLAGS		= $(CFLAGS) #${SANITIZE} #-ftemplate-backtrace-limit=0

LIBS		= -lTUTools++ -ltbb
LINKER		= $(CXX) #${SANITIZE}

BINDIR		= $(PREFIX)/bin
LIBDIR		= $(PREFIX)/lib
INCDIR		= $(PREFIX)/include

#########################
#  Macros set by mkmf	#
#########################
SUFFIX		= .cc:sC .cpp:sC .cu:sC
EXTHDRS		= ../../../TU/Array++.h \
		../../../TU/BoxFilter.h \
		../../../TU/Camera++.h \
		../../../TU/Filter2.h \
		../../../TU/Geometry++.h \
		../../../TU/GuidedFilter.h \
		../../../TU/Image++.h \
		../../../TU/Manip.h \
		../../../TU/Minimize.h \
		../../../TU/Profiler.h \
		../../../TU/Rectify.h \
		../../../TU/StereoBase.h \
		../../../TU/StereoUtility.h \
		../../../TU/TreeFilter.h \
		../../../TU/Vector++.h \
		../../../TU/Warp.h \
		../../../TU/algorithm.h \
		../../../TU/functional.h \
		../../../TU/io.h \
		../../../TU/iterator.h \
		../../../TU/pair.h \
		../../../TU/range.h \
		../../../TU/simd/Array++.h \
		../../../TU/simd/allocator.h \
		../../../TU/simd/arithmetic.h \
		../../../TU/simd/arm/arch.h \
		../../../TU/simd/arm/arithmetic.h \
		../../../TU/simd/arm/bit_shift.h \
		../../../TU/simd/arm/cast.h \
		../../../TU/simd/arm/compare.h \
		../../../TU/simd/arm/cvt.h \
		../../../TU/simd/arm/dup.h \
		../../../TU/simd/arm/insert_extract.h \
		../../../TU/simd/arm/load_store.h \
		../../../TU/simd/arm/logical.h \
		../../../TU/simd/arm/lookup.h \
		../../../TU/simd/arm/select.h \
		../../../TU/simd/arm/shift.h \
		../../../TU/simd/arm/type_traits.h \
		../../../TU/simd/arm/vec.h \
		../../../TU/simd/arm/zero.h \
		../../../TU/simd/bit_shift.h \
		../../../TU/simd/cast.h \
		../../../TU/simd/compare.h \
		../../../TU/simd/config.h \
		../../../TU/simd/cvt.h \
		../../../TU/simd/cvtdown_iterator.h \
		../../../TU/simd/cvtup_iterator.h \
		../../../TU/simd/dup.h \
		../../../TU/simd/insert_extract.h \
		../../../TU/simd/load_store.h \
		../../../TU/simd/load_store_iterator.h \
		../../../TU/simd/logical.h \
		../../../TU/simd/lookup.h \
		../../../TU/simd/misc.h \
		../../../TU/simd/select.h \
		../../../TU/simd/shift.h \
		../../../TU/simd/shift_iterator.h \
		../../../TU/simd/simd.h \
		../../../TU/simd/transform.h \
		../../../TU/simd/type_traits.h \
		../../../TU/simd/vec.h \
		../../../TU/simd/x86/arch.h \
		../../../TU/simd/x86/arithmetic.h \
		../../../TU/simd/x86/bit_shift.h \
		../../../TU/simd/x86/cast.h \
		../../../TU/simd/x86/compare.h \
		../../../TU/simd/x86/cvt.h \
		../../../TU/simd/x86/dup.h \
		../../../TU/simd/x86/insert_extract.h \
		../../../TU/simd/x86/load_store.h \
		../../../TU/simd/x86/logical.h \
		../../../TU/simd/x86/logical_base.h \
		../../../TU/simd/x86/lookup.h \
		../../../TU/simd/x86/select.h \
		../../../TU/simd/x86/shift.h \
		../../../TU/simd/x86/shuffle.h \
		../../../TU/simd/x86/svml.h \
		../../../TU/simd/x86/type_traits.h \
		../../../TU/simd/x86/unpack.h \
		../../../TU/simd/x86/vec.h \
		../../../TU/simd/x86/zero.h \
		../../../TU/simd/zero.h \
		../../../TU/tuple.h \
		../../../TU/type_traits.h
HDRS		=
SRCS		= Rectify.cc \
		main.cc
OBJS		= Rectify.o \
		main.o

#include $(PROJECT)/lib/rtc.mk		# IDLHDRS, IDLSRCS, CPPFLAGS, OBJS, LIBS
#include $(PROJECT)/lib/qt.mk		# MOCSRCS, OBJS
#include $(PROJECT)/lib/cnoid.mk	# CPPFLAGS, LIBS, LIBDIR
#include $(PROJECT)/lib/lib.mk		# PUBHDRS TARGHDRS
include $(PROJECT)/lib/common.mk
###
Rectify.o: ../../../TU/Rectify.h ../../../TU/Warp.h \
	../../../TU/simd/Array++.h ../../../TU/simd/simd.h \
	../../../TU/simd/config.h ../../../TU/simd/vec.h \
	../../../TU/simd/type_traits.h ../../../TU/simd/x86/type_traits.h \
	../../../TU/simd/arm/type_traits.h ../../../TU/simd/x86/vec.h \
	../../../TU/simd/x86/arch.h ../../../TU/simd/arm/vec.h \
	../../../TU/simd/arm/arch.h ../../../TU/simd/allocator.h \
	../../../TU/simd/load_store.h ../../../TU/simd/x86/load_store.h \
	../../../TU/simd/arm/load_store.h ../../../TU/simd/zero.h \
	../../../TU/simd/x86/zero.h ../../../TU/simd/arm/zero.h \
	../../../TU/simd/cast.h ../../../TU/simd/x86/cast.h \
	../../../TU/simd/arm/cast.h ../../../TU/simd/insert_extract.h \
	../../../TU/simd/x86/insert_extract.h \
	../../../TU/simd/arm/insert_extract.h ../../../TU/simd/shift.h \
	../../../TU/tuple.h ../../../TU/type_traits.h \
	../../../TU/simd/x86/shift.h ../../../TU/simd/arm/shift.h \
	../../../TU/simd/bit_shift.h ../../../TU/simd/x86/bit_shift.h \
	../../../TU/simd/arm/bit_shift.h ../../../TU/simd/dup.h \
	../../../TU/simd/cvt.h ../../../TU/iterator.h \
	../../../TU/simd/x86/cvt.h ../../../TU/simd/x86/unpack.h \
	../../../TU/simd/arm/cvt.h ../../../TU/simd/logical.h \
	../../../TU/simd/x86/logical.h ../../../TU/simd/x86/logical_base.h \
	../../../TU/simd/arm/logical.h ../../../TU/simd/x86/dup.h \
	../../../TU/simd/arm/dup.h ../../../TU/simd/compare.h \
	../../../TU/simd/x86/compare.h ../../../TU/simd/arm/compare.h \
	../../../TU/simd/select.h ../../../TU/simd/x86/select.h \
	../../../TU/simd/arm/select.h ../../../TU/simd/arithmetic.h \
	../../../TU/simd/x86/arithmetic.h ../../../TU/simd/arm/arithmetic.h \
	../../../TU/simd/misc.h ../../../TU/simd/x86/shuffle.h \
	../../../TU/simd/x86/svml.h ../../../TU/simd/transform.h \
	../../../TU/functional.h ../../../TU/simd/lookup.h \
	../../../TU/simd/x86/lookup.h ../../../TU/simd/arm/lookup.h \
	../../../TU/simd/load_store_iterator.h \
	../../../TU/simd/cvtdown_iterator.h ../../../TU/simd/cvtup_iterator.h \
	../../../TU/simd/shift_iterator.h ../../../TU/Array++.h \
	../../../TU/range.h ../../../TU/algorithm.h ../../../TU/Image++.h \
	../../../TU/pair.h ../../../TU/Manip.h ../../../TU/Camera++.h \
	../../../TU/Geometry++.h ../../../TU/Minimize.h \
	../../../TU/Vector++.h
main.o: ../../../TU/algorithm.h ../../../TU/io.h ../../../TU/Rectify.h \
	../../../TU/Warp.h ../../../TU/simd/Array++.h ../../../TU/simd/simd.h \
	../../../TU/simd/config.h ../../../TU/simd/vec.h \
	../../../TU/simd/type_traits.h ../../../TU/simd/x86/type_traits.h \
	../../../TU/simd/arm/type_traits.h ../../../TU/simd/x86/vec.h \
	../../../TU/simd/x86/arch.h ../../../TU/simd/arm/vec.h \
	../../../TU/simd/arm/arch.h ../../../TU/simd/allocator.h \
	../../../TU/simd/load_store.h ../../../TU/simd/x86/load_store.h \
	../../../TU/simd/arm/load_store.h ../../../TU/simd/zero.h \
	../../../TU/simd/x86/zero.h ../../../TU/simd/arm/zero.h \
	../../../TU/simd/cast.h ../../../TU/simd/x86/cast.h \
	../../../TU/simd/arm/cast.h ../../../TU/simd/insert_extract.h \
	../../../TU/simd/x86/insert_extract.h \
	../../../TU/simd/arm/insert_extract.h ../../../TU/simd/shift.h \
	../../../TU/tuple.h ../../../TU/type_traits.h \
	../../../TU/simd/x86/shift.h ../../../TU/simd/arm/shift.h \
	../../../TU/simd/bit_shift.h ../../../TU/simd/x86/bit_shift.h \
	../../../TU/simd/arm/bit_shift.h ../../../TU/simd/dup.h \
	../../../TU/simd/cvt.h ../../../TU/iterator.h \
	../../../TU/simd/x86/cvt.h ../../../TU/simd/x86/unpack.h \
	../../../TU/simd/arm/cvt.h ../../../TU/simd/logical.h \
	../../../TU/simd/x86/logical.h ../../../TU/simd/x86/logical_base.h \
	../../../TU/simd/arm/logical.h ../../../TU/simd/x86/dup.h \
	../../../TU/simd/arm/dup.h ../../../TU/simd/compare.h \
	../../../TU/simd/x86/compare.h ../../../TU/simd/arm/compare.h \
	../../../TU/simd/select.h ../../../TU/simd/x86/select.h \
	../../../TU/simd/arm/select.h ../../../TU/simd/arithmetic.h \
	../../../TU/simd/x86/arithmetic.h ../../../TU/simd/arm/arithmetic.h \
	../../../TU/simd/misc.h ../../../TU/simd/x86/shuffle.h \
	../../../TU/simd/x86/svml.h ../../../TU/simd/transform.h \
	../../../TU/functional.h ../../../TU/simd/lookup.h \
	../../../TU/simd/x86/lookup.h ../../../TU/simd/arm/lookup.h \
	../../../TU/simd/load_store_iterator.h \
	../../../TU/simd/cvtdown_iterator.h ../../../TU/simd/cvtup_iterator.h \
	../../../TU/simd/shift_iterator.h ../../../TU/Array++.h \
	../../../TU/range.h ../../../TU/Image++.h ../../../TU/pair.h \
	../../../TU/Manip.h ../../../TU/Camera++.h ../../../TU/Geometry++.h \
	../../../TU/Minimize.h ../../../TU/Vector++.h \
	../../../TU/GuidedFilter.h ../../../TU/BoxFilter.h \
	../../../TU/Filter2.h ../../../TU/TreeFilter.h ../../../TU/Profiler.h \
	../../../TU/StereoUtility.h ../../../TU/StereoBase.h
