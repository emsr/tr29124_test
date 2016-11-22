# -*- C++ -*-
# Integration utilities for the C++ library testsuite.
#
# Copyright (C) 2011, 2012 Free Software Foundation, Inc.
#
# This file is part of the GNU ISO C++ Library.  This library is free
# software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this library; see the file COPYING3.  If not see
# <http:#www.gnu.org/licenses/>.
#
#Written by Jason Dick

cmake_minimum_required(VERSION 2.8)

project(hermite_test)

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS " -g -Wall -Wno-missing-braces -std=c++17")
endif()

add_executable(hermite_test hermite_test.cpp)

install(TARGETS hermite_test DESTINATION bin)
