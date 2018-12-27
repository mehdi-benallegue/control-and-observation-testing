# Copyright 2013, CNRS
#
# Author: Mehdi Benallegue
#
# This file is part of state-observation.
# state-observation is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# state-observation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# state-observation. If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)



SET(PROJECT_NAME control-and-observation-testing)
SET(PROJECT_DESCRIPTION
  "Prototyping of simple controllers")
SET(PROJECT_URL "")

# Search for dependencies.
ADD_REQUIRED_DEPENDENCY("state-observation")
ADD_REQUIRED_DEPENDENCY("control-prototyping")
ADD_REQUIRED_DEPENDENCY("pinocchio")
ADD_REQUIRED_DEPENDENCY("yaml-cpp")
ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.0.0")

SET(CUSTOM_HEADER_DIR "${PROJECT_NAME}")

SET(${PROJECT_NAME}_HEADERS
  #include/control-and-observation-testing/file.hpp
)
SET(CXX_DISABLE_WERROR TRUE)

SETUP_PROJECT()

SET(BOOST_COMPONENTS filesystem system thread program_options unit_test_framework timer)

SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()

PKG_CONFIG_APPEND_LIBS(${PROJECT_NAME})

ADD_SUBDIRECTORY(src)
SETUP_PROJECT_FINALIZE()
