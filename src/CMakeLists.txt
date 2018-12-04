# Copyright 2008, 2009, 2010, 2013 Olivier Stasse, Mehdi Benallegue JRL, CNRS/AIST
#
# This file is part of control-prototyping-test.
# state-observation is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# state-observation is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with state-observation.  If not, see <http://www.gnu.org/licenses/>.

# prefix and suffix each element of list by ${prefix}elemnt${suffix}
macro(ADDPREFIX newlist prefix list_name)
    # create empty list - necessary?
    SET(${newlist})

     # prefix and suffix elements
    foreach(l ${${list_name}})
      list(APPEND ${newlist} ${prefix}${l} )
    endforeach()

endmacro(ADDPREFIX)

if(${_Eigen_VERSION} VERSION_LESS "3.2.0")
    add_definitions(-DEIGEN_VERSION_LESS_THAN_3_2)
endif()



ADDPREFIX(${PROJECT_NAME}_ABSOLUTE_HEADERS  "${CMAKE_SOURCE_DIR}/"  ${PROJECT_NAME}_HEADERS)

SET(LIBRARY_NAME ${PROJECT_NAME})



ADD_EXECUTABLE(${LIBRARY_NAME}
    integral-passivity-based-tests.cpp
  ${${PROJECT_NAME}_ABSOLUTE_HEADERS}
  )
  
  
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} state-observation)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} pinocchio)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} control-prototyping)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} yaml-cpp)

SET_TARGET_PROPERTIES(${LIBRARY_NAME}
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

INSTALL(TARGETS ${LIBRARY_NAME}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

TARGET_LINK_LIBRARIES(${LIBRARY_NAME} rt)
ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.0.0")
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})
SET(LINEAR_ALGEBRA_LIB_FOUND 1)
SET(_EIGEN_MATRIX_FOUND_ 1)


