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

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS}) 
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
ADDPREFIX(${PROJECT_NAME}_ABSOLUTE_HEADERS  "${CMAKE_SOURCE_DIR}/"  ${PROJECT_NAME}_HEADERS)


ADD_EXECUTABLE(integral-passivity-based-tests
    integral-passivity-based-tests.cpp
  ${${PROJECT_NAME}_ABSOLUTE_HEADERS}
  )
  
PKG_CONFIG_USE_DEPENDENCY(integral-passivity-based-tests state-observation)
PKG_CONFIG_USE_DEPENDENCY(integral-passivity-based-tests pinocchio)
PKG_CONFIG_USE_DEPENDENCY(integral-passivity-based-tests control-prototyping)
PKG_CONFIG_USE_DEPENDENCY(integral-passivity-based-tests yaml-cpp)



  
ADD_EXECUTABLE(model-based-flex-estimation-tester
    model-based-flex-estimation-tester.cpp
  ${${PROJECT_NAME}_ABSOLUTE_HEADERS}
  )


PKG_CONFIG_USE_DEPENDENCY(model-based-flex-estimation-tester state-observation)
PKG_CONFIG_USE_DEPENDENCY(model-based-flex-estimation-tester pinocchio)
PKG_CONFIG_USE_DEPENDENCY(model-based-flex-estimation-tester control-prototyping)
PKG_CONFIG_USE_DEPENDENCY(model-based-flex-estimation-tester yaml-cpp)


SET_TARGET_PROPERTIES(integral-passivity-based-tests 
#model-based-flex-estimation-tester
PROPERTIES SOVERSION ${PROJECT_VERSION}
            INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

INSTALL(TARGETS integral-passivity-based-tests 
#model-based-flex-estimation-tester  
DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)




