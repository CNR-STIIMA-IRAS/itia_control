cmake_minimum_required(VERSION 2.8.3)
project(itia_clik)
add_compile_options(-std=c++11 -funroll-loops -Wall -Ofast)
set(CMAKE_BUILD_TYPE Release)


find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  geometry_msgs
  sensor_msgs
  itia_mutils
  kdl_conversions
  roscpp
  tf_conversions
  kdl_parser
  itia_rutils
  dynamic_reconfigure 
  urdf
  rospy
)

find_package(orocos_kdl REQUIRED COMPONENTS system)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES itia_clik
  CATKIN_DEPENDS eigen_conversions geometry_msgs itia_mutils kdl_conversions roscpp tf_conversions kdl_parser itia_rutils sensor_msgs dynamic_reconfigure  rospy urdf
  DEPENDS  orocos_kdl 
)


include_directories(include
										${catkin_INCLUDE_DIRS}
										${orocos_kdl_INCLUDE_DIRS}
)

add_library(itia_clik
						src/${PROJECT_NAME}/itia_clik.cpp
)
add_dependencies(itia_clik 
								${${PROJECT_NAME}_EXPORTED_TARGETS} 
								${catkin_EXPORTED_TARGETS})


add_executable(itia_clik_test src/itia_clik_test.cpp)
add_dependencies(itia_clik_test itia_clik ) 
target_link_libraries(itia_clik_test
											${orocos_kdl_LIBRARIES}
											${catkin_LIBRARIES}
											itia_clik
)
