cmake_minimum_required(VERSION 2.8.3)
project(itia_controllers_and_filters)
add_compile_options(-std=c++11 -funroll-loops -Wall -Ofast)
set(CMAKE_BUILD_TYPE Release)


find_package(	catkin REQUIRED COMPONENTS
							roscpp
							cmake_modules
							std_msgs
							itia_msgs
							itia_futils
							itia_rutils
)

find_package(Eigen3 REQUIRED COMPONENTS Dense)
find_package(orocos_kdl REQUIRED )

catkin_package(
								INCLUDE_DIRS include
								LIBRARIES impedance_control
								CATKIN_DEPENDS roscpp cmake_modules std_msgs itia_msgs itia_futils itia_rutils
								DEPENDS orocos_kdl Eigen3
)

include_directories(
${catkin_INCLUDE_DIRS}
${Eigen3_INCLUDE_DIRS}
${orocos_kdl_INCLUDE_DIRS}
include
)


add_library(impedance_control
						src/${PROJECT_NAME}/c_impedance_control.cpp
)
add_dependencies(impedance_control ${catkin_EXPORTED_TARGETS})


add_executable(test_controller_node src/test_controller.cpp)
add_dependencies(test_controller_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(test_controller_node
											${catkin_LIBRARIES}
											${orocos_kdl_LIBRARIES}
											${Eigen3_LIBRARIES}
)
