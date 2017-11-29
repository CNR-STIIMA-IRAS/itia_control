cmake_minimum_required(VERSION 2.8.3)
project(control_models)
add_compile_options(-std=c++11 -funroll-loops -Wall -Ofast)
set(CMAKE_BUILD_TYPE Release)


## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  itia_futils
  itia_mutils
  itia_rutils
  kdl_conversions
  roscpp
  std_msgs
  tf_conversions
)
find_package(orocos_kdl REQUIRED)
find_package(Eigen3 REQUIRED COMPONENTS)


catkin_package(
								INCLUDE_DIRS include
								LIBRARIES control_models
								CATKIN_DEPENDS eigen_conversions itia_futils itia_mutils itia_rutils kdl_conversions roscpp std_msgs tf_conversions
								#  DEPENDS 
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen3_INCLUDE_DIRS}
)

add_library(control_models
  src/${PROJECT_NAME}/ss.cpp
)
add_dependencies(control_models ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_executable(control_models_node src/control_model_test.cpp)
add_dependencies(control_models_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(control_models_node
											control_models
											${catkin_LIBRARIES}
											${Eigen3_LIBRARIES}
)
