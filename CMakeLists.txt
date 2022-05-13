cmake_minimum_required(VERSION 3.11)
project(transformer)

## Find necessary packages
#find_package( Python3 REQUIRED COMPONENTS Interpreter Development)
##include_directories( ${PYTHON3_INCLUDE_DIRS} )

find_package(ROOT REQUIRED COMPONENTS RIO MathCore MathMore Matrix Tree Hist Gpad)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories( ./include )

add_library(AnalyticPhi23 SHARED src/AnalyticPhi23.cpp)
target_link_libraries(AnalyticPhi23 PUBLIC ROOT::MathMore)

add_library(CalDerivation SHARED src/CalDerivation.cpp)
target_link_libraries(CalDerivation PUBLIC AnalyticPhi23 ROOT::Matrix)

add_library(Base SHARED src/Base.cpp)
target_link_libraries(Base PUBLIC ROOT::MathMore CalDerivation)

#add_library(SetMatrix SHARED src/SetMatrix.cpp)
#target_link_libraries(SetMatrix PUBLIC CalDerivation ROOT::Matrix)

add_library(Contour SHARED src/Contour.cpp)
target_link_libraries(Contour PUBLIC Base)

add_library(Fitter SHARED src/Fitter.cpp)
target_link_libraries(Fitter PUBLIC ROOT::MathCore Base Contour)

add_executable(cont_test src/Contour_Test.cxx)
target_link_libraries(cont_test PUBLIC Fitter Base ROOT::RIO ROOT::MathCore ROOT::Tree ROOT::Hist ROOT::Gpad)

add_executable(chi2_test src/Chi2_Test.cxx)
#target_include_directories(${PROJECT_NAME} PUBLIC)
target_link_libraries(chi2_test Base ROOT::RIO ROOT::MathCore ROOT::Tree ROOT::Hist)

add_executable(newminos_test src/NewMinos_Test.cpp)
target_link_libraries(newminos_test Fitter Base ROOT::RIO ROOT::Tree)

#add_executable(transformer src/transformer.cpp)
#target_include_directories(${PROJECT_NAME} PRIVATE ${Python3_INCLUDE_DIRS})
#target_link_libraries(transformer SetMatrix ROOT::RIO ROOT::MathCore ROOT::Tree ROOT::Hist ${Python3_LIBRARIES})
