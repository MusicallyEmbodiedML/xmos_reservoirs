
cmake_minimum_required(VERSION 3.28)

cmake_policy(SET CMP0065 NEW) #lose -rdynamic

message("${CMAKE_EXE_LINKER_FLAGS_INIT}")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_COMPILER "xcc")
set(CMAKE_C_COMPILER_FORCED TRUE)
set(CMAKE_CXX_COMPILER "xcc")
set(CMAKE_CXX_COMPILER_FORCED TRUE)

project(XRES C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)


set(CMAKE_C_FLAGS "-Os -target=XK-EVK-XU316 -fdata-sections -ffunction-sections")
set(CMAKE_CXX_FLAGS -target=XK-EVK-XU316)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN -Os -report -fdata-sections -ffunction-sections")
SET(CMAKE_EXE_LINKER_FLAGS "-Wl --gc-sections")

project(XRES VERSION 0.1)
set(XRES 0.1)

# Collect all source files
#set(XRES_SOURCES ${XRES_SOURCES}  PARENT_SCOPE)

file(GLOB_RECURSE XRES_SOURCES main.cpp 
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/LinearAlgebra/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/SignalProcessing/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Miscellaneous/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Statistics/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/MachineLearning/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Hardware/*.c
#${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Optimization/*.c
)

#list(FILTER XRES_SOURCES EXCLUDE REGEX ".*quadprog.c$")

# Collect all library paths

set(XRES_HEADERS
#    ${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Headers
#    ${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/LinearAlgebra
#    ${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/SignalProcessing
#    ${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Statistics
#    ${CMAKE_CURRENT_LIST_DIR}/CControl/src/CControl/Sources/Optimization
)


add_executable(XRES ${XRES_SOURCES} ${XRES_HEADERS})

