cmake_minimum_required(VERSION 3.9)
set(PROJECT_NAME "DelaunuQuad")
project(${PROJECT_NAME} LANGUAGES CXX)

include_directories("include")

file(GLOB_RECURSE SRC_FILES "src/*.*")

add_executable(${PROJECT_NAME} ${INC_FILES} ${SRC_FILES})
set_property(TARGET ${PROJECT_NAME}  PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_NAME}  PROPERTY CXX_STANDARD_REQUIRED TRUE)
add_definitions(-DNOMINMAX)

include_directories("input/")

#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
