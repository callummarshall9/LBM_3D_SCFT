cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
project(cmdline_LBM)
include_directories("${CMAKE_SOURCE_DIR}/LBM" "${CMAKE_SOURCE_DIR}/LBM/headers" "${CMAKE_SOURCE_DIR}/lib/rapidjson/include/rapidjson" "${CMAKE_SOURCE_DIR}/lib/fast-cpp-csv-parser")
find_package (Threads)
# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/Catch2)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
add_executable(main_cmdline ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp SCFT/SCFT.cpp SCFT/SCFT.hpp)
target_sources(main_cmdline
PRIVATE
        ${CMAKE_SOURCE_DIR}/SCFT/qDaggerPropagator.cpp
        ${CMAKE_SOURCE_DIR}/SCFT/qDaggerPropagator.hpp
        ${CMAKE_SOURCE_DIR}/SCFT/qPropagator.cpp
        ${CMAKE_SOURCE_DIR}/SCFT/qPropagator.hpp
        ${CMAKE_SOURCE_DIR}/SCFT/SCFT.cpp
        ${CMAKE_SOURCE_DIR}/SCFT/SCFT.hpp
        ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp
        ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp
        ${CMAKE_SOURCE_DIR}/lib/rapidjson/include/rapidjson/document.h
        ${CMAKE_SOURCE_DIR}/lib/rapidjson/include/rapidjson/writer.h
        ${CMAKE_SOURCE_DIR}/lib/rapidjson/include/rapidjson/stringbuffer.h
        ${CMAKE_SOURCE_DIR}/lib/fast-cpp-csv-parser/csv.h
)
target_link_libraries(main_cmdline ${CMAKE_THREAD_LIBS_INIT} Catch)
