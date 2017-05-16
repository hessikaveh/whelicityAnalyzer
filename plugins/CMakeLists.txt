# Specify the minimum version for CMake

cmake_minimum_required(VERSION 2.8)
set (CMAKE_CXX_FLAGS "--std=gnu++11  ${CMAKE_CXX_FLAGS}")
# Project's name

project(whelicity)
include_directories("/cmsdata3/kaveh/workarea/cmssw/")
include_directories("/cmsdata3/kaveh/scratch/lhapdf-5.8.9/../local/include")
#include_directories("/home/hkaveh/root/include/")
include_directories("/opt/root-6.06.08/include/")

#include_directories("/home/hkaveh/cmsdata3/kaveh/workarea/cmssw/")

add_library(my MiniAnalyzer.cc ttamwtsolver.cpp MSETools.cc)

