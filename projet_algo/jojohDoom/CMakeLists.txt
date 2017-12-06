cmake_minimum_required(VERSION 3.8)
project(jojohDoom)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        fonctions_j.c
        fonctions_j.h
        main_j.c)

add_executable(jojohDoom ${SOURCE_FILES})