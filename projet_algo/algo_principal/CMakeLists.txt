cmake_minimum_required(VERSION 3.8)
project(algo_principal)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        fonctions.c
        fonctions.h
        fonctions_masques_kmers.c
        fonctions_PSSM.c
        fonctions_sequences.c
        main.c)

add_executable(algo_principal ${SOURCE_FILES})