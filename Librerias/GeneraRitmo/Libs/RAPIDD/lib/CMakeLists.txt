cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_FLAGS_INIT "-Werror -Wall -Wextra -Wpedantic")

list(APPEND CMAKE_REQUIRED_LIBRARIES m)

# set the project name
project(RAPIDD VERSION 1.0)
configure_file(RAPIDDConfig.h.in RAPIDDConfig.h)

add_library(RAPIDD SHARED pcubature.c hcubature.c source/halo.c source/effFormFact.c source/difCrosSec.c source/difRateGen.c source/WriteData.c source/coeffs_eft.c source/MenenFF.c source/binning_general.c source/difRateGen.c source/WriteData.c source/phys_consts.c source/efficiency_curve.c source/specific_exps.c)

find_package(GSL REQUIRED)
target_link_libraries(RAPIDD PUBLIC GSL::gsl GSL::gslcblas)
# add the executable
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(rpd main.c)
target_link_libraries(rpd PUBLIC RAPIDD GSL::gsl GSL::gslcblas m)

add_executable(rpd3 main.c)
target_link_libraries(rpd3 PUBLIC RAPIDD GSL::gsl GSL::gslcblas m)

#if(OpenMP_CXX_FOUND)
    #set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    #target_link_libraries(rpd PUBLIC OpenMP::OpenMP_CXX)
#endif()



#target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
