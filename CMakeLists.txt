cmake_minimum_required(VERSION 3.5)

# project name
project (HPCOT)

 

add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11")

include_directories(.)
INCLUDE_DIRECTORIES("/usr/include/eigen3")

find_package(Eigen3 REQUIRED)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")

set(CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA REQUIRED)

INCLUDE_DIRECTORIES(./${CUDA_INCLUDE_DIRS})
include_directories(.)
set(CUDA_NVCC_FLAGS  --fmad=true -Xptxas -dlcm=ca -O3
    -gencode=arch=compute_35,code=sm_35
   -gencode=arch=compute_50,code=sm_50
    -gencode=arch=compute_60,code=sm_60
   -gencode=arch=compute_70,code=sm_70
    -gencode=arch=compute_75,code=sm_75
   -gencode=arch=compute_75,code=compute_75 )

set(CUDA_NVCC_FLAGS
    --relocatable-device-code=true -Xcompiler
    -fPIC )

# set your path here
LINK_DIRECTORIES("/home")


INCLUDE_DIRECTORIES (LogSinkhorn)
#INCLUDE_DIRECTORIES (LP_Lemon)
INCLUDE_DIRECTORIES (Sinkhorn)
INCLUDE_DIRECTORIES (Common)
INCLUDE_DIRECTORIES (ShortCutSolver)


add_subdirectory(Common)
add_subdirectory(Sinkhorn)
#add_subdirectory(LP_Lemon)  // if the lemon is not support 
add_subdirectory(ShortCutSolver)
add_subdirectory(SparseSinkhorn)
add_subdirectory(LogSinkhorn)
#add_subdirectory(LP_CPLEX)   // if the cplex is not support 



add_executable(Multiscale  Multiscale.cpp  )
target_include_directories(Multiscale PUBLIC ${SINKHORN_INCLUDE_DIRECTORY})
target_compile_options(Multiscale PUBLIC ${SINKHORN_COMPILE_OPTIONS})
target_link_libraries(Multiscale
        Common
        Sinkhorn
        )


# add_executable(ShortCut_Lemon
#         ShortCut.cpp)

# target_compile_options(ShortCut_Lemon PUBLIC ${LEMON_COMPILE_OPTIONS})
# target_include_directories(ShortCut_Lemon PUBLIC ${LEMON_INCLUDE_DIRECTORY})
# target_link_libraries(ShortCut_Lemon
#         ${LEMON_LIBRARY}
#         ShortCutSolver
#         LP_Lemon
#         Common
#         )


add_executable(FastEMD   FastEMD.cpp)
target_link_libraries(FastEMD    Common)


# target_compile_options(ShortCut_Lemon PUBLIC ${LEMON_COMPILE_OPTIONS})
# target_include_directories(ShortCut_Lemon PUBLIC ${LEMON_INCLUDE_DIRECTORY})
# target_link_libraries(ShortCut_Lemon
#         ${LEMON_LIBRARY}
#         ShortCutSolver
#         LP_Lemon
#         Common
#         )



add_executable(SparseSK  SparseSK.cpp)
target_link_libraries(SparseSK  -lpthread  Common    SparseSinkhorn)

# log domain sinkhorn
# cuda implementation
add_executable(LogDomainSK  LogDomainSK.cpp)
target_link_libraries(LogDomainSK  LogSinkhorn)


# cplex


# set(CPLEX_LIBRARY /home/usr/cplex/cplexfile/cplex/lib/x86-64_linux/static_pic/libcplex.a
#     pthread )
# add_executable(TransportNetwork  CplexTransportNetwork.cpp)
# target_link_libraries(TransportNetwork
#     Common  ${CPLEX_LIBRARY}
#     -ldl LP_CPLEX   )



