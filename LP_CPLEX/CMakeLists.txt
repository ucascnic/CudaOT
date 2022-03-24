


set(CPLEX_COMPILE_OPTIONS -DNDEBUG -DIL_STD -DILOSTRICTPOD)
set(CPLEX_COMPILE_OPTIONS ${CPLEX_COMPILE_OPTIONS} -DUSE_CPLEX)
set(CPLEX_INCLUDE_DIRECTORY "/home/chenyidong/cplex/cplexfile/cplex/include")

add_library(LP_CPLEX STATIC TCPLEXNetSolver.cpp Interfaces-CPLEX.cpp)

target_compile_options(LP_CPLEX PUBLIC ${CPLEX_COMPILE_OPTIONS})
target_include_directories(LP_CPLEX PUBLIC ${CPLEX_INCLUDE_DIRECTORY})

install (TARGETS LP_CPLEX DESTINATION ${CMAKE_BINARY_DIR}/../bin)


install (FILES
        Interfaces-CPLEX.h
        TCPLEXNetSolver.h
        DESTINATION ${CMAKE_BINARY_DIR}/../include/LP_CPLEX)



