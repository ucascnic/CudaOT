


cuda_add_library(SparseSinkhorn STATIC image_ot_solver.h  kernel_manage.h solver.h
    image_ot_solver.cu  kernel_manage.cu )

 

target_link_libraries(SparseSinkhorn  -lpthread
    Common  Sinkhorn -lcusparse -lcublas)



