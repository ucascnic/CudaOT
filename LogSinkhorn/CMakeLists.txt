
 
cuda_add_library(LogSinkhorn STATIC log_sinkhorn.h  log_sinkhorn.cu )



target_link_libraries(LogSinkhorn  -lpthread
    Common   -lcublas)





