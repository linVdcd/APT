aux_source_directory(. DIR_LIB_SRCS)
include_directories(../algrithom)
include_directories(../include)
include_directories(../3partLib/include)

add_library(APtool ${DIR_LIB_SRCS})
target_link_libraries(APtool 3part)
