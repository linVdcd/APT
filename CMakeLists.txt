project(APtoolbox)

include_directories(include)
include_directories(algrithom)
include_directories(3partLib/include)
add_subdirectory(3partLib)
add_subdirectory(algrithom)
add_subdirectory(test)
