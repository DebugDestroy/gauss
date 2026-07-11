# Алгоритм RRT
add_library(rrt
    rrt.cpp
)

target_include_directories(rrt PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(rrt
    PUBLIC
        core
        gauss
        geometry
        path_common
)
