# Алгоритм RRT*
add_library(rrt_star
    rrt_star.cpp
)

target_include_directories(rrt_star PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(rrt_star
    PUBLIC
        core
        gauss
        geometry
        path_common
)
