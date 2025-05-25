# Алгоритм A*
add_library(a_star
    conditions.cpp
    graph.cpp
    heuristics.cpp
    path_finder.cpp
)

target_include_directories(a_star PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(a_star
    PUBLIC
        core
        gauss
        geometry
        kinematics
)
