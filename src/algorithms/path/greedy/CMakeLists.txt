# Алгоритм Greedy
add_library(greedy
    path_finder.cpp
)

target_include_directories(greedy PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(greedy
    PUBLIC
        core
        gauss
        geometry
        path_common
)
