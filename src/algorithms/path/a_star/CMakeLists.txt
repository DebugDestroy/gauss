add_library(a_star
    path_finder.cpp
)

target_include_directories(a_star PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(a_star
    PUBLIC
        path_common
)
