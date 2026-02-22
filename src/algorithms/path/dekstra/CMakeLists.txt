add_library(dekstra
    path_finder.cpp
)

target_include_directories(dekstra PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(dekstra
    PUBLIC
        path_common
)
