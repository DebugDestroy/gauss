# Алгоритмы smoothing
add_library(smoothing
    spline.cpp
    shortcut.cpp
)

target_include_directories(smoothing PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(smoothing
    PUBLIC
        core
        gauss
        geometry
        path_common
)
