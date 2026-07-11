add_library(path_common
    collision.cpp
    graph.cpp
    grid.cpp
    path_validator.cpp
)

target_include_directories(path_common PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(path_common
    PUBLIC
        core
        gauss
        geometry
        kinematics
)
