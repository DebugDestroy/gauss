# Геометрические структуры и алгоритмы
add_library(geometry
    bresenham_line.cpp
    math.cpp
    triangulator.cpp
    voronoi_diagram.cpp
)

target_include_directories(geometry PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(geometry
    PUBLIC
        core        # constants.hpp, logger.hpp
        gauss       # pole.hpp
)
