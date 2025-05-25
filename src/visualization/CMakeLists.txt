# Визуализация через Gnuplot и цветовая палитра
add_library(visualization
    color.cpp
    gnuplot.cpp
)

target_include_directories(visualization PUBLIC ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(visualization
    core
    command
    algorithms
)
