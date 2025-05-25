# Обработка BMP-файлов
add_library(io
    bmp_handler.cpp
)

target_include_directories(io PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(io
    core
    algorithms
)
