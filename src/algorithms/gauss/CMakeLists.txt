# Файлы для метода Гаусса и построения поля
add_library(gauss
    pole.cpp
    gauss_builder.cpp
)

target_include_directories(gauss PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(gauss
    PUBLIC
        core # logger.hpp, constants.hpp
)
