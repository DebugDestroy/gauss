# Обработка команд
add_library(command
    control.cpp
    interface.cpp
)

target_include_directories(command PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(command
    core
    io
    visualization
    algorithms
)
