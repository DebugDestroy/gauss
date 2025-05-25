# Базовые компоненты ядра
add_library(core
    config.cpp
    logger.cpp
)

target_include_directories(core PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)
