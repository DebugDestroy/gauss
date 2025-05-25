# Файлы компонентов (анализ связности, кластеризация и т.п.)
add_library(components
    binary.cpp
    cluster_service.cpp
    components_analysis.cpp
    copier_service.cpp
    find_components.cpp
    kmeans.cpp
)

target_include_directories(components PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(components
    PUBLIC
        core           # logger, constants, config
        gauss          # pole.hpp
        geometry       # geometry_structures.hpp
        visualization  # color.hpp
)
