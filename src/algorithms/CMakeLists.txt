# Поддиректории с алгоритмами
add_subdirectory(components)
add_subdirectory(gauss)
add_subdirectory(geometry)
add_subdirectory(kinematics)
add_subdirectory(path)

# Собери библиотеку algorithms из подблоков
add_library(algorithms INTERFACE)

# Подключи зависимости (в виде библиотек)
target_link_libraries(algorithms
    INTERFACE
        components
        gauss
        geometry
        kinematics
        path
)
