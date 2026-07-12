# Собираем библиотеку path из подмодулей
add_library(path INTERFACE)

# Добавляем подмодули
add_subdirectory(common)
add_subdirectory(a_star)
add_subdirectory(dekstra)
add_subdirectory(greedy)
add_subdirectory(rrt)
add_subdirectory(rrt_star)

# Прокидываем зависимости
target_link_libraries(path
    INTERFACE
        path_common
        a_star
        dekstra
        greedy
        rrt
        rrt_star
)
