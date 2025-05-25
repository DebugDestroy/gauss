# Основная папка path
# Собираем библиотеку path из подмодулей
add_library(path INTERFACE)

# Добавляем подмодули
add_subdirectory(a_star)

# Прокидываем зависимости
target_link_libraries(path
    INTERFACE
        a_star
)
