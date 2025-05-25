# Укажи исходники главного приложения
set(SRC_APP app/main.cpp)

# Создай главный исполняемый файл
add_executable(terrain_navigator ${SRC_APP})

# Заголовочные директории
target_include_directories(terrain_navigator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Добавь подпапки (каждая будет библиотекой)
add_subdirectory(algorithms)
add_subdirectory(command)
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(visualization)

# Линкуем библиотеки к исполняемому файлу
target_link_libraries(terrain_navigator
    algorithms
    command
    core
    io
    visualization
)

# Копирование конфига после сборки
add_custom_command(
    TARGET terrain_navigator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/config
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config
    COMMENT "Copying config/ to runtime output directory"
)
