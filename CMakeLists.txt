# ==============================================
# Terrain Navigator CMake Configuration
# ==============================================
# Доступные опции:
#   -DCMAKE_BUILD_TYPE=[Debug|Release] (по умолчанию Release)
#   -DENABLE_TESTING=ON/OFF (по умолчанию OFF)
#   -DINSTALL_SYSTEM=ON/OFF (установка в систему, по умолчанию OFF)
# 
# Примеры использования:
#   mkdir build && cd build
#   cmake .. -DCMAKE_BUILD_TYPE=Release
#   cmake --build . --parallel 4
# ==============================================

cmake_minimum_required(VERSION 3.12)
project(TerrainNavigator VERSION 1.0 LANGUAGES CXX)

# Опции конфигурации
option(ENABLE_TESTING "Enable testing" OFF)
option(INSTALL_SYSTEM "Install system-wide" OFF)

# Настройка выходных директорий
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Настройка стандарта C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Настройки компилятора
if(MSVC)
    add_compile_options(/W4 /WX /MP)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror -fPIC)
endif()

# Структура проекта
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src)

# Поиск исходников
file(GLOB_RECURSE SOURCES
    "${SOURCE_DIR}/app/*.cpp"
    "${SOURCE_DIR}/algorithms/*.cpp"
    "${SOURCE_DIR}/command/*.cpp"
    "${SOURCE_DIR}/core/*.cpp"
    "${SOURCE_DIR}/io/*.cpp"
    "${SOURCE_DIR}/utils/*.cpp"
    "${SOURCE_DIR}/visualization/*.cpp"
)

# Создание исполняемого файла
add_executable(terrain_navigator ${SOURCES})

# Подключение заголовков
target_include_directories(terrain_navigator PRIVATE
    ${INCLUDE_DIR}
    ${INCLUDE_DIR}/algorithms
    ${INCLUDE_DIR}/command
    ${INCLUDE_DIR}/core
    ${INCLUDE_DIR}/io
    ${INCLUDE_DIR}/utils
    ${INCLUDE_DIR}/visualization
)

# Проверка gnuplot
find_program(GNUPLOT_EXECUTABLE gnuplot)
if(NOT GNUPLOT_EXECUTABLE)
    message(WARNING "Gnuplot not found - visualization features will be limited")
endif()

# Установка (опционально)
if(INSTALL_SYSTEM)
    install(TARGETS terrain_navigator DESTINATION bin)
    install(DIRECTORY bin/etc DESTINATION share/terrain_navigator)
endif()

# Тестирование
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Режим сборки по умолчанию
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_BUILD_TYPE Release)
