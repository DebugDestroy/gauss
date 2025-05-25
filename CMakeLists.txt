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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_TESTING "Enable testing" OFF)
option(INSTALL_SYSTEM "Install system-wide" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Добавим подпапку src — пока только её!
add_subdirectory(src)

# Проверка gnuplot
find_program(GNUPLOT_EXECUTABLE gnuplot)
if(NOT GNUPLOT_EXECUTABLE)
    message(WARNING "Gnuplot not found - visualization features will be limited")
endif()

# Установка
if(INSTALL_SYSTEM)
    install(TARGETS terrain_navigator DESTINATION bin)
    install(DIRECTORY config DESTINATION share/terrain_navigator)
endif()

# Тесты
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Режим по умолчанию
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif()
