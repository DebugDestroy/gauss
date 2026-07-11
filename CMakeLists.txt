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

# Режим по умолчанию
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif()

if (MSVC)
    add_compile_options(
        /W4
        /permissive-
        /WX
    )
else()
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Werror

        -Wshadow
        -Wconversion
        -Wsign-conversion
        -Wfloat-equal
        -Wundef
        -Wcast-qual
        -Wold-style-cast
        -Woverloaded-virtual
        -Wnon-virtual-dtor
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Wunused
        -Wunused-parameter
        -Wmissing-declarations
    )
   
    add_link_options(
        -Wl,--no-undefined
        -Wl,--as-needed
        -Wl,--no-copy-dt-needed-entries
        -Wl,-z,relro
        -Wl,-z,now
    )
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_TESTING "Enable testing" OFF)
option(INSTALL_SYSTEM "Install system-wide" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        -fsanitize=address
        -fsanitize=undefined
    )

    add_link_options(
        -fsanitize=address
        -fsanitize=undefined
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-fdiagnostics-color=always)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wmaybe-uninitialized)
endif()

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

# Добавим подпапку src — пока только её!
add_subdirectory(src)
