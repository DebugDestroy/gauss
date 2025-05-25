# Модуль кинематики
add_library(kinematics
    incline_angle.cpp
)

target_include_directories(kinematics PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(kinematics
    PUBLIC
        core        # на всякий случай, если понадобится logger или constants
        geometry    # geometry_structures.hpp
        gauss       # pole.hpp
)
