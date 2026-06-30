add_library(statistics
    csv_writer.cpp
    statistics_manager.cpp
)

target_include_directories(statistics PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(statistics
    PUBLIC
        geometry
)
