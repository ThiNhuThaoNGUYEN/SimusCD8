# ============================================================================
# Generate config.h
# ============================================================================
configure_file(
        ${PROJECT_SOURCE_DIR}/src/config.h.in
        ${PROJECT_BINARY_DIR}/src/config.h
)
include_directories("${PROJECT_BINARY_DIR}/src ")

# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
#add_subdirectory(view)
add_subdirectory(simuscale-core)

# ============================================================================
# Add include directories
# ============================================================================
include_directories("." "simuscale-core")

# ============================================================================
# List all available plugins
# ============================================================================
#set(PLUGINS Cell_Spiky Cell_SyncClock  Lymphocyte )
set( PLUGINS Lymphocyte )

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(PLUGIN_DIR plugins)
foreach (PLUGIN IN LISTS PLUGINS)
    set(PLUGIN_SOURCES ${PLUGIN_SOURCES} ${PLUGIN_DIR}/${PLUGIN}.h ${PLUGIN_DIR}/${PLUGIN}.cpp)
endforeach(PLUGIN)


# ============================================================================
# Add target simuscale add its dependencies
# ============================================================================
add_executable(simuscale main.cpp ${PLUGIN_SOURCES})

target_link_libraries(simuscale PUBLIC simuscale-core GSL::gsl GSL::gslcblas)
# We use C++11
target_compile_options(simuscale PRIVATE "-std=c++11")
# Make STDC MACROS available (for fixed width integers)
target_compile_definitions(simuscale PUBLIC __STDC_FORMAT_MACROS
                                     PUBLIC __STDC_CONSTANT_MACROS)


# ============================================================================
# Set behaviour on make install
# ============================================================================
install(TARGETS simuscale RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )


# ============================================================================
# Include tests
# ============================================================================
add_subdirectory(tests EXCLUDE_FROM_ALL)
