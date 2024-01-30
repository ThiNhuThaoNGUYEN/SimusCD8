# ============================================================================
# Find packages
# ============================================================================
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

#  # It looks like there's a bug in CMake support for GLUT when libXi or libXmu are
#  # not present on the system. See
#  # http://public.kitware.com/Bug/view.php?id=14060 for more info and workaround
#  #
#  # Workaround (from the aforementioned webpage):
#  set(_glut_libraries)
#  foreach(_lib ${GLUT_LIBRARIES})
#    if(_lib)
#      list(APPEND _glut_libraries ${_lib})
#    endif()
#  endforeach()
#  set(GLUT_LIBRARIES ${_glut_libraries})



# ============================================================================
# Add target view add its dependencies
# ============================================================================
add_executable(view view.cpp)
#add_executable(viewGene1 viewGene1.cpp)
#add_executable(viewGene2 viewGene2.cpp)
#add_executable(viewM1 viewM1.cpp)
#add_executable(viewM2 viewM2.cpp)


# ============================================================================
# Use OpenGL and GLUT
# ============================================================================
target_link_libraries(view ${OPENGL_LIBRARIES} GLUT::GLUT )


# ============================================================================
# Set behaviour on make install
# ============================================================================
#install(TARGETS view DESTINATION ${CMAKE_INSTALL_BINDIR}/simuscale)
