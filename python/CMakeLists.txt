
# Find SWIG
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "-modern")
# Add include directories necessary for compilation of the new module
include_directories("${BEMPP_PYTHON_INCLUDE_DIRS}") # Python and NumPy headers
include_directories("${BEMPP_INCLUDE_DIR}")
include_directories("${BEMPP_INCLUDE_DIR}/bempp")
include_directories("${BEMPP_INCLUDE_DIR}/bempp/swig")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# Specify the source files and name of the new module
set_source_files_properties(core_periodic_laplace.i PROPERTIES CPLUSPLUS ON)
swig_add_module(core_periodic_laplace python core_periodic_laplace.i)

swig_link_libraries(core_periodic_laplace
    ${BEMPP_PYTHON_LIBRARY}
    periodic_laplace
    ${BEMPP_LIBRARY} ${BEMPP_TEUCHOSCORE_LIBRARY} ${BEMPP_TEUCHOSPARAMETERLIST_LIBRARY})

# Specify files to be installed

set(FILES_TO_INSTALL
    "${CMAKE_CURRENT_SOURCE_DIR}/periodic_laplace.py"
    "${CMAKE_CURRENT_BINARY_DIR}/core_periodic_laplace.py"
    "${CMAKE_CURRENT_BINARY_DIR}/_core_periodic_laplace.so")
install(FILES ${FILES_TO_INSTALL}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/python/bempp")
