
# Find all source and header files
file(GLOB_RECURSE LIB_SOURCES *.cpp)
file(GLOB_RECURSE LIB_HEADERS *.hpp)
include_directories("${BEMPP_INCLUDE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_library(periodic_laplace SHARED
  ${LIB_SOURCES}
  ${LIB_HEADERS}
)

# Add includes and target links

target_link_libraries(periodic_laplace
    ${BEMPP_PYTHON_LIBRARY}
    ${BEMPP_LIBRARY} ${BEMPP_TEUCHOSCORE_LIBRARY} ${BEMPP_TEUCHOSPARAMETERLIST_LIBRARY})

# Install library
install(TARGETS periodic_laplace
  EXPORT BemppLibraryDepends
  RUNTIME DESTINATION bempp/bin
  LIBRARY DESTINATION bempp/lib
  ARCHIVE DESTINATION bempp/lib)

# Install headers
install(DIRECTORY ${CMAKE_SOURCE_DIR}/lib/ DESTINATION bempp/include/bempp FILES_MATCHING PATTERN "*.hpp")

SET_TARGET_PROPERTIES(periodic_laplace
	PROPERTIES INSTALL_RPATH "${BEMPP_RPATH}")
