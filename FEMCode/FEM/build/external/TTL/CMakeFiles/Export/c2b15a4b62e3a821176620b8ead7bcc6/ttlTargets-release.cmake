#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ttl" for configuration "Release"
set_property(TARGET ttl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "/usr/local/lib/libttl.a"
  )

list(APPEND _cmake_import_check_targets ttl )
list(APPEND _cmake_import_check_files_for_ttl "/usr/local/lib/libttl.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
