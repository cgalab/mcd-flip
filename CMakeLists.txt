cmake_minimum_required (VERSION 3.4)
project (mcde)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

set (CMAKE_CXX_STANDARD_REQUIRED ON)
IF ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER "3.7")
  set (CMAKE_CXX_STANDARD 17)
ELSE ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
ENDIF ()

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall" )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -W -Wfloat-equal -Wundef -Wpointer-arith" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wwrite-strings" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wredundant-decls -Wchar-subscripts -Wcomment -Wformat=2" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wwrite-strings -Wmissing-declarations -Wredundant-decls" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wswitch-enum" )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Winit-self -Wmissing-field-initializers" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Waddress -Wmissing-noreturn -Wstrict-overflow=1" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wextra -Warray-bounds" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wshadow" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pedantic" )

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=switch-enum" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=shadow" )

### XXX REMOVE EVENTUALLY
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-parameter" )
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-unused-variable" )

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release."
       FORCE )
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" UPPERCASE_CMAKE_BUILD_TYPE)

if (NOT UPPERCASE_CMAKE_BUILD_TYPE MATCHES "^(DEBUG|RELEASE)$")
  message(FATAL_ERROR "Invalid value for CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

if(UPPERCASE_CMAKE_BUILD_TYPE MATCHES DEBUG)
  SET (DEBUG_OUTPUT_DEFAULT ON)
else()
  SET (DEBUG_OUTPUT_DEFAULT OFF)
endif()
if (NOT DEFINED DEBUG_OUTPUT)
  SET (DEBUG_OUTPUT ${DEBUG_OUTPUT_DEFAULT})
endif()

option(DEBUG_OUTPUT               "Include logging at debug level" ${DEBUG_OUTPUT})
option(DEBUG_OUTPUT_WITH_FILES    "Include filenames and line numbersin debug output" OFF)

# so we can run in valgrind
# add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON)

add_definitions(-DELPP_NO_DEFAULT_LOG_FILE)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  @ONLY
)

include_directories("${PROJECT_BINARY_DIR}")

subdirs (src)

enable_testing()
