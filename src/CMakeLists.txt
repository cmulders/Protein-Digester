cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CPP_SOURCES
	CPPDigester.cpp
	fastaReader.cpp
	fastaProtein.cpp
)

find_package(Boost 1.54.0 COMPONENTS program_options filesystem system)
if(NOT Boost_FOUND)
	message( FATAL_ERROR "Boost not found" )
endif()
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIR})

add_executable(CPPDigester ${CPP_SOURCES})

TARGET_LINK_LIBRARIES( CPPDigester ${Boost_LIBRARIES} )
 
install(TARGETS CPPDigester DESTINATION ${PROJECT_SOURCE_DIR}/bin)
