cmake_minimum_required (VERSION 3.1)

project(flux-analyze)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(PROG_NAME "flux-analyze")

set(CMAKE_CXX_FLAGS "-ggdb -Wall -O3")

find_package(ZLIB REQUIRED)
find_package(SQLite3 REQUIRED)

add_executable(${PROG_NAME}
	address_marks.cc
	main.cc
	rabin_karp/rabin_karp.cc
	flux_record.cc
	pulse_train.cc
	ordinal_search.cc
	sector_data.cc
	preambles.cc
	timeline.cc
	tools.cc)

add_executable(copy_flux
	copy_flux.cc)

# Add SQLite and zlib libraries; if they're not found, the REQUIRED
# should've erred out above, so we won't get any errors here.

include_directories( ${ZLIB_INCLUDE_DIRS} ${SQLite3_INCLUDE_DIRS} )
target_link_libraries( ${PROG_NAME} ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
target_link_libraries( copy_flux ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
