cmake_minimum_required (VERSION 3.1)

project(flux-analyze)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(PROG_NAME "flux-analyze")

set(CMAKE_CXX_FLAGS "-ggdb -Wall -O3")

add_executable(${PROG_NAME}
	address_marks.cc
	mfm_decoder.cc
	rabin_karp.cc
	flux_record.cc
	pulse_train.cc
	ordinal_search.cc
	tools.cc)

# TODO later: actually check that these libraries are installed

target_link_libraries(${PROG_NAME} -lz -lsqlite3)

