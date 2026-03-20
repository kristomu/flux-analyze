cmake_minimum_required (VERSION 3.1)

project(flux-analyze)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(PROG_NAME "flux-analyze")

set(CMAKE_CXX_FLAGS "-ggdb -Wall -O3 -Werror -Wextra -march=native -mtune=native")

find_package(ZLIB REQUIRED)
find_package(SQLite3 REQUIRED)
find_package(Protobuf REQUIRED)

add_subdirectory(protobuf)

add_executable(${PROG_NAME}
	pulse_train/alternating.cc
	pulse_train/classical_decoders.cc
	pulse_train/ewma.cc
	pulse_train/dp/one_mode_dp.cc
	pulse_train/offset_clock.cc
	rabin_karp/rabin_karp.cc
	address_marks.cc
	crc16.cc
	decoder.cc
	full_decoder.cc
	flux_record.cc
	main.cc
	ordinal_search.cc
	preambles.cc
	sector_data.cc
	timeline.cc
	tools.cc)

add_executable(grid-search
	pulse_train/alternating.cc
	pulse_train/classical_decoders.cc
	pulse_train/dp/one_mode_dp.cc
	pulse_train/ewma.cc
	pulse_train/kmedian_decoder.cc
	pulse_train/offset_clock.cc
	rabin_karp/rabin_karp.cc
	pulse_train/variable_offset.cc
	address_marks.cc
	crc16.cc
	decoder.cc
	full_decoder.cc
	flux_record.cc
	grid_search_main.cc
	ordinal_search.cc
	preambles.cc
	sector_data.cc
	stats/kmedian/kmedian.cc
	timeline.cc
	tools.cc)

add_executable(mfm_dp_test
	flux_record.cc
	dp/mfm_dp_test.cc)

add_executable(copy_flux
	copy_flux.cc)

# Add SQLite and zlib libraries, as well as the generated protobuf
# source; if they're not found, the REQUIRED should've erred out above,
# so we won't get any errors here.

# Protobuf solution from https://stackoverflow.com/a/20825093.

set(ProtobufIncludePath ${CMAKE_CURRENT_BINARY_DIR}
    CACHE INTERNAL "Path to generated protobuf files.")

include_directories( ${ZLIB_INCLUDE_DIRS} ${SQLite3_INCLUDE_DIRS} ${ProtobufIncludePath})
target_link_libraries( ${PROG_NAME} fl2_protobuf ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
target_link_libraries( grid_search fl2_protobuf ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
target_link_libraries( mfm_dp_test fl2_protobuf ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
target_link_libraries( copy_flux ${ZLIB_LIBRARIES} ${SQLite3_LIBRARIES} )
