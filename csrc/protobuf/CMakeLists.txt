file(GLOB ProtoFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.proto")
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders ${ProtoFiles})
add_library(fl2_protobuf STATIC ${ProtoSources} ${ProtoHeaders})
target_link_libraries(fl2_protobuf ${Boost_LIBRARIES} ${PROTOBUF_LIBRARY})