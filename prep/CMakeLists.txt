cmake_minimum_required(VERSION 3.0)
project(PREPINPUT)

SET(CMAKE_BUILD_TYPE Release)

add_subdirectory(lib)

include_directories (
	${PREPINPUT_SOURCE_DIR}/lib/bamtools/src
	${PREPINPUT_SOURCE_DIR}/lib/config
	${PREPINPUT_SOURCE_DIR}/lib/fastahack
	${PREPINPUT_SOURCE_DIR}/lib/genomedata
	${PREPINPUT_SOURCE_DIR}/lib/inputparser
	${PREPINPUT_SOURCE_DIR}/lib/kent/inc
	${PREPINPUT_SOURCE_DIR}/lib/matrix
	${PREPINPUT_SOURCE_DIR}/lib/mydefine
	${PREPINPUT_SOURCE_DIR}/lib/split
)

add_executable(prepInput src/prepInput.cpp)

target_link_libraries(prepInput mydefine)

SET_TARGET_PROPERTIES(prepInput PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PREPINPUT_BINARY_DIR}/bin)
