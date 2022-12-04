cmake_minimum_required(VERSION 3.0)

# Build the kent library
FILE(GLOB kent_src kent/lib/*.c)
include_directories(kent/inc)
add_library(kent ${kent_src})
target_link_libraries(kent m z)

# Build the Bamtools library
FILE(GLOB bamtools_api bamtools/src/api/*.cpp)
FILE(GLOB bamtools_internal bamtools/src/api/internal/*.cpp)

# create main BamTools API static library
include_directories(bamtools/src)
add_library(bamtools ${bamtools_api} ${bamtools_internal})
target_link_libraries(bamtools z)

# Build the matrix library
include_directories(matrix)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)

# Build the fastahack library
include_directories(fastahack)
add_library(fastahack fastahack/Fasta.cpp)
target_link_libraries(fastahack split)

# Build the inputparser library
include_directories(inputparser)
add_library(inputparser inputparser/InputParser.cpp)
target_link_libraries(inputparser mydefine)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config inputparser genomedata)

# Build the genomedata library
include_directories(genomedata)
add_library(genomedata genomedata/GenomeData.cpp)
target_link_libraries(genomedata config mydefine fastahack split bamtools kent)
