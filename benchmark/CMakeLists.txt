list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

if(RAMTOOLS_BUILD_TESTS OR RAMTOOLS_BUILD_BENCHMARKS)
    add_library(sam_generator STATIC
        generate_sam.cxx
    )

    target_include_directories(sam_generator PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

if(RAMTOOLS_BUILD_BENCHMARKS)
    include(GoogleBenchmark)

    add_library(benchmark_common STATIC
        bam_ops.cxx
        benchmark_config.cxx
        benchmark_utils.cxx
    )

    target_include_directories(benchmark_common PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )

    target_link_libraries(benchmark_common PUBLIC ramcore sam_generator PkgConfig::htslib)

    ROOT_EXECUTABLE(conversion_time_benchmark
        conversion_time_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
            benchmark_common
    )

    ROOT_EXECUTABLE(region_query_benchmark
        region_query_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
            benchmark_common
    )

    ROOT_EXECUTABLE(chromosome_split_benchmark
        chromosome_split_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            sam_generator
            benchmark_common
    )

    ROOT_EXECUTABLE(columnar_read_benchmark
        columnar_read_benchmark.cxx
        LIBRARIES
            benchmark::benchmark
            ramcore
            benchmark_common
    )

    install(TARGETS conversion_time_benchmark region_query_benchmark chromosome_split_benchmark columnar_read_benchmark
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    add_custom_target(benchmark
        COMMAND ${CMAKE_COMMAND} -E echo "=== Conversion Time Benchmark ==="
        COMMAND conversion_time_benchmark
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Region Query Benchmark ==="
        COMMAND region_query_benchmark
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Chromosome Split Benchmark ==="
        COMMAND chromosome_split_benchmark
        COMMAND ${CMAKE_COMMAND} -E echo ""
        COMMAND ${CMAKE_COMMAND} -E echo "=== Columnar Read Benchmark ==="
        COMMAND columnar_read_benchmark
        DEPENDS conversion_time_benchmark region_query_benchmark chromosome_split_benchmark columnar_read_benchmark
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running all RAM tools benchmarks"
    )
endif()
