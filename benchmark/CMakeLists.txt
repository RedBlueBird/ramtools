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

    # `cmake --build <dir> --target benchmark` runs the whole suite through the
    # runner, which gives each binary the flags it takes and writes the JSON to
    # <dir>/benchmark/results; RAMTOOLS_BENCH_ARGS holds the flags. Relative
    # paths in them are taken from the source directory.
    set(RAMTOOLS_BENCH_ARGS "" CACHE STRING
        "Flags for the benchmark target, e.g. --sam=reads.sam --threads=1,4")
    separate_arguments(ramtools_bench_args UNIX_COMMAND "${RAMTOOLS_BENCH_ARGS}")
    add_custom_target(benchmark
        COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/run_benchmark_suite.sh ${CMAKE_BINARY_DIR} -- ${ramtools_bench_args}
        DEPENDS conversion_time_benchmark region_query_benchmark chromosome_split_benchmark columnar_read_benchmark
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running all RAM tools benchmarks"
        VERBATIM
    )
endif()
