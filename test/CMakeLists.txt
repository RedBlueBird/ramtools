list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
include(RamtoolsGoogleTest)

function(add_ramcore_test test_name)
    ROOT_EXECUTABLE(${test_name}
        ${ARGN}
        LIBRARIES
            ramcore
            ROOT::Core
            ROOT::ROOTNTuple
            ROOT::RIO
            sam_generator
            gtest
            gtest_main
    )

    if(RAMTOOLS_BUILD_BENCHMARKS)
        target_link_libraries(${test_name} PUBLIC benchmark::benchmark)
    endif()

    target_include_directories(${test_name}
        PRIVATE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/inc>
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/benchmark>
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_ramcore_test(ramcoretests ramcoretests.cxx)
add_ramcore_test(chromosome_split_test chromosome_split_test.cxx)
add_ramcore_test(bam_conversion bam_conversion.cxx)
add_ramcore_test(parallel_conversion_test parallel_conversion_test.cxx)
add_ramcore_test(record_blocks_test record_blocks_test.cxx)

if(TARGET samtoramntuple)
    add_test(NAME samtoramntuple_rejects_bad_compression COMMAND samtoramntuple in.sam -compression asd)
    add_test(NAME bamtoramntuple_rejects_bad_compression COMMAND bamtoramntuple in.bam -compression 999)
    set_tests_properties(samtoramntuple_rejects_bad_compression bamtoramntuple_rejects_bad_compression
        PROPERTIES PASS_REGULAR_EXPRESSION "invalid -compression value")
    add_test(NAME samtoramntuple_rejects_bad_threads COMMAND samtoramntuple in.sam -threads 0)
    add_test(NAME bamtoramntuple_rejects_bad_threads COMMAND bamtoramntuple in.bam -threads x)
    # A number followed by garbage must not be silently truncated to the number.
    add_test(NAME samtoramntuple_rejects_trailing_threads COMMAND samtoramntuple in.sam -threads 3sdfd)
    add_test(NAME bamtoramntuple_rejects_trailing_threads COMMAND bamtoramntuple in.bam -threads 3sdfd)
    set_tests_properties(samtoramntuple_rejects_bad_threads bamtoramntuple_rejects_bad_threads
        samtoramntuple_rejects_trailing_threads bamtoramntuple_rejects_trailing_threads
        PROPERTIES PASS_REGULAR_EXPRESSION "invalid -threads value")
    add_test(NAME samtoramntuple_rejects_fqzcomp_with_illumina COMMAND samtoramntuple in.sam -fqzcomp -illumina)
    add_test(NAME bamtoramntuple_rejects_fqzcomp_with_dropqual COMMAND bamtoramntuple in.bam -dropqual -fqzcomp)
    set_tests_properties(samtoramntuple_rejects_fqzcomp_with_illumina bamtoramntuple_rejects_fqzcomp_with_dropqual
        PROPERTIES PASS_REGULAR_EXPRESSION "-fqzcomp compresses lossless qualities")
endif()

if(TARGET conversion_time_benchmark)
    # The benchmarks check their flags as the converters do, before anything is
    # measured: a typo must not give numbers for the default settings.
    add_test(NAME benchmark_rejects_bad_compression COMMAND conversion_time_benchmark --compression=999)
    add_test(NAME benchmark_rejects_trailing_reads COMMAND conversion_time_benchmark --reads=1000x)
    add_test(NAME benchmark_rejects_empty_threads COMMAND conversion_time_benchmark --threads=2,)
    add_test(NAME benchmark_rejects_bad_quality COMMAND conversion_time_benchmark --quality=lossy)
    set_tests_properties(benchmark_rejects_bad_compression benchmark_rejects_trailing_reads
        benchmark_rejects_empty_threads benchmark_rejects_bad_quality
        PROPERTIES PASS_REGULAR_EXPRESSION "invalid --[a-z]+ value")
    # A flag a binary does not take is left to Google Benchmark, which rejects it.
    add_test(NAME benchmark_rejects_flag_it_does_not_take COMMAND region_query_benchmark --compression=505)
    set_tests_properties(benchmark_rejects_flag_it_does_not_take
        PROPERTIES PASS_REGULAR_EXPRESSION "unrecognized command-line flag")
    add_test(NAME benchmark_sweeps_reads_and_threads
        COMMAND conversion_time_benchmark --reads=10,20 --threads=1,3 --benchmark_list_tests)
    set_tests_properties(benchmark_sweeps_reads_and_threads
        PROPERTIES PASS_REGULAR_EXPRESSION "reads:20/threads:3")
    add_test(NAME benchmark_takes_regions_from_environment
        COMMAND region_query_benchmark --reads=1000 --benchmark_list_tests)
    set_tests_properties(benchmark_takes_regions_from_environment PROPERTIES
        ENVIRONMENT "RAMTOOLS_BENCH_REGIONS=chr1:1-10,chr2:1-10,chr3:1-10"
        PASS_REGULAR_EXPRESSION "RegionQuery/r2/RAM"
        FAIL_REGULAR_EXPRESSION "RegionQuery/r3/")
    # Inputs: a path must exist, and inputs given two ways are refused.
    add_test(NAME benchmark_rejects_missing_sam COMMAND conversion_time_benchmark --sam=no_such_file.sam)
    set_tests_properties(benchmark_rejects_missing_sam PROPERTIES PASS_REGULAR_EXPRESSION "no such file")
    add_test(NAME benchmark_rejects_sam_with_reads
        COMMAND conversion_time_benchmark --sam=${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt --reads=10)
    set_tests_properties(benchmark_rejects_sam_with_reads PROPERTIES PASS_REGULAR_EXPRESSION "exclude each other")
    # CRAM is encoded against a reference, and is never measured without one.
    add_test(NAME benchmark_rejects_cram_without_ref COMMAND conversion_time_benchmark --formats=ram,cram)
    set_tests_properties(benchmark_rejects_cram_without_ref PROPERTIES PASS_REGULAR_EXPRESSION "needs --ref")
    add_test(NAME benchmark_rejects_repeated_format COMMAND conversion_time_benchmark --formats=bam,bam)
    set_tests_properties(benchmark_rejects_repeated_format PROPERTIES PASS_REGULAR_EXPRESSION "invalid --formats")
    # One run of each binary on a few synthetic records, in a temporary
    # directory: the default path works without any dataset and nothing fails.
    add_test(NAME benchmark_conversion_runs
        COMMAND conversion_time_benchmark --reads=2000 --threads=1,2 --benchmark_min_time=1x)
    add_test(NAME benchmark_region_query_runs
        COMMAND region_query_benchmark --reads=5000 --benchmark_min_time=1x)
    add_test(NAME benchmark_chromosome_split_runs
        COMMAND chromosome_split_benchmark --reads=2000 --threads=1,2 --benchmark_min_time=1x)
    add_test(NAME benchmark_columnar_read_runs
        COMMAND columnar_read_benchmark --reads=5000 --benchmark_min_time=1x)
    set_tests_properties(benchmark_conversion_runs benchmark_region_query_runs benchmark_chromosome_split_runs
        benchmark_columnar_read_runs PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR OCCURRED")
endif()

if(TARGET ramdump)
    add_ramcore_test(ramdump_test ramdump_test.cxx)
    add_dependencies(ramdump_test samtoramntuple ramdump ramstats)
    target_compile_definitions(ramdump_test PRIVATE
        SAMTORAMNTUPLE_BIN="$<TARGET_FILE:samtoramntuple>"
        RAMDUMP_BIN="$<TARGET_FILE:ramdump>"
        RAMSTATS_BIN="$<TARGET_FILE:ramstats>")
endif()

install(TARGETS ramcoretests chromosome_split_test bam_conversion parallel_conversion_test record_blocks_test
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
