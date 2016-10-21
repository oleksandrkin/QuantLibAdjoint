file(GLOB_RECURSE BOOST_SYSTEM "system/src/*.cpp")
file(GLOB_RECURSE BOOST_ARCHIVE "archive/impl/src/*.cpp")
file(GLOB_RECURSE BOOST_FILESYSTEM "filesystem/src/*.cpp")
file(GLOB_RECURSE BOOST_TEST
                                "test/impl/compiler_log_formatter.ipp"
                                "test/impl/debug.ipp"
                                "test/impl/exception_safety.ipp"
                                "test/impl/execution_monitor.ipp"
                                "test/impl/framework.ipp"
                                "test/impl/interaction_based.ipp"
                                "test/impl/logged_expectations.ipp"
                                "test/impl/plain_report_formatter.ipp"
                                "test/impl/progress_monitor.ipp"
                                "test/impl/results_collector.ipp"
                                "test/impl/results_reporter.ipp"
                                "test/impl/test_tools.ipp"
                                "test/impl/unit_test_log.ipp"
                                "test/impl/unit_test_main.ipp"
                                "test/impl/unit_test_monitor.ipp"
                                "test/impl/unit_test_parameters.ipp"
                                "test/impl/unit_test_suite.ipp"
                                "test/impl/xml_log_formatter.ipp"
                                "test/impl/xml_report_formatter.ipp"
    )
set_source_files_properties(${BOOST_TEST} PROPERTIES LANGUAGE CXX)
set_source_files_properties(${BOOST_TEST} PROPERTIES COMPILE_FLAGS "-x c++")
add_library(BoostLib SHARED ${BOOST_SYSTEM} ${BOOST_ARCHIVE} ${BOOST_FILESYSTEM} ${BOOST_TEST})
