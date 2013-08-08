#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Sect"

#include <boost/test/unit_test.hpp>

#include "sect/sect.hpp"

BOOST_AUTO_TEST_SUITE(SectTestSuite);
    BOOST_AUTO_TEST_CASE(FastaLoadTest) {
        BOOST_CHECK(true);
    }
BOOST_AUTO_TEST_SUITE_END();
