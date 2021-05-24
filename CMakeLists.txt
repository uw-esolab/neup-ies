cmake_minimum_required(VERSION 3.12)

Project(system_advisor_model VERSION 1.0.0)

add_subdirectory(lk)
add_subdirectory(wex)
add_subdirectory(ssc)
add_subdirectory(sam)

option(SAM_SKIP_TOOLS "Skips the sdktool and tcsconsole builds" ON)
option(SAM_SKIP_TESTS "Skips building tests" ON)
option(SAMAPI_EXPORT "Export of ssc binaries to the SAM_api directory; for Unix, compile ssc libraries for SAM_api" ON)
option(SAM_SKIP_AUTOGEN "Re-generate the SAMAPI files from export_config" OFF)
