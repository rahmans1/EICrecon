
# Check if catch2 exists and unit tests can be built
find_package(Catch2 3)
if(Catch2_FOUND)
    add_subdirectory(algorithms_test)
else()
    message(STATUS "Catch2 is not found. Skipping algorithms_test...")
endif()

add_subdirectory(pyjano)
add_subdirectory(BEMC_test)
add_subdirectory(reco_test)
add_subdirectory(tracking_test)
add_subdirectory(track_propagation_test)
add_subdirectory(track_seeding_test)
add_subdirectory(geometry_navigation_test)
add_subdirectory(podio_test)
