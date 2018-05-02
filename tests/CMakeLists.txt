set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories( ${CMAKE_SOURCE_DIR})
include_directories( ${CMAKE_SOURCE_DIR}/kratos)
include_directories( ${CMAKE_SOURCE_DIR}/applications)
include_directories( ${CMAKE_SOURCE_DIR}/applications/isogeometric_application)

###############################################################
set(name_list
    test_knot_array_1d
    test_knot_array_1d_test_2
    test_tsplines_1
    test_bezier_extraction_1d
    test_bezier_extraction_1d_test_2
    test_bezier_extraction_2d
    test_bezier_extraction_2d_test_2
    test_bezier_extraction_3d
    test_bezier_extraction_local_1d
    test_findspan_local_knots
    test_CreateRectangularControlPointGrid
)

foreach(str ${name_list})
    add_executable(${str} ${str}.cpp)
    target_link_libraries(${str} KratosCore)
    target_link_libraries(${str} KratosIsogeometricApplication)
    install(TARGETS ${str} DESTINATION libs )
endforeach()
