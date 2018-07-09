set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIsogeometricApplication ****")


include_directories( ${CMAKE_SOURCE_DIR})
include_directories( ${CMAKE_SOURCE_DIR}/kratos)
include_directories( ${CMAKE_SOURCE_DIR}/applications)

# add_subdirectory(custom_external_libraries/geomod)
add_subdirectory(custom_external_libraries/tetgen1.5.0)
add_definitions( -DISOGEOMETRIC_USE_TETGEN )
# add_definitions( -DENABLE_BEZIER_GEOMETRY ) # this was promoted to system level
# add_definitions( -DENABLE_PRECOMPUTE )

if(DEFINED $ENV{HDF5_ROOT})
    SET(HDF5_DIR $ENV{HDF5_ROOT}/share/cmake/hdf5)
    FIND_PACKAGE(HDF5 REQUIRED NAMES hdf5 hdf5_cpp)
    set(HDF5_PROVIDED TRUE)
elseif(DEFINED HDF5_ROOT)
    SET(HDF5_DIR ${HDF5_ROOT}/share/cmake/hdf5)
    FIND_PACKAGE(HDF5 REQUIRED NAMES hdf5 hdf5_cpp)
    set(HDF5_PROVIDED TRUE)
else()
    set(HDF5_PROVIDED FALSE)
endif()

set(ISOGEOMETRIC_USE_HDF5 FALSE)
if(${HDF5_PROVIDED} MATCHES TRUE)
    message("HDF5_INCLUDE_DIR: " ${HDF5_INCLUDE_DIR})
    message("HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
#    add_definitions( -DBOOST_FILESYSTEM_VERSION=2 )
    add_definitions(-DISOGEOMETRIC_USE_HDF5)
    set(ISOGEOMETRIC_USE_HDF5 TRUE)
endif()

set(ISOGEOMETRIC_USE_MPI FALSE)
if(DEFINED ISOGEOMETRIC_PARMETIS_ROOT)
    INCLUDE_DIRECTORIES(${ISOGEOMETRIC_PARMETIS_ROOT}/include)
    FIND_LIBRARY(METIS_FOR_IGA_LIBRARY metis "${ISOGEOMETRIC_PARMETIS_ROOT}/lib" NO_DEFAULT_PATH)
    FIND_LIBRARY(PARMETIS_FOR_IGA_LIBRARY parmetis "${ISOGEOMETRIC_PARMETIS_ROOT}/lib" NO_DEFAULT_PATH)
    set(ISOGEOMETRIC_USE_MPI TRUE)
    add_definitions(-DISOGEOMETRIC_USE_PARMETIS)
    add_definitions(-DISOGEOMETRIC_USE_MPI)
endif()

set(ISOGEOMETRIC_USE_GISMO FALSE)
if(DEFINED GISMO_ROOT)
    INCLUDE_DIRECTORIES(${GISMO_ROOT}/include)
    INCLUDE_DIRECTORIES(${GISMO_ROOT}/include/gismo)
    FIND_LIBRARY(GISMO_LIBRARY gismo PATHS ${GISMO_ROOT}/lib NO_DEFAULT_PATH)
    set(ISOGEOMETRIC_USE_GISMO TRUE)
    add_definitions(-DISOGEOMETRIC_USE_GISMO)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/kratos/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_deactivation.h)
    add_definitions(-D_residualbased_elimination_builder_and_solver_deactivation_existed_)
endif()

## generate variables with the sources
set( KRATOS_ISOGEOMETRIC_APPLICATION_SOURCES
#    ${PROJECT_SOURCE_DIR}/kratos/sources/io.cpp
#    ${PROJECT_SOURCE_DIR}/kratos/sources/model_part_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/bezier_model_part_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/bezier_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/bspline_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/bezier_post_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs/bsplines_fespace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs/cell_manager.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs/domain_manager_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/nurbs/domain_manager_3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/hbsplines/deprecated_hb_basis_function.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/hbsplines/deprecated_hb_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/tsplines/tsmesh_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/import_export/multi_nurbs_patch_geo_importer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dummy_isogeometric_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/isogeometric_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_backend_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_frontend_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_transformation_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_control_grids_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_nurbs_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_hbsplines_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_tsplines_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mesh_and_model_part_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_patches_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/isogeometric_python_application.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosIsogeometricApplication SHARED ${KRATOS_ISOGEOMETRIC_APPLICATION_SOURCES})

add_subdirectory(tests)

target_link_libraries(KratosIsogeometricApplication KratosCore)
#target_link_libraries(KratosIsogeometricApplication KratosStructuralApplication)
target_link_libraries(KratosIsogeometricApplication iga_tetgen150)
if(${ISOGEOMETRIC_USE_HDF5} MATCHES TRUE)
    target_link_libraries(KratosIsogeometricApplication ${HDF5_LIBRARIES})
endif()
if(${ISOGEOMETRIC_USE_MPI} MATCHES TRUE)
    target_link_libraries(KratosIsogeometricApplication ${METIS_FOR_IGA_LIBRARY})
    target_link_libraries(KratosIsogeometricApplication ${PARMETIS_FOR_IGA_LIBRARY})
    target_link_libraries(KratosIsogeometricApplication ${MPI_LIBRARIES})
endif()
if(${ISOGEOMETRIC_USE_GISMO} MATCHES TRUE)
    target_link_libraries(KratosIsogeometricApplication ${GISMO_LIBRARY})
endif()

set_target_properties(KratosIsogeometricApplication PROPERTIES PREFIX "")
install(TARGETS KratosIsogeometricApplication DESTINATION libs )


# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosIsogeometricApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/IsogeometricApplication.py" DESTINATION KratosMultiphysics)

