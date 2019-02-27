/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 2013-08-18 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_control_point_to_python.h"
#include "custom_python/add_control_grids_to_python.h"
#include "custom_python/add_transformation_to_python.h"
#include "custom_python/add_nurbs_to_python.h"
#include "custom_python/add_pbbsplines_to_python.h"
#include "custom_python/add_hbsplines_to_python.h"
#include "custom_python/add_tsplines_to_python.h"
#include "custom_python/add_fespace_to_python.hpp"
#include "custom_python/add_grid_functions_to_python.hpp"
#include "custom_python/add_patches_to_python.h"
#include "custom_python/add_mesh_and_model_part_to_python.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_python/add_io_to_python.h"
#include "custom_python/add_strategies_to_python.h"

////utilities

// Project includes
#include "includes/define.h"
#include "isogeometric_application.h"


namespace Kratos
{

namespace Python
{


BOOST_PYTHON_MODULE(KratosIsogeometricApplication)
{

    using namespace boost::python;

    class_<KratosIsogeometricApplication,
           KratosIsogeometricApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable > ("KratosIsogeometricApplication")
           ;

    IsogeometricApplication_AddBackendUtilitiesToPython();
    IsogeometricApplication_AddFrontendUtilitiesToPython();
    IsogeometricApplication_AddControlPointToPython();
    IsogeometricApplication_AddControlGridsToPython();
    IsogeometricApplication_AddTransformationToPython();
    IsogeometricApplication_AddFESpacesToPython<1>();
    IsogeometricApplication_AddFESpacesToPython<2>();
    IsogeometricApplication_AddFESpacesToPython<3>();
    IsogeometricApplication_AddGridFunctionsToPython<1>();
    IsogeometricApplication_AddGridFunctionsToPython<2>();
    IsogeometricApplication_AddGridFunctionsToPython<3>();
    IsogeometricApplication_AddPatchesToPython();
    IsogeometricApplication_AddNURBSToPython();
    IsogeometricApplication_AddPBBSplinesToPython();
    IsogeometricApplication_AddHBSplinesToPython();
    IsogeometricApplication_AddTSplinesToPython();
    IsogeometricApplication_AddMeshAndModelPartToPython();
    IsogeometricApplication_AddProcessesToPython();
    IsogeometricApplication_AddIOToPython();
    IsogeometricApplication_AddStrategiesToPython();

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NURBS_WEIGHT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NURBS_WEIGHTS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUM_DIVISION_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUM_DIVISION_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUM_DIVISION_3 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONTROL_POINT_COORDINATES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUM_IGA_INTEGRATION_METHOD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTROL_POINT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_LEFT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_RIGHT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_TOP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_BOTTOM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_FRONT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KNOT_BACK )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PATCH_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( HIERARCHICAL_LEVEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BASIS_FUNCTION_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EQUATION_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CELL_INDEX )

}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
