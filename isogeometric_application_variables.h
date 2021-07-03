//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 18, 2013 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_VARIABLES_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/control_point.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition
    KRATOS_DEFINE_VARIABLE( double, NURBS_WEIGHT )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_WEIGHTS )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_1 )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_2 )
    KRATOS_DEFINE_VARIABLE( Vector, NURBS_KNOTS_3 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_1 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_2 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DEGREE_3 )
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_1 ) //number of control points in 1st direction
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_2 ) //number of control points in 2nd direction
    KRATOS_DEFINE_VARIABLE( int, NURBS_DIMENSION_3 ) //number of control points in 3rd direction
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_1 ) //number of mesh points along 1st direction in post-processing
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_2 ) //number of mesh points along 2nd direction in post-processing
    KRATOS_DEFINE_VARIABLE( int, NUM_DIVISION_3 ) //number of mesh points along 3rd direction in post-processing
    KRATOS_DEFINE_VARIABLE( int, NUM_IGA_INTEGRATION_METHOD )
    KRATOS_DEFINE_VARIABLE( Matrix, EXTRACTION_OPERATOR )
    KRATOS_DEFINE_VARIABLE( Matrix, EXTRACTION_OPERATOR_MCSR )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_ROWPTR )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_COLIND )
    KRATOS_DEFINE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_VALUES )
    KRATOS_DEFINE_VARIABLE( ControlPoint<double>, CONTROL_POINT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_COORDINATES)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_COORDINATES)
    KRATOS_DEFINE_VARIABLE( double, KNOT_LEFT )
    KRATOS_DEFINE_VARIABLE( double, KNOT_RIGHT )
    KRATOS_DEFINE_VARIABLE( double, KNOT_TOP )
    KRATOS_DEFINE_VARIABLE( double, KNOT_BOTTOM )
    KRATOS_DEFINE_VARIABLE( double, KNOT_FRONT )
    KRATOS_DEFINE_VARIABLE( double, KNOT_BACK )
    KRATOS_DEFINE_VARIABLE( int, PATCH_INDEX )
    KRATOS_DEFINE_VARIABLE( int, HIERARCHICAL_LEVEL )
    KRATOS_DEFINE_VARIABLE( int, BASIS_FUNCTION_INDEX )
    KRATOS_DEFINE_VARIABLE( int, EQUATION_INDEX )
    KRATOS_DEFINE_VARIABLE( int, CELL_INDEX )

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_VARIABLES_H_INCLUDED  defined


