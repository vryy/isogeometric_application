//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2013 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_utilities/multipatch.h"
#include "isogeometric_application_variables.h"

namespace Kratos
{

typedef ControlPoint<double, double> ControlPointType;
typedef ControlPoint<std::complex<double>, double> ComplexControlPointType;

KRATOS_CREATE_VARIABLE( double, NURBS_WEIGHT )
KRATOS_CREATE_VARIABLE( Vector, NURBS_WEIGHTS )
KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_1 )
KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_2 )
KRATOS_CREATE_VARIABLE( Vector, NURBS_KNOTS_3 )
KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_1 )
KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_2 )
KRATOS_CREATE_VARIABLE( int, NURBS_DEGREE_3 )
KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_1 )
KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_2 )
KRATOS_CREATE_VARIABLE( int, NURBS_DIMENSION_3 )
KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_1 )
KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_2 )
KRATOS_CREATE_VARIABLE( int, NUM_DIVISION_3 )
KRATOS_CREATE_VARIABLE( int, NUM_IGA_INTEGRATION_METHOD )
KRATOS_CREATE_VARIABLE( Matrix, EXTRACTION_OPERATOR )
KRATOS_CREATE_VARIABLE( Matrix, EXTRACTION_OPERATOR_MCSR )
KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_ROWPTR )
KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_COLIND )
KRATOS_CREATE_VARIABLE( Vector, EXTRACTION_OPERATOR_CSR_VALUES )
KRATOS_CREATE_VARIABLE( ControlPointType, CONTROL_POINT )
KRATOS_CREATE_VARIABLE( ComplexControlPointType, COMPLEX_CONTROL_POINT )

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONTROL_POINT_COORDINATES )

KRATOS_CREATE_VARIABLE( double, KNOT_LEFT )
KRATOS_CREATE_VARIABLE( double, KNOT_RIGHT )
KRATOS_CREATE_VARIABLE( double, KNOT_TOP )
KRATOS_CREATE_VARIABLE( double, KNOT_BOTTOM )
KRATOS_CREATE_VARIABLE( double, KNOT_FRONT )
KRATOS_CREATE_VARIABLE( double, KNOT_BACK )

KRATOS_CREATE_VARIABLE( int, PATCH_INDEX )
KRATOS_CREATE_VARIABLE( int, HIERARCHICAL_LEVEL )
KRATOS_CREATE_VARIABLE( int, BASIS_FUNCTION_INDEX )
KRATOS_CREATE_VARIABLE( int, EQUATION_INDEX )
KRATOS_CREATE_VARIABLE( int, CELL_INDEX )

template<> PatchSelector<1>::RealPatch PatchSelector<1>::msRealPatch(1);
template<> PatchSelector<2>::RealPatch PatchSelector<2>::msRealPatch(1);
template<> PatchSelector<3>::RealPatch PatchSelector<3>::msRealPatch(1);
template<> PatchSelector<1>::ComplexPatch PatchSelector<1>::msComplexPatch(1);
template<> PatchSelector<2>::ComplexPatch PatchSelector<2>::msComplexPatch(1);
template<> PatchSelector<3>::ComplexPatch PatchSelector<3>::msComplexPatch(1);

PatchSelector<1> PatchSelector1DInstance;
PatchSelector<2> PatchSelector2DInstance;
PatchSelector<3> PatchSelector3DInstance;

}  // namespace Kratos.
