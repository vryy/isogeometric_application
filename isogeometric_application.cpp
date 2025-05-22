//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2013 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "isogeometric_application.h"
#include "isogeometric_application_variables.h"
#include "custom_geometries/geo_1d_bezier.h"
#include "custom_geometries/geo_2d_bezier.h"
#include "custom_geometries/geo_2d_bezier_3.h"
#include "custom_geometries/geo_3d_bezier.h"

namespace Kratos
{

typedef Element::GeometryType::PointType NodeType;

KratosIsogeometricApplication::KratosIsogeometricApplication()
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("KratosIsogeometricApplication")
#else
    : KratosApplication()
#endif
      // , mDummyElementBezier( 0, Element::GeometryType::Pointer( new Geo1dBezier<NodeType>() ) )
      // , mDummyElementBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<NodeType>() ) )
      // , mDummyElementBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<NodeType>() ) )
      // , mDummyConditionBezier( 0, Element::GeometryType::Pointer( new Geo1dBezier<NodeType>() ) )
      // , mDummyConditionBezier2D( 0, Element::GeometryType::Pointer( new Geo2dBezier<NodeType>() ) )
      // , mDummyConditionBezier2D3( 0, Element::GeometryType::Pointer( new Geo2dBezier3<NodeType>() ) )
      // , mDummyConditionBezier3D( 0, Element::GeometryType::Pointer( new Geo3dBezier<NodeType>() ) )
{}

void KratosIsogeometricApplication::Register()
{
    std::cout << "Initializing KratosIsogeometricApplication... " << std::endl;

    //register variables
    KRATOS_REGISTER_VARIABLE( NURBS_WEIGHT )
    KRATOS_REGISTER_VARIABLE( NURBS_WEIGHTS )
    KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_1 )
    KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_2 )
    KRATOS_REGISTER_VARIABLE( NURBS_KNOTS_3 )
    KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_1 )
    KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_2 )
    KRATOS_REGISTER_VARIABLE( NURBS_DEGREE_3 )
    KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_1 )
    KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_2 )
    KRATOS_REGISTER_VARIABLE( NURBS_DIMENSION_3 )
    KRATOS_REGISTER_VARIABLE( NUM_DIVISION_1 )
    KRATOS_REGISTER_VARIABLE( NUM_DIVISION_2 )
    KRATOS_REGISTER_VARIABLE( NUM_DIVISION_3 )
    KRATOS_REGISTER_VARIABLE( NUM_IGA_INTEGRATION_METHOD )
    KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR )
    KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_MCSR )
    KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_ROWPTR )
    KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_COLIND )
    KRATOS_REGISTER_VARIABLE( EXTRACTION_OPERATOR_CSR_VALUES )
    KRATOS_REGISTER_VARIABLE_IMPLEMENTATION( CONTROL_POINT )
    KRATOS_REGISTER_VARIABLE_IMPLEMENTATION( COMPLEX_CONTROL_POINT )
    KRATOS_REGISTER_VARIABLE( KNOT_LEFT )
    KRATOS_REGISTER_VARIABLE( KNOT_RIGHT )
    KRATOS_REGISTER_VARIABLE( KNOT_TOP )
    KRATOS_REGISTER_VARIABLE( KNOT_BOTTOM )
    KRATOS_REGISTER_VARIABLE( KNOT_FRONT )
    KRATOS_REGISTER_VARIABLE( KNOT_BACK )
    KRATOS_REGISTER_VARIABLE( PATCH_INDEX )
    KRATOS_REGISTER_VARIABLE( HIERARCHICAL_LEVEL )
    KRATOS_REGISTER_VARIABLE( BASIS_FUNCTION_INDEX )
    KRATOS_REGISTER_VARIABLE( EQUATION_INDEX )
    KRATOS_REGISTER_VARIABLE( CELL_INDEX )

    // to make sure the variable imported from other application is registered
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_COORDINATES )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTROL_POINT_COORDINATES )

#ifdef SD_APP_FORWARD_COMPATIBILITY
    // DO NOTHING
#else
    // register the geometries
    Geo1dBezier<NodeType> Geo1dBezierPrototype;
    Serializer::Register( "Geo1dBezier", Geo1dBezierPrototype );

    Geo2dBezier<NodeType> Geo2dBezierPrototype;
    Serializer::Register( "Geo2dBezier", Geo2dBezierPrototype );

    Geo2dBezier3<NodeType> Geo2dBezier3Prototype;
    Serializer::Register( "Geo2dBezier3", Geo2dBezier3Prototype );

    Geo3dBezier<NodeType> Geo3dBezierPrototype;
    Serializer::Register( "Geo3dBezier", Geo3dBezierPrototype );
#endif

    // // register elements
    // KRATOS_REGISTER_ELEMENT( "DummyElementBezier", mDummyElementBezier )
    // KRATOS_REGISTER_ELEMENT( "DummyElementBezier2D", mDummyElementBezier2D )
    // KRATOS_REGISTER_ELEMENT( "DummyElementBezier3D", mDummyElementBezier3D )

    // // register conditions
    // KRATOS_REGISTER_CONDITION( "DummyConditionBezier", mDummyConditionBezier )
    // KRATOS_REGISTER_CONDITION( "DummyConditionBezier2D", mDummyConditionBezier2D )
    // KRATOS_REGISTER_CONDITION( "DummyConditionBezier2D3", mDummyConditionBezier2D3 )
    // KRATOS_REGISTER_CONDITION( "DummyConditionBezier3D", mDummyConditionBezier3D )
}

}  // namespace Kratos.
