/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Nov 25 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_2D_BEZIER_3_H_INCLUDED )
#define  KRATOS_GEO_2D_BEZIER_3_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "integration/quadrature.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_geometries/geo_2d_bezier.h"
//#include "integration/quadrature.h"
//#include "integration/line_gauss_legendre_integration_points.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
// #define DEBUG_LEVEL3
#define ENABLE_PROFILING

namespace Kratos
{

/**
 * A geometry representing Bezier decomposition surface in 3D
 */
template<class TPointType>
class Geo2dBezier3 : public Geo2dBezier<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef Geo2dBezier<TPointType> BaseType;

    /**
     * Pointer definition of Geo2dBezier3
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo2dBezier3 );

    /**
     * Integration methods implemented in geometry.
     */
    typedef typename BaseType::IntegrationMethod IntegrationMethod;

    /**
     * A VectorType of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A VectorType of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A VectorType of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType
    ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType
    ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType
    ShapeFunctionsSecondDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef typename BaseType::MatrixType MatrixType;

    /**
     * Type of Vector
     */
    typedef typename BaseType::VectorType VectorType;

    /**
     * Type of values container
     */
    typedef typename BaseType::NormalType ValuesContainerType;

    /**
     * Life Cycle
     */

    Geo2dBezier3()
//    : BaseType( PointsArrayType(), &msGeometryData )
    : BaseType( PointsArrayType() )
    {}

    Geo2dBezier3( const PointsArrayType& ThisPoints )
//    : BaseType( ThisPoints, &msGeometryData )
    : BaseType( ThisPoints )
    {
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Geo2dBezier3( Geo2dBezier3 const& rOther )
    : BaseType( rOther )
    {}

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Geo2dBezier3( Geo2dBezier3<TOtherPointType> const& rOther )
    : BaseType( rOther )
    {}

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo2dBezier3()
    {}

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     * @see Clone
     * @see ClonePoints
     */
    Geo2dBezier3& operator=( const Geo2dBezier3& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Geo2dBezier3& operator=( Geo2dBezier3<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        Geo2dBezier3::Pointer pNewGeom = Geo2dBezier3::Pointer( new Geo2dBezier3( ThisPoints ) );
        ValuesContainerType DummyKnots;
        if (BaseType::mpBezierGeometryData != NULL)
        {
            pNewGeom->AssignGeometryData(DummyKnots, DummyKnots, DummyKnots,
                BaseType::mCtrlWeights, BaseType::mExtractionOperator, BaseType::mOrder1, BaseType::mOrder2, 0,
                static_cast<int>(BaseType::mpBezierGeometryData->DefaultIntegrationMethod()) + 1);
        }
        return pNewGeom;
    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
////        Geometry< Point<3> >::PointsArrayType NewPoints;
////        //making a copy of the nodes TO POINTS (not Nodes!!!)

////        for ( IndexType i = 0; i < this->Points().size(); ++i )
////        NewPoints.push_back( this->Points()[i] );

////        //creating a geometry with the new points
////        boost::shared_ptr< Geometry< Point<3> > >
////        p_clone( new Geo2dBezier3< Point<3> >( NewPoints ) );

////        p_clone->ClonePoints();

////        return p_clone;

//        KRATOS_THROW_ERROR(std::logic_error, "NURBS geometry does not support for Clone", *this)
//    }

    /**
     * Informations
     */

    virtual GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_NURBS;
    }

    virtual GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Bezier2D3;
    }

    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    // TODO

    /**
     * Jacobian
     */

    /**
     * Jacobians for given  method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return JacobiansType a VectorType of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        BaseType::CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Jacobians for given  method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return JacobiansType a VectorType of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        BaseType::CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Determinants of Jacobian for given  method.
     */
    virtual VectorType& DeterminantOfJacobian( VectorType& rResults,
            IntegrationMethod ThisMethod ) const
    {
        JacobiansType J;
        J = Jacobian( J, ThisMethod );

        if ( rResults.size() != J.size() )
        {
            rResults.resize(J.size(), false);
        }

        MatrixType JtJ(2, 2);
        for ( unsigned int pnt = 0; pnt < J.size(); ++pnt )
        {
            noalias(JtJ) = prod( trans(J[pnt]), J[pnt] );
            rResults[pnt] = sqrt(MathUtils<double>::Det(JtJ));
        }

        return rResults;
    }

    virtual JacobiansType& Jacobian0( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        BaseType::CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Compute the Bezier control points
     */
    virtual void ExtractLocalCoordinates(PointsArrayType& rPoints)
    {
        std::size_t number_of_points = this->PointsNumber();
        std::size_t number_of_local_points = BaseType::mExtractionOperator.size2();
        rPoints.clear();
        rPoints.reserve(number_of_local_points);

        // compute the Bezier weight
        VectorType bezier_weights = prod(trans(BaseType::mExtractionOperator), BaseType::mCtrlWeights);

        // compute the Bezier control points
        typedef typename PointType::Pointer PointPointerType;
        for(std::size_t i = 0; i < number_of_local_points; ++i)
        {
            PointPointerType pPoint = PointPointerType(new PointType(0, 0.0, 0.0, 0.0));
            for(std::size_t j = 0; j < number_of_points; ++j)
                noalias(*pPoint) += BaseType::mExtractionOperator(j, i) * this->GetPoint(j) * BaseType::mCtrlWeights[j] / bezier_weights[i];
            rPoints.push_back(pPoint);
        }
    }

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "2 dimensional Bezier decomposition surface in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << Info();
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
//        MatrixType jacobian;
//        Jacobian( jacobian, PointType() );
//        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * TO BE CALLED BY ELEMENT
     * TODO: optimized this by integrating pre-computed values at Gauss points
     */
    virtual void AssignGeometryData(
        const ValuesContainerType& Knots1, //not used
        const ValuesContainerType& Knots2, //not used
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        const int& Degree1,
        const int& Degree2,
        const int& Degree3, //not used
        const int& NumberOfIntegrationMethod
    )
    {
        BaseType::mCtrlWeights = Weights;
        BaseType::mOrder1 = Degree1;
        BaseType::mOrder2 = Degree2;
        BaseType::mNumber1 = BaseType::mOrder1 + 1;
        BaseType::mNumber2 = BaseType::mOrder2 + 1;

        // select the type of extraction operator to save in memory
        if(ExtractionOperator.size1() == 2 && ExtractionOperator.size2() != 2)
        // extraction operator is stored as compressed matrix
        {
            unsigned int size_ex_n = (unsigned int)(ExtractionOperator(0, 0) - 1);
            unsigned int size_ex_nz = ExtractionOperator.size2() - 1;
            if( ( (double)(size_ex_nz) ) / (size_ex_n * size_ex_n) < 0.2 )
            {
                BaseType::mExtractionOperator = IsogeometricMathUtils::MCSR2CSR(ExtractionOperator);
            }
            else
                BaseType::mExtractionOperator = IsogeometricMathUtils::MCSR2MAT(ExtractionOperator);
        }
        else if((ExtractionOperator.size1() != 2) && (ExtractionOperator.size1() == ExtractionOperator.size2()))
        // extraction operator is stored as full matrix
            BaseType::mExtractionOperator = ExtractionOperator;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid extraction operator", __FUNCTION__)

        // size checking
        if(BaseType::mNumber1 * BaseType::mNumber2 != this->size())
        {
            KRATOS_THROW_ERROR(std::logic_error, "The parametric parameters is not compatible.", __FUNCTION__)
        }

        // find the existing integration rule or create new one if not existed
        BezierUtils::RegisterIntegrationRule<2, 3, 2>(NumberOfIntegrationMethod, Degree1, Degree2);

        // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
        BaseType::mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<2, 3, 2>(NumberOfIntegrationMethod, Degree1, Degree2);
        BaseType::BaseType::mpGeometryData = &(*BaseType::mpBezierGeometryData);
    }

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */
//    static const GeometryData msGeometryData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    /**
     * Private Operations
     */

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo2dBezier3;

    /**
     * Un accessible methods
     */

};    // Class Geo2dBezier3

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo2dBezier3<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo2dBezier3<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/**
 * static initialisation for geometry_data
 */
// template<class TPointType>
// const GeometryData Geo2dBezier3<TPointType>::msGeometryData
// (
//    2,
//    3,
//    2,
//    GeometryData::GI_GAUSS_2,
//    IntegrationPointsContainerType(),
//    ShapeFunctionsValuesContainerType(),
//    ShapeFunctionsLocalGradientsContainerType()
// );

}    // namespace Kratos.

#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
#undef DEBUG_LEVEL3
#undef DEBUG_LEVEL4
#undef DEBUG_LEVEL5
#undef DEBUG_LEVEL6
#undef DEBUG_LEVEL7
#undef DEBUG_LEVEL8
#undef ENABLE_PROFILING

#endif

