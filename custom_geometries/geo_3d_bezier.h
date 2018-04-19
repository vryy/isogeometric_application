/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Feb 3 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_3D_BEZIER_H_INCLUDED )
#define  KRATOS_GEO_3D_BEZIER_H_INCLUDED

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
//#include "integration/quadrature.h"
//#include "integration/line_gauss_legendre_integration_points.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
// #define DEBUG_LEVEL3
#define ENABLE_PROFILING
#define ENABLE_CHECK_SIZE

namespace Kratos
{
/**
 * A geometry representing NURBS curve
 */

template<class TPointType>
class Geo3dBezier : public IsogeometricGeometry<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef IsogeometricGeometry<TPointType> BaseType;

    /**
     * Pointer definition of Geo3dBezier
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo3dBezier );

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
    typedef typename BaseType::PointType PointType;

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

    Geo3dBezier()
//    : BaseType( PointsArrayType(), &msGeometryData )
    : BaseType( PointsArrayType() ), mpBezierGeometryData(NULL)
    {}

    Geo3dBezier( const PointsArrayType& ThisPoints )
//    : BaseType( ThisPoints, &msGeometryData )
    : BaseType( ThisPoints ), mpBezierGeometryData(NULL)
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
    Geo3dBezier( Geo3dBezier const& rOther )
    : BaseType( rOther ), mpBezierGeometryData(NULL)
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
    template<class TOtherPointType> Geo3dBezier( Geo3dBezier<TOtherPointType> const& rOther )
    : BaseType( rOther ), mpBezierGeometryData(NULL)
    {}

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo3dBezier()
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
    Geo3dBezier& operator=( const Geo3dBezier& rOther )
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
    Geo3dBezier& operator=( Geo3dBezier<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        Geo3dBezier::Pointer pNewGeom = Geo3dBezier::Pointer( new Geo3dBezier( ThisPoints ) );
        if (mpBezierGeometryData != NULL)
        {
            ValuesContainerType DummyKnots;
            pNewGeom->AssignGeometryData(DummyKnots, DummyKnots, DummyKnots,
                mCtrlWeights, mExtractionOperator, mOrder1, mOrder2, mOrder3,
                static_cast<int>(mpBezierGeometryData->DefaultIntegrationMethod()) + 1);
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
////        p_clone( new Geo3dBezier< Point<3> >( NewPoints ) );

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
        return GeometryData::Kratos_Bezier3D;
    }

    /**
     * This method calculates and returns Length or characteristic
     * length of this geometry depending on it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */

    /**
     * Compute shape function values and local gradients at every integration points of an integration method.
     */
    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        IntegrationMethod ThisMethod
    ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

//        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber(ThisMethod);
        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

        shape_functions_values.resize(NumberOfIntegrationPoints, this->PointsNumber(), false);

        shape_functions_local_gradients.resize(NumberOfIntegrationPoints);
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType(this->PointsNumber(), 3));

        #ifdef DEBUG_LEVEL3
        KRATOS_WATCH(NumberOfIntegrationPoints)
        #endif

        const MatrixType& bezier_functions_values
            = mpBezierGeometryData->ShapeFunctionsValues( ThisMethod );

        const ShapeFunctionsGradientsType& bezier_functions_local_gradients
            = mpBezierGeometryData->ShapeFunctionsLocalGradients( ThisMethod );

        VectorType temp_bezier_values(bezier_functions_values.size2());
        VectorType bezier_weights(mExtractionOperator.size2());
        double denom, tmp1, tmp2, tmp3;
        VectorType tmp_gradients1(this->PointsNumber());
        VectorType tmp_gradients2(this->PointsNumber());
        VectorType tmp_gradients3(this->PointsNumber());
        for(IndexType i = 0; i < NumberOfIntegrationPoints; ++i)
        {
            noalias(temp_bezier_values) = row(bezier_functions_values, i);

            //compute the Bezier weight
            noalias(bezier_weights) = prod(trans(mExtractionOperator), mCtrlWeights);
            denom = inner_prod(temp_bezier_values, bezier_weights);

            //compute the shape function values
            VectorType temp_values = prod(mExtractionOperator, temp_bezier_values);
            for(IndexType j = 0; j < this->PointsNumber(); ++j)
                shape_functions_values(i, j) = (temp_values(j) * mCtrlWeights(j) / denom);

            //compute the shape function local gradients
//            shape_functions_local_gradients[i].resize(this->PointsNumber(), 3, false);
            tmp1 = inner_prod(row(bezier_functions_local_gradients[i], 0), bezier_weights);
            tmp2 = inner_prod(row(bezier_functions_local_gradients[i], 1), bezier_weights);
            tmp3 = inner_prod(row(bezier_functions_local_gradients[i], 2), bezier_weights);

            noalias(tmp_gradients1) = prod( mExtractionOperator,
                        (1 / denom) * row(bezier_functions_local_gradients[i], 0) - (tmp1 / pow(denom, 2)) * temp_bezier_values );

            noalias(tmp_gradients2) = prod(mExtractionOperator,
                        (1 / denom) * row(bezier_functions_local_gradients[i], 1) - (tmp2 / pow(denom, 2)) * temp_bezier_values );

            noalias(tmp_gradients3) = prod(mExtractionOperator,
                        (1 / denom) * row(bezier_functions_local_gradients[i], 2) - (tmp3 / pow(denom, 2)) * temp_bezier_values );

            for(IndexType j = 0; j < this->PointsNumber(); ++j)
            {
                shape_functions_local_gradients[i](j, 0) = tmp_gradients1(j) * mCtrlWeights(j);
                shape_functions_local_gradients[i](j, 1) = tmp_gradients2(j) * mCtrlWeights(j);
                shape_functions_local_gradients[i](j, 2) = tmp_gradients3(j) * mCtrlWeights(j);
            }
        }
    }

    /**
     * Compute Jacobian at every integration points for an integration method
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting local gradients of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

//        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );
        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 3 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 0, 2 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 2 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 2 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Compute Jacobian at every integration points for an integration method
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting local gradients of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

//        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );
        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 3 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 0, 2 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 2 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 2 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Compute Jacobian at an integration point for an integration method
     */
    virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting local gradients of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

        rResult.resize(3, 3, false);
        noalias(rResult) = ZeroMatrix( 3, 3 );

        //loop over all nodes
        for ( IndexType i = 0; i < this->PointsNumber(); ++i )
        {
            rResult( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 0, 2 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
            rResult( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 1, 2 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
            rResult( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 2, 1 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 2, 2 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
        }

        return rResult;
    }

    /**
     * Compute Jacobian at an integration point for an integration method
     */
    virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting local gradients of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

        rResult.resize(3, 3, false);
        noalias(rResult) = ZeroMatrix( 3, 3 );

        //loop over all nodes
        for ( IndexType i = 0; i < this->PointsNumber(); ++i )
        {
            rResult( 0, 0 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 0, 1 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 0, 2 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
            rResult( 1, 0 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 1, 1 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 1, 2 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
            rResult( 2, 0 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 2, 1 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 2, 2 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients[IntegrationPointIndex]( i, 2 ) );
        }

        return rResult;
    }

    /**
     * Compute Jacobian at a particular point
     */
    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
        VectorType shape_functions_values;
        MatrixType shape_functions_local_gradients;

        //getting local gradients of shape functions
        ShapeFunctionsValuesAndLocalGradients(shape_functions_values, shape_functions_local_gradients, rPoint);

        rResult.resize( 3, 3, false );
        noalias(rResult) = ZeroMatrix( 3, 3 );

        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            rResult( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_local_gradients( i, 2 ) );
        }

        return rResult;
    }

    /**
     * Compute Jacobian at a particular point
     */
    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint, Matrix& DeltaPosition ) const
    {
        VectorType shape_functions_values;
        MatrixType shape_functions_local_gradients;

        //getting local gradients of shape functions
        ShapeFunctionsValuesAndLocalGradients(shape_functions_values, shape_functions_local_gradients, rPoint);

        rResult.resize( 3, 3, false );
        noalias(rResult) = ZeroMatrix( 3, 3 );

        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            rResult( 0, 0 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( this->GetPoint( i ).Z() - DeltaPosition(i, 2) ) * ( shape_functions_local_gradients( i, 2 ) );
        }

        return rResult;
    }

    // REMARKS: Those InverseOfJacobian, DeterminantOfJacobian below is implemented in abstract Geometry class
//    /**
//     * Compute inverse of Jacobian at an integration point for an integration method
//     */
//    virtual Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
//    {
//        //current jacobian
//        MatrixType tempMatrix;
//        tempMatrix = Jacobian( tempMatrix, IntegrationPointIndex, ThisMethod );
//        double det = 0.0;
//        //inverse of jacobian
//        MathUtils<double>::InvertMatrix3( tempMatrix, rResult, det );
//        return rResult;
//    }
//
//    /**
//     * Compute determinant of Jacobian at an integration point for an integration method
//     */
//    virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
//    {
//        //current jacobian
//        MatrixType tempMatrix;
//        tempMatrix = Jacobian( tempMatrix, IntegrationPointIndex, ThisMethod );
//        return MathUtils<double>::Det( tempMatrix );
//    }
//
//    /**
//     * Compute inverse of Jacobian at every integration points for an integration method
//     */
//    virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
//    {
////        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );
//        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

//        if ( rResult.size() != NumberOfIntegrationPoints )
//        {
//            JacobiansType temp( NumberOfIntegrationPoints );
//            rResult.swap( temp );
//        }

//        //loop over all integration points
//        for ( unsigned int pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
//        {
//            Matrix tempMatrix = ZeroMatrix( 3, 3 );
//            rResult[pnt] = InverseOfJacobian( tempMatrix, pnt, ThisMethod );
//        }

//        return rResult;
//    }
//
//    /**
//     * Compute determinant of Jacobian at every integration points for an integration method
//     */
//    virtual VectorType& DeterminantOfJacobian( VectorType& rResult, IntegrationMethod ThisMethod ) const
//    {
////        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );
//        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

//        if ( rResult.size() != NumberOfIntegrationPoints )
//        {
//            rResult.resize( NumberOfIntegrationPoints );
//        }

//        //loop over all integration points
//        for ( unsigned int pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
//        {
//            rResult[pnt] = DeterminantOfJacobian( pnt, ThisMethod );
//        }

//        return rResult;
//    }
//
//    /**
//     * Compute inverse of Jacobian at a particular point
//     */
//    virtual Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
//    {
//        //current jacobian
//        Matrix tempMatrix = ZeroMatrix( 3, 3 );
//        tempMatrix = Jacobian( tempMatrix, rPoint );

//        //setting up result matrix
//        rResult.resize( 3, 3, false );
//        double det;
//        MathUtils<double>::InvertMatrix3( tempMatrix, rResult, det );
//        return rResult;
//    }
//
//    /**
//     * Compute determinant of Jacobian at a particular point
//     */
//    virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
//    {
//        //current jacobian
//        Matrix tempMatrix = ZeroMatrix( 3, 3 );
//        tempMatrix = Jacobian( tempMatrix, rPoint );
//        return MathUtils<double>::Det( tempMatrix );
//    }

    /**
     * Compute Jacobian at the reference configuration
     */
    virtual JacobiansType& Jacobian0( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            shape_functions_values,
            shape_functions_local_gradients,
            ThisMethod
        );

//        SizeType NumberOfIntegrationPoints = this->IntegrationPointsNumber( ThisMethod );
        SizeType NumberOfIntegrationPoints = mpBezierGeometryData->IntegrationPoints(ThisMethod).size();

        if ( rResult.size() != NumberOfIntegrationPoints )
        {
            JacobiansType temp( NumberOfIntegrationPoints );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( IndexType pnt = 0; pnt < NumberOfIntegrationPoints; ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 3, 3 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 0, 2 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 2 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( this->GetPoint( i ).Z0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 2, 2 ) += ( this->GetPoint( i ).Z0() ) * ( shape_functions_local_gradients[pnt]( i, 2 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Compute shape function values at a particular reference point. This function is kept to keep the backward compatibility. The function ShapeFunctionsValuesAndLocalGradients is more general and direct to use.
     */
    virtual Vector& ShapeFunctionsValues( Vector& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(mNumber1);
        VectorType bezier_functions_values2(mNumber2);
        VectorType bezier_functions_values3(mNumber3);
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        VectorType bezier_functions_derivatives3(mNumber3);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rCoordinates[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rCoordinates[1]);
        BezierUtils::bernstein(bezier_functions_values3, bezier_functions_derivatives3, mOrder3, rCoordinates[2]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;
                    bezier_functions_values(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                }
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        if(rResults.size() != this->PointsNumber())
            rResults.resize(this->PointsNumber(), false);
        noalias( rResults ) = prod(mExtractionOperator, bezier_functions_values);
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
            rResults(i) *= (mCtrlWeights(i) / denom);

        return rResults;
    }

    /**
     * Compute shape function local gradients at a particular reference point. This function is kept to keep the backward compatibility. The function ShapeFunctionsValuesAndLocalGradients is more general and direct to use.
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(mNumber1);
        VectorType bezier_functions_values2(mNumber2);
        VectorType bezier_functions_values3(mNumber3);
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        VectorType bezier_functions_derivatives3(mNumber3);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rCoordinates[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rCoordinates[1]);
        BezierUtils::bernstein(bezier_functions_values3, bezier_functions_derivatives3, mOrder3, rCoordinates[2]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;
                    bezier_functions_values(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                }
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        VectorType bezier_functions_local_derivatives1(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_derivatives2(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_derivatives3(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;

                    bezier_functions_local_derivatives1(index) =
                        bezier_functions_derivatives1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_derivatives2(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_derivatives2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_derivatives3(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_derivatives3(k);
                }
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local gradients
        rResults.resize(this->PointsNumber(), 3, false);
        double tmp1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double tmp2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
        double tmp3 = inner_prod(bezier_functions_local_derivatives3, bezier_weights);
        VectorType tmp_gradients1 =
            prod(
                mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives1 -
                        (tmp1 / pow(denom, 2)) * bezier_functions_values
            );
        VectorType tmp_gradients2 =
            prod(
                mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives2 -
                        (tmp2 / pow(denom, 2)) * bezier_functions_values
            );
        VectorType tmp_gradients3 =
            prod(
                mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives3 -
                        (tmp3 / pow(denom, 2)) * bezier_functions_values
            );
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults(i, 0) = tmp_gradients1(i) * mCtrlWeights(i);
            rResults(i, 1) = tmp_gradients2(i) * mCtrlWeights(i);
            rResults(i, 2) = tmp_gradients3(i) * mCtrlWeights(i);
        }

        return rResults;
    }

    /**
     * Compute shape function second derivatives at a particular reference point. This function is kept to keep the backward compatibility.
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rCoordinates ) const
    {
        KRATOS_THROW_ERROR( std::logic_error,
                      "Calling base class ShapeFunctionsSecondDerivatives at" , *this );
        return rResult;
    }

    /**
     * Compute the Bezier control points
     */
    virtual void ExtractLocalCoordinates(PointsArrayType& rPoints)
    {
        std::size_t number_of_points = this->PointsNumber();
        std::size_t number_of_local_points = mExtractionOperator.size2();
        rPoints.clear();
        rPoints.reserve(number_of_local_points);

        // compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);

        // compute the Bezier control points
        typedef typename PointType::Pointer PointPointerType;
        for(std::size_t i = 0; i < number_of_local_points; ++i)
        {
            PointPointerType pPoint = PointPointerType(new PointType(0, 0.0, 0.0, 0.0));
            for(std::size_t j = 0; j < number_of_points; ++j)
                noalias(*pPoint) += mExtractionOperator(j, i) * this->GetPoint(j) * mCtrlWeights[j] / bezier_weights[i];
            rPoints.push_back(pPoint);
        }
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        double tol = 1.0e-6;
        if ( (rResult[0] > -tol) && (rResult[0] < 1 + tol) )
            if ( (rResult[1] > -tol) && (rResult[1] < 1 + tol) )
                if ( (rResult[2] > -tol) && (rResult[2] < 1 + tol) )
                    return true;

        return false;
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
        return "3 dimensional Bezier decomposition volume in 3D space";
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
        rOStream << std::endl;
    }

    virtual void AssignGeometryData(
        const ValuesContainerType& Knots1, //not used
        const ValuesContainerType& Knots2, //not used
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        const int& Degree1,
        const int& Degree2,
        const int& Degree3,
        const int& NumberOfIntegrationMethod
    )
    {
        mCtrlWeights = Weights;
        mOrder1 = Degree1;
        mOrder2 = Degree2;
        mOrder3 = Degree3;
        mNumber1 = mOrder1 + 1;
        mNumber2 = mOrder2 + 1;
        mNumber3 = mOrder3 + 1;
        mExtractionOperator = ExtractionOperator;

        // size checking
        if(mExtractionOperator.size1() != this->PointsNumber())
        {
            KRATOS_WATCH(this->PointsNumber())
            KRATOS_WATCH(mExtractionOperator)
            KRATOS_THROW_ERROR(std::logic_error, "The number of row of extraction operator must be equal to number of nodes", __FUNCTION__)
        }
        if(mExtractionOperator.size2() != (mOrder1 + 1) * (mOrder2 + 1) * (mOrder3 + 1))
        {
            KRATOS_WATCH(mExtractionOperator)
            KRATOS_WATCH(mOrder1)
            KRATOS_WATCH(mOrder2)
            KRATOS_WATCH(mOrder3)
            KRATOS_THROW_ERROR(std::logic_error, "The number of column of extraction operator must be equal to (p_u+1) * (p_v+1) * (p_w+1), error at", __FUNCTION__)
        }

        if(NumberOfIntegrationMethod > 0)
        {
            // find the existing integration rule or create new one if not existed
            BezierUtils::RegisterIntegrationRule<3, 3, 3>(NumberOfIntegrationMethod, Degree1, Degree2, Degree3);

            // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
            mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<3, 3, 3>(NumberOfIntegrationMethod, Degree1, Degree2, Degree3);
            #ifndef ENABLE_PRECOMPUTE
            BaseType::mpGeometryData = &(*mpBezierGeometryData);
            #else
            // precompute the values at each integration points; note that it can generate a LOT of data
            IntegrationPointsContainerType all_integration_points
                    = BezierUtils::AllIntegrationPoints(NumberOfIntegrationMethod, mOrder1, mOrder2, mOrder3);

            ShapeFunctionsValuesContainerType shape_functions_values;
            ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

            for(IndexType i = 0; i < NumberOfIntegrationMethod; ++i)
            {
                CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                    shape_functions_values[i],
                    shape_functions_local_gradients[i],
                    (GeometryData::IntegrationMethod)i
                );
            }

            mpGeometryData = GeometryData::Pointer(
                new GeometryData(
                        3,
                        3,
                        3,
                        GeometryData::GI_GAUSS_1,           //ThisDefaultMethod
                        all_integration_points,             //ThisIntegrationPoints
                        shape_functions_values,             //ThisShapeFunctionsValues
                        shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                    )
                );

            BaseType::mpGeometryData = &(*mpGeometryData);
            #endif
        }
    }

protected:

    /**
     * there are no protected class members
     */

    GeometryData::Pointer mpBezierGeometryData;
    #ifdef ENABLE_PRECOMPUTE
    GeometryData::Pointer mpGeometryData;
    #endif

    MatrixType mExtractionOperator;

    ValuesContainerType mCtrlWeights; //weight of control points

    int mOrder1; //order of the surface at parametric direction 1
    int mOrder2; //order of the surface at parametric direction 2
    int mOrder3; //order of the surface at parametric direction 3

    int mNumber1; //number of bezier shape functions define the surface on parametric direction 1
    int mNumber2; //number of bezier shape functions define the surface on parametric direction 2
    int mNumber3; //number of bezier shape functions define the surface on parametric direction 3

private:

    /**
     * Static Member Variables
     */
//    static const GeometryData msGeometryData; // see COMMENTS below

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
     * Calculate shape function values and local gradient at a particular point
     */
    virtual void ShapeFunctionsValuesAndLocalGradients(
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(mNumber1);
        VectorType bezier_functions_values2(mNumber2);
        VectorType bezier_functions_values3(mNumber3);
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        VectorType bezier_functions_derivatives3(mNumber3);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rPoint[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rPoint[1]);
        BezierUtils::bernstein(bezier_functions_values3, bezier_functions_derivatives3, mOrder3, rPoint[2]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;
                    bezier_functions_values(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                }
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        VectorType bezier_functions_local_derivatives1(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_derivatives2(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_derivatives3(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;

                    bezier_functions_local_derivatives1(index) =
                        bezier_functions_derivatives1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_derivatives2(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_derivatives2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_derivatives3(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_derivatives3(k);
                }
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        if(shape_functions_values.size() != this->PointsNumber())
            shape_functions_values.resize(this->PointsNumber(), false);
        noalias( shape_functions_values ) = prod(mExtractionOperator, bezier_functions_values);
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
            shape_functions_values(i) *= (mCtrlWeights(i) / denom);

        //compute the shape function local gradients
        if(shape_functions_local_gradients.size1() != this->PointsNumber()
            || shape_functions_local_gradients.size2() != 3)
            shape_functions_local_gradients.resize(this->PointsNumber(), 3, false);
        double tmp1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double tmp2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
        double tmp3 = inner_prod(bezier_functions_local_derivatives3, bezier_weights);
        VectorType tmp_gradients1 = prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives1 - (tmp1 / pow(denom, 2)) * bezier_functions_values );
        VectorType tmp_gradients2 = prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives2 - (tmp2 / pow(denom, 2)) * bezier_functions_values );
        VectorType tmp_gradients3 = prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives3 - (tmp3 / pow(denom, 2)) * bezier_functions_values );
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            shape_functions_local_gradients(i, 0) = tmp_gradients1(i) * mCtrlWeights(i);
            shape_functions_local_gradients(i, 1) = tmp_gradients2(i) * mCtrlWeights(i);
            shape_functions_local_gradients(i, 2) = tmp_gradients3(i) * mCtrlWeights(i);
        }
    }

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo3dBezier;

    /**
     * Un accessible methods
     */

};    // Class Geo3dBezier

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo3dBezier<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo3dBezier<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/**
 * static initialisation for geometry_data
 * TODO: to be deleted; the geometry_data is obtained by register with the BezierUtils
 * the static initialization is not applicable for Bezier geometry since each geometry are not the same
 */
// template<class TPointType>
// const GeometryData Geo3dBezier<TPointType>::msGeometryData
// (
//    3,
//    3,
//    3,
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
#undef ENABLE_PRECOMPUTE

#endif

