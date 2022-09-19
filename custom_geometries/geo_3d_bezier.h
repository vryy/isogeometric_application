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
 * Geometry representing Bezier decomposition of a NURBS volume
 */
template<class TPointType>
class Geo3dBezier : public IsogeometricGeometry<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Pointer definition of Geo3dBezier
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo3dBezier );

    /**
     * IsogeometricGeometry as base class.
     */
    typedef IsogeometricGeometry<TPointType> BaseType;

    /**
     * The original geometry type
     */
    typedef typename BaseType::GeometryType GeometryType;

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
     * A first order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A second order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /**
     * A first order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A second order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /** A fourth order tensor to hold shape functions' local third order derivatives
     */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

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
    typedef boost::numeric::ublas::compressed_matrix<typename MatrixType::value_type> CompressedMatrixType;

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
    : BaseType( PointsArrayType() ), mpBezierGeometryData(NULL)
    {}

    Geo3dBezier( const PointsArrayType& ThisPoints )
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
    : BaseType( rOther )
    , mpBezierGeometryData(rOther.mpBezierGeometryData)
    , mOrder1(rOther.mOrder1)
    , mOrder2(rOther.mOrder2)
    , mOrder3(rOther.mOrder3)
    , mNumber1(rOther.mNumber1)
    , mNumber2(rOther.mNumber2)
    , mNumber3(rOther.mNumber3)
    , mExtractionOperator(rOther.mExtractionOperator)
    , mCtrlWeights(rOther.mCtrlWeights)
    {
        GeometryType::mpGeometryData = &(*mpBezierGeometryData);
    }

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
    : IsogeometricGeometry<TOtherPointType>( rOther )
    , mpBezierGeometryData(rOther.mpBezierGeometryData)
    , mOrder1(rOther.mOrder1)
    , mOrder2(rOther.mOrder2)
    , mOrder3(rOther.mOrder3)
    , mNumber1(rOther.mNumber1)
    , mNumber2(rOther.mNumber2)
    , mNumber3(rOther.mNumber3)
    , mExtractionOperator(rOther.mExtractionOperator)
    , mCtrlWeights(rOther.mCtrlWeights)
    {
        Geometry<TOtherPointType>::mpGeometryData = &(*mpBezierGeometryData);
    }

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
        this->mpBezierGeometryData = rOther.mpBezierGeometryData;
        GeometryType::mpGeometryData = &(*(this->mpBezierGeometryData));
        this->mOrder1 = rOther.mOrder1;
        this->mOrder2 = rOther.mOrder2;
        this->mOrder3 = rOther.mOrder3;
        this->mNumber1 = rOther.mNumber1;
        this->mNumber2 = rOther.mNumber2;
        this->mNumber3 = rOther.mNumber3;
        this->mExtractionOperator = rOther.mExtractionOperator;
        this->mCtrlWeights = rOther.mCtrlWeights;
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
        IsogeometricGeometry<TOtherPointType>::operator=( rOther );
        this->mpBezierGeometryData = rOther.mpBezierGeometryData;
        Geometry<TOtherPointType>::mpGeometryData = &(*(this->mpBezierGeometryData));
        this->mOrder1 = rOther.mOrder1;
        this->mOrder2 = rOther.mOrder2;
        this->mOrder3 = rOther.mOrder3;
        this->mNumber1 = rOther.mNumber1;
        this->mNumber2 = rOther.mNumber2;
        this->mNumber3 = rOther.mNumber3;
        this->mExtractionOperator = rOther.mExtractionOperator;
        this->mCtrlWeights = rOther.mCtrlWeights;
        return *this;
    }

    /**
     * Operations
     */

    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const override
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

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        #ifdef SD_APP_FORWARD_COMPATIBILITY
        return static_cast<GeometryData::KratosGeometryType>(IsogeometricGeometryData::Kratos_Bezier3D);
        #else
        return GeometryData::Kratos_Bezier3D;
        #endif
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
    void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        IntegrationMethod ThisMethod
    ) const final
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
        VectorType bezier_weights(mNumber1 * mNumber2 * mNumber3);
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
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const final
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
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const final
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
    Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const final
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
    Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const final
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
    Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const final
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
    Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint, Matrix& DeltaPosition ) const final
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
//    Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const final
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
//    double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const final
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
//    JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const final
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
//    VectorType& DeterminantOfJacobian( VectorType& rResult, IntegrationMethod ThisMethod ) const final
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
//    Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const final
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
//    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const final
//    {
//        //current jacobian
//        Matrix tempMatrix = ZeroMatrix( 3, 3 );
//        tempMatrix = Jacobian( tempMatrix, rPoint );
//        return MathUtils<double>::Det( tempMatrix );
//    }

    /**
     * Compute Jacobian at the reference configuration
     */
    JacobiansType& Jacobian0( JacobiansType& rResult, IntegrationMethod ThisMethod ) const final
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
    Vector& ShapeFunctionsValues( Vector& rResults, const CoordinatesArrayType& rCoordinates ) const final
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
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResults, const CoordinatesArrayType& rCoordinates ) const final
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
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResults, const CoordinatesArrayType& rCoordinates ) const final
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
        VectorType bezier_functions_second_derivatives1(mNumber1);
        VectorType bezier_functions_second_derivatives2(mNumber2);
        VectorType bezier_functions_second_derivatives3(mNumber3);
        BezierUtils::bernstein(bezier_functions_values1,
                               bezier_functions_derivatives1,
                               bezier_functions_second_derivatives1,
                               mOrder1,
                               rCoordinates[0]);
        BezierUtils::bernstein(bezier_functions_values2,
                               bezier_functions_derivatives2,
                               bezier_functions_second_derivatives2,
                               mOrder2,
                               rCoordinates[1]);
        BezierUtils::bernstein(bezier_functions_values3,
                               bezier_functions_derivatives3,
                               bezier_functions_second_derivatives3,
                               mOrder3,
                               rCoordinates[2]);

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

        //compute trivariate Bezier shape functions second derivatives w.r.t local coordinates
        VectorType bezier_functions_local_second_derivatives11(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_second_derivatives12(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_second_derivatives13(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_second_derivatives22(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_second_derivatives23(mNumber1 * mNumber2 * mNumber3);
        VectorType bezier_functions_local_second_derivatives33(mNumber1 * mNumber2 * mNumber3);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                for(IndexType k = 0; k < mNumber3; ++k)
                {
                    IndexType index = k + (j + i * mNumber2) * mNumber3;

                    bezier_functions_local_second_derivatives11(index) =
                        bezier_functions_second_derivatives1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_second_derivatives12(index) =
                        bezier_functions_derivatives1(i) *
                        bezier_functions_derivatives2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_second_derivatives13(index) =
                        bezier_functions_derivatives1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_derivatives3(k);

                    bezier_functions_local_second_derivatives22(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_second_derivatives2(j) *
                        bezier_functions_values3(k);
                    bezier_functions_local_second_derivatives23(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_derivatives2(j) *
                        bezier_functions_derivatives3(k);

                    bezier_functions_local_second_derivatives33(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_second_derivatives3(k);
                }
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local second gradients
        rResults.resize(this->PointsNumber(), false);
        double aux1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double aux2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
        double aux3 = inner_prod(bezier_functions_local_derivatives3, bezier_weights);
        double auxs11 = inner_prod(bezier_functions_local_second_derivatives11, bezier_weights);
        double auxs12 = inner_prod(bezier_functions_local_second_derivatives12, bezier_weights);
        double auxs13 = inner_prod(bezier_functions_local_second_derivatives13, bezier_weights);
        double auxs22 = inner_prod(bezier_functions_local_second_derivatives22, bezier_weights);
        double auxs23 = inner_prod(bezier_functions_local_second_derivatives23, bezier_weights);
        double auxs33 = inner_prod(bezier_functions_local_second_derivatives33, bezier_weights);
        VectorType tmp_gradients11 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives11
                    - (aux1 / pow(denom, 2)) * bezier_functions_local_derivatives1 * 2
                    - (auxs11 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * pow(aux1, 2) / pow(denom, 3) * bezier_functions_values
            );
        VectorType tmp_gradients12 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives12
                    - ((aux1 + aux2) / pow(denom, 2)) * bezier_functions_local_derivatives1
                    - (auxs12 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * aux1 * aux2 / pow(denom, 3) * bezier_functions_values
            );
        VectorType tmp_gradients13 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives13
                    - ((aux1 + aux3) / pow(denom, 2)) * bezier_functions_local_derivatives1
                    - (auxs13 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * aux1 * aux3 / pow(denom, 3) * bezier_functions_values
            );
        VectorType tmp_gradients22 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives22
                    - (aux2 / pow(denom, 2)) * bezier_functions_local_derivatives2 * 2
                    - (auxs22 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * pow(aux2, 2) / pow(denom, 3) * bezier_functions_values
            );
        VectorType tmp_gradients23 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives23
                    - ((aux2 + aux3) / pow(denom, 2)) * bezier_functions_local_derivatives2
                    - (auxs23 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * aux2 * aux3 / pow(denom, 3) * bezier_functions_values
            );
        VectorType tmp_gradients33 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives33
                    - (aux3 / pow(denom, 2)) * bezier_functions_local_derivatives3 * 2
                    - (auxs33 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * pow(aux3, 2) / pow(denom, 3) * bezier_functions_values
            );
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults[i].resize(3, 3, false);
            rResults[i](0, 0) = tmp_gradients11(i) * mCtrlWeights(i);
            rResults[i](0, 1) = tmp_gradients12(i) * mCtrlWeights(i);
            rResults[i](0, 2) = tmp_gradients13(i) * mCtrlWeights(i);
            rResults[i](1, 0) = rResults[i](0, 1);
            rResults[i](1, 1) = tmp_gradients22(i) * mCtrlWeights(i);
            rResults[i](1, 2) = tmp_gradients23(i) * mCtrlWeights(i);
            rResults[i](2, 0) = rResults[i](0, 2);
            rResults[i](2, 1) = rResults[i](1, 2);
            rResults[i](2, 2) = tmp_gradients33(i) * mCtrlWeights(i);
        }

        return rResults;
    }

    /**
     * Compute shape function third derivatives at a particular reference point.
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResults, const CoordinatesArrayType& rPoint ) const final
    {
        // TODO

        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not yet implemented")

        rResults.resize(this->PointsNumber(), false);

        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults[i].resize(3, false);
            for(IndexType j = 0; j < 3; ++j)
            {
                rResults[i][j].resize(3, 3, false);
                noalias(rResults[i][j]) = ZeroMatrix(3, 3);
            }
        }

        return rResults;
    }

    /**
     * Compute the Bezier control points
     */
    void ExtractControlPoints(PointsArrayType& rPoints) final
    {
        std::size_t number_of_points = this->PointsNumber();
        std::size_t number_of_local_points = mNumber1 * mNumber2 * mNumber3;
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
                noalias(*pPoint) += mExtractionOperator(j, i) * this->GetPoint(j).GetInitialPosition() * mCtrlWeights[j] / bezier_weights[i];
            pPoint->SetInitialPosition(*pPoint);
            pPoint->SetSolutionStepVariablesList(this->GetPoint(0).pGetVariablesList());
            pPoint->SetBufferSize(this->GetPoint(0).GetBufferSize());
            rPoints.push_back(pPoint);
        }
    }

    /**
     * Sampling the points on NURBS/Bezier geometry
     */
    void ExtractPoints(PointsArrayType& rPoints, const std::vector<int>& sampling_size) final
    {
        CoordinatesArrayType p_ref;
        CoordinatesArrayType p;

        // create and add nodes
        typedef typename PointType::Pointer PointPointerType;
        for(int i = 0; i <= sampling_size[0]; ++i)
        {
            p_ref[0] = ((double) i) / sampling_size[0];
            for(int j = 0; j <= sampling_size[1]; ++j)
            {
                p_ref[1] = ((double) j) / sampling_size[1];
                for(int k = 0; k <= sampling_size[2]; ++k)
                {
                    p_ref[2] = ((double) k) / sampling_size[2];
                    p = BaseType::GlobalCoordinates0(p, p_ref);
                    PointPointerType pPoint = PointPointerType(new PointType(0, p));
                    pPoint->SetSolutionStepVariablesList(this->GetPoint(0).pGetVariablesList());
                    pPoint->SetBufferSize(this->GetPoint(0).GetBufferSize());
                    rPoints.push_back(pPoint);
                }
            }
        }
    }

    /**
     * Extract the control values from NURBS/Bezier geometry
     */
    void ExtractControlValues(const Variable<double>& rVariable, std::vector<double>& rValues) const final
    {
        this->template ExtractControlValues_<double>(rVariable, rValues);
    }

    /**
     * Extract the control values from NURBS/Bezier geometry
     */
    void ExtractControlValues(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues) const final
    {
        this->template ExtractControlValues_<array_1d<double, 3> >(rVariable, rValues);
    }

    /**
     * Sampling the values on NURBS/Bezier geometry
     */
    void ExtractValues(const Variable<double>& rVariable, std::vector<double>& rValues, const std::vector<int>& sampling_size) const final
    {
        this->template ExtractValues_<double>(rVariable, rValues, sampling_size);
    }

    /**
     * Sampling the values on NURBS/Bezier geometry
     */
    void ExtractValues(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const std::vector<int>& sampling_size) const final
    {
        this->template ExtractValues_<array_1d<double, 3> >(rVariable, rValues, sampling_size);
    }

    /**
     * Returns whether given local point is inside the Geometry
     */
    bool IsInside( const CoordinatesArrayType& rPoint ) const final
    {
        const double tol = 1.0e-8;
        if ( (rPoint[0] > -tol) && (rPoint[0] < 1.0 + tol) )
            if ( (rPoint[1] > -tol) && (rPoint[1] < 1.0 + tol) )
                if ( (rPoint[2] > -tol) && (rPoint[2] < 1.0 + tol) )
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
    std::string Info() const override
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
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Geo3dBezier";
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
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        rOStream << std::endl;
        rOStream << "    Control Weights: " << mCtrlWeights << std::endl;
        rOStream << "    Order: " << mOrder1 << " " << mOrder2 << " " << mOrder3 << std::endl;
        rOStream << "    Number: " << mNumber1 << " " << mNumber2 << " " << mNumber3 << std::endl;
        rOStream << "    Extraction Operator: " << mExtractionOperator << std::endl;
    }

    void AssignGeometryData(
        const ValuesContainerType& Knots1, //not used
        const ValuesContainerType& Knots2, //not used
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        const int& Degree1,
        const int& Degree2,
        const int& Degree3,
        const int& NumberOfIntegrationMethod
    ) final
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
        if(mExtractionOperator.size2() != mNumber1 * mNumber2 * mNumber3)
        {
            KRATOS_WATCH(mExtractionOperator)
            KRATOS_WATCH(mOrder1)
            KRATOS_WATCH(mOrder2)
            KRATOS_WATCH(mOrder3)
            KRATOS_THROW_ERROR(std::logic_error, "The number of column of extraction operator must be equal to (p_u+1) * (p_v+1) * (p_w+1), error at", __FUNCTION__)
        }
        if(mCtrlWeights.size() != this->PointsNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The number of weights must be equal to number of nodes", __FUNCTION__)

        if(NumberOfIntegrationMethod > 0)
        {
            // find the existing integration rule or create new one if not existed
            BezierUtils::RegisterIntegrationRule<3, 3, 3>(NumberOfIntegrationMethod, Degree1, Degree2, Degree3);

            // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
            mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<3, 3, 3>(NumberOfIntegrationMethod, Degree1, Degree2, Degree3);
            #ifndef ENABLE_PRECOMPUTE
                #ifdef SD_APP_FORWARD_COMPATIBILITY
                BaseType::SetGeometryData(&(*mpBezierGeometryData));
                #else
                BaseType::mpGeometryData = &(*mpBezierGeometryData);
                #endif
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

                #ifdef SD_APP_FORWARD_COMPATIBILITY
                BaseType::SetGeometryData(&(*mpGeometryData));
                #else
                BaseType::mpGeometryData = &(*mpGeometryData);
                #endif
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

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    /**
     * Private Operations
     */

    /**
     * Calculate shape function values and local gradient at a particular point
     */
    void ShapeFunctionsValuesAndLocalGradients(
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    ) const final
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

    template<typename TDataType>
    void ExtractValues_(const Variable<TDataType>& rVariable, std::vector<TDataType>& rValues, const std::vector<int>& sampling_size) const
    {
        Vector shape_functions_values;

        CoordinatesArrayType p_ref;
        CoordinatesArrayType p;

        typedef typename PointType::Pointer PointPointerType;
        for(int i = 0; i <= sampling_size[0]; ++i)
        {
            p_ref[0] = ((double) i) / sampling_size[0];
            for(int j = 0; j <= sampling_size[1]; ++j)
            {
                p_ref[1] = ((double) j) / sampling_size[1];
                for(int k = 0; k <= sampling_size[2]; ++k)
                {
                    p_ref[2] = ((double) k) / sampling_size[2];

                    ShapeFunctionsValues(shape_functions_values, p_ref);

                    TDataType rResult = shape_functions_values( 0 ) * this->GetPoint( 0 ).GetSolutionStepValue(rVariable);
                    for ( IndexType i = 1 ; i < this->size() ; ++i )
                    {
                        rResult += shape_functions_values( i ) * this->GetPoint( i ).GetSolutionStepValue(rVariable);
                    }

                    rValues.push_back(rResult);
                }
            }
        }
    }

    template<typename TDataType>
    void ExtractControlValues_(const Variable<TDataType>& rVariable, std::vector<TDataType>& rValues) const
    {
        std::size_t number_of_points = this->PointsNumber();
        std::size_t number_of_local_points = mNumber1 * mNumber2 * mNumber3;
        if (rValues.size() != number_of_local_points)
            rValues.resize(number_of_local_points);

        // compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);

        // compute the Bezier control points
        for(std::size_t i = 0; i < number_of_local_points; ++i)
        {
            rValues[i] = TDataType(0.0);
            for(std::size_t j = 0; j < number_of_points; ++j)
                rValues[i] += mExtractionOperator(j, i) * this->GetPoint(j).GetSolutionStepValue(rVariable) * mCtrlWeights[j] / bezier_weights[i];
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

