/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Jan 29 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_2D_BEZIER_H_INCLUDED )
#define  KRATOS_GEO_2D_BEZIER_H_INCLUDED

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
//#include "integration/quadrature.h"
//#include "integration/line_gauss_legendre_integration_points.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
// #define DEBUG_LEVEL3
#define ENABLE_PROFILING

namespace Kratos
{

/**
 * A geometry representing Bezier decomposition of NURBS surface in 2D. In this implementation, surface XY is considerred. For a surface in 3D space, used Geo2dBezier3 instead
 */

template<class TPointType>
class Geo2dBezier : public IsogeometricGeometry<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * IsogeometricGeometry as base class.
     */
    typedef IsogeometricGeometry<TPointType> BaseType;

    /**
     * Pointer definition of Geo2dBezier
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo2dBezier );

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

    Geo2dBezier()
//    : BaseType( PointsArrayType(), &msGeometryData )
    : BaseType( PointsArrayType() ), mpBezierGeometryData(NULL)
    {}

    Geo2dBezier( const PointsArrayType& ThisPoints )
//    : BaseType( ThisPoints, &msGeometryData )
    : BaseType( ThisPoints ), mpBezierGeometryData(NULL)
    {}

//    Geo2dBezier( const PointsArrayType& ThisPoints, const GeometryData* pGeometryData )
//    : BaseType( ThisPoints, pGeometryData )
//    {}

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Geo2dBezier( Geo2dBezier const& rOther )
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
    template<class TOtherPointType> Geo2dBezier( Geo2dBezier<TOtherPointType> const& rOther )
    : BaseType( rOther ), mpBezierGeometryData(NULL)
    {}

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo2dBezier()
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
    Geo2dBezier& operator=( const Geo2dBezier& rOther )
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
    Geo2dBezier& operator=( Geo2dBezier<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        Geo2dBezier::Pointer pNewGeom = Geo2dBezier::Pointer( new Geo2dBezier( ThisPoints ) );
        ValuesContainerType DummyKnots;
        if (mpBezierGeometryData != NULL)
        {
            pNewGeom->AssignGeometryData(DummyKnots, DummyKnots, DummyKnots,
                mCtrlWeights, mExtractionOperator, mOrder1, mOrder2, 0,
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
////        p_clone( new Geo2dBezier< Point<3> >( NewPoints ) );

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
        return GeometryData::Kratos_Bezier2D;
    }

    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        IntegrationMethod ThisMethod
    ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        IndexType NumberOfIntegrationPoints = this->IntegrationPointsNumber(ThisMethod);

        shape_functions_values.resize(NumberOfIntegrationPoints, this->PointsNumber(), false);

        shape_functions_local_gradients.resize(NumberOfIntegrationPoints);
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType(this->PointsNumber(), 2));

        #ifdef DEBUG_LEVEL3
        KRATOS_WATCH(NumberOfIntegrationPoints)
        KRATOS_WATCH(mCtrlWeights)
        KRATOS_WATCH(mExtractionOperator)
        KRATOS_WATCH(mNumber1)
        KRATOS_WATCH(mNumber2)
        KRATOS_WATCH(this->PointsNumber())
        #endif

        const MatrixType& bezier_functions_values
//                = this->ShapeFunctionsValues(ThisMethod); // this is correct but dangerous
            = mpBezierGeometryData->ShapeFunctionsValues( ThisMethod );

        const ShapeFunctionsGradientsType& bezier_functions_local_gradients
//                = this->ShapeFunctionsLocalGradients(ThisMethod); // this is correct but dangerous
            = mpBezierGeometryData->ShapeFunctionsLocalGradients( ThisMethod );

        VectorType temp_bezier_values(bezier_functions_values.size2());
        VectorType bezier_weights(mExtractionOperator.size2());
        double denom, tmp1, tmp2;
        VectorType tmp_gradients1(this->PointsNumber());
        VectorType tmp_gradients2(this->PointsNumber());
        for(IndexType i = 0; i < NumberOfIntegrationPoints; ++i)
        {
            noalias(temp_bezier_values) = row(bezier_functions_values, i);

            //compute the Bezier weight
            noalias(bezier_weights) = prod(trans(mExtractionOperator), mCtrlWeights);
            denom = inner_prod(temp_bezier_values, bezier_weights);

            //compute the shape function values
            VectorType temp_values = prod(mExtractionOperator, temp_bezier_values);
            for(IndexType j = 0; j < this->PointsNumber(); ++j)
                shape_functions_values(i, j) = (temp_values(j) * mCtrlWeights(j)) / denom;

            //compute the shape function local gradients
//            shape_functions_local_gradients[i].resize(this->PointsNumber(), 2, false); // is not necessary when fill is used above
            tmp1 = inner_prod(row(bezier_functions_local_gradients[i], 0), bezier_weights);
            tmp2 = inner_prod(row(bezier_functions_local_gradients[i], 1), bezier_weights);

            noalias(tmp_gradients1) = prod(mExtractionOperator,
                    (1 / denom) * row(bezier_functions_local_gradients[i], 0) - (tmp1 / pow(denom, 2)) * temp_bezier_values );

            noalias(tmp_gradients2) = prod(mExtractionOperator,
                    (1 / denom) * row(bezier_functions_local_gradients[i], 1) - (tmp2 / pow(denom, 2)) * temp_bezier_values );

            for(IndexType j = 0; j < this->PointsNumber(); ++j)
            {
                shape_functions_local_gradients[i](j, 0) = tmp_gradients1(j) * mCtrlWeights(j);
                shape_functions_local_gradients[i](j, 1) = tmp_gradients2(j) * mCtrlWeights(j);
            }
        }
    }

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
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
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
            MatrixType jacobian = ZeroMatrix( 2, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
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
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
    {
        MatrixType shape_functions_values;
        ShapeFunctionsGradientsType shape_functions_local_gradients;

        //getting derivatives of shape functions
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
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
            MatrixType jacobian = ZeroMatrix( 2, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() - DeltaPosition(i, 0) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() - DeltaPosition(i, 1) ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    // REMARKS: Those InverseOfJacobian, DeterminantOfJacobian below is implemented in abstract Geometry class
//    /**
//     * TODO
//     */
//    virtual VectorType& DeterminantOfJacobian( VectorType& rResults,
//            IntegrationMethod ThisMethod ) const
//    {
//        JacobiansType J;
//        J = Jacobian( J, ThisMethod );
//
//        if ( rResults.size() != J.size() )
//        {
//            rResults.resize(J.size(), false);
//        }

//        for ( unsigned int pnt = 0; pnt < J.size(); ++pnt )
//        {
//            rResults[pnt] = MathUtils<double>::Det(J[pnt]);
//        }

//        return rResults;
//    }

//    /**
//     * TODO
//     */
//    virtual VectorType& DeterminantOfJacobian( VectorType& rResults,
//            IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
//    {
//        JacobiansType J;
//        J = Jacobian( J, ThisMethod, DeltaPosition );
//
//        if ( rResults.size() != J.size() )
//        {
//            rResults.resize(J.size(), false);
//        }

//        for ( unsigned int pnt = 0; pnt < J.size(); ++pnt )
//        {
//            rResults[pnt] = MathUtils<double>::Det(J[pnt]);
//        }

//        return rResults;
//    }

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
            MatrixType jacobian = ZeroMatrix( 2, 2 );

            //loop over all nodes
            for ( IndexType i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y0() ) * ( shape_functions_local_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Shape Function
     */

    /**
     * Compute shape function values at a particular reference point. This function is kept to keep the backward compatibility. The function ShapeFunctionsValuesAndLocalGradients is more general and direct to use.
     */
    virtual Vector& ShapeFunctionsValues( Vector& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(mNumber1);
        VectorType bezier_functions_values2(mNumber2);
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rCoordinates[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rCoordinates[1]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;
                bezier_functions_values(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_values2(j);
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
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rCoordinates[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rCoordinates[1]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;
                bezier_functions_values(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_values2(j);
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        VectorType bezier_functions_local_derivatives1(mNumber1 * mNumber2);
        VectorType bezier_functions_local_derivatives2(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;

                bezier_functions_local_derivatives1(index) =
                    bezier_functions_derivatives1(i) *
                    bezier_functions_values2(j);
                bezier_functions_local_derivatives2(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_derivatives2(j);
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local gradients
        rResults.resize(this->PointsNumber(), 2, false);
        double tmp1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double tmp2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
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
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults(i, 0) = tmp_gradients1(i) * mCtrlWeights(i);
            rResults(i, 1) = tmp_gradients2(i) * mCtrlWeights(i);
        }

        return rResults;
    }

    /**
     * Compute shape function second derivatives at a particular reference point. This function is kept to keep the backward compatibility.
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(mNumber1);
        VectorType bezier_functions_values2(mNumber2);
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        VectorType bezier_functions_second_derivatives1(mNumber1);
        VectorType bezier_functions_second_derivatives2(mNumber2);
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

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;
                bezier_functions_values(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_values2(j);
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        VectorType bezier_functions_local_derivatives1(mNumber1 * mNumber2);
        VectorType bezier_functions_local_derivatives2(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;

                bezier_functions_local_derivatives1(index) =
                    bezier_functions_derivatives1(i) *
                    bezier_functions_values2(j);
                bezier_functions_local_derivatives2(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_derivatives2(j);
            }
        }

        //compute trivariate Bezier shape functions second derivatives w.r.t local coordinates
        VectorType bezier_functions_local_second_derivatives11(mNumber1 * mNumber2);
        VectorType bezier_functions_local_second_derivatives12(mNumber1 * mNumber2);
        VectorType bezier_functions_local_second_derivatives22(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;

                bezier_functions_local_second_derivatives11(index) =
                    bezier_functions_second_derivatives1(i) *
                    bezier_functions_values2(j);
                bezier_functions_local_second_derivatives12(index) =
                    bezier_functions_derivatives1(i) *
                    bezier_functions_derivatives2(j);
                bezier_functions_local_second_derivatives22(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_second_derivatives2(j);
            }
        }

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local second gradients
        rResults.resize(this->PointsNumber(), false);
        double aux1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double aux2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
        double auxs11 = inner_prod(bezier_functions_local_second_derivatives11, bezier_weights);
        double auxs12 = inner_prod(bezier_functions_local_second_derivatives12, bezier_weights);
        double auxs22 = inner_prod(bezier_functions_local_second_derivatives22, bezier_weights);
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
        VectorType tmp_gradients22 =
            prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_second_derivatives22
                    - (aux2 / pow(denom, 2)) * bezier_functions_local_derivatives2 * 2
                    - (auxs22 / pow(denom, 2)) * bezier_functions_values
                    + 2.0 * pow(aux2, 2) / pow(denom, 3) * bezier_functions_values
            );
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults[i].resize(2, 2, false);
            rResults[i](0, 0) = tmp_gradients11(i) * mCtrlWeights(i);
            rResults[i](0, 1) = tmp_gradients12(i) * mCtrlWeights(i);
            rResults[i](1, 0) = rResults[i](0, 1);
            rResults[i](1, 1) = tmp_gradients22(i) * mCtrlWeights(i);
        }

        return rResults;
    }

    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        double tol = 1.0e-6;
        if ( (rResult[0] > -tol) && (rResult[0] < 1 + tol) )
            if ( (rResult[1] > -tol) && (rResult[1] < 1 + tol) )
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
        return "2 dimensional Bezier decomposition surface in 2D space";
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
        mCtrlWeights = Weights;
        mOrder1 = Degree1;
        mOrder2 = Degree2;
        mNumber1 = mOrder1 + 1;
        mNumber2 = mOrder2 + 1;
        // TODO we have to check here if the compressed_matrix copy is called or not. Otherwise, maybe the full matrix is populated.
        mExtractionOperator = ExtractionOperator;

        // size checking
        if(mExtractionOperator.size1() != this->PointsNumber())
            KRATOS_THROW_ERROR(std::logic_error, "The number of row of extraction operator must be equal to number of nodes, mExtractionOperator.size1() =", mExtractionOperator.size1())
        if(mExtractionOperator.size2() != (mOrder1 + 1) * (mOrder2 + 1))
            KRATOS_THROW_ERROR(std::logic_error, "The number of column of extraction operator must be equal to (p_u+1) * (p_v+1), mExtractionOperator.size2() =", mExtractionOperator.size2())

        if(NumberOfIntegrationMethod > 0)
        {
            // find the existing integration rule or create new one if not existed
            BezierUtils::RegisterIntegrationRule<2, 2, 2>(NumberOfIntegrationMethod, Degree1, Degree2);

            // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
            mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<2, 2, 2>(NumberOfIntegrationMethod, Degree1, Degree2);
            BaseType::mpGeometryData = &(*mpBezierGeometryData);
        }
    }

protected:

//    static const GeometryData msGeometryData;
    GeometryData::Pointer mpBezierGeometryData;

    MatrixType mExtractionOperator;
    // CompressedMatrixType mExtractionOperator;

    ValuesContainerType mCtrlWeights; //weight of control points

    int mOrder1; //order of the surface at parametric direction 1
    int mOrder2; //order of the surface at parametric direction 2

    int mNumber1; //number of bezier shape functions define the surface on parametric direction 1
    int mNumber2; //number of bezier shape functions define the surface on parametric direction 2

private:

    /**
     * Static Member Variables
     */

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
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
        VectorType bezier_functions_derivatives1(mNumber1);
        VectorType bezier_functions_derivatives2(mNumber2);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, mOrder1, rPoint[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, mOrder2, rPoint[1]);

        //compute trivariate Bezier shape functions values
        VectorType bezier_functions_values(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;
                bezier_functions_values(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_values2(j);
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        VectorType bezier_functions_local_derivatives1(mNumber1 * mNumber2);
        VectorType bezier_functions_local_derivatives2(mNumber1 * mNumber2);
        for(IndexType i = 0; i < mNumber1; ++i)
        {
            for(IndexType j = 0; j < mNumber2; ++j)
            {
                IndexType index = j + i * mNumber2;

                bezier_functions_local_derivatives1(index) =
                    bezier_functions_derivatives1(i) *
                    bezier_functions_values2(j);
                bezier_functions_local_derivatives2(index) =
                    bezier_functions_values1(i) *
                    bezier_functions_derivatives2(j);
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
            || shape_functions_local_gradients.size2() != 2)
            shape_functions_local_gradients.resize(this->PointsNumber(), 2, false);
        double tmp1 = inner_prod(bezier_functions_local_derivatives1, bezier_weights);
        double tmp2 = inner_prod(bezier_functions_local_derivatives2, bezier_weights);
        VectorType tmp_gradients1 = prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives1 - (tmp1 / pow(denom, 2)) * bezier_functions_values );
        VectorType tmp_gradients2 = prod(mExtractionOperator,
                    (1 / denom) * bezier_functions_local_derivatives2 - (tmp2 / pow(denom, 2)) * bezier_functions_values );
        for(IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            shape_functions_local_gradients(i, 0) = tmp_gradients1(i) * mCtrlWeights(i);
            shape_functions_local_gradients(i, 1) = tmp_gradients2(i) * mCtrlWeights(i);
        }
    }

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo2dBezier;

    /**
     * Un accessible methods
     */

};    // Class Geo2dBezier

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo2dBezier<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo2dBezier<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/**
 * static initialisation for geometry_data
 * TODO: to be deleted; the geometry_data is obtained by register with the BezierUtils
 */
// template<class TPointType>
// const GeometryData Geo2dBezier<TPointType>::msGeometryData
// (
//    2,
//    2,
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

