/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Jan 28 $
//
//
#if !defined(KRATOS_GEO_1D_BEZIER_H_INCLUDED )
#define  KRATOS_GEO_1D_BEZIER_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "integration/quadrature.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "custom_utilities/bezier_utils.h"

namespace Kratos
{

/**
 * A geometry representing Bezier decomposition of a NURBS curve
 */
template<class TPointType>
class Geo1dBezier: public IsogeometricGeometry<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Pointer definition of Geo1dBezier
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo1dBezier );

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
     * Type used for indexing in geometry class. IndexType used for indexing
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

    Geo1dBezier() : BaseType( PointsArrayType() )
    {}

    Geo1dBezier(const PointsArrayType& ThisPoints)
        : BaseType( ThisPoints )
    {}

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Geo1dBezier( Geo1dBezier const& rOther )
        : BaseType( rOther )
        , mpBezierGeometryData(rOther.mpBezierGeometryData)
        , mOrder(rOther.mOrder)
        , mNumber(rOther.mNumber)
        , mExtractionOperator(rOther.mExtractionOperator)
        , mCtrlWeights(rOther.mCtrlWeights)
    {
        GeometryType::SetGeometryData(mpBezierGeometryData.get());
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
    template<class TOtherPointType> Geo1dBezier( Geo1dBezier<TOtherPointType> const& rOther )
        : IsogeometricGeometry<TOtherPointType>( rOther )
        , mpBezierGeometryData(rOther.mpBezierGeometryData)
        , mOrder(rOther.mOrder)
        , mNumber(rOther.mNumber)
        , mExtractionOperator(rOther.mExtractionOperator)
        , mCtrlWeights(rOther.mCtrlWeights)
    {
        Geometry<TOtherPointType>::SetGeometryData(mpBezierGeometryData.get());
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Geo1dBezier() override
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
    Geo1dBezier& operator=( const Geo1dBezier& rOther )
    {
        BaseType::operator=( rOther );
        this->mpBezierGeometryData = rOther.mpBezierGeometryData;
        GeometryType::mpGeometryData = this->mpBezierGeometryData.get();
        this->mOrder = rOther.mOrder;
        this->mNumber = rOther.mNumber;
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
    Geo1dBezier& operator=( Geo1dBezier<TOtherPointType> const & rOther )
    {
        IsogeometricGeometry<TOtherPointType>::operator=( rOther );
        this->mpBezierGeometryData = rOther.mpBezierGeometryData;
        Geometry<TOtherPointType>::mpGeometryData = this->mpBezierGeometryData.get();
        this->mOrder = rOther.mOrder;
        this->mNumber = rOther.mNumber;
        this->mExtractionOperator = rOther.mExtractionOperator;
        this->mCtrlWeights = rOther.mCtrlWeights;
        return *this;
    }

    /**
     * Operations
     */

    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        Geo1dBezier::Pointer pNewGeom = Geo1dBezier::Pointer( new Geo1dBezier( ThisPoints ) );
        ValuesContainerType DummyKnots;
        if (mpBezierGeometryData != nullptr)
        {
            pNewGeom->AssignGeometryData(DummyKnots, DummyKnots, DummyKnots,
                                         mCtrlWeights, mExtractionOperator, mOrder, 0, 0,
                                         static_cast<int>(mpBezierGeometryData->DefaultIntegrationMethod()) + 1);
        }
        return pNewGeom;
    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
////        Geometry< Point<3> >::PointsArrayType NewPoints;
////        //making a copy of the nodes TO POINTS (not Nodes!!!)

////        for ( IndexType i = 0; i < this->Points().size(); ++i )
////            #if defined(KRATOS_SD_REF_NUMBER_2)
////            NewPoints.push_back( this->Points()[i] );
////            #elif defined(KRATOS_SD_REF_NUMBER_3)
////            NewPoints.push_back(boost::make_shared< Point<3> >(( *this )[i]));
////            #endif

////        //creating a geometry with the new points
////        boost::shared_ptr< Geometry< Point<3> > >
////        p_clone( new Geo1dBezier< Point<3> >( NewPoints ) );

////        p_clone->ClonePoints();

////        return p_clone;

//        KRATOS_ERROR << "NURBS geometry does not support for Clone";
//    }

    /**
     * Informations
     */

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        return static_cast<GeometryData::KratosGeometryType>(IsogeometricGeometryData::KratosIsogeometricGeometryType::Kratos_Bezier1D);
#else
        return GeometryData::KratosGeometryType::Kratos_Bezier1D;
#endif
    }

    // /**
    //  * This method calculates and returns Length or charactereistic
    //  * length of this geometry depending on it's dimension. For one
    //  * dimensional geometry for example Line it returns length of it
    //  * and for the other geometries it gives Characteristic length
    //  * otherwise.
    //  *
    //  * @return double value contains length or Characteristic
    //  * length
    //  * @see Area()
    //  * @see Volume()
    //  * @see DomainSize()
    //  */
    // double Length() const final
    // {
    //     //TODO: reimplement to account for nurbs curve
    //     KRATOS_ERROR <<  __FUNCTION__ << " is not yet implemented";
    //     return 0.0;
    // }

    // /**
    //  * This method calculates and returns area or surface area of
    //  * this geometry depending to it's dimension. For one dimensional
    //  * geometry it returns zero, for two dimensional it gives area
    //  * and for three dimensional geometries it gives surface area.
    //  *
    //  * @return double value contains area or surface area.
    //  * @see Length()
    //  * @see Volume()
    //  * @see DomainSize()
    //  */
    // double Area() const final
    // {
    //     // Area is not relevant for 1d geometry
    //     return 0.0;
    // }

    // double Volume() const final
    // {
    //     // Volume is not relevant for 1d geometry
    //     return 0.0;
    // }

    /**
     * This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    // double DomainSize() const final
    // {
    //     return Length();
    // }

    /**
     * Returns whether given local point is inside the Geometry
     */
    bool IsInside( const CoordinatesArrayType& rPoint ) const final
    {
        const double tol = 1.0e-8;
        if ( (rPoint[0] > -tol) && (rPoint[0] < 1.0 + tol) )
        {
            return true;
        }

        return false;
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
    JacobiansType& Jacobian( JacobiansType& rResult,
                             IntegrationMethod ThisMethod ) const override
    {
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //getting values of shape functions
        MatrixType shape_functions_values =
            CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 1, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
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
     * @param DeltaPosition MatrixType with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian( JacobiansType& rResult,
                             IntegrationMethod ThisMethod,
                             const MatrixType& DeltaPosition ) const override
    {
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //getting values of shape functions
        MatrixType shape_functions_values =
            CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 1, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() + DeltaPosition(i, 0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Jacobian in specific integration point of given integration
     * method. This method calculate jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return MatrixType(double) Jacobian matrix \f$ J_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    MatrixType& Jacobian( MatrixType& rResult,
                          IndexType IntegrationPointIndex,
                          IntegrationMethod ThisMethod ) const override
    {
        //setting up size of jacobian matrix
        rResult.resize( 1, 1 );
        //derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        MatrixType ShapeFunctionsGradientInIntegrationPoint =
            shape_functions_gradients( IntegrationPointIndex );

        //values of shape functions in integration points
        vector<double> ShapeFunctionValuesInIntegrationPoint = ZeroVector( 3 );
        ShapeFunctionValuesInIntegrationPoint = row( CalculateShapeFunctionsIntegrationPointsValues( ThisMethod ),
                                                IntegrationPointIndex );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;
        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;

        return rResult;
    }

    /**
     * Jacobian in given point. This method calculate jacobian
     * matrix in given point.
     *
     * @param rPoint point which jacobians has to
     * be calculated in it.
     *
     * @return MatrixType of double which is jacobian matrix \f$ J \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    MatrixType& Jacobian( MatrixType& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        //setting up size of jacobian matrix
        rResult.resize( 1, 1 );

        //derivatives of shape functions
        MatrixType shape_functions_gradients;
        shape_functions_gradients = ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;

        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;

        return rResult;
    }

    /**
     * Determinant of jacobians for given integration method.
     * This method calculate determinant of jacobian in all
     * integrations points of given integration method.
     *
     * @return VectorType of double which is vector of determinants of
     * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    // VectorType& DeterminantOfJacobian( VectorType& rResult,
    //         IntegrationMethod ThisMethod ) const final
    // {
    //     // TODO
    //     return rResult;
    // }

    /**
     * Determinant of jacobian in specific integration point of
     * given integration method. This method calculate determinant
     * of jacobian in given integration point of given integration
     * method.
     *
     * @param IntegrationPointIndex index of integration point which
     * jacobians has to be calculated in it.
     *
     * @param IntegrationPointIndex index of integration point
     * which determinant of jacobians has to be calculated in it.
     *
     * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    // double DeterminantOfJacobian( IndexType IntegrationPointIndex,
    //         IntegrationMethod ThisMethod ) const final
    // {
    //     // TODO
    //     return 0.0;
    // }

    /**
     * Determinant of jacobian in given point.
     * This method calculate determinant of jacobian
     * matrix in given point.
     *
     * @param rPoint point which determinant of jacobians has to
     * be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
     * point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: PointType needed for proper functionality
     * KLUDGE: works only with explicitly generated MatrixType object
     */
    // double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const final
    // {
    //     // TODO
    //     return 0.0;
    // }

    /**
     * Inverse of jacobians for given integration method.
     * This method calculate inverse of jacobians matrices in all
     * integrations points of given integration method.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     * @return Inverse of jacobian
     * matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     *
     * KLUDGE: works only with explicitly generated MatrixType object
     */
    // JacobiansType& InverseOfJacobian( JacobiansType& rResult,
    //         IntegrationMethod ThisMethod ) const final
    // {
    //     // TODO
    //     return rResult;
    // }

    /**
     * Inverse of jacobian in specific integration point of given integration
     * method.
     * This method calculate Inverse of jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which
     * inverse of jacobians has to be calculated in it.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     *
     * @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     *
     * KLUDGE: works only with explicitly generated MatrixType object
     */
    // MatrixType& InverseOfJacobian( MatrixType& rResult,
    //         IndexType IntegrationPointIndex,
    //         IntegrationMethod ThisMethod ) const final
    // {
    //     // TODO
    //     return rResult;
    // }

    /**
     * Inverse of jacobian in given point.
     * This method calculate inverse of jacobian matrix in given point.
     *
     * @param rPoint point which inverse of jacobians has to
     * be calculated in it.
     * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: works only with explicitly generated MatrixType object
     */
    // MatrixType& InverseOfJacobian( MatrixType& rResult, const CoordinatesArrayType& rPoint ) const final
    // {
    //     // TODO
    //     return rResult;
    // }

    /**
     * Shape Function
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     */
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                               const CoordinatesArrayType& rPoint ) const final
    {
        //compute all Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        BezierUtils::bernstein(bezier_functions_values, mOrder, rPoint[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        VectorType shape_functions_values(this->PointsNumber());
        noalias( shape_functions_values ) = prod(mExtractionOperator, bezier_functions_values);

        return shape_functions_values(ShapeFunctionIndex) *
               mCtrlWeights(ShapeFunctionIndex) / denom;
    }

    /**
     * Calculates the value of all shape functions at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     */
    VectorType& ShapeFunctionsValues( VectorType& rResults,
                                      const CoordinatesArrayType& rPoint ) const final
    {
        //compute all Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        BezierUtils::bernstein(bezier_functions_values, mOrder, rPoint[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        rResults.resize(this->PointsNumber(), false);
        noalias( rResults ) = prod(mExtractionOperator, bezier_functions_values);

        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults(i) *= (mCtrlWeights(i) / denom);
        }

        return rResults;
    }

    /**
     * Calculates the local gradients at a given point
     */
    MatrixType& ShapeFunctionsLocalGradients( MatrixType& rResult,
            const CoordinatesArrayType& rPoint ) const final
    {
        //compute all Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        VectorType bezier_functions_derivatives(mNumber);
        BezierUtils::bernstein(bezier_functions_values, bezier_functions_derivatives, mOrder, rPoint[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        VectorType shape_functions_values(this->PointsNumber());
        noalias(shape_functions_values) = prod(mExtractionOperator, bezier_functions_values);
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            shape_functions_values(i) *= (mCtrlWeights(i) / denom);
        }

        //compute the shape function local gradients
        rResult.resize(this->PointsNumber(), 1, false);
        double tmp = inner_prod(bezier_functions_derivatives, bezier_weights);
        VectorType tmp_gradients =
            prod(
                mExtractionOperator,
                (1 / denom) * bezier_functions_derivatives -
                (tmp / pow(denom, 2)) * bezier_functions_values
            );
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResult(i, 0) = tmp_gradients(i) * mCtrlWeights(i);
        }

        return rResult;
    }

    void ShapeFunctionsValuesAndLocalGradients( VectorType& shape_functions_values,
            MatrixType& shape_functions_local_gradients, const CoordinatesArrayType& rPoint ) const
    {
        //compute all Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        VectorType bezier_functions_derivatives(mNumber);
        BezierUtils::bernstein(bezier_functions_values, bezier_functions_derivatives, mOrder, rPoint[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function values
        shape_functions_values.resize(this->PointsNumber(), false);
        noalias( shape_functions_values ) = prod(mExtractionOperator, bezier_functions_values);
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            shape_functions_values(i) *= (mCtrlWeights(i) / denom);
        }

        //compute the shape function local gradients
        shape_functions_local_gradients.resize(this->PointsNumber(), 1, false);
        double tmp = inner_prod(bezier_functions_derivatives, bezier_weights);
        VectorType tmp_gradients =
            prod(
                mExtractionOperator,
                (1 / denom) * bezier_functions_derivatives -
                (tmp / pow(denom, 2)) * bezier_functions_values
            );
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            shape_functions_local_gradients(i, 0) = tmp_gradients(i) * mCtrlWeights(i);
        }
    }

    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to the global
     * coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     * @return the gradients of all shape functions with regard to the global coordinates
     *
     * KLUDGE: method call only works with explicit JacobiansType rather than creating
     * JacobiansType within argument list
     */
    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, IntegrationMethod ThisMethod ) const final
    {
        KRATOS_ERROR << __FUNCTION__ << " is not yet implemented.";

        return rResult;
    }

    /**
     * Compute shape function second derivatives at a particular reference point.
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResults, const CoordinatesArrayType& rCoordinates ) const final
    {
#ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
#endif

        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        VectorType bezier_functions_derivatives(mNumber);
        VectorType bezier_functions_second_derivatives(mNumber);
        BezierUtils::bernstein(bezier_functions_values,
                               bezier_functions_derivatives,
                               bezier_functions_second_derivatives,
                               mOrder,
                               rCoordinates[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local second gradients
        rResults.resize(this->PointsNumber(), false);
        double aux = inner_prod(bezier_functions_derivatives, bezier_weights);
        double aux2 = inner_prod(bezier_functions_second_derivatives, bezier_weights);
        VectorType tmp_gradients =
            prod(mExtractionOperator,
                 (1 / denom) * bezier_functions_second_derivatives
                 - 2.0 * (aux / pow(denom, 2)) * bezier_functions_derivatives
                 - (aux2 / pow(denom, 2)) * bezier_functions_values
                 + 2.0 * pow(aux, 2) / pow(denom, 3) * bezier_functions_values
                );
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults[i].resize(1, 1, false);
            rResults[i](0, 0) = tmp_gradients(i) * mCtrlWeights(i);
        }

        return rResults;
    }

    /**
     * Compute shape function third derivatives at a particular reference point.
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResults, const CoordinatesArrayType& rCoordinates ) const final
    {
        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values(mNumber);
        VectorType bezier_functions_derivatives(mNumber);
        VectorType bezier_functions_second_derivatives(mNumber);
        VectorType bezier_functions_third_derivatives(mNumber);
        BezierUtils::bernstein(bezier_functions_values,
                               bezier_functions_derivatives,
                               bezier_functions_second_derivatives,
                               bezier_functions_third_derivatives,
                               mOrder,
                               rCoordinates[0]);

        //compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);
        double denom = inner_prod(bezier_functions_values, bezier_weights);

        //compute the shape function local third gradients
        rResults.resize(this->PointsNumber(), false);
        double aux = inner_prod(bezier_functions_derivatives, bezier_weights);
        double aux2 = inner_prod(bezier_functions_second_derivatives, bezier_weights);
        double aux3 = inner_prod(bezier_functions_third_derivatives, bezier_weights);
        VectorType tmp_gradients =
            prod(mExtractionOperator,
                 (1 / denom) * bezier_functions_third_derivatives
                 - 3.0 * (aux / pow(denom, 2)) * bezier_functions_second_derivatives
                 + ( - 3.0 * (aux2 / pow(denom, 2))
                     + 6.0 * (pow(aux, 2) / pow(denom, 3)) ) * bezier_functions_derivatives
                 + ( -(aux3 / pow(denom, 2))
                     + 6.0 * (aux2 * aux / pow(denom, 3))
                     - 6.0 * (pow(aux, 3) / pow(denom, 4)) ) * bezier_functions_values
                );
        for (IndexType i = 0; i < this->PointsNumber(); ++i)
        {
            rResults[i].resize(1, false);
            for (IndexType j = 0; j < 1; ++j)
            {
                rResults[i][j].resize(1, 1, false);
                rResults[i][j](0, 0) = tmp_gradients(i) * mCtrlWeights(i);
            }
        }

        return rResults;
    }

    /**
     * Compute the Bezier control points
     */
    void ExtractControlPoints(PointsArrayType& rPoints) final
    {
        SizeType number_of_points = this->PointsNumber();
        SizeType number_of_local_points = mNumber;
        rPoints.clear();
        rPoints.reserve(number_of_local_points);

        // compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);

        // compute the Bezier control points
        typedef typename PointType::Pointer PointPointerType;
        for (IndexType i = 0; i < number_of_local_points; ++i)
        {
            PointPointerType pPoint = PointPointerType(new PointType(0, 0.0, 0.0, 0.0));
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                noalias(*pPoint) += mExtractionOperator(j, i) * this->GetPoint(j).GetInitialPosition() * mCtrlWeights[j] / bezier_weights[i];
            }
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
        p_ref[1] = 0.0;
        p_ref[2] = 0.0;
        typedef typename PointType::Pointer PointPointerType;
        for (int i = 0; i <= sampling_size[0]; ++i)
        {
            p_ref[0] = ((double) i) / sampling_size[0];
            p = BaseType::GlobalCoordinates0(p, p_ref);
            PointPointerType pPoint = PointPointerType(new PointType(0, p));
            pPoint->SetSolutionStepVariablesList(this->GetPoint(0).pGetVariablesList());
            pPoint->SetBufferSize(this->GetPoint(0).GetBufferSize());
            rPoints.push_back(pPoint);
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
        return "1 dimensional Bezier extraction from NURBS curve in 1D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const final
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
    void PrintData( std::ostream& rOStream ) const final
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        MatrixType jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     *
     */
    void AssignGeometryData(
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2,
        const ValuesContainerType& Knots3,
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        int Degree1,
        int Degree2,
        int Degree3,
        int NumberOfIntegrationMethod
    ) override
    {
        mCtrlWeights = Weights;
        mOrder = Degree1;
        mNumber = mOrder + 1;
        mExtractionOperator = ExtractionOperator;

        // size checking
        if (mExtractionOperator.size1() != this->PointsNumber())
        {
            KRATOS_ERROR << "The number of row of extraction operator must be equal to number of nodes";
        }
        if (mExtractionOperator.size2() != mNumber)
        {
            KRATOS_ERROR << "The number of column of extraction operator must be equal to (p_u+1)";
        }
        if (mCtrlWeights.size() != this->PointsNumber())
        {
            KRATOS_ERROR << "The number of weights must be equal to number of nodes";
        }

        // find the existing integration rule or create new one if not existed
        BezierUtils::RegisterIntegrationRule<1, 1, 1>(NumberOfIntegrationMethod, Degree1);

        // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
        mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<1, 1, 1>(NumberOfIntegrationMethod, Degree1);
        BaseType::SetGeometryData(mpBezierGeometryData.get());
    }

protected:

    /**
     * Member Variables
     */

    GeometryData::Pointer mpBezierGeometryData;

    MatrixType mExtractionOperator;

    ValuesContainerType mCtrlWeights; // weight of control points

    int mOrder; // order of the curve

    int mNumber; // number of Bezier shape functions

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // method to build GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     * KLUDGE: values are hard-coded!
     */
    MatrixType CalculateShapeFunctionsIntegrationPointsValues(IntegrationMethod ThisMethod ) const
    {
        const IntegrationPointsArrayType& integration_points = BaseType::IntegrationPoints(ThisMethod);

        //number of integration points
        const int integration_points_number = integration_points.size();

        //setting up return matrix
        MatrixType shape_functions_values( integration_points_number, this->PointsNumber() );

        //loop over all integration points
        for (IndexType pnt = 0; pnt < integration_points_number; ++pnt)
        {
            for (IndexType node = 0; node < this->PointsNumber(); ++node)
            {
                shape_functions_values( pnt, node ) = ShapeFunctionValue(node, integration_points[pnt]);
            }
        }

        return shape_functions_values;
    }

    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
     * KLUGDE: gradients are hard-coded!
     */
    ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(IntegrationMethod ThisMethod ) const
    {
        const IntegrationPointsArrayType& integration_points = BaseType::IntegrationPoints(ThisMethod);
        ShapeFunctionsGradientsType DN_De( integration_points.size() );
        std::fill( DN_De.begin(), DN_De.end(), MatrixType( this->PointsNumber(), 1 ) );

        for (IndexType it_gp = 0; it_gp < integration_points.size(); ++it_gp )
        {
            ShapeFunctionsLocalGradients(DN_De[it_gp], integration_points[it_gp]);
        }

        return DN_De;
    }

    void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        const IntegrationPointsArrayType& integration_points
    ) const final
    {
        if (shape_functions_local_gradients.size() != integration_points.size())
        {
            shape_functions_local_gradients.resize(integration_points.size());
        }
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType());
        if (shape_functions_values.size1() != integration_points.size() || shape_functions_values.size2() != this->PointsNumber())
        {
            shape_functions_values.resize(integration_points.size(), this->PointsNumber(), false);
        }

        for (IndexType it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            VectorType temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                temp_values,
                shape_functions_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for (IndexType node = 0; node < this->PointsNumber(); ++node)
            {
                shape_functions_values( it_gp, node ) = temp_values(node);
            }
        }
    }

    void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        IntegrationMethod ThisMethod
    ) const final
    {
        const IntegrationPointsArrayType& integration_points = BaseType::IntegrationPoints(ThisMethod);
        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(shape_functions_values, shape_functions_local_gradients, integration_points);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // end of method to build to GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<typename TDataType>
    void ExtractValues_(const Variable<TDataType>& rVariable, std::vector<TDataType>& rValues, const std::vector<int>& sampling_size) const
    {
        Vector shape_functions_values;

        CoordinatesArrayType p_ref;
        CoordinatesArrayType p;

        p_ref[1] = 0.0;
        p_ref[2] = 0.0;
        typedef typename PointType::Pointer PointPointerType;
        for (int i = 0; i <= sampling_size[0]; ++i)
        {
            p_ref[0] = ((double) i) / sampling_size[0];

            ShapeFunctionsValues(shape_functions_values, p_ref);

            TDataType rResult = shape_functions_values( 0 ) * this->GetPoint( 0 ).GetSolutionStepValue(rVariable);
            for ( IndexType i = 1 ; i < this->size() ; ++i )
            {
                rResult += shape_functions_values( i ) * this->GetPoint( i ).GetSolutionStepValue(rVariable);
            }

            rValues.push_back(rResult);
        }
    }

    template<typename TDataType>
    void ExtractControlValues_(const Variable<TDataType>& rVariable, std::vector<TDataType>& rValues) const
    {
        SizeType number_of_points = this->PointsNumber();
        SizeType number_of_local_points = mNumber;
        if (rValues.size() != number_of_local_points)
        {
            rValues.resize(number_of_local_points);
        }

        // compute the Bezier weight
        VectorType bezier_weights = prod(trans(mExtractionOperator), mCtrlWeights);

        // compute the Bezier control points
        for (IndexType i = 0; i < number_of_local_points; ++i)
        {
            rValues[i] = TDataType(0.0);
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                rValues[i] += mExtractionOperator(j, i) * this->GetPoint(j).GetSolutionStepValue(rVariable) * mCtrlWeights[j] / bezier_weights[i];
            }
        }
    }

private:

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
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo1dBezier;

    /**
     * Un accessible methods
     */

};    // Class Geo1dBezier

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
    std::istream& rIStream, Geo1dBezier<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
    std::ostream& rOStream, const Geo1dBezier<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}    // namespace Kratos.

#endif
