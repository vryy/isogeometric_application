/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2023 Aug 21 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_1D_BEZIER_2_H_INCLUDED )
#define  KRATOS_GEO_1D_BEZIER_2_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "custom_geometries/geo_1d_bezier.h"


namespace Kratos
{

/**
 * A geometry representing Bezier decomposition of a NURBS curve in 2D space
 */
template<class TPointType>
class Geo1dBezier2: public Geo1dBezier<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Pointer definition of Geo1dBezier2
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo1dBezier2 );

    /**
     * Geo1dBezier as base class.
     */
    typedef Geo1dBezier<TPointType> BaseType;

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

    Geo1dBezier2() : BaseType( PointsArrayType() )
    {}

    Geo1dBezier2(const PointsArrayType& ThisPoints)
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
    Geo1dBezier2( Geo1dBezier2 const& rOther )
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
    template<class TOtherPointType> Geo1dBezier2( Geo1dBezier2<TOtherPointType> const& rOther )
    : BaseType( rOther )
    {}

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo1dBezier2()
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
    Geo1dBezier2& operator=( const Geo1dBezier2& rOther )
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
    Geo1dBezier2& operator=( Geo1dBezier2<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Operations
     */

    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const final
    {
        Geo1dBezier2::Pointer pNewGeom = Geo1dBezier2::Pointer( new Geo1dBezier2( ThisPoints ) );
        ValuesContainerType DummyKnots;
        if (BaseType::mpBezierGeometryData != NULL)
        {
            pNewGeom->AssignGeometryData(DummyKnots, DummyKnots, DummyKnots,
                BaseType::mCtrlWeights, BaseType::mExtractionOperator, BaseType::mOrder, 0, 0,
                static_cast<int>(BaseType::mpBezierGeometryData->DefaultIntegrationMethod()) + 1);
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
////        p_clone( new Geo1dBezier2< Point<3> >( NewPoints ) );

////        p_clone->ClonePoints();

////        return p_clone;

//        KRATOS_THROW_ERROR(std::logic_error, "NURBS geometry does not support for Clone", *this)
//    }

    /**
     * Informations
     */

    GeometryData::KratosGeometryType GetGeometryType() const final
    {
        #ifdef SD_APP_FORWARD_COMPATIBILITY
        return static_cast<GeometryData::KratosGeometryType>(IsogeometricGeometryData::KratosIsogeometricGeometryType::Kratos_Bezier1D2);
        #else
        return GeometryData::KratosGeometryType::Kratos_Bezier1D2;
        #endif
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
            IntegrationMethod ThisMethod ) const final
    {
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            BaseType::CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //getting values of shape functions
        MatrixType shape_functions_values =
            BaseType::CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 2, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
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
            MatrixType& DeltaPosition ) const final
    {
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            BaseType::CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //getting values of shape functions
        MatrixType shape_functions_values =
            BaseType::CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            MatrixType jacobian = ZeroMatrix( 2, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
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
            IntegrationMethod ThisMethod ) const final
    {
        //setting up size of jacobian matrix
        rResult.resize( 2, 1 );
        //derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            BaseType::CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        MatrixType ShapeFunctionsGradientInIntegrationPoint =
            shape_functions_gradients( IntegrationPointIndex );

        //values of shape functions in integration points
        VectorType ShapeFunctionValuesInIntegrationPoint = ZeroVector( 3 );
        ShapeFunctionValuesInIntegrationPoint = row( BaseType::CalculateShapeFunctionsIntegrationPointsValues( ThisMethod ),
                IntegrationPointIndex );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;
        double j1 = 0.0;
        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            j1 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;
        rResult( 1, 0 ) = j1;

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
    MatrixType& Jacobian( MatrixType& rResult, const CoordinatesArrayType& rPoint ) const final
    {
        //setting up size of jacobian matrix
        rResult.resize( 2, 1 );

        //derivatives of shape functions
        MatrixType shape_functions_gradients;
        shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;
        double j1 = 0.0;

        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
            j1 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;
        rResult( 1, 0 ) = j1;

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
    VectorType& DeterminantOfJacobian( VectorType& rResult,
            IntegrationMethod ThisMethod ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return rResult;
    }

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
    double DeterminantOfJacobian( IndexType IntegrationPointIndex,
            IntegrationMethod ThisMethod ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return 0.0;
    }

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
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return 0.0;
    }

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
    JacobiansType& InverseOfJacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return rResult;
    }

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
    MatrixType& InverseOfJacobian( MatrixType& rResult,
            IndexType IntegrationPointIndex,
            IntegrationMethod ThisMethod ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return rResult;
    }

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
    MatrixType& InverseOfJacobian( MatrixType& rResult, const CoordinatesArrayType& rPoint ) const final
    {
        KRATOS_ERROR << "Jacobian is not square";
        return rResult;
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
    std::string Info() const final
    {
        return "1 dimensional Bezier extraction from NURBS curve in 2D space";
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
    ) final
    {
        BaseType::mCtrlWeights = Weights;
        BaseType::mOrder = Degree1;
        BaseType::mNumber = BaseType::mOrder + 1;
        BaseType::mExtractionOperator = ExtractionOperator;

        // size checking
        if(BaseType::mExtractionOperator.size1() != this->PointsNumber())
            KRATOS_ERROR << "The number of row of extraction operator must be equal to number of nodes";
        if(BaseType::mExtractionOperator.size2() != BaseType::mNumber)
            KRATOS_ERROR << "The number of column of extraction operator must be equal to (p_u+1)";
        if(BaseType::mCtrlWeights.size() != this->PointsNumber())
            KRATOS_ERROR << "The number of weights must be equal to number of nodes";

        // find the existing integration rule or create new one if not existed
        BezierUtils::RegisterIntegrationRule<1, 2, 1>(NumberOfIntegrationMethod, Degree1);

        // get the geometry_data according to integration rule. Note that this is a static geometry_data of a reference Bezier element, not the real Bezier element.
        BaseType::mpBezierGeometryData = BezierUtils::RetrieveIntegrationRule<1, 2, 1>(NumberOfIntegrationMethod, Degree1);
        #ifdef SD_APP_FORWARD_COMPATIBILITY
        BaseType::SetGeometryData(&(*BaseType::mpBezierGeometryData));
        #else
        BaseType::mpGeometryData = &(*BaseType::mpBezierGeometryData);
        #endif
    }

private:

    /**
     * Member Variables
     */

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

    template<class TOtherPointType> friend class Geo1dBezier2;

    /**
     * Un accessible methods
     */

};    // Class Geo1dBezier2

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo1dBezier2<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo1dBezier2<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}    // namespace Kratos.

#endif // KRATOS_GEO_1D_BEZIER_2_H_INCLUDED
