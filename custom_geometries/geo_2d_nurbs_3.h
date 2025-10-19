//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2016 Jun 1 $
//
//
#if !defined(KRATOS_GEO_2D_NURBS_3_H_INCLUDED )
#define  KRATOS_GEO_2D_NURBS_3_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/openmp_utils.h"
#include "geometries/geometry.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "custom_geometries/geo_2d_nurbs.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_LEVEL3
#define ENABLE_PROFILING

namespace Kratos
{

/**
 * A geometry representing NURBS surface in 3D
 */
template<class TPointType>
class Geo2dNURBS3: public Geo2dNURBS<TPointType>
{
public:

    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef IsogeometricGeometry<TPointType> SuperType;
    typedef Geo2dNURBS<TPointType> BaseType;

    /**
     * Pointer definition of Geo2dNURBS3
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo2dNURBS3 );

    /**
     * Integration methods implemented in geometry.
     */
    typedef typename BaseType::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries. Used for
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
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
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

    Geo2dNURBS3(): BaseType( PointsArrayType() )
    {}

    Geo2dNURBS3( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints )
    {
//        KRATOS_WATCH("at Geo2dNURBS3 constructor")
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
    Geo2dNURBS3( Geo2dNURBS3 const& rOther )
        : BaseType( rOther )
    {
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
    template<class TOtherPointType> Geo2dNURBS3( Geo2dNURBS3<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Geo2dNURBS3() override
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
    Geo2dNURBS3& operator=( const Geo2dNURBS3& rOther )
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
    Geo2dNURBS3& operator=( Geo2dNURBS3<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::BaseType::BaseType::Pointer( new Geo2dNURBS3( ThisPoints ) );
    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
//        Geometry< Point<3> >::PointsArrayType NewPoints;
//        //making a copy of the nodes TO POINTS (not Nodes!!!)

//        for ( IndexType i = 0; i < this->Points().size(); ++i )
//        NewPoints.push_back( this->Points()[i] );

//        //creating a geometry with the new points
//        boost::shared_ptr< Geometry< Point<3> > >
//        p_clone( new Geo2dNURBS3< Point<3> >( NewPoints ) );

//        p_clone->ClonePoints();

//        return p_clone;
//    }

    /**
     * TODO: TO BE CHECKED!!!!!!!!!!!
     */
    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        //TODO
        return rResult;
    }

    /**
     * Informations
     */

    /**
     * Jacobian
     */

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Jacobians for given  method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod ) const
    {
//        //getting derivatives of shape functions
//        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
//        //getting values of shape functions
//        const Matrix& shape_functions_values = BaseType::ShapeFunctionsValues( ThisMethod );

//        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
//        {
//            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
//            rResult.swap( temp );
//        }

//        //loop over all integration points
//        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
//        {
//            //defining single jacobian matrix
//            Matrix jacobian = ZeroMatrix( 2, 2 );
//            //loop over all nodes

//            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
//            {
//                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
//                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
//                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
//                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
//            }

//            rResult[pnt] = jacobian;
//        } //end of loop over all integration points

//        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Jacobians for given  method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @param DeltaPosition Matrix with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian( JacobiansType& rResult,
                             IntegrationMethod ThisMethod,
                             const Matrix& DeltaPosition ) const override
    {
//        //getting derivatives of shape functions
//        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
//        //getting values of shape functions
//        const Matrix& shape_functions_values = BaseType::ShapeFunctionsValues( ThisMethod );

//        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
//        {
//            // KLUDGE: While there is a bug in ublas
//            // vector resize, I have to put this beside resizing!!
//            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
//            rResult.swap( temp );
//        }

//        //loop over all integration points
//        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
//        {
//            //defining single jacobian matrix
//            Matrix jacobian = ZeroMatrix( 2, 2 );
//            //loop over all nodes
//            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
//            {
//                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
//                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
//                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
//                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
//            }

//            rResult[pnt] = jacobian;
//        } //end of loop over all integration points

//        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Jacobian in specific integration point of given integration
     * method. This method calculate jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return Matrix(double) Jacobian matrix \f$ J_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual Matrix& Jacobian( Matrix& rResult,
                              IndexType IntegrationPointIndex,
                              IntegrationMethod ThisMethod ) const
    {
//        //setting up size of jacobian matrix
//        rResult.resize( 2, 2 );
//        //derivatives of shape functions
//        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
//        const Matrix& ShapeFunctionsGradientInIntegrationPoint = shape_functions_gradients( IntegrationPointIndex );

//        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
//        //loop over all nodes
//        double j00 = 0.0;
//        double j01 = 0.0;
//        double j10 = 0.0;
//        double j11 = 0.0;
//        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
//        {
//            j00 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
//            j01 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
//            j10 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
//            j11 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
//        }

//        rResult( 0, 0 ) = j00;
//        rResult( 0, 1 ) = j01;
//        rResult( 1, 0 ) = j10;
//        rResult( 1, 1 ) = j11;
//
//        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Jacobian in given point. This method calculate jacobian
     * matrix in given point.
     *
     * @param rPoint point which jacobians has to
     * be calculated in it.
     *
     * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
//        //setting up size of jacobian matrix
//        rResult.resize( 2, 2 );
//        //derivatives of shape functions
//        Matrix shape_functions_gradients;
//        shape_functions_gradients = this->ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );
//        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
//        //loop over all nodes

//        double j00 = 0.0;
//        double j01 = 0.0;
//        double j10 = 0.0;
//        double j11 = 0.0;
//        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
//        {
//            j00 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
//            j01 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 1 ) );
//            j10 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
//            j11 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 1 ) );
//        }
//
//        rResult( 0, 0 ) = j00;
//        rResult( 0, 1 ) = j01;
//        rResult( 1, 0 ) = j10;
//        rResult( 1, 1 ) = j11;

//        return rResult;
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
        return "2 dimensional NURBS surface in 3D space";
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
//        Matrix jacobian;
//        Jacobian( jacobian, PointType() );
//        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

protected:

    /**
     * there are no protected class members
     */

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

    template<class TOtherPointType> friend class Geo2dNURBS3;

    /**
     * Un accessible methods
     */

};    // Class Geo2dNURBS3

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
    std::istream& rIStream, Geo2dNURBS3<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
    std::ostream& rOStream, const Geo2dNURBS3<TPointType>& rThis)
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

#endif
