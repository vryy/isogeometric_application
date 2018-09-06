/*
 ==============================================================================
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 CIMNE (International Center for Numerical Methods in Engineering),
 Gran Capita' s/n, 08034 Barcelona, Spain

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 ==============================================================================
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013 Sep 14 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_2D_NURBS_H_INCLUDED )
#define  KRATOS_GEO_2D_NURBS_H_INCLUDED

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
#include "integration/quadrature.h"
#include "custom_utilities/bspline_utils.h"
#include "integration/quadrature.h"
#include "integration/line_gauss_legendre_integration_points.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_LEVEL3
#define ENABLE_PROFILING

namespace Kratos
{

/**
 * A geometry representing NURBS surface
 */
template<class TPointType>
class Geo2dNURBS: public IsogeometricGeometry<TPointType>
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
     * Pointer definition of Geo2dNURBS
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo2dNURBS );

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

    Geo2dNURBS(): BaseType( PointsArrayType() )
    {}

    Geo2dNURBS( const PointsArrayType& ThisPoints )
    : BaseType( ThisPoints )
    {
//        KRATOS_WATCH("at Geo2dNURBS constructor")
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
    Geo2dNURBS( Geo2dNURBS const& rOther )
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
    template<class TOtherPointType> Geo2dNURBS( Geo2dNURBS<TOtherPointType> const& rOther )
    : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo2dNURBS()
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
    Geo2dNURBS& operator=( const Geo2dNURBS& rOther )
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
    Geo2dNURBS& operator=( Geo2dNURBS<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::BaseType::Pointer( new Geo2dNURBS( ThisPoints ) );
    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
//        Geometry< Point<3> >::PointsArrayType NewPoints;
//        //making a copy of the nodes TO POINTS (not Nodes!!!)

//        for ( IndexType i = 0; i < this->Points().size(); ++i )
//        NewPoints.push_back( this->Points()[i] );

//        //creating a geometry with the new points
//        boost::shared_ptr< Geometry< Point<3> > >
//        p_clone( new Geo2dNURBS< Point<3> >( NewPoints ) );

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
     * This method calculates and returns Length or charactereistic
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
     *
     * :TODO: might need to be changed to be useful!
     */
    virtual double Length() const
    {
        //this can be considered as characteristic length of this element
        double length = 0.000;
        length = sqrt( fabs( Area() ) );
        return length;
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
     *
     * :TODO: might need to be changed to be useful!
     */
    virtual double Area() const
    {
        Vector temp;
        DeterminantOfJacobian( temp, mpGeometryData->DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( mpGeometryData->DefaultIntegrationMethod() );
        double Area = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); ++i )
        {
            Area += temp[i] * integration_points[i].Weight();
        }
        return Area;
    }

    virtual double Volume() const
    {
        // Volume is not relevant for 2d geometry
        return 0.0;
    }

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
     *
     * :TODO: might need to be changed to be useful!
     */
    virtual double DomainSize() const
    {
        return Area();
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( fabs( rResult[0] ) <= 1.0 )
        {
            if ( fabs( rResult[1] ) <= 1.0 )
            {
                return true;
            }
        }

        return false;
    }

    /**
    * Returns the local coordinates of a given arbitrary point
    */
    virtual CoordinatesArrayType& PointLocalCoordinates( CoordinatesArrayType& rResult,
            const CoordinatesArrayType& rPoint )
    {
        #ifdef DEBUG_LEVEL2
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        int LocalDim = BaseType::LocalSpaceDimension();

        Matrix J = ZeroMatrix( LocalDim, LocalDim );

        if ( rResult.size() != LocalDim )
            rResult.resize( LocalDim );

        //starting with xi = 0
        rResult = ZeroVector( LocalDim );

        Vector DeltaXi = ZeroVector( LocalDim );

        CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );

        //Newton iteration:
        double tol = 1.0e-9;

        int maxiter = 300;

        for ( int k = 0; k < maxiter; ++k )
        {
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(k)
            #endif
            noalias( CurrentGlobalCoords ) = ZeroVector( 3 );
            BaseType::GlobalCoordinates( CurrentGlobalCoords, rResult );
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(CurrentGlobalCoords)
            #endif
            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(CurrentGlobalCoords)
            #endif
            InverseOfJacobian( J, rResult );
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(J)
            #endif
            noalias( DeltaXi ) = prod( J, CurrentGlobalCoords );
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(DeltaXi)
            #endif
            noalias( rResult ) += DeltaXi;
            #ifdef DEBUG_LEVEL2
            KRATOS_WATCH(rResult)
            KRATOS_WATCH(MathUtils<double>::Norm( DeltaXi ))
            #endif

            if ( MathUtils<double>::Norm( DeltaXi ) < tol )
            {
                break;
            }
        }

        #ifdef DEBUG_LEVEL2
        KRATOS_WATCH(rResult)
        #endif

        return( rResult );
    }

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
        //getting derivatives of shape functions
        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
        //getting values of shape functions
        const Matrix& shape_functions_values = BaseType::ShapeFunctionsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 2, 2 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;

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
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod,
            Matrix & DeltaPosition ) const
    {
        //getting derivatives of shape functions
        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
        //getting values of shape functions
        const Matrix& shape_functions_values = BaseType::ShapeFunctionsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            //defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 2, 2 );
            //loop over all nodes
            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
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
        //setting up size of jacobian matrix
        rResult.resize( 2, 2 );
        //derivatives of shape functions
        const ShapeFunctionsGradientsType& shape_functions_gradients = BaseType::ShapeFunctionsLocalGradients( ThisMethod );
        const Matrix& ShapeFunctionsGradientInIntegrationPoint = shape_functions_gradients( IntegrationPointIndex );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        //loop over all nodes
        double j00 = 0.0;
        double j01 = 0.0;
        double j10 = 0.0;
        double j11 = 0.0;
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j00 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            j01 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
            j10 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            j11 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
        }

        rResult( 0, 0 ) = j00;
        rResult( 0, 1 ) = j01;
        rResult( 1, 0 ) = j10;
        rResult( 1, 1 ) = j11;

        return rResult;
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
        //setting up size of jacobian matrix
        rResult.resize( 2, 2 );
        //derivatives of shape functions
        Matrix shape_functions_gradients;
        shape_functions_gradients = this->ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );
        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        //loop over all nodes

        double j00 = 0.0;
        double j01 = 0.0;
        double j10 = 0.0;
        double j11 = 0.0;
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j00 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
            j01 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 1 ) );
            j10 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
            j11 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 1 ) );
        }

        rResult( 0, 0 ) = j00;
        rResult( 0, 1 ) = j01;
        rResult( 1, 0 ) = j10;
        rResult( 1, 1 ) = j11;

        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Determinant of jacobians for given integration method.
     * This method calculate determinant of jacobian in all
     * integrations points of given integration method.
     *
     * @return Vector of double which is vector of determinants of
     * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    virtual Vector& DeterminantOfJacobian( Vector& rResult,
            IntegrationMethod ThisMethod ) const
    {
        //workaround by riccardo
        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            //                Vector temp(this->IntegrationPointsNumber(ThisMethod));
            //                rResult.swap(temp);
        }

        //for all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            rResult[pnt] = DeterminantOfJacobian( pnt, ThisMethod );
        }

        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
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
    virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex,
            IntegrationMethod ThisMethod ) const
    {
        Matrix jacobian = ZeroMatrix( 2, 2 );
        jacobian = Jacobian( jacobian, IntegrationPointIndex, ThisMethod );
        return(( jacobian( 0, 0 ) * jacobian( 1, 1 ) ) - ( jacobian( 0, 1 ) * jacobian( 1, 0 ) ) );
    }

    /**
     * :TODO: TO BE TESTED
     */
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
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
    {
        Matrix jacobian = ZeroMatrix( 2, 2 );
        jacobian = Jacobian( jacobian, rPoint );
        return(( jacobian( 0, 0 ) * jacobian( 1, 1 ) ) - ( jacobian( 0, 1 ) * jacobian( 1, 0 ) ) );
    }

    /**
     * :TODO: TO BE TESTED
     */
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
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod ) const
    {
        //workaround by riccardo
        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); ++pnt )
        {
            Matrix tempMatrix( 2, 2 );
            rResult[pnt] = InverseOfJacobian( tempMatrix, pnt, ThisMethod );
        }

        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
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
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual Matrix& InverseOfJacobian( Matrix& rResult,
            IndexType IntegrationPointIndex,
            IntegrationMethod ThisMethod ) const
    {
        //current jacobian
        Matrix tempMatrix = ZeroMatrix( 2, 2 );
        tempMatrix = Jacobian( tempMatrix, IntegrationPointIndex, ThisMethod );
        //determinant of jacobian
        double det_j = DeterminantOfJacobian( IntegrationPointIndex, ThisMethod );
        //checking for singularity

        if ( det_j == 0.00 )
        {
            KRATOS_THROW_ERROR( std::runtime_error,
                    "Zero determinant of jacobian during inversion of matrix!" ,
                    *this );
        }

        //setting up result matrix
        rResult.resize( 2, 2 );

        //filling matrix
        rResult( 0, 0 ) = ( tempMatrix( 1, 1 ) ) / ( det_j );

        rResult( 1, 0 ) = -( tempMatrix( 1, 0 ) ) / ( det_j );

        rResult( 0, 1 ) = -( tempMatrix( 0, 1 ) ) / ( det_j );

        rResult( 1, 1 ) = ( tempMatrix( 0, 0 ) ) / ( det_j );

        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
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
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
        //current jacobian
        Matrix tempMatrix = ZeroMatrix( 2, 2 );
        tempMatrix = Jacobian( tempMatrix, rPoint );
        //deteminant of Jacobian
        double det_j = DeterminantOfJacobian( rPoint );
        //checking for singularity

        if ( det_j == 0.00 )
        {
            KRATOS_THROW_ERROR( std::runtime_error,
                    "Zero determinant of jacobian during inversion of matrix!",
                    *this );
        }

        //setting up result matrix
        rResult.resize( 2, 2 );

        //filling matrix
        rResult( 0, 0 ) = ( tempMatrix( 1, 1 ) ) / ( det_j );

        rResult( 1, 0 ) = -( tempMatrix( 1, 0 ) ) / ( det_j );

        rResult( 0, 1 ) = -( tempMatrix( 0, 1 ) ) / ( det_j );

        rResult( 1, 1 ) = ( tempMatrix( 0, 0 ) ) / ( det_j );

        return rResult;
    }

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
     * TODO: implemented but not yet tested
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
            const CoordinatesArrayType& rPoint ) const
    {
        #ifdef DEBUG_LEVEL1
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        int Index1 = ShapeFunctionIndex / mNumber2;
        int Index2 = ShapeFunctionIndex % mNumber2;

        int Span1 = BSplineUtils::FindSpan(mNumber1, mOrder1, rPoint[0], mKnots1);
        int Span2 = BSplineUtils::FindSpan(mNumber2, mOrder2, rPoint[1], mKnots2);

        #ifdef DEBUG_LEVEL1
        KRATOS_WATCH(Span1)
        KRATOS_WATCH(Span2)
        #endif

        int Start1 = Span1 - mOrder1;
        int Start2 = Span2 - mOrder2;

        // bound checking
        if((Index1 - Start1 > mOrder1) || (Index1 - Start1 < 0))
        {
            return 0.0;
        }
        if((Index2 - Start2 > mOrder2) || (Index2 - Start2 < 0))
        {
            return 0.0;
        }

        ValuesContainerType ShapeFunctionValues1(mOrder1 + 1);
        ValuesContainerType ShapeFunctionValues2(mOrder2 + 1);

        BSplineUtils::BasisFuns(ShapeFunctionValues1, Span1, rPoint[0], mOrder1, mKnots1);
        BSplineUtils::BasisFuns(ShapeFunctionValues2, Span2, rPoint[1], mOrder2, mKnots2);

        #ifdef DEBUG_LEVEL1
        KRATOS_WATCH(ShapeFunctionValues1)
        KRATOS_WATCH(ShapeFunctionValues2)
        KRATOS_WATCH(Start1)
        KRATOS_WATCH(Start2)
        KRATOS_WATCH(mCtrlWeights)
        KRATOS_WATCH(mNumber2)
        #endif

        double Denom = 0.0;
        for(unsigned int i = Start1; i <= Span1; ++i) //i \in {0, ... , mNumber1 -1}
        {
            for(unsigned int j = Start2; j <= Span2; ++j) //j \in {0, ... , mNumber2 -1}
            {
                Denom += mCtrlWeights[i * mNumber2 + j] * ShapeFunctionValues1[i - Start1] * ShapeFunctionValues2[j - Start2];
            }
        }

        #ifdef DEBUG_LEVEL1
        KRATOS_WATCH(Denom)
        KRATOS_WATCH(Index1)
        KRATOS_WATCH(Index2)
        #endif

        return ShapeFunctionValues1[Index1 - Start1] *
                ShapeFunctionValues2[Index2 - Start2] *
                 ( mCtrlWeights[ShapeFunctionIndex] / Denom );
    }

    /**
     * Calculates the value of all shape functions at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     */
    virtual Vector& ShapeFunctionValues( Vector& rResults,
        const CoordinatesArrayType& rCoordinates ) const
    {
        rResults.resize( mNumber1 * mNumber2 );
        noalias( rResults ) = ZeroVector( mNumber1 * mNumber2 );

        //compute the b-spline shape functions
        ValuesContainerType ShapeFunctionValues1(mOrder1 + 1);
        ValuesContainerType ShapeFunctionValues2(mOrder2 + 1);

        int Span1 = BSplineUtils::FindSpan(mNumber1, mOrder1, rCoordinates[0], mKnots1);
        int Span2 = BSplineUtils::FindSpan(mNumber2, mOrder2, rCoordinates[1], mKnots2);

        BSplineUtils::BasisFuns(ShapeFunctionValues1, Span1, rCoordinates[0], mOrder1, mKnots1);
        BSplineUtils::BasisFuns(ShapeFunctionValues2, Span2, rCoordinates[1], mOrder2, mKnots2);

        int Start1 = Span1 - mOrder1;
        int Start2 = Span2 - mOrder2;

        double Denom = 0.0;
        double N1;
        double N2;
        double W;

        unsigned int i, j, Index;
        for(i = Start1; i <= Span1; ++i)
        {
            for(j = Start2; j <= Span2; ++j)
            {
                Index = i * mNumber2 + j; // TODO check this again. I'm skeptical.

                W = mCtrlWeights[Index];
                N1 = ShapeFunctionValues1(i - Start1);
                N2 = ShapeFunctionValues2(j - Start2);

                rResults(Index) = W * N1 * N2;

                Denom += rResults(Index);
            }
        }

        rResults *= (1.0 / Denom);

        return rResults;
    }

    /**
     * Calculates the local gradients at a given point
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const
    {
        // setting up result matrix
        rResult.resize( mNumber1 * mNumber2, 2 );
        noalias( rResult ) = ZeroMatrix( mNumber1 * mNumber2, 2 );

        // compute the b-spline shape functions & first derivatives
        const int NumberOfDerivatives = 1;
        Matrix ShapeFunctionsValuesAndDerivatives1(NumberOfDerivatives + 1, mOrder1 + 1);
        Matrix ShapeFunctionsValuesAndDerivatives2(NumberOfDerivatives + 1, mOrder2 + 1);
        int Span1 = BSplineUtils::FindSpan(mNumber1, mOrder1, rPoint[0], mKnots1);
        int Span2 = BSplineUtils::FindSpan(mNumber2, mOrder2, rPoint[1], mKnots2);
        int Start1 = Span1 - mOrder1;
        int Start2 = Span2 - mOrder2;
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives1, Span1, rPoint[0], mOrder1, mKnots1, NumberOfDerivatives, BSplineUtils::MatrixOp());
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives2, Span2, rPoint[1], mOrder2, mKnots2, NumberOfDerivatives, BSplineUtils::MatrixOp());
        double Denom = 0.0;
        double Denom_der1 = 0.0;
        double Denom_der2 = 0.0;
        double N1, dN1;
        double N2, dN2;
        double W;

        unsigned int i, j, Index;
        for(i = Start1; i <= Span1; ++i)
        {
            for(j = Start2; j <= Span2; ++j)
            {
                Index = i * mNumber2 + j;

                W = mCtrlWeights[Index];
                N1 = ShapeFunctionsValuesAndDerivatives1(0, i - Start1);
                dN1 = ShapeFunctionsValuesAndDerivatives1(1, i - Start1);
                N2 = ShapeFunctionsValuesAndDerivatives2(0, j - Start2);
                dN2 = ShapeFunctionsValuesAndDerivatives2(1, j - Start2);

                Denom += W * N1 * N2;
                Denom_der1 += W * dN1 * N2;
                Denom_der2 += W * N1 * dN2;
            }
        }

        for(i = Start1; i <= Span1; ++i)
        {
            for(j = Start2; j <= Span2; ++j)
            {
                Index = i * mNumber2 + j;

                W = mCtrlWeights[Index];
                N1 = ShapeFunctionsValuesAndDerivatives1(0, i - Start1);
                dN1 = ShapeFunctionsValuesAndDerivatives1(1, i - Start1);
                N2 = ShapeFunctionsValuesAndDerivatives2(0, j - Start2);
                dN2 = ShapeFunctionsValuesAndDerivatives2(1, j - Start2);

                rResult(Index, 0) = W * N2 * (dN1 * Denom - Denom_der1 * N1) / pow(Denom, 2);
                rResult(Index, 1) = W * N1 * (dN2 * Denom - Denom_der2 * N2) / pow(Denom, 2);
            }
        }

        return( rResult );
    }

    void ShapeFunctionsValuesAndLocalGradients
    (
        Vector& shape_functions_values,
        Matrix& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        //setting up result matrix
        shape_functions_local_gradients.resize( mNumber1 * mNumber2, 2 );
        noalias( shape_functions_local_gradients ) = ZeroMatrix( mNumber1 * mNumber2, 2 );
        shape_functions_values.resize( mNumber1 * mNumber2 );
        noalias( shape_functions_values ) = ZeroVector( mNumber1 * mNumber2 );

        #ifdef DEBUG_LEVEL3
        KRATOS_WATCH(mNumber1)
        KRATOS_WATCH(mNumber2)
        KRATOS_WATCH(mOrder1)
        KRATOS_WATCH(mOrder2)
        KRATOS_WATCH(mKnots1)
        KRATOS_WATCH(mKnots2)
        KRATOS_WATCH(rPoint)
        #endif

        //compute the b-spline shape functions & first derivatives
        const int NumberOfDerivatives = 1;
        Matrix ShapeFunctionsValuesAndDerivatives1(NumberOfDerivatives + 1, mOrder1 + 1);
        Matrix ShapeFunctionsValuesAndDerivatives2(NumberOfDerivatives + 1, mOrder2 + 1);
        int Span1 = BSplineUtils::FindSpan(mNumber1, mOrder1, rPoint[0], mKnots1);
        int Span2 = BSplineUtils::FindSpan(mNumber2, mOrder2, rPoint[1], mKnots2);
        int Start1 = Span1 - mOrder1;
        int Start2 = Span2 - mOrder2;

        #ifdef DEBUG_LEVEL3
        KRATOS_WATCH(Span1)
        KRATOS_WATCH(Span2)
        #endif

        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives1, Span1, rPoint[0], mOrder1, mKnots1, NumberOfDerivatives, BSplineUtils::MatrixOp());
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives2, Span2, rPoint[1], mOrder2, mKnots2, NumberOfDerivatives, BSplineUtils::MatrixOp());
        double Denom = 0.0;
        double Denom_der1 = 0.0;
        double Denom_der2 = 0.0;
        double N1, dN1;
        double N2, dN2;
        double W;

        unsigned int i, j, Index;
        for(i = Start1; i <= Span1; ++i)
        {
            for(j = Start2; j <= Span2; ++j)
            {
                Index = i * mNumber2 + j;

                W = mCtrlWeights[Index];
                N1 = ShapeFunctionsValuesAndDerivatives1(0, i - Start1);
                dN1 = ShapeFunctionsValuesAndDerivatives1(1, i - Start1);
                N2 = ShapeFunctionsValuesAndDerivatives2(0, j - Start2);
                dN2 = ShapeFunctionsValuesAndDerivatives2(1, j - Start2);

                Denom += W * N1 * N2;
                Denom_der1 += W * dN1 * N2;
                Denom_der2 += W * N1 * dN2;
            }
        }

        for(i = Start1; i <= Span1; ++i)
        {
            for(j = Start2; j <= Span2; ++j)
            {
                Index = i * mNumber2 + j;

                W = mCtrlWeights[Index];
                N1 = ShapeFunctionsValuesAndDerivatives1(0, i - Start1);
                dN1 = ShapeFunctionsValuesAndDerivatives1(1, i - Start1);
                N2 = ShapeFunctionsValuesAndDerivatives2(0, j - Start2);
                dN2 = ShapeFunctionsValuesAndDerivatives2(1, j - Start2);

                shape_functions_values(Index) = W * N1 * N2 / Denom;
                shape_functions_local_gradients(Index, 0) = W * N2 * (dN1 * Denom - Denom_der1 * N1) / pow(Denom, 2);
                shape_functions_local_gradients(Index, 1) = W * N1 * (dN2 * Denom - Denom_der2 * N2) / pow(Denom, 2);
            }
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
     *
     * :TODO: TESTING!!!
     */
    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
            ShapeFunctionsGradientsType& rResult,
            IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number = mpGeometryData->IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
        {
            KRATOS_THROW_ERROR( std::logic_error, "This integration method is not supported" , *this );
        }

        //workaround by riccardo
        if ( rResult.size() != integration_points_number )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( integration_points_number );
            rResult.swap( temp );
        }

        //calculating the local gradients
        const ShapeFunctionsGradientsType& locG = BaseType::ShapeFunctionsLocalGradients( ThisMethod );

        //getting the inverse jacobian matrices
        JacobiansType temp( integration_points_number );

        JacobiansType invJ = InverseOfJacobian( temp, ThisMethod );

        unsigned int i, j, pnt;
        // loop over all integration points
        for ( pnt = 0; pnt < integration_points_number; ++pnt )
        {
            rResult[pnt].resize( mNumber1 * mNumber2, 2 );

            for ( i = 0; i < mNumber1 * mNumber2; ++i )
            {
                for ( j = 0; j < 2; ++j )
                {
                    rResult[pnt]( i, j ) = ( locG[pnt]( i, 0 ) * invJ[pnt]( j, 0 ) ) + ( locG[pnt]( i, 1 ) * invJ[pnt]( j, 1 ) );
                }
            }
        }// end of loop over integration points

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
    virtual std::string Info() const
    {
        return "2 dimensional NURBS surface in 2D space";
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

    /**
     * TO BE CALLED BY ELEMENT
     */
    virtual void GenerateGeometryData(
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2,
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator, //not used
        const int& Degree1,
        const int& Degree2,
        const int& Degree3, //not used
        const int& NumberOfIntegrationMethod
    )
    {
        mKnots1 = Knots1;
        mKnots2 = Knots2;
        mCtrlWeights = Weights;
        mOrder1 = Degree1;
        mOrder2 = Degree2;
        mNumber1 = Knots1.size() - Degree1 - 1;
        mNumber2 = Knots2.size() - Degree2 - 1;

        if(mNumber1 * mNumber2 != this->size())
        {
            KRATOS_WATCH(Knots1.size())
            KRATOS_WATCH(Knots2.size())
            KRATOS_WATCH(Degree1)
            KRATOS_WATCH(Degree2)
            KRATOS_THROW_ERROR(std::logic_error, "The parametric parameters is not compatible, knots.length != n+p+1.", __FUNCTION__)
        }

        #ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        #endif

        //generate all integration points
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(NumberOfIntegrationMethod);

        #ifdef DEBUG_LEVEL3
        std::cout << "Generate integration_points completed." << std::endl;
        KRATOS_WATCH(all_integration_points.size())
//        for(int i = 0; i < all_integration_points.size(); i++)
//            KRATOS_WATCH(all_integration_points[i])
        #endif

        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Generate integration points: " << end_compute - start_compute << " s" << std::endl;
        start_compute = end_compute;
        #endif

        //generate all shape function values and derivatives
        ShapeFunctionsValuesContainerType shape_functions_values;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        for (unsigned int i = 0; i < NumberOfIntegrationMethod; ++i)
        {
            CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients
            (
                shape_functions_values[i],
                shape_functions_local_gradients[i],
                all_integration_points[i]
            );
        }

        #ifdef ENABLE_PROFILING
        end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Generate shape function values and derivatives: " << end_compute - start_compute << " s" << std::endl;
        start_compute = end_compute;
        #endif

        #ifdef DEBUG_LEVEL3
        std::cout << "Generate shape functions and local gradients completed." << std::endl;
        for(int i = 0; i < shape_functions_values.size(); i++)
            KRATOS_WATCH(shape_functions_values[i])
        for(int i = 0; i < shape_functions_local_gradients.size(); i++)
            KRATOS_WATCH(shape_functions_local_gradients[i])
        #endif

        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
                new GeometryData(
                        2,        //ThisDimension
                        2,//ThisWorkingSpaceDimension
                        2,//ThisLocalSpaceDimension
                        GeometryData::GI_GAUSS_1,//ThisDefaultMethod
                        all_integration_points,//ThisIntegrationPoints
                        shape_functions_values,//ThisShapeFunctionsValues
                        shape_functions_local_gradients//ThisShapeFunctionsLocalGradients
                )
        );

        #ifdef ENABLE_PROFILING
        end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Initialize pNewGeometryData: " << end_compute - start_compute << " s" << std::endl;
        start_compute = end_compute;

//        for ( unsigned int i = 0 ; i < GeometryData::NumberOfIntegrationMethods ; i++ )
//        {
//            boost::numeric::ublas::vector<Matrix> temp( all_integration_points[i].size() );

//            for ( unsigned int j = 0 ; j < all_integration_points[i].size() ; j++ )
//                temp[j] = outer_prod( row( shape_functions_values[i], j ), row( shape_functions_values[i], j ) );
//        }
//
//        end_compute = OpenMPUtils::GetCurrentTime();
//        std::cout << "Dummy calculation: " << end_compute - start_compute << " s" << std::endl;
//        start_compute = end_compute;
        #endif


        //generate an empty GeometryData
//        IntegrationPointsContainerType integration_points =
//        {};
//        ShapeFunctionsValuesContainerType shape_functions_values =
//        {};
//        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
//        {};
//        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
//                new GeometryData(
//                        3,
//                        3,
//                        1,
//                        GeometryData::GI_GAUSS_1,
//                        integration_points,
//                        shape_functions_values,
//                        shape_functions_local_gradients
//                )
//        );

        mpGeometryData.swap(pNewGeometryData);
        BaseType::mpGeometryData = &(*mpGeometryData);
    }

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */

    GeometryData::Pointer mpGeometryData;

    ValuesContainerType mKnots1; //knots vector
    ValuesContainerType mKnots2;//knots vector

    ValuesContainerType mCtrlWeights;//weight of control points

    int mOrder1;//order of the surface at parametric direction 1
    int mOrder2;//order of the surface at parametric direction 2

    int mNumber1;//number of shape functions define the surface on parametric direction 1
    int mNumber2;//number of shape functions define the surface on parametric direction 2

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // method to build GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     * KLUDGE: values are hard-coded!
     */
    Matrix CalculateShapeFunctionsIntegrationPointsValues(IntegrationMethod ThisMethod ) const
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(ThisMethod);
        const IntegrationPointsArrayType& integration_points = all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, mNumber1 * mNumber2 );
        //loop over all integration points

        //TODO: this can be optimized
        for(unsigned int pnt = 0; pnt < integration_points_number; ++pnt)
        {
            for(unsigned int node = 0; node < mNumber1 * mNumber2; ++node)
            {
                shape_function_values( pnt, node ) = ShapeFunctionValue(node, integration_points[pnt]);
            }
        }

        return shape_function_values;
    }

    /**
     * TODO: implemented but not yet tested
     */
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
    ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(IntegrationMethod ThisMethod ) const
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(ThisMethod);
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
        ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
        std::fill( DN_De.begin(), DN_De.end(), Matrix() );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            ShapeFunctionsLocalGradients(DN_De[it_gp], IntegrationPoints[it_gp]);
        }

        return DN_De;
    }

    /**
     * TODO
     */
    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        Matrix& shape_function_values,
        ShapeFunctionsGradientsType& shape_function_local_gradients,
        const IntegrationPointsArrayType& integration_points
    ) const
    {
        #ifdef DEBUG_LEVEL3
        std::cout << typeid(*this).name() << "::" << __FUNCTION__ << std::endl;
        #endif

        shape_function_local_gradients.resize(integration_points.size());
        std::fill(shape_function_local_gradients.begin(), shape_function_local_gradients.end(), Matrix());
        shape_function_values.resize(integration_points.size(), mNumber1 * mNumber2);

        #ifdef DEBUG_LEVEL3
        KRATOS_WATCH(integration_points.size())
        #endif

        for (unsigned int it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            Vector temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                temp_values,
                shape_function_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for(unsigned int node = 0; node < mNumber1 * mNumber2; ++node)
            {
                shape_function_values( it_gp, node ) = temp_values(node);
            }
        }
    }

    void FilterUniqueKnots(ValuesContainerType& UnrepeatedKnots, const ValuesContainerType& Knots) const
    {
        UnrepeatedKnots.resize(Knots.size(), false);

        UnrepeatedKnots[0] = Knots[0];
        double knot = UnrepeatedKnots[0];
        int cnt = 1;
        for(unsigned int i = 1; i < Knots.size(); ++i)
        {
            if(Knots[i] > knot)
            {
                knot = Knots[i];
                UnrepeatedKnots[cnt++] = knot;
            }
        }

        UnrepeatedKnots.resize(cnt, true);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // end of method to build to GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // method to construct GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * TODO: TO BE VERIFIED
     */
    virtual IntegrationPointsContainerType AllIntegrationPoints(int NumberOfIntegrationMethod) const
    {
        //check the size of relevant integration method
        if( (NumberOfIntegrationMethod > GeometryData::NumberOfIntegrationMethods - mOrder1 + 1) ||
                (NumberOfIntegrationMethod > GeometryData::NumberOfIntegrationMethods - mOrder2 + 1) )
        {
            KRATOS_WATCH(GeometryData::NumberOfIntegrationMethods)
            KRATOS_WATCH(mOrder1)
            KRATOS_WATCH(mOrder2)
            KRATOS_THROW_ERROR(std::logic_error, "Number of integration methods exceeds the allowance defined by boost.array", __FUNCTION__)
        }

        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

        //define the base integration rule for 1st parametric dimension. This ensures that the number of integration points in each direction equal to order + 1
        std::vector<IntegrationPointsArrayType> BaseRule;
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints6, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints7, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints8, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints9, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints10, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());

        ///////////////////////////////////////////////////////////////
        // Remarks: this current implementation supports integration with order up to 9

        ValuesContainerType UnrepeatedKnots1;
        ValuesContainerType UnrepeatedKnots2;
        this->FilterUniqueKnots(UnrepeatedKnots1, mKnots1);
        this->FilterUniqueKnots(UnrepeatedKnots2, mKnots2);

        KRATOS_WATCH(UnrepeatedKnots1)
        KRATOS_WATCH(UnrepeatedKnots2)

        unsigned int i1, i2, j1, j2, k, offset1, offset2;
        double a1, b1, left1, right1;
        double a2, b2, left2, right2;

        //TODO: parallelize this process
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset1 = k + mOrder1; // this allows for one more order than the order of the curve
            offset2 = k + mOrder2; // this allows for one more order than the order of the curve

            IntegrationPointsArrayType TempGaussRule;

            for(i1 = 0; i1 < UnrepeatedKnots1.size() - 1; ++i1)
            {
                left1 = UnrepeatedKnots1[i1];
                right1 = UnrepeatedKnots1[i1 + 1];

                a1 = (right1 - left1) / 2;
                b1 = (right1 + left1) / 2;

                for(i2 = 0; i2 < UnrepeatedKnots2.size() - 1; ++i2)
                {
                    left2 = UnrepeatedKnots2[i2];
                    right2 = UnrepeatedKnots2[i2 + 1];

                    a2 = (right2 - left2) / 2;
                    b2 = (right2 + left2) / 2;

                    for(j1 = 0; j1 < BaseRule[offset1].size(); ++j1)
                    {
                        IntegrationPointType& temp1 = BaseRule[offset1][j1];

                        for(j2 = 0; j2 < BaseRule[offset2].size(); ++j2)
                        {
                            // this allows for one more order than the order of the curve
                            IntegrationPointType& temp2 = BaseRule[offset2][j2];

                            IntegrationPointType temp;

                            temp.X() = a1 * temp1.X() + b1;
                            temp.Y() = a2 * temp2.X() + b2;
                            temp.Weight() = temp1.Weight() * temp2.Weight() * a1 * a2;

                            TempGaussRule.push_back(temp);
                        }
                    }
                }
            }
            GaussRule.push_back(TempGaussRule);
        }

        KRATOS_WATCH(NumberOfIntegrationMethod)

        IntegrationPointsContainerType integration_points;
        for (unsigned int k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            integration_points[k] = GaussRule[k];

//            std::cout << "GaussRule[" << k << "]: " << std::endl;
//            for(int i = 0; i < GaussRule[k].size(); ++i)
//            {
//                std::cout << " " << GaussRule[k][i] << std::endl;
//            }
//            std::cout << std::endl;
        }

        return integration_points;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // end of method to construct GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo2dNURBS;

    /**
     * Un accessible methods
     */

};    // Class Geo2dNURBS

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo2dNURBS<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo2dNURBS<TPointType>& rThis)
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

