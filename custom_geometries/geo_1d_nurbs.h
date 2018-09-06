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
//   Date:                $Date: 2013 Sep 7 $
//   Revision:            $Revision: 1.1 $
//
//
#if !defined(KRATOS_GEO_1D_NURBS_H_INCLUDED )
#define  KRATOS_GEO_1D_NURBS_H_INCLUDED

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
#include "custom_utilities/bspline_utils.h"
#include "integration/quadrature.h"
#include "integration/line_gauss_legendre_integration_points.h"

namespace Kratos
{
/**
 * A geometry representing NURBS curve
 */

template<class TPointType>
class Geo1dNURBS: public IsogeometricGeometry<TPointType>
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
     * Pointer definition of Geo1dNURBS
     */
    KRATOS_CLASS_POINTER_DEFINITION( Geo1dNURBS );

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
     * Type of values container
     */
    typedef typename BaseType::NormalType ValuesContainerType;

    /**
     * Type of Matrix
     */
    typedef typename BaseType::MatrixType MatrixType;

    /**
     * Type of Vector
     */
    typedef typename BaseType::VectorType VectorType;

    /**
     * Life Cycle
     */

    Geo1dNURBS(): BaseType( PointsArrayType() )
    {}

    Geo1dNURBS(
            const PointsArrayType& ThisPoints
    )
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
    Geo1dNURBS( Geo1dNURBS const& rOther )
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
    template<class TOtherPointType> Geo1dNURBS( Geo1dNURBS<TOtherPointType> const& rOther )
    : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Geo1dNURBS()
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
    Geo1dNURBS& operator=( const Geo1dNURBS& rOther )
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
    Geo1dNURBS& operator=( Geo1dNURBS<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */

    typename BaseType::BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::BaseType::Pointer( new Geo1dNURBS( ThisPoints ) );
    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
//        Geometry< Point<3> >::PointsArrayType NewPoints;
//        //making a copy of the nodes TO POINTS (not Nodes!!!)

//        for ( IndexType i = 0; i < this->Points().size(); ++i )
//        NewPoints.push_back( this->Points()[i] );

//        //creating a geometry with the new points
//        boost::shared_ptr< Geometry< Point<3> > >
//        p_clone( new Geo1dNURBS< Point<3> >( NewPoints ) );

//        p_clone->ClonePoints();

//        return p_clone;
//    }

    /**
     * TODO: TO BE CHECKED!!!!!!!!!!!
     */
    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        rResult.resize( 3, false );
        rResult[0] = 0.25;
        rResult[2] = 0.5;
        rResult[1] = 0.25;
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
        //TODO: reimplement to account for nurbs curve
        return 0.0;
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
        // Area is not relevant for 1d geometry
        return 0.0;
    }

    virtual double Volume() const
    {
        // Volume is not relevant for 1d geometry
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
        return Length();
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( fabs( rResult[0] ) < 1 + 1.0e-8 )
        return true;

        return false;
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
        ShapeFunctionsGradientsType shape_functions_gradients =
        CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        //getting values of shape functions
        Matrix shape_functions_values =
        CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            //defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 3, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
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
        ShapeFunctionsGradientsType shape_functions_gradients =
        CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        //getting values of shape functions
        Matrix shape_functions_values =
        CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            //defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 3, 1 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
            {
                jacobian( 0, 0 ) += ( this->GetPoint( i ).X() + DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() + DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 0 ) += ( this->GetPoint( i ).Z() + DeltaPosition(i,2) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
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
        rResult.resize( 3, 1 );
        //derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
        CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        Matrix ShapeFunctionsGradientInIntegrationPoint =
        shape_functions_gradients( IntegrationPointIndex );
        //values of shape functions in integration points
        vector<double> ShapeFunctionValuesInIntegrationPoint = ZeroVector( 3 );
        ShapeFunctionValuesInIntegrationPoint = row( CalculateShapeFunctionsIntegrationPointsValues( ThisMethod ),
                IntegrationPointIndex );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;
        double j1 = 0.0;
        double j2 = 0.0;
        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            j1 += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            j2 += ( this->GetPoint( i ).Z() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;
        rResult( 1, 0 ) = j1;
        rResult( 2, 0 ) = j2;

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
        rResult.resize( 3, 1 );
        //derivatives of shape functions
        Matrix shape_functions_gradients;
        shape_functions_gradients = ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        double j0 = 0.0;
        double j1 = 0.0;
        double j2 = 0.0;
        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            j0 += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
            j1 += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
            j2 += ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients( i, 0 ) );
        }

        rResult( 0, 0 ) = j0;
        rResult( 1, 0 ) = j1;
        rResult( 2, 0 ) = j2;

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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
        return 0.0;
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
        return 0.0;
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
        int span = BSplineUtils::FindSpan(mNumber, mOrder, rPoint[0], mKnots);
        int start = span - mOrder;

        // bound checking
        if(ShapeFunctionIndex - start > mOrder || ShapeFunctionIndex - start < 0)
        {
            return 0.0;
        }

        ValuesContainerType ShapeFunctionValues(mOrder + 1);

        BSplineUtils::BasisFuns(ShapeFunctionValues, span, rPoint[0], mOrder, mKnots);

        double denom = 0.0;
        for(unsigned int i = start; i <= span; ++i)
        {
            denom += mCtrlWeights[i] * ShapeFunctionValues[i - start];
        }

        return ShapeFunctionValues[ShapeFunctionIndex - start] *
                ( mCtrlWeights[ShapeFunctionIndex] / denom );
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
        rResults.resize( mNumber );
        noalias( rResults ) = ZeroVector( mNumber );

        //compute the b-spline shape functions
        ValuesContainerType ShapeFunctionValues1(mOrder + 1);

        int Span = BSplineUtils::FindSpan(mNumber, mOrder, rCoordinates[0], mKnots);

        BSplineUtils::BasisFuns(ShapeFunctionValues1, Span, rCoordinates[0], mOrder, mKnots);

        int Start = Span - mOrder;

        double Denom = 0.0;
        double N;
        double W;

        unsigned int i, Index;
        for(Index = Start; Index <= Span; ++Index)
        {
            W = mCtrlWeights[Index];
            N = ShapeFunctionValues1(Index - Start);

            rResults(Index) = W * N;

            Denom += rResults(Index);
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
        //setting up result matrix
        rResult.resize( mNumber, 1 );
        noalias( rResult ) = ZeroMatrix( mNumber, 1 );

        //compute the b-spline shape functions & first derivatives
        const int NumberOfDerivatives = 1;
        Matrix ShapeFunctionsValuesAndDerivatives(NumberOfDerivatives + 1, mOrder + 1);
        int span = BSplineUtils::FindSpan(mNumber, mOrder, rPoint[0], mKnots);
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives, span, rPoint[0], mOrder, mKnots, NumberOfDerivatives, BSplineUtils::MatrixOp());
        double denom = 0.0;
        double denom_der = 0.0;
        int start = span - mOrder;
        double N, dN, W;

        unsigned int i;
        for(i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            denom += W * N;
            denom_der += W * dN;
        }

        for(i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            rResult(i, 0) = W * (dN * denom - N * denom_der) / pow(denom, 2);
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
        //setting up result matrix
        shape_functions_local_gradients.resize(mNumber, 1);
        noalias( shape_functions_local_gradients ) = ZeroMatrix( mNumber, 1 );
        shape_functions_values.resize(mNumber);
        noalias( shape_functions_values ) = ZeroVector(mNumber);

        //compute the b-spline shape functions & first derivatives
        const int NumberOfDerivatives = 1;
        Matrix ShapeFunctionsValuesAndDerivatives(NumberOfDerivatives + 1, mOrder + 1);
        int span = BSplineUtils::FindSpan(mNumber, mOrder, rPoint[0], mKnots);
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives, span, rPoint[0], mOrder, mKnots, NumberOfDerivatives, BSplineUtils::MatrixOp());
        double denom = 0.0;
        double denom_der = 0.0;
        int start = span - mOrder;
        double N, dN, W;

        unsigned int i;
        for(i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            denom += W * N;
            denom_der += W * dN;
        }

        for(i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            shape_functions_local_gradients(i, 0) = W * (dN * denom - N * denom_der) / pow(denom, 2);
            shape_functions_values(i) = W * N / denom;
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
        KRATOS_THROW_ERROR( std::logic_error, "Jacobian is not square" , __FUNCTION__);

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
        return "1 dimensional NURBS curve in 3D space";
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
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * TO BE CALLED BY ELEMENT
     */
    virtual void GenerateGeometryData(
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2, //not used
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator, //not used
        const int& Degree1,
        const int& Degree2, //not used
        const int& Degree3, //not used
        const int& NumberOfIntegrationMethod
    )
    {
        mKnots = Knots1;
        mCtrlWeights = Weights;
        mOrder = Degree1;
        mNumber = Knots1.size() - Degree1 - 1;

        if(mNumber != this->size())
        {
            KRATOS_THROW_ERROR(std::logic_error, "The parametric parameters is not compatible, knots.length != n+p+1.", __FUNCTION__)
        }

        //TODO: calculate everything related to geometry data here
        //generate all integration points
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(NumberOfIntegrationMethod);

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

        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
            new GeometryData(
                    3,        //ThisDimension
                    3,//ThisWorkingSpaceDimension
                    1,//ThisLocalSpaceDimension
                    GeometryData::GI_GAUSS_1,//ThisDefaultMethod
                    all_integration_points,//ThisIntegrationPoints
                    shape_functions_values,//ThisShapeFunctionsValues
                    shape_functions_local_gradients//ThisShapeFunctionsLocalGradients
            )
        );

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

    ValuesContainerType mKnots; //knot vector

    ValuesContainerType mCtrlWeights;//weight of control points

    int mOrder;//order of the curve

    int mNumber;//number of shape functions define the curve

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
        Matrix shape_function_values( integration_points_number, mNumber );
        //loop over all integration points

        //TODO: this can be optimized
        for(unsigned int pnt = 0; pnt < integration_points_number; ++pnt)
        {
            for(unsigned int node = 0; node < mNumber; node++)
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
        std::fill( DN_De.begin(), DN_De.end(), Matrix( mNumber, 1 ) );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            ShapeFunctionsLocalGradients(DN_De[it_gp], IntegrationPoints[it_gp]);
        }

        return DN_De;
    }

    /**
     * TODO
     */
    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients
    (
        Matrix& shape_function_values,
        ShapeFunctionsGradientsType& shape_function_local_gradients,
        const IntegrationPointsArrayType& integration_points
    ) const
    {
        shape_function_local_gradients.resize(integration_points.size());
        std::fill(shape_function_local_gradients.begin(), shape_function_local_gradients.end(), Matrix());
        shape_function_values.resize(integration_points.size(), mNumber);

        for (unsigned int it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            Vector temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                temp_values,
                shape_function_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for(unsigned int node = 0; node < mNumber; ++node)
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
        if(NumberOfIntegrationMethod > GeometryData::NumberOfIntegrationMethods - mOrder + 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Number of integration methods exceeds the allowance defined by boost.array", __FUNCTION__)
        }

        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

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

        ValuesContainerType UnrepeatedKnots;
        this->FilterUniqueKnots(UnrepeatedKnots, mKnots);
        unsigned int i, j, k, offset;

        double a, b, left, right;
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset = k + mOrder;

            IntegrationPointsArrayType TempGaussRule;

            for(i = 0; i < UnrepeatedKnots.size() - 1; ++i)
            {
                left = UnrepeatedKnots[i];
                right = UnrepeatedKnots[i + 1];

                a = (right - left) / 2;
                b = (right + left) / 2;

                for(j = 0; j < BaseRule[offset].size(); ++j)
                {
                    IntegrationPointType temp = BaseRule[offset][j];
                    temp.X() = a * temp.X() + b;
                    temp.Weight() *= a;

                    TempGaussRule.push_back(temp);
                }
            }
            GaussRule.push_back(TempGaussRule);
        }

        IntegrationPointsContainerType integration_points;
        for (unsigned int k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            integration_points[k] = GaussRule[k];
        }

        return integration_points;
    }

    //This 2 methods works but invoking repetitive calculation of shape functions. This can be optimized by merging these functions in one call. These functions are kept as reference but never be called to compute geometry_data
    /**
     * TODO: TO BE VERIFIED
     */
    ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 )
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 )
            }
        };
        return shape_functions_local_gradients;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // end of method to construct GeometryData
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Geo1dNURBS;

    /**
     * Un accessible methods
     */

};    // Class Geo1dNURBS

/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >>(
        std::istream& rIStream, Geo1dNURBS<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator <<(
        std::ostream& rOStream, const Geo1dNURBS<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}    // namespace Kratos.

#endif

