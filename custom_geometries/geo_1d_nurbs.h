/*
see isogeometric_application/LICENSE.txt
 */

//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013 Sep 7 $
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
     * The original geometry type
     */
    typedef typename BaseType::GeometryType GeometryType;

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

    /** Type used to represent the coordinate, althought it is called DataType
     * It should be called CoordinateType, but many downstream geometries depend
     * on this type, I will leave it for now.
     */
    typedef typename PointType::DataType DataType;

    /** Type used to represent the real values, like shape function
     */
    typedef typename BaseType::ValueType ValueType;

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
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

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

    /** This type used for representing the local coordinates of
    an integration point
    */
    typedef typename BaseType::LocalCoordinatesArrayType LocalCoordinatesArrayType;

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
    ~Geo1dNURBS() override
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

    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename GeometryType::Pointer( new Geo1dNURBS( ThisPoints ) );
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
    Vector& LumpingFactors( Vector& rResult ) const override
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
    double Length() const override
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
    double Area() const override
    {
        // Area is not relevant for 1d geometry
        return 0.0;
    }

    double Volume() const override
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
    double DomainSize() const override
    {
        return Length();
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    bool IsInside( const LocalCoordinatesArrayType& rPoint, const ValueType tol ) const final
    {
        if ( fabs( rPoint[0] ) < 1 + tol )
        {
            return true;
        }

        return false;
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
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                               const CoordinatesArrayType& rPoint ) const override
    {
        int span = BSplineUtils::FindSpan(mNumber, mOrder, rPoint[0], mKnots);
        int start = span - mOrder;

        // bound checking
        if (ShapeFunctionIndex - start > mOrder || ShapeFunctionIndex - start < 0)
        {
            return 0.0;
        }

        ValuesContainerType ShapeFunctionValues(mOrder + 1);

        BSplineUtils::BasisFuns(ShapeFunctionValues, span, rPoint[0], mOrder, mKnots);

        double denom = 0.0;
        for (unsigned int i = start; i <= span; ++i)
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
    Vector& ShapeFunctionsValues( Vector& rResults,
                                  const LocalCoordinatesArrayType& rCoordinates ) const override
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
        for (Index = Start; Index <= Span; ++Index)
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
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
                                          const LocalCoordinatesArrayType& rPoint ) const override
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
        for (i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            denom += W * N;
            denom_der += W * dN;
        }

        for (i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            rResult(i, 0) = W * (dN * denom - N * denom_der) / pow(denom, 2);
        }

        return ( rResult );
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
        for (i = start; i <= span; ++i)
        {
            W = mCtrlWeights[i];
            N = ShapeFunctionsValuesAndDerivatives(0, i - start);
            dN = ShapeFunctionsValuesAndDerivatives(1, i - start);
            denom += W * N;
            denom_der += W * dN;
        }

        for (i = start; i <= span; ++i)
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
    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "To be implemented";

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
    std::string Info() const override
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
    void PrintInfo( std::ostream& rOStream ) const override
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
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        this->Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * TO BE CALLED BY ELEMENT
     */
    void GenerateGeometryData(
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2, //not used
        const ValuesContainerType& Knots3, //not used
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator, //not used
        int Degree1,
        int Degree2, //not used
        int Degree3, //not used
        int NumberOfIntegrationMethod
    ) override
    {
        mKnots = Knots1;
        mCtrlWeights = Weights;
        mOrder = Degree1;
        mNumber = Knots1.size() - Degree1 - 1;

        if (mNumber != this->size())
        {
            KRATOS_ERROR << "The parametric parameters is not compatible, knots.length != n+p+1.";
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
                        GeometryData::IntegrationMethod::GI_GAUSS_1,//ThisDefaultMethod
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
//                        GeometryData::IntegrationMethod::GI_GAUSS_1,
//                        integration_points,
//                        shape_functions_values,
//                        shape_functions_local_gradients
//                )
//        );

        mpGeometryData.swap(pNewGeometryData);
        GeometryType::SetGeometryData(mpGeometryData.get());
    }

private:

    /**
     * Member Variables
     */

    GeometryData::Pointer mpGeometryData;

    ValuesContainerType mKnots; // knot vector

    ValuesContainerType mCtrlWeights; // weight of control points

    int mOrder; // order of the curve

    int mNumber; // number of shape functions define the curve

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

    ///@}

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
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(static_cast<int>(ThisMethod));
        const IntegrationPointsArrayType& integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, mNumber );
        //loop over all integration points

        //TODO: this can be optimized
        for (unsigned int pnt = 0; pnt < integration_points_number; ++pnt)
        {
            for (unsigned int node = 0; node < mNumber; node++)
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
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints(static_cast<int>(ThisMethod));
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[static_cast<int>(ThisMethod)];
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
            for (unsigned int node = 0; node < mNumber; ++node)
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
        for (unsigned int i = 1; i < Knots.size(); ++i)
        {
            if (Knots[i] > knot)
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
        if (NumberOfIntegrationMethod > static_cast<int>(GeometryData::IntegrationMethod::NumberOfIntegrationMethods) - mOrder + 1)
        {
            KRATOS_ERROR << "Number of integration methods too big";
        }

        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

        std::vector<IntegrationPointsArrayType> BaseRule;
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<1>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<2>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<3>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<4>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<5>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<6>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<7>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<8>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<9>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        BaseRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints<10>, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());

        ValuesContainerType UnrepeatedKnots;
        this->FilterUniqueKnots(UnrepeatedKnots, mKnots);
        unsigned int i, j, k, offset;

        double a, b, left, right;
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset = k + mOrder;

            IntegrationPointsArrayType TempGaussRule;

            for (i = 0; i < UnrepeatedKnots.size() - 1; ++i)
            {
                left = UnrepeatedKnots[i];
                right = UnrepeatedKnots[i + 1];

                a = (right - left) / 2;
                b = (right + left) / 2;

                for (j = 0; j < BaseRule[offset].size(); ++j)
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
                CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_1 )
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
                CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_1 )
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
