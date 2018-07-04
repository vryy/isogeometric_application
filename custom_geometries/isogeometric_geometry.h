
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
//   Date:                $Date: 18 Feb 2014 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_GEOMETRY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_GEOMETRY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes
#include <boost/array.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/pointer_vector.h"
#include "integration/integration_point.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/math_utils.h"
#include "integration/quadrature.h"
#include "integration/line_gauss_legendre_integration_points.h"


namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

///IsogeometricGeometry base class.
/** As a base class IsogeometricGeometry has all the common
 * interface of Kratos' geometries. Also it contains array of
 * pointers to its points, reference to shape functions values in
 * all integrations points and also local gradients of shape
 * functions evaluated in all integrations points.
 *
 * IsogeometricGeometry is a template class with just one template parameter:
 * - TPointType which reperesent the type of the point this geometry
 * type contain and build on.
 *
 * @see Point
 * @see Node
 * @see Formulation
 * @see GeometryAndFormulationElement
 */
template<class TPointType>
class IsogeometricGeometry : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// This IsogeometricGeometry type.
    typedef IsogeometricGeometry<TPointType> GeometryType;

    /// Pointer definition of IsogeometricGeometry
    KRATOS_CLASS_POINTER_DEFINITION( IsogeometricGeometry );

    /** Base type for geometry.
    */
    typedef Geometry<TPointType> BaseType;


    /** The bounding box */
    /*typedef BoundingBox<TPointType, GeometryType>  BoundingBoxType; */

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef PointerVector<TPointType> PointsArrayType;

    /** Integration methods implemented in geometry.
    */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
    returning edges of the geometry.
     */
    typedef PointerVector<GeometryType> GeometriesArrayType;

    /** Redefinition of geometry template parameter TPointType as this geometry point type.
     */
    typedef typename BaseType::PointType PointType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef std::size_t IndexType;


    /** This typed used to return size or dimension in
    geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef std::size_t SizeType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /** This type used for representing an integration point in
    geometry. This integration point is a point with an
    additional weight component.
    */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry.
    */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
    integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef typename BaseType::JacobiansType JacobiansType;

    /** A third order tensor to hold shape functions'  gradients.
    ShapefunctionsGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** A third order tensor to hold shape functions' local second derivatives.
    ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /** A fourth order tensor to hold shape functions' local third order derivatives
     */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    typedef Vector ValuesContainerType;

    /**
     * Type of iterators
     */
    typedef typename BaseType::iterator                     iterator;
    typedef typename BaseType::const_iterator               const_iterator;
    typedef typename BaseType::reverse_iterator             reverse_iterator;
    typedef typename BaseType::const_reverse_iterator       const_reverse_iterator;
    typedef typename BaseType::ptr_iterator                 ptr_iterator;
    typedef typename BaseType::ptr_const_iterator           ptr_const_iterator;
    typedef typename BaseType::ptr_reverse_iterator         ptr_reverse_iterator;
    typedef typename BaseType::ptr_const_reverse_iterator   ptr_const_reverse_iterator;

    /**
     * Type of Matrix
     */
    typedef Matrix MatrixType;

    /**
     * Type of Vector
     */
    typedef Vector VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    IsogeometricGeometry() : BaseType()
    #ifndef ENABLE_PRECOMPUTE
    , mIsInitialized(false)
    #endif
    {
    }

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it usefull
    in the case of constructing new geometry without mapping and
    integrating properties.

    @param ThisPoints Vector of pointers to points which this
    geometry constructing on them. Points must have dimension
    equal or greater than working space dimension though there
    is no control on it.

    TODO
    */
    IsogeometricGeometry( const PointsArrayType& ThisPoints,
              GeometryData const* pThisGeometryData = 0 )
    : BaseType( ThisPoints, pThisGeometryData )
    #ifndef ENABLE_PRECOMPUTE
    , mIsInitialized(false)
    #endif
    {
    }

    /** Copy constructor.
    Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    IsogeometricGeometry( const IsogeometricGeometry& rOther )
    : BaseType( rOther )
    #ifndef ENABLE_PRECOMPUTE
    , mIsInitialized(false)
    #endif
    {
    }

    /** Copy constructor from a geometry with other point type.
    Construct this geometry as a copy of given geometry which
    has different type of points. The given goemetry's
    TOtherPointType* must be implicity convertible to this
    geometry PointType.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    template<class TOtherPointType> IsogeometricGeometry( IsogeometricGeometry<TOtherPointType> const & rOther )
    : BaseType( rOther.begin(), rOther.end() )
    #ifndef ENABLE_PRECOMPUTE
    , mIsInitialized(false)
    #endif
    {
    }

    /// Destructor. Do nothing!!!
    virtual ~IsogeometricGeometry() {}

    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    IsogeometricGeometry& operator=( const IsogeometricGeometry& rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /** Assignment operator for geometries with different point type.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    template<class TOtherPointType>
    IsogeometricGeometry& operator=( IsogeometricGeometry<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_generic_family;
    }

    virtual GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_generic_type;
    }

//    virtual typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
//    {
//        return BaseType::Pointer( new IsogeometricGeometry( ThisPoints, BaseType::mpGeometryData ) );
//    }

//    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
//    {
//        IsogeometricGeometry< Point<3> >::PointsArrayType NewPoints;

//        //making a copy of the nodes TO POINTS (not Nodes!!!)

//        for ( IndexType i = 0 ; i < this->size() ; i++ )
//            #if defined(KRATOS_SD_REF_NUMBER_2)
//            NewPoints.push_back(( *this )[i] );
//            #elif defined(KRATOS_SD_REF_NUMBER_3)
//            NewPoints.push_back(boost::make_shared< Point<3> >(( *this )[i]));
//            #endif

//        //NewPoints[i] = typename Point<3>::Pointer(new Point<3>(*mPoints[i]));

//        //creating a geometry with the new points
//        boost::shared_ptr< IsogeometricGeometry< Point<3> > > p_clone( new IsogeometricGeometry< Point<3> >( NewPoints ) );

//        p_clone->ClonePoints();

//        return p_clone;
//    }

//    virtual inline SizeType NURBS_Dimension() const
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
//    }

//    virtual inline SizeType NURBS_WorkingSpaceDimension() const
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
//    }

//    virtual inline SizeType NURBS_LocalSpaceDimension() const
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
//    }

    /**
     * TODO to be removed
     */
    virtual void GenerateGeometryData
    (
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2,
        const ValuesContainerType& Knots3,
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        const int& Degree1,
        const int& Degree2,
        const int& Degree3,
        const int& NumberOfIntegrationMethod
    )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling IsogeometricGeometry base class function", __FUNCTION__)
    }

    /**
     * Subroutine to pass in the data to the Bezier element. This subroutine shall be called from the element/condition
     */
    virtual void AssignGeometryData
    (
        const ValuesContainerType& Knots1,
        const ValuesContainerType& Knots2,
        const ValuesContainerType& Knots3,
        const ValuesContainerType& Weights,
        const MatrixType& ExtractionOperator,
        const int& Degree1,
        const int& Degree2,
        const int& Degree3,
        const int& NumberOfIntegrationMethod
    )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling IsogeometricGeometry base class function", __FUNCTION__)
    }

    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        const IntegrationPointsArrayType& integration_points
    ) const
    {
        if (shape_functions_values.size1() != integration_points.size()
            || shape_functions_values.size1() != this->PointsNumber())
        shape_functions_values.resize(integration_points.size(), this->PointsNumber(), false);

        if (shape_functions_local_gradients.size() != integration_points.size())
            shape_functions_local_gradients.resize(integration_points.size());

        VectorType tmp_values(this->PointsNumber());
        MatrixType tmp_local_gradients(this->PointsNumber(), this->WorkingSpaceDimension());
        for (std::size_t PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            this->ShapeFunctionsValuesAndLocalGradients(tmp_values, tmp_local_gradients, integration_points[PointNumber]);
            noalias(row(shape_functions_values, PointNumber)) = tmp_values;
            shape_functions_local_gradients[PointNumber] = tmp_local_gradients;
        }
    }

    virtual void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        IntegrationMethod ThisMethod
    ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling IsogeometricGeometry base class function", __FUNCTION__)
    }

    /**
     * TODO
     */
    virtual JacobiansType& Jacobian0( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling IsogeometricGeometry base class function", __FUNCTION__)
        return rResult;
    }

    virtual CoordinatesArrayType& GlobalCoordinates( CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates )
    {
        noalias( rResult ) = ZeroVector( 3 );

        Vector shape_functions_values;

        ShapeFunctionsValues(shape_functions_values, LocalCoordinates);

        for ( IndexType i = 0 ; i < this->size() ; ++i )
        {
            noalias( rResult ) += shape_functions_values( i ) * this->GetPoint( i );
        }

        return rResult;
    }

    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "This function is not available for isogeometric geometry" , "" );
        return rResult;
    }

    virtual void ExtractLocalCoordinates(PointsArrayType& rPoints)
    {
        KRATOS_THROW_ERROR( std::logic_error, "Calling base class function" , __FUNCTION__ );
    }


    /******************************************************
        OVERRIDE FROM GEOMETRY
    *******************************************************/
    virtual void Initialize(IntegrationMethod ThisMethod)
    {
        #ifndef ENABLE_PRECOMPUTE
        if(!mIsInitialized)
        {
            mpInternal_Ncontainer = boost::shared_ptr<Matrix>(new Matrix());
            mpInternal_DN_De = boost::shared_ptr<ShapeFunctionsGradientsType>(new ShapeFunctionsGradientsType());
            this->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(*mpInternal_Ncontainer, *mpInternal_DN_De, ThisMethod);
            mIsInitialized = true;
        }
        #endif
    }

    virtual void Initialize(const IntegrationPointsArrayType& integration_points)
    {
        #ifndef ENABLE_PRECOMPUTE
        if(!mIsInitialized)
        {
            mpInternal_Ncontainer = boost::shared_ptr<Matrix>(new Matrix());
            mpInternal_DN_De = boost::shared_ptr<ShapeFunctionsGradientsType>(new ShapeFunctionsGradientsType());
            this->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(*mpInternal_Ncontainer, *mpInternal_DN_De, integration_points);
            mIsInitialized = true;
        }
        #else
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not available in the PRECOMPUTE mode")
        #endif
    }

    virtual void Clean()
    {
        #ifndef ENABLE_PRECOMPUTE
        if (mIsInitialized)
        {
            mpInternal_DN_De.reset();
            mpInternal_Ncontainer.reset();
            mIsInitialized = false;
        }
        #endif
    }

    #ifndef ENABLE_PRECOMPUTE
    virtual const Matrix& ShapeFunctionsValues( IntegrationMethod ThisMethod )  const
    {
        return *mpInternal_Ncontainer;
    }

    virtual const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients( IntegrationMethod ThisMethod ) const
    {
        return *mpInternal_DN_De;
    }
    #else
    virtual const Matrix& ShapeFunctionsValues( IntegrationMethod ThisMethod )  const
    {
        return BaseType::ShapeFunctionsValues( ThisMethod );
    }

    virtual const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients( IntegrationMethod ThisMethod ) const
    {
        return BaseType::ShapeFunctionsLocalGradients( ThisMethod );
    }
    #endif

    virtual Vector& ShapeFunctionsValues( Vector& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Calling base class ShapeFunctionsValues method instead of derived class one. Please check the definition of derived class.", *this );
        return rResults;
    }

    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResults, const CoordinatesArrayType& rCoordinates ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Calling base class ShapeFunctionsLocalGradients method instead of derived class one. Please check the definition of derived class.", *this );
        return rResults;
    }

    /******************************************************
        END OVERRIDE FROM GEOMETRY
    *******************************************************/


    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << BaseType::Dimension()
               << " dimensional isogeometric geometry in "
               << BaseType::WorkingSpaceDimension()
               << "D space";
        return buffer.str();
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << BaseType::Dimension()
                 << " dimensional isogeometric geometry in "
                 << BaseType::WorkingSpaceDimension()
                 << "D space";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@}
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    /** Protected Constructor.
    Avoids object to be created Except for derived classes
    */


    ///@}



private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    #ifndef ENABLE_PRECOMPUTE
    bool mIsInitialized;
    boost::shared_ptr<ShapeFunctionsGradientsType> mpInternal_DN_De;
    boost::shared_ptr<Matrix> mpInternal_Ncontainer;
    #endif

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
//        rSerializer.save( "IsogeometricGeometry Data", mpGeometryData );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
//        rSerializer.load( "IsogeometricGeometry Data", const_cast<GeometryData*>( mpGeometryData ) );
    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    virtual void ShapeFunctionsValuesAndLocalGradients(
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling IsogeometricGeometry base class function", __FUNCTION__)
    }

    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class IsogeometricGeometry;

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class IsogeometricGeometry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   IsogeometricGeometry<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const IsogeometricGeometry<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_GEOMETRY_H_INCLUDED  defined


