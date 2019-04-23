// see isogeometric_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014-01-28 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BEZIER_UTILS_H_INCLUDED )
#define  KRATOS_BEZIER_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "integration/quadrature.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "custom_utilities/isogeometric_math_utils.h"

#define ENABLE_PROFILING
#define USE_EQUAL_ORDER_INTEGRATION_IN_ALL_DIRECTION

namespace Kratos
{

///@{

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

/**
 * An implementation of multikey for choosing the correct integration rule
 * Reference: http://cnx.org/content/m35767/latest/
 */
class BezierGeometryDataKey
{
public:

    typedef std::size_t SizeType;

    BezierGeometryDataKey(SizeType NumberOfIntegrationMethod,
                          SizeType Order1 = 0,
                          SizeType Order2 = 0,
                          SizeType Order3 = 0,
                          SizeType Dimension = 3,
                          SizeType WorkingSpaceDimension = 3,
                          SizeType LocalSpaceDimension = 3)
    {
        mNumberOfIntegrationMethod = NumberOfIntegrationMethod;
        mOrder1 = Order1;
        mOrder2 = Order2;
        mOrder3 = Order3;
        mDimension = Dimension;
        mWorkingSpaceDimension = WorkingSpaceDimension;
        mLocalSpaceDimension = LocalSpaceDimension;
    }

    bool operator<(const BezierGeometryDataKey& Key) const
    {
        if ( mNumberOfIntegrationMethod == Key.mNumberOfIntegrationMethod )
        {
            if ( mOrder1 == Key.mOrder1 )
            {
                if ( mOrder2 == Key.mOrder2 )
                {
                    if ( mOrder3 ==  Key.mOrder3 )
                    {
                        if ( mDimension == Key.mDimension )
                        {
                            if ( mWorkingSpaceDimension == Key.mWorkingSpaceDimension )
                            {
                                return mLocalSpaceDimension < Key.mLocalSpaceDimension;
                            }
                            else
                                return mWorkingSpaceDimension < Key.mWorkingSpaceDimension;
                        }
                        else
                            return mDimension < Key.mDimension;
                    }
                    else
                        return mOrder3 < Key.mOrder3;
                }
                else
                {
                    return mOrder2 < Key.mOrder2;
                }
            }
            else
            {
                return mOrder1 < Key.mOrder1;
            }
        }
        else
        {
            return mNumberOfIntegrationMethod < Key.mNumberOfIntegrationMethod;
        }
    }

    void PrintInfo(std::ostream& os) const
    {
        os << "(NumberOfIntegrationMethod = " << mNumberOfIntegrationMethod
           << ", Order1 = " << mOrder1
           << ", Order2 = " << mOrder2
           << ", Order3 = " << mOrder3
           << ", Dimension = " << mDimension
           << ", WorkingSpaceDimension = " << mWorkingSpaceDimension
           << ", LocalSpaceDimension = " << mLocalSpaceDimension
           << ")";
    }

private:
    SizeType mNumberOfIntegrationMethod;
    SizeType mOrder1;
    SizeType mOrder2;
    SizeType mOrder3;
    SizeType mDimension;
    SizeType mWorkingSpaceDimension;
    SizeType mLocalSpaceDimension;
};

//Output function for BezierGeometryDataKey
inline std::ostream& operator<<(std::ostream& os, const BezierGeometryDataKey& Key)
{
    Key.PrintInfo(os);
    return os;
}

/// Short class definition.
/**
 * Detail class definition.
 */
class BezierUtils
{
public:

    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;
    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;
    typedef std::size_t IndexType;
    typedef std::map<BezierGeometryDataKey, GeometryData::Pointer> MapType;
    typedef std::pair<BezierGeometryDataKey, GeometryData::Pointer> PairType;
    typedef typename Element::GeometryType GeometryType;
    typedef typename GeometryType::PointType PointType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointType IntegrationPointType;
    typedef typename GeometryType::IntegrationMethod IntegrationMethod;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
    typedef typename ModelPart::ElementsContainerType ElementsArrayType;
    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;

    /// Pointer definition of BezierUtils
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(BezierUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BezierUtils()
    {}

    /// Destructor.
    virtual ~BezierUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /********************************************************
            Fundamental Bezier handling functions
     ********************************************************/

    /**
     * Computes Bernstein basis function B(i, p)(x) on [0, 1]
     * Remark when debugging hierarchical NURBS, it may be not correct
     */
    static inline double bernstein(const int& i, const int& p, const double& x)
    {
        //bound checking
        if(i < 0 || i > p)
            return 0.0;

        double a = x;
        double b = 1 - x;
        double coeff;

        if(p == 0)
        {
            return 1.0;
        }
        else if(p == 1)
        {
            coeff = 1.0;
        }
        else
        {
            int cnt;
            for(int j = 0; j < p - 2; ++j)
            {
                cnt += (j + 2) / 2;
            }
            coeff = msBernsteinCoefs[cnt - 1];
        }

        return coeff * pow(a, i) * pow(b, p - i);
    }

    /**
     * Computes Bernstein basis function B(i, p)(x) on [0, 1] using recursive iteration
     */
    static double bernstein2(const int& i, const int& p, const double& x)
    {
        if(i < 0 || i > p)
            return 0.0;

        if(p == 0)
        {
            return 1.0;
        }

        double a = x;
        double b = 1 - x;

        double tmp1 = bernstein2(i, p - 1, x);
        double tmp2 = bernstein2(i - 1, p - 1, x);
        return b * tmp1 + a * tmp2;
    }

    /**
     * Computes Bernstein basis function value & derivative B(i, p)(x) on [0, 1]
     */
    static inline void bernstein(double& v, double& d, const int& i, const int& p, const double& x)
    {
        if(i < 0 || i > p)
        {
            v = 0.0;
            d = 0.0;
            return;
        }

        if(p == 0)
        {
            v = 1.0;
            d = 0.0;
            return;
        }

        double a = x;
        double b = 1 - x;

        double tmp1 = bernstein2(i, p - 1, x);
        double tmp2 = bernstein2(i - 1, p - 1, x);

        v = b * tmp1 + a * tmp2;
        d = p * (tmp2 - tmp1);
    }

    /**
     * Computes Bernstein basis function value & derivative & second derivative B(i, p)(x) on [0, 1]
     */
    static inline void bernstein(double& v, double& d, double& d2, const int& i, const int& p, const double& x)
    {
        if(i < 0 || i > p)
        {
            v = 0.0;
            d = 0.0;
            d2 = 0.0;
            return;
        }

        if(p == 0)
        {
            v = 1.0;
            d = 0.0;
            d2 = 0.0;
            return;
        }

        if(p == 1)
        {
            bernstein(v, d, i, p, x);
            d2 = 0.0;
            return;
        }

        double a = x;
        double b = 1 - x;

        double tmp1 = bernstein2(i    , p - 1, x);
        double tmp2 = bernstein2(i - 1, p - 1, x);
        double tmp3 = bernstein2(i    , p - 2, x);
        double tmp4 = bernstein2(i - 1, p - 2, x);
        double tmp5 = bernstein2(i - 2, p - 2, x);

        v = b * tmp1 + a * tmp2;
        d = p * (tmp2 - tmp1);
        d2 = p * (p-1) * (tmp3 - 2 * tmp4 + tmp5);
    }

    template<class ValuesContainerType>
    static inline void bernstein(ValuesContainerType& rS, const int& p, const double& x)
    {
        double a = x;
        double b = 1 - x;

        if(p == 0)
        {
            rS[0] = 1.0;
            return;
        }
        else if(p == 1)
        {
            rS[0] = b;
            rS[1] = a;
            return;
        }
        else
        {
            rS[0] = 1.0;
            rS[p] = 1.0;
            int cnt = 0;

            for(int j = 0; j < p - 2; ++j)
            {
                cnt += (j + 2) / 2;
            }

            for(int j = 0; j < p / 2; ++j)
            {
                rS[j + 1] = msBernsteinCoefs[cnt + j];
            }

            for(int j = 0; j < (p - 1) / 2; ++j)
            {
                rS[p - j - 1] = rS[j + 1];
            }

            for(int j = 0; j < p + 1; ++j)
            {
                rS[j] *= (pow(a, j) * pow(b, p - j));
            }

            return;
        }
    }

    template<class ValuesContainerType>
    static inline void bernstein(ValuesContainerType& rS,
                                 ValuesContainerType& rD,
                                 const int& p,
                                 const double& x)
    {
        for(int i = 0; i < p + 1; ++i)
        {
            bernstein(rS[i], rD[i], i, p, x);
        }
    }

    template<class ValuesContainerType>
    static inline void bernstein(ValuesContainerType& rS,
                                 ValuesContainerType& rD,
                                 ValuesContainerType& rD2,
                                 const int& p,
                                 const double& x)
    {
        for(int i = 0; i < p + 1; ++i)
        {
            bernstein(rS[i], rD[i], rD2[i], i, p, x);
        }
    }

    /********************************************************
            End of Fundamental Bezier handling functions
     ********************************************************/

    /********************************************************
            Bezier integration utilities
     ********************************************************/
    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static void RegisterIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order
    )
    {
        //define the key
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order, 0, 0, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);

        //find the key in existing map
        if(mIntegrationMethods.find(Key) != mIntegrationMethods.end())
        //found the key
        {
//            std::cout << "Key " << Key << " existed and ready to be used" << std::endl;
        }
        else
        {
            //created integration rule and insert the key
            //define the integration rule
            IntegrationPointsContainerType all_integration_points
                = AllIntegrationPoints(NumberOfIntegrationMethod, Order);

            ShapeFunctionsValuesContainerType shape_functions_values;
            ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

            for (IndexType i = 0; i < NumberOfIntegrationMethod; ++i)
            {
                CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                    Order,
                    shape_functions_values[i],
                    shape_functions_local_gradients[i],
                    all_integration_points[i]
                );
            }

            //create the geometry_data pointer
            GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
                new GeometryData(
                    TDimension,
                    TWorkingSpaceDimension,
                    TLocalSpaceDimension,
                    GeometryData::GI_GAUSS_2,           //ThisDefaultMethod
                    all_integration_points,             //ThisIntegrationPoints
                    shape_functions_values,             //ThisShapeFunctionsValues
                    shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                )
            );

            //insert value to map
            mIntegrationMethods.insert(PairType(Key, pNewGeometryData));

            std::cout << "Registered BezierGeometryData " << Key << " successfully" << std::endl;
        }
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static void RegisterIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2
    )
    {
        //define the key
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order1, Order2, 0, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);

        //find the key in existing map
        if(mIntegrationMethods.find(Key) != mIntegrationMethods.end())
        //found the key
        {
//            std::cout << "Key " << Key << " existed and ready to be used" << std::endl;
        }
        else
        {
            IntegrationPointsContainerType all_integration_points
                = AllIntegrationPoints(NumberOfIntegrationMethod, Order1, Order2);

            ShapeFunctionsValuesContainerType shape_functions_values;
            ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

            for (IndexType i = 0; i < NumberOfIntegrationMethod; ++i)
            {
                CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                    Order1,
                    Order2,
                    shape_functions_values[i],
                    shape_functions_local_gradients[i],
                    all_integration_points[i]
                );
            }

            GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
                new GeometryData(
                    TDimension,
                    TWorkingSpaceDimension,
                    TLocalSpaceDimension,
                    GeometryData::GI_GAUSS_2,           //ThisDefaultMethod
                    all_integration_points,             //ThisIntegrationPoints
                    shape_functions_values,             //ThisShapeFunctionsValues
                    shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                )
            );

            //insert value to map
            mIntegrationMethods.insert(PairType(Key, pNewGeometryData));

            std::cout << "Registered BezierGeometryData " << Key << " successfully" << std::endl;
        }
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static void RegisterIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2,
        unsigned int Order3
    )
    {
        //define the key
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order1, Order2, Order3, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);

        //find the key in existing map
        if(mIntegrationMethods.find(Key) != mIntegrationMethods.end())
        //found the key
        {
//            std::cout << "Key " << Key << " existed and ready to be used" << std::endl;
            return;
        }
        else
        {
            IntegrationPointsContainerType all_integration_points
                = AllIntegrationPoints(NumberOfIntegrationMethod, Order1, Order2, Order3);

            ShapeFunctionsValuesContainerType shape_functions_values;
            ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

            for (IndexType i = 0; i < NumberOfIntegrationMethod; ++i)
            {
                CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                    Order1,
                    Order2,
                    Order3,
                    shape_functions_values[i],
                    shape_functions_local_gradients[i],
                    all_integration_points[i]
                );
            }

            GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
                new GeometryData(
                    TDimension,
                    TWorkingSpaceDimension,
                    TLocalSpaceDimension,
                    GeometryData::GI_GAUSS_2,           //ThisDefaultMethod
                    all_integration_points,             //ThisIntegrationPoints
                    shape_functions_values,             //ThisShapeFunctionsValues
                    shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
                )
            );

            //insert value to map
            mIntegrationMethods.insert(PairType(Key, pNewGeometryData));

            std::cout << "Registered BezierGeometryData " << Key << " successfully" << std::endl;
        }
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer RetrieveIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1
    )
    {
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order1, 0, 0, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);
        return mIntegrationMethods[Key];
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer RetrieveIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2
    )
    {
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order1, Order2, 0, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);
        return mIntegrationMethods[Key];
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer RetrieveIntegrationRule(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2,
        unsigned int Order3
    )
    {
        BezierGeometryDataKey Key(NumberOfIntegrationMethod, Order1, Order2, Order3, TDimension, TWorkingSpaceDimension, TLocalSpaceDimension);
        return mIntegrationMethods[Key];
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer CreateIntegrationRule(
        const GeometryData::IntegrationMethod ThisIntegrationMethod,
        unsigned int Order1,
        const IntegrationPointsArrayType& integration_points
    )
    {
        IntegrationPointsContainerType all_integration_points;
        ShapeFunctionsValuesContainerType shape_functions_values;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        all_integration_points[ThisIntegrationMethod] = integration_points;

        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Order1,
            shape_functions_values[ThisIntegrationMethod],
            shape_functions_local_gradients[ThisIntegrationMethod],
            all_integration_points[ThisIntegrationMethod]
        );

        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
            new GeometryData(
                TDimension,
                TWorkingSpaceDimension,
                TLocalSpaceDimension,
                ThisIntegrationMethod,              //ThisDefaultMethod
                all_integration_points,             //ThisIntegrationPoints
                shape_functions_values,             //ThisShapeFunctionsValues
                shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
            )
        );

        std::cout << "Create BezierGeometryData successfully for " << integration_points.size() << " integration points" << std::endl;

        return pNewGeometryData;
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer CreateIntegrationRule(
        const GeometryData::IntegrationMethod ThisIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2,
        const IntegrationPointsArrayType& integration_points
    )
    {
        IntegrationPointsContainerType all_integration_points;
        ShapeFunctionsValuesContainerType shape_functions_values;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        all_integration_points[ThisIntegrationMethod] = integration_points;

        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Order1,
            Order2,
            shape_functions_values[ThisIntegrationMethod],
            shape_functions_local_gradients[ThisIntegrationMethod],
            all_integration_points[ThisIntegrationMethod]
        );

        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
            new GeometryData(
                TDimension,
                TWorkingSpaceDimension,
                TLocalSpaceDimension,
                ThisIntegrationMethod,              //ThisDefaultMethod
                all_integration_points,             //ThisIntegrationPoints
                shape_functions_values,             //ThisShapeFunctionsValues
                shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
            )
        );

        std::cout << "Create BezierGeometryData successfully for " << integration_points.size() << " integration points" << std::endl;

        return pNewGeometryData;
    }

    template<std::size_t TDimension, std::size_t TWorkingSpaceDimension, std::size_t TLocalSpaceDimension>
    static GeometryData::Pointer CreateIntegrationRule(
        const GeometryData::IntegrationMethod ThisIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2,
        unsigned int Order3,
        const IntegrationPointsArrayType& integration_points
    )
    {
        IntegrationPointsContainerType all_integration_points;
        ShapeFunctionsValuesContainerType shape_functions_values;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;

        all_integration_points[ThisIntegrationMethod] = integration_points;

        CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
            Order1,
            Order2,
            Order3,
            shape_functions_values[ThisIntegrationMethod],
            shape_functions_local_gradients[ThisIntegrationMethod],
            all_integration_points[ThisIntegrationMethod]
        );

        GeometryData::Pointer pNewGeometryData = GeometryData::Pointer(
            new GeometryData(
                TDimension,
                TWorkingSpaceDimension,
                TLocalSpaceDimension,
                ThisIntegrationMethod,              //ThisDefaultMethod
                all_integration_points,             //ThisIntegrationPoints
                shape_functions_values,             //ThisShapeFunctionsValues
                shape_functions_local_gradients     //ThisShapeFunctionsLocalGradients
            )
        );

        std::cout << "Create BezierGeometryData successfully for " << integration_points.size() << " integration points" << std::endl;

        return pNewGeometryData;
    }

    /********************************************************
            End of Bezier integration utilities
     ********************************************************/

    ///@}
    ///@name Access
    ///@{
    /********************************************************
            Geometry handling routines
     ********************************************************/
    /*Important remarks: this function only works correctly with Geo2dBezier and Geo3dBezier */
    template<class TElementType>
    static PointType& ComputeCentroid(typename TElementType::Pointer pElem, PointType& P)
    {
//        KRATOS_WATCH(typeid(*pElem).name())
//        GeometryType& Geom = *(pElem->pGetGeometry());
//        GeometryType& Geom = pElem->GetGeometry();
        IsogeometricGeometryType& Geom = dynamic_cast<IsogeometricGeometryType&>(pElem->GetGeometry());
//        KRATOS_WATCH(typeid(Geom).name())
        unsigned int dim = Geom.WorkingSpaceDimension();

        //number of integration points used
        IntegrationMethod ThisIntegrationMethod = Geom.GetDefaultIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints(ThisIntegrationMethod);

        //initializing the Jacobian, the inverse Jacobian and Jacobian determinant in the reference configuration
        GeometryType::JacobiansType J0(integration_points.size());

        //calculating the Jacobian
        J0 = Geom.Jacobian0(J0, ThisIntegrationMethod);

        //calculating contributions
        noalias(P) = ZeroVector(3);
        CoordinatesArrayType p_ref;
        CoordinatesArrayType p;
        Matrix InvJ0(dim, dim); // just a dummy variable
        double DetJ0;
        double TotalDomainInitialSize = 0.0;
        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix(J0[PointNumber], InvJ0, DetJ0);
            //calculating the total area
            double dS = DetJ0 * IntegrationWeight;
            TotalDomainInitialSize += dS;
            //calculating the global coordinates
            p_ref[0] = integration_points[PointNumber].X();
            p_ref[1] = integration_points[PointNumber].Y();
            p_ref[2] = integration_points[PointNumber].Z();
//            KRATOS_WATCH(p_ref)
            GlobalCoordinates(Geom, p, p_ref);
//            KRATOS_WATCH(p)
            P[0] += p[0] * dS;
            P[1] += p[1] * dS;
            P[2] += p[2] * dS;
        }
        P *= (1.0 / TotalDomainInitialSize);
        return P;
    }

    /*Important remarks: this function works with Geo2dBezier, Geo2dBezier3 and Geo3dBezier, but rather slow, don't know why */
    // TODO: to be profiled
//    template<class T>
//    static PointType& ComputeCentroid(typename T::Pointer& pElem, PointType& P)
//    {
//        IsogeometricGeometryType& Geom = dynamic_cast<IsogeometricGeometryType&>(pElem->GetGeometry());
//        unsigned int dim = Geom.WorkingSpaceDimension();
//
//        //number of integration points used
//        IntegrationMethod ThisIntegrationMethod = Geom.GetDefaultIntegrationMethod();
//        const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints(ThisIntegrationMethod);
//
//        //calculating determinants of Jacobian
//        Vector DetJ;
//        DetJ = Geom.DeterminantOfJacobian(DetJ, ThisIntegrationMethod);
//
//        //calculating contributions
//        noalias(P) = ZeroVector(3);
//        CoordinatesArrayType p_ref;
//        CoordinatesArrayType p;
//        double TotalDomainInitialSize = 0.0;
//        for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
//        {
//            //getting informations for integration
//            double IntegrationWeight = integration_points[PointNumber].Weight();
//
//            //calculating the area fraction
//            double dS = DetJ[PointNumber] * IntegrationWeight;
//            TotalDomainInitialSize += dS;
//
//            //getting the local coodinates of integration point
//            p_ref[0] = integration_points[PointNumber].X();
//            p_ref[1] = integration_points[PointNumber].Y();
//            p_ref[2] = integration_points[PointNumber].Z();
//
//            //calculating the global coordinates of integration point
//            GlobalCoordinates(Geom, p, p_ref);
//
//            //contributing to the global center of mass
//            P[0] += p[0] * dS;
//            P[1] += p[1] * dS;
//            P[2] += p[2] * dS;
//        }
// //        KRATOS_WATCH(TotalDomainInitialSize)
//        P *= (1.0 / TotalDomainInitialSize);
//        return P;
//    }

    /********************************************************
            End of Geometry handling routines
     ********************************************************/

    /********************************************************
            Dumping utilities
     ********************************************************/
    static void DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        ModelPart& rModelPart,
        std::string FileName
    )
    {
        #ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        #endif

        IntegrationMethod ThisMethod = GeometryData::GI_GAUSS_1;

        ElementsArrayType& pElements = rModelPart.Elements();

        std::ofstream Out;
        Out.open(FileName.c_str());

        Out << "{" << std::endl;
        for(typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin();
            it != pElements.ptr_end();
            ++it)
        {
            MatrixType shape_functions_values;
            ShapeFunctionsGradientsType shape_functions_local_gradients;

            IsogeometricGeometryType& rGeometry = dynamic_cast<IsogeometricGeometryType&>((*it)->GetGeometry());

            rGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                shape_functions_values,
                shape_functions_local_gradients,
                ThisMethod
            );

            unsigned int integration_points_number = shape_functions_values.size1();
            unsigned int number_of_nodes = shape_functions_values.size2();


            Out << "{" << std::endl;
            for(unsigned int i = 0; i < integration_points_number; ++i)
            {
                Out << "{";
                for(unsigned int j = 0; j < number_of_nodes - 1; ++j)
                {
                    Out << shape_functions_values(i, j) << ", ";
                }
                if(i != integration_points_number - 1)
                {
                    Out << shape_functions_values(i, number_of_nodes - 1) << "}," << std::endl;
                }
                else
                    Out << shape_functions_values(i, number_of_nodes - 1) << "}" << std::endl;;
            }

            if(it != pElements.ptr_end() - 1)
            {
                Out << "}," << std::endl;
            }
            else
                Out << "}" << std::endl;
        }
        Out << "};" << std::endl;
        Out.close();

        #ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Dumping shape functions values and local gradients completed: "
                  << (end_compute - start_compute) << " s" << std::endl;
        #endif
    }

    /********************************************************
            End of Dumping utilities
     ********************************************************/

    /********************************************************
            INTEGRATION POINT INTERFACE
     ********************************************************/

    static IntegrationPointsContainerType AllIntegrationPoints(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order
    )
    {
        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

        //define the base integration rule for 1st parametric dimension. This ensures that the number of integration points in each direction equal to order + 1
        std::vector<IntegrationPointsArrayType> BaseRule;
        BaseRule = GenerateBaseIntegrationRule(BaseRule);

        ///////////////////////////////////////////////////////////////
        // Remarks: this current implementation supports integration with order up to 9
        IndexType k, j1, offset1;
        IndexType base_offset = 1 + Order / 2;

        IntegrationPointsContainerType integration_points;
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset1 = k + base_offset;

            if(offset1 >= BaseRule.size())
                KRATOS_THROW_ERROR(std::logic_error, "There are not enough Gauss point to support for integration", __FUNCTION__)

            for(j1 = 0; j1 < BaseRule[offset1].size(); ++j1)
            {
                const IntegrationPointType& temp1 = BaseRule[offset1][j1];

                IntegrationPointType temp;

                temp.X() = 0.5 * (temp1.X() + 1);
                temp.Weight() = 0.5 * temp1.Weight();

                integration_points[k].push_back(temp);
            }
        }

        return integration_points;
    }


    static IntegrationPointsContainerType AllIntegrationPoints(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2
    )
    {
        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

        //define the base integration rule for 1st parametric dimension. This ensures that the number of integration points in each direction equal to order + 1
        std::vector<IntegrationPointsArrayType> BaseRule;
        BaseRule = GenerateBaseIntegrationRule(BaseRule);

        ///////////////////////////////////////////////////////////////
        // Remarks: this current implementation supports integration with order up to 9
        IndexType k, j1, j2, offset1, offset2;
        #ifdef USE_EQUAL_ORDER_INTEGRATION_IN_ALL_DIRECTION
        unsigned int IntOrder = std::max(Order1, Order2);
        IndexType base_offset1 = 1 + IntOrder / 2;
        IndexType base_offset2 = 1 + IntOrder / 2;
        #else
        IndexType base_offset1 = 1 + Order1 / 2;
        IndexType base_offset2 = 1 + Order2 / 2;
        #endif

        IntegrationPointsContainerType integration_points;
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset1 = k + base_offset1;
            offset2 = k + base_offset2;

            if(offset1 >= BaseRule.size() || offset2 >= BaseRule.size())
                KRATOS_THROW_ERROR(std::logic_error, "There are not enough Gauss point to support for integration", __FUNCTION__)

            for(j1 = 0; j1 < BaseRule[offset1].size(); ++j1)
            {
                const IntegrationPointType& temp1 = BaseRule[offset1][j1];

                for(j2 = 0; j2 < BaseRule[offset2].size(); ++j2)
                {
                    const IntegrationPointType& temp2 = BaseRule[offset2][j2];

                    IntegrationPointType temp;

                    temp.X() = 0.5 * (temp1.X() + 1);
                    temp.Y() = 0.5 * (temp2.X() + 1);
                    temp.Weight() = 0.25 * temp1.Weight() * temp2.Weight();

                    integration_points[k].push_back(temp);
                }
            }
        }

        return integration_points;
    }


    static IntegrationPointsContainerType AllIntegrationPoints(
        unsigned int NumberOfIntegrationMethod,
        unsigned int Order1,
        unsigned int Order2,
        unsigned int Order3
    )
    {
        //generate integration points for GI_GAUSS_1 rule. Note that GI_GAUSS_1 is the default rule which take the minimum order of integration in each direction
        std::vector<IntegrationPointsArrayType> GaussRule;

        //define the base integration rule for 1st parametric dimension. This ensures that the number of integration points in each direction equal to order + 1
        std::vector<IntegrationPointsArrayType> BaseRule;
        BaseRule = GenerateBaseIntegrationRule(BaseRule);

        ///////////////////////////////////////////////////////////////
        // Remarks: this current implementation supports integration with order up to 9
        IndexType k, j1, j2, j3, offset1, offset2, offset3;
        #ifdef USE_EQUAL_ORDER_INTEGRATION_IN_ALL_DIRECTION
        unsigned int IntOrder = std::max(std::max(Order1, Order2), Order3);
        IndexType base_offset1 = 1 + IntOrder / 2;
        IndexType base_offset2 = 1 + IntOrder / 2;
        IndexType base_offset3 = 1 + IntOrder / 2;
        #else
//        IndexType base_offset1 = 1 + Order1 / 2;
//        IndexType base_offset2 = 1 + Order2 / 2;
//        IndexType base_offset3 = 1 + Order3 / 2;
        IndexType base_offset1 = Order1 / 2;
        IndexType base_offset2 = Order2 / 2;
        IndexType base_offset3 = Order3 / 2;
//        IndexType base_offset1 = (Order1 / 2 >= 1) ? (Order1 / 2 - 1) : 0;
//        IndexType base_offset2 = (Order2 / 2 >= 1) ? (Order2 / 2 - 1) : 0;
//        IndexType base_offset3 = (Order3 / 2 >= 1) ? (Order3 / 2 - 1) : 0;
        #endif

        IntegrationPointsContainerType integration_points;
        for (k = 0; k < NumberOfIntegrationMethod; ++k)
        {
            offset1 = k + base_offset1; // this allows for one more order than the order of the curve
            offset2 = k + base_offset2; // this allows for one more order than the order of the curve
            offset3 = k + base_offset3; // this allows for one more order than the order of the curve

            if( offset1 >= BaseRule.size() ||
                offset2 >= BaseRule.size() ||
                offset3 >= BaseRule.size() )
                KRATOS_THROW_ERROR(std::logic_error, "There are not enough Gauss point to support for integration", __FUNCTION__)

            for(j1 = 0; j1 < BaseRule[offset1].size(); ++j1)
            {
                const IntegrationPointType& temp1 = BaseRule[offset1][j1];

                for(j2 = 0; j2 < BaseRule[offset2].size(); ++j2)
                {
                    const IntegrationPointType& temp2 = BaseRule[offset2][j2];

                    for(j3 = 0; j3 < BaseRule[offset3].size(); ++j3)
                    {
                        const IntegrationPointType& temp3 = BaseRule[offset3][j3];

                        IntegrationPointType temp;

                        temp.X() = 0.5 * (temp1.X() + 1);
                        temp.Y() = 0.5 * (temp2.X() + 1);
                        temp.Z() = 0.5 * (temp3.X() + 1);
                        temp.Weight() = 0.125 * temp1.Weight() * temp2.Weight() * temp3.Weight();

                        integration_points[k].push_back(temp);
                    }
                }
            }

//            KRATOS_WATCH(base_offset1)
//            KRATOS_WATCH(base_offset2)
//            KRATOS_WATCH(base_offset3)

//            KRATOS_WATCH(offset1)
//            KRATOS_WATCH(offset2)
//            KRATOS_WATCH(offset3)

//            KRATOS_WATCH(BaseRule[offset1].size())
//            KRATOS_WATCH(BaseRule[offset2].size())
//            KRATOS_WATCH(BaseRule[offset3].size())

            std::cout << BaseRule[offset1].size() * BaseRule[offset2].size() * BaseRule[offset3].size() << " integration points are generated" << std::endl;
        }
        return integration_points;
    }


    /********************************************************
            END OF INTEGRATION POINT INTERFACE
     ********************************************************/

    /********************************************************
            Bezier extraction subroutines
     ********************************************************/

    /**
    * Compute Bezier extraction for NURBS in 1D
    */
    template<class TValuesContainerType>
    static void bezier_extraction_1d(
        std::vector<Matrix>& C,
        int& nb, // number of elements
        const TValuesContainerType& U,
        const int p)
    {
        int m = U.size() - p - 1;
        int a = p + 1;
        int b = a + 1;
        nb = 1;
        int i, j, k, l, r, s, save, multiplicity;
        double alpha, numerator;
        std::vector<double> alphas(p+1);

        C.push_back(IdentityMatrix(p+1));

        while (b <= m)
        {
            if (nb+1 > C.size())
                C.push_back(IdentityMatrix(p+1));
            else
                noalias(C[nb]) = IdentityMatrix(p+1);
            i = b;
            while (b<=m && U[b]==U[b-1]) ++b;

            multiplicity = b - i + 1;
            if (multiplicity < p)
            {
                numerator = U[b-1] - U[a-1];
                for (j = p; j>=multiplicity+1; --j) alphas[j-multiplicity-1] = numerator/(U[a+j-1]-U[a-1]);
                r = p - multiplicity;
                for (j = 1; j <= r; ++j)
                {
                    save = r-j+1;
                    s = multiplicity + j;
                    for (k = p+1; k >= s+1; --k)
                    {
                        alpha = alphas[k-s-1];
                        for (l = 1; l <= p+1; ++l)
                            C[nb-1](l-1,k-1) = alpha*C[nb-1](l-1,k-1) + (1-alpha)*C[nb-1](l-1,k-2);
                    }
                    if (b <= m)
                    {
                        for (l = 0; l <= j; ++l)
                            C[nb](save+l-1, save-1) = C[nb-1](p-j+l, p);
                    }
                }
                nb = nb+1;
                if (b <= m)
                {
                    a = b;
                    b = b+1;
                }
            }
            else if (multiplicity == p)
            {
                if (b <= m)
                {
                    nb = nb+1;
                    a = b;
                    b = b+1;
                }
            }
        }

        assert(nb == C.size());
    }

    /**
    * Compute Bezier extraction for NURBS in 2D
    */
    template<class TValuesContainerType>
    static void bezier_extraction_2d(
        std::vector<Matrix>& C,
        int& nb1, // number of elements in u-direction
        int& nb2, // number of elements in v-direction
        const TValuesContainerType& U,
        const TValuesContainerType& V,
        const int p,
        const int q)
    {
        std::vector<Matrix> Cxi, Cet;

        BezierUtils::bezier_extraction_1d(Cxi, nb1, U, p);
        BezierUtils::bezier_extraction_1d(Cet, nb2, V, q);

        C.resize(nb1*nb2);

        int eta, xi, e, row, col, ird, jrd, icd, jcd, i, j;
        for (eta = 0; eta < nb2; ++eta)
        {
            for (xi = 0; xi < nb1; ++xi)
            {
                e = eta*nb1 + xi;
                C[e].resize((p+1)*(q+1), (p+1)*(q+1), false);
                for (row = 0; row < q+1; ++row)
                {
                    ird = row*(p+1);
                    jrd = (row+1)*(p+1)-1;
                    for (col = 0; col < q+1; ++col)
                    {
                        icd = col*(p+1);
                        jcd = (col+1)*(p+1)-1;
                        for (i = 0; i < p+1; ++i)
                        {
                            for (j = 0; j < p+1; ++j)
                            {
                                C[e](ird+i, icd+j) = Cet[eta](row, col) * Cxi[xi](i, j);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
    * Compute Bezier extraction for NURBS in 3D
    */
    template<class TValuesContainerType>
    static void bezier_extraction_3d(
        std::vector<Matrix>& C,
        int& nb1, // number of elements in u-direction
        int& nb2, // number of elements in v-direction
        int& nb3, // number of elements in w-direction
        const TValuesContainerType& U,
        const TValuesContainerType& V,
        const TValuesContainerType& W,
        const int p,
        const int q,
        const int r)
    {
        std::vector<Matrix> Cxi, Cet, Cze;
        Matrix C_et_xi;

        bezier_extraction_1d(Cxi, nb1, U, p);
        bezier_extraction_1d(Cet, nb2, V, q);
        bezier_extraction_1d(Cze, nb3, W, r);

        C_et_xi.resize((p+1)*(q+1), (p+1)*(q+1), false);
        C.resize(nb1*nb2*nb3);

        int eta, xi,zeta, e, row, col, ird, jrd, icd, jcd, i, j;
        for (eta = 0; eta < nb2; ++eta)
        {
            for (xi = 0; xi < nb1; ++xi)
            {
                for (row = 0; row < q+1; ++row)
                {
                    ird = row*(p+1);
                    jrd = (row+1)*(p+1)-1;
                    for (col = 0; col < q+1; ++col)
                    {
                        icd = col*(p+1);
                        jcd = (col+1)*(p+1)-1;
                        for (i = 0; i < p+1; ++i)
                        {
                            for (j = 0; j < p+1; ++j)
                            {
                                C_et_xi(ird+i, icd+j) = Cet[eta](row, col) * Cxi[xi](i, j);
                            }
                        }
                    }
                }

                for (zeta = 0; zeta < nb3; ++zeta)
                {
                    e = (zeta * nb2 + eta) * nb1 + xi;
                    C[e].resize((p+1)*(q+1)*(r+1), (p+1)*(q+1)*(r+1));
                    for (row = 0; row < r+1; ++row)
                    {
                        ird = row*(p+1)*(q+1);
                        jrd = (row+1)*(p+1)*(q+1)-1;
                        for (col = 0; col < r+1; ++col)
                        {
                            icd = col*(p+1)*(q+1);
                            jcd = (col+1)*(p+1)*(q+1)-1;
                            for (i = 0; i < (p+1)*(q+1); ++i)
                            {
                                for (j = 0; j < (p+1)*(q+1); ++j)
                                {
                                    C[e](ird+i, icd+j) = Cze[zeta](row, col) * C_et_xi(i, j);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
        Compute extended knot vector given the local knot vector
     */
//    static void compute_extended_knot_vector(
//        Vector& Ubar,
//        int& nt,
//        const std::vector<double>& Xi,
//        const int p);

    /**
        Compute the Bezier extraction of T-splines basis function on knot spans
        This is 1D version of the code
     */
    static void bezier_extraction_tsplines_1d(
        std::vector<Vector>& Crows,
        int& nb,
        Vector& Ubar,
        const std::vector<double>& Xi,
        const std::vector<double>& U,
        const std::vector<int>& spans,
        const int p);

    /**
        Compute the Bezier extraction of T-splines basis function on knot spans
        This is a not efficient but more stable version of the bezier_extraction_tsplines_1d
     */
    static void bezier_extraction_local_1d(
        std::vector<Vector>& Crows,
        int& nb,
        Vector& Ubar,
        const std::vector<double>& Xi,
        const std::vector<double>& U,
        const int p);

    /**
        Compute the Bezier extraction of T-splines basis function on knot spans
        This is 2D version of the code
     */
    static void bezier_extraction_tsplines_2d(
        std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const std::vector<int>& spans_xi,
        const std::vector<int>& spans_eta,
        const int p,
        const int q);

    /**
        Compute the Bezier extraction of T-splines basis function on knot spans
        This is a not efficient but more stable version of the bezier_extraction_tsplines_2d
     */
    static void bezier_extraction_local_2d(
        std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const int p,
        const int q);

    /**
        Compute the Bezier extraction of T-splines/NURBS basis function on knot spans
     */
    static void bezier_extraction_local_3d(
        std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        int& nb_zeta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        Vector& Ubar_zeta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Zeta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const std::vector<double>& Uzeta,
        const int p,
        const int q,
        const int r);

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "BezierUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BezierUtils";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
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

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const int msBernsteinCoefs[];

    static MapType mIntegrationMethods;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Calculate global coodinates w.r.t initial configuration
     */
    static CoordinatesArrayType& GlobalCoordinates(
        GeometryType& rGeometry,
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates
    )
    {
        noalias( rResult ) = ZeroVector( 3 );

        Vector ShapesFunctionValues;

        rGeometry.ShapeFunctionsValues(ShapesFunctionValues, LocalCoordinates);

        for ( IndexType i = 0 ; i < rGeometry.size() ; ++i )
        {
            noalias( rResult ) += ShapesFunctionValues( i ) * rGeometry.GetPoint( i ).GetInitialPosition();
        }

        return rResult;
    }

    /**************************************************************************
                1st order NURBS
     **************************************************************************/
    static void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        IndexType Order,
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        const IntegrationPointsArrayType& integration_points
    )
    {
        shape_functions_local_gradients.resize(integration_points.size());
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType());
        shape_functions_values.resize(integration_points.size(), Order + 1);

        for (unsigned int it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            VectorType temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                Order,
                temp_values,
                shape_functions_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for(unsigned int i = 0; i < Order + 1; ++i)
            {
                shape_functions_values( it_gp, i ) = temp_values(i);
            }
        }
    }

    static void ShapeFunctionsValuesAndLocalGradients(
        IndexType Order,
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    )
    {
        //compute all Bezier shape functions & derivatives at rPoint
        shape_functions_values.resize(Order + 1);
        VectorType shape_functions_derivatives(Order + 1);
        bernstein(shape_functions_values, shape_functions_derivatives, Order, rPoint[0]);

        shape_functions_local_gradients.resize(1, Order + 1);
        for(int i = 0; i < Order + 1; ++i)
            shape_functions_local_gradients(0, i) = shape_functions_derivatives(i);
    }

    /**************************************************************************
                2nd order NURBS
     **************************************************************************/
    static void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        IndexType Order1,
        IndexType Order2,
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        const IntegrationPointsArrayType& integration_points
    )
    {
        shape_functions_local_gradients.resize(integration_points.size());
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType());
        shape_functions_values.resize(integration_points.size(), (Order1 + 1) * (Order2 + 1));

        for (unsigned int it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            VectorType temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                Order1,
                Order2,
                temp_values,
                shape_functions_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for(unsigned int i = 0; i < (Order1 + 1) * (Order2 + 1); ++i)
            {
                shape_functions_values( it_gp, i ) = temp_values(i);
            }
        }
    }

    static void ShapeFunctionsValuesAndLocalGradients(
        IndexType Order1,
        IndexType Order2,
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    )
    {
        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(Order1 + 1);
        VectorType bezier_functions_values2(Order2 + 1);
        VectorType bezier_functions_derivatives1(Order1 + 1);
        VectorType bezier_functions_derivatives2(Order2 + 1);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, Order1, rPoint[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, Order2, rPoint[1]);

        //compute bivariate Bezier shape functions values
        shape_functions_values.resize((Order1 + 1) * (Order2 + 1));
        for(IndexType i = 0; i < Order1 + 1; ++i)
        {
            for(IndexType j = 0; j < Order2 + 1; ++j)
            {
                IndexType index = j + i * (Order2 + 1);

                shape_functions_values(index) =
                    bezier_functions_values1(i) * bezier_functions_values2(j);
            }
        }

        //compute bivariate Bezier shape functions derivatives w.r.t local coordinates
        shape_functions_local_gradients.resize(2, (Order1 + 1) * (Order2 + 1));
        for(IndexType i = 0; i < Order1 + 1; ++i)
        {
            for(IndexType j = 0; j < Order2 + 1; ++j)
            {
                IndexType index = j + i * (Order2 + 1);

                shape_functions_local_gradients(0, index) =
                    bezier_functions_derivatives1(i) * bezier_functions_values2(j);
                shape_functions_local_gradients(1, index) =
                    bezier_functions_values1(i) * bezier_functions_derivatives2(j);
            }
        }
    }

    /**************************************************************************
                3rd order NURBS
     **************************************************************************/
    static void CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
        IndexType Order1,
        IndexType Order2,
        IndexType Order3,
        MatrixType& shape_functions_values,
        ShapeFunctionsGradientsType& shape_functions_local_gradients,
        const IntegrationPointsArrayType& integration_points
    )
    {
        shape_functions_local_gradients.resize(integration_points.size());
        std::fill(shape_functions_local_gradients.begin(), shape_functions_local_gradients.end(), MatrixType());
        shape_functions_values.resize(integration_points.size(), (Order1 + 1) * (Order2 + 1) * (Order3 + 1));

        for (unsigned int it_gp = 0; it_gp < integration_points.size(); ++it_gp)
        {
            VectorType temp_values;
            ShapeFunctionsValuesAndLocalGradients(
                Order1,
                Order2,
                Order3,
                temp_values,
                shape_functions_local_gradients[it_gp],
                integration_points[it_gp]
            );
            for(unsigned int i = 0; i < (Order1 + 1) * (Order2 + 1) * (Order3 + 1); ++i)
            {
                shape_functions_values( it_gp, i ) = temp_values(i);
            }
        }
    }

    static void ShapeFunctionsValuesAndLocalGradients(
        IndexType Order1,
        IndexType Order2,
        IndexType Order3,
        VectorType& shape_functions_values,
        MatrixType& shape_functions_local_gradients,
        const CoordinatesArrayType& rPoint
    )
    {
        //compute all univariate Bezier shape functions & derivatives at rPoint
        VectorType bezier_functions_values1(Order1 + 1);
        VectorType bezier_functions_values2(Order2 + 1);
        VectorType bezier_functions_values3(Order3 + 1);
        VectorType bezier_functions_derivatives1(Order1 + 1);
        VectorType bezier_functions_derivatives2(Order2 + 1);
        VectorType bezier_functions_derivatives3(Order3 + 1);
        BezierUtils::bernstein(bezier_functions_values1, bezier_functions_derivatives1, Order1, rPoint[0]);
        BezierUtils::bernstein(bezier_functions_values2, bezier_functions_derivatives2, Order2, rPoint[1]);
        BezierUtils::bernstein(bezier_functions_values3, bezier_functions_derivatives3, Order3, rPoint[2]);

        //compute trivariate Bezier shape functions values
        shape_functions_values.resize((Order1 + 1) * (Order2 + 1) * (Order3 + 1));
        for(IndexType i = 0; i < (Order1 + 1); ++i)
        {
            for(IndexType j = 0; j < (Order2 + 1); ++j)
            {
                for(IndexType k = 0; k < (Order3 + 1); ++k)
                {
                    IndexType index = k + (j + i * (Order2 + 1)) * (Order3 + 1);

                    shape_functions_values(index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                }
            }
        }

        //compute trivariate Bezier shape functions derivatives w.r.t local coordinates
        shape_functions_local_gradients.resize(3, (Order1 + 1) * (Order2 + 1) * (Order3 + 1));
        for(IndexType i = 0; i < Order1 + 1; ++i)
        {
            for(IndexType j = 0; j < Order2 + 1; ++j)
            {
                for(IndexType k = 0; k < Order3 + 1; ++k)
                {
                    IndexType index = k + (j + i * (Order2 + 1)) * (Order3 + 1);

                    shape_functions_local_gradients(0, index) =
                        bezier_functions_derivatives1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_values3(k);
                    shape_functions_local_gradients(1, index) =
                        bezier_functions_values1(i) *
                        bezier_functions_derivatives2(j) *
                        bezier_functions_values3(k);
                    shape_functions_local_gradients(2, index) =
                        bezier_functions_values1(i) *
                        bezier_functions_values2(j) *
                        bezier_functions_derivatives3(k);
                }
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    static std::vector<IntegrationPointsArrayType>& GenerateBaseIntegrationRule(
        std::vector<IntegrationPointsArrayType>& rRule
    )
    {
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints6, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints7, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints8, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints9, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());
        rRule.push_back(Quadrature<LineGaussLegendreIntegrationPoints10, 1, IntegrationPoint<3> >::GenerateIntegrationPoints());

        return rRule;
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BezierUtils& operator=(BezierUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BezierUtils(BezierUtils const& rOther)
    {
    }

    ///@}

}; // Class BezierUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BezierUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const BezierUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef ENABLE_PROFILING
#undef USE_EQUAL_ORDER_INTEGRATION_IN_ALL_DIRECTION

#endif // KRATOS_BEZIER_UTILS_H_INCLUDED  defined

