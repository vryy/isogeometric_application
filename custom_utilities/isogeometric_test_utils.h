//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013 Sep 12 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_TEST_UTILS_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_TEST_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

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

/// Short class definition.
/** Detail class definition.
 */
class IsogeometricTestUtils
{
public:
    ///@name Type Definitions
    ///@{

//      typedef std::vector<double>       ValueContainerType;
    typedef boost::numeric::ublas::vector<double> ValueContainerType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename Element::GeometryType::PointType::CoordinatesArrayType CoordinatesArrayType;

    /// Pointer definition of IsogeometricTestUtils
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricTestUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricTestUtils()
    {}

    /// Destructor.
    virtual ~IsogeometricTestUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Test1(ModelPart& r_model_part, std::size_t NumTestPoints)
    {
        std::cout << std::setprecision(10);

        ElementsContainerType& pElements = r_model_part.Elements();
        Element::Pointer pElement = (*(pElements.ptr_begin()));

        std::vector<CoordinatesArrayType> Points;
        for(unsigned int i = 0; i < NumTestPoints + 1; ++i)
        {
            CoordinatesArrayType Point;
            Point(0) = (double) i / NumTestPoints;
            Point(1) = 0.0;
            Point(2) = 0.0;
            Points.push_back(Point);
        }
        KRATOS_WATCH(Points.size())

        std::cout << "List of test points:" << std::endl;
        for (unsigned int j = 0; j < Points.size(); j++)
            std::cout << "xi[" << j  << "] = " << Points[j] << std::endl;

        int NumberOfNodes = pElement->GetGeometry().size();

        std::cout << "Inspecting shape function values at test points" << std::endl;
        Matrix Results(NumberOfNodes, Points.size());
        for(unsigned int i = 0; i < Points.size(); ++i)
        {
            std::cout << "N(xi[" << i << "])(...) = (";
            for(unsigned int j = 0; j < NumberOfNodes; ++j)
            {
                double temp = pElement->GetGeometry().ShapeFunctionValue(j, Points[i]);
                std::cout << temp << " ";
            }
            std::cout << ")" << std::endl;
        }
        std::cout << "------------------------------------------" << std::endl;

        //checking derivatives
        CoordinatesArrayType Probe;
        Probe(0) = 0.5;
        Probe(1) = 0.0;
        Probe(2) = 0.0;
        Matrix temp;
        temp = pElement->GetGeometry().ShapeFunctionsLocalGradients(temp, Probe);

        CoordinatesArrayType Probe1;
        Probe1(0) = 0.5 + 1e-6;
        Probe1(1) = 0.0;
        Probe1(2) = 0.0;

        int node = 2;
        std::cout << "Inspecting derivative at node #" << node << ":" << std::endl;
        double v1 = pElement->GetGeometry().ShapeFunctionValue(node, Probe);
        double v2 = pElement->GetGeometry().ShapeFunctionValue(node, Probe1);
        std::cout << "exact derivative = " << temp(2, 0) << std::endl;
        std::cout << "approximate derivative = " << (v2 - v1) / 1e-6 << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting shape function derivatives at test points" << std::endl;
        for(unsigned int i = 0; i < Points.size(); ++i)
        {
            std::cout << "dN(xi[" << i << "])(...) = ";
            Matrix temp;
            temp = pElement->GetGeometry().ShapeFunctionsLocalGradients(temp, Points[i]);
            std::cout << temp << std::endl;
        }
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting shape function second derivatives at test points" << std::endl;
        for(unsigned int i = 0; i < Points.size(); ++i)
        {
            std::cout << "d2N(xi[" << i << "])(...) = ";
            GeometryType::ShapeFunctionsSecondDerivativesType temp;
            temp = pElement->GetGeometry().ShapeFunctionsSecondDerivatives(temp, Points[i]);
            std::cout << temp << std::endl;
        }
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Calculating global coordinates along the curve" << std::endl;
        std::ofstream LogFile;
        LogFile.open("curve_data.txt");
        for (unsigned int i = 0; i < Points.size(); ++i)
        {
            CoordinatesArrayType p;
            p = pElement->GetGeometry().GlobalCoordinates(p, Points[i]);

            std::cout << "Global coordinate at " << Points[i] << ": " << p << std::endl;
            LogFile << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        LogFile.close();
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Calculating Jacobian along the curve" << std::endl;
        LogFile.open("curve_data.txt", std::ofstream::out | std::ofstream::app);
        Matrix N;
        for (unsigned int i = 0; i < Points.size(); ++i)
        {
            N = pElement->GetGeometry().ShapeFunctionsLocalGradients(temp, Points[i]);
            CoordinatesArrayType t;
            noalias(t) = ZeroVector(3);
            for(unsigned int j = 0; j < NumberOfNodes; ++j)
                noalias(t) += N(j, 0)*pElement->GetGeometry()[j].GetInitialPosition();

            std::cout << "Tangent at " << Points[i] << ": " << t << std::endl;
            LogFile << t[0] << " " << t[1] << " " << t[2] << std::endl;
        }
        LogFile.close();
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Calculating Hessian along the curve" << std::endl;
        LogFile.open("curve_data.txt", std::ofstream::out | std::ofstream::app);
        GeometryType::ShapeFunctionsSecondDerivativesType H;
        for (unsigned int i = 0; i < Points.size(); ++i)
        {
            H = pElement->GetGeometry().ShapeFunctionsSecondDerivatives(H, Points[i]);
            CoordinatesArrayType h;
            noalias(h) = ZeroVector(3);
            for(unsigned int j = 0; j < NumberOfNodes; ++j)
                noalias(h) += H[j](0, 0)*pElement->GetGeometry()[j].GetInitialPosition();

            std::cout << "Hessian at " << Points[i] << ": " << h << std::endl;
            LogFile << h[0] << " " << h[1] << " " << h[2] << std::endl;
        }
        LogFile.close();
        std::cout << "------------------------------------------" << std::endl;

        // profiling basisfuns
//        std::cout << "Calculating global coordinates along the curve" << std::endl;
//        double start = OpenMPUtils::GetCurrentTime();
//        for (unsigned int i = 0; i < Points.size(); ++i)
//        {
//            CoordinatesArrayType p;
//            p = pElement->GetGeometry().GlobalCoordinates(p, Points[i]);
//        }
//        std::cout << "Calculating time = "<< OpenMPUtils::GetCurrentTime() - start << std::endl;

    }

    void Test2(ModelPart& r_model_part)
    {
        ElementsContainerType& pElements = r_model_part.Elements();
        Element::Pointer pElement = (*(pElements.ptr_begin()));

        std::cout << "Inspecting all integration points:" << std::endl;
        const GeometryType::IntegrationPointsArrayType& integration_points =
        pElement->GetGeometry().IntegrationPoints( (GeometryData::IntegrationMethod)0 );

        for (unsigned int i = 0; i < integration_points.size(); ++i)
        {
            KRATOS_WATCH(integration_points[i])
        }
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting all shape function values at all integration points:" << std::endl;
        const Matrix& Ncontainer = pElement->GetGeometry().ShapeFunctionsValues( (GeometryData::IntegrationMethod)0 );
        KRATOS_WATCH(Ncontainer)
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting all shape function local gradients at all integration points:" << std::endl;
        const GeometryType::ShapeFunctionsGradientsType& DN_De = pElement->GetGeometry().ShapeFunctionsLocalGradients( (GeometryData::IntegrationMethod)0 );
        for (unsigned int i = 0; i < DN_De.size(); ++i)
        {
            KRATOS_WATCH(DN_De[i])
        }
        std::cout << "------------------------------------------" << std::endl;

    }

    void ProbeGlobalCoordinates(Element::Pointer& pElement, double X, double Y, double Z = 0.0)
    {
        CoordinatesArrayType p;
        p[0] = X;
        p[1] = Y;
        p[2] = Z;

        CoordinatesArrayType Result = ZeroVector( 3 );

        pElement->GetGeometry().GlobalCoordinates(Result, p);

        std::cout << "Global coordinates at " << p << ": " << Result << std::endl;
    }

    void ProbeShapeFunctionValues(Element::Pointer& pElement, double X, double Y, double Z = 0.0)
    {
        CoordinatesArrayType p;
        p[0] = X;
        p[1] = Y;
        p[2] = Z;

        Vector Result;
        Result = pElement->GetGeometry().ShapeFunctionsValues(Result, p);

        for(int i = 0; i < pElement->GetGeometry().size(); ++i)
        {
            std::cout << "Shape function value " << (i+1) << " at " << p << ": " << Result(i) << std::endl;
        }
    }

    void ProbeShapeFunctionDerivatives(Element::Pointer& pElement, double X, double Y, double Z = 0.0)
    {
        CoordinatesArrayType p;
        p[0] = X;
        p[1] = Y;
        p[2] = Z;

        Matrix Result;

        pElement->GetGeometry().ShapeFunctionsLocalGradients(Result, p);

        std::cout << "Shape function local gradients at " << p << ":\n" << Result << std::endl;
    }

    void ProbeJacobian(Element::Pointer& pElement, double X, double Y, double Z = 0.0)
    {
        CoordinatesArrayType p;
        p[0] = X;
        p[1] = Y;
        p[2] = Z;

        Matrix Result;

        pElement->GetGeometry().Jacobian(Result, p);

        std::cout << "Jacobian at " << p << ":\n" << Result << std::endl;
    }

    template<class TDataType>
    void DumpNodalValues(
        Variable<TDataType>& rVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        std::cout << "Dumping nodal results " << rVariable.Name() << ": " << std::endl;
        for(NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end();  ++it)
        {
            std::cout << "Node " << (*it)->Id() << ": " << (*it)->GetSolutionStepValue(rVariable) << std::endl;
        }
    }

    ///@}
    ///@name Access
    ///@{

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
        buffer << "IsogeometricTestUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricTestUtils";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

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

    ///@}
    ///@name Member Variables
    ///@{

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

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IsogeometricTestUtils& operator=(IsogeometricTestUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricTestUtils(IsogeometricTestUtils const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricTestUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, IsogeometricTestUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const IsogeometricTestUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_TEST_UTILS_H_INCLUDED  defined
