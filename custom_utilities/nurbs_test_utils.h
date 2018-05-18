//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013 Sep 12 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_NURBS_TEST_UTILS_H_INCLUDED )
#define  KRATOS_NURBS_TEST_UTILS_H_INCLUDED

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
class NURBSTestUtils
{
public:
    ///@name Type Definitions
    ///@{

//      typedef std::vector<double>       ValueContainerType;
    typedef boost::numeric::ublas::vector<double> ValueContainerType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename Element::GeometryType::PointType::CoordinatesArrayType CoordinatesArrayType;

    /// Pointer definition of NURBSTestUtils
    KRATOS_CLASS_POINTER_DEFINITION(NURBSTestUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NURBSTestUtils()
    {}

    /// Destructor.
    virtual ~NURBSTestUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Test1(ModelPart& r_model_part, int NumTestPoints)
    {
        ElementsArrayType& pElements = r_model_part.Elements();

        std::vector<CoordinatesArrayType> Points;
        for(unsigned int i = 0; i < NumTestPoints + 1; ++i)
        {
            CoordinatesArrayType Point;
            Point(0) = (double) i / NumTestPoints;
            Point(1) = 0.0;
            Point(2) = 0.0;
            Points.push_back(Point);
        }

        int NumberOfNodes = (*(pElements.ptr_begin()))->GetGeometry().size();

        std::cout << "Inspecting shape function values at points" << std::endl;
        Matrix Results(NumberOfNodes, Points.size());
        for(unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            std::cout << "N[" << i << "](...) = ( ";
            for (unsigned int j = 0; j < Points.size(); j++)
            {
                Results(i, j) = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionValue(i, Points[j]);
                std::cout << Results(i, j) << " ";
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
        temp = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionsLocalGradients(temp, Probe);

        CoordinatesArrayType Probe1;
        Probe1(0) = 0.5 + 1e-6;
        Probe1(1) = 0.0;
        Probe1(2) = 0.0;

        int node = 2;
        std::cout << "Inspecting derivative at node 2:" << std::endl;
        double v1 = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionValue(node, Probe);
        double v2 = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionValue(node, Probe1);
        std::cout << "exact derivative = " << temp(2, 0) << std::endl;
        std::cout << "approximate derivative = " << (v2 - v1) / 1e-6 << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting shape function derivatives at points" << std::endl;
        for(unsigned int i = 0; i < Points.size(); ++i)
        {
            std::cout << "dN(xi[" << i << "])(...) = ";
            Matrix temp;
            temp = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionsLocalGradients(temp, Points[i]);
            std::cout << temp << std::endl;
        }
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Calculating global coordinates along the curve" << std::endl;
        std::ofstream LogFile;
        LogFile.open("curve_data.txt");
        for (unsigned int i = 0; i < Points.size(); ++i)
        {
            CoordinatesArrayType p;
            p = (*(pElements.ptr_begin()))->GetGeometry().GlobalCoordinates(p, Points[i]);

            std::cout << "Global coordinate at " << Points[i] << ": " << p << std::endl;
            LogFile << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        LogFile.close();
        std::cout << "------------------------------------------" << std::endl;

        // profiling basisfuns
//        std::cout << "Calculating global coordinates along the curve" << std::endl;
//        double start = OpenMPUtils::GetCurrentTime();
//        for (unsigned int i = 0; i < Points.size(); ++i)
//        {
//            CoordinatesArrayType p;
//            p = (*(pElements.ptr_begin()))->GetGeometry().GlobalCoordinates(p, Points[i]);
//        }
//        std::cout << "Calculating time = "<< OpenMPUtils::GetCurrentTime() - start << std::endl;

    }

    void Test2(ModelPart& r_model_part)
    {
        ElementsArrayType& pElements = r_model_part.Elements();

        std::cout << "Inspecting all integration points:" << std::endl;
        const GeometryType::IntegrationPointsArrayType& integration_points =
        (*(pElements.ptr_begin()))->GetGeometry().IntegrationPoints( (GeometryData::IntegrationMethod)0 );

        for (unsigned int i = 0; i < integration_points.size(); ++i)
        {
            KRATOS_WATCH(integration_points[i])
        }
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting all shape function values at all integration points:" << std::endl;
        const Matrix& Ncontainer = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionsValues( (GeometryData::IntegrationMethod)0 );
        KRATOS_WATCH(Ncontainer)
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Inspecting all shape function local gradients at all integration points:" << std::endl;
        const GeometryType::ShapeFunctionsGradientsType& DN_De = (*(pElements.ptr_begin()))->GetGeometry().ShapeFunctionsLocalGradients( (GeometryData::IntegrationMethod)0 );
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
        buffer << "NURBSTestUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NURBSTestUtils";
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
    NURBSTestUtils& operator=(NURBSTestUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    NURBSTestUtils(NURBSTestUtils const& rOther)
    {
    }

    ///@}

}; // Class NURBSTestUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, NURBSTestUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const NURBSTestUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_BSPLINE_UTILS_H_INCLUDED  defined
