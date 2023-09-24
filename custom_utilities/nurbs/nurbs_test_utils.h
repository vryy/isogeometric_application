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

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"


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

template<int TDim>
class NURBSTestUtils;

template<>
class NURBSTestUtils<1>
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValueContainerType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename Element::GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;

    typedef typename GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

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
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function computes and the values, derivatives and second derivatives at the integration of the conditions
     * and then compares it with the values computed from the patch
     */
    static void ProbeAndTestValuesOnPatch(typename Patch<1>::Pointer pPatch, ConditionsContainerType& rConditions,
        const IntegrationMethod& ThisIntegrationMethod, double tol)
    {
        // get the FESpace
        typename BSplinesFESpace<1>::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpace<1> >(pPatch->pFESpace());
        if (pFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace<1> is failed.", "")

        // extract the control point grid function
        // auto& rControlGridFunction = pPatch->ControlPointGridFunction();
        auto pControlGridFunction = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

        Vector shape_values;
        Matrix shape_local_gradients;
        ShapeFunctionsSecondDerivativesType shape_second_derivatives;
        std::vector<double> xi(1);
        PointType C, T, dT;

        for (typename ConditionsContainerType::ptr_iterator it = rConditions.ptr_begin();
                it != rConditions.ptr_end(); ++it)
        {
            std::cout << "Probing condition " << (*it)->Id() << std::endl;

            // get the knot span
            double kleft = (*it)->GetValue(KNOT_LEFT);
            double kright = (*it)->GetValue(KNOT_RIGHT);
            KRATOS_WATCH(kleft)
            KRATOS_WATCH(kright)
            // obtain the integration points

            const GeometryType::IntegrationPointsArrayType& integration_points = (*it)->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
            for (auto point = integration_points.begin(); point != integration_points.end(); ++point)
            {
                KRATOS_WATCH(*point)

                // compute the shape function values from the condition
                shape_values = (*it)->GetGeometry().ShapeFunctionsValues(shape_values, *point);
                KRATOS_WATCH(shape_values)

                // compute the point
                noalias(C) = ZeroVector(3);
                for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
                    noalias(C) += shape_values[i] * (*it)->GetGeometry()[i].GetInitialPosition();
                KRATOS_WATCH(C)

                // compute the shape function local gradients from the condition
                shape_local_gradients = (*it)->GetGeometry().ShapeFunctionsLocalGradients(shape_local_gradients, *point);
                // shape_local_gradients *= 1.0 / (kright - kleft);
                KRATOS_WATCH(shape_local_gradients)

                // compute the tangent
                noalias(T) = ZeroVector(3);
                for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
                    noalias(T) += shape_local_gradients(i, 0) * (*it)->GetGeometry()[i].GetInitialPosition();
                KRATOS_WATCH(T)

                // compute the shape function second derivatives from the condition
                shape_second_derivatives = (*it)->GetGeometry().ShapeFunctionsSecondDerivatives(shape_second_derivatives, *point);
                KRATOS_WATCH(shape_second_derivatives)

                // compute the derivatives of tangent
                noalias(dT) = ZeroVector(3);
                for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
                    noalias(dT) += shape_second_derivatives[i](0, 0) * (*it)->GetGeometry()[i].GetInitialPosition();
                KRATOS_WATCH(dT)

                //////////////////////////////

                // compute the knot on the parameter space
                xi[0] = (1.0 - point->X())*kleft + point->X()*kright;
                KRATOS_WATCH(xi[0])

                // compute the point from the grid function
                auto Cref = pControlGridFunction->GetValue(xi);
                KRATOS_WATCH(Cref)
                std::vector<double> shape_values_vec;
                pFESpace->GetValues(shape_values_vec, xi);
                std::cout << "shape_values_vec:";
                for (std::size_t i = 0; i < shape_values_vec.size(); ++i)
                    std::cout << " " << shape_values_vec[i];
                std::cout << std::endl;

                // compute the tangent from the grid function
                auto Tref = pControlGridFunction->GetDerivative(xi);
                // std::cout << "Tref: " << Tref[0] << ", " << Tref[1] << std::endl;
                KRATOS_WATCH(Tref.size())
                KRATOS_WATCH(Tref[0])

                std::vector<std::vector<double> > shape_local_gradients_vec;
                pFESpace->GetDerivatives(shape_local_gradients_vec, xi);
                std::cout << "shape_local_gradients_vec:";
                for (std::size_t i = 0; i < shape_local_gradients_vec.size(); ++i)
                {
                    std::cout << " (";
                    for (std::size_t j = 0; j < shape_local_gradients_vec[i].size(); ++j)
                        std::cout << " " << shape_local_gradients_vec[i][j];
                    std::cout << ")";
                }
                std::cout << std::endl;

                //////////////////////////////

                double Cdiff = std::abs(C[0] - Cref[0]);
                if (Cdiff > tol)
                    KRATOS_THROW_ERROR(std::logic_error, "Error computing the point", "")

                // double Tdiff = std::abs(T[0] - Tref[0]);
                // if (Tdiff > tol)
                //     KRATOS_THROW_ERROR(std::logic_error, "Error computing the tangent", "")

                std::cout << "--------------" << std::endl;
            }

            std::cout << "------------------------------" << std::endl;
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
template<int TDim>
inline std::istream& operator >>(std::istream& rIStream, NURBSTestUtils<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSTestUtils<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_NURBS_TEST_UTILS_H_INCLUDED  defined
