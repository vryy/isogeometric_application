//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ISOGEOMETRIC_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ISOGEOMETRIC_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/multipatch.h"

namespace Kratos
{

/**
 * Abstract class for all isogeometric utility.
 */
class IsogeometricUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricUtility);

    /// Type definition

    /// Default constructor
    IsogeometricUtility() {}

    /// Destructor
    virtual ~IsogeometricUtility() {}

    /// Get the integration method from the integration order
    static GeometryData::IntegrationMethod GetIntegrationMethod(int integration_order)
    {
        if (integration_order == 1)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if (integration_order == 2)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if (integration_order == 3)
        {
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else
            KRATOS_ERROR << "Invalid integration order " << integration_order;
        }


    /// Find the element in the KRATOS container with specific key
    template<class TContainerType, class TKeyType>
    static typename TContainerType::const_iterator FindKey(const TContainerType& ThisContainer,
            TKeyType ThisKey, std::string ComponentName)
    {
        typename TContainerType::const_iterator i_result;
        if ((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            KRATOS_ERROR << ComponentName << " #" << ThisKey << " is not found.";
        }

        return i_result;
    }

    /// Create a condition taking the same geometry as the parent element
    static Condition::Pointer CreateConditionFromElement(const std::string& sample_condition_name,
            std::size_t& lastConditionId, Element::Pointer pElement, Properties::Pointer pProperties )
    {
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

        // REMARK: when creating the condition here, the integration rule is not passed. Instead the default integration rule of this element_type is applied, which may be not the same as the original element.
        Condition::Pointer pNewCondition = r_clone_condition.Create(++lastConditionId, pElement->pGetGeometry(), pProperties);

        std::cout << "1 condition of type " << sample_condition_name << " is created" << std::endl;

        return pNewCondition;
    }

    /// List the nodes, elements, conditions of a model_part
    template<class TModelPartType>
    static void ListModelPart(TModelPartType& r_model_part)
    {
        for (auto it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            std::cout << "Node #" << it->Id() << ": ("
                      << it->X0() << ", " << it->Y0() << ", " << it->Z0() << ")" << std::endl;
        }

        for (auto it = r_model_part.Elements().begin(); it != r_model_part.Elements().end(); ++it)
        {
            std::cout << typeid(*it).name() << ": " << it->Id() << std::endl;
        }

        for (auto it = r_model_part.Conditions().begin(); it != r_model_part.Conditions().end(); ++it)
        {
            std::cout << typeid(*it).name() << ": " << it->Id() << std::endl;
        }
    }

    /// Get the equation id of a dof associated with a node
    template<typename TVariableType>
    static std::size_t GetEquationId(ModelPart::NodeType& rNode, const TVariableType& rVariable)
    {
        return rNode.GetDof(rVariable).EquationId();
    }

    /// Compute the normal vector of a surface patch
    template<typename TCoordinatesType>
    static array_1d<double, 3> ComputeNormal(const Patch<2>& rPatch, const TCoordinatesType& xi, const bool& normalize = true)
    {
        auto tangents = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES)->GetDerivative(xi);

        const auto& a = tangents[0];
        const auto& b = tangents[1];
        array_1d<double, 3> normal;
        normal[0] = a[1] * b[2] - a[2] * b[1];
        normal[1] = a[2] * b[0] - a[0] * b[2];
        normal[2] = a[0] * b[1] - a[1] * b[0];

        if (normalize)
        {
            double norm = norm_2(normal);
            normal *= 1.0 / norm;
        }

        return normal;
    }

    /// Check if normal on the side of the patch is pointing outward
    /// This is done by comparing the vector connecting the patch inner point and face center point
    /// and the face normal
    static bool IsNormalPointingOutward(const Patch<3>& rPatch, const BoundarySide& side, double delta = 1.0e-2)
    {
        Patch<2>::Pointer pBoundaryPatch = rPatch.ConstructBoundaryPatch(side);

        const std::vector<double> xic = {0.5, 0.5, 0.5};

        const std::vector<std::vector<double> > xic_boundary =
        {
            {0.0, 0.5, 0.5}, {1.0, 0.5, 0.5},
            {0.5, 0.0, 0.5}, {0.5, 1.0, 0.5},
            {0.5, 0.5, 0.0}, {0.5, 0.5, 1.0}
        };

        array_1d<double, 3> face_center = pBoundaryPatch->pGetGridFunction(CONTROL_POINT_COORDINATES)->GetValue(xic);
        array_1d<double, 3> face_normal = ComputeNormal(*pBoundaryPatch, xic, true);
        // KRATOS_WATCH(face_center)
        // KRATOS_WATCH(face_normal)

        // we find the local coordinates of the face center in the local frame of the patch by look-up technique.
        // Since we know that the face center can only have values contained in xic_boundary
        array_1d<double, 3> ref_point;
        array_1d<double, 3> ref_local_point;
        bool found = false;
        for (std::size_t i = 0; i < xic_boundary.size(); ++i)
        {
            noalias(ref_point) = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES)->GetValue(xic_boundary[i]);
            if (norm_2(ref_point - face_center) < 1.0e-10)
            {
                for (int j = 0; j < 3; ++j)
                {
                    ref_local_point[j] = xic_boundary[i][j];
                }
                found = true;
                break;
            }
        }
        if (found)
        {
            // if we find the local coordinates, we adjust the other coordinate (0.0 or 1.0) with a small pertubation
            for (int j = 0; j < 3; ++j)
            {
                if (std::abs(ref_local_point[j]) < 1.0e-10)
                {
                    ref_local_point[j] += delta;
                }
                else if (std::abs(ref_local_point[j] - 1.0) < 1.0e-10)
                {
                    ref_local_point[j] -= delta;
                }
            }
        }
        else
        {
            // if we can't find it, which is not likely to happen, we take the center point of the patch as a reference point
            for (int j = 0; j < 3; ++j)
            {
                ref_local_point[j] = xic[j];
            }
        }
        // KRATOS_WATCH(ref_local_point)
        noalias(ref_point) = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES)->GetValue(ref_local_point);
        // KRATOS_WATCH(ref_point)

        array_1d<double, 3> v = face_center - ref_point;

        return (inner_prod(v, face_normal) > 1.0e-13);
    }

    /// Information
    template<class TClassType>
    static void PrintAddress(std::ostream& rOStream, typename TClassType::Pointer pInstance)
    {
        rOStream << pInstance << std::endl;
    }

    virtual std::string Info() const
    {
        return "IsogeometricUtility";
    }

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

}; // end class IsogeometricUtility

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const IsogeometricUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ISOGEOMETRIC_UTILITY_H_INCLUDED defined
