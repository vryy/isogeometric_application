//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Oct 2019 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_LAGRANGE_MESH_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_LAGRANGE_MESH_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/isogeometric_post_utility.h"

// #define DEBUG_MESH_GENERATION

namespace Kratos
{

/**
Construct the standard FEM mesh based on Lagrange basis functions from isogeometric patch.
At the end, the resulting model_part will have nodal values interpolated from patch. This class is useful for post-processing all types of isogeometric patches, including NURBS, hierarchical B-Splines and T-Splines.
 */
template<int TDim>
class PatchLagrangeMesh
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PatchLagrangeMesh);

    /// Type definition
    typedef typename Element::GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename Element::GeometryType::PointType NodeType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef std::size_t IndexType;

    /// Default constructor
    PatchLagrangeMesh()
    {}

    /// Destructor
    virtual ~PatchLagrangeMesh() {}

    /// Append to model_part, the quad/hex element from patches
    template<typename TEntityType, typename TEntityContainerType>
    static void WriteEntities(ModelPart& r_model_part, TEntityContainerType& r_entities,
        typename Patch<TDim>::Pointer pPatch, const TEntityType& r_clone_entity,
        const std::vector<std::size_t>& num_divisions,
        std::size_t& last_node_id, std::size_t& last_entity_id,
        Properties::Pointer pProperties, int echo_level = 0)
    {
        if (echo_level > 0)
            std::cout << "invoking PatchLagrangeMesh::" << __FUNCTION__ << std::endl;

        if (num_divisions.size() < TDim)
            KRATOS_THROW_ERROR(std::logic_error, "Insufficient number of division", "")

        // generate nodes and elements for each patch

        if (echo_level > 1)
            std::cout << "Elements/Conditions will be created on patch " << pPatch->Id() << std::endl;

        // generate the connectivities
        std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;

        if (TDim == 2)
        {
            // create new nodes and elements
            if (echo_level > 1)
                std::cout << "Divisioning for patch " << pPatch->Id() << ": " << num_divisions[0]
                          << " " << num_divisions[1] << std::endl;

            std::vector<array_1d<double, 3> > corners(4);

            IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                    corners[2], corners[3], ++last_node_id, num_divisions[0], num_divisions[1]);
        }
        else if (TDim == 3)
        {
            // create new nodes and elements
            if (echo_level > 1)
                std::cout << "Divisioning for patch " << pPatch->Id() << ": " << num_divisions[0]
                          << " " << num_divisions[1] << " " << num_divisions[2] << std::endl;

            std::vector<array_1d<double, 3> > corners(8);

            IsogeometricPostUtility::GenerateBox(corners, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateHexGrid(corners[0], corners[1], corners[2], corners[3],
                    corners[4], corners[5], corners[6], corners[7], ++last_node_id, num_divisions[0], num_divisions[1], num_divisions[2]);
        }

        // create nodes
        for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
            IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *pPatch, r_model_part, last_node_id++);

        // create elements
        const std::string NodeKey = std::string("Node");
        TEntityContainerType pNewEntities = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, TEntityType, TEntityContainerType>(
            points_and_connectivities.second, r_model_part, r_clone_entity, last_entity_id, pProperties, NodeKey);

        for (typename TEntityContainerType::ptr_iterator it2 = pNewEntities.ptr_begin(); it2 != pNewEntities.ptr_end(); ++it2)
        {
            r_entities.push_back(*it2);
            if (echo_level > 2)
            {
                std::cout << "Element/Condition " << (*it2)->Id() << " is created with connectivity:";
                for (std::size_t n = 0; n < (*it2)->GetGeometry().size(); ++n)
                    std::cout << " " << (*it2)->GetGeometry()[n].Id();
                std::cout << std::endl;
            }
        }

        // just to make sure everything is organized properly
        r_entities.Unique();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PatchLagrangeMesh<" << TDim << ">";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const PatchLagrangeMesh<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_MESH_GENERATION

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_LAGRANGE_MESH_H_INCLUDED defined

