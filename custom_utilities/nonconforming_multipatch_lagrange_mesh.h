//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/isogeometric_post_utility.h"

// #define DEBUG_MESH_GENERATION

namespace Kratos
{

/**
 * Construct the standard FEM mesh based on Lagrange basis functions from isogeometric multipatch.
 * Each patch can have different division and is non-conformed at the boundary.
 * The principle is that each patch will be sampled based on number of divisions defined by user.
 * Therefore user is not able to see the knot density.
 * At the end, the resulting model_part will have nodal values interpolated from patch.
 * This class is useful for post-processing all types of isogeometric patches, including NURBS,
 * hierarchical B-Splines and T-Splines.
 */
template<int TDim>
class NonConformingMultipatchLagrangeMesh : public IsogeometricEcho
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NonConformingMultipatchLagrangeMesh);

    /// Type definition
    typedef typename Element::GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename Element::GeometryType::PointType NodeType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef std::size_t IndexType;

    /// Default constructor
    NonConformingMultipatchLagrangeMesh(typename MultiPatch<TDim>::Pointer pMultiPatch)
        : mpMultiPatch(pMultiPatch), mLastNodeId(0), mLastElemId(0)
    {}

    /// Destructor
    virtual ~NonConformingMultipatchLagrangeMesh() {}

    /// Set the division for all the patches the same number of division in each dimension
    /// Note that if the division is changed, the post_model_part must be generated again
    void SetUniformDivision(IndexType num_division)
    {
        typedef typename MultiPatch<TDim>::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            for (IndexType dim = 0; dim < TDim; ++dim)
            {
                mNumDivision[it->Id()][dim] = num_division;
            }
        }
    }

    /// Set the division for the patch at specific dimension
    /// Note that if the division is changed, the post_model_part must be generated again
    void SetDivision(IndexType patch_id, int dim, IndexType num_division)
    {
        if (dim >= TDim)
        {
            KRATOS_ERROR << "The dimension " << dim << " is invalid";
        }

        if (mpMultiPatch->Patches().find(patch_id) == mpMultiPatch->end())
        {
            KRATOS_ERROR << "Patch " << patch_id << " is not found in the multipatch";
        }

        mNumDivision[patch_id][dim] = num_division;
    }

    /// Set the base element name
    void SetBaseElementName(const std::string& BaseElementName) {mBaseElementName = BaseElementName;}

    /// Set the last node index
    void SetLastNodeId(IndexType lastNodeId) {mLastNodeId = lastNodeId;}

    /// Set the last element index
    void SetLastElemId(IndexType lastElemId) {mLastElemId = lastElemId;}

    /// Append to model_part, the quad/hex element from patches
    void WriteModelPart(ModelPart& r_model_part) const
    {
        if (GetEchoLevel() > 0)
        {
            std::cout << "NonConformingMultipatchLagrangeMesh::" << __FUNCTION__ << " - started" << std::endl;
        }

        // get the sample element
        std::string element_name = mBaseElementName;
        if constexpr (TDim == 2)
        {
            element_name = element_name + "2D4N";
        }
        else if constexpr (TDim == 3)
        {
            element_name = element_name + "3D8N";
        }

        const std::string NodeKey = std::string("Node");

        if (!KratosComponents<Element>::Has(element_name))
        {
            KRATOS_ERROR << "Element " << element_name << " is not registered in Kratos."
                         << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
        }

        Element const& rCloneElement = KratosComponents<Element>::Get(element_name);

        // generate nodes and elements for each patch
        IndexType NodeCounter = mLastNodeId + 1;
        IndexType ElementCounter = mLastElemId;
        typedef typename MultiPatch<TDim>::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            if (!it->Is(ACTIVE))
            {
                continue;
            }

            if (GetEchoLevel() > 1)
            {
                std::cout << "Elements will be created on patch " << it->Id() << std::endl;
            }

            // create new properties and add to model_part
            if (it->LayerIndex() < 0)
            {
                KRATOS_WATCH(it->LayerIndex())
                KRATOS_WATCH(it->Id())
                KRATOS_WATCH(it->Prefix())
                KRATOS_ERROR << "Patch " << it->Id() << " has invalid layer index " << it->LayerIndex();
            }
            Properties::Pointer pNewProperties = r_model_part.pGetProperties(it->LayerIndex());

            // generate the connectivities
            std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;

            if constexpr (TDim == 2)
            {
                // create new nodes and elements
                typename std::map<IndexType, boost::array<IndexType, TDim> >::const_iterator it_num = mNumDivision.find(it->Id());
                if (it_num == mNumDivision.end())
                {
                    KRATOS_ERROR << "NumDivision is not set for patch " << it->Id();
                }

                IndexType NumDivision1 = it_num->second[0];
                IndexType NumDivision2 = it_num->second[1];
                if (GetEchoLevel() > 1)
                    std::cout << "Divisioning for patch " << it->Id() << ": " << NumDivision1
                              << " " << NumDivision2 << std::endl;

                std::vector<array_1d<double, 3> > corners(4);

                IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

                points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                                            corners[2], corners[3], NodeCounter, NumDivision1, NumDivision2);
            }
            else if constexpr (TDim == 3)
            {
                // create new nodes and elements
                typename std::map<IndexType, boost::array<IndexType, TDim> >::const_iterator it_num = mNumDivision.find(it->Id());
                if (it_num == mNumDivision.end())
                {
                    KRATOS_ERROR << "NumDivision is not set for patch " << it->Id();
                }

                IndexType NumDivision1 = it_num->second[0];
                IndexType NumDivision2 = it_num->second[1];
                IndexType NumDivision3 = it_num->second[2];
                if (GetEchoLevel() > 1)
                    std::cout << "Divisioning for patch " << it->Id() << ": " << NumDivision1
                              << " " << NumDivision2 << " " << NumDivision3 << std::endl;

                std::vector<array_1d<double, 3> > corners(8);

                IsogeometricPostUtility::GenerateBox(corners, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

                points_and_connectivities = IsogeometricPostUtility::GenerateHexGrid(corners[0], corners[1], corners[2], corners[3],
                                            corners[4], corners[5], corners[6], corners[7], NodeCounter, NumDivision1, NumDivision2, NumDivision3);
            }

            // create nodes and transfer values to nodes
            for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
            {
                IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *it, r_model_part, NodeCounter++);
            }

            // create elements
            ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                                                 points_and_connectivities.second, r_model_part, rCloneElement, ElementCounter, pNewProperties, NodeKey);

            for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
            {
                r_model_part.AddElement(*it2);
                if (GetEchoLevel() > 2)
                {
                    std::cout << "Element " << (*it2)->Id() << " is created with connectivity:";
                    for (std::size_t n = 0; n < (*it2)->GetGeometry().size(); ++n)
                    {
                        std::cout << " " << (*it2)->GetGeometry()[n].Id();
                    }
                    std::cout << std::endl;
                }
            }

            // just to make sure everything is organized properly
            r_model_part.Elements().Unique();

            // create and add conditions on the boundary
            // TODO
        }

        if (GetEchoLevel() > 0)
        {
            std::cout << "NonConformingMultipatchLagrangeMesh::" << __FUNCTION__ << " - completed" << std::endl;
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NonConformingMultipatchLagrangeMesh<" << TDim << ">";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    typename MultiPatch<TDim>::Pointer mpMultiPatch;

    std::map<IndexType, boost::array<IndexType, TDim> > mNumDivision;

    std::string mBaseElementName;
    IndexType mLastNodeId;
    IndexType mLastElemId;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const NonConformingMultipatchLagrangeMesh<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_MESH_GENERATION

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED defined
