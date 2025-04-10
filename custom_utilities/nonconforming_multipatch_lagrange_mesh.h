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
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;
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

    /// Set the base condition name
    void SetBaseConditionName(const std::string& BaseConditionName) {mBaseConditionName = BaseConditionName;}

    /// Set the last node index
    void SetLastNodeId(IndexType lastNodeId) {mLastNodeId = lastNodeId;}

    /// Set the last element index
    void SetLastElemId(IndexType lastElemId) {mLastElemId = lastElemId;}

    /// Set the last condition index
    void SetLastCondId(IndexType lastCondId) {mLastCondId = lastCondId;}

    /// Mark the side of the patch where the conditions are created
    void MarkConditionFace(IndexType patch_id, const BoundarySide side, const IndexType prop_id)
    {
        mConditionMap.push_back(std::make_pair(patch_id, std::make_pair(side, prop_id)));
    }

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

        // get the sample condition
        std::string condition_name = mBaseConditionName;
        if constexpr (TDim == 2)
        {
            condition_name = condition_name + "2D2N";
        }
        else if constexpr (TDim == 3)
        {
            condition_name = condition_name + "3D4N";
        }

        if (!KratosComponents<Condition>::Has(condition_name))
        {
            KRATOS_ERROR << "Condition " << condition_name << " is not registered in Kratos."
                         << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
        }

        Condition const& rCloneCondition = KratosComponents<Condition>::Get(condition_name);

        // generate nodes and elements for each patch
        IndexType NodeCounter = mLastNodeId + 1;
        IndexType ElementCounter = mLastElemId + 1;
        IndexType ConditionCounter = mLastCondId + 1;
        for (auto it = mpMultiPatch->ptr_begin(); it != mpMultiPatch->ptr_end(); ++it)
        {
            if (!(*it)->Is(ACTIVE))
            {
                continue;
            }

            // create new properties and add to model_part
            if ((*it)->LayerIndex() < 0)
            {
                KRATOS_WATCH((*it)->LayerIndex())
                KRATOS_WATCH((*it)->Id())
                KRATOS_WATCH((*it)->Prefix())
                KRATOS_ERROR << "Patch " << (*it)->Id() << " has invalid layer index " << (*it)->LayerIndex();
            }
            Properties::Pointer pNewProperties = r_model_part.pGetProperties((*it)->LayerIndex());

            this->AddElements(r_model_part, *it, rCloneElement, NodeCounter, ElementCounter, pNewProperties);
        }

        // create and add conditions on the boundary
        for (std::size_t i = 0; i < mConditionMap.size(); ++i)
        {
            typename Patch<TDim>::Pointer cond_patch = mpMultiPatch->pGetPatch(mConditionMap[i].first);
            const BoundarySide side = mConditionMap[i].second.first;
            const IndexType prop_id = mConditionMap[i].second.second;

            Properties::Pointer pProperties = r_model_part.pGetProperties(prop_id);

            this->AddConditions(r_model_part, cond_patch, side, rCloneCondition, NodeCounter, ConditionCounter, pProperties);
        }

        if (GetEchoLevel() > 0)
        {
            std::cout << "NonConformingMultipatchLagrangeMesh::" << __FUNCTION__ << " - completed" << std::endl;
        }
    }

    /// create the elements out from the patch and add to the model_part
    ElementsContainerType AddElements(ModelPart& r_model_part, typename Patch<TDim>::Pointer pPatch, Element const& rCloneElement,
            std::size_t& starting_node_id, std::size_t& starting_elem_id, Properties::Pointer pProperties) const
    {
        const std::string NodeKey = std::string("Node");

        // generate nodes and elements for the patch
        IndexType& NodeCounter = starting_node_id;
        IndexType& ElementCounter = starting_elem_id;
        if (!pPatch->Is(ACTIVE))
        {
            ElementsContainerType Dummy;
            return Dummy;
        }

        if (GetEchoLevel() > 1)
        {
            std::cout << "Elements will be created on patch " << pPatch->Id() << std::endl;
        }

        // generate the points and connectivities
        std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;

        if constexpr (TDim == 2)
        {
            // create new nodes and elements
            auto it_num = mNumDivision.find(pPatch->Id());
            if (it_num == mNumDivision.end())
            {
                KRATOS_ERROR << "NumDivision is not set for patch " << pPatch->Id();
            }

            IndexType NumDivision1 = it_num->second[0];
            IndexType NumDivision2 = it_num->second[1];
            if (GetEchoLevel() > 1)
                std::cout << "Divisioning for patch " << pPatch->Id() << ": " << NumDivision1
                          << " " << NumDivision2 << std::endl;

            std::vector<array_1d<double, 3> > corners(4);

            IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                                        corners[2], corners[3], NodeCounter, NumDivision1, NumDivision2);
        }
        else if constexpr (TDim == 3)
        {
            // create new nodes and elements
            const auto it_num = mNumDivision.find(pPatch->Id());
            if (it_num == mNumDivision.end())
            {
                KRATOS_ERROR << "NumDivision is not set for patch " << pPatch->Id();
            }

            IndexType NumDivision1 = it_num->second[0];
            IndexType NumDivision2 = it_num->second[1];
            IndexType NumDivision3 = it_num->second[2];
            if (GetEchoLevel() > 1)
                std::cout << "Divisioning for patch " << pPatch->Id() << ": " << NumDivision1
                          << " " << NumDivision2 << " " << NumDivision3 << std::endl;

            std::vector<array_1d<double, 3> > corners(8);

            IsogeometricPostUtility::GenerateBox(corners, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateHexGrid(corners[0], corners[1], corners[2], corners[3],
                                        corners[4], corners[5], corners[6], corners[7], NodeCounter, NumDivision1, NumDivision2, NumDivision3);
        }

        // create nodes and transfer values to nodes
        for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
        {
            IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *pPatch, r_model_part, NodeCounter++);
        }

        // create elements
        ElementsContainerType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsContainerType, 1>(
                                             points_and_connectivities.second, r_model_part, rCloneElement, ElementCounter, pProperties, NodeKey);

        for (typename ElementsContainerType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
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

        return pNewElements;
    }

    /// create the conditions out from the boundary patch and add to the model_part
    ConditionsContainerType AddConditions(ModelPart& r_model_part, typename Patch<TDim>::Pointer pPatch, const BoundarySide side,
            Condition const& rCloneCondition, std::size_t& starting_node_id, std::size_t& starting_cond_id, Properties::Pointer pProperties) const
    {
        // construct the boundary patch
        const auto pBoundaryPatch = pPatch->ConstructBoundaryPatch(side);

        // get the sampling numbers
        auto it_num = mNumDivision.find(pPatch->Id());
        if (it_num == mNumDivision.end())
        {
            KRATOS_ERROR << "NumDivision is not set for patch " << pPatch->Id();
        }

        boost::array<double, TDim-1> div;
        if constexpr (TDim == 2)
        {
            const auto pdir = ParameterDirection<TDim>::Get(side);
            div[0] = it_num->second[pdir[0]];
        }
        else if constexpr (TDim == 3)
        {
            const auto pdir = ParameterDirection<TDim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            div[1] = it_num->second[pdir[1]];
        }

        // create conditions
        return this->AddConditions(r_model_part, pBoundaryPatch, div, rCloneCondition, starting_node_id, starting_cond_id, pProperties);
    }

    /// create the conditions out from the slice patch and add to the model_part
    ConditionsContainerType AddConditions(ModelPart& r_model_part, typename Patch<TDim>::Pointer pPatch, const int idir, const double xi,
            Condition const& rCloneCondition, std::size_t& starting_node_id, std::size_t& starting_cond_id, Properties::Pointer pProperties) const
    {
        // extract the slice patch
        const auto pSlicedPatch = pPatch->ConstructSlicedPatch(idir, xi);

        // get the sampling numbers
        auto it_num = mNumDivision.find(pPatch->Id());
        if (it_num == mNumDivision.end())
        {
            KRATOS_ERROR << "NumDivision is not set for patch " << pPatch->Id();
        }

        boost::array<double, TDim-1> div;
        const BoundarySide side = ParameterDirection<TDim>::GetSide(idir);
        if constexpr (TDim == 2)
        {
            const auto pdir = ParameterDirection<TDim>::Get(side);
            div[0] = it_num->second[pdir[0]];
        }
        else if constexpr (TDim == 3)
        {
            const auto pdir = ParameterDirection<TDim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            div[1] = it_num->second[pdir[1]];
        }

        // create conditions
        return this->AddConditions(r_model_part, pSlicedPatch, div, rCloneCondition, starting_node_id, starting_cond_id, pProperties);
    }

    /// create the conditions out from the boudnary patch and add to the model_part
    ConditionsContainerType AddConditions(ModelPart& r_model_part, typename Patch<TDim-1>::Pointer pBoundaryPatch,
            const boost::array<double, TDim-1>& nsampling,
            Condition const& rCloneCondition,
            std::size_t& starting_node_id, std::size_t& starting_cond_id, Properties::Pointer pProperties) const
    {
        std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;
        IndexType& NodeCounter = starting_node_id;
        IndexType& ConditionCounter = starting_cond_id;

        // generate points and connectivities
        if constexpr (TDim == 2)
        {
            std::vector<array_1d<double, 3> > corners(2);

            IsogeometricPostUtility::GenerateLine(corners, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateLineGrid(corners[0], corners[1], NodeCounter, nsampling[0]);
        }
        else if constexpr (TDim == 3)
        {
            std::vector<array_1d<double, 3> > corners(4);

            IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                                                corners[2], corners[3], NodeCounter, nsampling[0], nsampling[1]);
        }

        // create nodes and transfer values to nodes
        for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
        {
            IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *pBoundaryPatch, r_model_part, NodeCounter++);
        }

        // create conditions
        const std::string NodeKey = std::string("Node");
        ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Condition, ConditionsContainerType, 1>(
                                             points_and_connectivities.second, r_model_part, rCloneCondition, ConditionCounter, pProperties, NodeKey);

        for (typename ConditionsContainerType::ptr_iterator it2 = pNewConditions.ptr_begin(); it2 != pNewConditions.ptr_end(); ++it2)
        {
            r_model_part.AddCondition(*it2);
            if (GetEchoLevel() > 2)
            {
                std::cout << "Condition " << (*it2)->Id() << " is created with connectivity:";
                for (std::size_t n = 0; n < (*it2)->GetGeometry().size(); ++n)
                {
                    std::cout << " " << (*it2)->GetGeometry()[n].Id();
                }
                std::cout << std::endl;
            }
        }

        // just to make sure everything is organized properly
        r_model_part.Conditions().Unique();

        return pNewConditions;
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

    std::vector<std::pair<IndexType, std::pair<BoundarySide, IndexType> > > mConditionMap;

    std::string mBaseElementName;
    std::string mBaseConditionName;
    IndexType mLastNodeId;
    IndexType mLastElemId;
    IndexType mLastCondId;
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
