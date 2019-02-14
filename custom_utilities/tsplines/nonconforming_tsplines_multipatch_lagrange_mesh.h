//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_TSPLINES_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_TSPLINES_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED

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
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/nonconforming_multipatch_lagrange_mesh.h"

// #define DEBUG_MESH_GENERATION

namespace Kratos
{

/**
Construct the standard FEM mesh based on Lagrange basis functions from isogeometric T-Splines multipatch. Each patch can have different division and is non-conformed at the boundary.
The principle is that each patch will be sampled based on number of divisions defined by used. Therefore user is not able to see the knot density.
At the end, the resulting model_part will have nodal values interpolated from patch. This class is useful for post-processing all types of isogeometric patches, including NURBS, hierarchical B-Splines and T-Splines.
 */
template<class TFESpaceType>
class NonConformingTSplinesMultipatchLagrangeMesh
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NonConformingTSplinesMultipatchLagrangeMesh);

    /// Type definition
    typedef typename Element::GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename Element::GeometryType::PointType NodeType;

    /// Default constructor
    NonConformingTSplinesMultipatchLagrangeMesh(typename MultiPatch<TFESpaceType::Dim()>::Pointer pMultiPatch) : mpMultiPatch(pMultiPatch) {}

    /// Destructor
    virtual ~NonConformingTSplinesMultipatchLagrangeMesh() {}

    /// Set the division for all the patches the same number of division in each dimension
    /// Note that if the division is changed, the post_model_part must be generated again
    void SetUniformDivision(const std::size_t& num_division)
    {
        typedef typename MultiPatch<TFESpaceType::Dim()>::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            for (std::size_t dim = 0; dim < TFESpaceType::Dim(); ++dim)
                mNumDivision[it->Id()][dim] = num_division;
        }

    }

    /// Set the division for the patch at specific dimension
    /// Note that if the division is changed, the post_model_part must be generated again
    void SetDivision(const std::size_t& patch_id, const int& dim, const std::size_t& num_division)
    {
        if (mpMultiPatch->Patches().find(patch_id) == mpMultiPatch->end())
        {
            std::stringstream ss;
            ss << "Patch " << patch_id << " is not found in the multipatch";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

        mNumDivision[patch_id][dim] = num_division;
    }

    /// Set the base element name
    void SetBaseElementName(const std::string& BaseElementName) {mBaseElementName = BaseElementName;}

    /// Set the last node index
    void SetLastNodeId(const std::size_t& lastNodeId) {mLastNodeId = lastNodeId;}

    /// Set the last element index
    void SetLastElemId(const std::size_t& lastElemId) {mLastElemId = lastElemId;}

    /// Set the last properties index
    void SetLastPropId(const std::size_t& lastPropId) {mLastPropId = lastPropId;}

    /// Append to model_part, the quad/hex element from patches
    void WriteModelPart(ModelPart& r_model_part) const
    {
        // get the sample element
        std::string element_name = mBaseElementName;
        if (TFESpaceType::Dim() == 2)
            element_name = element_name + "2D4N";
        else if (TFESpaceType::Dim() == 3)
            element_name = element_name + "3D8N";

        std::string NodeKey = std::string("Node");

        if(!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
            KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
        }

        Element const& rCloneElement = KratosComponents<Element>::Get(element_name);

        // generate nodes and elements for each patch
        std::size_t NodeCounter = mLastNodeId;
        std::size_t NodeCounter_old = NodeCounter;
        std::size_t ElementCounter = mLastElemId;
        // std::size_t PropertiesCounter = mLastPropId;
        std::vector<double> p_ref(TFESpaceType::Dim());
        typedef typename MultiPatch<TFESpaceType::Dim()>::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            // create new properties and add to model_part
            // Properties::Pointer pNewProperties = Properties::Pointer(new Properties(PropertiesCounter++));
            Properties::Pointer pNewProperties = Properties::Pointer(new Properties(it->Id()));
            r_model_part.AddProperties(pNewProperties);

            typename TFESpaceType::cell_container_t::Pointer pFaceManager = boost::dynamic_pointer_cast<TFESpaceType>(it->pFESpace())->pFaceManager();

            if (TFESpaceType::Dim() == 2)
            {
                typename std::map<std::size_t, boost::array<std::size_t, TFESpaceType::Dim()> >::const_iterator it_num = mNumDivision.find(it->Id());
                if (it_num == mNumDivision.end())
                    KRATOS_THROW_ERROR(std::logic_error, "NumDivision is not set for patch", it->Id())

                std::size_t NumDivision1 = it_num->second[0];
                std::size_t NumDivision2 = it_num->second[1];
                #ifdef DEBUG_MESH_GENERATION
                KRATOS_WATCH(NumDivision1)
                KRATOS_WATCH(NumDivision2)
                #endif

                for (typename TFESpaceType::cell_iterator it_cell = pFaceManager->begin(); it_cell != pFaceManager->end(); ++it_cell)
                {
                    double xi_min = (*it_cell)->XiMin();
                    double xi_max = (*it_cell)->XiMax();
                    double eta_min = (*it_cell)->EtaMin();
                    double eta_max = (*it_cell)->EtaMax();

                    // create new nodes for each face
                    for (std::size_t i = 0; i <= NumDivision1; ++i)
                    {
                        p_ref[0] = ((double) i) / NumDivision1 * (xi_max - xi_min) + xi_min;
                        for (std::size_t j = 0; j <= NumDivision2; ++j)
                        {
                            p_ref[1] = ((double) j) / NumDivision2 * (eta_max - eta_min) + eta_min;

                            #ifdef DEBUG_MESH_GENERATION
                            std::cout << "p_ref: " << p_ref[0] << " " << p_ref[1] << std::endl;
                            #endif

                            NonConformingMultipatchLagrangeMesh<TFESpaceType::Dim()>::CreateNode(p_ref, *it, r_model_part, NodeCounter);
                            ++NodeCounter;
                        }
                    }

                    // create and add element
                    Element::NodesArrayType temp_element_nodes;
                    for (std::size_t i = 0; i < NumDivision1; ++i)
                    {
                        for (std::size_t j = 0; j < NumDivision2; ++j)
                        {
                            std::size_t Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j;
                            std::size_t Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
                            std::size_t Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j;
                            std::size_t Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;

                            // TODO: check if jacobian checking is necessary
                            temp_element_nodes.clear();
                            temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node1, NodeKey).base()));
                            temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node2, NodeKey).base()));
                            temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node4, NodeKey).base()));
                            temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node3, NodeKey).base()));

                            Element::Pointer pNewElement = rCloneElement.Create(ElementCounter++, temp_element_nodes, pNewProperties);
                            r_model_part.AddElement(pNewElement);
                            #ifdef DEBUG_MESH_GENERATION
                            std::cout << "Element " << pNewElement->Id() << " is created with connectivity:";
                            for (std::size_t n = 0; n < pNewElement->GetGeometry().size(); ++n)
                                std::cout << " " << pNewElement->GetGeometry()[n].Id();
                            std::cout << std::endl;
                            #endif
                        }
                    }

                    // create and add conditions on the boundary
                    // TODO

                    // update the node counter
                    NodeCounter_old = NodeCounter;
                }

                // just to make sure everything is organized properly
                r_model_part.Elements().Unique();
            }
            else if (TFESpaceType::Dim() == 3)
            {
                // create new nodes
                typename std::map<std::size_t, boost::array<std::size_t, TFESpaceType::Dim()> >::const_iterator it_num = mNumDivision.find(it->Id());
                if (it_num == mNumDivision.end())
                    KRATOS_THROW_ERROR(std::logic_error, "NumDivision is not set for patch", it->Id())

                std::size_t NumDivision1 = it_num->second[0];
                std::size_t NumDivision2 = it_num->second[1];
                std::size_t NumDivision3 = it_num->second[2];
                #ifdef DEBUG_MESH_GENERATION
                KRATOS_WATCH(NumDivision1)
                KRATOS_WATCH(NumDivision2)
                KRATOS_WATCH(NumDivision3)
                #endif

                for (typename TFESpaceType::cell_iterator it_cell = pFaceManager->begin(); it_cell != pFaceManager->end(); ++it_cell)
                {
                    double xi_min = (*it_cell)->XiMin();
                    double xi_max = (*it_cell)->XiMax();
                    double eta_min = (*it_cell)->EtaMin();
                    double eta_max = (*it_cell)->EtaMax();
                    double zeta_min = (*it_cell)->ZetaMin();
                    double zeta_max = (*it_cell)->ZetaMax();

                    for (std::size_t i = 0; i <= NumDivision1; ++i)
                    {
                        p_ref[0] = ((double) i) / NumDivision1 * (xi_max - xi_min) + xi_min;
                        for (std::size_t j = 0; j <= NumDivision2; ++j)
                        {
                            p_ref[1] = ((double) j) / NumDivision2 * (eta_max - eta_min) + eta_min;
                            for (std::size_t k = 0; k <= NumDivision3; ++k)
                            {
                                p_ref[2] = ((double) k) / NumDivision3 * (zeta_max - zeta_min) + zeta_min;

                                NonConformingMultipatchLagrangeMesh<TFESpaceType::Dim()>::CreateNode(p_ref, *it, r_model_part, NodeCounter);
                                ++NodeCounter;
                            }
                        }
                    }

                    // create and add element
                    Element::NodesArrayType temp_element_nodes;
                    for (std::size_t i = 0; i < NumDivision1; ++i)
                    {
                        for (std::size_t j = 0; j < NumDivision2; ++j)
                        {
                            for (std::size_t k = 0; k < NumDivision3; ++k)
                            {
                                std::size_t Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k;
                                std::size_t Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k;
                                std::size_t Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k;
                                std::size_t Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k;
                                std::size_t Node5 = Node1 + 1;
                                std::size_t Node6 = Node2 + 1;
                                std::size_t Node7 = Node3 + 1;
                                std::size_t Node8 = Node4 + 1;

                                // TODO: check if jacobian checking is necessary
                                temp_element_nodes.clear();
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node1, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node2, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node4, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node3, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node5, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node6, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node8, NodeKey).base()));
                                temp_element_nodes.push_back(*(MultiPatchUtility::FindKey(r_model_part.Nodes(), Node7, NodeKey).base()));

                                Element::Pointer pNewElement = rCloneElement.Create(ElementCounter++, temp_element_nodes, pNewProperties);
                                r_model_part.AddElement(pNewElement);
                                #ifdef DEBUG_MESH_GENERATION
                                std::cout << "Element " << pNewElement->Id() << " is created with connectivity:";
                                for (std::size_t n = 0; n < pNewElement->GetGeometry().size(); ++n)
                                    std::cout << " " << pNewElement->GetGeometry()[n].Id();
                                std::cout << std::endl;
                                #endif
                            }
                        }
                    }

                    // create and add conditions on the boundary
                    // TODO

                    // update the node counter
                    NodeCounter_old = NodeCounter;
                }

                // just to make sure everything is organized properly
                r_model_part.Elements().Unique();
            }
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NonConformingTSplinesMultipatchLagrangeMesh";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    typename MultiPatch<TFESpaceType::Dim()>::Pointer mpMultiPatch;

    std::map<std::size_t, boost::array<std::size_t, TFESpaceType::Dim()> > mNumDivision;

    std::string mBaseElementName;
    std::size_t mLastNodeId;
    std::size_t mLastElemId;
    std::size_t mLastPropId;

};

/// output stream function
template<class TFESpaceType>
inline std::ostream& operator <<(std::ostream& rOStream, const NonConformingTSplinesMultipatchLagrangeMesh<TFESpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_MESH_GENERATION

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NONCONFORMING_MULTIPATCH_LAGRANGE_MESH_H_INCLUDED defined

