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
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef std::size_t IndexType;

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

                    // create nodes and elements
                    std::vector<array_1d<double, 3> > corners(4);

                    IsogeometricPostUtility::GenerateRectangle(corners, xi_min, xi_max, eta_min, eta_max);

                    std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities
                        = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1], corners[2], corners[3],
                            NodeCounter, NumDivision1, NumDivision2);

                    for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
                    {
                        IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *it, r_model_part, NodeCounter++);
                    }

                    ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                        points_and_connectivities.second, r_model_part, rCloneElement, ElementCounter, pNewProperties, NodeKey);

                    for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
                    {
                        r_model_part.AddElement(*it2);
                        #ifdef DEBUG_MESH_GENERATION
                        std::cout << "Element " << (*it2)->Id() << " is created with connectivity:";
                        for (std::size_t n = 0; n < (*it2)->GetGeometry().size(); ++n)
                            std::cout << " " << (*it2)->GetGeometry()[n].Id();
                        std::cout << std::endl;
                        #endif
                    }

                    // create and add conditions on the boundary
                    // TODO
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

                    // create nodes and elements
                    std::vector<array_1d<double, 3> > corners(8);

                    IsogeometricPostUtility::GenerateBox(corners, xi_min, xi_max, eta_min, eta_max, zeta_min, zeta_max);

                    std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities
                        = IsogeometricPostUtility::GenerateHexGrid(corners[0], corners[1], corners[2], corners[3],
                            corners[4], corners[5], corners[6], corners[7], NodeCounter, NumDivision1, NumDivision2, NumDivision3);

                    for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
                    {
                        IsogeometricPostUtility::CreateNodeAndTransferValues(points_and_connectivities.first[i], *it, r_model_part, NodeCounter++);
                    }

                    ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                        points_and_connectivities.second, r_model_part, rCloneElement, ElementCounter, pNewProperties, NodeKey);

                    for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
                    {
                        r_model_part.AddElement(*it2);
                        #ifdef DEBUG_MESH_GENERATION
                        std::cout << "Element " << (*it2)->Id() << " is created with connectivity:";
                        for (std::size_t n = 0; n < (*it2)->GetGeometry().size(); ++n)
                            std::cout << " " << (*it2)->GetGeometry()[n].Id();
                        std::cout << std::endl;
                        #endif
                    }

                    // create and add conditions on the boundary
                    // TODO
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

