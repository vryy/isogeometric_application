//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2021 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_LAGRANGE_CONTROL_MESH_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_LAGRANGE_CONTROL_MESH_H_INCLUDED

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
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/structured_control_grid.h"

// #define DEBUG_MESH_GENERATION
// #define USE_FUNCTION_ID_FOR_NODE_ID // this option will allow duplicated node with the same id in the post model_part
                // However, it could cause segmentation fault if the model_part cleans the duplicated nodes

namespace Kratos
{

/**
Construct the standard FEM mesh for the control point mesh. The nodal values are taken from the control point values.
 */
template<int TDim>
class MultipatchLagrangeControlMesh : public IsogeometricEcho
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultipatchLagrangeControlMesh);

    /// Type definition
    typedef typename Element::GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename Element::GeometryType::PointType NodeType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef std::size_t IndexType;

    typedef typename Patch<TDim>::ControlPointType ControlPointType;

    /// Default constructor
    MultipatchLagrangeControlMesh(typename MultiPatch<TDim>::Pointer pMultiPatch)
        : mpMultiPatch(pMultiPatch), mEchoLevel(1)
    {}

    /// Destructor
    virtual ~MultipatchLagrangeControlMesh() {}

    /// Set the base element name
    void SetBaseElementName(const std::string& Name) {mBaseElementName = Name;}

    /// Set the base condition name
    void SetBaseConditionName(const std::string& Name) { /*NO NOTHING*/ }

    /// Set the last node index
    void SetLastNodeId(IndexType lastId) {mLastNodeId = lastId;}

    /// Set the last element index
    void SetLastElemId(IndexType lastId) {mLastElemId = lastId;}

    /// Set the last condition index
    void SetLastCondId(IndexType lastId) { /*DO NOTHING*/ }

    /// Append to model_part, the quad/hex element from patches
    void WriteModelPart(ModelPart& r_model_part) const
    {
        if (mEchoLevel > 0)
        {
            std::cout << "invoking MultipatchLagrangeControlMesh::" << __FUNCTION__ << std::endl;
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
        IndexType NodeCounter = mLastNodeId;
        IndexType ElementCounter = mLastElemId;
        typedef typename MultiPatch<TDim>::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            if (!it->Is(ACTIVE))
            {
                continue;
            }

            if (mEchoLevel > 1)
            {
                std::cout << "Elements will be created on patch " << it->Id() << std::endl;
            }

            // create new properties and add to model_part
            if (it->LayerIndex() < 0)
            {
                KRATOS_WATCH(it->LayerIndex())
                KRATOS_WATCH(it->Id())
                KRATOS_WATCH(it->Prefix())
                KRATOS_ERROR << "Invalid layer index " << it->LayerIndex();
            }
            Properties::Pointer pNewProperties = r_model_part.pGetProperties(it->LayerIndex());

            // generate the connectivities
            std::vector<std::vector<IndexType> > connectivities;
            IndexType offset = NodeCounter + 1;

            if (it->pFESpace()->Type() == BSplinesFESpace<TDim>::StaticType())
            {
                typename StructuredControlGrid<TDim, ControlPointType>::Pointer pcontrol_grid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, ControlPointType> >(it->ControlPointGridFunction().pControlGrid());

                const auto func_ids = it->pFESpace()->FunctionIndices();

                #ifdef USE_FUNCTION_ID_FOR_NODE_ID
                std::map<IndexType, IndexType> true_node_ids;
                #endif

                for (std::size_t i = 0; i < pcontrol_grid->size(); ++i)
                {
                    const ControlPointType& cp = pcontrol_grid->operator[](i);
                    #ifdef USE_FUNCTION_ID_FOR_NODE_ID
                    ++NodeCounter;
                    typename NodeType::Pointer pNewNode = r_model_part.CreateNewNode(func_ids[i]+1, cp.X(), cp.Y(), cp.Z());
                    true_node_ids[NodeCounter] = pNewNode->Id();
                    #else
                    typename NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++NodeCounter, cp.X(), cp.Y(), cp.Z());
                    #endif

                    pNewNode->GetSolutionStepValue(BASIS_FUNCTION_INDEX) = func_ids[i] + 1;
                    // TODO transfer the nodal values

                    #ifdef DEBUG_MESH_GENERATION
                    std::cout << "Node " << pNewNode->Id() << " is added to model_part " << r_model_part.Name() << std::endl;
                    #endif
                }

                #ifdef USE_FUNCTION_ID_FOR_NODE_ID
                // obtain the raw connectivity
                std::vector<std::vector<IndexType> > tmp_connectivities;
                pcontrol_grid->CreateConnectivity(offset, tmp_connectivities);

                // modify the connectivity
                connectivities = tmp_connectivities;
                for (std::size_t i = 0; i < tmp_connectivities.size(); ++i)
                {
                    for (std::size_t j = 0; j < tmp_connectivities[i].size(); ++j)
                    {
                        connectivities[i][j] = true_node_ids[tmp_connectivities[i][j]];
                    }
                }
                #else
                pcontrol_grid->CreateConnectivity(offset, connectivities);
                #endif

                #ifdef DEBUG_MESH_GENERATION
                for (std::size_t i = 0; i < connectivities.size(); ++i)
                {
                    std::cout << "connectivities[" << i << "]:";
                    for (std::size_t j = 0; j < connectivities[i].size(); ++j)
                    {
                        std::cout << " " << connectivities[i][j];
                    }
                    std::cout << std::endl;
                }
                #endif
            }
            else
            {
                KRATOS_ERROR << "FESpace of type " << it->pFESpace()->Type() << " is not yet supported";
            }

            // create elements
            ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<ModelPart, std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                                                 connectivities, r_model_part, rCloneElement, ElementCounter, pNewProperties, NodeKey);

            for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
            {
                r_model_part.AddElement(*it2);
                if (mEchoLevel > 2)
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
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultipatchLagrangeControlMesh<" << TDim << ">";
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

    int mEchoLevel;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultipatchLagrangeControlMesh<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_MESH_GENERATION
#undef USE_FUNCTION_ID_FOR_NODE_ID

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_LAGRANGE_CONTROL_MESH_H_INCLUDED defined
