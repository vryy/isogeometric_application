//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Dec 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_MULTIPATCH_MODEL_PART_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_MULTIPATCH_MODEL_PART_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
#include "includes/model_part.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "containers/model.h"
#endif
#include "utilities/openmp_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/multipatch_model_part.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "isogeometric_application_variables.h"

#define ENABLE_PROFILING

namespace Kratos
{

/**
Coupling between KRATOS model_part and multiple multipatch structure. THis is useful for simulation involving mixed elements.
All the multipatch must have the same underlying space.
 */
template<class TMultiPatchType, class TModelPartType>
class MultiMultiPatchModelPart : public IsogeometricEcho
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiMultiPatchModelPart);

    /// Constants
    static constexpr int Dim = TMultiPatchType::Dim;

    /// Type definition
    typedef TMultiPatchType MultiPatchType;
    typedef typename MultiPatchType::PatchType PatchType;
    typedef typename PatchType::BoundaryPatchType BoundaryPatchType;
    typedef typename PatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename PatchType::ControlPointType ControlPointType;

    typedef FESpace<Dim, LocalCoordinateType> FESpaceType;
    typedef FESpace<Dim-1, LocalCoordinateType> BoundaryFESpaceType;

    typedef TModelPartType ModelPartType;
    typedef typename ModelPartType::NodeType NodeType;
    typedef typename ModelPartType::ElementType ElementType;
    typedef typename ModelPartType::ConditionType ConditionType;
    typedef typename ModelPartType::NodesContainerType NodesContainerType;
    typedef typename ModelPartType::ConditionsContainerType ConditionsContainerType;
    typedef typename ModelPartType::ElementsContainerType ElementsContainerType;

    typedef IsogeometricGeometry<NodeType> IsogeometricGeometryType;

    /// Default constructor
    MultiMultiPatchModelPart() : mIsModelPartReady(false), mNodeOffset(0)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        mpModel = Model::Pointer(new Model());
#else
        mpModelPart = typename ModelPartType::Pointer(new ModelPartType("MultiMultiPatch"));
#endif
    }

    /// Destructor
    virtual ~MultiMultiPatchModelPart() {}

    /// Get the underlying model_part
#ifdef SD_APP_FORWARD_COMPATIBILITY
    ModelPartType& GetModelPart() {return static_cast<ModelPartType&>(mpModel->GetModelPart("MultiMultiPatch"));}
#else
    ModelPartType& GetModelPart() {return *mpModelPart;}
#endif

    /// Get the underlying model_part
#ifdef SD_APP_FORWARD_COMPATIBILITY
    const ModelPartType& GetModelPart() const {return static_cast<const ModelPartType&>(mpModel->GetModelPart("MultiMultiPatch"));}
#else
    const ModelPartType& GetModelPart() const {return *mpModelPart;}
#endif

    /// Add the multipatch to the list
    void AddMultiPatch(typename MultiPatchType::Pointer pMultiPatch)
    {
        mpMultiPatches.push_back(pMultiPatch);
    }

    /// Get the underlying multipatch pointer
    typename MultiPatchType::Pointer pMultiPatch(std::size_t ip) {return mpMultiPatches[ip];}

    /// Get the underlying multipatch pointer
    typename MultiPatchType::ConstPointer pMultiPatch(std::size_t ip) const {return mpMultiPatches[ip];}

    /// Check if the multipatch model_part ready for transferring/transmitting data
    bool IsReady() const
    {
        bool is_ready = mIsModelPartReady;
        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            is_ready = is_ready && mpMultiPatches[ip]->IsEnumerated();
        }
        return is_ready;
    }

    /// Start the process to create and fill new model_part. This function will first create the new model_part instance and add in the nodes (which are the control points in the multipatch)
    void BeginModelPart()
    {
        mIsModelPartReady = false;
        mNodeOffset = 0;

        // always enumerate the multipatch first
        std::size_t EquationSystemSize = 0;
        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            EquationSystemSize = mpMultiPatches[ip]->Enumerate(EquationSystemSize);
            // KRATOS_WATCH(EquationSystemSize)
        }

#ifdef SD_APP_FORWARD_COMPATIBILITY
        mpModel->DeleteModelPart("MultiMultiPatch");
        mpModel->CreateModelPart("MultiMultiPatch");
#else
        // create new model_part
        typename ModelPartType::Pointer pNewModelPart = typename ModelPartType::Pointer(new ModelPartType(mpModelPart->Name()));

        // swap the internal model_part with new model_part
        mpModelPart.swap(pNewModelPart);
#endif
    }

    /// Start the process to fill existing model_part. This function will first create the new model_part instance and add in the nodes (which are the control points in the multipatch)
    void BeginModelPart(typename ModelPartType::Pointer pModelPart)
    {
        mIsModelPartReady = false;
        mNodeOffset = pModelPart->GetLastNodeId() + 1;

        // always enumerate the multipatch first
        std::size_t EquationSystemSize = 0;
        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            EquationSystemSize = mpMultiPatches[ip]->Enumerate(EquationSystemSize);
            // KRATOS_WATCH(EquationSystemSize)
        }

        // store the model_part
        mpModelPart = pModelPart;
    }

    /// create the nodes from the control points and add to the model_part
    void CreateNodes()
    {
#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        std::size_t node_counter = 0;
        std::size_t cnt = 0;

        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            if (!mpMultiPatches[ip]->IsEnumerated())
            {
                KRATOS_ERROR << "The multipatch " << ip << " is not yet enumerated";
            }

            // create new nodes from control points
            for (std::size_t i = 0; i < mpMultiPatches[ip]->EquationSystemSize(); ++i)
            {
                std::tuple<std::size_t, std::size_t> loc = mpMultiPatches[ip]->EquationIdLocation(cnt++);

                std::size_t patch_id = std::get<0>(loc);
                std::size_t local_id = std::get<1>(loc);
                // KRATOS_WATCH(patch_id)
                // KRATOS_WATCH(local_id)

                const ControlPointType& point = mpMultiPatches[ip]->pGetPatch(patch_id)->pControlPointGridFunction()->pControlGrid()->GetData(local_id);
                // KRATOS_WATCH(point)

                auto pNewNode = this->GetModelPart().CreateNewNode(CONVERT_INDEX_IGA_TO_KRATOS(node_counter) + mNodeOffset, point.X(), point.Y(), point.Z());
                ++node_counter;
            }
        }

        if (this->GetEchoLevel() > 0)
        {
#ifdef ENABLE_PROFILING
            std::cout << "+++ " << __FUNCTION__ << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
#else
            std::cout << __FUNCTION__ << " completed" << std::endl;
#endif
        }
    }

    /// create the elements out from the patches and add to the model_part
    ElementsContainerType AddElements(std::vector<typename PatchType::Pointer> pPatches,
            const std::string& element_name,
            std::size_t starting_id, Properties::Pointer pProperties)
    {
        if (IsReady()) { return ElementsContainerType(); } // call BeginModelPart first before adding elements

#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        // get the list of FESpaces and control grids
        std::vector<typename FESpaceType::ConstPointer> pFESpaces;
        std::vector<typename ControlGrid<ControlPointType>::ConstPointer> pControlGrids;

        for (std::size_t i = 0; i < pPatches.size(); ++i)
        {
            pFESpaces.push_back(pPatches[i]->pFESpace());
            const auto& rControlPointGridFunction = pPatches[i]->ControlPointGridFunction();
            pControlGrids.push_back(rControlPointGridFunction.pControlGrid());
        }

        // create new elements and add to the model_part
        ElementsContainerType pNewElements = CreateEntitiesFromFESpace<ElementType, FESpaceType,
                ControlGrid<ControlPointType>, NodesContainerType>(
            pFESpaces, pControlGrids, this->GetModelPart().Nodes(),
            element_name, starting_id, mNodeOffset,
            pProperties, this->GetEchoLevel());

        for (auto it = pNewElements.ptr_begin(); it != pNewElements.ptr_end(); ++it)
        {
            this->GetModelPart().Elements().push_back(*it);
        }

        // sort the element container and make it consistent
        this->GetModelPart().Elements().Unique();

        if (this->GetEchoLevel() > 0)
        {
#ifdef ENABLE_PROFILING
            std::cout << "+++ " << __FUNCTION__ << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s, " << pNewElements.size() << " elements of type " << element_name << " are generated" << std::endl;
#else
            std::cout << __FUNCTION__ << " completed" << std::endl;
#endif
        }

        return pNewElements;
    }

    /// create the elements out from the patches and add to the model_part
    ConditionsContainerType AddConditions(std::vector<typename PatchType::Pointer> pPatches,
            const std::string& condition_name,
            std::size_t starting_id, Properties::Pointer pProperties)
    {
        if (IsReady()) { return ConditionsContainerType(); } // call BeginModelPart first before adding elements

#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        // get the list of FESpaces and control grids
        std::vector<typename FESpaceType::ConstPointer> pFESpaces;
        std::vector<typename ControlGrid<ControlPointType>::ConstPointer> pControlGrids;

        for (std::size_t i = 0; i < pPatches.size(); ++i)
        {
            pFESpaces.push_back(pPatches[i]->pFESpace());
            const auto& rControlPointGridFunction = pPatches[i]->ControlPointGridFunction();
            pControlGrids.push_back(rControlPointGridFunction.pControlGrid());
        }

        // create new elements and add to the model_part
        ConditionsContainerType pNewConditions = CreateEntitiesFromFESpace<ConditionType, FESpaceType,
                ControlGrid<ControlPointType>, NodesContainerType>(
            pFESpaces, pControlGrids, this->GetModelPart().Nodes(),
            condition_name, starting_id, mNodeOffset,
            pProperties, this->GetEchoLevel());

        for (auto it = pNewConditions.ptr_begin(); it != pNewConditions.ptr_end(); ++it)
        {
            this->GetModelPart().Conditions().push_back(*it);
        }

        // sort the element container and make it consistent
        this->GetModelPart().Conditions().Unique();

        if (this->GetEchoLevel() > 0)
        {
#ifdef ENABLE_PROFILING
            std::cout << "+++ " << __FUNCTION__ << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s, " << pNewConditions.size() << " conditions of type " << condition_name << " are generated" << std::endl;
#else
            std::cout << __FUNCTION__ << " completed" << std::endl;
#endif
        }

        return pNewConditions;
    }

    /// create the conditions out from the boundary of the patch and add to the model_part
    ConditionsContainerType AddConditions(typename PatchType::Pointer pPatch, const BoundarySide& side,
            const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
    {
        if (IsReady()) { return ConditionsContainerType(); } // call BeginModelPart first before adding conditions

        // construct the boundary patch
        typename BoundaryPatchType::Pointer pBoundaryPatch = pPatch->ConstructBoundaryPatch(side);

        return AddConditions(pBoundaryPatch, condition_name, starting_id, pProperties);
    }

    /// create the conditions out from the boundary of the patch and add to the model_part
    ConditionsContainerType AddConditions(typename BoundaryPatchType::Pointer pBoundaryPatch,
            const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
    {
        if (IsReady()) { return ConditionsContainerType(); } // call BeginModelPart first before adding conditions

#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        // KRATOS_WATCH(*pBoundaryPatch)

        // get the grid function for control points
        const auto& rControlPointGridFunction = pBoundaryPatch->ControlPointGridFunction();

        // create new conditions and add to the model_part
        ConditionsContainerType pNewConditions = MultiPatchModelPart<TMultiPatchType, ModelPartType>::template CreateEntitiesFromFESpace<ConditionType, BoundaryFESpaceType,
                ControlGrid<ControlPointType>, NodesContainerType>(
            pBoundaryPatch->pFESpace(), rControlPointGridFunction.pControlGrid(),
            this->GetModelPart().Nodes(), condition_name,
            starting_id, mNodeOffset,
            pProperties, this->GetEchoLevel());

        for (auto it = pNewConditions.ptr_begin(); it != pNewConditions.ptr_end(); ++it)
        {
            this->GetModelPart().Conditions().push_back(*it);
        }

        // sort the condition container and make it consistent
        this->GetModelPart().Conditions().Unique();

        if (this->GetEchoLevel() > 0)
        {
#ifdef ENABLE_PROFILING
            std::cout << "+++ " << __FUNCTION__ << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s, " << pNewConditions.size() << " conditions of type " << condition_name << " are generated" << std::endl;
#else
            std::cout << __FUNCTION__ << " completed" << std::endl;
#endif
        }

        return pNewConditions;
    }

    /// Finalize the model_part creation process
    void EndModelPart()
    {
        if (IsReady()) { return; }
        mIsModelPartReady = true;
    }

    /// Synchronize from multipatch to model_part
    template<class TVariableType>
    void SynchronizeForward(const TVariableType& rVariable)
    {
        if (!IsReady()) { return; }

        // transfer data from from control points to nodes
        std::size_t cnt = 0;
        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            if (!mpMultiPatches[ip]->IsEnumerated())
            {
                KRATOS_ERROR << "The multipatch " << ip << " is not yet enumerated";
            }

            for (std::size_t i = 0; i < mpMultiPatches[ip]->EquationSystemSize(); ++i)
            {
                std::tuple<std::size_t, std::size_t> loc = mpMultiPatches[ip]->EquationIdLocation(i);

                std::size_t patch_id = std::get<0>(loc);
                std::size_t local_id = std::get<1>(loc);
                // KRATOS_WATCH(patch_id)
                // KRATOS_WATCH(local_id)

                const typename TVariableType::Type& value = mpMultiPatches[ip]->pGetPatch(patch_id)->pGetGridFunction(rVariable)->pControlGrid()->GetData(local_id);
                // KRATOS_WATCH(value)

                auto pNode = this->GetModelPart().pGetNode(CONVERT_INDEX_IGA_TO_KRATOS(cnt++) + mNodeOffset);

                pNode->GetSolutionStepValue(rVariable) = value;
            }
        }
    }

    /// Synchronize from model_part to the multipatch
    template<class TVariableType>
    void SynchronizeBackward(std::size_t ip, const TVariableType& rVariable)
    {
        if (!IsReady()) { return; }

        // loop through each patch, we construct a map from each function id to the patch id
        typedef typename MultiPatchType::patch_iterator patch_iterator;
        for (patch_iterator it = mpMultiPatches[ip]->begin(); it != mpMultiPatches[ip]->end(); ++it)
        {
            std::vector<std::size_t> func_ids = it->pFESpace()->FunctionIndices();

            // check if the grid function existed in the patch
            if (!it->template HasGridFunction<TVariableType>(rVariable))
            {
                // if not then create the new grid function
                typename ControlGrid<typename TVariableType::Type>::Pointer pNewControlGrid = UnstructuredControlGrid<typename TVariableType::Type>::Create(it->pFESpace()->TotalNumber());
                it->template CreateGridFunction<TVariableType>(rVariable, pNewControlGrid);
            }

            // get the control grid
            typename ControlGrid<typename TVariableType::Type>::Pointer pControlGrid = it->pGetGridFunction(rVariable)->pControlGrid();

            // set the data for the control grid
            for (std::size_t i = 0; i < pControlGrid->size(); ++i)
            {
                std::size_t global_id = func_ids[i];
                std::size_t node_id = CONVERT_INDEX_IGA_TO_KRATOS(global_id) + mNodeOffset;

                pControlGrid->SetData(i, this->GetModelPart().Nodes()[node_id].GetSolutionStepValue(rVariable));
            }
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiMultiPatchModelPart";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "+++ModelPart:" << std::endl;
        rOStream << this->GetModelPart() << std::endl;
        for (std::size_t ip = 0; ip < mpMultiPatches.size(); ++ip)
        {
            rOStream << "+++MultiPatch " << ip << std::endl;
            rOStream << *mpMultiPatches[ip] << std::endl;
        }
    }

private:

    bool mIsModelPartReady;
    std::size_t mNodeOffset;

#ifdef SD_APP_FORWARD_COMPATIBILITY
    Model::Pointer mpModel;
#else
    typename ModelPartType::Pointer mpModelPart;
#endif
    std::vector<typename MultiPatchType::Pointer> mpMultiPatches;

    /// Create entities (elements/conditions) from FESpaces
    /// @param pFESpaces the list of finite element space to provide the cell manager
    /// @param pControlGris control grids to provide control points
    /// @param rNodes model_part Nodes to look up for when creating elements
    /// @param element_name name of the sample element
    /// @param starting_id the first id of the newly created entities, from there the id is incremental
    /// @param p_temp_properties the Properties to create new entities
    template<class TEntityType, class TFESpace, class TControlGridType, class TNodeContainerType>
    static PointerVectorSet<TEntityType, IndexedObject> CreateEntitiesFromFESpace(
        std::vector<typename TFESpace::ConstPointer> pFESpaces,
        std::vector<typename TControlGridType::ConstPointer> pControlGrids,
        TNodeContainerType& rNodes, const std::string& element_name,
        const std::size_t starting_id, const std::size_t node_offset,
        Properties::Pointer p_temp_properties, int echo_level)
    {
#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        // construct the cell manager out from the FESpaces
        typedef typename TFESpace::cell_container_t cell_container_t;

        std::vector<typename cell_container_t::Pointer> pCellManagers;
        for (std::size_t ip = 0; ip < pFESpaces.size(); ++ip)
        {
            pCellManagers.push_back(pFESpaces[ip]->ConstructCellManager());
        }

        // TODO compare all cell managers to make sure they are matching
        for (std::size_t ip = 1; ip < pCellManagers.size(); ++ip)
        {
            if ((*pCellManagers[ip]) != (*pCellManagers[0]))
                KRATOS_ERROR << "The cell manager does not match at index " << ip;
            }

        if (echo_level > 0)
        {
#ifdef ENABLE_PROFILING
            std::cout << "  ++ ConstructCellManager: " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
            start = OpenMPUtils::GetCurrentTime();
#endif
        }

        // container for newly created elements
        PointerVectorSet<TEntityType, IndexedObject> pNewElements;

        // get the sample element
        if (!KratosComponents<TEntityType>::Has(element_name))
        {
            KRATOS_ERROR << "Entity (Element/Condition) " << element_name << " is not registered in Kratos.";
            return pNewElements;
        }

        TEntityType const& r_clone_element = KratosComponents<TEntityType>::Get(element_name);

        // loop through each cell in the space
        typename TEntityType::NodesArrayType temp_element_nodes;
        std::size_t cnt = starting_id;
        Vector dummy;
        int max_integration_method = 1;
        if (p_temp_properties->Has(NUM_IGA_INTEGRATION_METHOD))
        {
            max_integration_method = (*p_temp_properties)[NUM_IGA_INTEGRATION_METHOD];
        }

        std::size_t ic = 0; // this is to mark the location of the iterator
        for (typename cell_container_t::iterator it_dummy = pCellManagers[0]->begin(); it_dummy != pCellManagers[0]->end(); ++it_dummy)
        {
            std::vector<typename ElementType::GeometryType::Pointer> p_temp_geometries;

            // fill the vector of geometries
            for (std::size_t ip = 0; ip < pFESpaces.size(); ++ip)
            {
                typename cell_container_t::iterator it_cell = pCellManagers[ip]->begin();
                std::advance(it_cell, ic);
                typename cell_container_t::cell_t pcell = *it_cell;
                // KRATOS_WATCH(*pcell)

                // get new nodes
                temp_element_nodes.clear();

                const std::vector<std::size_t>& anchors = pcell->GetSupportedAnchors();
                Vector weights(anchors.size());
                for (std::size_t i = 0; i < anchors.size(); ++i)
                {
                    temp_element_nodes.push_back(( *(MultiPatchUtility::FindKey(rNodes, CONVERT_INDEX_IGA_TO_KRATOS(anchors[i]) + node_offset, "Node").base())));
                    weights[i] = pControlGrids[ip]->GetData(pFESpaces[ip]->LocalId(anchors[i])).W();
                }

                if (echo_level > 1)
                {
                    std::cout << "nodes:";
                    for (std::size_t i = 0; i < anchors.size(); ++i)
                    {
                        std::cout << " " << CONVERT_INDEX_IGA_TO_KRATOS(anchors[i]) + node_offset;
                    }
                    std::cout << std::endl;
                    KRATOS_WATCH(weights)
                    // KRATOS_WATCH(pcell->GetExtractionOperator())
                    KRATOS_WATCH(pcell->GetCompressedExtractionOperator())
                    KRATOS_WATCH(pFESpaces[ip]->Order(0))
                    KRATOS_WATCH(pFESpaces[ip]->Order(1))
                    KRATOS_WATCH(pFESpaces[ip]->Order(2))
                }

                // create the geometry
                typename IsogeometricGeometryType::Pointer p_temp_geometry
                    = iga::dynamic_pointer_cast<IsogeometricGeometryType>(r_clone_element.GetGeometry().Create(temp_element_nodes));
                if (p_temp_geometry == NULL)
                    KRATOS_ERROR << "The cast to IsogeometricGeometry is failed.";

                    p_temp_geometry->AssignGeometryData(dummy,
                                                        dummy,
                                                        dummy,
                                                        weights,
                                                        // pcell->GetExtractionOperator(),
                                                        pcell->GetCompressedExtractionOperator(),
                                                        static_cast<int>(pFESpaces[ip]->Order(0)),
                                                        static_cast<int>(pFESpaces[ip]->Order(1)),
                                                        static_cast<int>(pFESpaces[ip]->Order(2)),
                                                        max_integration_method);

                p_temp_geometries.push_back(p_temp_geometry);
            }

            ++ic;

            // create the element and add to the list
            typename TEntityType::Pointer pNewElement = r_clone_element.Create(cnt++, p_temp_geometries, p_temp_properties);
            pNewElement->SetValue(ACTIVATION_LEVEL, 0);
            pNewElement->SetValue(IS_INACTIVE, false);
            pNewElement->Set(ACTIVE, true);
            pNewElements.push_back(pNewElement);
        }

        if (echo_level > 0)
        {
            std::cout << "  ++ generate entities: " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
            start = OpenMPUtils::GetCurrentTime();
        }

        return pNewElements;
    }

};

/// output stream function
template<class TMultiPatchType, class TModelPartType>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiMultiPatchModelPart<TMultiPatchType, TModelPartType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_GEN_ENTITY
#undef ENABLE_PROFILING

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_MULTIPATCH_MODEL_PART_H_INCLUDED
