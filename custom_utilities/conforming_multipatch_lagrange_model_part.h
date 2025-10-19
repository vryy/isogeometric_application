//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 30 Oct 2024 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONFORMING_MULTIPATCH_LAGRANGE_MODEL_PART_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONFORMING_MULTIPATCH_LAGRANGE_MODEL_PART_H_INCLUDED

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
#include "layer_application/custom_utilities/auto_collapse_spatial_binning.h"

// #define DEBUG_MESH_GENERATION

namespace Kratos
{

/**
 * Construct the standard FEM mesh based on Lagrange basis functions from isogeometric multipatch.
 * Each patch can have different divisio. The total mesh is collapsed to conform at the boundary.
 * The principle is that each patch can be sampled differently.
 * This class is useful for pre-processing, i.e., constructing finite element mesh, of all types of
 * isogeometric patches, including NURBS, hierarchical B-Splines and T-Splines.
 */
template<class TMultiPatchType, class TModelPartType>
class ConformingMultipatchLagrangeModelPart : public IsogeometricEcho
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConformingMultipatchLagrangeModelPart);

    /// Constants
    static constexpr int Dim = TMultiPatchType::Dim;

    /// Type definition
    typedef TMultiPatchType MultiPatchType;
    typedef typename MultiPatchType::PatchType PatchType;
    typedef typename PatchType::BoundaryPatchType BoundaryPatchType;
    typedef typename PatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename PatchType::ControlPointType ControlPointType;

    typedef TModelPartType ModelPartType;
    typedef typename ModelPartType::IndexType IndexType;
    typedef typename ModelPartType::NodeType NodeType;
    typedef typename ModelPartType::ElementType ElementType;
    typedef typename ModelPartType::ConditionType ConditionType;
    typedef typename ModelPartType::NodesContainerType NodesContainerType;
    typedef typename ModelPartType::ElementsContainerType ElementsContainerType;
    typedef typename ModelPartType::ConditionsContainerType ConditionsContainerType;

    typedef AutoCollapseSpatialBinning<IndexType, double> BinningType;
    typedef std::map<IndexType, std::vector<std::vector<IndexType> > > connectivity_t;

    /// Default constructor
    ConformingMultipatchLagrangeModelPart(typename TMultiPatchType::Pointer pMultiPatch)
        : mpMultiPatch(pMultiPatch), mIsModelPartReady(false), mLastNodeId(0)
        , mName("MultiPatch")
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        mpModel = Model::Pointer(new Model());
#else
        mpModelPart = typename ModelPartType::Pointer(new ModelPartType(mName));
#endif
        mpBinning = typename BinningType::Pointer(new BinningType(0.0, 0.0, 0.0, 1e-3, 1e-3, 1e-3, 1e-10));
    }

    /// Destructor
    virtual ~ConformingMultipatchLagrangeModelPart() {}

    /// Get the underlying model_part
#ifdef SD_APP_FORWARD_COMPATIBILITY
    ModelPartType& GetModelPart() {return static_cast<ModelPartType&>(mpModel->GetModelPart(mName));}
#else
    ModelPartType& GetModelPart() {return static_cast<ModelPartType&>(*mpModelPart);}
#endif

    /// Get the underlying model_part
#ifdef SD_APP_FORWARD_COMPATIBILITY
    const ModelPartType& GetModelPart() const {return static_cast<const ModelPartType&>(mpModel->GetModelPart(mName));}
#else
    const ModelPartType& GetModelPart() const {return *mpModelPart;}
#endif

    /// Access the underlying binning
    typename BinningType::Pointer pGetBinning() const
    {
        return mpBinning;
    }

    /// Set the uniform sampling for the patch at specific dimension
    void SetDivision(IndexType patch_id, int dim, const IndexType nsampling)
    {
        if (dim >= Dim)
        {
            KRATOS_ERROR << "The dimension " << dim << " is invalid";
        }

        if (mpMultiPatch->Patches().find(patch_id) == mpMultiPatch->end())
        {
            KRATOS_ERROR << "Patch " << patch_id << " is not found in the multipatch";
        }

        std::vector<double> sampling;
        for (std::size_t i = 0; i <= nsampling; ++i)
            sampling.push_back((double) i / nsampling);

        mPatchSampling[patch_id][dim] = sampling;
    }

    /// Set the sampling for the patch at specific dimension
    void SetSampling(IndexType patch_id, int dim, const std::vector<double>& sampling)
    {
        if (dim >= Dim)
        {
            KRATOS_ERROR << "The dimension " << dim << " is invalid";
        }

        if (mpMultiPatch->Patches().find(patch_id) == mpMultiPatch->end())
        {
            KRATOS_ERROR << "Patch " << patch_id << " is not found in the multipatch";
        }

        mPatchSampling[patch_id][dim] = sampling;
    }

    /// Get the underlying multipatch pointer
    typename TMultiPatchType::Pointer pMultiPatch() {return mpMultiPatch;}

    /// Get the underlying multipatch pointer
    typename TMultiPatchType::ConstPointer pMultiPatch() const {return mpMultiPatch;}

    /// Check if the multipatch model_part ready for transferring/transmitting data
    bool IsReady() const {return mpMultiPatch->IsEnumerated() && mIsModelPartReady;}

    /// Assign the appropriate starting node and element index
    void BeginModelPart()
    {
        mIsModelPartReady = false;
        mLastNodeId = 0;
        mName = "MultiPatch";

#ifdef SD_APP_FORWARD_COMPATIBILITY
        mpModel->DeleteModelPart(mName);
        mpModel->CreateModelPart(mName);
#else
        // create new model_part
        typename ModelPartType::Pointer pNewModelPart = typename ModelPartType::Pointer(new ModelPartType(mpModelPart->Name()));

        // swap the internal model_part with new model_part
        mpModelPart.swap(pNewModelPart);
#endif

        // initialize the mesh
        this->Initialize(*mpBinning, mconnectivities);
    }

    /// Start the process to fill existing model_part. This function will first create the new model_part instance and add in the nodes (which are the control points in the multipatch)
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    void BeginModelPart(const std::string& Name)
    {
        mIsModelPartReady = false;
        mLastNodeId = MultiPatchUtility::GetLastNodeId(mpModel->GetModelPart(Name));

        // store the model_part name
        mName = Name;

        // initialize the mesh
        this->Initialize(*mpBinning, mconnectivities);
    }
    #else
    void BeginModelPart(typename ModelPartType::Pointer pModelPart)
    {
        mIsModelPartReady = false;
        mLastNodeId = pModelPart->GetLastNodeId();

        // store the model_part
        mpModelPart = pModelPart;

        // initialize the mesh
        this->Initialize(*mpBinning, mconnectivities);
    }
    #endif

    /// create the nodes from the points and add to the model_part
    NodesContainerType CreateNodes()
    {
        NodesContainerType pNodes;

        // create new nodes from points
        for (std::size_t i = 0; i < mpBinning->NumberOfNodes(); ++i)
        {
            auto pNewNode = this->GetModelPart().CreateNewNode(i + 1 + mLastNodeId, mpBinning->GetX(i+1), mpBinning->GetY(i+1), mpBinning->GetZ(i+1));
            pNodes.push_back(pNewNode);
        }

        if (this->GetEchoLevel() > 0)
        {
            std::cout << Info() << "::" << __FUNCTION__ << " completed"
                      << ", " << pNodes.size() << " new nodes are added to the model_part " << this->GetModelPart().Name()
                      << std::endl;
        }

        return pNodes;
    }

    /// create the elements out from the patch and add to the model_part
    ElementsContainerType AddElements(typename PatchType::Pointer pPatch, const std::string& element_name,
            std::size_t starting_id, Properties::Pointer pProperties)
    {
        const auto it = mconnectivities.find(pPatch->Id());
        if (it == mconnectivities.end())
            return ElementsContainerType();

        const std::string NodeKey = std::string("Node");

        if (!KratosComponents<Element>::Has(element_name))
        {
            KRATOS_ERROR << "Element " << element_name << " is not registered in Kratos."
                         << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
        }

        ElementType const& rCloneElement = KratosComponents<ElementType>::Get(element_name);

        // create elements
        ElementsContainerType pNewElements = IsogeometricPostUtility::CreateEntities<ModelPartType, std::vector<std::vector<IndexType> >, ElementType, ElementsContainerType, 1>(
                                             it->second, GetModelPart(), rCloneElement, starting_id, pProperties, NodeKey);

        for (auto it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
        {
            this->GetModelPart().AddElement(*it2);
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
        this->GetModelPart().Elements().Unique();

        return pNewElements;
    }

    /// create the conditions out from the boundary patch and add to the model_part
    ConditionsContainerType AddConditions(typename PatchType::Pointer pPatch, const BoundarySide side,
            const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
    {
        // construct the boundary patch
        const auto pBoundaryPatch = pPatch->ConstructBoundaryPatch(side);

        // get the sampling numbers
        const auto it_num = mPatchSampling.find(pPatch->Id());
        if (it_num == mPatchSampling.end())
        {
            KRATOS_ERROR << "Number of sampling is not set for patch " << pPatch->Id();
        }

        boost::array<std::vector<double>, Dim-1> div;
        if constexpr (Dim == 2)
        {
            const auto pdir = ParameterDirection<Dim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            if (GetEchoLevel() > 1)
            {
                std::cout << "Divisioning for patch " << pPatch->Id() << ":" << std::endl;
                KRATOS_WATCH_STD_CON(div[0])
            }
        }
        else if constexpr (Dim == 3)
        {
            const auto pdir = ParameterDirection<Dim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            div[1] = it_num->second[pdir[1]];
            if (GetEchoLevel() > 1)
            {
                std::cout << "Divisioning for patch " << pPatch->Id() << ":" << std::endl;
                KRATOS_WATCH_STD_CON(div[0])
                KRATOS_WATCH_STD_CON(div[1])
            }
        }

        // create conditions
        return this->AddConditions(pBoundaryPatch, div, condition_name, starting_id, pProperties);
    }

    /// create the conditions out from the slice patch and add to the model_part
    ConditionsContainerType AddConditions(typename PatchType::Pointer pPatch, const int idir, const double xi,
            const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
    {
        // extract the slice patch
        const auto pSlicedPatch = pPatch->ConstructSlicedPatch(idir, xi);

        // get the sampling numbers
        const auto it_num = mPatchSampling.find(pPatch->Id());
        if (it_num == mPatchSampling.end())
        {
            KRATOS_ERROR << "Number of sampling is not set for patch " << pPatch->Id();
        }

        boost::array<std::vector<double>, Dim-1> div;
        const BoundarySide side = ParameterDirection<Dim>::GetSide(idir);
        if constexpr (Dim == 2)
        {
            const auto pdir = ParameterDirection<Dim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            if (GetEchoLevel() > 1)
            {
                std::cout << "Divisioning for patch " << pPatch->Id() << ":" << std::endl;
                KRATOS_WATCH_STD_CON(div[0])
            }
        }
        else if constexpr (Dim == 3)
        {
            const auto pdir = ParameterDirection<Dim>::Get(side);
            div[0] = it_num->second[pdir[0]];
            div[1] = it_num->second[pdir[1]];
            if (GetEchoLevel() > 1)
            {
                std::cout << "Divisioning for patch " << pPatch->Id() << ":" << std::endl;
                KRATOS_WATCH_STD_CON(div[0])
                KRATOS_WATCH_STD_CON(div[1])
            }
        }

        // create conditions
        return this->AddConditions(pSlicedPatch, div, condition_name, starting_id, pProperties);
    }

    /// create the conditions out from the boudnary patch and add to the model_part
    ConditionsContainerType AddConditions(typename BoundaryPatchType::Pointer pBoundaryPatch,
            const boost::array<std::vector<double>, Dim-1>& nsampling,
            const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
    {
        if (!KratosComponents<Condition>::Has(condition_name))
        {
            KRATOS_ERROR << "Condition " << condition_name << " is not registered in Kratos."
                         << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
        }

        ConditionType const& rCloneCondition = KratosComponents<ConditionType>::Get(condition_name);

        std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;
        IndexType NodeCounter = 1;

        // generate points and connectivities
        if constexpr (Dim == 2)
        {
            std::vector<array_1d<double, 3> > corners(2);

            IsogeometricPostUtility::GenerateLine(corners, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateLineGrid(corners[0], corners[1], NodeCounter, nsampling[0]);
        }
        else if constexpr (Dim == 3)
        {
            std::vector<array_1d<double, 3> > corners(4);

            IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

            points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                                                corners[2], corners[3], NodeCounter, nsampling[0], nsampling[1]);
        }

        // add the points to the binning
        std::map<IndexType, IndexType> old_to_new;
        IndexType cnt = 0;
        if constexpr ((Dim == 2) || (Dim == 3))
        {
            for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
            {
                const auto& p = points_and_connectivities.first[i];

                const auto cp = IsogeometricPostUtility::CreatePoint(p, *pBoundaryPatch);

                ++cnt;
                const IndexType id = mpBinning->GetNode(cp.X(), cp.Y(), cp.Z());
                old_to_new[cnt] = id + mLastNodeId;
            }
        }

        // create new conditions
        ConditionsContainerType pNewConditions;
        const std::string NodeKey = std::string("Node");
        for (std::size_t i = 0; i < points_and_connectivities.second.size(); ++i)
        {
            const auto& c = points_and_connectivities.second[i];

            std::vector<IndexType> newc;
            for (std::size_t j = 0; j < c.size(); ++j)
            {
                newc.push_back(old_to_new[c[j]]);
            }

            // create condition
            typename ConditionType::Pointer pNewCondition = IsogeometricPostUtility::CreateEntity<ModelPartType, std::vector<IndexType>, ConditionType>(
                                                        newc, GetModelPart(), rCloneCondition, i + starting_id, pProperties, NodeKey);
            this->GetModelPart().AddCondition(pNewCondition);
            pNewConditions.push_back(pNewCondition);
        }

        // just to make sure everything is organized properly
        this->GetModelPart().Conditions().Unique();

        return pNewConditions;
    }

    /// Finalize the model_part creation
    void EndModelPart()
    {
        if (IsReady()) { return; }
        mIsModelPartReady = true;
    }

    std::string Info() const
    {
        std::stringstream ss;
        ss << "ConformingMultipatchLagrangeModelPart<" << Dim << ">";
        return ss.str();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConformingMultipatchLagrangeModelPart<" << Dim << ">";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    void Initialize(BinningType& rBinning, connectivity_t& connectivities) const
    {
        IndexType NodeCounter;

        typedef typename TMultiPatchType::patch_const_iterator patch_const_iterator;
        for (patch_const_iterator it = mpMultiPatch->begin(); it != mpMultiPatch->end(); ++it)
        {
            if (!it->Is(ACTIVE))
            {
                continue;
            }

            // generate the connectivities
            std::pair<std::vector<array_1d<double, 3> >, std::vector<std::vector<IndexType> > > points_and_connectivities;

            NodeCounter = 1;
            if constexpr (Dim == 2)
            {
                // create new nodes and elements
                auto it_num = mPatchSampling.find(it->Id());
                if (it_num == mPatchSampling.end())
                {
                    KRATOS_ERROR << "Number of sampling is not set for patch " << it->Id();
                }

                const std::vector<double>& div_1 = it_num->second[0];
                const std::vector<double>& div_2 = it_num->second[1];
                if (GetEchoLevel() > 1)
                {
                    std::cout << "Divisioning for patch " << it->Id() << ":" << std::endl;
                    KRATOS_WATCH_STD_CON(div_1)
                    KRATOS_WATCH_STD_CON(div_2)
                }

                std::vector<array_1d<double, 3> > corners(4);

                IsogeometricPostUtility::GenerateRectangle(corners, 0.0, 1.0, 0.0, 1.0);

                points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(corners[0], corners[1],
                                            corners[2], corners[3], NodeCounter, div_1, div_2);
            }
            else if constexpr (Dim == 3)
            {
                // create new nodes and elements
                auto it_num = mPatchSampling.find(it->Id());
                if (it_num == mPatchSampling.end())
                {
                    KRATOS_ERROR << "Number of sampling is not set for patch " << it->Id();
                }

                const std::vector<double>& div_1 = it_num->second[0];
                const std::vector<double>& div_2 = it_num->second[1];
                const std::vector<double>& div_3 = it_num->second[2];
                if (GetEchoLevel() > 1)
                {
                    std::cout << "Divisioning for patch " << it->Id() << ":" << std::endl;
                    KRATOS_WATCH_STD_CON(div_1)
                    KRATOS_WATCH_STD_CON(div_2)
                    KRATOS_WATCH_STD_CON(div_3)
                }

                std::vector<array_1d<double, 3> > corners(8);

                IsogeometricPostUtility::GenerateBox(corners, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

                points_and_connectivities = IsogeometricPostUtility::GenerateHexGrid(corners[0], corners[1], corners[2], corners[3],
                                            corners[4], corners[5], corners[6], corners[7], NodeCounter, div_1, div_2, div_3);
            }

            // if (this->GetEchoLevel() > 1)
            // {
            //     for (std::size_t i = 0; i < points_and_connectivities.second.size(); ++i)
            //     {
            //         const auto& c = points_and_connectivities.second[i];
            //         std::cout << "conn[" << i << "]:";
            //         for (std::size_t j = 0; j < c.size(); ++j)
            //         {
            //             std::cout << " " << c[j];
            //         }
            //         std::cout << std::endl;
            //     }
            // }

            // add points to the binning
            std::map<IndexType, IndexType> old_to_new;
            IndexType cnt = 0;
            for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
            {
                const auto& p = points_and_connectivities.first[i];

                const auto cp = IsogeometricPostUtility::CreatePoint(p, *it);

                ++cnt;
                const IndexType id = rBinning.AddNode(cp.X(), cp.Y(), cp.Z());
                old_to_new[cnt] = id + mLastNodeId;
            }

            // add the connectivities
            auto& conn = connectivities[it->Id()];
            for (std::size_t i = 0; i < points_and_connectivities.second.size(); ++i)
            {
                const auto& c = points_and_connectivities.second[i];

                std::vector<IndexType> newc;
                for (std::size_t j = 0; j < c.size(); ++j)
                {
                    newc.push_back(old_to_new[c[j]]);
                }
                conn.push_back(newc);
            }
        }

        if (this->GetEchoLevel() > 1)
        {
            for (auto it = connectivities.begin(); it != connectivities.end(); ++it)
            {
                for (std::size_t i = 0; i < it->second.size(); ++i)
                {
                    const auto& c = it->second[i];
                    std::cout << "conn[" << i << "]:";
                    for (std::size_t j = 0; j < c.size(); ++j)
                    {
                        std::cout << " " << c[j];
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    bool mIsModelPartReady;
    std::string mName;
    IndexType mLastNodeId;

#ifdef SD_APP_FORWARD_COMPATIBILITY
    Model::Pointer mpModel;
#else
    typename ModelPartType::Pointer mpModelPart;
#endif
    typename TMultiPatchType::Pointer mpMultiPatch;

    std::map<IndexType, boost::array<std::vector<double>, Dim> > mPatchSampling;

    typename BinningType::Pointer mpBinning;
    connectivity_t mconnectivities;
};

/// output stream function
template<class TMultiPatchType, class TModelPartType>
inline std::ostream& operator <<(std::ostream& rOStream, const ConformingMultipatchLagrangeModelPart<TMultiPatchType, TModelPartType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_MESH_GENERATION

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONFORMING_MULTIPATCH_LAGRANGE_MODEL_PART_H_INCLUDED defined
