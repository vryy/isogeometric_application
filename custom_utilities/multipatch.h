//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED

#include "custom_utilities/patch.h"
#include "custom_utilities/bending_strip_patch.h"

namespace Kratos
{

/**
This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
 */
template<int TDim>
class MultiPatch : public boost::enable_shared_from_this<MultiPatch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);

    /// Type definition
    typedef Patch<TDim> PatchType;
    typedef PointerVectorSet<PatchType, IndexedObject> PatchContainerType;

    typedef typename Patch<TDim>::vertex_t vertex_t;
    typedef typename Patch<TDim>::edge_t edge_t;
    typedef typename Patch<TDim>::face_t face_t;
    typedef typename Patch<TDim>::volume_t volume_t;

    /// Default constructor
    MultiPatch() : mIsEnumerated(false) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPatch(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
        mIsEnumerated = false;
    }

    /// Reset Id for all the patches
    void ResetId()
    {
        std::size_t Id = 0;
        for (typename PatchContainerType::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->SetId(++Id);
        }
    }

    /// Check the enumeration flag
    const bool& IsEnumerated() const {return mIsEnumerated;}

    /// Get the equation system size
    std::size_t EquationSystemSize() const {return mEquationSystemSize;}

    /// Locate the patch of the global equation id and the corresponding local id to determine the control value
    std::tuple<std::size_t, std::size_t> EquationIdLocation(const std::size_t& global_id) const
    {
        if (!IsEnumerated())
            KRATOS_THROW_ERROR(std::logic_error, "The multipatch is not enumerated", "")

        std::map<std::size_t, std::size_t>::const_iterator it = mGlobalToPatch.find(global_id);
        if (it == mGlobalToPatch.end())
        {
            KRATOS_WATCH(global_id)
            KRATOS_WATCH(mEquationSystemSize)
            std::cout << "global_to_patch map:" << std::endl;
            for (std::map<std::size_t, std::size_t>::const_iterator it2 = mGlobalToPatch.begin(); it2 != mGlobalToPatch.end(); ++it2)
                std::cout << " " << it2->first << ": " << it2->second << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "The global id does not exist in the global_to_patch map.", "")
        }

        const std::size_t& patch_id = it->second;
        const std::size_t& local_id = pGetPatch(it->second)->pFESpace()->LocalId(global_id);

        return std::make_tuple(patch_id, local_id);
    }

    /// Validate the MultiPatch
    virtual bool Validate() const
    {
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            bool check = it->Validate();
            if (!check)
                return false;
        }

        return true;
    }

    /// Get the patch with specific Id
    typename PatchType::Pointer pGetPatch(const std::size_t& Id)
    {
        typename PatchContainerType::ptr_iterator it_patch = mpPatches.find(Id).base();
        if(it_patch == mpPatches.ptr_end())
        {
            std::stringstream ss;
            ss << "The patch " << Id << " does not exist in the multipatch";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
        return *it_patch;
    }

    /// Get the patch with specific Id
    typename PatchType::ConstPointer pGetPatch(const std::size_t& Id) const
    {
        typename PatchContainerType::ptr_const_iterator it_patch = mpPatches.find(Id).base();
        if(it_patch == mpPatches.ptr_end())
        {
            std::stringstream ss;
            ss << "The patch " << Id << " does not exist in the multipatch";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
        return typename PatchType::ConstPointer(*it_patch);
    }

    /// Access the underlying list of patches
    /// WARNING!!! be careful with this routine
    PatchContainerType& Patches() {return mpPatches;}

    /// Access the underlying list of patches
    const PatchContainerType& Patches() const {return mpPatches;}

    /// iterators
    typename PatchContainerType::iterator begin() {return mpPatches.begin();}
    typename PatchContainerType::const_iterator begin() const {return mpPatches.begin();}
    typename PatchContainerType::iterator end() {return mpPatches.end();}
    typename PatchContainerType::const_iterator end() const {return mpPatches.end();}

    /// Get the number of patches
    std::size_t size() const {return mpPatches.size();}

    /// Enumerate all the patches, starting at 0
    std::size_t Enumerate()
    {
        return this->Enumerate(0);
    }

    /// Enumerate all the patches, with the given starting id
    std::size_t Enumerate(const std::size_t& start)
    {
        // reset global ids for each patch
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->ResetFunctionIndices();
        }

        // enumerate each patch
        std::size_t last = start;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsBendingStrip() == false)
            {
                last = (*it)->pFESpace()->Enumerate(last);
                // KRATOS_WATCH(last)

                // transfer the enumeration to neighbor boundary
                for (int i = _LEFT_; i <= _BACK_; ++i)
                {
                    BoundarySide side = static_cast<BoundarySide>(i);

                    if ((*it)->pNeighbor(side) != NULL)
                    {
                        // find the side of the other neighbor
                        BoundarySide other_side = (*it)->pNeighbor(side)->FindBoundarySide(*it);

                        if (other_side == _NUMBER_OF_BOUNDARY_SIDE)
                            KRATOS_THROW_ERROR(std::logic_error, "No neighbor of the neighbor is the same as this. Error setting the neighbor.", "")

                        // check the boundary compatibility again
                        if (!(*it)->CheckBoundaryCompatibility(side, *((*it)->pNeighbor(side)), other_side))
                        {
                            KRATOS_WATCH(side)
                            KRATOS_WATCH(other_side)
                            KRATOS_THROW_ERROR(std::logic_error, "The boundary compatibility with the neighbor is not satisfied", "")
                        }
                        else
                        {
                            std::vector<std::size_t> func_indices = (*it)->pFESpace()->ExtractBoundaryFunctionIndices(side);
                            (*it)->pNeighbor(side)->pFESpace()->AssignBoundaryFunctionIndices(other_side, func_indices);
                        }
                    }
                }
            }
        }

        // check if a patch is a bending strip, then that patch must be enumerated again using the enumeration info from the parent patches
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsBendingStrip() == true)
            {
                typename BendingStripPatch<TDim>::Pointer pBendPatch = boost::dynamic_pointer_cast<BendingStripPatch<TDim> >(*it);

                std::vector<std::size_t> patch_indices = pBendPatch->GetIndicesFromParent();
                pBendPatch->pFESpace()->ResetFunctionIndices(patch_indices);
            }
        }

        // collect all the enumerated numbers and reassign with new to make it consecutive
        std::set<std::size_t> all_indices;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            all_indices.insert((*it)->pFESpace()->FunctionIndices().begin(), (*it)->pFESpace()->FunctionIndices().end());
        }
        mEquationSystemSize = all_indices.size();

        std::map<std::size_t, std::size_t> new_indices;
        std::size_t cnt = start;
        for (std::set<std::size_t>::iterator it = all_indices.begin(); it != all_indices.end(); ++it)
        {
            new_indices[*it] = cnt++;
        }

        // reassign the new indices to each patch
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->UpdateFunctionIndices(new_indices);
        }

        // rebuild the global to patch map
        mGlobalToPatch.clear();
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            const std::vector<std::size_t>& global_indices = (*it)->pFESpace()->FunctionIndices();
            for (std::size_t i = 0; i < global_indices.size(); ++i)
                mGlobalToPatch[global_indices[i]] = (*it)->Id();
        }

        // turn on the enumerated flag
        mIsEnumerated = true;

        return start + mEquationSystemSize;
    }

    /// Make the two patches neighbor. This requires that two patches are conformed at the interface.
    static void MakeNeighbor(typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2)
    {
        typename FESpace<TDim-1>::Pointer pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);
        typename FESpace<TDim-1>::Pointer pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2);

        if( (*pBFESpace1) == (*pBFESpace2) )
        {
            // set the neighbor information
            pPatch1->pSetNeighbor(side1, pPatch2);
            pPatch2->pSetNeighbor(side2, pPatch1);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The two patch's boundaries are not conformed", "")
    }

    /// Information
    void PrintAddress(typename PatchType::Pointer pPatch)
    {
        std::cout << pPatch << std::endl;
    }

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch overview: Number of patches = " << mpPatches.size();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch details:" << std::endl;
        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
            rOStream << (*it) << std::endl;
    }

private:

    PatchContainerType mpPatches; // container for all the patches
    bool mIsEnumerated;
    std::size_t mEquationSystemSize; // this is the number of equation id in this multipatch
    std::map<std::size_t, std::size_t> mGlobalToPatch; // this is to map each global id to a patch id

};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatch<TDim>& rThis)
{
    rOStream << ">>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    rOStream << "-------------Begin MultiPatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << "-------------End MultiPatchInfo-------------" << std::endl;
    rOStream << ">>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    return rOStream;
}

} // end namespace Kratos

#endif
