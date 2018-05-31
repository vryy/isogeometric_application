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
#include "custom_utilities/patch_interface.h"

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
    MultiPatch() : mEquationSystemSize(0) {}

    /// Destructor
    virtual ~MultiPatch() {}

    /// Add the patch
    void AddPatch(typename Patch<TDim>::Pointer pPatch)
    {
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
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
    bool IsEnumerated() const
    {
        return (this->EquationSystemSize() != 0) && ((this->GetLastEquationId() - this->GetFirstEquationId() + 1) == this->EquationSystemSize());
    }

    /// Get the equation system size
    std::size_t EquationSystemSize() const {return mEquationSystemSize;}

    /// Locate the patch of the global equation id and the corresponding local id to determine the control value
    /// IMPORTANT: user must make sure that the multipatch is fully enumerated by checking IsEnumerated()
    std::tuple<std::size_t, std::size_t> EquationIdLocation(const std::size_t& global_id) const
    {
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

            for (std::size_t i = 0; i < it->NumberOfInterfaces(); ++i)
            {
                bool check2 = it->pInterface(i)->Validate();
                if (!check2)
                    return false;
            }
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

    /// Get the first equation_id accross all patches
    std::size_t GetFirstEquationId() const
    {
        std::size_t first_id;

        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            std::size_t patch_first_id = it->pFESpace()->GetFirstEquationId();

            if (patch_first_id == -1)
            {
                return -1;
            }
            else
            {
                if (it == this->begin())
                {
                    first_id = patch_first_id;
                }
                else
                {
                    if (patch_first_id < first_id)
                        first_id = patch_first_id;
                }
            }
        }

        return first_id;
    }

    /// Get the last equation_id accross all patches
    std::size_t GetLastEquationId() const
    {
        std::size_t last_id = -1;
        bool hit = false;

        for (typename PatchContainerType::const_iterator it = this->begin(); it != this->end(); ++it)
        {
            std::size_t patch_last_id = it->pFESpace()->GetLastEquationId();

            if (patch_last_id != -1)
            {
                if (!hit)
                {
                    last_id = patch_last_id;
                    hit = true;
                }
                else
                {
                    if (patch_last_id > last_id)
                        last_id = patch_last_id;
                }
            }
        }

        return last_id;
    }

    /// Reset global ids for each patch
    /// In principle, it initializes all the equation_id to -1
    void ResetFunctionIndices()
    {
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->ResetFunctionIndices();
        }
    }

    /// Enumerate all the patches, starting at 0
    std::size_t Enumerate()
    {
        return this->Enumerate(0);
    }

    /// Enumerate all the patches, with the given starting id
    std::size_t Enumerate(const std::size_t& start)
    {
        // enumerate each patch
        std::size_t last = start;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsPrimary() == true)
            {
                last = (*it)->pFESpace()->Enumerate(last);
                std::cout << "At " << __FUNCTION__ << ", last: " << last << std::endl;

                // enumerate the interface
                for (std::size_t i = 0; i < (*it)->NumberOfInterfaces(); ++i)
                    (*it)->pInterface(i)->Enumerate();
            }
        }

        // check if a patch is not a primary patch, then that patch must be enumerated again using the enumeration info from the other patches
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsPrimary() == false)
            {
                (*it)->Enumerate();
            }
        }

        // collect all the enumerated numbers and reassign with new to make it consecutive
        std::set<std::size_t> all_indices;
        for (typename PatchContainerType::ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            std::vector<std::size_t> func_indices = (*it)->pFESpace()->FunctionIndices();
            all_indices.insert(func_indices.begin(), func_indices.end());
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
            std::vector<std::size_t> global_indices = (*it)->pFESpace()->FunctionIndices();
            for (std::size_t i = 0; i < global_indices.size(); ++i)
                mGlobalToPatch[global_indices[i]] = (*it)->Id();
        }

        return start + mEquationSystemSize;
    }

    /// Information
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

