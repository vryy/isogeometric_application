//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED

#include "custom_utilities/patch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/patch_interface_utility.h"

namespace Kratos
{

/**
 * Base class for all multipatch
 */
class BaseMultiPatch
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BaseMultiPatch);

    /// Default constructor
    BaseMultiPatch()
    {}

    virtual ~BaseMultiPatch()
    {}

    /// Provide the information on dimension
    virtual int Dimension() const {return 0;}

    /// Provide the identifier for local coordinate type
    virtual const char* LocalCoordinateTypeStr() const {return "";}

    /// Provide the identifier for physical coordinate type
    virtual const char* CoordinateTypeStr() const {return "";}

    /// Provide the identifier for underlying data type
    virtual const char* DataTypeStr() const {return "";}

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {}

    virtual void load(Serializer& rSerializer)
    {}
    ///@}
};

/**
 * This class represents an isogeometric multipatch in parametric coordinates. An isogeometric multipatch comprises a list of similar type patches, i.e NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
 */
template<int TDim, typename TLocalCoordinateType = double, typename TCoordinateType = double, typename TDataType = double>
class MultiPatch
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : public std::enable_shared_from_this<MultiPatch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#else
    : public boost::enable_shared_from_this<MultiPatch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#endif
    , public BaseMultiPatch
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatch);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const MultiPatch> ConstPointer;
#endif
    /// Type definition
    typedef Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;
    typedef PointerVectorSet<PatchType, IndexedObject> PatchContainerType;

    typedef typename PatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename PatchType::CoordinateType CoordinateType;
    typedef typename PatchType::DataType DataType;
    typedef typename PatchType::CoordinateValueType CoordinateValueType;

    typedef typename PatchType::vertex_t vertex_t;
    typedef typename PatchType::edge_t edge_t;
    typedef typename PatchType::face_t face_t;
    typedef typename PatchType::volume_t volume_t;

    typedef typename PatchContainerType::iterator patch_iterator;
    typedef typename PatchContainerType::const_iterator patch_const_iterator;
    typedef typename PatchContainerType::ptr_iterator patch_ptr_iterator;
    typedef typename PatchContainerType::ptr_const_iterator patch_ptr_const_iterator;

    typedef typename PatchType::PatchInterfaceType PatchInterfaceType;
    typedef typename PatchType::interface_iterator interface_iterator;
    typedef typename PatchType::interface_const_iterator interface_const_iterator;

    /// Constants
    static constexpr int Dim = TDim;

    /// Default constructor
    MultiPatch() : mEquationSystemSize(0), mpPatches() {}

    /// The copy constructor does not make sense for MultiPatch because during adding patches,
    /// a shared_ptr to this patch is needed, which is not yet available. This will lead to
    /// the error RuntimeError: tr1::bad_weak_ptr.
    /// Therefore, the copy constructor is disabled completely.
    MultiPatch(const MultiPatch& rOther) = delete;

    /// Destructor
    ~MultiPatch() override {}

    /// Clone the multipatch
    MultiPatch::Pointer Clone() const
    {
        MultiPatch::Pointer pNewMultiPatch = MultiPatch::Pointer(new MultiPatch());
        *pNewMultiPatch = *this;
        return pNewMultiPatch;
    }

    /// (derived!)
    int Dimension() const override {return TDim;}

    /// (derived!)
    const char* LocalCoordinateTypeStr() const override {return DataTypeToString<TLocalCoordinateType>::Get();}

    /// (derived!)
    const char* CoordinateTypeStr() const override {return DataTypeToString<TCoordinateType>::Get();}

    /// (derived!)
    const char* DataTypeStr() const override {return DataTypeToString<TDataType>::Get();}

    /// Create an empty multipatch instance
    static MultiPatch::Pointer Create()
    {
        return MultiPatch::Pointer(new MultiPatch());
    }

    /// Cast a BaseMultiPatch to this multipatch type
    static MultiPatch::Pointer Cast(BaseMultiPatch::Pointer pMultiPatch)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        return std::dynamic_pointer_cast<MultiPatch>(pMultiPatch);
#else
        return boost::dynamic_pointer_cast<MultiPatch>(pMultiPatch);
#endif
    }

    /// Add the patch
    void AddPatch(typename PatchType::Pointer pPatch)
    {
        const auto it = mpPatches.find(pPatch->Id());
        if (it != mpPatches.end())
            KRATOS_ERROR << "WARNING!!! Patch " << pPatch->Id() << " exists";
        mpPatches.push_back(pPatch);
        pPatch->pSetParentMultiPatch(this->shared_from_this());
    }

    /// Add the patch
    void RemovePatch(typename PatchType::Pointer pPatch)
    {
        // look for all interfaces and remove the corresponding interfaces in other patches
        for (interface_iterator it = pPatch->InterfaceBegin(); it != pPatch->InterfaceEnd(); ++it)
        {
            typename PatchType::Pointer pNeighborPatch = (*it)->pPatch2();

            for (interface_iterator it2 = pNeighborPatch->InterfaceBegin(); it2 != pNeighborPatch->InterfaceEnd(); ++it2)
            {
                if ((*it2)->pPatch2() == pPatch)
                {
                    pNeighborPatch->RemoveInterface(*it2);
                    break;
                }
            }
        }

        // remove the patch
        mpPatches.erase(pPatch->Id());

        // modify the global to patch map data
        for (std::map<std::size_t, std::size_t>::iterator it = mGlobalIdToPatchId.begin(); it != mGlobalIdToPatchId.end();)
        {
            if (it->second == pPatch->Id())
            {
                mGlobalIdToPatchId.erase(it++);
            }
            else
            {
                ++it;
            }
        }
    }

    /// Reset Id for all the patches
    void ResetId()
    {
        std::size_t Id = 0;
        for (patch_iterator it = this->begin(); it != this->end(); ++it)
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
    std::tuple<std::size_t, std::size_t> EquationIdLocation(std::size_t global_id) const
    {
        std::map<std::size_t, std::size_t>::const_iterator it = mGlobalIdToPatchId.find(global_id);
        if (it == mGlobalIdToPatchId.end())
        {
            KRATOS_WATCH(global_id)
            KRATOS_WATCH(mEquationSystemSize)
            std::cout << "global_to_patch map:" << std::endl;
            for (std::map<std::size_t, std::size_t>::const_iterator it2 = mGlobalIdToPatchId.begin(); it2 != mGlobalIdToPatchId.end(); ++it2)
            {
                std::cout << " " << it2->first << ": " << it2->second << std::endl;
            }
            KRATOS_ERROR << "The global id " << global_id << " does not exist in the global_to_patch map.";
        }

        std::size_t patch_id = it->second;
        std::size_t local_id = pGetPatch(it->second)->pFESpace()->LocalId(global_id);

        return std::make_tuple(patch_id, local_id);
    }

    /// Validate the MultiPatch
    virtual bool Validate() const
    {
        for (patch_const_iterator it = this->begin(); it != this->end(); ++it)
        {
            bool check = it->Validate();
            if (!check)
            {
                return false;
            }

            for (interface_const_iterator it2 = it->InterfaceBegin(); it2 != it->InterfaceEnd(); ++it2)
            {
                bool check2 = (*it2)->Validate(false, PatchType::DISTANCE_TOLERANCE);
                if (!check2)
                {
                    return false;
                }
            }
        }

        return true;
    }

    /// Get the patch with specific Id
    typename PatchType::Pointer pGetPatch(std::size_t Id)
    {
        patch_ptr_iterator it_patch = mpPatches.find(Id).base();
        if (it_patch == mpPatches.ptr_end())
        {
            KRATOS_ERROR << "The patch with Id = " << Id << " does not exist in the multipatch";
        }
        return *it_patch;
    }

    /// Get the patch with specific Id
    typename PatchType::ConstPointer pGetPatch(std::size_t Id) const
    {
        typename PatchContainerType::ptr_const_iterator it_patch = mpPatches.find(Id).base();
        if (it_patch == mpPatches.ptr_end())
        {
            KRATOS_ERROR << "The patch with Id = " << Id << " does not exist in the multipatch";
        }
        return typename PatchType::ConstPointer(*it_patch);
    }

    /// Access the underlying list of patches
    /// WARNING!!! be careful with this routine
    PatchContainerType& Patches() {return mpPatches;}

    /// Access the underlying list of patches
    const PatchContainerType& Patches() const {return mpPatches;}

    /// iterators
    patch_iterator begin() {return mpPatches.begin();}
    patch_ptr_iterator ptr_begin() {return mpPatches.ptr_begin();}
    patch_const_iterator begin() const {return mpPatches.begin();}
    patch_ptr_const_iterator ptr_begin() const {return mpPatches.ptr_begin();}
    patch_iterator end() {return mpPatches.end();}
    patch_ptr_iterator ptr_end() {return mpPatches.ptr_end();}
    patch_const_iterator end() const {return mpPatches.end();}
    patch_ptr_const_iterator ptr_end() const {return mpPatches.ptr_end();}

    /// Get the number of patches
    std::size_t size() const {return mpPatches.size();}

    /// Get the first equation_id accross all patches
    std::size_t GetFirstEquationId() const
    {
        std::size_t first_id;

        for (patch_const_iterator it = this->begin(); it != this->end(); ++it)
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
                    {
                        first_id = patch_first_id;
                    }
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

        for (patch_const_iterator it = this->begin(); it != this->end(); ++it)
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
                    {
                        last_id = patch_last_id;
                    }
                }
            }
        }

        return last_id;
    }

    /// Reset global ids for each patch
    /// In principle, it initializes all the equation_id to -1
    void ResetFunctionIndices()
    {
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->ResetFunctionIndices();
        }
    }

    /// Enumerate all the patches, starting at 0
    std::size_t Enumerate()
    {
        this->ResetFunctionIndices();
        return this->Enumerate(0);
    }

    /// Enumerate all the patches, with the given starting id
    std::size_t Enumerate(std::size_t start)
    {
        // enumerate each patch
        std::size_t last = start;
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsPrimary() == true)
            {
                last = (*it)->pFESpace()->Enumerate(last);

                // enumerate the interface
                for (interface_iterator it2 = (*it)->InterfaceBegin(); it2 != (*it)->InterfaceEnd(); ++it2)
                {
                    (*it2)->Enumerate();
                }
            }
        }

        // check if a patch is not a primary patch, then that patch must be enumerated again using the enumeration info from the other patches
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            if ((*it)->IsPrimary() == false)
            {
                (*it)->Enumerate();
            }
        }

        // collect all the enumerated numbers and reassign with new to make it consecutive
        std::set<std::size_t> all_indices;
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
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
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            (*it)->pFESpace()->UpdateFunctionIndices(new_indices);
        }

        // rebuild the global to patch map
        mGlobalIdToPatchId.clear();
        for (patch_ptr_iterator it = Patches().ptr_begin(); it != Patches().ptr_end(); ++it)
        {
            std::vector<std::size_t> global_indices = (*it)->pFESpace()->FunctionIndices();
            for (std::size_t i = 0; i < global_indices.size(); ++i)
            {
                mGlobalIdToPatchId[global_indices[i]] = (*it)->Id();
            }
        }

        return start + mEquationSystemSize;
    }

    /// Compute the local coordinates of a point. On output returns the index of the patch containing the point if found, -1 otherwise
    int LocalCoordinates(const array_1d<CoordinateType, 3>& point, array_1d<LocalCoordinateType, 3>& xi) const
    {
        int error_code;

        const array_1d<LocalCoordinateType, 3> xi0 = xi;
        for (patch_const_iterator it = this->begin(); it != this->end(); ++it)
        {
            noalias(xi) = xi0;
            error_code = it->LocalCoordinates(point, xi);
            if (error_code == 0)
            {
                return it->Id();
            }
        }

        return -1;
    }

    /// Overload assignment operator
    MultiPatch& operator=(const MultiPatch& rOther)
    {
        this->mEquationSystemSize = rOther.mEquationSystemSize;
        this->mGlobalIdToPatchId = rOther.mGlobalIdToPatchId;

        /* clone the patches */
        for (patch_const_iterator itp = rOther.begin(); itp != rOther.end(); ++itp)
        {
            this->AddPatch(itp->Clone());
        }

        /* clone the interfaces */
        for (patch_const_iterator itp = rOther.begin(); itp != rOther.end(); ++itp)
        {
            typename PatchType::Pointer pNewPatch = this->pGetPatch(itp->Id());

            for (interface_const_iterator iti = itp->InterfaceBegin(); iti != itp->InterfaceEnd(); ++iti)
            {
                typename PatchType::Pointer pOldNeighborPatch = (*iti)->pPatch2();
                typename PatchType::Pointer pNewNeighborPatch = this->pGetPatch(pOldNeighborPatch->Id());

                auto pNewInterface = (*iti)->Create(pNewPatch, pNewNeighborPatch);
                pNewPatch->AddInterface(pNewInterface);

                // search and assign the other interface
                for (interface_const_iterator iti2 = pNewNeighborPatch->InterfaceBegin(); iti2 != pNewNeighborPatch->InterfaceEnd(); ++iti2)
                {
                    if ((*iti2)->Side1() == pNewInterface->Side2()
                     && (*iti2)->pPatch1()->Id() == pNewInterface->pPatch2()->Id())
                    {
                        pNewInterface->SetOtherInterface(*iti2);
                        (*iti2)->SetOtherInterface(pNewInterface);
                    }
                }
            }
        }

        return *this;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch overview: Number of patches = " << mpPatches.size();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "MultiPatch details:" << std::endl;
        for (patch_const_iterator it = this->begin(); it != this->end(); ++it)
        {
            rOStream << (*it) << std::endl;
        }
    }

private:

    PatchContainerType mpPatches; // container for all the patches
    std::size_t mEquationSystemSize; // this is the number of equation id in this multipatch
    std::map<std::size_t, std::size_t> mGlobalIdToPatchId; // this is to map each global id to a patch id

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        rSerializer.save("size", this->size());

        std::size_t ip = 0;
        for (patch_const_iterator itp = this->begin(); itp != this->end(); ++itp, ++ip)
        {
            std::stringstream ss_name;
            ss_name << "Patch_" << ip;
            rSerializer.save(ss_name.str(), *itp);
        }

        ip = 0;
        for (patch_const_iterator itp = this->begin(); itp != this->end(); ++itp, ++ip)
        {
            std::stringstream ss_isize;
            ss_isize << "Patch_" << ip << "_Interface_size";
            rSerializer.save(ss_isize.str(), itp->NumberOfInterfaces());

            std::size_t ii = 0;
            for (interface_const_iterator iti = itp->InterfaceBegin(); iti != itp->InterfaceEnd(); ++iti, ++ii)
            {
                std::stringstream ss_type_name;
                ss_type_name << "Patch_" << ip << "_Interface_" << ii << "_type";
                rSerializer.save(ss_type_name.str(), (*iti)->Type());

                std::stringstream ss_patch1_Id;
                ss_patch1_Id << "Patch_" << ip << "_Interface_" << ii << "_patch1_id";
                if ((*iti)->pPatch1() != nullptr)
                    rSerializer.save(ss_patch1_Id.str(), (*iti)->pPatch1()->Id());
                else
                    rSerializer.save(ss_patch1_Id.str(), 0);

                std::stringstream ss_patch2_Id;
                ss_patch2_Id << "Patch_" << ip << "_Interface_" << ii << "_patch2_id";
                if ((*iti)->pPatch2() != nullptr)
                    rSerializer.save(ss_patch2_Id.str(), (*iti)->pPatch2()->Id());
                else
                    rSerializer.save(ss_patch2_Id.str(), 0);

                std::stringstream ss_iname;
                ss_iname << "Patch_" << ip << "_Interface_" << ii;
                rSerializer.save(ss_iname.str(), **iti);
            }
        }

        rSerializer.save("EquationSystemSize", mEquationSystemSize);
        rSerializer.save("GlobalIdToPatchId", mGlobalIdToPatchId);
    }

    void load(Serializer& rSerializer) override
    {
        std::size_t size;
        rSerializer.load("size", size);

        for (std::size_t ip = 0; ip < size; ++ip)
        {
            typename PatchType::Pointer pPatch = typename PatchType::Pointer(new PatchType(0));

            std::stringstream ss_name;
            ss_name << "Patch_" << ip;
            rSerializer.load(ss_name.str(), *pPatch);

            this->AddPatch(pPatch);
            pPatch->pSetParentMultiPatch(this->shared_from_this());
        }

        for (std::size_t ip = 0; ip < size; ++ip)
        {
            std::stringstream ss_isize;
            ss_isize << "Patch_" << ip << "_Interface_size";
            std::size_t isize;
            rSerializer.load(ss_isize.str(), isize);

            for (std::size_t ii = 0; ii < isize; ++ii)
            {
                std::stringstream ss_type_name;
                ss_type_name << "Patch_" << ip << "_Interface_" << ii << "_type";
                std::string itype;
                rSerializer.load(ss_type_name.str(), itype);

                std::stringstream ss_patch1_Id;
                ss_patch1_Id << "Patch_" << ip << "_Interface_" << ii << "_patch1_id";
                std::size_t patch1_Id;
                rSerializer.load(ss_patch1_Id.str(), patch1_Id);

                std::stringstream ss_patch2_Id;
                ss_patch2_Id << "Patch_" << ip << "_Interface_" << ii << "_patch2_id";
                std::size_t patch2_Id;
                rSerializer.load(ss_patch2_Id.str(), patch2_Id);

                auto pNewInterface = PatchInterfaceUtility::template CreateEmptyPatchInterface<PatchInterfaceType>(itype);

                std::stringstream ss_iname;
                ss_iname << "Patch_" << ip << "_Interface_" << ii;
                rSerializer.load(ss_iname.str(), *pNewInterface);

                if (patch1_Id > 0 && patch2_Id > 0)
                {
                    typename PatchType::Pointer pPatch = this->pGetPatch(patch1_Id);
                    typename PatchType::Pointer pOtherPatch = this->pGetPatch(patch2_Id);

                    pNewInterface->SetPatch1(pPatch);
                    pNewInterface->SetPatch2(pOtherPatch);

                    pPatch->AddInterface(pNewInterface);

                    // search and assign the other interface
                    for (interface_const_iterator iti2 = pOtherPatch->InterfaceBegin(); iti2 != pOtherPatch->InterfaceEnd(); ++iti2)
                    {
                        if ((*iti2)->Side1() == pNewInterface->Side2()
                         && (*iti2)->pPatch1()->Id() == pNewInterface->pPatch2()->Id())
                        {
                            pNewInterface->SetOtherInterface(*iti2);
                            (*iti2)->SetOtherInterface(pNewInterface);
                        }
                    }
                }
            }
        }

        rSerializer.load("EquationSystemSize", mEquationSystemSize);
        rSerializer.load("GlobalIdToPatchId", mGlobalIdToPatchId);
    }
    ///@}

}; // class MultiPatch

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const MultiPatch<TDim, TLocalCoordinateType, TCoordinateType, TDataType>& rThis)
{
    rOStream << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    rOStream << "-------------Begin MultiPatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << "-------------End MultiPatchInfo-------------" << std::endl;
    rOStream << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    return rOStream;
}

/// Selector for patch based on dimension
template<int TDim>
class PatchSelector
{
public:
    typedef Patch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE> RealPatch;
    typedef Patch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE> ComplexPatch;
    typedef Patch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE> GComplexPatch;

    typedef PatchInterface<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE> RealPatchInterface;
    typedef PatchInterface<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE> ComplexPatchInterface;
    typedef PatchInterface<TDim, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE> GComplexPatchInterface;

    typedef MultiPatch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE> RealMultiPatch;
    typedef MultiPatch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE> ComplexMultiPatch;
    typedef MultiPatch<TDim, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE> GComplexMultiPatch;

    const RealPatch& GetRealPatch() const {return msRealPatch;}
    const ComplexPatch& GetComplexPatch() const {return msComplexPatch;}
    const GComplexPatch& GetGComplexPatch() const {return msGComplexPatch;}

private:
    static RealPatch msRealPatch;
    static ComplexPatch msComplexPatch;
    static GComplexPatch msGComplexPatch;
};

extern PatchSelector<1> PatchSelector1DInstance;
extern PatchSelector<2> PatchSelector2DInstance;
extern PatchSelector<3> PatchSelector3DInstance;

} // end namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_H_INCLUDED
