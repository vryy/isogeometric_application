//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Apr 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/patch_interface.h"

#define DEBUG_DESTROY

namespace Kratos
{

/**
 * This class represents an interface connecting two patches.
 * The interface is designed using half-edge philosophy, in which the interface keeps a pointer to the interface of the other side.
 */
template<int TDim>
class BSplinesPatchInterface : public PatchInterface<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<TDim> BaseType;
    typedef typename BaseType::PatchType PatchType;

    /// Empty Constructor, be careful when using
    BSplinesPatchInterface() : BaseType()
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide& side1,
        typename PatchType::Pointer pPatch2, const BoundarySide& side2)
    : BaseType(pPatch1, side1, pPatch2, side2)
    {
        if (TDim == 2)
        {
            mLocalParameterMap[0] = 0;
            mDirections[0] = _FORWARD_;
        }
        else if (TDim == 3)
        {
            mLocalParameterMap[0] = 0;
            mLocalParameterMap[1] = 1;
            mDirections[0] = _FORWARD_;
            mDirections[1] = _FORWARD_;
        }
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 2D.
    /// The direction variable is to indicate if the local parameters shall be reversed.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide& side1,
        typename PatchType::Pointer pPatch2, const BoundarySide& side2,
        const BoundaryDirection& direction)
    : BaseType(pPatch1, side1, pPatch2, side2)
    {
        if (TDim == 3)
            KRATOS_THROW_ERROR(std::logic_error, "This constructor is not designed for 3D", "")

        mLocalParameterMap[0] = 0;
        mDirections[0] = direction;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 3D.
    /// if uv_or_vu is true, the local configuration {(u:ub) - (v:vb)} is assumed; otherwise, it is {(u:vb) - (v:ub)}
    /// The direction indicates if the respective local parameters shall be reversed or not
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide& side1,
        typename PatchType::Pointer pPatch2, const BoundarySide& side2,
        const bool& uv_or_vu,
        const BoundaryDirection& direction1,
        const BoundaryDirection& direction2)
    : BaseType(pPatch1, side1, pPatch2, side2)
    {
        if (TDim == 2)
            KRATOS_THROW_ERROR(std::logic_error, "This constructor is not designed for 2D", "")

        if (uv_or_vu == true)
        {
            mLocalParameterMap[0] = 0;
            mLocalParameterMap[1] = 1;
        }
        else
        {
            mLocalParameterMap[0] = 1;
            mLocalParameterMap[1] = 0;
        }

        mDirections[0] = direction1;
        mDirections[1] = direction2;
    }

    /// Destructor
    virtual ~BSplinesPatchInterface()
    {
        #ifdef DEBUG_DESTROY
        std::cout << "BSplinesPatchInterface" << TDim << "D, Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Create a clone of this interface
    virtual typename BaseType::Pointer Clone() const
    {
        if (TDim == 2)
        {
            return typename BaseType::Pointer(new BSplinesPatchInterface<TDim>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), this->Direction(0)));
        }
        else if (TDim == 3)
        {
            if (this->LocalParameterMapping(0) == 0)
                return typename BaseType::Pointer(new BSplinesPatchInterface<TDim>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), true, this->Direction(0), this->Direction(1)));
            else if (this->LocalParameterMapping(0) == 1)
                return typename BaseType::Pointer(new BSplinesPatchInterface<TDim>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), false, this->Direction(0), this->Direction(1)));
        }
    }

    /// Get the local parameter space mapping
    std::size_t LocalParameterMapping(const std::size_t& dim) const
    {
        std::map<std::size_t, std::size_t>::const_iterator it = mLocalParameterMap.find(dim);
        if (it != mLocalParameterMap.end())
            return it->second;
        else
            KRATOS_THROW_ERROR(std::logic_error, "The dimension is invalid", "")
    }

    /// Get the relative direction between the local parameter space
    const BoundaryDirection& Direction(const std::size_t& dim) const {return mDirections[dim];}

    /// Validate the compatibility of two patches on the interface
    virtual bool Validate() const
    {
        // TODO validate the parameter direction

        typename Patch<TDim-1>::Pointer BPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename Patch<TDim-1>::Pointer BPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

        return (*BPatch1) == (*BPatch2);
    }

    /// Enumerate on the interface, i.e. to make sure that the enumeration on the two patch interfaces are compatible
    virtual void Enumerate()
    {
        std::vector<std::size_t> func_indices = this->pPatch1()->pFESpace()->ExtractBoundaryFunctionIndices(this->Side1());
        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }

private:

    std::map<std::size_t, std::size_t> mLocalParameterMap;
    // this variable stores the information about the local mapping at the interface. In 2D, it is simply 0:0.
    // In 3D, it shall have the values {0:0, 1:1} or {0:1, 1:0}. The latter values corresponds to the {(u:vb) - (v:ub)} configuration.

    boost::array<BoundaryDirection, TDim> mDirections;
    // this variable indicates if the relative local direction shall be reversed of each other.

};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesPatchInterface<TDim>& rThis)
{
    rOStream << "-------------Begin BSplinesPatchInterfaceInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End BSplinesPatchInterfaceInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED defined

