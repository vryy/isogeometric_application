//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Apr 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/patch.h"

#define DEBUG_DESTROY

namespace Kratos
{

/**
 * This class represents an interface connecting two patches.
 */
template<int TDim>
class PatchInterface
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PatchInterface);

    typedef Patch<TDim> PatchType;

    /// Empty Constructor, be careful when using
    PatchInterface()
    : mSide1(_NUMBER_OF_BOUNDARY_SIDE), mSide2(_NUMBER_OF_BOUNDARY_SIDE)
    {}

    /// Full Constructor
    PatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide& side1,
        typename PatchType::Pointer pPatch2, const BoundarySide& side2)
    : mpPatch1(pPatch1->shared_from_this()), mpPatch2(pPatch2->shared_from_this()), mSide1(side1), mSide2(side2)
    {}

    /// Destructor
    virtual ~PatchInterface()
    {
        #ifdef DEBUG_DESTROY
        std::cout << "PatchInterface" << TDim << "D, Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Get the first patch
    typename PatchType::Pointer pPatch1() const {return mpPatch1.lock();}

    /// Get the second patch
    typename PatchType::Pointer pPatch2() const {return mpPatch2.lock();}

    /// Get the side of the first patch where the strip locates
    const BoundarySide& Side1() const {return mSide1;}

    /// Get the side of the second patch where the strip locates
    const BoundarySide& Side2() const {return mSide2;}

    /// Get the indices of the control points on this strip from the parent patch
    virtual std::vector<std::size_t> GetIndicesFromParent() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PatchInterface" << TDim << "D, Addr = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        if (this->pPatch1() == NULL)
            rOStream << "patch 1 is not specified" << std::endl;
        else
            rOStream << ">> patch 1: " << *(this->pPatch1()) << std::endl;
        if (this->pPatch2() == NULL)
            rOStream << "patch 2 is not specified" << std::endl;
        else
            rOStream << ">> patch 2: " << *(this->pPatch2()) << std::endl;
        rOStream << ">> side 1: " << mSide1 << std::endl;
        rOStream << ">> side 2: " << mSide2 << std::endl;
    }

private:

    BoundarySide mSide1;
    BoundarySide mSide2;

    typename PatchType::WeakPointer mpPatch1;
    typename PatchType::WeakPointer mpPatch2;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const PatchInterface<TDim>& rThis)
{
    rOStream << "-------------Begin PatchInterfaceInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End PatchInterfaceInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_H_INCLUDED defined

