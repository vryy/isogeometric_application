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

namespace Kratos
{

/**
 * This class represents an interface connecting two patches.
 * The interface is designed using half-edge philosophy, in which the interface keeps a pointer to the interface of the other side.
 */
template<int TDim, typename TLocalCoordinateType = double, typename TCoordinateType = double, typename TDataType = double>
class PatchInterface
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : public std::enable_shared_from_this<PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#else
    : public boost::enable_shared_from_this<PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#endif
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PatchInterface);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const PatchInterface> ConstPointer;
#endif

    typedef Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;
    typedef PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchInterfaceType;

    /// Empty Constructor, be careful when using
    PatchInterface()
        : mSide1(_NUMBER_OF_BOUNDARY_SIDE), mSide2(_NUMBER_OF_BOUNDARY_SIDE)
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    PatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                   typename PatchType::Pointer pPatch2, const BoundarySide side2)
        : mpPatch1(pPatch1->shared_from_this()), mpPatch2(pPatch2->shared_from_this())
        , mSide1(side1), mSide2(side2)
    {}

    /// Destructor
    virtual ~PatchInterface()
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "PatchInterface" << TDim << "D, Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Create a clone of this interface
    virtual typename PatchInterfaceType::Pointer Clone() const
    {
        return typename PatchInterfaceType::Pointer(new PatchInterfaceType(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2()));
    }

    /// Get/Set the other half interface
    void SetOtherInterface(typename PatchInterfaceType::Pointer pOther) {mpOtherInterface = pOther->shared_from_this();}
    typename PatchInterfaceType::Pointer pOtherInterface() const {return mpOtherInterface.lock();}

    /// Get/Set the first patch
    void SetPatch1(typename PatchType::Pointer pPatch) {mpPatch1 = pPatch->shared_from_this();}
    typename PatchType::Pointer pPatch1() {return mpPatch1.lock();}
    typename PatchType::Pointer pPatch1() const {return mpPatch1.lock();}
    // typename PatchType::ConstPointer pPatch1() const {return mpPatch1.lock();}

    /// Get/Set the second patch
    void SetPatch2(typename PatchType::Pointer pPatch) {mpPatch2 = pPatch->shared_from_this();}
    typename PatchType::Pointer pPatch2() {return mpPatch2.lock();}
    typename PatchType::Pointer pPatch2() const {return mpPatch2.lock();}
    // typename PatchType::ConstPointer pPatch2() const {return mpPatch2.lock();}

    /// Get the side of the first patch where the strip locates
    BoundarySide Side1() const {return mSide1;}

    /// Set the side 1
    void SetSide1(BoundarySide side) {mSide1 = side;}

    /// Flip the boundary side 1
    void FlipSide1() {mSide1 = ReversedBoundarySide::Get(mSide1);}

    /// Get the side of the second patch where the strip locates
    BoundarySide Side2() const {return mSide2;}

    /// Set the side 2
    void SetSide2(BoundarySide side) {mSide2 = side;}

    /// Flip the boundary side 2
    void FlipSide2() {mSide2 = ReversedBoundarySide::Get(mSide2);}

    /// Get the indices of the control points on this strip from the parent patch
    virtual std::vector<std::size_t> GetIndicesFromParent() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Validate the compatibility of two patches on the interface
    virtual bool Validate(const bool debug, const double dist_tol) const
    {
        return true;
    }

    /// Enumerate on the interface, i.e. to make sure that the enumeration on the two patch interfaces are compatible
    virtual void Enumerate()
    {
        std::vector<std::size_t> func_indices = this->pPatch1()->pFESpace()->ExtractBoundaryFunctionIndices(this->Side1());
        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices, false);
    }

    /// Overload operator ==
    virtual bool operator==(const PatchInterfaceType& rOther) const
    {
        if (this->pPatch1() != rOther.pPatch1()) { return false; }
        if (this->pPatch2() != rOther.pPatch2()) { return false; }
        if (this->Side1() != rOther.Side1()) { return false; }
        if (this->Side2() != rOther.Side2()) { return false; }
        return true;
    }

    /// Return information as string
    virtual std::string Info() const
    {
        if constexpr (TDim == 1)
            return "PatchInterface1D";
        else if constexpr (TDim == 2)
            return "PatchInterface2D";
        else if constexpr (TDim == 3)
            return "PatchInterface3D";
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << ", Addr = " << this;
        rOStream << ", Patch ";
        if (this->pPatch1() == NULL)
        {
            rOStream << "null";
        }
        else
        {
            rOStream << this->pPatch1()->Id() << "(" << this->pPatch1() << ")"
                     << ": " << this->Side1();
        }
        rOStream << " - Patch ";
        if (this->pPatch2() == NULL)
        {
            rOStream << "null";
        }
        else
        {
            rOStream << this->pPatch2()->Id() << "(" << this->pPatch2() << ")"
                     << ": " << this->Side2();
        }
        rOStream << ", Opposite Interface: ";
        if (this->pOtherInterface() == NULL) { rOStream << "null"; }
        else { rOStream << this->pOtherInterface(); }
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        if (this->pPatch1() == NULL)
        {
            rOStream << "patch 1 is not specified" << std::endl;
        }
        else
        {
            rOStream << ">> patch 1: " << *(this->pPatch1()) << std::endl;
        }
        if (this->pPatch2() == NULL)
        {
            rOStream << "patch 2 is not specified" << std::endl;
        }
        else
        {
            rOStream << ">> patch 2: " << *(this->pPatch2()) << std::endl;
        }
        rOStream << ">> side 2: " << mSide2 << std::endl;
    }

private:

    BoundarySide mSide1;
    BoundarySide mSide2;

    typename PatchType::WeakPointer mpPatch1;
    typename PatchType::WeakPointer mpPatch2;

    typename PatchInterfaceType::WeakPointer mpOtherInterface;
}; // class PatchInterface

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType>& rThis)
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_H_INCLUDED defined
