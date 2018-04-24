//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Apr 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_PATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_PATCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/patch.h"


namespace Kratos
{

/**
This class represents an isogeometric bending strip patch connecting two patches.
Idea is taken from Kiendl et al, The bending strip method for isogeometric analysis of Kirchhoff–Love shell structures comprised of multiple patches.
 */
template<int TDim>
class BendingStripPatch : public Patch<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BendingStripPatch);

    typedef Patch<TDim> BaseType;
    typedef typename BaseType::ControlPointType ControlPointType;

    /// Default Constructor
    BendingStripPatch(const std::size_t& Id, const int& Order) : BaseType(Id), mOrder(Order)
    {
    }

    /// Full Constructor
    BendingStripPatch(const std::size_t& Id,
        typename BaseType::Pointer pPatch1, const BoundarySide& side1,
        typename BaseType::Pointer pPatch2, const BoundarySide& side2,
        const int& Order)
    : BaseType(Id), mOrder(Order), mpPatch1(pPatch1), mpPatch2(pPatch2), mSide1(side1), mSide2(side2)
    {
    }

    /// Destructor
    virtual ~BendingStripPatch()
    {
        #ifdef DEBUG_DESTROY
        std::cout << Type() << ", Id = " << Id() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Return true if this patch is a bending strip patch
    virtual bool IsBendingStrip() const
    {
        return true;
    }

    /// Get the order of the strip patch in the orthogonal direction
    const int& Order() const {return mOrder;}

    /// Get the first patch
    typename BaseType::Pointer pPatch1() const {return mpPatch1;}

    /// Get the second patch
    typename BaseType::Pointer pPatch2() const {return mpPatch2;}

    /// Get the side of the first patch where the strip locates
    const BoundarySide& Side1() const {return mSide1;}

    /// Get the side of the second patch where the strip locates
    const BoundarySide& Side2() const {return mSide2;}

    /// Get the string representing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "BendingStripPatch" << TDim << "D";
        return ss.str();
    }

    /// Get the indices of the control points on this strip from the parent patch
    virtual std::vector<std::size_t> GetIndicesFromParent() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

private:
    int mOrder; // this is the bending strip order in the normal direction to the boundary

    BoundarySide mSide1;
    BoundarySide mSide2;

    typename BaseType::Pointer mpPatch1;
    typename BaseType::Pointer mpPatch2;
};

//template<>
//class BendingStripPatch<1> : public Patch<1>
//{
//public:
//    /// Pointer definition
//    KRATOS_CLASS_POINTER_DEFINITION(BendingStripPatch);

//    typedef Patch<1> BaseType;
//    typedef typename BaseType::ControlPointType ControlPointType;

//    /// Default Constructor
//    BendingStripPatch(const std::size_t& Id, const int& Order) : BaseType(Id), mOrder(Order)
//    {
//    }

//    /// Full Constructor
//    BendingStripPatch(const std::size_t& Id,
//        typename BaseType::Pointer pPatch1, const BoundarySide& side1,
//        typename BaseType::Pointer pPatch2, const BoundarySide& side2,
//        const int& Order) : BaseType(Id), mOrder(Order), mpPatch1(pPatch1), mpPatch2(pPatch2)
//    {
//    }

//    /// Destructor
//    virtual ~BendingStripPatch()
//    {
//        #ifdef DEBUG_DESTROY
//        std::cout << Type() << ", Id = " << Id() << ", Addr = " << this << " is destroyed" << std::endl;
//        #endif
//    }
//};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const BendingStripPatch<TDim>& rThis)
{
    rOStream << "-------------Begin BendingStripPatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End BendingStripPatchInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_PATCH_H_INCLUDED defined

