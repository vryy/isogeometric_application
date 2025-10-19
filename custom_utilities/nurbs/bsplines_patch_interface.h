//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Apr 2018 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"

namespace Kratos
{

/**
 * This class represents an interface connecting two patches.
 * The interface is designed using half-edge philosophy, in which the interface keeps a pointer to the interface of the other side.
 */
template<int TDim, typename TLocalCoordinateType = double, typename TCoordinateType = double, typename TDataType = double>
class BSplinesPatchInterface : public PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType>
{
};

/**
 * Template specialization for 1D
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class BSplinesPatchInterface<1, TLocalCoordinateType, TCoordinateType, TDataType> : public PatchInterface<1, TLocalCoordinateType, TCoordinateType, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<1, TLocalCoordinateType, TCoordinateType, TDataType> BaseType;
    typedef BSplinesPatchInterface<1, TLocalCoordinateType, TCoordinateType, TDataType> ThisType;
    typedef typename BaseType::PatchType PatchType;

    /// Empty Constructor
    BSplinesPatchInterface() : BaseType()
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
    }

    /// Create a clone of this interface
    typename BaseType::Pointer Clone() const override
    {
        return typename BaseType::Pointer(new ThisType(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2()));
    }

    /// Get the local parameter space mapping
    std::size_t LocalParameterMapping(std::size_t dim) const
    {
        return 0;
    }

    /// Get the relative direction between the local parameter space
    BoundaryDirection Direction() const {return BoundaryDirection::_UNDEFINED_DIR_;}
    BoundaryDirection Direction(std::size_t dim) const {return BoundaryDirection::_UNDEFINED_DIR_;}

    /// Validate the compatibility of two patches on the interface
    bool Validate(const bool debug, const double dist_tol) const override
    {
        typename PatchType::BoundaryPatchType::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename PatchType::BoundaryPatchType::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

        if (debug)
        {
            KRATOS_WATCH(*pBPatch1)
            KRATOS_WATCH(*pBPatch2)
        }

        bool is_valid = true;
        if (dist_tol == 0.0)
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0));
        else
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0, dist_tol));

        return is_valid;
    }

    /// Return information as string
    std::string Info() const override
    {
        return "BSplinesPatchInterface1D";
    }
};

/**
 * Template specialization for 2D
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class BSplinesPatchInterface<2, TLocalCoordinateType, TCoordinateType, TDataType> : public PatchInterface<2, TLocalCoordinateType, TCoordinateType, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<2, TLocalCoordinateType, TCoordinateType, TDataType> BaseType;
    typedef BSplinesPatchInterface<2, TLocalCoordinateType, TCoordinateType, TDataType> ThisType;
    typedef typename BaseType::PatchType PatchType;

    /// Empty Constructor
    BSplinesPatchInterface() : BaseType()
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mDirection = BoundaryDirection::_FORWARD_;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 2D.
    /// The direction variable is to indicate if the local parameters shall be reversed.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const BoundaryDirection direction)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mDirection = direction;
    }

    /// Destructor
    ~BSplinesPatchInterface() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "BSplinesPatchInterface2D, Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Create a clone of this interface
    typename BaseType::Pointer Clone() const override
    {
        return typename BaseType::Pointer(new ThisType(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), this->Direction()));
    }

    /// Get the local parameter space mapping
    std::size_t LocalParameterMapping(std::size_t dim) const
    {
        return 0;
    }

    /// Get the relative direction between the local parameter space
    BoundaryDirection Direction() const {return mDirection;}
    BoundaryDirection Direction(std::size_t dim) const
    {
        if (dim == 0)
        {
            return mDirection;
        }
        else
        {
            return BoundaryDirection::_UNDEFINED_DIR_;
        }
    }

    /// Validate the compatibility of two patches on the interface
    bool Validate(const bool debug, const double dist_tol) const override
    {
        typename PatchType::BoundaryPatchType::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename PatchType::BoundaryPatchType::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

        if (mDirection == BoundaryDirection::_REVERSED_)
        {
            BSplinesPatchUtility::Reverse(pBPatch2, 0);
        }

        if (debug)
        {
            KRATOS_WATCH(*pBPatch1)
            KRATOS_WATCH(*pBPatch2)
        }

        bool is_valid = true;
        if (dist_tol == 0.0)
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0));
        else
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0, dist_tol));

        return is_valid;
    }

    /// Enumerate on the interface, i.e. to make sure that the enumeration on the two patch interfaces are compatible
    void Enumerate() override
    {
        std::vector<std::size_t> func_indices = this->pPatch1()->pFESpace()->ExtractBoundaryFunctionIndices(this->Side1());
        if (mDirection == BoundaryDirection::_REVERSED_)
        {
            std::reverse(func_indices.begin(), func_indices.end());
        }
        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices, false);
    }

    /// Return information as string
    std::string Info() const override
    {
        return "BSplinesPatchInterface2D";
    }

private:

    BoundaryDirection mDirection;
    // this variable indicates if the relative local direction shall be reversed of each other.
};

/**
 * Template specialization for 3D
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class BSplinesPatchInterface<3, TLocalCoordinateType, TCoordinateType, TDataType> : public PatchInterface<3, TLocalCoordinateType, TCoordinateType, TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<3, TLocalCoordinateType, TCoordinateType, TDataType> BaseType;
    typedef BSplinesPatchInterface<3, TLocalCoordinateType, TCoordinateType, TDataType> ThisType;
    typedef typename BaseType::PatchType PatchType;

    /// Empty Constructor
    BSplinesPatchInterface() : BaseType()
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mLocalParameterMap[0] = 0;
        mLocalParameterMap[1] = 1;
        mDirections[0] = BoundaryDirection::_FORWARD_;
        mDirections[1] = BoundaryDirection::_FORWARD_;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 2D.
    /// The direction variable is to indicate if the local parameters shall be reversed.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const BoundaryDirection direction)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mLocalParameterMap[0] = 0;
        mDirections[0] = direction;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 3D.
    /// if uv_or_vu is true, the local configuration {(u:ub) - (v:vb)} is assumed; otherwise, it is {(u:vb) - (v:ub)}
    /// The direction indicates if the respective local parameters shall be reversed
    /// In principal, the boundary face in 3D has local parameters u,v. The local direction u (or v)
    /// corresponds to which parameter direction of the 3D patch depends on how the boundary patch
    /// is constructed. The local parameter direction (uv) of one patch may align with those of
    /// the other patch. In that case uv_or_vu == true, otherwise uv_or_vu == false. The correct
    /// rotation is then determined by the direction information in each local parameter direction.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const bool uv_or_vu,
                           const BoundaryDirection direction1,
                           const BoundaryDirection direction2)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
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
    ~BSplinesPatchInterface() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "BSplinesPatchInterface3D, Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Create a clone of this interface
    typename BaseType::Pointer Clone() const override
    {
        const bool uv_or_vu = (this->LocalParameterMapping(0) == 0);
        if (uv_or_vu)
        {
            return typename BaseType::Pointer(new ThisType(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), true, this->Direction(0), this->Direction(1)));
        }
        else
        {
            return typename BaseType::Pointer(new ThisType(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), false, this->Direction(0), this->Direction(1)));
        }
    }

    /// Get the local parameter space mapping
    std::size_t LocalParameterMapping(std::size_t dim) const
    {
        std::map<std::size_t, std::size_t>::const_iterator it = mLocalParameterMap.find(dim);
        if (it != mLocalParameterMap.end())
        {
            return it->second;
        }
        else
            KRATOS_ERROR << "The dimension " << dim << " is invalid";
        }

    /// Get the relative direction between the local parameter space
    BoundaryDirection Direction(std::size_t dim) const {return mDirections[dim];}

    /// Validate the compatibility of two patches on the interface
    bool Validate(const bool debug, const double dist_tol) const override
    {
        typename PatchType::BoundaryPatchType::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());

        std::vector<BoundaryDirection> directions = {mDirections[LocalParameterMapping(0)], mDirections[LocalParameterMapping(1)]};

        typename PatchType::BoundaryPatchType::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2(), mLocalParameterMap, directions);

        if (debug)
        {
            auto pBFESpace1 = this->pPatch1()->pFESpace()->ConstructBoundaryFESpace(this->Side1());
            auto pBFESpace2 = this->pPatch2()->pFESpace()->ConstructBoundaryFESpace(this->Side2(), mLocalParameterMap, directions);
            KRATOS_WATCH(*pBFESpace1)
            KRATOS_WATCH(*pBFESpace2)

            std::cout << "Comparing two boundary patches on interface " << this << std::endl;
            KRATOS_WATCH(*pBPatch1)
            KRATOS_WATCH(*pBPatch2)
        }

        bool is_valid = true;
        if (dist_tol == 0.0)
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0));
        else
            is_valid = is_valid && (pBPatch1->IsEquivalent(*pBPatch2, debug ? 1 : 0, dist_tol));

        return is_valid;
    }

    /// Enumerate on the interface, i.e. to make sure that the enumeration on the two patch interfaces are compatible
    void Enumerate() override
    {
        std::vector<std::size_t> size_info;
        std::vector<std::size_t> func_indices = this->pPatch1()->pFESpace()->ExtractBoundaryFunctionIndices(size_info, this->Side1());

        const bool uv_or_vu = (this->LocalParameterMapping(0) == 0);

        // KRATOS_WATCH_STD_CON(func_indices)
        // BSplinesIndexingUtility::Transform(func_indices, size_info, uv_or_vu, mDirections[0], mDirections[1]);
        BSplinesIndexingUtility::Transform(func_indices, size_info, uv_or_vu, mDirections[this->LocalParameterMapping(0)], mDirections[this->LocalParameterMapping(1)]);

        // std::vector<std::size_t> size_info2;
        // std::vector<std::size_t> other_func_indices = this->pPatch2()->pFESpace()->ExtractBoundaryFunctionIndices(size_info2, this->Side2());

        // std::cout << "Patch " << this->pPatch1()->Id() << " transfers boundary function indices to Patch " << this->pPatch2()->Id() << ":";
        // KRATOS_WATCH_STD_CON(func_indices)

        // std::cout << "Current boundary function indices of Patch " << this->pPatch2()->Id() << ":";
        // KRATOS_WATCH_STD_CON(other_func_indices)

        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices, false);

        // other_func_indices = this->pPatch2()->pFESpace()->ExtractBoundaryFunctionIndices(size_info2, this->Side2());
        // std::cout << "New boundary function indices of Patch " << this->pPatch2()->Id() << ":";
        // KRATOS_WATCH_STD_CON(other_func_indices)
        // KRATOS_WATCH(uv_or_vu)
        // KRATOS_WATCH(mDirections[0])
        // KRATOS_WATCH(mDirections[1])
        // const auto pBFESpace1 = this->pPatch1()->pFESpace()->ConstructBoundaryFESpace(this->Side1());
        // KRATOS_WATCH(*pBFESpace1)
        // std::vector<BoundaryDirection> directions = {mDirections[LocalParameterMapping(0)], mDirections[LocalParameterMapping(1)]};
        // const auto pBFESpace2 = this->pPatch2()->pFESpace()->ConstructBoundaryFESpace(this->Side2(), mLocalParameterMap, directions);
        // KRATOS_WATCH(*pBFESpace2);
        // KRATOS_WATCH("-------------------")
        // KRATOS_WATCH("-------------------")
        // KRATOS_WATCH("-------------------")
        // KRATOS_WATCH("")
        // KRATOS_WATCH("")
    }

    /// Return information as string
    std::string Info() const override
    {
        return "BSplinesPatchInterface3D";
    }

private:

    std::map<std::size_t, std::size_t> mLocalParameterMap;
    // this variable stores the information about the local mapping at the interface. In 2D, it is simply 0:0.
    // In 3D, it shall have the values {0:0, 1:1} or {0:1, 1:0}. The latter values corresponds to the {(u:vb) - (v:ub)} configuration.

    boost::array<BoundaryDirection, 2> mDirections;
    // this variable indicates if the relative local direction shall be reversed of each other.
};

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED defined
