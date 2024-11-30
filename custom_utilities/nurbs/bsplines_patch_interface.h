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
#include "custom_utilities/nurbs/bsplines_patch_utility.h"

namespace Kratos
{

/**
 * This class represents an interface connecting two patches.
 * The interface is designed using half-edge philosophy, in which the interface keeps a pointer to the interface of the other side.
 */
template<int TDim>
class BSplinesPatchInterface : public PatchInterface<TDim>
{
};

/**
 * Template specialization for 1D
 */
template<>
class BSplinesPatchInterface<1> : public PatchInterface<1>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<1> BaseType;
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
        return typename BaseType::Pointer(new BSplinesPatchInterface<1>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2()));
    }

    /// Get the local parameter space mapping
    std::size_t LocalParameterMapping(std::size_t dim) const
    {
        return 0;
    }

    /// Get the relative direction between the local parameter space
    BoundaryDirection Direction() const {return _UNDEFINED_DIR_;}
    BoundaryDirection Direction(std::size_t dim) const {return _UNDEFINED_DIR_;}

    /// Validate the compatibility of two patches on the interface
    bool Validate(const bool debug, const double dist_tol) const override
    {
        typename Patch<0>::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename Patch<0>::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

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
};

/**
 * Template specialization for 2D
 */
template<>
class BSplinesPatchInterface<2> : public PatchInterface<2>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<2> BaseType;
    typedef typename BaseType::PatchType PatchType;

    /// Empty Constructor
    BSplinesPatchInterface() : BaseType()
    {}

    /// Full Constructor, the default local configuration is assumed and there is no reverse.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mDirection = _FORWARD_;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 2D.
    /// The direction variable is to indicate if the local parameters shall be reversed.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const BoundaryDirection& direction)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mDirection = direction;
    }

    /// Destructor
    virtual ~BSplinesPatchInterface()
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "BSplinesPatchInterface2D, Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Create a clone of this interface
    typename BaseType::Pointer Clone() const override
    {
        return typename BaseType::Pointer(new BSplinesPatchInterface<2>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), this->Direction()));
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
            return _UNDEFINED_DIR_;
        }
    }

    /// Validate the compatibility of two patches on the interface
    bool Validate(const bool debug, const double dist_tol) const override
    {
        typename Patch<1>::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename Patch<1>::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

        if (mDirection == _REVERSED_)
        {
            BSplinesPatchUtility::Reverse<1>(pBPatch2, 0);
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
        if (mDirection == _REVERSED_)
        {
            std::reverse(func_indices.begin(), func_indices.end());
        }
        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices);
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        BaseType::PrintInfo(rOStream);
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    BoundaryDirection mDirection;
    // this variable indicates if the relative local direction shall be reversed of each other.
};

/**
 * Template specialization for 3D
 */
template<>
class BSplinesPatchInterface<3> : public PatchInterface<3>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchInterface);

    typedef PatchInterface<3> BaseType;
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
        mDirections[0] = _FORWARD_;
        mDirections[1] = _FORWARD_;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 2D.
    /// The direction variable is to indicate if the local parameters shall be reversed.
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const BoundaryDirection& direction)
        : BaseType(pPatch1, side1, pPatch2, side2)
    {
        mLocalParameterMap[0] = 0;
        mDirections[0] = direction;
    }

    /// Full Constructor with the relative direction. This constructor is only valid in 3D.
    /// if uv_or_vu is true, the local configuration {(u:ub) - (v:vb)} is assumed; otherwise, it is {(u:vb) - (v:ub)}
    /// The direction indicates if the respective local parameters shall be reversed or not
    BSplinesPatchInterface(typename PatchType::Pointer pPatch1, const BoundarySide side1,
                           typename PatchType::Pointer pPatch2, const BoundarySide side2,
                           const bool& uv_or_vu,
                           const BoundaryDirection& direction1,
                           const BoundaryDirection& direction2)
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
    virtual ~BSplinesPatchInterface()
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "BSplinesPatchInterface3D, Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Create a clone of this interface
    typename BaseType::Pointer Clone() const override
    {
        if (this->LocalParameterMapping(0) == 0)
        {
            return typename BaseType::Pointer(new BSplinesPatchInterface<3>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), true, this->Direction(0), this->Direction(1)));
        }
        else if (this->LocalParameterMapping(0) == 1)
        {
            return typename BaseType::Pointer(new BSplinesPatchInterface<3>(this->pPatch1(), this->Side1(), this->pPatch2(), this->Side2(), false, this->Direction(0), this->Direction(1)));
        }
        else
            KRATOS_ERROR << "Invalid local parameter mapping";
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
        typename Patch<2>::Pointer pBPatch1 = this->pPatch1()->ConstructBoundaryPatch(this->Side1());
        typename Patch<2>::Pointer pBPatch2 = this->pPatch2()->ConstructBoundaryPatch(this->Side2());

        if (mDirections[0] == _REVERSED_)
        {
            BSplinesPatchUtility::Reverse<2>(pBPatch2, 0);
        }
        if (mDirections[1] == _REVERSED_)
        {
            BSplinesPatchUtility::Reverse<2>(pBPatch2, 1);
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
        std::vector<std::size_t> size_info;
        std::vector<std::size_t> func_indices = this->pPatch1()->pFESpace()->ExtractBoundaryFunctionIndices(size_info, this->Side1());

        if (mDirections[0] == _REVERSED_)
        {
            for (std::size_t i = 0; i < size_info[0]; ++i)
            {
                std::vector<std::size_t> tmp(size_info[1]);
                for (std::size_t j = 0; j < size_info[1]; ++j)
                {
                    tmp[j] = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, size_info[1] - j, size_info[0], size_info[1])];
                }
                for (std::size_t j = 0; j < size_info[1]; ++j)
                {
                    func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, size_info[0], size_info[1])] = tmp[j];
                }
            }
        }

        if (mDirections[1] == _REVERSED_)
        {
            for (std::size_t j = 0; j < size_info[1]; ++j)
            {
                std::vector<std::size_t> tmp(size_info[0]);
                for (std::size_t i = 0; i < size_info[0]; ++i)
                {
                    tmp[i] = func_indices[BSplinesIndexingUtility_Helper::Index2D(size_info[0] - i, j + 1, size_info[0], size_info[1])];
                }
                for (std::size_t i = 0; i < size_info[0]; ++i)
                {
                    func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, size_info[0], size_info[1])] = tmp[i];
                }
            }
        }

        this->pPatch2()->pFESpace()->AssignBoundaryFunctionIndices(this->Side2(), func_indices);
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        BaseType::PrintInfo(rOStream);
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    std::map<std::size_t, std::size_t> mLocalParameterMap;
    // this variable stores the information about the local mapping at the interface. In 2D, it is simply 0:0.
    // In 3D, it shall have the values {0:0, 1:1} or {0:1, 1:0}. The latter values corresponds to the {(u:vb) - (v:ub)} configuration.

    boost::array<BoundaryDirection, 2> mDirections;
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_INTERFACE_H_INCLUDED defined
