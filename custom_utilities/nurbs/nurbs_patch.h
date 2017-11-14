//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/control_point.h"
#include "custom_utilities/nurbs/grid_function.h"
#include "custom_utilities/nurbs/patch.h"

namespace Kratos
{

/**
This class represents a single NURBS patch in parametric coordinates.
 */
template<std::size_t TDim>
class NURBSPatch : public Patch<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSPatch);

    /// Type definition
    typedef Patch<TDim> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef typename BaseType::ControlPointType ControlPointType;

    /// Default constructor
    NURBSPatch() : BaseType() {}

    /// Constructor with id
    NURBSPatch(const std::size_t& Id) : BaseType(Id) {}

    /// Destructor
    virtual ~NURBSPatch() {}

    /// Get the order of the NURBS patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return mOrders[i];}

    /// Get the number of control points of the NURBSPatch in specific direction
    const std::size_t Number(const std::size_t& i) const {return mNumbers[i];}

    /// Get the number of basis functions defined over the NURBS NURBSPatch
    virtual const std::size_t TotalNumber() const
    {
        std::size_t Number = 1;
        for (std::size_t i = 0; i < TDim; ++i)
            Number *= mNumbers[i];
        return Number;
    }

    /// Get the string representing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "NURBSPatch" << TDim << "D";
        return ss.str();
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {
        mpKnotVectors[i] = p_knot_vector;
    }

    /// Create and set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const std::vector<double>& values)
    {
        if (i >= TDim)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
        }
        else
        {
            for (std::size_t j = 0; j < values.size(); ++j)
                mpKnotVectors[i].pCreateKnot(values[j]);
        }
    }

    /// Get the knot vector in i-direction
    const knot_container_t& KnotVector(const std::size_t& i) const {return mpKnotVectors[i];}

    /// Set the NURBS information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {
        mOrders[i] = Order;
        mNumbers[i] = Number;
    }

    /// Validate the NURBSPatch before using
    virtual bool Validate() const
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (mpKnotVectors[i].size() != mNumbers[i] + mOrders[i] + 1)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The knot vector is incompatible at dimension", i)
                return false;
            }
        }

        return BaseType::Validate();
    }

    /// Compare between two NURBS patches in terms of parametric information
    virtual bool IsCompatible(const Patch<TDim>& rOtherPatch) const
    {
        if (rOtherPatch.Type() != Type())
        {
            KRATOS_WATCH(rOtherPatch.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const NURBSPatch<TDim>& rOtherNURBSPatch = dynamic_cast<const NURBSPatch<TDim>&>(rOtherPatch);

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Number(i)) == rOtherNURBSPatch.Number(i))
                return false;
            if (!(this->Order(i)) == rOtherNURBSPatch.Order(i))
                return false;
            if (!(this->KnotVector(i) == rOtherNURBSPatch.KnotVector(i)))
                return false;
        }

        return true;
    }

    /// Construct the boundary patch based on side
    virtual typename Patch<TDim-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    {
        typename NURBSPatch<TDim-1>::Pointer pBPatch = typename NURBSPatch<TDim-1>::Pointer(new NURBSPatch<TDim-1>(-1));

        if (TDim == 2)
        {
            if (side == _LEFT_)
            {
                pBPatch->SetKnotVector(0, KnotVector(1));
                pBPatch->SetInfo(0, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _RIGHT_)
            {
                pBPatch->SetKnotVector(0, KnotVector(1));
                pBPatch->SetInfo(0, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _TOP_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetInfo(0, Number(0), Order(0));
                // TODO set the control points and grid functions
            }
            else if (side == _BOTTOM_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetInfo(0, Number(0), Order(0));
                // TODO set the control points and grid functions
            }
        }
        else if (TDim == 3)
        {
            if (side == _LEFT_)
            {
                pBPatch->SetKnotVector(0, KnotVector(1));
                pBPatch->SetKnotVector(1, KnotVector(2));
                pBPatch->SetInfo(0, Number(1), Order(1));
                pBPatch->SetInfo(1, Number(2), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _RIGHT_)
            {
                pBPatch->SetKnotVector(0, KnotVector(1));
                pBPatch->SetKnotVector(1, KnotVector(2));
                pBPatch->SetInfo(0, Number(1), Order(1));
                pBPatch->SetInfo(1, Number(2), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _TOP_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetKnotVector(1, KnotVector(1));
                pBPatch->SetInfo(0, Number(0), Order(0));
                pBPatch->SetInfo(1, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _BOTTOM_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetKnotVector(1, KnotVector(1));
                pBPatch->SetInfo(0, Number(0), Order(0));
                pBPatch->SetInfo(1, Number(1), Order(1));
                // TODO set the control points and grid functions
            }
            else if (side == _FRONT_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetKnotVector(1, KnotVector(2));
                pBPatch->SetInfo(0, Number(0), Order(0));
                pBPatch->SetInfo(1, Number(1), Order(2));
                // TODO set the control points and grid functions
            }
            else if (side == _BACK_)
            {
                pBPatch->SetKnotVector(0, KnotVector(0));
                pBPatch->SetKnotVector(1, KnotVector(2));
                pBPatch->SetInfo(0, Number(0), Order(0));
                pBPatch->SetInfo(1, Number(1), Order(2));
                // TODO set the control points and grid functions
            }
        }

        return pBPatch;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Id = " << BaseType::Id() << ", Add = " << this << ", n = (";
        for (std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << this->Number(i);
        rOStream << "), p = (";
        for (std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << this->Order(i);
        rOStream << ")";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            rOStream << " knot vector " << i << ":";
            for (std::size_t j = 0; j < mpKnotVectors[i].size(); ++j)
                rOStream << " " << mpKnotVectors[i].pKnotAt(j)->Value();
            rOStream << std::endl;
        }
        if (BaseType::ControlPointGrid() != NULL)
            rOStream << *(BaseType::ControlPointGrid());
        BaseType::PrintData(rOStream);
    }

private:

    /**
     * internal data to construct the shape functions on the NURBS
     */
    boost::array<std::size_t, TDim> mOrders;
    boost::array<std::size_t, TDim> mNumbers;
    boost::array<knot_container_t, TDim> mpKnotVectors;
};

/**
 * Template specific instantiation for null-D NURBS patch to terminate the compilation
 */
template<>
class NURBSPatch<0> : public Patch<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(NURBSPatch);

    /// Type definition
    typedef Patch<0> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    NURBSPatch() : BaseType() {}

    /// Constructor with id
    NURBSPatch(const std::size_t& Id) : BaseType(Id) {}

    /// Destructor
    virtual ~NURBSPatch() {}

    /// Get the order of the NURBS patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return 0;}

    /// Get the number of basis functions defined over the NURBS NURBSPatch
    virtual const std::size_t Number() const {return 0;}

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "NURBSPatch0D";
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(const std::size_t& i, const knot_container_t& p_knot_vector)
    {}

    /// Set the NURBS information in the direction i
    void SetInfo(const std::size_t& i, const std::size_t& Number, const std::size_t& Order)
    {}

    /// Validate the NURBSPatch before using
    virtual bool Validate() const
    {
        return BaseType::Validate();
    }

    /// Compare between two NURBS patches in terms of parametric information
    virtual bool IsCompatible(const Patch<0>& rOtherPatch) const
    {
        if (rOtherPatch.Type() != Type())
        {
            KRATOS_WATCH(rOtherPatch.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        return true;
    }
};

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const NURBSPatch<TDim>& rThis)
{
    rOStream << "-------------Begin NURBS PatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << "-------------End NURBS PatchInfo-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_NURBS_PATCH_H_INCLUDED defined

