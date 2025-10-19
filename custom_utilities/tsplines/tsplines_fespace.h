//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Nov 2018 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/nurbs/pbbsplines_fespace.h"

/*#define DEBUG_GEN_CELL*/

namespace Kratos
{

/**
 * Abstract class represents the FESpace for a single T-Splines patch.
 */
template<int TDim, typename TLocalCoordinateType, typename TBasisFunctionType, typename TCellManagerType>
class TSplinesFESpace : public PBBSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TSplinesFESpace);

    /// Type definition
    typedef PBBSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType> BaseType;
    typedef TSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType> ThisType;
    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;
    typedef typename BaseType::BasisFunctionType BasisFunctionType;
    typedef typename BaseType::bf_t bf_t;
    typedef typename BaseType::bf_container_t bf_container_t;
    typedef typename BaseType::bf_iterator bf_iterator;
    typedef typename BaseType::bf_const_iterator bf_const_iterator;

    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_container_t cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::knot_container_t knot_container_t;
    typedef typename BaseType::knot_t knot_t;

    typedef typename BaseType::function_map_t function_map_t;

    /// Default constructor
    TSplinesFESpace() : BaseType()
    {
        mpFaceManager = typename cell_container_t::Pointer(new cell_container_t());
    }

    /// Destructor
    ~TSplinesFESpace() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << Type() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper to create new TSplinesFESpace pointer
    static typename ThisType::Pointer Create()
    {
        return typename ThisType::Pointer(new ThisType());
    }

    /// Get the string representing the type of the patch
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "TSplinesFESpace" << TDim << "D";
        return ss.str();
    }

    /// Get the underlying face manager
    typename cell_container_t::Pointer pFaceManager() {return mpFaceManager;}

    /// Get the underlying face manager
    typename cell_container_t::ConstPointer pFaceManager() const {return mpFaceManager;}

    /// Validate the TSplinesFESpace
    bool Validate() const override
    {
        // TODO validate more
        return BaseType::Validate();
    }

    /// Compare between two BSplines patches in terms of parametric information
    bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const override
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const ThisType* pOtherTSplinesFESpace = dynamic_cast<const TSplinesFESpace*>(&rOtherFESpace);
        if (pOtherTSplinesFESpace == nullptr)
            return false;

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Order(i)) == pOtherTSplinesFESpace->Order(i))
            {
                return false;
            }
        }

        return true;
    }

    /// Overload assignment operator
    TSplinesFESpace& operator=(const TSplinesFESpace& rOther)
    {
        // TODO copy more
        KRATOS_ERROR << "The assignment operator is not yet fully implemented";
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    typename FESpaceType::Pointer Clone() const override
    {
        typename TSplinesFESpace::Pointer pNewFESpace = typename TSplinesFESpace::Pointer(new TSplinesFESpace());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Type() << ", Addr = " << this << ", n = " << this->TotalNumber();
        rOStream << ", p = (";
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            rOStream << " " << this->Order(dim);
        }
        rOStream << ")";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "###################" << std::endl;

        // print the basis functions
        rOStream << "Basis functions:" << std::endl;
        for (bf_const_iterator it = BaseType::bf_begin(); it != BaseType::bf_end(); ++it)
        {
            rOStream << " >> " << *(*it) << std::endl;
        }

        // print the cells
        rOStream << "Cells:" << std::endl;
        for (typename cell_container_t::iterator it = BaseType::mpCellManager->begin(); it != BaseType::mpCellManager->end(); ++it)
        {
            rOStream << " >> " << *(*it) << std::endl;
        }
    }

protected:

    typename cell_container_t::Pointer mpFaceManager;

};

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TBasisFunctionType, typename TCellManagerType>
inline std::ostream& operator <<(std::ostream& rOStream, const TSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType>& rThis)
{
    rOStream << "-------------Begin TSplinesFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End TSplinesFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_GEN_CELL    /// Get the underlying cell manager

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_FESPACE_H_INCLUDED defined
