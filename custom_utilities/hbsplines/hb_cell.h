//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 May 2015 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/nurbs/cell.h"

namespace Kratos
{

/**
    Represent a cell in hierarchical B-Splines topology
 */
template<class TBasisFuncType>
class HBCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBCell);

    /// Type definitions
    typedef Knot<double>::Pointer knot_t;

    typedef Cell BaseType;

    typedef typename TBasisFuncType::Pointer bf_t;
    typedef std::set<bf_t> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    /// Constructor for 1D
    HBCell(const std::size_t& Id, knot_t pLeft, knot_t pRight)
    : BaseType(Id, pLeft, pRight)
    {}

    /// Constructor for 2D
    HBCell(const std::size_t& Id, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp)
    : BaseType(Id, pLeft, pRight, pDown, pUp)
    {}

    /// Constructor for 3D
    HBCell(const std::size_t& Id, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp, knot_t pBelow, knot_t pAbove)
    : BaseType(Id, pLeft, pRight, pDown, pUp, pBelow, pAbove)
    {}

    /// Destructor
    virtual ~HBCell()
    {}

    /// Set the level for this cell
    void SetLevel(const std::size_t& Level) {mLevel = Level;}

    /// Get the level of this cell
    const std::size_t& Level() const {return mLevel;}

    /// Add the basis function to the set. If it does exist in the set, return the internal one.
    /// Typically one shall add the basis function that has support domain covering this cell.
    bf_t AddBf(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            if(*it == p_bf)
                return *it;
        mpBasisFuncs.insert(p_bf);
        return p_bf;
    }

    /// Remove basis function from the set
    void RemoveBf(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if(*it == p_bf)
            {
                mpBasisFuncs.erase(it);
                break;
            }
        }
    }

    /// Absorb the information from the other cell
    virtual void Absorb(BaseType::Pointer pOther)
    {
        BaseType::Absorb(pOther);

        try
        {
            typedef HBCell<TBasisFuncType> HBCellType;

            typename HBCellType::Pointer pOtherCell = boost::dynamic_pointer_cast<HBCellType>(pOther);
            if (pOtherCell == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBCell is failed.", "")
            for(typename HBCellType::bf_iterator it_bf = pOtherCell->bf_begin(); it_bf != pOtherCell->bf_end(); ++it_bf)
            {
                this->AddBf(*it_bf);
            }
        }
        catch (std::exception const& e)
        {
            // DO NOTHING
        }
    }

    /// Iterators
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_const_iterator bf_begin() const {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_const_iterator bf_end() const {return mpBasisFuncs.end();}

    /// Number of basis functions that covers this cell
    std::size_t size() const {return mpBasisFuncs.size();}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ", supporting basis functions: (";
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            rOStream << " " << (*it)->Id();
        rOStream << ")";
        BaseType::PrintData(rOStream);
    }

private:

    std::size_t mLevel;
    bf_container_t mpBasisFuncs; // list of basis functions contain this cell in its support
};

/// output stream function
template<class TBasisFuncType>
inline std::ostream& operator <<(std::ostream& rOStream, const HBCell<TBasisFuncType>& rThis)
{
    rOStream << "hbcell ";
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HB_CELL_H_INCLUDED

