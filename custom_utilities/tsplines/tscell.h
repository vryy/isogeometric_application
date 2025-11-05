//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 Nov 2018 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TSCELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TSCELL_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/cell.h"


namespace Kratos
{

class TsCell;

template<int TDim>
struct TsCell_Helper
{
    static bool IsCovered(const TsCell& this_cell, const TsCell& other_cell);
};

/**
 * Abstract class for a cell in NURBS/hierarchical B-Splines/T-splines mesh topology.
 * A cell is the smaller unit in the isogeometric topology mesh (e.g. it represents an element, or Bezier element of the T-splines basis function).
 * A cell is determined by its topology index of its vertices.
 */
class TsCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsCell);
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const TsCell> ConstPointer;
    #endif

    /// Type definitions
    typedef Cell BaseType;
    typedef double knot_t;
    typedef std::vector<knot_t> knot_container_t;
    typedef BaseType::IndexType IndexType;

    /// Constructor with knots
    TsCell(const IndexType Id, knot_t XiMin, knot_t XiMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(0.0), mEtaMin(0.0), mZetaMax(0.0), mZetaMin(0.0)
    {}

    /// Constructor with knots
    TsCell(const IndexType Id, knot_t XiMin, knot_t XiMax, knot_t EtaMin, knot_t EtaMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(EtaMax), mEtaMin(EtaMin), mZetaMax(0.0), mZetaMin(0.0)
    {}

    /// Constructor with knots
    TsCell(const IndexType Id, knot_t XiMin, knot_t XiMax, knot_t EtaMin, knot_t EtaMax, knot_t ZetaMin, knot_t ZetaMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(EtaMax), mEtaMin(EtaMin), mZetaMax(ZetaMax), mZetaMin(ZetaMin)
    {}

    /// Destructor
    ~TsCell() override {}

    /// Get the level of this cell
    IndexType Level() const override
    {
        return 1;
    }

    /// Get the coordinates
    knot_t XiMin() const {return mXiMin;}
    knot_t XiMinValue() const {return mXiMin;}

    knot_t XiMax() const {return mXiMax;}
    knot_t XiMaxValue() const {return mXiMax;}

    knot_t EtaMin() const {return mEtaMin;}
    knot_t EtaMinValue() const {return mEtaMin;}

    knot_t EtaMax() const {return mEtaMax;}
    knot_t EtaMaxValue() const {return mEtaMax;}

    knot_t ZetaMin() const {return mZetaMin;}
    knot_t ZetaMinValue() const {return mZetaMin;}

    knot_t ZetaMax() const {return mZetaMax;}
    knot_t ZetaMaxValue() const {return mZetaMax;}

    /// Check if this cell is covered by another cell
    template<int TDim>
    bool IsCovered(TsCell::ConstPointer p_cell) const
    {
        return TsCell_Helper<TDim>::IsCovered(*this, *p_cell);
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(const knot_t rXi, const knot_t rEta) const
    {
        if(    XiMinValue()  <= rXi  && XiMaxValue()  >= rXi
            && EtaMinValue() <= rEta && EtaMaxValue() >= rEta )
            return true;
        return false;
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(const knot_t rXi, const knot_t rEta, const knot_t rZeta) const
    {
        if(    XiMinValue()   <= rXi   && XiMaxValue()   >= rXi
            && EtaMinValue()  <= rEta  && EtaMaxValue()  >= rEta
            && ZetaMinValue() <= rZeta && ZetaMaxValue() >= rZeta )
            return true;
        return false;
    }

    /// check if this cell is the same as the reference cell. Two cells are the same if it has the same bounding knot values.
    bool IsSame(const TsCell::Pointer p_cell, const knot_t tol) const
    {
        if(    std::abs( XiMinValue()   - p_cell->XiMinValue()   ) < tol
            && std::abs( XiMaxValue()   - p_cell->XiMaxValue()   ) < tol
            && std::abs( EtaMinValue()  - p_cell->EtaMinValue()  ) < tol
            && std::abs( EtaMaxValue()  - p_cell->EtaMaxValue()  ) < tol
            && std::abs( ZetaMinValue() - p_cell->ZetaMinValue() ) < tol
            && std::abs( ZetaMaxValue() - p_cell->ZetaMaxValue() ) < tol )
                return true;
        return false;
    }

    /// Helper function to get out the value of the knot
    static double GetValue(knot_t knot)
    {
        return knot;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "{Id:" << Id() << ",range:([" << XiMinValue() << " " << XiMaxValue() << "];[" << EtaMinValue() << " " << EtaMaxValue() << "];[" << ZetaMinValue() << " " << ZetaMaxValue() << "])}";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:

    knot_t mXiMin;
    knot_t mXiMax;
    knot_t mEtaMax;
    knot_t mEtaMin;
    knot_t mZetaMax;
    knot_t mZetaMin;
};

template<>
inline bool TsCell_Helper<1>::IsCovered(const TsCell& this_cell, const TsCell& other_cell)
{
    if(    this_cell.XiMinValue() >= other_cell.XiMinValue() && this_cell.XiMaxValue() <= other_cell.XiMaxValue() )
        return true;
    return false;
}

template<>
inline bool TsCell_Helper<2>::IsCovered(const TsCell& this_cell, const TsCell& other_cell)
{
    if(    this_cell.XiMinValue()  >= other_cell.XiMinValue()  && this_cell.XiMaxValue()  <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue() >= other_cell.EtaMinValue() && this_cell.EtaMaxValue() <= other_cell.EtaMaxValue() )
        return true;
    return false;
}

template<>
inline bool TsCell_Helper<3>::IsCovered(const TsCell& this_cell, const TsCell& other_cell)
{
    if(    this_cell.XiMinValue()   >= other_cell.XiMinValue()   && this_cell.XiMaxValue()   <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue()  >= other_cell.EtaMinValue()  && this_cell.EtaMaxValue()  <= other_cell.EtaMaxValue()
        && this_cell.ZetaMinValue() >= other_cell.ZetaMinValue() && this_cell.ZetaMaxValue() <= other_cell.ZetaMaxValue() )
        return true;
    return false;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TSCELL_H_INCLUDED
