//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 Nov 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TCELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TCELL_H_INCLUDED

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

class TCell;

template<int TDim>
struct TCell_Helper
{
    static bool IsCovered(const TCell& this_cell, const TCell& other_cell);
};

/**
 * Abstract class for a cell in NURBS/hierarchical B-Splines/T-splines mesh topology.
 * A cell is the smaller unit in the isogeometric topology mesh (e.g. it represents an element, or Bezier element of the T-splines basis function).
 * A cell is determined by its topology index of its vertices.
 */
class TCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TCell);

    /// Type definitions
    typedef Cell BaseType;
    typedef double knot_t;
    typedef std::vector<double> knot_container_t;

    /// Constructor with knots
    TCell(const std::size_t& Id, double XiMin, double XiMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(0.0), mEtaMin(0.0), mZetaMax(0.0), mZetaMin(0.0)
    {}

    /// Constructor with knots
    TCell(const std::size_t& Id, double XiMin, double XiMax, double EtaMin, double EtaMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(EtaMax), mEtaMin(EtaMin), mZetaMax(0.0), mZetaMin(0.0)
    {}

    /// Constructor with knots
    TCell(const std::size_t& Id, double XiMin, double XiMax, double EtaMin, double EtaMax, double ZetaMin, double ZetaMax)
    : BaseType(Id), mXiMin(XiMin), mXiMax(XiMax), mEtaMax(EtaMax), mEtaMin(EtaMin), mZetaMax(ZetaMax), mZetaMin(ZetaMin)
    {}

    /// Destructor
    virtual ~TCell() {}

    /// Get the coordinates
    double XiMin() const {return mXiMin;}
    double XiMinValue() const {return mXiMin;}

    double XiMax() const {return mXiMax;}
    double XiMaxValue() const {return mXiMax;}

    double EtaMin() const {return mEtaMin;}
    double EtaMinValue() const {return mEtaMin;}

    double EtaMax() const {return mEtaMax;}
    double EtaMaxValue() const {return mEtaMax;}

    double ZetaMin() const {return mZetaMin;}
    double ZetaMinValue() const {return mZetaMin;}

    double ZetaMax() const {return mZetaMax;}
    double ZetaMaxValue() const {return mZetaMax;}

    /// Check if this cell is covered by another cell
    template<int TDim>
    bool IsCovered(TCell::ConstPointer p_cell) const
    {
        return TCell_Helper<TDim>::IsCovered(*this, *p_cell);
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(const double& rXi, const double& rEta) const
    {
        if(    XiMinValue()  <= rXi  && XiMaxValue()  >= rXi
            && EtaMinValue() <= rEta && EtaMaxValue() >= rEta )
            return true;
        return false;
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(const double& rXi, const double& rEta, const double& rZeta) const
    {
        if(    XiMinValue()   <= rXi   && XiMaxValue()   >= rXi
            && EtaMinValue()  <= rEta  && EtaMaxValue()  >= rEta
            && ZetaMinValue() <= rZeta && ZetaMaxValue() >= rZeta )
            return true;
        return false;
    }

    /// check if this cell is the same as the reference cell. Two cells are the same if it has the same bounding knot values.
    bool IsSame(const TCell::Pointer p_cell, const double& tol) const
    {
        if(    fabs( XiMinValue()   - p_cell->XiMinValue()   ) < tol
            && fabs( XiMaxValue()   - p_cell->XiMaxValue()   ) < tol
            && fabs( EtaMinValue()  - p_cell->EtaMinValue()  ) < tol
            && fabs( EtaMaxValue()  - p_cell->EtaMaxValue()  ) < tol
            && fabs( ZetaMinValue() - p_cell->ZetaMinValue() ) < tol
            && fabs( ZetaMaxValue() - p_cell->ZetaMaxValue() ) < tol )
                return true;
        return false;
    }

    /// Helper function to get out the value of the knot
    static double GetValue(knot_t knot)
    {
        return knot;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "{Id:" << Id() << ",range:([" << XiMinValue() << " " << XiMaxValue() << "];[" << EtaMinValue() << " " << EtaMaxValue() << "];[" << ZetaMinValue() << " " << ZetaMaxValue() << "])}";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }

private:

    double mXiMin;
    double mXiMax;
    double mEtaMax;
    double mEtaMin;
    double mZetaMax;
    double mZetaMin;
};

template<>
inline bool TCell_Helper<1>::IsCovered(const TCell& this_cell, const TCell& other_cell)
{
    if(    this_cell.XiMinValue() >= other_cell.XiMinValue() && this_cell.XiMaxValue() <= other_cell.XiMaxValue() )
        return true;
    return false;
}

template<>
inline bool TCell_Helper<2>::IsCovered(const TCell& this_cell, const TCell& other_cell)
{
    if(    this_cell.XiMinValue()  >= other_cell.XiMinValue()  && this_cell.XiMaxValue()  <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue() >= other_cell.EtaMinValue() && this_cell.EtaMaxValue() <= other_cell.EtaMaxValue() )
        return true;
    return false;
}

template<>
inline bool TCell_Helper<3>::IsCovered(const TCell& this_cell, const TCell& other_cell)
{
    if(    this_cell.XiMinValue()   >= other_cell.XiMinValue()   && this_cell.XiMaxValue()   <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue()  >= other_cell.EtaMinValue()  && this_cell.EtaMaxValue()  <= other_cell.EtaMaxValue()
        && this_cell.ZetaMinValue() >= other_cell.ZetaMinValue() && this_cell.ZetaMaxValue() <= other_cell.ZetaMaxValue() )
        return true;
    return false;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TCell& rThis)
{
    rOStream << "cell ";
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TCELL_H_INCLUDED

