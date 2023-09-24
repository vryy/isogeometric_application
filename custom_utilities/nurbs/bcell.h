//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/cell.h"
#include "custom_utilities/nurbs/knot_array_1d.h"

namespace Kratos
{

class BCell;

template<int TDim>
struct BCell_Helper
{
    static bool IsCovered(const BCell& this_cell, const BCell& other_cell);
};

/**
 * Abstract class for a cell in NURBS/hierarchical B-Splines/T-splines mesh topology.
 * A cell is the smaller unit in the isogeometric topology mesh (e.g. it represents an element, or Bezier element of the T-splines basis function).
 * A cell is determined by its topology index of its vertices.
 */
class BCell : public Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BCell);
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const BCell> ConstPointer;
    #endif

    /// Type definitions
    typedef Cell BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef typename knot_container_t::KnotType KnotType;

    /// Constructor with knots
    BCell(std::size_t Id, knot_t pXiMin, knot_t pXiMax)
    : BaseType(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(new KnotType(0.0)), mpEtaMin(new KnotType(0.0)), mpZetaMax(new KnotType(0.0)), mpZetaMin(new KnotType(0.0))
    {}

    /// Constructor with knots
    BCell(std::size_t Id, knot_t pXiMin, knot_t pXiMax, knot_t pEtaMin, knot_t pEtaMax)
    : BaseType(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(pEtaMax), mpEtaMin(pEtaMin), mpZetaMax(new KnotType(0.0)), mpZetaMin(new KnotType(0.0))
    {}

    /// Constructor with knots
    BCell(std::size_t Id, knot_t pXiMin, knot_t pXiMax, knot_t pEtaMin, knot_t pEtaMax, knot_t pZetaMin, knot_t pZetaMax)
    : BaseType(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(pEtaMax), mpEtaMin(pEtaMin), mpZetaMax(pZetaMax), mpZetaMin(pZetaMin)
    {}

    /// Destructor
    virtual ~BCell() {}

    /// Get the level of this cell
    virtual std::size_t Level() const {return 1;}

    /// Get the coordinates
    knot_t XiMin() const {return mpXiMin;}
    int XiMinIndex() const {return mpXiMin->Index();}
    double XiMinValue() const {return mpXiMin->Value();}

    knot_t XiMax() const {return mpXiMax;}
    int XiMaxIndex() const {return mpXiMax->Index();}
    double XiMaxValue() const {return mpXiMax->Value();}

    knot_t EtaMax() const {return mpEtaMax;}
    int EtaMaxIndex() const {return mpEtaMax->Index();}
    double EtaMaxValue() const {return mpEtaMax->Value();}

    knot_t EtaMin() const {return mpEtaMin;}
    int EtaMinIndex() const {return mpEtaMin->Index();}
    double EtaMinValue() const {return mpEtaMin->Value();}

    knot_t ZetaMax() const {return mpZetaMax;}
    int ZetaMaxIndex() const {return mpZetaMax->Index();}
    double ZetaMaxValue() const {return mpZetaMax->Value();}

    knot_t ZetaMin() const {return mpZetaMin;}
    int ZetaMinIndex() const {return mpZetaMin->Index();}
    double ZetaMinValue() const {return mpZetaMin->Value();}

    /// Check if the cell is covered by knot spans; the comparison is based on indexing, so the knot vectors must be sorted a priori
    template<typename TIndexType>
    bool IsCovered(const std::vector<TIndexType>& rKnotsIndex1) const
    {
        TIndexType anchor_cover_xi_min  = *std::min_element(rKnotsIndex1.begin(), rKnotsIndex1.end());
        TIndexType anchor_cover_xi_max  = *std::max_element(rKnotsIndex1.begin(), rKnotsIndex1.end());

        if(     XiMinIndex()  >= anchor_cover_xi_min
            && XiMaxIndex() <= anchor_cover_xi_max  )
        {
            return true;
        }
        return false;
    }

    /// Check if the cell is covered by knot spans; the comparison is based on indexing, so the knot vectors must be sorted a priori
    template<typename TIndexType>
    bool IsCovered(const std::vector<TIndexType>& rKnotsIndex1, const std::vector<TIndexType>& rKnotsIndex2) const
    {
        TIndexType anchor_cover_xi_min  = *std::min_element(rKnotsIndex1.begin(), rKnotsIndex1.end());
        TIndexType anchor_cover_xi_max  = *std::max_element(rKnotsIndex1.begin(), rKnotsIndex1.end());
        TIndexType anchor_cover_eta_min = *std::min_element(rKnotsIndex2.begin(), rKnotsIndex2.end());
        TIndexType anchor_cover_eta_max = *std::max_element(rKnotsIndex2.begin(), rKnotsIndex2.end());

        if(     XiMinIndex()  >= anchor_cover_xi_min
            && XiMaxIndex() <= anchor_cover_xi_max
            && EtaMinIndex()  >= anchor_cover_eta_min
            && EtaMaxIndex()    <= anchor_cover_eta_max )
        {
            return true;
        }
        return false;
    }

    /// Check if the cell is covered by knot spans; the comparison is based on indexing, so the knot vectors must be sorted a priori
    template<typename TIndexType>
    bool IsCovered(const std::vector<TIndexType>& rKnotsIndex1, const std::vector<TIndexType>& rKnotsIndex2,
            const std::vector<TIndexType>& rKnotsIndex3) const
    {
        TIndexType anchor_cover_xi_min  = *std::min_element(rKnotsIndex1.begin(), rKnotsIndex1.end());
        TIndexType anchor_cover_xi_max  = *std::max_element(rKnotsIndex1.begin(), rKnotsIndex1.end());
        TIndexType anchor_cover_eta_min = *std::min_element(rKnotsIndex2.begin(), rKnotsIndex2.end());
        TIndexType anchor_cover_eta_max = *std::max_element(rKnotsIndex2.begin(), rKnotsIndex2.end());
        TIndexType anchor_cover_zeta_min = *std::min_element(rKnotsIndex3.begin(), rKnotsIndex3.end());
        TIndexType anchor_cover_zeta_max = *std::max_element(rKnotsIndex3.begin(), rKnotsIndex3.end());

        if(     XiMinIndex()  >= anchor_cover_xi_min
            && XiMaxIndex() <= anchor_cover_xi_max
            && EtaMinIndex()  >= anchor_cover_eta_min
            && EtaMaxIndex()    <= anchor_cover_eta_max
            && ZetaMinIndex()  >= anchor_cover_zeta_min
            && ZetaMaxIndex()    <= anchor_cover_zeta_max )
        {
            return true;
        }
        return false;
    }

    /// Check if this cell is covered by another cell
    template<int TDim>
    bool IsCovered(BCell::ConstPointer p_cell) const
    {
        return BCell_Helper<TDim>::IsCovered(*this, *p_cell);
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(double rXi, double rEta) const
    {
        if(    XiMinValue()  <= rXi  && XiMaxValue()  >= rXi
            && EtaMinValue() <= rEta && EtaMaxValue() >= rEta )
            return true;
        return false;
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(double rXi, double rEta, double rZeta) const
    {
        if(    XiMinValue()   <= rXi   && XiMaxValue()   >= rXi
            && EtaMinValue()  <= rEta  && EtaMaxValue()  >= rEta
            && ZetaMinValue() <= rZeta && ZetaMaxValue() >= rZeta )
            return true;
        return false;
    }

    /// check if this cell is the same as the reference cell. Two cells are the same if it has the same bounding knot values.
    bool IsSame(BCell::ConstPointer p_cell, double tol) const
    {
        if(    fabs( XiMinValue()   - p_cell->XiMinValue() )   < tol
            && fabs( XiMaxValue()   - p_cell->XiMaxValue() )   < tol
            && fabs( EtaMinValue()  - p_cell->EtaMinValue()  ) < tol
            && fabs( EtaMaxValue()  - p_cell->EtaMaxValue()  ) < tol
            && fabs( ZetaMinValue() - p_cell->ZetaMinValue() ) < tol
            && fabs( ZetaMaxValue() - p_cell->ZetaMaxValue() ) < tol )
                return true;
        return false;
    }

    /// Helper function to get out the value of the knot
    static double GetValue(knot_t p_knot)
    {
        return p_knot->Value();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "{Id:" << Id() << ",range:([" << XiMinIndex() << " " << XiMaxIndex() << "];[" << EtaMinIndex() << " " << EtaMaxIndex() << "];[" << ZetaMinIndex() << " " << ZetaMaxIndex() << "])";
        rOStream << "<=>([" << XiMinValue() << " " << XiMaxValue() << "];[" << EtaMinValue() << " " << EtaMaxValue() << "];[" << ZetaMinValue() << " " << ZetaMaxValue() << "])}";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }

private:

    knot_t mpXiMin;
    knot_t mpXiMax;
    knot_t mpEtaMax;
    knot_t mpEtaMin;
    knot_t mpZetaMax;
    knot_t mpZetaMin;
};

template<>
inline bool BCell_Helper<1>::IsCovered(const BCell& this_cell, const BCell& other_cell)
{
    if(    this_cell.XiMinValue() >= other_cell.XiMinValue() && this_cell.XiMaxValue() <= other_cell.XiMaxValue() )
        return true;
    return false;
}

template<>
inline bool BCell_Helper<2>::IsCovered(const BCell& this_cell, const BCell& other_cell)
{
    if(    this_cell.XiMinValue()  >= other_cell.XiMinValue()  && this_cell.XiMaxValue()  <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue() >= other_cell.EtaMinValue() && this_cell.EtaMaxValue() <= other_cell.EtaMaxValue() )
        return true;
    return false;
}

template<>
inline bool BCell_Helper<3>::IsCovered(const BCell& this_cell, const BCell& other_cell)
{
    if(    this_cell.XiMinValue()   >= other_cell.XiMinValue()   && this_cell.XiMaxValue()   <= other_cell.XiMaxValue()
        && this_cell.EtaMinValue()  >= other_cell.EtaMinValue()  && this_cell.EtaMaxValue()  <= other_cell.EtaMaxValue()
        && this_cell.ZetaMinValue() >= other_cell.ZetaMinValue() && this_cell.ZetaMaxValue() <= other_cell.ZetaMaxValue() )
        return true;
    return false;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BCell& rThis)
{
    rOStream << "bcell ";
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_H_INCLUDED

