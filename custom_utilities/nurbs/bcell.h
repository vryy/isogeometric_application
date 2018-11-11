//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <boost/numeric/ublas/vector_sparse.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"

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
class BCell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BCell);

    /// Type definitions
    typedef Knot<double> KnotType;
    typedef KnotType::Pointer knot_t;
    typedef boost::numeric::ublas::mapped_vector<double> SparseVectorType;
    // typedef boost::numeric::ublas::vector<double> SparseVectorType;

    /// Constructor with knots
    BCell(const std::size_t& Id, knot_t pXiMin, knot_t pXiMax)
    : mId(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(new KnotType(0.0)), mpEtaMin(new KnotType(0.0)), mpZetaMax(new KnotType(0.0)), mpZetaMin(new KnotType(0.0))
    {}

    /// Constructor with knots
    BCell(const std::size_t& Id, knot_t pXiMin, knot_t pXiMax, knot_t pEtaMin, knot_t pEtaMax)
    : mId(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(pEtaMax), mpEtaMin(pEtaMin), mpZetaMax(new KnotType(0.0)), mpZetaMin(new KnotType(0.0))
    {}

    /// Constructor with knots
    BCell(const std::size_t& Id, knot_t pXiMin, knot_t pXiMax, knot_t pEtaMin, knot_t pEtaMax, knot_t pZetaMin, knot_t pZetaMax)
    : mId(Id), mpXiMin(pXiMin), mpXiMax(pXiMax), mpEtaMax(pEtaMax), mpEtaMin(pEtaMin), mpZetaMax(pZetaMax), mpZetaMin(pZetaMin)
    {}

    /// Destructor
    virtual ~BCell() {}

    /// Get the Id
    std::size_t Id() const {return mId;}

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
    bool IsCoverage(const double& rXi, const double& rEta) const
    {
        if(    XiMinValue()  <= rXi  && XiMaxValue() >= rXi
            && EtaMinValue()  <= rEta && EtaMaxValue()    >= rEta )
            return true;
        return false;
    }

    /// Check if this cell cover a point in knot space
    bool IsCoverage(const double& rXi, const double& rEta, const double& rZeta) const
    {
        if(    XiMinValue()  <= rXi   && XiMaxValue() >= rXi
            && EtaMinValue()  <= rEta  && EtaMaxValue()    >= rEta
            && ZetaMinValue() <= rZeta && ZetaMaxValue() >= rZeta )
            return true;
        return false;
    }

    /// check if this cell is the same as the reference cell. Two cells are the same if it has the same bounding knot values.
    bool IsSame(const BCell::Pointer p_cell, const double& tol) const
    {
        if(    fabs( XiMinValue()  - p_cell->XiMinValue()  ) < tol
            && fabs( XiMaxValue() - p_cell->XiMaxValue() ) < tol
            && fabs( EtaMinValue()  - p_cell->EtaMinValue()  ) < tol
            && fabs( EtaMaxValue()    - p_cell->EtaMaxValue()    ) < tol
            && fabs( ZetaMinValue() - p_cell->ZetaMinValue() ) < tol
            && fabs( ZetaMaxValue() - p_cell->ZetaMaxValue() ) < tol )
                return true;
        return false;
    }

    /// Clear internal data of this cell
    void Reset()
    {
        mSupportedAnchors.clear();
        mAnchorWeights.clear();
        mCrows.clear();
    }

    /// Add supported anchor and the respective extraction operator of this cell to the anchor
    void AddAnchor(const std::size_t& Id, const double& W, const Vector& Crow)
    {
        mSupportedAnchors.push_back(Id);
        mAnchorWeights.push_back(W);

        std::vector<std::size_t> inz;
        inz.reserve(Crow.size());
        for (std::size_t i = 0; i < Crow.size(); ++i)
            if (Crow[i] != 0.0)
                inz.push_back(i);
        SparseVectorType Crow_sparse(Crow.size(), inz.size());
        for (std::size_t i = 0; i < inz.size(); ++i)
            Crow_sparse[inz[i]] = Crow[inz[i]];
        mCrows.push_back(Crow_sparse);
    }

    /// Absorb the information from the other cell
    virtual void Absorb(BCell::Pointer pOther)
    {
        for (std::size_t i = 0; i < pOther->NumberOfAnchors(); ++i)
        {
            if (std::find(mSupportedAnchors.begin(), mSupportedAnchors.end(), pOther->GetSupportedAnchors()[i]) == mSupportedAnchors.end())
            {
                this->AddAnchor(pOther->GetSupportedAnchors()[i], pOther->GetAnchorWeights()[i], pOther->GetCrows()[i]);
            }
        }
    }

    /// This action is called when the cell is removed from the cell manager, see e.g. hb_cell.
    virtual void ClearTrace()
    {
        // DO NOTHING
    }

    /// Get the number of supported anchors of this cell. In other words, it is the number of basis functions that the support domain includes this cell.
    std::size_t NumberOfAnchors() const {return mSupportedAnchors.size();}

    /// Get the supported anchors of this cell
    const std::vector<std::size_t>& GetSupportedAnchors() const {return mSupportedAnchors;}

    /// Get the weights of all the supported anchors
    const std::vector<double>& GetAnchorWeights() const {return mAnchorWeights;}
    void GetAnchorWeights(Vector& rWeights) const
    {
        if(rWeights.size() != mAnchorWeights.size())
            rWeights.resize(mAnchorWeights.size(), false);
        std::copy(mAnchorWeights.begin(), mAnchorWeights.end(), rWeights.begin());
    }

    /// Get the internal data of row of the extraction operator
    const std::vector<SparseVectorType>& GetCrows() const {return mCrows;}

    /// Get the extraction operator matrix
    Matrix GetExtractionOperator() const
    {
        Matrix M(mCrows.size(), mCrows[0].size());
        for(std::size_t i = 0; i < mCrows.size(); ++i)
            noalias(row(M, i)) = mCrows[i];
        return M;
    }

    /// Get the extraction as compressed matrix
    CompressedMatrix GetCompressedExtractionOperator() const
    {
        CompressedMatrix M(mCrows.size(), mCrows[0].size());
        // here we just naively copy the data. However we shall initialize the compressed matrix properly as in the kernel (TODO)
        for(std::size_t i = 0; i < mCrows.size(); ++i)
            noalias(row(M, i)) = mCrows[i];
        M.complete_index1_data();
        return M;
    }

    /// Get the extraction operator as CSR triplet
    void GetExtractionOperator(std::vector<int>& rowPtr, std::vector<int>& colInd, std::vector<double>& values) const
    {
        int cnt = 0;
        rowPtr.push_back(cnt);
        for(std::size_t i = 0; i < mCrows.size(); ++i)
        {
            for(std::size_t j = 0; j < mCrows[i].size(); ++j)
            {
                if(mCrows[i](j) != 0)
                {
                    colInd.push_back(j);
                    values.push_back(mCrows[i](j));
                    ++cnt;
                }
            }
            rowPtr.push_back(cnt);
        }
    }

    /// Implement relational operator for automatic arrangement in container
    inline bool operator==(const BCell& rA) const
    {
        return this->Id() == rA.Id();
    }

    inline bool operator<(const BCell& rA) const
    {
        return this->Id() < rA.Id();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "{Id:" << Id() << ",range:([" << XiMinIndex() << " " << XiMaxIndex() << "];[" << EtaMinIndex() << " " << EtaMaxIndex() << "];[" << ZetaMinIndex() << " " << ZetaMaxIndex() << "])";
        rOStream << "<=>([" << XiMinValue() << " " << XiMaxValue() << "];[" << EtaMinValue() << " " << EtaMaxValue() << "];[" << ZetaMinValue() << " " << ZetaMaxValue() << "])}";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ", supporting anchors: ";
        rOStream << "(";
        for(std::vector<std::size_t>::const_iterator it = mSupportedAnchors.begin(); it != mSupportedAnchors.end(); ++it)
            rOStream << " " << (*it);
        rOStream << ")";
    }

private:

    std::size_t mId;
    knot_t mpXiMin;
    knot_t mpXiMax;
    knot_t mpEtaMax;
    knot_t mpEtaMin;
    knot_t mpZetaMax;
    knot_t mpZetaMin;
    std::vector<std::size_t> mSupportedAnchors;
    std::vector<double> mAnchorWeights; // weight of the anchor
    std::vector<SparseVectorType> mCrows; // bezier extraction operator row to each anchor
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
    rOStream << "cell ";
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_H_INCLUDED

