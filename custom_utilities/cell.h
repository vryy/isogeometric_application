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

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

// External includes
#include <boost/numeric/ublas/vector_sparse.hpp>


namespace Kratos
{

/**
 * Abstract class for a cell
 */
class Cell
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Cell);

    /// Type definitions
    typedef boost::numeric::ublas::mapped_vector<double> SparseVectorType;
    // typedef boost::numeric::ublas::vector<double> SparseVectorType;

    /// Default constructor
    Cell(const std::size_t& Id) : mId(Id)
    {}

    /// Destructor
    virtual ~Cell()
    {}

    /// Get the Id
    std::size_t Id() const {return mId;}

    /// Clear internal data of this cell
    virtual void Reset()
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
    virtual void Absorb(Cell::Pointer pOther)
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
    inline bool operator==(const Cell& rA) const
    {
        return this->Id() == rA.Id();
    }

    inline bool operator<(const Cell& rA) const
    {
        return this->Id() < rA.Id();
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "{Id:" << Id() << "}";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << ", supporting anchors: ";
        rOStream << "(";
        for(std::vector<std::size_t>::const_iterator it = mSupportedAnchors.begin(); it != mSupportedAnchors.end(); ++it)
            rOStream << " " << (*it);
        rOStream << ")";
    }

protected:

    std::size_t mId;
    std::vector<std::size_t> mSupportedAnchors;
    std::vector<double> mAnchorWeights; // weight of the anchor
    std::vector<SparseVectorType> mCrows; // bezier extraction operator row to each anchor
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const Cell& rThis)
{
    rOStream << "cell ";
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_H_INCLUDED

