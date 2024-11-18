//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 16 Oct 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_TREE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_TREE_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"

namespace Kratos
{

bool CellTree_RtreeSearchCallback(std::size_t id, void* arg);

/**
  TODO
 */
template<class TCellType>
class CellTree
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(CellTree);

    /// Type definitions
    typedef Knot<double> KnotType;
    typedef KnotType::Pointer knot_t;

    typedef typename TCellType::Pointer cell_t;
    struct cell_compare
    {
        bool operator() (const cell_t& lhs, const cell_t& rhs) const {return lhs->Id() < rhs->Id();}
    };
    typedef std::set<cell_t, cell_compare> cell_container_t;
    typedef std::map<std::size_t, cell_t> map_t;
    typedef typename cell_container_t::iterator iterator;

    /// Default constructor
    CellTree() : mTol(1.0e-10), mLastId(0)
    {}

    /// Destructor
    virtual ~CellTree()
    {}

    /// Set the tolerance for the internal searching algorithm
    void SetTolerance(double Tol) {mTol = Tol;}

    /// Get the tolerance for the internal searching algorithm
    double GetTolerance() const {return mTol;}

    /// Get the current level of the cell tree top node
    int Level() const {return mLevel;}

    /// Check if the cell exists in the list; ortherwise create new cell and return
    virtual cell_t CreateCell(unsigned int Level, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Check if the cell exists in the list; ortherwise create new cell and return
    virtual cell_t CreateCell(unsigned int Level, knot_t pLeft, knot_t pRight, knot_t pDown, knot_t pUp, knot_t pBelow, knot_t pAbove)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Iterators
    iterator begin() {return mpCells.begin();}
    iterator begin() const {return mpCells.begin();}
    iterator end() {return mpCells.end();}
    iterator end() const {return mpCells.end();}

    /// Get the number of cells of this manager
    std::size_t size() const {return mpCells.size();}

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Get a cell based on its Id
    cell_t get(std::size_t Id)
    {
        // create the index map if it's not created yet
        if (!cell_map_is_created)
        {
            CreateCellsMap();
        }

        // return the bf if its Id exist in the list
        typename map_t::iterator it = mCellsMap.find(Id);
        if (it != mCellsMap.end())
        {
            return it->second;
        }
        else
            KRATOS_ERROR << "Index " << Id << " is not found";
    }

    /// Overload operator[]
    cell_t operator[](std::size_t Id)
    {
        return get(Id);
    }

    /// Search the cells coverred in another cell. In return p_cell covers all the cells of std::vector<cell_t>
    virtual std::vector<cell_t> GetCells(cell_t p_cell)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Reset all the Id of all the basis functions. Remarks: use it with care, you have to be responsible to the old indexing data of the basis functions before calling this function
    /// Disable this function for temporary
//    std::size_t ReIndexing()
//    {
//        mLastId = 0;
//        for(iterator it = mpCells.begin(); it != mpCells.end(); ++it)
//            (*it)->SetId(++mLastId);
//        return mLastId;
//    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

protected:
    int mLevel;
    cell_container_t mpCells;
    mutable map_t mCellsMap; // map from cell id to the basis function. It's mainly used to search for the cell quickly. But it needs to be re-initialized whenever new cells are added to the set
    bool cell_map_is_created;
    std::size_t mLastId;
    std::vector<CellTree::Pointer> mpSubTrees;

private:
    double mTol;

    void CreateCellsMap()
    {
        mCellsMap.clear();
        for (iterator it = mpCells.begin(); it != mpCells.end(); ++it)
        {
            mCellsMap[(*it)->Id()] = *it;
        }
        cell_map_is_created = true;
    }
};

/// output stream function
template<class TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const CellTree<TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_TREE_H_INCLUDED
