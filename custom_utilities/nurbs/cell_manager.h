//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_H_INCLUDED

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

bool CellManager_RtreeSearchCallback(std::size_t id, void* arg);

/**
A cell manager is a collection of cells. It provides facility to search for cells, or obtain cells in the consistent manner.
 */
template<class TCellType>
class CellManager
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(CellManager);

    /// Type definitions
    typedef Knot<double> KnotType;
    typedef KnotType::Pointer knot_t;

    typedef TCellType CellType;
    typedef typename CellType::Pointer cell_t;
    struct cell_compare
    {
        bool operator() (const cell_t& lhs, const cell_t& rhs) const {return lhs->Id() < rhs->Id();}
    };
    typedef std::set<cell_t, cell_compare> cell_container_t;
    typedef std::map<std::size_t, cell_t> map_t;
    typedef typename cell_container_t::iterator iterator;
    typedef typename cell_container_t::const_iterator const_iterator;

    /// Default constructor
    CellManager() : mTol(1.0e-10), mLastId(0)
    {}

    /// Destructor
    virtual ~CellManager()
    {}

    /// Set the tolerance for the internal searching algorithm
    void SetTolerance(const double& Tol) {mTol = Tol;}

    /// Get the tolerance for the internal searching algorithm
    const double& GetTolerance() const {return mTol;}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the virtual function", __FUNCTION__)
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the virtual function", __FUNCTION__)
    }

    /// Iterators
    iterator begin() {return mpCells.begin();}
    const_iterator begin() const {return mpCells.begin();}
    iterator end() {return mpCells.end();}
    const_iterator end() const {return mpCells.end();}

    /// Get the number of cells of this manager
    std::size_t size() const {return mpCells.size();}

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the virtual function", __FUNCTION__)
    }

    /// Get a cell based on its Id
    cell_t get(const std::size_t& Id)
    {
        // create the index map if it's not created yet
        if(!cell_map_is_created)
            this->CreateCellsMap();

        // return the bf if its Id exist in the list
        typename map_t::iterator it = mCellsMap.find(Id);
        if(it != mCellsMap.end())
            return it->second;
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Access index is not found:", Id)
    }

    /// Overload operator[]
    cell_t operator[](const std::size_t& Id)
    {
        return get(Id);
    }

    /// Overload comparison operator
    bool operator==(const CellManager<TCellType>& rOther)
    {
        if (this->size() != rOther.size())
            return false;

        iterator it_this = this->begin();
        iterator it_other = rOther.begin();

        while ((it_this != this->end()) && (it_other != rOther.end()))
        {
            if (!(*it_this)->IsSame(*it_other, 1.0e-10))
                return false;
            ++it_this;
            ++it_other;
        }

        return true;
    }

    /// Overload comparison operator
    bool operator!=(const CellManager<TCellType>& rOther)
    {
        return !(*this == rOther);
    }

    /// Search the cells covered in another cell. In return, std::vector<cell_t> are all the cells covered by p_cell.
    virtual std::vector<cell_t> GetCells(cell_t p_cell)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the virtual function", __FUNCTION__)
    }

    /// Collapse the overlapping cells
    void CollapseCells()
    {
        bool hit;
        do
        {
            iterator it_cell = this->begin();
            hit = this->CollapseCells(it_cell, this->end());
            if (hit) this->erase(*it_cell);
        } while (hit);
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

    cell_container_t mpCells;
    mutable map_t mCellsMap; // map from cell id to the basis function. It's mainly used to search for the cell quickly. But it needs to be re-initialized whenever new cell is added to the set
    bool cell_map_is_created;
    std::size_t mLastId;

private:

    double mTol;

    void CreateCellsMap()
    {
        mCellsMap.clear();
        for(iterator it = mpCells.begin(); it != mpCells.end(); ++it)
            mCellsMap[(*it)->Id()] = *it;
        cell_map_is_created = true;
    }

    /// Collapse the first found overlapping cells
    bool CollapseCells(iterator& it_cell, const_iterator cell_end)
    {
        do
        {
            std::vector<cell_t> inner_cells = this->GetCells(*it_cell);

            if (inner_cells.size() != 0)
            {
                for (std::size_t i = 0; i < inner_cells.size(); ++i)
                {
                    inner_cells[i]->Absorb(*it_cell);
                }
                return true;
            }
            else
                ++it_cell;
        } while (it_cell != cell_end);

        return false;
    }
};

/// output stream function
template<class TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const CellManager<TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_H_INCLUDED

