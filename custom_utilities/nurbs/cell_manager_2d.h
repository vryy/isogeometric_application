//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_2D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_2D_H_INCLUDED

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
#include "custom_utilities/nurbs/cell_manager.h"

// #define USE_BRUTE_FORCE_TO_SEARCH_FOR_CELLS
#define USE_R_TREE_TO_SEARCH_FOR_CELLS

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
#include "custom_external_libraries/RTree.h"
#endif

namespace Kratos
{

/**
 * Cell Manager in 2D
 */
template<class TCellType>
class CellManager2D : public CellManager<TCellType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(CellManager2D);

    /// Type definitions
    typedef CellManager<TCellType> BaseType;
    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::knot_t knot_t;
    typedef typename BaseType::iterator iterator;

    /// Default constructor
    CellManager2D() : BaseType()
    {}

    /// Destructor
    virtual ~CellManager2D()
    {}

    /// Helper function to create new instance of cell manager
    static typename BaseType::Pointer Create() {return typename BaseType::Pointer(new CellManager2D<CellType>());}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
    {
        assert(pKnots.size() == 4);

        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future.
        for(iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if( (*it)->Left()  == pKnots[0] // left
             && (*it)->Right() == pKnots[1] // right
             && (*it)->Down()  == pKnots[2] // down
             && (*it)->Up()    == pKnots[3] ) // up
                return *it;
        }

        // otherwise create new cell
        cell_t p_cell = cell_t(new TCellType(++BaseType::mLastId, pKnots[0], pKnots[1], pKnots[2], pKnots[3]));
        BaseType::mpCells.insert(p_cell);
        BaseType::cell_map_is_created = false;

        #ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {pKnots[0]->Value(), pKnots[2]->Value()};
        double cmax[] = {pKnots[1]->Value(), pKnots[3]->Value()};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
        #endif

        return p_cell;
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for(iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if(*it == p_cell)
                return it;

        // otherwise insert new cell
        iterator it = BaseType::mpCells.insert(p_cell).first;
        BaseType::cell_map_is_created = false;

        #ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {p_cell->LeftValue(), p_cell->DownValue()};
        double cmax[] = {p_cell->RightValue(), p_cell->UpValue()};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
        #endif

        return it;
    }

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        for(iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if(*it == p_cell)
            {
                BaseType::mpCells.erase(it);

                #ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
                // update the r-tree
                double cmin[] = {p_cell->LeftValue(), p_cell->DownValue(), p_cell->BelowValue()};
                double cmax[] = {p_cell->RightValue(), p_cell->UpValue(), p_cell->AboveValue()};
                rtree_cells.Remove(cmin, cmax, p_cell->Id());
                #endif

                break;
            }
        }
    }

    /// Search the cells covered in another cell. In return p_cell covers all the cells of std::vector<cell_t>
    virtual std::vector<cell_t> GetCells(cell_t p_cell)
    {
        std::vector<cell_t> p_cells;

        #ifdef USE_BRUTE_FORCE_TO_SEARCH_FOR_CELLS
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for(iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if(*it != p_cell)
                if((*it)->template IsCovered<2>(p_cell))
                    p_cells.push_back(*it);
        #endif

        #ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // determine the overlapping cells; for now, this only works in 3D
        std::vector<std::size_t> OverlappingCells;
        double cmin[] = {p_cell->LeftValue(), p_cell->DownValue()};
        double cmax[] = {p_cell->RightValue(), p_cell->UpValue()};
        int nhits = rtree_cells.Search(cmin, cmax, CellManager_RtreeSearchCallback, (void*)(&OverlappingCells));
//        printf("Search resulted in %d hits\n", nhits);

        // check within overlapping cells the one covered in p_cell
        for(std::size_t i = 0; i < OverlappingCells.size(); ++i)
        {
            cell_t pthis_cell = this->get(OverlappingCells[i]);
            if(pthis_cell != p_cell)
                if(pthis_cell->template IsCovered<2>(p_cell))
                    p_cells.push_back(pthis_cell);
        }
        #endif

        return p_cells;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    #ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
    RTree<std::size_t, double, 2, double> rtree_cells;
    #endif
};

/// output stream function
template<class TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const CellManager2D<TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#undef USE_BRUTE_FORCE_TO_SEARCH_FOR_CELLS
#undef USE_R_TREE_TO_SEARCH_FOR_CELLS

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_MANAGER_2D_H_INCLUDED

