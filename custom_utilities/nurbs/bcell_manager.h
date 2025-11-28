//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_MANAGER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_MANAGER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/cell_container.h"

// #define USE_BRUTE_FORCE_TO_SEARCH_FOR_CELLS
#define USE_R_TREE_TO_SEARCH_FOR_CELLS

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
#include "custom_external_libraries/RTree.h"
#endif

namespace Kratos
{

bool BCellManager_RtreeSearchCallback(std::size_t id, void* arg);

struct BCellManager_Helper
{
    template<typename knot_t>
    static inline double GetKnotValue(const knot_t& p_knot)
    {
        return 0.0;
    }
};

template<>
inline double BCellManager_Helper::GetKnotValue<double>(const double& r_knot)
{
    return r_knot;
}

template<>
inline double BCellManager_Helper::GetKnotValue<Knot<double> >(const Knot<double>& r_knot)
{
    return r_knot.Value();
}

template<>
inline double BCellManager_Helper::GetKnotValue<typename Knot<double>::Pointer>(const typename Knot<double>::Pointer& p_knot)
{
    return p_knot->Value();
}

/**
 * Abstract cell manager for management of collection of b-cells. It provides facility to search for cells, or obtain cells in the consistent manner.
 * TCellType must be sub-class of BCell
 */
template<class TCellType>
class KRATOS_API(ISOGEOMETRIC_APPLICATION) BaseBCellManager : public CellContainer
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BaseBCellManager);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const BaseBCellManager> ConstPointer;
#endif

    /// Type definitions
    typedef CellContainer BaseType;
    typedef TCellType CellType;
    typedef typename CellType::knot_t knot_t;
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::IndexType IndexType;
    struct cell_compare
    {
        bool operator() (const cell_t& lhs, const cell_t& rhs) const {return lhs->Id() < rhs->Id();}
    };
    typedef std::set<cell_t, cell_compare> cell_container_t;
    typedef std::map<std::size_t, cell_t> map_t;
    typedef typename cell_container_t::iterator iterator;
    typedef typename cell_container_t::const_iterator const_iterator;

    /// Default constructor
    BaseBCellManager() : mTol(1.0e-10), mLastId(0)
    {}

    /// Destructor
    ~BaseBCellManager() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        this->PrintInfo(std::cout); std::cout << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Set the tolerance for the internal searching algorithm
    void SetTolerance(double Tol) {mTol = Tol;}

    /// Get the tolerance for the internal searching algorithm
    double GetTolerance() const {return mTol;}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
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
    iterator end() {return mpCells.end();}
    const_iterator cbegin() const {return mpCells.begin();}
    const_iterator cend() const {return mpCells.end();}

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
            this->CreateCellsMap();
        }

        // return the bf if its Id exist in the list
        typename map_t::iterator it = mCellsMap.find(Id);
        if (it != mCellsMap.end())
        {
            return it->second;
        }
        else
            KRATOS_ERROR << "Access index " << Id << " is not found";
        }

    /// Overload operator[]
    cell_t operator[](std::size_t Id)
    {
        return get(Id);
    }

    /// Overload comparison operator
    bool operator==(const BaseBCellManager<TCellType>& rOther)
    {
        if (this->size() != rOther.size())
        {
            return false;
        }

        iterator it_this = this->begin();
        iterator it_other = rOther.begin();

        while ((it_this != this->end()) && (it_other != rOther.end()))
        {
            if (!(*it_this)->IsSame(*it_other, 1.0e-10))
            {
                return false;
            }
            ++it_this;
            ++it_other;
        }

        return true;
    }

    /// Overload comparison operator
    bool operator!=(const BaseBCellManager<TCellType>& rOther)
    {
        return !(*this == rOther);
    }

    /// Search the cells covered in another cell. In return, std::vector<cell_t> are all the cells covered by p_cell.
    virtual std::vector<cell_t> GetCells(cell_t p_cell)
    {
        KRATOS_ERROR << "Calling the base class function";
    }

    /// Collapse the overlapping cells
    void CollapseCells()
    {
        bool hit;
        do
        {
            iterator it_cell = this->begin();
            hit = this->CollapseCells(it_cell, this->end());
            if (hit)
            {
                (*it_cell)->ClearTrace();
                this->erase(*it_cell);
            }
        }
        while (hit);
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
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BaseBCellManager";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:

    cell_container_t mpCells;
    mutable map_t mCellsMap; // map from cell id to the basis function. It's mainly used to search for the cell quickly. But it needs to be re-initialized whenever new cell is added to the set
    bool cell_map_is_created;
    IndexType mLastId;

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
            {
                ++it_cell;
            }
        }
        while (it_cell != cell_end);

        return false;
    }
};

/**
 * Abstract BCell Manager
 */
template<int TDim, class TCellType>
class KRATOS_API(ISOGEOMETRIC_APPLICATION) BCellManager : public BaseBCellManager<TCellType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BCellManager);
};

/**
 * BCell Manager in 1D
 */
template<class TCellType>
class BCellManager<1, TCellType> : public BaseBCellManager<TCellType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BCellManager);

    /// Type definitions
    typedef BaseBCellManager<TCellType> BaseType;
    typedef typename BaseType::BaseType SuperType;
    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::knot_t knot_t;
    typedef typename BaseType::iterator iterator;

    /// Default constructor
    BCellManager() : BaseType()
    {}

    /// Destructor
    ~BCellManager() override
    {}

    /// Helper function to create new instance of cell manager
    static typename BaseType::Pointer Create() {return typename BaseType::Pointer(new BCellManager<1, CellType>());}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
    {
        assert(pKnots.size() == 2);

        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future.
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if ( (*it)->XiMin() == pKnots[0] // left
                    && (*it)->XiMax() == pKnots[1] ) // right
            {
                return *it;
            }
        }

        // otherwise create new cell
        cell_t p_cell = cell_t(new TCellType(++BaseType::mLastId, pKnots[0], pKnots[1]));
        BaseType::mpCells.insert(p_cell);
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {BCellManager_Helper::GetKnotValue(pKnots[0])};
        double cmax[] = {BCellManager_Helper::GetKnotValue(pKnots[1])};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return p_cell;
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it == p_cell)
            {
                return it;
            }

        // otherwise insert new cell
        iterator it = BaseType::mpCells.insert(p_cell).first;
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {p_cell->XiMinValue()};
        double cmax[] = {p_cell->XiMaxValue()};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return it;
    }

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if (*it == p_cell)
            {
                BaseType::mpCells.erase(it);
                SuperType::erase(&(**it));

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
                // update the r-tree
                double cmin[] = {p_cell->XiMinValue()};
                double cmax[] = {p_cell->XiMaxValue()};
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
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it != p_cell)
                if ((*it)->template IsCovered<1>(p_cell))
                {
                    p_cells.push_back(*it);
                }
#endif

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // determine the overlapping cells; for now, this only works in 3D
        std::vector<std::size_t> OverlappingCells;
        double cmin[] = {p_cell->XiMinValue()};
        double cmax[] = {p_cell->XiMaxValue()};
        int nhits = rtree_cells.Search(cmin, cmax, BCellManager_RtreeSearchCallback, (void*)(&OverlappingCells));
//        printf("Search resulted in %d hits\n", nhits);

        // check within overlapping cells the one covered in p_cell
        for (std::size_t i = 0; i < OverlappingCells.size(); ++i)
        {
            cell_t pthis_cell = this->get(OverlappingCells[i]);
            if (pthis_cell != p_cell)
                if (pthis_cell->template IsCovered<1>(p_cell))
                {
                    p_cells.push_back(pthis_cell);
                }
        }
#endif

        return p_cells;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BCellManager1D";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
    RTree<std::size_t, double, 1, double> rtree_cells;
#endif
};

/**
 * BCell Manager in 2D
 */
template<class TCellType>
class BCellManager<2, TCellType> : public BaseBCellManager<TCellType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BCellManager);

    /// Type definitions
    typedef BaseBCellManager<TCellType> BaseType;
    typedef typename BaseType::BaseType SuperType;
    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::knot_t knot_t;
    typedef typename BaseType::iterator iterator;

    /// Default constructor
    BCellManager() : BaseType()
    {}

    /// Destructor
    ~BCellManager() override
    {}

    /// Helper function to create new instance of cell manager
    static typename BaseType::Pointer Create() {return typename BaseType::Pointer(new BCellManager<2, CellType>());}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
    {
        assert(pKnots.size() == 4);

        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future.
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if ( (*it)->XiMin()  == pKnots[0] // left
                    && (*it)->XiMax()  == pKnots[1] // right
                    && (*it)->EtaMin() == pKnots[2] // down
                    && (*it)->EtaMax() == pKnots[3] ) // up
            {
                return *it;
            }
        }

        // otherwise create new cell
        cell_t p_cell = cell_t(new TCellType(++BaseType::mLastId, pKnots[0], pKnots[1], pKnots[2], pKnots[3]));
        BaseType::mpCells.insert(p_cell);
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {BCellManager_Helper::GetKnotValue(pKnots[0]), BCellManager_Helper::GetKnotValue(pKnots[2])};
        double cmax[] = {BCellManager_Helper::GetKnotValue(pKnots[1]), BCellManager_Helper::GetKnotValue(pKnots[3])};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return p_cell;
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it == p_cell)
            {
                return it;
            }

        // otherwise insert new cell
        iterator it = BaseType::mpCells.insert(p_cell).first;
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue()};
        double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue()};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return it;
    }

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if (*it == p_cell)
            {
                BaseType::mpCells.erase(it);
                SuperType::erase(&(**it));

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
                // update the r-tree
                double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue(), p_cell->ZetaMinValue()};
                double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue(), p_cell->ZetaMaxValue()};
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
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it != p_cell)
                if ((*it)->template IsCovered<2>(p_cell))
                {
                    p_cells.push_back(*it);
                }
#endif

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // determine the overlapping cells; for now, this only works in 3D
        std::vector<std::size_t> OverlappingCells;
        double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue()};
        double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue()};
        int nhits = rtree_cells.Search(cmin, cmax, BCellManager_RtreeSearchCallback, (void*)(&OverlappingCells));
//        printf("Search resulted in %d hits\n", nhits);

        // check within overlapping cells the one covered in p_cell
        for (std::size_t i = 0; i < OverlappingCells.size(); ++i)
        {
            cell_t pthis_cell = this->get(OverlappingCells[i]);
            if (pthis_cell != p_cell)
                if (pthis_cell->template IsCovered<2>(p_cell))
                {
                    p_cells.push_back(pthis_cell);
                }
        }
#endif

        return p_cells;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BCellManager2D";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
    RTree<std::size_t, double, 2, double> rtree_cells;
#endif
};

/**
 * BCell Manager in 3D
 */
template<class TCellType>
class BCellManager<3, TCellType> : public BaseBCellManager<TCellType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BCellManager);

    /// Type definitions
    typedef BaseBCellManager<TCellType> BaseType;
    typedef typename BaseType::BaseType SuperType;
    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::knot_t knot_t;
    typedef typename BaseType::iterator iterator;

    /// Default constructor
    BCellManager() : BaseType()
    {}

    /// Destructor
    ~BCellManager() override
    {}

    /// Helper function to create new instance of cell manager
    static typename BaseType::Pointer Create() {return typename BaseType::Pointer(new BCellManager<3, CellType>());}

    /// Check if the cell exists in the list; otherwise create new cell and return
    virtual cell_t CreateCell(const std::vector<knot_t>& pKnots)
    {
        assert(pKnots.size() == 6);

        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future.
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if ( (*it)->XiMin()  == pKnots[0] // left
                    && (*it)->XiMax()  == pKnots[1] // right
                    && (*it)->EtaMin() == pKnots[2] // down
                    && (*it)->EtaMax() == pKnots[3] // up
                    && (*it)->ZetaMin() == pKnots[4] // below
                    && (*it)->ZetaMax() == pKnots[5] ) // above
            {
                return *it;
            }
        }

        // otherwise create new cell
        cell_t p_cell = cell_t(new TCellType(++BaseType::mLastId, pKnots[0], pKnots[1], pKnots[2], pKnots[3], pKnots[4], pKnots[5]));
        BaseType::mpCells.insert(p_cell);
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {BCellManager_Helper::GetKnotValue(pKnots[0]), BCellManager_Helper::GetKnotValue(pKnots[2]), BCellManager_Helper::GetKnotValue(pKnots[4])};
        double cmax[] = {BCellManager_Helper::GetKnotValue(pKnots[1]), BCellManager_Helper::GetKnotValue(pKnots[3]), BCellManager_Helper::GetKnotValue(pKnots[5])};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return p_cell;
    }

    /// Insert a cell to the container. If the cell is existed in the container, the iterator of the existed one will be returned.
    virtual iterator insert(cell_t p_cell)
    {
        // search in the list of cell if any cell has the same knot span
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it == p_cell)
            {
                return it;
            }

        // otherwise insert new cell
        iterator it = BaseType::mpCells.insert(p_cell).first;
        SuperType::insert(&(*p_cell));
        BaseType::cell_map_is_created = false;

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // update the r-tree
        double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue(), p_cell->ZetaMinValue()};
        double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue(), p_cell->ZetaMaxValue()};
        rtree_cells.Insert(cmin, cmax, p_cell->Id());
#endif

        return it;
    }

    /// Remove a cell by its Id from the set
    virtual void erase(cell_t p_cell)
    {
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
        {
            if (*it == p_cell)
            {
                BaseType::mpCells.erase(it);
                SuperType::erase(&(**it));

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
                // update the r-tree
                double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue(), p_cell->ZetaMinValue()};
                double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue(), p_cell->ZetaMaxValue()};
                rtree_cells.Remove(cmin, cmax, p_cell->Id());
#endif

                break;
            }
        }
    }

    /// Search the cells coverred in another cell. In return p_cell covers all the cells of std::vector<cell_t>
    virtual std::vector<cell_t> GetCells(cell_t p_cell)
    {
        std::vector<cell_t> p_cells;

#ifdef USE_BRUTE_FORCE_TO_SEARCH_FOR_CELLS
        // Currently I use the brute-force approach. I know it is not efficient. I will improve it in the future. TODO
        for (iterator it = BaseType::mpCells.begin(); it != BaseType::mpCells.end(); ++it)
            if (*it != p_cell)
                if ((*it)->template IsCovered<3>(p_cell))
                {
                    p_cells.push_back(*it);
                }
#endif

#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
        // determine the overlapping cells; for now, this's only working in 3D
        std::vector<std::size_t> OverlappingCells;
        double cmin[] = {p_cell->XiMinValue(), p_cell->EtaMinValue(), p_cell->ZetaMinValue()};
        double cmax[] = {p_cell->XiMaxValue(), p_cell->EtaMaxValue(), p_cell->ZetaMaxValue()};
        int nhits = rtree_cells.Search(cmin, cmax, BCellManager_RtreeSearchCallback, (void*)(&OverlappingCells));
//        printf("Search resulted in %d hits\n", nhits);

        // check within overlapping cells the one covered in p_cell
        for (std::size_t i = 0; i < OverlappingCells.size(); ++i)
        {
            cell_t pthis_cell = this->get(OverlappingCells[i]);
            if (pthis_cell != p_cell)
                if (pthis_cell->template IsCovered<3>(p_cell))
                {
                    p_cells.push_back(pthis_cell);
                }
        }
#endif

        return p_cells;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BCellManager3D";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
#ifdef USE_R_TREE_TO_SEARCH_FOR_CELLS
    RTree<std::size_t, double, 3, double> rtree_cells;
#endif
};

/// output stream function
template<int TDim, class TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const BCellManager<TDim, TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BCELL_MANAGER_H_INCLUDED
