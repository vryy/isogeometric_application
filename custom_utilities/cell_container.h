//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED

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
#include "custom_utilities/cell.h"


namespace Kratos
{

/**
 * CellContainer provides general container to stores Cell's. It provides the iterator interface and the methods to insert or remove cells.
 */
class CellContainer
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(CellContainer);

    /// Type definitions
    typedef Cell CellType;
    typedef CellType* cell_t;
    struct cell_compare
    {
        bool operator() (const cell_t& lhs, const cell_t& rhs) const {return lhs->Id() < rhs->Id();}
    };
    typedef std::set<cell_t, cell_compare> cell_container_t;
    typedef typename cell_container_t::iterator iterator;
    typedef typename cell_container_t::const_iterator const_iterator;

    /// Default constructor
    CellContainer()
    {}

    /// Destructor
    virtual ~CellContainer()
    {
        #ifdef ISOGEOMETRIC_DEBUG_DESTROY
        this->PrintInfo(std::cout); std::cout << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Insert a cell to the container.
    void insert(cell_t p_cell)
    {
        mpCells.insert(p_cell);
    }

    /// Iterators
    iterator begin() {return mpCells.begin();}
    const_iterator begin() const {return mpCells.begin();}
    iterator end() {return mpCells.end();}
    const_iterator end() const {return mpCells.end();}

    /// Get the number of cells of this container
    std::size_t size() const {return mpCells.size();}

    /// Remove a cell from the container
    void erase(cell_t p_cell)
    {
        mpCells.erase(p_cell);
    }

    /// Overload comparison operator
    bool operator==(const CellContainer& rOther)
    {
        if (this->size() != rOther.size())
            return false;

        // TODO check more

        return true;
    }

    /// Overload comparison operator
    bool operator!=(const CellContainer& rOther)
    {
        return !(*this == rOther);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CellContainer";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    cell_container_t mpCells;

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const CellContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED

