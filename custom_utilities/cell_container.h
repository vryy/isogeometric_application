//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED

// System includes
#include <set>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/cell.h"

namespace Kratos
{

/**
 * CellContainer provides general container to stores Cell's. It provides the iterator interface to iterator through cell.
 * However, the CellContainer cannot modify itself from outside. The reason for that is the CellContainer is the base for
 * all kind of cell container with any kind of cell, to avoid lock in a specific type of cell via template. It shall only
 * be modified from the subclass.
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

    /// Iterators
    const_iterator cbegin() const {return mpCells.begin();}
    const_iterator cend() const {return mpCells.end();}

    /// Get the number of cells of this container
    std::size_t size() const {return mpCells.size();}

    /// Overload comparison operator
    bool operator==(const CellContainer& rOther) const
    {
        // check size
        if (this->size() != rOther.size())
        {
            return false;
        }

        // check specific ids
        const_iterator it1, it2;
        for (it1 = this->cbegin(), it2 = rOther.cbegin(); it1 != this->cend(); ++it1, ++it2)
        {
            if ((*it1)->Id() != (*it2)->Id())
                return false;
        }

        return true;
    }

    /// Overload comparison operator
    bool operator!=(const CellContainer& rOther) const
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

protected:

    iterator begin() {return mpCells.begin();}
    iterator end() {return mpCells.end();}

    /// Insert a cell to the container.
    void insert(cell_t p_cell)
    {
        mpCells.insert(p_cell);
    }

    /// Remove a cell from the container
    void erase(cell_t p_cell)
    {
        mpCells.erase(p_cell);
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

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CELL_CONTAINER_H_INCLUDED
