/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 5, 2025 $
//
//

// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define_python.h"
#include "custom_utilities/cell.h"
#include "custom_utilities/cell_container.h"
#include "custom_utilities/nurbs/bcell.h"
#include "custom_utilities/nurbs/bcell_manager.h"
#include "custom_python/add_cell_container_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

boost::python::list CellContainerToList(const CellContainer& c)
{
    boost::python::list l;
    for (auto it = c.cbegin(); it != c.cend(); ++it)
        l.append(boost::python::ptr(*it));  // wrap Cell* into Python object
    return l;
}

boost::python::list Cell_GetSupportedAnchors(Cell& c)
{
    boost::python::list result;

    const auto& anchors = c.GetSupportedAnchors();
    for (auto i : anchors)
        result.append(i);

    return result;
}

void IsogeometricApplication_AddCellContainerToPython()
{

    class_<Cell, Cell::Pointer, boost::noncopyable>
    ("Cell", no_init)
    .add_property("Id", &Cell::Id)
    .add_property("Level", &Cell::Level)
    .def("GetSupportedAnchors", &Cell_GetSupportedAnchors)
    .def(self_ns::str(self))
    ;

    class_<BCell, BCell::Pointer, bases<Cell>, boost::noncopyable>
    ("BCell", no_init)
    .def(self_ns::str(self))
    ;

    class_<CellContainer, CellContainer::Pointer, boost::noncopyable>
    ("CellContainer", init<>())
    .def("__len__", &CellContainer::size)
    .def("ToList", &CellContainerToList)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos
