/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 24, 2017 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes
#include <string>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "python/pointer_vector_set_python_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/pbsplines_basis_function.h"
#include "custom_utilities/pbsplines_fespace.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_point_based_control_grid_to_python.h"
#include "custom_python/add_import_export_to_python.h"
#include "custom_python/add_control_grids_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<int TDim>
void IsogeometricApplication_AddPBSplinesSpaceToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "PBSplinesBasisFunction" << TDim << "D";
    typedef PBSplinesBasisFunction<TDim, Cell> PBSplinesBasisFunctionType;
    class_<PBSplinesBasisFunctionType, typename PBSplinesBasisFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<const std::size_t&>())
    .add_property("Id", Isogeometric_GetId<PBSplinesBasisFunctionType>, Isogeometric_SetId<PBSplinesBasisFunctionType>)
    .add_property("EquationId", Isogeometric_GetEquationId<PBSplinesBasisFunctionType>, Isogeometric_SetEquationId<PBSplinesBasisFunctionType>)
    .def("Weight", &PBSplinesBasisFunctionType::Weight)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "PBSplinesFESpace" << TDim << "D";
    typedef FESpace<TDim> FESpaceType;
    typedef PBSplinesFESpace<TDim, PBSplinesBasisFunctionType> PBSplinesFESpaceType;
    class_<PBSplinesFESpaceType, typename PBSplinesFESpaceType::Pointer, bases<FESpaceType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__getitem__", &FESpace_GetItem<PBSplinesFESpaceType>)
    .def("UpdateCells", &PBSplinesFESpaceType::UpdateCells)
    .def(self_ns::str(self))
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, PBSplinesFESpaceType>::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, PBSplinesFESpaceType>::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, PBSplinesFESpaceType>::Execute();

}

////////////////////////////////////////

void IsogeometricApplication_AddPBSplinesToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////Point-based BSplines//////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddPBSplinesSpaceToPython<1>();
    IsogeometricApplication_AddPBSplinesSpaceToPython<2>();
    IsogeometricApplication_AddPBSplinesSpaceToPython<3>();

    class_<ControlGridUtility, ControlGridUtility::Pointer, boost::noncopyable>
    ("PointBasedControlGridUtility", init<>())
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<1, PBSplinesBasisFunction<1, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<2, PBSplinesBasisFunction<2, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<3, PBSplinesBasisFunction<3, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<1, PBSplinesBasisFunction<1, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<2, PBSplinesBasisFunction<2, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<3, PBSplinesBasisFunction<3, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<1, PBSplinesBasisFunction<1, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<2, PBSplinesBasisFunction<2, Cell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<3, PBSplinesBasisFunction<3, Cell> > >)
    ;

}

}  // namespace Python.

} // Namespace Kratos

