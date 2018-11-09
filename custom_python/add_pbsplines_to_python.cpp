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
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/pbsplines_basis_function.h"
#include "custom_utilities/pbsplines_fespace.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_point_based_control_grid_to_python.h"
#include "custom_python/add_import_export_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim>
std::size_t PBSplinesBasisFunction_GetId(PBSplinesBasisFunction<TDim>& rDummy)
{
    return rDummy.Id();
}

template<int TDim>
void PBSplinesBasisFunction_SetId(PBSplinesBasisFunction<TDim>& rDummy, std::size_t Id)
{
    // DO NOTHING
}

template<int TDim>
std::size_t PBSplinesBasisFunction_GetEquationId(PBSplinesBasisFunction<TDim>& rDummy)
{
    return rDummy.EquationId();
}

template<int TDim>
void PBSplinesBasisFunction_SetEquationId(PBSplinesBasisFunction<TDim>& rDummy, std::size_t EquationId)
{
    rDummy.SetEquationId(EquationId);
}

template<int TDim>
boost::python::list PBSplinesFESpace_ExtractBoundaryBfsByFlag(PBSplinesFESpace<TDim>& rDummy, std::size_t boundary_id)
{
    typedef typename PBSplinesFESpace<TDim>::bf_t bf_t;

    std::vector<bf_t> bf_list = rDummy.ExtractBoundaryBfsByFlag(boundary_id);

    boost::python::list Output;
    for (std::size_t i = 0; i < bf_list.size(); ++i)
        Output.append(bf_list[i]);

    return Output;
}

template<int TDim>
typename PBSplinesBasisFunction<TDim>::Pointer PBSplinesFESpace_GetItem(PBSplinesFESpace<TDim>& rDummy, std::size_t i)
{
    return rDummy[i];
}

////////////////////////////////////////

template<int TDim>
void PBSplinesPatchUtility_ListBoundaryBfs(PBSplinesPatchUtility& rDummy,
    typename PBSplinesFESpace<TDim>::Pointer pFESpace, BoundarySide side)
{
    rDummy.ListBoundaryBfs<TDim>(std::cout, pFESpace, side);
}

template<typename TDataType, class TFESpaceType>
typename ControlGrid<TDataType>::Pointer ControlGridUtility_CreatePointBasedControlGrid(
        ControlGridUtility& rDummy,
        const Variable<TDataType>& rVariable, typename TFESpaceType::Pointer pFESpace)
{
    return rDummy.CreatePointBasedControlGrid<TDataType, TFESpaceType>(rVariable, pFESpace);
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddPBSplinesSpaceToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "PBSplinesBasisFunction" << TDim << "D";
    class_<PBSplinesBasisFunction<TDim>, typename PBSplinesBasisFunction<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<const std::size_t&, const std::size_t&>())
    .add_property("Id", PBSplinesBasisFunction_GetId<TDim>, PBSplinesBasisFunction_SetId<TDim>)
    .add_property("EquationId", PBSplinesBasisFunction_GetEquationId<TDim>, PBSplinesBasisFunction_SetEquationId<TDim>)
    .def("Weight", &PBSplinesBasisFunction<TDim>::Weight)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "PBSplinesFESpace" << TDim << "D";
    typename FESpace<TDim-1>::Pointer(PBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace1)(const BoundarySide& side) const = &PBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    // typename FESpace<TDim-1>::Pointer(PBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace2)(const BoundarySide& side, const BoundaryRotation& rotation) const = &PBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    class_<PBSplinesFESpace<TDim>, typename PBSplinesFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__getitem__", &PBSplinesFESpace_GetItem<TDim>)
    .def("GetBoundaryBfs", &PBSplinesFESpace_ExtractBoundaryBfsByFlag<TDim>) // deprecated
    .def("ExtractBoundaryBfsByFlag", &PBSplinesFESpace_ExtractBoundaryBfsByFlag<TDim>)
    .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace1)
    // .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace2)
    .def("UpdateCells", &PBSplinesFESpace<TDim>::UpdateCells)
    .def(self_ns::str(self))
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, PBSplinesFESpace<TDim> >::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, PBSplinesFESpace<TDim> >::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, PBSplinesFESpace<TDim> >::Execute();

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
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBSplinesFESpace<3> >)
    ;

}

}  // namespace Python.

} // Namespace Kratos

