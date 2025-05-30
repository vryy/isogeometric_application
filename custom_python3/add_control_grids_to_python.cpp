/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/control_grid_library.h"
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/tsplines/tcell.h"
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/nurbs/pbbsplines_fespace.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_python3/add_control_grids_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateLinearControlPointGrid(
    ControlGridLibrary& rDummy,
    double start_x, double start_y, double start_z,
    std::size_t n_points_u,
    double end_x, double end_y, double end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateStructuredControlPointGrid<1>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid1(
    ControlGridLibrary& rDummy,
    double start_x, double start_y,
    std::size_t n_points_u, std::size_t n_points_v,
    double end_x, double end_y)
{
    std::vector<double> start(2);
    start[0] = start_x;
    start[1] = start_y;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> end(2);
    end[0] = end_x;
    end[1] = end_y;

    return rDummy.CreateStructuredControlPointGrid<2>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid2(
    ControlGridLibrary& rDummy,
    double start_x, double start_y, double start_z,
    std::size_t n_points_u, std::size_t n_points_v,
    double space1_x, double space1_y, double space1_z,
    double space2_x, double space2_y, double space2_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> space1(3);
    space1[0] = space1_x;
    space1[1] = space1_y;
    space1[2] = space1_z;

    std::vector<double> space2(3);
    space2[0] = space2_x;
    space2[1] = space2_y;
    space2[2] = space2_z;

    std::vector<std::vector<double> > spacing_vectors(2);
    spacing_vectors[0] = space1;
    spacing_vectors[1] = space2;

    return rDummy.CreateStructuredControlPointGrid<2>(start, ngrid, spacing_vectors);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid1(
    ControlGridLibrary& rDummy,
    double start_x, double start_y, double start_z,
    std::size_t n_points_u, std::size_t n_points_v, std::size_t n_points_w,
    double end_x, double end_y, double end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateStructuredControlPointGrid<3>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid2(
    ControlGridLibrary& rDummy,
    double start_x, double start_y, double start_z,
    std::size_t n_points_u, std::size_t n_points_v, std::size_t n_points_w,
    const pybind11::list& spacing_vectors_data)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<std::vector<double> > spacing_vectors;
    std::size_t cnt1 = 0, cnt2 = 0;
    for (auto vect : spacing_vectors_data)
    {
        std::vector<double> space_vect;
        for (auto v : vect.cast<pybind11::list>())
        {
            space_vect.push_back(v.cast<double>());
        }
        spacing_vectors.push_back(space_vect);
    }

    return rDummy.CreateStructuredControlPointGrid<3>(start, ngrid, spacing_vectors);
}

////////////////////////////////////////

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateLinearZeroControlGridWithVariable(
    ControlGridLibrary& rDummy,
    const TVariableType& rVariable,
    std::size_t n_points_u)
{
    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;
    return rDummy.CreateStructuredZeroControlGrid<1, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateRectangularZeroControlGridWithVariable(
    ControlGridLibrary& rDummy,
    const TVariableType& rVariable,
    std::size_t n_points_u, std::size_t n_points_v)
{
    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    return rDummy.CreateStructuredZeroControlGrid<2, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateCubicZeroControlGridWithVariable(
    ControlGridLibrary& rDummy,
    const TVariableType& rVariable,
    std::size_t n_points_u, std::size_t n_points_v, std::size_t n_points_w)
{
    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;
    return rDummy.CreateStructuredZeroControlGrid<3, TVariableType>(rVariable, ngrid);
}

template<typename TDataType, class TFESpaceType>
inline typename ControlGrid<TDataType>::Pointer ControlGridUtility_CreatePointBasedControlGrid(
    ControlGridUtility& rDummy,
    const Variable<TDataType>& rVariable, typename TFESpaceType::Pointer pFESpace)
{
    return rDummy.CreatePointBasedControlGrid<TDataType, TFESpaceType>(rVariable, pFESpace);
}

////////////////////////////////////////

void IsogeometricApplication_AddControlGridsToPython(pybind11::module& m)
{
    /////////////////////////////////////////////////////////////////
    ///////////////////////CONTROL GRIDS/////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<ControlGrid<ControlPoint<double> >, ControlGrid<ControlPoint<double> >::Pointer>
    (m, "ControlPointControlGrid")
    .def(init<>())
    .def("Size", &ControlGrid<ControlPoint<double> >::Size)
    .def("size", &ControlGrid<ControlPoint<double> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<ControlPoint<double> >)
    .def("__getitem__", &ControlGrid_GetItem<ControlPoint<double> >)
    .def("__str__", &PrintObject<ControlGrid<ControlPoint<double> > >)
    ;

    class_<ControlGrid<double>, ControlGrid<double>::Pointer>
    (m, "DoubleControlGrid")
    .def(init<>())
    .def("Size", &ControlGrid<double>::Size)
    .def("size", &ControlGrid<double>::Size)
    .def("__setitem__", &ControlGrid_SetItem<double>)
    .def("__getitem__", &ControlGrid_GetItem<double>)
    .def("__str__", &PrintObject<ControlGrid<double> >)
    ;

    class_<ControlGrid<array_1d<double, 3> >, ControlGrid<array_1d<double, 3> >::Pointer>
    (m, "Array1DControlGrid")
    .def(init<>())
    .def("Size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<array_1d<double, 3> >)
    .def("__getitem__", &ControlGrid_GetItem<array_1d<double, 3> >)
    .def("__str__", &PrintObject<ControlGrid<array_1d<double, 3> > >)
    ;

    class_<ControlGrid<Vector>, ControlGrid<Vector>::Pointer>
    (m, "VectorControlGrid")
    .def(init<>())
    .def("Size", &ControlGrid<Vector>::Size)
    .def("size", &ControlGrid<Vector>::Size)
    .def("__setitem__", &ControlGrid_SetItem<Vector>)
    .def("__getitem__", &ControlGrid_GetItem<Vector>)
    .def("__str__", &PrintObject<ControlGrid<Vector> >)
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<UnstructuredControlGrid<ControlPoint<double> >, UnstructuredControlGrid<ControlPoint<double> >::Pointer, ControlGrid<ControlPoint<double> > >
    (m, "UnstructuredControlPointGrid")
    .def(init<std::size_t>())
    .def("__str__", &PrintObject<UnstructuredControlGrid<ControlPoint<double> > >)
    ;

    class_<UnstructuredControlGrid<double>, UnstructuredControlGrid<double>::Pointer, ControlGrid<double> >
    (m, "UnstructuredDoubleControlGrid")
    .def(init<std::size_t>())
    .def("__str__", &PrintObject<UnstructuredControlGrid<double> >)
    ;

    class_<UnstructuredControlGrid<array_1d<double, 3> >, UnstructuredControlGrid<array_1d<double, 3> >::Pointer, ControlGrid<array_1d<double, 3> > >
    (m, "UnstructuredArray1DControlGrid")
    .def(init<std::size_t>())
    .def("__str__", &PrintObject<UnstructuredControlGrid<array_1d<double, 3> > >)
    ;

    class_<UnstructuredControlGrid<Vector>, UnstructuredControlGrid<Vector>::Pointer, ControlGrid<Vector> >
    (m, "UnstructuredVectorControlGrid")
    .def(init<std::size_t>())
    .def("__str__", &PrintObject<UnstructuredControlGrid<Vector> >)
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<ControlGridLibrary, ControlGridLibrary::Pointer>
    (m, "ControlGridLibrary")
    .def(init<>())
    .def("CreateLinearControlPointGrid", &ControlGridLibrary_CreateLinearControlPointGrid)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid1)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid2)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid1)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid2)
    .def("CreateLinearZeroDoubleControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<double> >)
    .def("CreateRectangularZeroDoubleControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<double> >)
    .def("CreateCubicZeroDoubleControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<double> >)
    .def("CreateLinearZeroArray1DControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateRectangularZeroArray1DControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateCubicZeroArray1DControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    ///////////////////////CONTROL GRID UTILITIY/////////////////////
    /////////////////////////////////////////////////////////////////

    class_<ControlGridUtility, ControlGridUtility::Pointer>
    (m, "PointBasedControlGridUtility")
    .def(init<>())
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBBSplinesFESpace<1, PBBSplinesBasisFunction<1, TCell>, BCellManager<1, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBBSplinesFESpace<2, PBBSplinesBasisFunction<2, TCell>, BCellManager<2, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, PBBSplinesFESpace<3, PBBSplinesBasisFunction<3, TCell>, BCellManager<3, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBBSplinesFESpace<1, PBBSplinesBasisFunction<1, TCell>, BCellManager<1, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBBSplinesFESpace<2, PBBSplinesBasisFunction<2, TCell>, BCellManager<2, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, PBBSplinesFESpace<3, PBBSplinesBasisFunction<3, TCell>, BCellManager<3, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBBSplinesFESpace<1, PBBSplinesBasisFunction<1, TCell>, BCellManager<1, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBBSplinesFESpace<2, PBBSplinesBasisFunction<2, TCell>, BCellManager<2, TCell> > >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, PBBSplinesFESpace<3, PBBSplinesBasisFunction<3, TCell>, BCellManager<3, TCell> > >)
    ////
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<3> >)
    ;
}

}  // namespace Python.

} // Namespace Kratos
