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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/grid_function.h"


namespace Kratos
{

namespace Python
{

template<class TGridFrunctionType>
typename TGridFrunctionType::FESpaceType::Pointer GridFunction_GetFESpace(TGridFrunctionType& rDummy)
{
    return rDummy.pFESpace();
}

template<class TGridFrunctionType>
void GridFunction_SetFESpace(TGridFrunctionType& rDummy, typename TGridFrunctionType::FESpaceType::Pointer pNewFESpace)
{
    rDummy.SetFESpace(pNewFESpace);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::ControlGridType::Pointer GridFunction_GetControlGrid(TGridFrunctionType& rDummy)
{
    return rDummy.pControlGrid();
}

template<class TGridFrunctionType>
void GridFunction_SetControlGrid(TGridFrunctionType& rDummy, typename TGridFrunctionType::ControlGridType::Pointer pNewControlGrid)
{
    rDummy.SetControlGrid(pNewControlGrid);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::DataType GridFunction_GetValue1(TGridFrunctionType& rDummy, const pybind11::list& py_xi)
{
    std::vector<double> xi_vec;
    for (auto v : py_xi)
        xi_vec.push_back(v.cast<double>());

    return rDummy.GetValue(xi_vec);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::DataType GridFunction_GetValue2(TGridFrunctionType& rDummy, const array_1d<double, 3>& xi)
{
    return rDummy.GetValue(xi);
}

template<class TGridFrunctionType>
pybind11::list GridFunction_GetDerivative1(TGridFrunctionType& rDummy, const pybind11::list& py_xi)
{
    std::vector<double> xi_vec;
    for (auto v : py_xi)
        xi_vec.push_back(v.cast<double>());

    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetDerivative(derivatives, xi_vec);

    pybind11::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType>
pybind11::list GridFunction_GetDerivative2(TGridFrunctionType& rDummy, const array_1d<double, 3>& xi)
{
    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetDerivative(derivatives, xi);

    pybind11::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType, typename TCoordinatesType>
pybind11::list GridFunction_LocalCoordinates(TGridFrunctionType& rDummy,
        const typename TGridFrunctionType::DataType& v, const TCoordinatesType& xi0)
{
    TCoordinatesType xi = xi0;
    int stat = rDummy.LocalCoordinates(v, xi);
    // KRATOS_WATCH(v)
    // KRATOS_WATCH(stat)
    // KRATOS_WATCH(xi)
    pybind11::list output;
    output.append(stat);
    output.append(xi);
    return output;
}

///////////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddGridFunctionsToPython(pybind11::module& m)
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "ControlPointGridFunction" << TDim << "D";
    typedef GridFunction<TDim, ControlPoint<double> > ControlPointGridFunctionType;
    pybind11::class_<ControlPointGridFunctionType, typename ControlPointGridFunctionType::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename FESpace<TDim>::Pointer, typename ControlGrid<ControlPoint<double> >::Pointer>())
    .def_property("FESpace", GridFunction_GetFESpace<ControlPointGridFunctionType>, GridFunction_SetFESpace<ControlPointGridFunctionType>)
    .def_property("ControlGrid", GridFunction_GetControlGrid<ControlPointGridFunctionType>, GridFunction_SetControlGrid<ControlPointGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<ControlPointGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<ControlPointGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<ControlPointGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<ControlPointGridFunctionType>)
    .def("__str__", &PrintObject<ControlPointGridFunctionType>)
    ;

    ss.str(std::string());
    ss << "DoubleGridFunction" << TDim << "D";
    typedef GridFunction<TDim, double> DoubleGridFunctionType;
    pybind11::class_<DoubleGridFunctionType, typename DoubleGridFunctionType::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename FESpace<TDim>::Pointer, typename ControlGrid<double>::Pointer>())
    .def_property("FESpace", GridFunction_GetFESpace<DoubleGridFunctionType>, GridFunction_SetFESpace<DoubleGridFunctionType>)
    .def_property("ControlGrid", GridFunction_GetControlGrid<DoubleGridFunctionType>, GridFunction_SetControlGrid<DoubleGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<DoubleGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<DoubleGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<DoubleGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<DoubleGridFunctionType>)
    .def("__str__", &PrintObject<DoubleGridFunctionType>)
    ;

    ss.str(std::string());
    ss << "Array1DGridFunction" << TDim << "D";
    typedef GridFunction<TDim, array_1d<double, 3> > Array1DGridFunctionType;
    pybind11::class_<Array1DGridFunctionType, typename Array1DGridFunctionType::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename FESpace<TDim>::Pointer, typename ControlGrid<array_1d<double, 3> >::Pointer>())
    .def_property("FESpace", GridFunction_GetFESpace<Array1DGridFunctionType>, GridFunction_SetFESpace<Array1DGridFunctionType>)
    .def_property("ControlGrid", GridFunction_GetControlGrid<Array1DGridFunctionType>, GridFunction_SetControlGrid<Array1DGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<Array1DGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<Array1DGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<Array1DGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<Array1DGridFunctionType>)
    .def("LocalCoordinates", &GridFunction_LocalCoordinates<Array1DGridFunctionType, array_1d<double, 3> >)
    .def("__str__", &PrintObject<Array1DGridFunctionType>)
    ;

    ss.str(std::string());
    ss << "VectorGridFunction" << TDim << "D";
    typedef GridFunction<TDim, Vector> VectorGridFunctionType;
    pybind11::class_<VectorGridFunctionType, typename VectorGridFunctionType::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename FESpace<TDim>::Pointer, typename ControlGrid<Vector>::Pointer>())
    .def_property("FESpace", GridFunction_GetFESpace<VectorGridFunctionType>, GridFunction_SetFESpace<VectorGridFunctionType>)
    .def_property("ControlGrid", GridFunction_GetControlGrid<VectorGridFunctionType>, GridFunction_SetControlGrid<VectorGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<VectorGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<VectorGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<VectorGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<VectorGridFunctionType>)
    .def("LocalCoordinates", &GridFunction_LocalCoordinates<VectorGridFunctionType, array_1d<double, 3> >)
    .def("__str__", &PrintObject<VectorGridFunctionType>)
    ;
}

}  // namespace Python.

} // Namespace Kratos

