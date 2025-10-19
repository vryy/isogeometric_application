/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "python/python_utils.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/grid_function.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

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
typename TGridFrunctionType::DataType GridFunction_GetValue1(TGridFrunctionType& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    PythonUtils::Unpack<double, double>(xi, xi_vec);
    return rDummy.GetValue(xi_vec);
}

template<class TGridFrunctionType>
typename TGridFrunctionType::DataType GridFunction_GetValue2(TGridFrunctionType& rDummy, const array_1d<double, 3>& xi)
{
    return rDummy.GetValue(xi);
}

template<class TGridFrunctionType>
boost::python::list GridFunction_GetDerivative1(TGridFrunctionType& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    PythonUtils::Unpack<double, double>(xi, xi_vec);

    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetDerivative(derivatives, xi_vec);

    boost::python::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType>
boost::python::list GridFunction_GetDerivative2(TGridFrunctionType& rDummy, const array_1d<double, 3>& xi)
{
    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetDerivative(derivatives, xi);

    boost::python::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType>
boost::python::list GridFunction_GetSecondDerivative1(TGridFrunctionType& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    PythonUtils::Unpack<double, double>(xi, xi_vec);

    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetSecondDerivative(derivatives, xi_vec);

    boost::python::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType>
boost::python::list GridFunction_GetSecondDerivative2(TGridFrunctionType& rDummy, const array_1d<double, 3>& xi)
{
    std::vector<typename TGridFrunctionType::DataType> derivatives;
    rDummy.GetSecondDerivative(derivatives, xi);

    boost::python::list results;
    for (std::size_t i = 0; i < derivatives.size(); ++i)
        results.append(derivatives[i]);

    return results;
}

template<class TGridFrunctionType, typename TLocalCoordinatesType>
boost::python::list GridFunction_LocalCoordinates(TGridFrunctionType& rDummy,
        const typename TGridFrunctionType::DataType& v, const TLocalCoordinatesType& xi0)
{
    TLocalCoordinatesType xi = xi0;
    int stat = rDummy.LocalCoordinates(v, xi);
    // KRATOS_WATCH(v)
    // KRATOS_WATCH(stat)
    // KRATOS_WATCH(xi)
    boost::python::list output;
    output.append(stat);
    output.append(xi);
    return output;
}

///////////////////////////////////////////////////////

template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TWeightType>
void IsogeometricApplication_AddControlPointGridFunctionToPython(const std::string& Prefix)
{
    std::stringstream ss;

    typedef ControlPoint<TCoordinateType, TWeightType> ControlPointType;
    typedef GridFunction<TDim, TLocalCoordinateType, ControlPointType> ControlPointGridFunctionType;
    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;

    ss.str(std::string());
    ss << Prefix << "ControlPointGridFunction" << TDim << "D";
    class_<ControlPointGridFunctionType, typename ControlPointGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpaceType::Pointer, typename ControlGrid<ControlPointType>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<ControlPointGridFunctionType>, GridFunction_SetFESpace<ControlPointGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<ControlPointGridFunctionType>, GridFunction_SetControlGrid<ControlPointGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<ControlPointGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<ControlPointGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<ControlPointGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<ControlPointGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative1<ControlPointGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative2<ControlPointGridFunctionType>)
    .def(self_ns::str(self))
    ;
}

template<int TDim, typename TLocalCoordinateType, typename TDataType>
void IsogeometricApplication_AddGridFunctionsToPython(const std::string& Prefix)
{
    typedef typename MatrixVectorTypeSelector<TDataType>::VectorType VectorType;

    typedef GridFunction<TDim, TLocalCoordinateType, TDataType> DoubleGridFunctionType;
    typedef GridFunction<TDim, TLocalCoordinateType, array_1d<TDataType, 3> > Array1DGridFunctionType;
    typedef GridFunction<TDim, TLocalCoordinateType, VectorType> VectorGridFunctionType;

    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;

    std::stringstream ss;

    ss.str(std::string());
    ss << Prefix << "DoubleGridFunction" << TDim << "D";
    class_<DoubleGridFunctionType, typename DoubleGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpaceType::Pointer, typename ControlGrid<TDataType>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<DoubleGridFunctionType>, GridFunction_SetFESpace<DoubleGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<DoubleGridFunctionType>, GridFunction_SetControlGrid<DoubleGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<DoubleGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<DoubleGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<DoubleGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<DoubleGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative1<DoubleGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative2<DoubleGridFunctionType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << Prefix << "Array1DGridFunction" << TDim << "D";
    class_<Array1DGridFunctionType, typename Array1DGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpaceType::Pointer, typename ControlGrid<array_1d<TDataType, 3> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<Array1DGridFunctionType>, GridFunction_SetFESpace<Array1DGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<Array1DGridFunctionType>, GridFunction_SetControlGrid<Array1DGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<Array1DGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<Array1DGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<Array1DGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<Array1DGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative1<Array1DGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative2<Array1DGridFunctionType>)
    .def("LocalCoordinates", &GridFunction_LocalCoordinates<Array1DGridFunctionType, array_1d<TLocalCoordinateType, 3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << Prefix << "VectorGridFunction" << TDim << "D";
    class_<VectorGridFunctionType, typename VectorGridFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpaceType::Pointer, typename ControlGrid<VectorType>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<VectorGridFunctionType>, GridFunction_SetFESpace<VectorGridFunctionType>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<VectorGridFunctionType>, GridFunction_SetControlGrid<VectorGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue1<VectorGridFunctionType>)
    .def("GetValue", &GridFunction_GetValue2<VectorGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative1<VectorGridFunctionType>)
    .def("GetDerivative", &GridFunction_GetDerivative2<VectorGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative1<VectorGridFunctionType>)
    .def("GetSecondDerivative", &GridFunction_GetSecondDerivative2<VectorGridFunctionType>)
    .def("LocalCoordinates", &GridFunction_LocalCoordinates<VectorGridFunctionType, array_1d<TLocalCoordinateType, 3> >)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos
