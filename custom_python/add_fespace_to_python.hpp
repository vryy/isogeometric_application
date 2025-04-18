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
#include <boost/python.hpp>

// Project includes
#include "python/python_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/weighted_fespace.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim>
void FESpace_ResetFunctionIndices(FESpace<TDim>& rDummy)
{
    rDummy.ResetFunctionIndices();
}

template<int TDim>
std::size_t FESpace_Enumerate(FESpace<TDim>& rDummy)
{
    std::size_t eq_size = 0;
    eq_size = rDummy.Enumerate(eq_size);
    return eq_size;
}

template<int TDim>
boost::python::list FESpace_GetValue(FESpace<TDim>& rDummy, boost::python::list xi_list)
{
    std::vector<double> xi;
    PythonUtils::Unpack<double, double>(xi_list, xi);

    std::vector<double> values = rDummy.GetValues(xi);

    boost::python::list values_list;
    for (std::size_t i = 0; i < values.size(); ++i)
        values_list.append(values[i]);

    return values_list;
}

template<int TDim>
bool FESpace_IsInside(FESpace<TDim>& rDummy, boost::python::list xi_list)
{
    std::vector<double> xi;
    PythonUtils::Unpack<double, double>(xi_list, xi);

    return rDummy.IsInside(xi);
}

template<int TDim>
boost::python::list FESpace_FunctionIndices(FESpace<TDim>& rDummy)
{
    boost::python::list indices;
    std::vector<std::size_t> all_indices = rDummy.FunctionIndices();
    for (std::size_t i = 0; i < all_indices.size(); ++i)
        indices.append<int>(static_cast<int>(all_indices[i]));
    return indices;
}

template<int TDim>
boost::python::list FESpace_BoundaryFunctionIndices(FESpace<TDim>& rDummy, int iside)
{
    boost::python::list indices;
    BoundarySide side = static_cast<BoundarySide>(iside);
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i]));
    return indices;
}

template<int TDim>
boost::python::list FESpace_BoundaryFunctionIndicesByFlag(FESpace<TDim>& rDummy, int boundary_id)
{
    boost::python::list indices;
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndicesByFlag(boundary_id);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i]));
    return indices;
}

template<int TDim>
boost::python::list FESpace_BoundaryShiftedFunctionIndices(FESpace<TDim>& rDummy, int iside)
{
    boost::python::list indices;
    BoundarySide side = static_cast<BoundarySide>(iside);
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i] + 1));
    return indices;
}

template<int TDim>
typename FESpace<TDim-1>::Pointer FESpace_ConstructBoundaryFESpace(FESpace<TDim>& rDummy, int iside)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.ConstructBoundaryFESpace(side);
}

template<int TDim>
void IsogeometricApplication_AddFESpacesToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "FESpace" << TDim << "D";
    class_<FESpace<TDim>, typename FESpace<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Order", &FESpace<TDim>::Order)
    .def("TotalNumber", &FESpace<TDim>::TotalNumber)
    .def("GetValue", &FESpace_GetValue<TDim>)
    .def("IsInside", &FESpace_IsInside<TDim>)
    .def("ResetFunctionIndices", &FESpace_ResetFunctionIndices<TDim>)
    .def("Reverse", &FESpace<TDim>::Reverse)
    .def("Enumerate", &FESpace_Enumerate<TDim>)
    .def("FunctionIndices", &FESpace_FunctionIndices<TDim>)
    .def("BoundaryFunctionIndices", &FESpace_BoundaryFunctionIndices<TDim>)
    .def("BoundaryFunctionIndicesByFlag", &FESpace_BoundaryFunctionIndicesByFlag<TDim>)
    .def("BoundaryShiftedFunctionIndices", &FESpace_BoundaryShiftedFunctionIndices<TDim>)
    .def("ConstructBoundaryFESpace", &FESpace_ConstructBoundaryFESpace<TDim>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "WeightedFESpace" << TDim << "D";
    class_<WeightedFESpace<TDim>, typename WeightedFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, const std::vector<double>&>())
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

