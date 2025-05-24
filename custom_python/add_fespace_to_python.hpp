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

template<class TFESpaceType>
void FESpace_ResetFunctionIndices(TFESpaceType& rDummy)
{
    rDummy.ResetFunctionIndices();
}

template<class TFESpaceType>
std::size_t FESpace_Enumerate(TFESpaceType& rDummy)
{
    std::size_t eq_size = 0;
    eq_size = rDummy.Enumerate(eq_size);
    return eq_size;
}

template<class TFESpaceType>
boost::python::list FESpace_GetValue(TFESpaceType& rDummy, const boost::python::list& xi_list)
{
    std::vector<typename TFESpaceType::LocalCoordinateType> xi;
    PythonUtils::Unpack<double, typename TFESpaceType::LocalCoordinateType>(xi_list, xi);

    std::vector<double> values = rDummy.GetValues(xi);

    boost::python::list values_list;
    for (std::size_t i = 0; i < values.size(); ++i)
        values_list.append(values[i]);

    return values_list;
}

template<class TFESpaceType>
bool FESpace_IsInside(TFESpaceType& rDummy, const boost::python::list& xi_list)
{
    std::vector<typename TFESpaceType::LocalCoordinateType> xi;
    PythonUtils::Unpack<double, typename TFESpaceType::LocalCoordinateType>(xi_list, xi);

    return rDummy.IsInside(xi);
}

template<class TFESpaceType>
boost::python::list FESpace_FunctionIndices(TFESpaceType& rDummy)
{
    boost::python::list indices;
    std::vector<std::size_t> all_indices = rDummy.FunctionIndices();
    for (std::size_t i = 0; i < all_indices.size(); ++i)
        indices.append<int>(static_cast<int>(all_indices[i]));
    return indices;
}

template<class TFESpaceType>
boost::python::list FESpace_BoundaryFunctionIndices(TFESpaceType& rDummy, int iside)
{
    boost::python::list indices;
    BoundarySide side = static_cast<BoundarySide>(iside);
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i]));
    return indices;
}

template<class TFESpaceType>
boost::python::list FESpace_BoundaryFunctionIndicesByFlag(TFESpaceType& rDummy, int boundary_id)
{
    boost::python::list indices;
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndicesByFlag(boundary_id);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i]));
    return indices;
}

template<class TFESpaceType>
boost::python::list FESpace_BoundaryShiftedFunctionIndices(TFESpaceType& rDummy, int iside)
{
    boost::python::list indices;
    BoundarySide side = static_cast<BoundarySide>(iside);
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i] + 1));
    return indices;
}

template<class TFESpaceType>
typename TFESpaceType::BoundaryFESpaceType::Pointer FESpace_ConstructBoundaryFESpace(TFESpaceType& rDummy, int iside)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.ConstructBoundaryFESpace(side);
}

template<int TDim, typename TLocalCoordinateType>
void IsogeometricApplication_AddFESpacesToPython()
{

    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;
    typedef WeightedFESpace<TDim, TLocalCoordinateType> WeightedFESpaceType;

    std::stringstream ss;

    ss.str(std::string());
    ss << "FESpace" << TDim << "D";
    class_<FESpaceType, typename FESpaceType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Order", &FESpaceType::Order)
    .def("TotalNumber", &FESpaceType::TotalNumber)
    .def("GetValue", &FESpace_GetValue<FESpaceType>)
    .def("IsInside", &FESpace_IsInside<FESpaceType>)
    .def("ResetFunctionIndices", &FESpace_ResetFunctionIndices<FESpaceType>)
    .def("Reverse", &FESpaceType::Reverse)
    .def("Enumerate", &FESpace_Enumerate<FESpaceType>)
    .def("FunctionIndices", &FESpace_FunctionIndices<FESpaceType>)
    .def("BoundaryFunctionIndices", &FESpace_BoundaryFunctionIndices<FESpaceType>)
    .def("BoundaryFunctionIndicesByFlag", &FESpace_BoundaryFunctionIndicesByFlag<FESpaceType>)
    .def("BoundaryShiftedFunctionIndices", &FESpace_BoundaryShiftedFunctionIndices<FESpaceType>)
    .def("ConstructBoundaryFESpace", &FESpace_ConstructBoundaryFESpace<FESpaceType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "WeightedFESpace" << TDim << "D";
    class_<WeightedFESpaceType, typename WeightedFESpaceType::Pointer, bases<FESpaceType>, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpaceType::Pointer, const std::vector<double>&>())
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

