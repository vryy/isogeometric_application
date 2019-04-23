/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 21, 2017 $
//   Revision:            $Revision: 1.0 $
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
#include "includes/model_part.h"
#include "custom_utilities/nonconforming_multipatch_lagrange_mesh.h"
#include "custom_utilities/nonconforming_variable_multipatch_lagrange_mesh.h"
#include "custom_utilities/multipatch_model_part.h"
// #include "custom_utilities/multi_multipatch_model_part.h"
#include "custom_python/add_mesh_and_model_part_to_python.h"


namespace Kratos
{

namespace Python
{

////////////////////////////////////////

template<class T>
ModelPart& MultiPatchModelPart_GetModelPart(T& rDummy)
{
    return rDummy.GetModelPart();
}

template<class T>
typename T::MultiPatchType& MultiPatchModelPart_GetMultiPatch(T& rDummy)
{
    return *(rDummy.pMultiPatch());
}

template<class T>
typename T::MultiPatchType& MultiPatchModelPart_GetMultiPatch2(T& rDummy, const std::size_t& i)
{
    return *(rDummy.pMultiPatch(i));
}

template<int TDim>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions(MultiPatchModelPart<TDim>& rDummy,
    typename Patch<TDim>::Pointer pPatch,
    const std::string& condition_name, const std::size_t& starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pPatch, condition_name, starting_id, pProperties);
}

template<int TDim>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnBoundary(MultiPatchModelPart<TDim>& rDummy,
    typename Patch<TDim>::Pointer pPatch, const int& iside,
    const std::string& condition_name, const std::size_t& starting_id, Properties::Pointer pProperties)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
}

////////////////////////////////////////

// template<class T>
// ModelPart::ElementsContainerType MultiMultiPatchModelPart_AddElements(T& rDummy, pybind11::list py_patch_list,
//     const std::string& element_name, const std::size_t& starting_id, Properties::Pointer pProperties)
// {
//     std::vector<typename T::PatchType::Pointer> pPatches;
//     for (auto v : py_patch_list)
//         pPatches.push_back(v.cast<typename T::PatchType::Pointer>());

//     return rDummy.AddElements(pPatches, element_name, starting_id, pProperties);
// }

// template<class T>
// ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions(T& rDummy, pybind11::list py_patch_list,
//     const std::string& condition_name, const std::size_t& starting_id, Properties::Pointer pProperties)
// {
//     std::vector<typename T::PatchType::Pointer> pPatches;
//     for (auto v : py_patch_list)
//         pPatches.push_back(v.cast<typename T::PatchType::Pointer>());

//     return rDummy.AddConditions(pPatches, condition_name, starting_id, pProperties);
// }

// template<int TDim>
// ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions_OnBoundary(MultiMultiPatchModelPart<TDim>& rDummy,
//     typename Patch<TDim>::Pointer pPatch, const int& iside,
//     const std::string& condition_name, const std::size_t& starting_id, Properties::Pointer pProperties)
// {
//     BoundarySide side = static_cast<BoundarySide>(iside);
//     return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
// }

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddMeshToPython(pybind11::module& m)
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "NonConformingMultipatchLagrangeMesh" << TDim << "D";
    pybind11::class_<NonConformingMultipatchLagrangeMesh<TDim>, typename NonConformingMultipatchLagrangeMesh<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename MultiPatch<TDim>::Pointer>())
    .def("SetBaseElementName", &NonConformingMultipatchLagrangeMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastElemId)
    .def("SetLastPropId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastPropId)
    .def("SetDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetDivision)
    .def("SetUniformDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetUniformDivision)
    .def("WriteModelPart", &NonConformingMultipatchLagrangeMesh<TDim>::WriteModelPart)
    .def("__str__", &PrintObject<NonConformingMultipatchLagrangeMesh<TDim> >)
    ;

    ss.str(std::string());
    ss << "NonConformingVariableMultipatchLagrangeMesh" << TDim << "D";
    pybind11::class_<NonConformingVariableMultipatchLagrangeMesh<TDim>, typename NonConformingVariableMultipatchLagrangeMesh<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename MultiPatch<TDim>::Pointer, ModelPart&>())
    .def("SetBaseElementName", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetLastElemId)
    .def("SetLastPropId", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetLastPropId)
    .def("SetDivision", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetDivision)
    .def("SetUniformDivision", &NonConformingVariableMultipatchLagrangeMesh<TDim>::SetUniformDivision)
    .def("WriteModelPart", &NonConformingVariableMultipatchLagrangeMesh<TDim>::WriteModelPart)
    .def("TransferVariables", &NonConformingVariableMultipatchLagrangeMesh<TDim>::template TransferVariables<Variable<double> >)
    .def("TransferVariables", &NonConformingVariableMultipatchLagrangeMesh<TDim>::template TransferVariables<Variable<array_1d<double, 3> > >)
    .def("TransferVariables", &NonConformingVariableMultipatchLagrangeMesh<TDim>::template TransferVariables<Variable<Vector> >)
    .def("__str__", &PrintObject<NonConformingVariableMultipatchLagrangeMesh<TDim> >)
    ;
}

template<int TDim>
void IsogeometricApplication_AddModelPartToPython(pybind11::module& m)
{
    std::stringstream ss;

    typedef MultiPatchModelPart<TDim> MultiPatchModelPartType;
    ss.str(std::string());
    ss << "MultiPatchModelPart" << TDim << "D";
    pybind11::class_<MultiPatchModelPartType, typename MultiPatchModelPartType::Pointer>
    (m, ss.str().c_str())
    .def(pybind11::init<typename MultiPatch<TDim>::Pointer>())
    .def("BeginModelPart", &MultiPatchModelPartType::BeginModelPart)
    .def("CreateNodes", &MultiPatchModelPartType::CreateNodes)
    .def("AddElements", &MultiPatchModelPartType::AddElements)
    .def("AddConditions", &MultiPatchModelPart_AddConditions<TDim>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary<TDim>)
    .def("EndModelPart", &MultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiPatchModelPartType>, pybind11::return_value_policy::reference)
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch<MultiPatchModelPartType>, pybind11::return_value_policy::reference)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    .def("__str__", &PrintObject<MultiPatchModelPartType>)
    ;

    // typedef MultiMultiPatchModelPart<TDim> MultiMultiPatchModelPartType;
    // ss.str(std::string());
    // ss << "MultiMultiPatchModelPart" << TDim << "D";
    // pybind11::class_<MultiMultiPatchModelPartType, typename MultiMultiPatchModelPartType::Pointer>
    // (m, ss.str().c_str())
    // .def(pybind11::init<>())
    // .def("AddMultiPatch", &MultiMultiPatchModelPartType::AddMultiPatch)
    // .def("BeginModelPart", &MultiMultiPatchModelPartType::BeginModelPart)
    // .def("CreateNodes", &MultiMultiPatchModelPartType::CreateNodes)
    // .def("AddElements", &MultiMultiPatchModelPart_AddElements<MultiMultiPatchModelPartType>)
    // .def("AddConditions", &MultiMultiPatchModelPart_AddConditions<MultiMultiPatchModelPartType>)
    // .def("AddConditions", &MultiMultiPatchModelPart_AddConditions_OnBoundary<TDim>)
    // .def("EndModelPart", &MultiMultiPatchModelPartType::EndModelPart)
    // .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiMultiPatchModelPartType>, pybind11::return_value_policy::reference)
    // .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch2<MultiMultiPatchModelPartType>, pybind11::return_value_policy::reference)
    // .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    // .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    // .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    // .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    // .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    // .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    // .def("__str__", &PrintObject<MultiMultiPatchModelPartType>)
    // ;
}


void IsogeometricApplication_AddMeshAndModelPartToPython(pybind11::module& m)
{

    IsogeometricApplication_AddMeshToPython<2>(m);
    IsogeometricApplication_AddMeshToPython<3>(m);

    IsogeometricApplication_AddModelPartToPython<2>(m);
    IsogeometricApplication_AddModelPartToPython<3>(m);

}

}  // pybind11 Python.

} // pybind11 Kratos

