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

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/patch_lagrange_mesh.h"
#include "custom_utilities/nonconforming_multipatch_lagrange_mesh.h"
#include "custom_utilities/nonconforming_variable_multipatch_lagrange_mesh.h"
#include "custom_utilities/multipatch_model_part.h"
//#include "custom_utilities/multi_multipatch_model_part.h"
#include "custom_python3/add_mesh_and_model_part_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

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
typename T::MultiPatchType& MultiPatchModelPart_GetMultiPatch2(T& rDummy, std::size_t i)
{
    return *(rDummy.pMultiPatch(i));
}

template<int TDim>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions(MultiPatchModelPart<TDim>& rDummy,
    typename Patch<TDim>::Pointer pPatch,
    const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pPatch, condition_name, starting_id, pProperties);
}

template<int TDim>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnBoundary(MultiPatchModelPart<TDim>& rDummy,
    typename Patch<TDim>::Pointer pPatch, int iside,
    const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
}

template<int TDim>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnBoundary2(MultiPatchModelPart<TDim>& rDummy,
    typename Patch<TDim-1>::Pointer pBoundaryPatch,
    const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pBoundaryPatch, condition_name, starting_id, pProperties);
}

// template<int TDim>
// ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions_OnBoundary(MultiMultiPatchModelPart<TDim>& rDummy,
//     typename Patch<TDim>::Pointer pPatch, int iside,
//     const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
// {
//     BoundarySide side = static_cast<BoundarySide>(iside);
//     return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
// }

// template<int TDim>
// ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions_OnBoundary2(MultiMultiPatchModelPart<TDim>& rDummy,
//     typename Patch<TDim-1>::Pointer pBoundaryPatch,
//     const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
// {
//     return rDummy.AddConditions(pBoundaryPatch, condition_name, starting_id, pProperties);
// }

// ////////////////////////////////////////

// template<class T>
// ModelPart::ElementsContainerType MultiMultiPatchModelPart_AddElements(T& rDummy, pybind11::list& patch_list,
//     const std::string& element_name, std::size_t starting_id, Properties::Pointer pProperties)
// {
//     std::vector<typename T::PatchType::Pointer> pPatches;

//     for (auto v : patch_list)
//         pPatches.push_back(v.cast<typename T::PatchType::Pointer>());

//     return rDummy.AddElements(pPatches, element_name, starting_id, pProperties);
// }

// template<class T>
// ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions(T& rDummy, pybind11::list& patch_list,
//     const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
// {
//     std::vector<typename T::PatchType::Pointer> pPatches;

//     for (auto v : patch_list)
//         pPatches.push_back(v.cast<typename T::PatchType::Pointer>());

//     return rDummy.AddConditions(pPatches, condition_name, starting_id, pProperties);
// }

////////////////////////////////////////

template<int TDim>
pybind11::list PatchLagrangeMesh_WriteElements(PatchLagrangeMesh<TDim>& rDummy, ModelPart& r_model_part,
    typename Patch<TDim>::Pointer pPatch, const std::string& sample_element_name, pybind11::list& list_divs,
    std::size_t last_node_id, std::size_t last_elem_id,
    Properties::Pointer pProperties, int echo_level)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << sample_element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
        KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    std::vector<std::size_t> num_divisions;
    for (auto v : list_divs)
        num_divisions.push_back(v.cast<int>());

    std::size_t my_last_node_id = last_node_id;
    std::size_t my_last_elem_id = last_elem_id;

    rDummy.WriteEntities(r_model_part, r_model_part.Elements(),
        pPatch, r_clone_element, num_divisions,
        my_last_node_id, my_last_elem_id,
        pProperties, echo_level);

    pybind11::list output;
    output.append(my_last_node_id);
    output.append(my_last_elem_id);
    return output;
}

template<int TDim>
pybind11::list PatchLagrangeMesh_WriteConditions(PatchLagrangeMesh<TDim>& rDummy, ModelPart& r_model_part,
    typename Patch<TDim>::Pointer pPatch, const std::string& sample_condition_name, pybind11::list& list_divs,
    std::size_t last_node_id, std::size_t last_cond_id,
    Properties::Pointer pProperties, int echo_level)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
    {
        std::stringstream buffer;
        buffer << "Condition " << sample_condition_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
        KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

    std::vector<std::size_t> num_divisions;
    for (auto v : list_divs)
        num_divisions.push_back(v.cast<int>());

    std::size_t my_last_node_id = last_node_id;
    std::size_t my_last_cond_id = last_cond_id;

    rDummy.WriteEntities(r_model_part, r_model_part.Conditions(),
        pPatch, r_clone_condition, num_divisions,
        my_last_node_id, my_last_cond_id,
        pProperties, echo_level);

    pybind11::list output;
    output.append(my_last_node_id);
    output.append(my_last_cond_id);
    return output;
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddMeshToPython(pybind11::module& m)
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "PatchLagrangeMesh" << TDim << "D";
    class_<PatchLagrangeMesh<TDim>, typename PatchLagrangeMesh<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<>())
    .def("WriteElements", &PatchLagrangeMesh_WriteElements<TDim>)
    .def("WriteConditions", &PatchLagrangeMesh_WriteConditions<TDim>)
    .def("__str__", &PrintObject<PatchLagrangeMesh<TDim>>)
    ;

    ss.str(std::string());
    ss << "NonConformingMultipatchLagrangeMesh" << TDim << "D";
    class_<NonConformingMultipatchLagrangeMesh<TDim>, typename NonConformingMultipatchLagrangeMesh<TDim>::Pointer, IsogeometricEcho>
    (m, ss.str().c_str())
    .def(init<typename MultiPatch<TDim>::Pointer>())
    .def("SetEchoLevel", &NonConformingMultipatchLagrangeMesh<TDim>::SetEchoLevel)
    .def("SetBaseElementName", &NonConformingMultipatchLagrangeMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastElemId)
    .def("SetDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetDivision)
    .def("SetUniformDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetUniformDivision)
    .def("WriteModelPart", &NonConformingMultipatchLagrangeMesh<TDim>::WriteModelPart)
    .def("__str__", &PrintObject<NonConformingMultipatchLagrangeMesh<TDim> >)
    ;

    ss.str(std::string());
    ss << "NonConformingVariableMultipatchLagrangeMesh" << TDim << "D";
    class_<NonConformingVariableMultipatchLagrangeMesh<TDim>, typename NonConformingVariableMultipatchLagrangeMesh<TDim>::Pointer, IsogeometricEcho>
    (m, ss.str().c_str())
    .def(init<typename MultiPatch<TDim>::Pointer, ModelPart&>())
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
    class_<MultiPatchModelPartType, typename MultiPatchModelPartType::Pointer, IsogeometricEcho>
    (m, ss.str().c_str())
    .def(init<typename MultiPatch<TDim>::Pointer>())
    .def("BeginModelPart", &MultiPatchModelPartType::BeginModelPart)
    .def("CreateNodes", &MultiPatchModelPartType::CreateNodes)
    .def("AddElements", &MultiPatchModelPartType::AddElements)
    .def("AddConditions", &MultiPatchModelPart_AddConditions<TDim>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary<TDim>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary2<TDim>)
    .def("EndModelPart", &MultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiPatchModelPartType>, return_value_policy::reference_internal)
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch<MultiPatchModelPartType>, return_value_policy::reference_internal)
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
    // class_<MultiMultiPatchModelPartType, typename MultiMultiPatchModelPartType::Pointer, IsogeometricEcho>
    // (m, ss.str().c_str())
    // .def(init<>())
    // .def("AddMultiPatch", &MultiMultiPatchModelPartType::AddMultiPatch)
    // .def("BeginModelPart", &MultiMultiPatchModelPartType::BeginModelPart)
    // .def("CreateNodes", &MultiMultiPatchModelPartType::CreateNodes)
    // .def("AddElements", &MultiMultiPatchModelPart_AddElements<MultiMultiPatchModelPartType>)
    // .def("AddConditions", &MultiMultiPatchModelPart_AddConditions<MultiMultiPatchModelPartType>)
    // .def("AddConditions", &MultiMultiPatchModelPart_AddConditions_OnBoundary<TDim>)
    // .def("AddConditions", &MultiMultiPatchModelPart_AddConditions_OnBoundary2<TDim>)
    // .def("EndModelPart", &MultiMultiPatchModelPartType::EndModelPart)
    // .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiMultiPatchModelPartType>, return_value_policy::reference_internal)
    // .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch2<MultiMultiPatchModelPartType>, return_value_policy::reference_internal)
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

    IsogeometricApplication_AddModelPartToPython<1>(m);
    IsogeometricApplication_AddModelPartToPython<2>(m);
    IsogeometricApplication_AddModelPartToPython<3>(m);

}

}  // namespace Python.

} // Namespace Kratos

