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
#include "custom_utilities/multipatch_lagrange_control_mesh.h"
#include "custom_utilities/multipatch_model_part.h"
#include "custom_utilities/multi_multipatch_model_part.h"
#include "custom_utilities/conforming_multipatch_lagrange_model_part.h"
#include "custom_python/iga_python_utils.h"
#include "custom_python/add_mesh_and_model_part_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

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

template<class T>
void MultiPatchModelPart_BeginModelPart1(T& rDummy)
{
    rDummy.BeginModelPart();
}

template<class T>
void MultiPatchModelPart_BeginModelPart2(T& rDummy, ModelPart::Pointer pModelPart)
{
    rDummy.BeginModelPart(pModelPart);
}

template<class T>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions(T& rDummy,
        typename T::PatchType::Pointer pPatch,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pPatch, condition_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnBoundary(T& rDummy,
        typename T::PatchType::Pointer pPatch, int iside,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnBoundary2(T& rDummy,
        typename T::BoundaryPatchType::Pointer pBoundaryPatch,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pBoundaryPatch, condition_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiPatchModelPart_AddConditions_OnSlice(T& rDummy,
        typename T::PatchType::Pointer pPatch, int idir, double xi,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pPatch, idir, xi, condition_name, starting_id, pProperties);
}

template<class T, typename TVariableType>
void MultiPatchModelPart_SynchronizeBackward(T& rDummy,
        const TVariableType& rVariable, const boost::python::dict& rPatchNodalValues)
{
    std::map<std::size_t, std::map<std::size_t, typename TVariableType::Type> > patch_nodal_values;
    IsogeometricPythonUtils::Unpack(rPatchNodalValues, patch_nodal_values);

    rDummy.SynchronizeBackwardFromData(rVariable, patch_nodal_values);
}

////////////////////////////////////////

template<class T>
ModelPart::ElementsContainerType MultiMultiPatchModelPart_AddElements(T& rDummy, const boost::python::list& patch_list,
        const std::string& element_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    std::vector<typename T::PatchType::Pointer> pPatches;
    IsogeometricPythonUtils::Unpack<T>(patch_list, pPatches);

    return rDummy.AddElements(pPatches, element_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions(T& rDummy, const boost::python::list& patch_list,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    std::vector<typename T::PatchType::Pointer> pPatches;
    IsogeometricPythonUtils::Unpack<T>(patch_list, pPatches);

    return rDummy.AddConditions(pPatches, condition_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions_OnBoundary(T& rDummy,
        typename T::PatchType::Pointer pPatch, int iside,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.AddConditions(pPatch, side, condition_name, starting_id, pProperties);
}

template<class T>
ModelPart::ConditionsContainerType MultiMultiPatchModelPart_AddConditions_OnBoundary2(T& rDummy,
        typename T::BoundaryPatchType::Pointer pBoundaryPatch,
        const std::string& condition_name, std::size_t starting_id, Properties::Pointer pProperties)
{
    return rDummy.AddConditions(pBoundaryPatch, condition_name, starting_id, pProperties);
}

////////////////////////////////////////

template<int TDim>
boost::python::list PatchLagrangeMesh_WriteElements(PatchLagrangeMesh<TDim>& rDummy, ModelPart& r_model_part,
        typename Patch<TDim>::Pointer pPatch, const std::string& sample_element_name, const boost::python::list& list_divs,
        std::size_t last_node_id, std::size_t last_elem_id,
        Properties::Pointer pProperties, int echo_level)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
    {
        KRATOS_ERROR << "Element " << sample_element_name << " is not registered in Kratos."
                     << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);

    std::vector<std::size_t> num_divisions;
    IsogeometricPythonUtils::Unpack<int, std::size_t>(list_divs, num_divisions);

    std::size_t my_last_node_id = last_node_id;
    std::size_t my_last_elem_id = last_elem_id;

    rDummy.WriteEntities(r_model_part, r_model_part.Elements(),
                         pPatch, r_clone_element, num_divisions,
                         my_last_node_id, my_last_elem_id,
                         pProperties, echo_level);

    boost::python::list output;
    output.append(my_last_node_id);
    output.append(my_last_elem_id);
    return output;
}

template<int TDim>
boost::python::list PatchLagrangeMesh_WriteConditions(PatchLagrangeMesh<TDim>& rDummy, ModelPart& r_model_part,
        typename Patch<TDim>::Pointer pPatch, const std::string& sample_condition_name, const boost::python::list& list_divs,
        std::size_t last_node_id, std::size_t last_cond_id,
        Properties::Pointer pProperties, int echo_level)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
    {
        KRATOS_ERROR << "Condition " << sample_condition_name << " is not registered in Kratos."
                     << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

    std::vector<std::size_t> num_divisions;
    IsogeometricPythonUtils::Unpack<int, std::size_t>(list_divs, num_divisions);

    std::size_t my_last_node_id = last_node_id;
    std::size_t my_last_cond_id = last_cond_id;

    rDummy.WriteEntities(r_model_part, r_model_part.Conditions(),
                         pPatch, r_clone_condition, num_divisions,
                         my_last_node_id, my_last_cond_id,
                         pProperties, echo_level);

    boost::python::list output;
    output.append(my_last_node_id);
    output.append(my_last_cond_id);
    return output;
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddMeshToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "PatchLagrangeMesh" << TDim << "D";
    class_<PatchLagrangeMesh<TDim>, typename PatchLagrangeMesh<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("WriteElements", &PatchLagrangeMesh_WriteElements<TDim>)
    .def("WriteConditions", &PatchLagrangeMesh_WriteConditions<TDim>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "NonConformingMultipatchLagrangeMesh" << TDim << "D";
    class_<NonConformingMultipatchLagrangeMesh<TDim>, typename NonConformingMultipatchLagrangeMesh<TDim>::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("SetBaseElementName", &NonConformingMultipatchLagrangeMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &NonConformingMultipatchLagrangeMesh<TDim>::SetLastElemId)
    .def("SetDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetDivision)
    .def("SetUniformDivision", &NonConformingMultipatchLagrangeMesh<TDim>::SetUniformDivision)
    .def("WriteModelPart", &NonConformingMultipatchLagrangeMesh<TDim>::WriteModelPart)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "NonConformingVariableMultipatchLagrangeMesh" << TDim << "D";
    class_<NonConformingVariableMultipatchLagrangeMesh<TDim>, typename NonConformingVariableMultipatchLagrangeMesh<TDim>::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer, ModelPart&>())
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
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultipatchLagrangeControlMesh" << TDim << "D";
    class_<MultipatchLagrangeControlMesh<TDim>, typename MultipatchLagrangeControlMesh<TDim>::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("SetBaseElementName", &MultipatchLagrangeControlMesh<TDim>::SetBaseElementName)
    .def("SetLastNodeId", &MultipatchLagrangeControlMesh<TDim>::SetLastNodeId)
    .def("SetLastElemId", &MultipatchLagrangeControlMesh<TDim>::SetLastElemId)
    .def("WriteModelPart", &MultipatchLagrangeControlMesh<TDim>::WriteModelPart)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddModelPartToPython()
{
    std::stringstream ss;

    typedef MultiPatchModelPart<TDim> MultiPatchModelPartType;
    ss.str(std::string());
    ss << "MultiPatchModelPart" << TDim << "D";
    class_<MultiPatchModelPartType, typename MultiPatchModelPartType::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart1<MultiPatchModelPartType>)
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart2<MultiPatchModelPartType>)
    .def("CreateNodes", &MultiPatchModelPartType::CreateNodes)
    .def("AddElements", &MultiPatchModelPartType::AddElements)
    .def("AddConditions", &MultiPatchModelPart_AddConditions<MultiPatchModelPartType>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary<MultiPatchModelPartType>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary2<MultiPatchModelPartType>)
    .def("EndModelPart", &MultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiPatchModelPartType>, return_internal_reference<>())
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch<MultiPatchModelPartType>, return_internal_reference<>())
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeForward", &MultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiPatchModelPart_SynchronizeBackward<MultiPatchModelPartType, Variable<double> >)
    .def("SynchronizeBackward", &MultiPatchModelPart_SynchronizeBackward<MultiPatchModelPartType, Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiPatchModelPart_SynchronizeBackward<MultiPatchModelPartType, Variable<Vector> >)
    .def(self_ns::str(self))
    ;

    typedef MultiMultiPatchModelPart<TDim> MultiMultiPatchModelPartType;
    ss.str(std::string());
    ss << "MultiMultiPatchModelPart" << TDim << "D";
    class_<MultiMultiPatchModelPartType, typename MultiMultiPatchModelPartType::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("AddMultiPatch", &MultiMultiPatchModelPartType::AddMultiPatch)
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart1<MultiMultiPatchModelPartType>)
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart2<MultiMultiPatchModelPartType>)
    .def("CreateNodes", &MultiMultiPatchModelPartType::CreateNodes)
    .def("AddElements", &MultiMultiPatchModelPart_AddElements<MultiMultiPatchModelPartType>)
    .def("AddConditions", &MultiMultiPatchModelPart_AddConditions<MultiMultiPatchModelPartType>)
    .def("AddConditions", &MultiMultiPatchModelPart_AddConditions_OnBoundary<MultiMultiPatchModelPartType>)
    .def("AddConditions", &MultiMultiPatchModelPart_AddConditions_OnBoundary2<MultiMultiPatchModelPartType>)
    .def("EndModelPart", &MultiMultiPatchModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<MultiMultiPatchModelPartType>, return_internal_reference<>())
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch2<MultiMultiPatchModelPartType>, return_internal_reference<>())
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<double> >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<double> >)
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<array_1d<double, 3> > >)
    .def("SynchronizeForward", &MultiMultiPatchModelPartType::template SynchronizeForward<Variable<Vector> >)
    .def("SynchronizeBackward", &MultiMultiPatchModelPartType::template SynchronizeBackward<Variable<Vector> >)
    .def(self_ns::str(self))
    ;

    typedef ConformingMultipatchLagrangeModelPart<TDim> ConformingMultipatchLagrangeModelPartType;
    ss.str(std::string());
    ss << "ConformingMultipatchLagrangeModelPart" << TDim << "D";
    class_<ConformingMultipatchLagrangeModelPartType, typename ConformingMultipatchLagrangeModelPartType::Pointer, bases<IsogeometricEcho>, boost::noncopyable>
    (ss.str().c_str(), init<typename MultiPatch<TDim>::Pointer>())
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart1<ConformingMultipatchLagrangeModelPartType>)
    .def("BeginModelPart", &MultiPatchModelPart_BeginModelPart2<ConformingMultipatchLagrangeModelPartType>)
    .def("CreateNodes", &ConformingMultipatchLagrangeModelPartType::CreateNodes)
    .def("AddElements", &ConformingMultipatchLagrangeModelPartType::AddElements)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnBoundary<ConformingMultipatchLagrangeModelPartType>)
    .def("AddConditions", &MultiPatchModelPart_AddConditions_OnSlice<ConformingMultipatchLagrangeModelPartType>)
    .def("EndModelPart", &ConformingMultipatchLagrangeModelPartType::EndModelPart)
    .def("GetModelPart", &MultiPatchModelPart_GetModelPart<ConformingMultipatchLagrangeModelPartType>, return_internal_reference<>())
    .def("GetMultiPatch", &MultiPatchModelPart_GetMultiPatch<ConformingMultipatchLagrangeModelPartType>, return_internal_reference<>())
    .def("SetUniformSampling", &ConformingMultipatchLagrangeModelPartType::SetUniformSampling)
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddMeshAndModelPartToPython()
{

    IsogeometricApplication_AddMeshToPython<2>();
    IsogeometricApplication_AddMeshToPython<3>();

    IsogeometricApplication_AddModelPartToPython<1>();
    IsogeometricApplication_AddModelPartToPython<2>();
    IsogeometricApplication_AddModelPartToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos
