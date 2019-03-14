/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 9, 2013 $
//   Revision:            $Revision: 1.1 $
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
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/isogeometric_post_utility.h"
#include "custom_utilities/bezier_classical_post_utility.h"
#include "custom_utilities/bezier_post_utility.h"
#include "custom_utilities/nurbs_test_utils.h"
#include "custom_utilities/bezier_test_utils.h"
#include "custom_utilities/isogeometric_merge_utility.h"

#ifdef ISOGEOMETRIC_USE_HDF5
#include "custom_utilities/hdf5_post_utility.h"
#endif

#ifdef ISOGEOMETRIC_USE_GISMO
#include "custom_utilities/gismo/gismo_mesh.h"
#endif

namespace Kratos
{

namespace Python
{

using namespace boost::python;

int BSplineUtils_FindSpan(
    BSplineUtils& dummy,
    const int rN,
    const int rP,
    const double rXi,
    const Vector& rU
)
{
    return dummy.FindSpan(rN, rP, rXi, rU);
}

void BSplineUtils_BasisFuns(
    BSplineUtils& dummy,
    Vector& rS,
    const int rI,
    const double rXi,
    const int rP,
    const Vector& rU
)
{
    dummy.BasisFuns(rS, rI, rXi, rP, rU);
}

void BezierUtils_Bernstein(
    BezierUtils& dummy,
    BezierUtils::ValuesContainerType& rS,
    const int p,
    const double x
)
{
    dummy.bernstein(rS, p, x);
}

void BezierUtils_Bernstein_der(
    BezierUtils& dummy,
    BezierUtils::ValuesContainerType& rS,
    BezierUtils::ValuesContainerType& rD,
    const int p,
    const double x
)
{
    dummy.bernstein(rS, rD, p, x);
}

void BezierUtils_DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(
    BezierUtils& dummy,
    ModelPart::Pointer pModelPart,
    std::string FileName
)
{
    dummy.DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(pModelPart, FileName);
}

template<class T>
void BezierUtils_ComputeCentroid(
    BezierUtils& dummy,
    typename T::Pointer& pElem,
    typename T::GeometryType::PointType& P
)
{
    dummy.ComputeCentroid<T>(pElem, P);
}

void NURBSTestUtils_ProbeGlobalCoordinates2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeGlobalCoordinates(pElement, X, Y);
}

void NURBSTestUtils_ProbeGlobalCoordinates3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeGlobalCoordinates(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeShapeFunctionValues2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeShapeFunctionValues(pElement, X, Y);
}

void NURBSTestUtils_ProbeShapeFunctionValues3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeShapeFunctionValues(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeShapeFunctionDerivatives2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement, X, Y);
}

void NURBSTestUtils_ProbeShapeFunctionDerivatives3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement, X, Y, Z);
}

void NURBSTestUtils_ProbeJacobian2(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y
)
{
    dummy.ProbeJacobian(pElement, X, Y);
}

void NURBSTestUtils_ProbeJacobian3(
    NURBSTestUtils& dummy,
    Element::Pointer& pElement,
    double X,
    double Y,
    double Z
)
{
    dummy.ProbeJacobian(pElement, X, Y, Z);
}

//////////////////////////////////////////////////////////

ModelPart::ElementsContainerType IsogeometricPostUtility_TransferElements(IsogeometricPostUtility& rDummy, ModelPart::ElementsContainerType& pElements,
    ModelPart& r_other_model_part, const std::string& sample_element_name, Properties::Pointer pProperties, const bool& retain_prop_id)
{
    std::size_t last_element_id = IsogeometricPostUtility::GetLastElementId(r_other_model_part);
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_element_name, "is not registered to the Kratos kernel")
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    ModelPart::ElementsContainerType pNewElements = IsogeometricPostUtility::CreateEntities(pElements, r_other_model_part, r_clone_element, last_element_id, pProperties, retain_prop_id);

    for(ModelPart::ElementsContainerType::ptr_iterator it = pNewElements.ptr_begin(); it != pNewElements.ptr_end(); ++it)
        r_other_model_part.Elements().push_back(*it);

    std::cout << "Transfer mesh completed, "
              << pNewElements.size() << " elements of type " << sample_element_name
              << " was added to model_part" << r_other_model_part.Name() << std::endl;

    return pNewElements;
}

ModelPart::ConditionsContainerType IsogeometricPostUtility_TransferConditions(IsogeometricPostUtility& rDummy, ModelPart::ConditionsContainerType& pConditions,
    ModelPart& r_other_model_part, const std::string& sample_condition_name, Properties::Pointer pProperties, const bool& retain_prop_id)
{
    std::size_t last_condition_id = IsogeometricPostUtility::GetLastConditionId(r_other_model_part);
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_condition_name, "is not registered to the Kratos kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities(pConditions, r_other_model_part, r_clone_condition, last_condition_id, pProperties, retain_prop_id);

    for(ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin(); it != pNewConditions.ptr_end(); ++it)
        r_other_model_part.Conditions().push_back(*it);

    std::cout << "Transfer mesh completed, "
              << pNewConditions.size() << " conditions of type " << sample_condition_name
              << " was added to model_part" << r_other_model_part.Name() << std::endl;

    return pNewConditions;
}

ModelPart::ElementsContainerType IsogeometricPostUtility_FindElements(IsogeometricPostUtility& rDummy,
    ModelPart::ElementsContainerType& pElements, const std::string& sample_element_name)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_element_name, "is not registered to the Kratos kernel")
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    return IsogeometricPostUtility::FindEntities(pElements, r_clone_element);
}

ModelPart::ConditionsContainerType IsogeometricPostUtility_FindConditions(IsogeometricPostUtility& rDummy,
    ModelPart::ConditionsContainerType& pConditions, const std::string& sample_condition_name)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_condition_name, "is not registered to the Kratos kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    return IsogeometricPostUtility::FindEntities(pConditions, r_clone_condition);
}

template<typename TCoordinatesType, typename TPatchType>
boost::python::list IsogeometricPostUtility_CreateConditionsByTriangulation(IsogeometricPostUtility& rDummy,
    const boost::python::list& list_physical_points,
    const Vector& center, const Vector& normal, const Vector& t1, const Vector& t2,
    const boost::python::list& list_local_points, const std::size_t& nrefine,
    typename TPatchType::Pointer pPatch, ModelPart& r_model_part,
    const std::string& sample_condition_name,
    const std::size_t& last_node_id, const std::size_t& last_condition_id,
    Properties::Pointer pProperties)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_condition_name, "is not registered to the Kratos kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

    typedef boost::python::stl_input_iterator<array_1d<double, 3> > iterator_value_type;

    std::vector<TCoordinatesType> physical_points;
    BOOST_FOREACH(const iterator_value_type::value_type& p,
                std::make_pair(iterator_value_type(list_physical_points), // begin
                iterator_value_type() ) ) // end
    {
        physical_points.push_back(p);
    }

    std::vector<TCoordinatesType> local_points;
    BOOST_FOREACH(const iterator_value_type::value_type& p,
                std::make_pair(iterator_value_type(list_local_points), // begin
                iterator_value_type() ) ) // end
    {
        local_points.push_back(p);
    }

    std::size_t offset = last_node_id + 1;
    std::pair<std::vector<TCoordinatesType>, std::vector<std::vector<std::size_t> > >
    points_and_connectivities = IsogeometricPostUtility::GenerateTriangleGrid(physical_points, center, normal, t1, t2, local_points, offset, nrefine);

    boost::python::list new_nodes;
    boost::python::list new_local_points;
    std::size_t starting_node_id = last_node_id + 1;
    for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
    {
        new_local_points.append(points_and_connectivities.first[i]);
        ModelPart::NodeType::Pointer pNewNode = IsogeometricPostUtility::CreateNode(points_and_connectivities.first[i], *pPatch, r_model_part, starting_node_id++);
        new_nodes.append(pNewNode);
        // std::cout << "node " << pNewNode->Id() << " (" << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << ") is created at " << points_and_connectivities.first[i] << std::endl;
    }

    // std::cout << "connectivities:" << std::endl;
    // for (std::size_t i = 0; i < points_and_connectivities.second.size(); ++i)
    // {
    //     std::cout << " ";
    //     for (std::size_t j = 0; j < points_and_connectivities.second[i].size(); ++j)
    //         std::cout << " " << points_and_connectivities.second[i][j];
    //     std::cout << std::endl;
    // }

    std::size_t last_condition_id_new = last_condition_id;
    const std::string NodeKey = std::string("Node");
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<std::size_t> >, Condition, ModelPart::ConditionsContainerType>(
        points_and_connectivities.second, r_model_part, r_clone_condition, last_condition_id_new, pProperties, NodeKey);

    boost::python::list output;
    output.append(new_local_points);
    output.append(new_nodes);
    output.append(pNewConditions);
    return output;
}

template<typename TCoordinatesType, typename TPatchType>
boost::python::list IsogeometricPostUtility_CreateConditionsByQuadrilateralization(IsogeometricPostUtility& rDummy,
    const TCoordinatesType& p1, const TCoordinatesType& p2,
    const TCoordinatesType& p3, const TCoordinatesType& p4,
    const std::size_t& num_div_1, const std::size_t& num_div_2,
    typename TPatchType::Pointer pPatch, ModelPart& r_model_part,
    const std::string& sample_condition_name,
    const std::size_t& last_node_id, const std::size_t& last_condition_id,
    Properties::Pointer pProperties)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_condition_name, "is not registered to the Kratos kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

    std::size_t starting_node_id = last_node_id + 1;
    std::pair<std::vector<TCoordinatesType>, std::vector<std::vector<std::size_t> > >
    points_and_connectivities = IsogeometricPostUtility::GenerateQuadGrid(p1, p2, p3, p4, starting_node_id, num_div_1, num_div_2);

    boost::python::list new_nodes;
    boost::python::list new_local_points;
    for (std::size_t i = 0; i < points_and_connectivities.first.size(); ++i)
    {
        new_local_points.append(points_and_connectivities.first[i]);
        ModelPart::NodeType::Pointer pNewNode = IsogeometricPostUtility::CreateNode(points_and_connectivities.first[i], *pPatch, r_model_part, starting_node_id++);
        new_nodes.append(pNewNode);
        // std::cout << "node " << pNewNode->Id() << " (" << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << ") is created at " << points_and_connectivities.first[i] << std::endl;
    }

    // std::cout << "connectivities:" << std::endl;
    // for (std::size_t i = 0; i < points_and_connectivities.second.size(); ++i)
    // {
    //     std::cout << " ";
    //     for (std::size_t j = 0; j < points_and_connectivities.second[i].size(); ++j)
    //         std::cout << " " << points_and_connectivities.second[i][j];
    //     std::cout << std::endl;
    // }

    std::size_t last_condition_id_new = last_condition_id;
    const std::string NodeKey = std::string("Node");
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<std::size_t> >, Condition, ModelPart::ConditionsContainerType>(
        points_and_connectivities.second, r_model_part, r_clone_condition, last_condition_id_new, pProperties, NodeKey);

    boost::python::list output;
    output.append(new_local_points);
    output.append(new_nodes);
    output.append(pNewConditions);
    return output;
}

//////////////////////////////////////////////////////////

void BezierClassicalPostUtility_GenerateConditions(BezierClassicalPostUtility& dummy,
        ModelPart& rModelPart,
        Condition& rCondition,
        const std::string& sample_condition_name,
        const std::size_t& starting_node_id,
        const std::size_t& starting_condition_id)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_THROW_ERROR(std::logic_error, sample_condition_name, "is not registered to the Kratos kernel")
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    std::size_t NodeCounter = starting_node_id;
    std::size_t NodeCounter_old = NodeCounter;
    std::size_t ConditionCounter = starting_condition_id;
    const std::string NodeKey = std::string("Node");
    dummy.GenerateForOneEntity<Condition, ModelPart::ConditionsContainerType, 2>(rModelPart, rCondition,
            r_clone_condition, NodeCounter_old, NodeCounter, ConditionCounter, NodeKey);
}

void BezierClassicalPostUtility_GenerateModelPart2WithCondition(BezierClassicalPostUtility& dummy, ModelPart::Pointer pModelPartPost)
{
    dummy.GenerateModelPart2(pModelPartPost, true);
}

void BezierClassicalPostUtility_GenerateModelPart2(BezierClassicalPostUtility& dummy, ModelPart::Pointer pModelPartPost, const bool& generate_for_condition)
{
    dummy.GenerateModelPart2(pModelPartPost, generate_for_condition);
}

//////////////////////////////////////////////////////////

template<typename TVariableType>
void BezierPostUtility_TransferVariablesToNodes_ModelPart(BezierPostUtility& rDummy,
    const TVariableType& rThisVariable,
    ModelPart& r_model_part,
    BezierPostUtility::LinearSolverType::Pointer pSolver)
{
    rDummy.TransferVariablesToNodes(rThisVariable, r_model_part, pSolver);
}

template<typename TVariableType>
void BezierPostUtility_TransferVariablesToNodes_Elements(BezierPostUtility& rDummy,
    const TVariableType& rThisVariable,
    ModelPart& r_model_part, BezierPostUtility::ElementsArrayType& ElementsArray,
    BezierPostUtility::LinearSolverType::Pointer pSolver)
{
    rDummy.TransferVariablesToNodes(rThisVariable, r_model_part, ElementsArray, pSolver);
}

//////////////////////////////////////////////////////////

void IsogeometricApplication_AddBackendUtilitiesToPython()
{
    enum_<PostElementType>("PostElementType")
    .value("Triangle", _TRIANGLE_)
    .value("Quadrilateral", _QUADRILATERAL_)
    .value("Tetrahedra", _TETRAHEDRA_)
    .value("Hexahedra", _HEXAHEDRA_)
    ;

    class_<BSplineUtils, BSplineUtils::Pointer, boost::noncopyable>("BSplineUtils", init<>())
    .def("FindSpan", BSplineUtils_FindSpan)
    .def("BasisFuns", BSplineUtils_BasisFuns)
    .def("test_ComputeBsplinesKnotInsertionCoefficients1DLocal", &BSplineUtils::test_ComputeBsplinesKnotInsertionCoefficients1DLocal)
    ;

    class_<BezierUtils, BezierUtils::Pointer, boost::noncopyable>("BezierUtils", init<>())
    .def("Bernstein", BezierUtils_Bernstein)
    .def("BernsteinDerivative", BezierUtils_Bernstein_der)
    .def("DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients", BezierUtils_DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients)
    .def("ComputeCentroid", BezierUtils_ComputeCentroid<Element>)
    .def("ComputeCentroid", BezierUtils_ComputeCentroid<Condition>)
//    .def("compute_extended_knot_vector", &BezierUtils::compute_extended_knot_vector)
//    .def("bezier_extraction_tsplines_1d", &BezierUtils::bezier_extraction_tsplines_1d)
    ;

    class_<IsogeometricPostUtility,IsogeometricPostUtility::Pointer, boost::noncopyable>("IsogeometricPostUtility", init<>())
    .def("TransferElements", &IsogeometricPostUtility_TransferElements)
    .def("TransferConditions", &IsogeometricPostUtility_TransferConditions)
    .def("FindElements", &IsogeometricPostUtility_FindElements)
    .def("FindConditions", &IsogeometricPostUtility_FindConditions)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByTriangulation<array_1d<double, 3>, Patch<3> >)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByQuadrilateralization<Vector, Patch<3> >)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByQuadrilateralization<array_1d<double, 3>, Patch<3> >)
    ;

    class_<BezierClassicalPostUtility, BezierClassicalPostUtility::Pointer, boost::noncopyable>("BezierClassicalPostUtility", init<ModelPart::Pointer>())
    .def("GenerateConditions", &BezierClassicalPostUtility_GenerateConditions)
    .def("GenerateModelPart", &BezierClassicalPostUtility::GenerateModelPart)
    .def("GenerateModelPart2", &BezierClassicalPostUtility_GenerateModelPart2WithCondition)
    .def("GenerateModelPart2", &BezierClassicalPostUtility_GenerateModelPart2)
    .def("GenerateModelPart2AutoCollapse", &BezierClassicalPostUtility::GenerateModelPart2AutoCollapse)
    .def("TransferNodalResults", &BezierClassicalPostUtility::TransferNodalResults<Variable<double> >)
    .def("TransferNodalResults", &BezierClassicalPostUtility::TransferNodalResults<Variable<Vector> >)
    .def("TransferNodalResults", &BezierClassicalPostUtility::TransferNodalResults<Variable<array_1d<double, 3> > >)
    .def("TransferIntegrationPointResults", &BezierClassicalPostUtility::TransferIntegrationPointResults<Variable<double> >)
    .def("TransferIntegrationPointResults", &BezierClassicalPostUtility::TransferIntegrationPointResults<Variable<Vector> >)
    .def("SynchronizeActivation", &BezierClassicalPostUtility::SynchronizeActivation)
    .def("TransferElementalData", &BezierClassicalPostUtility::TransferElementalData<Variable<bool> >)
    .def("TransferConditionalData", &BezierClassicalPostUtility::TransferConditionalData<Variable<bool> >)
    .def("TransferVariablesToNodes", &BezierClassicalPostUtility::TransferVariablesToNodes<Variable<double> >)
    .def("TransferVariablesToNodes", &BezierClassicalPostUtility::TransferVariablesToNodes<Variable<Vector> >)
    .def("GlobalNodalRenumbering", &BezierClassicalPostUtility::GlobalNodalRenumbering)
    ;

    class_<BezierPostUtility, BezierPostUtility::Pointer, boost::noncopyable>("BezierPostUtility", init<>())
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<double> >)
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<Vector> >)
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<array_1d<double, 3> > >)
    .def("TransferIntegrationPointResults", &BezierPostUtility::TransferIntegrationPointResults<Variable<double> >)
    .def("TransferIntegrationPointResults", &BezierPostUtility::TransferIntegrationPointResults<Variable<Vector> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_ModelPart<Variable<double> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_ModelPart<Variable<Vector> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements<Variable<double> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements<Variable<Vector> >)
    ;

    #ifdef ISOGEOMETRIC_USE_HDF5
    class_<HDF5PostUtility, HDF5PostUtility::Pointer, boost::noncopyable>("HDF5PostUtility", init<const std::string>())
    .def(init<const std::string, const std::string>())
    .def("WriteNodes", &HDF5PostUtility::WriteNodes)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<double>)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<array_1d<double, 3> >)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<Vector>)
    .def("WriteElementalData", &HDF5PostUtility::WriteElementalData<bool>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<double>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<array_1d<double, 3> >)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<Vector>)
    .def("ReadElementalData", &HDF5PostUtility::ReadElementalData<bool>)
    ;
    #endif

    class_<NURBSTestUtils, NURBSTestUtils::Pointer, boost::noncopyable>("NURBSTestUtils", init<>())
    .def("Test1", &NURBSTestUtils::Test1)
    .def("Test2", &NURBSTestUtils::Test2)
    .def("ProbeGlobalCoordinates", &NURBSTestUtils_ProbeGlobalCoordinates2)
    .def("ProbeGlobalCoordinates", &NURBSTestUtils_ProbeGlobalCoordinates3)
    .def("ProbeShapeFunctionValues", &NURBSTestUtils_ProbeShapeFunctionValues2)
    .def("ProbeShapeFunctionValues", &NURBSTestUtils_ProbeShapeFunctionValues3)
    .def("ProbeShapeFunctionDerivatives", &NURBSTestUtils_ProbeShapeFunctionDerivatives2)
    .def("ProbeShapeFunctionDerivatives", &NURBSTestUtils_ProbeShapeFunctionDerivatives3)
    .def("ProbeJacobian", &NURBSTestUtils_ProbeJacobian2)
    .def("ProbeJacobian", &NURBSTestUtils_ProbeJacobian3)
    .def("DumpNodalValues", &NURBSTestUtils::DumpNodalValues<double>)
    .def("DumpNodalValues", &NURBSTestUtils::DumpNodalValues<array_1d<double, 3> >)
    ;

    class_<IsogeometricMergeUtility, IsogeometricMergeUtility::Pointer, boost::noncopyable>(
        "IsogeometricMergeUtility", init<>())
    .def("Add", &IsogeometricMergeUtility::Add)
    .def("Export", &IsogeometricMergeUtility::Export)
    .def("DumpNodalVariablesList", &IsogeometricMergeUtility::DumpNodalVariablesList)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////GISMO/////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    #ifdef ISOGEOMETRIC_USE_GISMO
    class_<GismoMesh, GismoMesh::Pointer, boost::noncopyable>
    ("GismoMesh", init<std::string>())
    .def("SetEchoLevel", &GismoMesh::SetEchoLevel)
    .def("ReadMesh", &GismoMesh::ReadMesh)
    .def(self_ns::str(self))
    ;
    #endif

}

}  // namespace Python.

} // Namespace Kratos

