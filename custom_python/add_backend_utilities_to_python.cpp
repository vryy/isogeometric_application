/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 9, 2013 $
//
//

// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "python/python_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/isogeometric_utility.h"
#include "custom_utilities/isogeometric_post_utility.h"
#include "custom_utilities/bezier_classical_post_utility.h"
#include "custom_utilities/bezier_post_utility.h"
#include "custom_utilities/isogeometric_test_utils.h"
#include "custom_utilities/bezier_test_utils.h"
#include "custom_utilities/isogeometric_merge_utility.h"
#include "custom_utilities/isogeometric_math_utils.h"
#include "custom_utilities/isogeometric_projection_utility.h"
#include "custom_python/add_utilities_to_python.h"

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
    ModelPart& r_model_part,
    std::string FileName
)
{
    dummy.DumpShapeFunctionsIntegrationPointsValuesAndLocalGradients(r_model_part, FileName);
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

template<typename TEntityType>
void IsogeometricTestUtils_ProbeGlobalCoordinates1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeGlobalCoordinates(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeGlobalCoordinates2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeGlobalCoordinates(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeGlobalCoordinates3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeGlobalCoordinates(pElement->GetGeometry(), X, Y, Z);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionValues1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeShapeFunctionValues(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionValues2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeShapeFunctionValues(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionValues3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeShapeFunctionValues(pElement->GetGeometry(), X, Y, Z);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionDerivatives1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionDerivatives2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionDerivatives3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeShapeFunctionDerivatives(pElement->GetGeometry(), X, Y, Z);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeJacobian1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeJacobian(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeJacobian2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeJacobian(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeJacobian3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeJacobian(pElement->GetGeometry(), X, Y, Z);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeShapeFunctionSecondDerivatives(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeShapeFunctionSecondDerivatives(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeShapeFunctionSecondDerivatives(pElement->GetGeometry(), X, Y, Z);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives1(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X
)
{
    dummy.ProbeShapeFunctionThirdDerivatives(pElement->GetGeometry(), X, 0.0, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives2(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y
)
{
    dummy.ProbeShapeFunctionThirdDerivatives(pElement->GetGeometry(), X, Y, 0.0);
}

template<typename TEntityType>
void IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives3(
    IsogeometricTestUtils& dummy,
    typename TEntityType::Pointer pElement,
    double X, double Y, double Z
)
{
    dummy.ProbeShapeFunctionThirdDerivatives(pElement->GetGeometry(), X, Y, Z);
}

//////////////////////////////////////////////////////////

ModelPart::ElementsContainerType IsogeometricPostUtility_TransferElements(IsogeometricPostUtility& rDummy, ModelPart::ElementsContainerType& pElements,
        ModelPart& r_other_model_part, const std::string& sample_element_name, Properties::Pointer pProperties, const bool& retain_prop_id)
{
    std::size_t last_element_id = r_other_model_part.GetLastElementId();
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_ERROR << sample_element_name << " is not registered to the Kratos kernel";
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    ModelPart::ElementsContainerType pNewElements = IsogeometricPostUtility::CreateEntities(pElements, r_other_model_part, r_clone_element, last_element_id, pProperties, retain_prop_id);

    for (ModelPart::ElementsContainerType::ptr_iterator it = pNewElements.ptr_begin(); it != pNewElements.ptr_end(); ++it)
    {
        r_other_model_part.Elements().push_back(*it);
    }

    std::cout << "Transfer mesh completed, "
              << pNewElements.size() << " elements of type " << sample_element_name
              << " was added to model_part" << r_other_model_part.Name() << std::endl;

    return pNewElements;
}

ModelPart::ConditionsContainerType IsogeometricPostUtility_TransferConditions(IsogeometricPostUtility& rDummy, ModelPart::ConditionsContainerType& pConditions,
        ModelPart& r_other_model_part, const std::string& sample_condition_name, Properties::Pointer pProperties, const bool& retain_prop_id)
{
    std::size_t last_condition_id = r_other_model_part.GetLastConditionId();
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities(pConditions, r_other_model_part, r_clone_condition, last_condition_id, pProperties, retain_prop_id);

    for (ModelPart::ConditionsContainerType::ptr_iterator it = pNewConditions.ptr_begin(); it != pNewConditions.ptr_end(); ++it)
    {
        r_other_model_part.Conditions().push_back(*it);
    }

    std::cout << "Transfer mesh completed, "
              << pNewConditions.size() << " conditions of type " << sample_condition_name
              << " was added to model_part" << r_other_model_part.Name() << std::endl;

    return pNewConditions;
}

ModelPart::ElementsContainerType IsogeometricPostUtility_FindElements(IsogeometricPostUtility& rDummy,
        ModelPart::ElementsContainerType& pElements, const std::string& sample_element_name)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_ERROR << sample_element_name << " is not registered to the Kratos kernel";
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    return IsogeometricPostUtility::FindEntities(pElements, r_clone_element);
}

ModelPart::ConditionsContainerType IsogeometricPostUtility_FindConditions(IsogeometricPostUtility& rDummy,
        ModelPart::ConditionsContainerType& pConditions, const std::string& sample_condition_name)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    return IsogeometricPostUtility::FindEntities(pConditions, r_clone_condition);
}

template<typename TCoordinatesType, typename TPatchType>
boost::python::list IsogeometricPostUtility_CreateConditionsByTriangulation(IsogeometricPostUtility& rDummy,
        const boost::python::list& list_physical_points,
        const Vector& center, const Vector& normal, const Vector& t1, const Vector& t2,
        const boost::python::list& list_local_points, std::size_t nrefine,
        typename TPatchType::Pointer pPatch, ModelPart& r_model_part,
        const std::string& sample_condition_name,
        std::size_t last_node_id, std::size_t last_condition_id,
        Properties::Pointer pProperties)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);

    std::vector<TCoordinatesType> physical_points;
    PythonUtils::Unpack<array_1d<double, 3>, TCoordinatesType>(list_physical_points, physical_points);

    std::vector<TCoordinatesType> local_points;
    PythonUtils::Unpack<array_1d<double, 3>, TCoordinatesType>(list_local_points, local_points);

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
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities<ModelPart, std::vector<std::vector<std::size_t> >, Condition, ModelPart::ConditionsContainerType>(
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
        std::size_t num_div_1, std::size_t num_div_2,
        typename TPatchType::Pointer pPatch, ModelPart& r_model_part,
        const std::string& sample_condition_name,
        std::size_t last_node_id, std::size_t last_condition_id,
        Properties::Pointer pProperties)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
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
    ModelPart::ConditionsContainerType pNewConditions = IsogeometricPostUtility::CreateEntities<ModelPart, std::vector<std::vector<std::size_t> >, Condition, ModelPart::ConditionsContainerType>(
                points_and_connectivities.second, r_model_part, r_clone_condition, last_condition_id_new, pProperties, NodeKey);

    boost::python::list output;
    output.append(new_local_points);
    output.append(new_nodes);
    output.append(pNewConditions);
    return output;
}

template<class TPatchType>
void IsogeometricPostUtility_TransferValuesToNodes(IsogeometricPostUtility& rDummy, Element::GeometryType::PointType& rNode, const TPatchType& rPatch)
{
    rDummy.TransferValuesToNodes(rNode, rPatch);
}

template<class TEntityType, typename TVariableType, class TPatchType>
void IsogeometricPostUtility_TransferValuesToGaussPoints(IsogeometricPostUtility& rDummy, TEntityType& rElement,
        const TVariableType& rVariable, const TPatchType& rPatch, const ProcessInfo& rProcessInfo)
{
    rDummy.TransferValuesToGaussPoints(rElement, rVariable, rPatch, rProcessInfo);
}

//////////////////////////////////////////////////////////

boost::python::list BezierClassicalPostUtility_GenerateConditions(BezierClassicalPostUtility& dummy,
        ModelPart& rModelPart,
        const Condition& rCondition,
        const std::string& sample_condition_name,
        std::size_t starting_node_id,
        std::size_t starting_condition_id)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    std::size_t NodeCounter = starting_node_id;
    std::size_t NodeCounter_old = NodeCounter;
    std::size_t ConditionCounter = starting_condition_id;
    const std::string NodeKey = std::string("Node");
    std::vector<std::size_t> node_ids;
    std::vector<std::size_t> element_ids;
    dummy.GenerateForOneEntity<Condition, ModelPart::ConditionsContainerType, 2>(rModelPart, rCondition,
            r_clone_condition, NodeCounter_old, NodeCounter, ConditionCounter, NodeKey, false,
            node_ids, element_ids, true);

    boost::python::list list_nodes;
    boost::python::list list_elements;
    for (std::size_t i = 0; i < node_ids.size(); ++i) { list_nodes.append(node_ids[i]); }
    for (std::size_t i = 0; i < element_ids.size(); ++i) { list_elements.append(element_ids[i]); }
    boost::python::list Output;
    Output.append(list_nodes);
    Output.append(list_elements);
    return Output;
}

boost::python::list BezierClassicalPostUtility_GenerateConditionsWithNodalVariables(BezierClassicalPostUtility& dummy,
        ModelPart& rModelPart,
        const Condition& rCondition,
        const std::string& sample_condition_name,
        std::size_t starting_node_id,
        std::size_t starting_condition_id)
{
    if (!KratosComponents<Condition>::Has(sample_condition_name))
        KRATOS_ERROR << sample_condition_name << " is not registered to the Kratos kernel";
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(sample_condition_name);
    std::size_t NodeCounter = starting_node_id;
    std::size_t NodeCounter_old = NodeCounter;
    std::size_t ConditionCounter = starting_condition_id;
    const std::string NodeKey = std::string("Node");
    std::vector<std::size_t> node_ids;
    std::vector<std::size_t> element_ids;
    dummy.GenerateForOneEntity<Condition, ModelPart::ConditionsContainerType, 2>(rModelPart, rCondition,
            r_clone_condition, NodeCounter_old, NodeCounter, ConditionCounter, NodeKey, true,
            node_ids, element_ids, true);

    boost::python::list list_nodes;
    boost::python::list list_elements;
    for (std::size_t i = 0; i < node_ids.size(); ++i) { list_nodes.append(node_ids[i]); }
    for (std::size_t i = 0; i < element_ids.size(); ++i) { list_elements.append(element_ids[i]); }
    boost::python::list Output;
    Output.append(list_nodes);
    Output.append(list_elements);
    return Output;
}

boost::python::list BezierClassicalPostUtility_GenerateElements(BezierClassicalPostUtility& dummy,
        ModelPart& rModelPart,
        const Element& rElement,
        const std::string& sample_element_name,
        std::size_t starting_node_id,
        std::size_t starting_element_id)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_ERROR << sample_element_name << " is not registered to the Kratos kernel";
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    std::size_t NodeCounter = starting_node_id;
    std::size_t NodeCounter_old = NodeCounter;
    std::size_t ElementCounter = starting_element_id;
    const std::string NodeKey = std::string("Node");
    std::vector<std::size_t> node_ids;
    std::vector<std::size_t> element_ids;
    dummy.GenerateForOneEntity<Element, ModelPart::ElementsContainerType, 2>(rModelPart, rElement,
            r_clone_element, NodeCounter_old, NodeCounter, ElementCounter, NodeKey, false,
            node_ids, element_ids, true);

    boost::python::list list_nodes;
    boost::python::list list_elements;
    for (std::size_t i = 0; i < node_ids.size(); ++i) { list_nodes.append(node_ids[i]); }
    for (std::size_t i = 0; i < element_ids.size(); ++i) { list_elements.append(element_ids[i]); }
    boost::python::list Output;
    Output.append(list_nodes);
    Output.append(list_elements);
    return Output;
}

boost::python::list BezierClassicalPostUtility_GenerateElementsWithNodalVariables(BezierClassicalPostUtility& dummy,
        ModelPart& rModelPart,
        const Element& rElement,
        const std::string& sample_element_name,
        std::size_t starting_node_id,
        std::size_t starting_element_id)
{
    if (!KratosComponents<Element>::Has(sample_element_name))
        KRATOS_ERROR << sample_element_name << " is not registered to the Kratos kernel";
    Element const& r_clone_element = KratosComponents<Element>::Get(sample_element_name);
    std::size_t NodeCounter = starting_node_id;
    std::size_t NodeCounter_old = NodeCounter;
    std::size_t ElementCounter = starting_element_id;
    const std::string NodeKey = std::string("Node");
    std::vector<std::size_t> node_ids;
    std::vector<std::size_t> element_ids;
    dummy.GenerateForOneEntity<Element, ModelPart::ElementsContainerType, 2>(rModelPart, rElement,
            r_clone_element, NodeCounter_old, NodeCounter, ElementCounter, NodeKey, true,
            node_ids, element_ids, true);

    boost::python::list list_nodes;
    boost::python::list list_elements;
    for (std::size_t i = 0; i < node_ids.size(); ++i) { list_nodes.append(node_ids[i]); }
    for (std::size_t i = 0; i < element_ids.size(); ++i) { list_elements.append(element_ids[i]); }
    boost::python::list Output;
    Output.append(list_nodes);
    Output.append(list_elements);
    return Output;
}

void BezierClassicalPostUtility_GenerateModelPart2WithCondition(BezierClassicalPostUtility& dummy, ModelPart& rModelPartPost)
{
    dummy.GenerateModelPart2(rModelPartPost, true);
}

void BezierClassicalPostUtility_GenerateModelPart2(BezierClassicalPostUtility& dummy, ModelPart& rModelPartPost, bool generate_for_condition)
{
    dummy.GenerateModelPart2(rModelPartPost, generate_for_condition);
}

//////////////////////////////////////////////////////////

bool IsogeometricUtility_IsNormalPointingOutward(IsogeometricUtility& rDummy, const Patch<3>& rPatch, std::size_t iside)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return IsogeometricUtility::IsNormalPointingOutward(rPatch, side);
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
        ModelPart& r_model_part, const BezierPostUtility::ElementsContainerType& ElementsArray,
        BezierPostUtility::LinearSolverType::Pointer pSolver)
{
    rDummy.TransferVariablesToNodes(rThisVariable, r_model_part, ElementsArray, pSolver);
}

void BezierPostUtility_TransferVariablesToNodes_Elements_Vector(BezierPostUtility& rDummy,
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part, const BezierPostUtility::ElementsContainerType& ElementsArray,
        BezierPostUtility::LinearSolverType::Pointer pSolver, std::size_t ncomponents)
{
    rDummy.TransferVariablesToNodes(rThisVariable, r_model_part, ElementsArray, pSolver, ncomponents);
}

boost::python::dict BezierPostUtility_TransferVariablesToNodalArray_Elements_Double(BezierPostUtility& rDummy,
        const Variable<double>& rThisVariable,
        const ModelPart& r_model_part, const BezierPostUtility::ElementsContainerType& ElementsArray,
        BezierPostUtility::LinearSolverType::Pointer pSolver)
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    BezierPostUtility::SerialSparseSpaceType::VectorType g;

    rDummy.TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part, ElementsArray, rThisVariable);

    boost::python::dict Output;

    for (auto it = active_nodes.begin(); it != active_nodes.end(); ++it)
    {
        std::size_t row = node_row_id[*it];
        Output[*it] = g[row];
    }

    return Output;
}

boost::python::dict BezierPostUtility_TransferVariablesToNodalArray_Elements_Array1D(BezierPostUtility& rDummy,
        const Variable<array_1d<double, 3> >& rThisVariable,
        const ModelPart& r_model_part, const BezierPostUtility::ElementsContainerType& ElementsArray,
        BezierPostUtility::LinearSolverType::Pointer pSolver)
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    BezierPostUtility::SerialDenseSpaceType::MatrixType g;

    rDummy.TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part,
                                         ElementsArray, rThisVariable);

    boost::python::dict Output;

    array_1d<double, 3> aux;
    for (auto it = active_nodes.begin(); it != active_nodes.end(); ++it)
    {
        std::size_t this_row = node_row_id[*it];
        noalias(aux) = row(g, this_row);
        Output[*it] = aux;
    }

    return Output;
}

boost::python::dict BezierPostUtility_TransferVariablesToNodalArray_Elements_Vector(BezierPostUtility& rDummy,
        const Variable<Vector>& rThisVariable,
        const ModelPart& r_model_part, const BezierPostUtility::ElementsContainerType& ElementsArray,
        BezierPostUtility::LinearSolverType::Pointer pSolver, std::size_t ncomponents)
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    BezierPostUtility::SerialDenseSpaceType::MatrixType g;

    rDummy.TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part, ElementsArray, rThisVariable, ncomponents);

    boost::python::dict Output;

    Vector aux(ncomponents);
    for (auto it = active_nodes.begin(); it != active_nodes.end(); ++it)
    {
        std::size_t this_row = node_row_id[*it];
        noalias(aux) = row(g, this_row);
        Output[*it] = aux;
    }

    return Output;
}

//////////////////////////////////////////////////////////

template<typename TPointType>
boost::python::list IsogeometricMathUtils_ComputeIntersection2D(IsogeometricMathUtils<double>& rDummy,
        const TPointType& P0, const TPointType& P1,
        const TPointType& P2, const TPointType& P3,
        const double TOL)
{
    TPointType P;
    int error_code = IsogeometricMathUtils<double>::ComputeIntersection2D(P0, P1, P2, P3, P, TOL);

    boost::python::list output;
    output.append(error_code);
    output.append(P);
    return std::move(output);
}

template<typename TPointType>
boost::python::list IsogeometricMathUtils_ComputeProjection2D(IsogeometricMathUtils<double>& rDummy,
        const TPointType& P0,
        const TPointType& P1, const TPointType& P2,
        const double TOL)
{
    TPointType P;
    int error_code = IsogeometricMathUtils<double>::ComputeProjection2D(P0, P1, P2, P, TOL);

    boost::python::list output;
    output.append(error_code);
    output.append(P);
    return std::move(output);
}

//////////////////////////////////////////////////////////

/// Python wrapper for IsogeometricProjectionUtility
struct IsogeometricProjectionUtilityPythonInterface
{
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricProjectionUtilityPythonInterface);
};

template<typename TPointType, int TDim>
boost::python::list IsogeometricProjectionUtility_ComputeRayProjection(IsogeometricProjectionUtilityPythonInterface& rDummy,
        const TPointType& rPoint, const TPointType& rDirection, typename MultiPatch<TDim>::Pointer pMultiPatch,
        double TOL, int max_iters, const boost::python::list& list_nsamplings, const int echo_level)
{
    std::vector<double> LocalPoint(TDim);
    TPointType GlobalPoint;
    int patch_id;
    std::array<unsigned int, TDim> nsamplings;
    PythonUtils::Unpack<int, unsigned int, TDim>(list_nsamplings, nsamplings);
    int error_code = IsogeometricProjectionUtility<TPointType, TDim>::ComputeRayProjection(rPoint, rDirection, LocalPoint, GlobalPoint, patch_id,
                                        pMultiPatch, TOL, max_iters, nsamplings, echo_level);

    boost::python::list output;
    output.append(error_code);
    output.append(GlobalPoint);

    return std::move(output);
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

    class_<IsogeometricEcho, boost::noncopyable>("IsogeometricEcho", init<>())
    .def("SetEchoLevel", &IsogeometricEcho::SetEchoLevel)
    // .def("GetEchoLevel", IsogeometricEcho_GetEchoLevel)
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

    class_<IsogeometricPostUtility, IsogeometricPostUtility::Pointer, boost::noncopyable>("IsogeometricPostUtility", init<>())
    .def("TransferElements", &IsogeometricPostUtility_TransferElements)
    .def("TransferConditions", &IsogeometricPostUtility_TransferConditions)
    .def("TransferValuesToNodes", &IsogeometricPostUtility_TransferValuesToNodes<Patch<2> >)
    .def("TransferValuesToNodes", &IsogeometricPostUtility_TransferValuesToNodes<Patch<2> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<double>, Patch<2> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<array_1d<double, 3> >, Patch<2> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<Vector>, Patch<2> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<double>, Patch<3> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<array_1d<double, 3> >, Patch<3> >)
    .def("TransferValuesToGaussPoints", &IsogeometricPostUtility_TransferValuesToGaussPoints<Element, Variable<Vector>, Patch<3> >)
    .def("FindElements", &IsogeometricPostUtility_FindElements)
    .def("FindConditions", &IsogeometricPostUtility_FindConditions)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByTriangulation<array_1d<double, 3>, Patch<3> >)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByQuadrilateralization<Vector, Patch<3> >)
    .def("CreateConditions", &IsogeometricPostUtility_CreateConditionsByQuadrilateralization<array_1d<double, 3>, Patch<3> >)
    ;

    class_<BezierClassicalPostUtility, BezierClassicalPostUtility::Pointer, boost::noncopyable>("BezierClassicalPostUtility", init<ModelPart&>())
    .def("GenerateConditions", &BezierClassicalPostUtility_GenerateConditions)
    .def("GenerateConditionsWithNodalVariables", &BezierClassicalPostUtility_GenerateConditionsWithNodalVariables)
    .def("GenerateElements", &BezierClassicalPostUtility_GenerateElements)
    .def("GenerateElementsWithNodalVariables", &BezierClassicalPostUtility_GenerateElementsWithNodalVariables)
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
    .def(init<const bool>())
    .def("SetCheckActive", &BezierPostUtility::SetCheckActive)
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<double> >)
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<Vector> >)
    .def("TransferNodalResults", &BezierPostUtility::TransferNodalResults<Variable<array_1d<double, 3> > >)
    .def("TransferIntegrationPointResults", &BezierPostUtility::TransferIntegrationPointResults<Variable<double> >)
    .def("TransferIntegrationPointResults", &BezierPostUtility::TransferIntegrationPointResults<Variable<Vector> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_ModelPart<Variable<double> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_ModelPart<Variable<Vector> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_ModelPart<Variable<array_1d<double, 3> > >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements<Variable<double> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements<Variable<Vector> >)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements_Vector)
    .def("TransferVariablesToNodes", &BezierPostUtility_TransferVariablesToNodes_Elements<Variable<array_1d<double, 3> > >)
    .def("TransferVariablesToNodalArray", &BezierPostUtility_TransferVariablesToNodalArray_Elements_Double)
    .def("TransferVariablesToNodalArray", &BezierPostUtility_TransferVariablesToNodalArray_Elements_Vector)
    .def("TransferVariablesToNodalArray", &BezierPostUtility_TransferVariablesToNodalArray_Elements_Array1D)
    ;

    class_<IsogeometricTestUtils, IsogeometricTestUtils::Pointer, boost::noncopyable>("IsogeometricTestUtils", init<>())
    .def("Test1", &IsogeometricTestUtils::Test1)
    .def("Test2", &IsogeometricTestUtils::Test2)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates1<Element>)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates1<Condition>)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates2<Element>)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates2<Condition>)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates3<Element>)
    .def("ProbeGlobalCoordinates", &IsogeometricTestUtils_ProbeGlobalCoordinates3<Condition>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues1<Element>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues1<Condition>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues2<Element>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues2<Condition>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues3<Element>)
    .def("ProbeShapeFunctionValues", &IsogeometricTestUtils_ProbeShapeFunctionValues3<Condition>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives1<Element>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives1<Condition>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives2<Element>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives2<Condition>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives3<Element>)
    .def("ProbeShapeFunctionDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionDerivatives3<Condition>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives1<Element>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives1<Condition>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives2<Element>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives2<Condition>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives3<Element>)
    .def("ProbeShapeFunctionSecondDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionSecondDerivatives3<Condition>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives1<Element>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives1<Condition>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives2<Element>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives2<Condition>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives3<Element>)
    .def("ProbeShapeFunctionThirdDerivatives", &IsogeometricTestUtils_ProbeShapeFunctionThirdDerivatives3<Condition>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian1<Element>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian1<Condition>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian2<Element>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian2<Condition>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian3<Element>)
    .def("ProbeJacobian", &IsogeometricTestUtils_ProbeJacobian3<Condition>)
    .def("DumpNodalValues", &IsogeometricTestUtils::DumpNodalValues<double>)
    .def("DumpNodalValues", &IsogeometricTestUtils::DumpNodalValues<array_1d<double, 3> >)
    .def("DumpNodalValues", &IsogeometricTestUtils::DumpNodalValues<Vector>)
    ;

    class_<IsogeometricMergeUtility, IsogeometricMergeUtility::Pointer, boost::noncopyable>(
        "IsogeometricMergeUtility", init<>())
    .def("Add", &IsogeometricMergeUtility::Add)
    .def("Export", &IsogeometricMergeUtility::Export)
    .def("DumpNodalVariablesList", &IsogeometricMergeUtility::DumpNodalVariablesList)
    ;

    class_<IsogeometricUtility, IsogeometricUtility::Pointer, boost::noncopyable>(
        "IsogeometricUtility", init<>())
    // .def("IsNormalPointingOutward", (bool (IsogeometricUtility::*)(const Patch<3>&, const BoundarySide&)&IsogeometricUtility::IsNormalPointingOutward)).staticmethod("IsNormalPointingOutward")
    .def("IsNormalPointingOutward", &IsogeometricUtility_IsNormalPointingOutward);

    class_<IsogeometricMathUtils<double>, IsogeometricMathUtils<double>::Pointer, boost::noncopyable>(
        "IsogeometricMathUtils", init<>())
    .def("ComputeIntersection2D", &IsogeometricMathUtils_ComputeIntersection2D<array_1d<double, 3> >)
    .def("ComputeProjection2D", &IsogeometricMathUtils_ComputeProjection2D<array_1d<double, 3> >)
    ;

    class_<IsogeometricProjectionUtilityPythonInterface, IsogeometricProjectionUtilityPythonInterface::Pointer, boost::noncopyable>(
        "IsogeometricProjectionUtility", init<>())
    .def("ComputeRayProjection", &IsogeometricProjectionUtility_ComputeRayProjection<array_1d<double, 3>, 1>)
    .def("ComputeRayProjection", &IsogeometricProjectionUtility_ComputeRayProjection<Element::GeometryType::PointType::PointType, 1>)
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
