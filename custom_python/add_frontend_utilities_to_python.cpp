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
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "python/pointer_vector_set_python_interface.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/multipatch_refinement_utility.h"
#include "custom_utilities/bending_strip_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim>
typename Patch<TDim>::Pointer MultiPatchUtility_CreatePatchPointer(MultiPatchUtility& rDummy, const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
{
    return rDummy.template CreatePatchPointer<TDim>(Id, pFESpace);
}

template<int TDim>
void MultiPatchUtility_MakeInterface(MultiPatchUtility& rDummy, typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
    typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2)
{
    rDummy.MakeInterface<TDim>(pPatch1, side1, pPatch2, side2);
}

std::size_t MultiPatchUtility_GetLastNodeId(MultiPatchUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t MultiPatchUtility_GetLastElementId(MultiPatchUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t MultiPatchUtility_GetLastConditionId(MultiPatchUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

Condition::Pointer MultiPatchUtility_CreateConditionFromElement(MultiPatchUtility& rDummy,
        const std::string& sample_condition_name, std::size_t lastConditionId,
        Element::Pointer pElement, Properties::Pointer pProperties)
{
    return rDummy.CreateConditionFromElement(sample_condition_name, lastConditionId, pElement, pProperties);
}

void MultiPatchUtility_ListModelPart(MultiPatchUtility& rDummy, ModelPart& r_model_part)
{
    rDummy.ListModelPart(r_model_part);
}

template<typename TVariableType>
std::size_t MultiPatchUtility_GetEquationId(MultiPatchUtility& rDummy, ModelPart::NodeType& rNode, const TVariableType& rVariable)
{
    return rDummy.GetEquationId(rNode, rVariable);
}

template<class TClassType>
void MultiPatchUtility_PrintAddress(MultiPatchUtility& rDummy, typename TClassType::Pointer pInstance)
{
    rDummy.PrintAddress<TClassType>(std::cout, pInstance);
}

//////////////////////////////////////////////////

template<int TDim>
void MultiPatchRefinementUtility_InsertKnots(MultiPatchRefinementUtility& rDummy,
       typename Patch<TDim>::Pointer& pPatch,
       boost::python::list ins_knots)
{
    std::vector<std::vector<double> > ins_knots_array(TDim);
    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& ins_knots_x,
                std::make_pair(iterator_value_type(ins_knots), // begin
                iterator_value_type() ) ) // end
    {
        std::vector<double> knots;

        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        BOOST_FOREACH(const iterator_value_type2::value_type& knot,
                    std::make_pair(iterator_value_type2(ins_knots_x), // begin
                    iterator_value_type2() ) ) // end
        {
            knots.push_back(knot);
        }

        ins_knots_array[dim++] = knots;
        if (dim == TDim)
            break;
    }

    if (dim != TDim)
        KRATOS_THROW_ERROR(std::logic_error, "insufficient dimension", "")

    rDummy.InsertKnots<TDim>(pPatch, ins_knots_array);
}

template<int TDim>
boost::python::dict MultiPatchRefinementUtility_InsertKnots2(MultiPatchRefinementUtility& rDummy,
       typename Patch<TDim>::Pointer& pPatch,
       boost::python::list ins_knots)
{
    std::vector<std::vector<double> > ins_knots_array(TDim);
    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& ins_knots_x,
                std::make_pair(iterator_value_type(ins_knots), // begin
                iterator_value_type() ) ) // end
    {
        std::vector<double> knots;

        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        BOOST_FOREACH(const iterator_value_type2::value_type& knot,
                    std::make_pair(iterator_value_type2(ins_knots_x), // begin
                    iterator_value_type2() ) ) // end
        {
            knots.push_back(knot);
        }

        ins_knots_array[dim++] = knots;
        if (dim == TDim)
            break;
    }

    if (dim != TDim)
        KRATOS_THROW_ERROR(std::logic_error, "insufficient dimension", "")

    std::map<std::size_t, Matrix> trans_mats;
    rDummy.InsertKnots<TDim>(pPatch, ins_knots_array, trans_mats);
    // KRATOS_WATCH(trans_mats.size())

    boost::python::dict res;
    for (std::map<std::size_t, Matrix>::iterator it = trans_mats.begin(); it != trans_mats.end(); ++it)
        res[it->first] = it->second;

    return res;
}

template<int TDim>
void MultiPatchRefinementUtility_DegreeElevate(MultiPatchRefinementUtility& rDummy,
       typename Patch<TDim>::Pointer& pPatch,
       boost::python::list order_increment)
{
    std::vector<std::size_t> order_incr_array(TDim);
    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& t,
               std::make_pair(iterator_value_type(order_increment), // begin
               iterator_value_type() ) ) // end
    {
       order_incr_array[dim++] = static_cast<std::size_t>(t);
       if (dim == TDim)
            break;
    }

    if (dim != TDim)
        KRATOS_THROW_ERROR(std::logic_error, "insufficient dimension", "")

    rDummy.DegreeElevate<TDim>(pPatch, order_incr_array);
}

//////////////////////////////////////////////////

template<int TDim>
typename Patch<TDim>::Pointer BSplinesPatchUtility_CreateConnectedPatch(BSplinesPatchUtility& dummy,
        typename Patch<TDim-1>::Pointer pPatch1, typename Patch<TDim-1>::Pointer pPatch2)
{
    return BSplinesPatchUtility::CreateConnectedPatch<TDim>(pPatch1, pPatch2);
}

template<int TDim>
typename Patch<TDim>::Pointer BSplinesPatchUtility_CreateConnectedPatch2(BSplinesPatchUtility& dummy,
        boost::python::list patch_list, const int& order)
{
    std::vector<typename Patch<TDim-1>::Pointer> pPatches;

    typedef boost::python::stl_input_iterator<typename Patch<TDim-1>::Pointer> iterator_value_type;
    BOOST_FOREACH(const typename iterator_value_type::value_type& pPatch,
                std::make_pair(iterator_value_type(patch_list), // begin
                iterator_value_type() ) ) // end
    {
        pPatches.push_back(pPatch);
    }

    return BSplinesPatchUtility::CreateConnectedPatch<TDim>(pPatches, order);
}

boost::python::list BSplinesPatchUtility_CreatePatchFromGeo(BSplinesPatchUtility& dummy,
        const std::string& filename)
{
    int Dim = BSplinesPatchUtility::GetDimensionOfGeo(filename);
    boost::python::list patches;
    if (Dim == 2)
        patches.append(BSplinesPatchUtility::CreatePatchFromGeo<2>(filename));
    else if (Dim == 3)
        patches.append(BSplinesPatchUtility::CreatePatchFromGeo<3>(filename));
    else
        KRATOS_THROW_ERROR(std::logic_error, "The dimension of the patch is invalid", "")
    return patches;
}

void BSplinesPatchUtility_MakeInterface2D(BSplinesPatchUtility& rDummy,
    typename Patch<2>::Pointer pPatch1, const int& iside1,
    typename Patch<2>::Pointer pPatch2, const int& iside2,
    const BoundaryDirection& direction)
{
    BoundarySide side1 = static_cast<BoundarySide>(iside1);
    BoundarySide side2 = static_cast<BoundarySide>(iside2);
    rDummy.MakeInterface2D(pPatch1, side1, pPatch2, side2, direction);
}

void BSplinesPatchUtility_MakeInterface3D(BSplinesPatchUtility& rDummy,
    typename Patch<3>::Pointer pPatch1, const int& iside1,
    typename Patch<3>::Pointer pPatch2, const int& iside2,
    const bool& uv_or_vu,
    const BoundaryDirection& direction1, const BoundaryDirection& direction2)
{
    BoundarySide side1 = static_cast<BoundarySide>(iside1);
    BoundarySide side2 = static_cast<BoundarySide>(iside2);
    rDummy.MakeInterface3D(pPatch1, side1, pPatch2, side2, uv_or_vu, direction1, direction2);
}

//////////////////////////////////////////////////

template<int TDim>
typename Patch<TDim>::Pointer BendingStripUtility_CreateBendingStripNURBSPatch1(
        BendingStripUtility& rDummy,
        const std::size_t& Id,
        typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
        typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2,
        const int& Order)
{
    return rDummy.CreateBendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Order);
}

template<int TDim>
typename Patch<TDim>::Pointer BendingStripUtility_CreateBendingStripNURBSPatch2(
        BendingStripUtility& rDummy,
        const std::size_t& Id,
        typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
        typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2,
        boost::python::list order_list)
{
    std::vector<int> Orders(TDim);

    std::size_t dim = 0;

    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& t,
                std::make_pair(iterator_value_type(order_list), // begin
                iterator_value_type() ) ) // end
    {
        Orders[dim++] = static_cast<int>(t);
        if (dim == TDim)
            break;
    }

    return rDummy.CreateBendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Orders);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void IsogeometricApplication_AddFrontendUtilitiesToPython()
{

    class_<MultiPatchUtility, MultiPatchUtility::Pointer, boost::noncopyable>
    ("MultiPatchUtility", init<>())
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<1>)
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<2>)
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<3>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<1>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<2>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<3>)
    .def("GetLastNodeId", &MultiPatchUtility_GetLastNodeId)
    .def("GetLastElementId", &MultiPatchUtility_GetLastElementId)
    .def("GetLastConditionId", &MultiPatchUtility_GetLastConditionId)
    .def("CreateConditionFromElement", &MultiPatchUtility_CreateConditionFromElement)
    .def("ListModelPart", &MultiPatchUtility_ListModelPart)
    .def("GetEquationId", &MultiPatchUtility_GetEquationId<Variable<double> >)
    .def("GetEquationId", &MultiPatchUtility_GetEquationId<Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<Patch<1> >)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<Patch<2> >)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<Patch<3> >)
    ;

    class_<MultiPatchRefinementUtility, MultiPatchRefinementUtility::Pointer, boost::noncopyable>
    ("MultiPatchRefinementUtility", init<>())
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<1>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<2>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<3>)
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<1>)
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<2>)
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<3>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<1>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<2>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<3>)
    ;

    class_<BSplinesPatchUtility, BSplinesPatchUtility::Pointer, boost::noncopyable>
    ("BSplinesPatchUtility", init<>())
    .def("CreateConnectedPatch", &BSplinesPatchUtility_CreateConnectedPatch<2>)
    .def("CreateConnectedPatch", &BSplinesPatchUtility_CreateConnectedPatch<3>)
    .def("CreateConnectedPatch", &BSplinesPatchUtility_CreateConnectedPatch2<2>)
    .def("CreateConnectedPatch", &BSplinesPatchUtility_CreateConnectedPatch2<3>)
    .def("CreatePatchFromGeo", &BSplinesPatchUtility_CreatePatchFromGeo)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface2D)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface3D)
    ;

    class_<BendingStripUtility, BendingStripUtility::Pointer, boost::noncopyable>
    ("BendingStripUtility", init<>())
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch1<2>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch1<3>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch2<2>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch2<3>)
    ;

}

}  // namespace Python.

} // Namespace Kratos

