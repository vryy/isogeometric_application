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
#include "includes/model_part.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/python_utils.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"
#include "custom_utilities/multipatch_utility.h"
#include "custom_utilities/multipatch_refinement_utility.h"
#include "custom_utilities/bending_strip_utility.h"
#include "custom_utilities/trim/isogeometric_intersection_utility.h"
#include "custom_python/add_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<class TPatchType>
typename TPatchType::Pointer MultiPatchUtility_CreatePatchPointer(MultiPatchUtility& rDummy, std::size_t Id, typename TPatchType::FESpaceType::Pointer pFESpace)
{
    return rDummy.template CreatePatchPointer<TPatchType>(Id, pFESpace);
}

template<class TPatchType>
typename TPatchType::Pointer MultiPatchUtility_CreateEmptyPatchPointer(MultiPatchUtility& rDummy, std::size_t Id)
{
    return rDummy.template CreateEmptyPatchPointer<TPatchType>(Id);
}

template<class TPatchType>
void MultiPatchUtility_MakeInterface(MultiPatchUtility& rDummy, typename TPatchType::Pointer pPatch1, const BoundarySide side1,
                                     typename TPatchType::Pointer pPatch2, const BoundarySide side2)
{
    rDummy.MakeInterface<TPatchType>(pPatch1, side1, pPatch2, side2);
}

template<class TMultiPatchType>
void MultiPatchUtility_CheckInterfaces(MultiPatchUtility& rDummy, typename TMultiPatchType::Pointer pMultiPatch)
{
    rDummy.CheckInterfaces<TMultiPatchType>(*pMultiPatch);
}

template<class TMultiPatchType>
void MultiPatchUtility_CheckInterfacesWithDebug(MultiPatchUtility& rDummy, typename TMultiPatchType::Pointer pMultiPatch, const bool debug)
{
    rDummy.CheckInterfaces<TMultiPatchType>(*pMultiPatch, debug);
}

template<class TMultiPatchType>
void MultiPatchUtility_CheckInterfacesFull(MultiPatchUtility& rDummy, typename TMultiPatchType::Pointer pMultiPatch, const bool debug, const double dist_tol)
{
    rDummy.CheckInterfaces<TMultiPatchType>(*pMultiPatch, debug, dist_tol);
}

template<class TMultiPatchType>
boost::python::list MultiPatchUtility_LocalCoordinates1(MultiPatchUtility& rDummy, typename TMultiPatchType::Pointer pMultiPatch,
        const boost::python::list& P, const boost::python::list& list_nsampling, const int echo_level)
{
    std::vector<double> P_vec;
    PythonUtils::Unpack<double, double>(P, P_vec);

    array_1d<double, 3> point, xi;
    noalias(point) = ZeroVector(3);
    noalias(xi) = ZeroVector(3);

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    std::vector<int> nsampling;
    PythonUtils::Unpack<int, int>(list_nsampling, nsampling);

    int patch_id = rDummy.LocalCoordinates(*pMultiPatch, point, xi, nsampling, echo_level);

    boost::python::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    boost::python::list output;
    output.append(patch_id);
    output.append(out_point);
    return output;
}

template<class TMultiPatchType>
boost::python::list MultiPatchUtility_LocalCoordinates(MultiPatchUtility& rDummy, typename TMultiPatchType::Pointer pMultiPatch,
        const boost::python::list& P, const boost::python::list& list_nsampling)
{
    return MultiPatchUtility_LocalCoordinates1<TMultiPatchType>(rDummy, pMultiPatch, P, list_nsampling, 0);
}

template<int TDim>
Matrix MultiPatchUtility_ComputeSpatialDerivatives(MultiPatchUtility& rDummy,
        const GridFunction<TDim, double, array_1d<double, 3> >& rControlPointGridFunction,
        const GridFunction<TDim, double, array_1d<double, 3> >& rControlValueGridFunction,
        const boost::python::list& xi_list)
{
    std::vector<double> xi;
    PythonUtils::Unpack<double, double>(xi_list, xi);

    return rDummy.ComputeSpatialDerivatives(rControlPointGridFunction,
                                            rControlValueGridFunction, xi);
}

std::size_t MultiPatchUtility_GetLastNodeId(MultiPatchUtility& rDummy, const ModelPart& r_model_part)
{
    return r_model_part.GetLastNodeId();
}

std::size_t MultiPatchUtility_GetLastElementId(MultiPatchUtility& rDummy, const ModelPart& r_model_part)
{
    return r_model_part.GetLastElementId();
}

std::size_t MultiPatchUtility_GetLastConditionId(MultiPatchUtility& rDummy, const ModelPart& r_model_part)
{
    return r_model_part.GetLastConditionId();
}

std::size_t MultiPatchUtility_GetLastConstraintId(MultiPatchUtility& rDummy, const ModelPart& r_model_part)
{
    return r_model_part.GetLastConstraintId();
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

std::size_t MultiPatchUtility_BoundaryFlag(MultiPatchUtility& rDummy, const BoundarySide side)
{
    return BOUNDARY_FLAG(side);
}

std::size_t MultiPatchUtility_BoundaryFlag2D(MultiPatchUtility& rDummy, const BoundarySide2D side)
{
    return BOUNDARY_FLAG(side);
}

std::size_t MultiPatchUtility_BoundaryFlag3D(MultiPatchUtility& rDummy, const BoundarySide3D side)
{
    return BOUNDARY_FLAG(side);
}

template<class TClassType>
void MultiPatchUtility_PrintAddress(MultiPatchUtility& rDummy, typename TClassType::Pointer pInstance)
{
    rDummy.PrintAddress<TClassType>(std::cout, pInstance);
}

//////////////////////////////////////////////////

template<class TPatchType>
void MultiPatchRefinementUtility_InsertKnots(typename TPatchType::Pointer& pPatch,
        const boost::python::list& ins_knots)
{
    std::vector<std::vector<typename TPatchType::LocalCoordinateType> > ins_knots_array(TPatchType::Dim);
    std::size_t dim = PythonUtils::Unpack<double, typename TPatchType::LocalCoordinateType>(ins_knots, ins_knots_array, TPatchType::Dim);

    if (dim != TPatchType::Dim)
        KRATOS_ERROR << "invalid dimension " << dim;

    MultiPatchRefinementUtility::InsertKnots<TPatchType>(pPatch, ins_knots_array);
}

template<class TPatchType>
boost::python::dict MultiPatchRefinementUtility_InsertKnots2(typename TPatchType::Pointer& pPatch,
        const boost::python::list& ins_knots)
{
    std::vector<std::vector<typename TPatchType::LocalCoordinateType> > ins_knots_array(TPatchType::Dim);
    std::size_t dim = PythonUtils::Unpack<double, typename TPatchType::LocalCoordinateType>(ins_knots, ins_knots_array, TPatchType::Dim);

    if (dim != TPatchType::Dim)
        KRATOS_ERROR << "Invalid dimension " << dim;

    std::map<std::size_t, Matrix> trans_mats;
    MultiPatchRefinementUtility::InsertKnots<TPatchType>(pPatch, ins_knots_array, trans_mats);
    // KRATOS_WATCH(trans_mats.size())

    boost::python::dict res;
    for (std::map<std::size_t, Matrix>::iterator it = trans_mats.begin(); it != trans_mats.end(); ++it)
    {
        res[it->first] = it->second;
    }

    return res;
}

template<class TPatchType>
void MultiPatchRefinementUtility_DegreeElevate(typename TPatchType::Pointer& pPatch,
        const boost::python::list& order_increment)
{
    std::vector<std::size_t> order_incr_array(TPatchType::Dim);
    std::size_t dim = PythonUtils::Unpack<int, std::size_t>(order_increment, order_incr_array, TPatchType::Dim);

    if (dim != TPatchType::Dim)
        KRATOS_ERROR << "Invalid dimension " << dim;

    MultiPatchRefinementUtility::DegreeElevate<TPatchType>(pPatch, order_incr_array);
}

//////////////////////////////////////////////////

template<class TPatchType>
typename TPatchType::Pointer BSplinesPatchUtility_CreateLoftPatch(BSplinesPatchUtility& dummy,
        typename TPatchType::BoundaryPatchType::Pointer pPatch1, typename TPatchType::BoundaryPatchType::Pointer pPatch2)
{
    return BSplinesPatchUtility::CreateLoftPatch<TPatchType>(pPatch1, pPatch2);
}

template<class TPatchType>
typename TPatchType::Pointer BSplinesPatchUtility_CreateLoftPatchFromList(BSplinesPatchUtility& dummy,
        const boost::python::list& patch_list, int order)
{
    typedef typename TPatchType::BoundaryPatchType::Pointer TBoundaryPatchPointerType;
    std::vector<TBoundaryPatchPointerType> pPatches;
    PythonUtils::Unpack<TBoundaryPatchPointerType, TBoundaryPatchPointerType>(patch_list, pPatches);

    return BSplinesPatchUtility::CreateLoftPatch<TPatchType>(pPatches, order);
}

boost::python::list BSplinesPatchUtility_CreatePatchFromGeo(BSplinesPatchUtility& dummy,
        const std::string& filename)
{
    const int Dim = BSplinesPatchUtility::GetDimensionOfGeo(filename);
    boost::python::list patches;
    if (Dim == 2)
    {
        patches.append(BSplinesPatchUtility::CreatePatchFromGeo<2>(filename));
    }
    else if (Dim == 3)
    {
        patches.append(BSplinesPatchUtility::CreatePatchFromGeo<3>(filename));
    }
    else
        KRATOS_ERROR << "Invalid dimension " << Dim;
    return patches;
}

template<typename TPatchPointerType>
void BSplinesPatchUtility_MakeInterface1D(BSplinesPatchUtility& rDummy,
        TPatchPointerType pPatch1, int iside1,
        TPatchPointerType pPatch2, int iside2)
{
    BoundarySide side1 = static_cast<BoundarySide>(iside1);
    BoundarySide side2 = static_cast<BoundarySide>(iside2);
    rDummy.MakeInterface1D(pPatch1, side1, pPatch2, side2);
}

template<typename TPatchPointerType>
void BSplinesPatchUtility_MakeInterface2D(BSplinesPatchUtility& rDummy,
        TPatchPointerType pPatch1, int iside1,
        TPatchPointerType pPatch2, int iside2,
        const BoundaryDirection direction)
{
    BoundarySide side1 = static_cast<BoundarySide>(iside1);
    BoundarySide side2 = static_cast<BoundarySide>(iside2);
    rDummy.MakeInterface2D(pPatch1, side1, pPatch2, side2, direction);
}

template<typename TPatchPointerType>
void BSplinesPatchUtility_MakeInterface3D(BSplinesPatchUtility& rDummy,
        TPatchPointerType pPatch1, int iside1,
        TPatchPointerType pPatch2, int iside2,
        const bool uv_or_vu,
        const BoundaryDirection direction1, const BoundaryDirection direction2)
{
    BoundarySide side1 = static_cast<BoundarySide>(iside1);
    BoundarySide side2 = static_cast<BoundarySide>(iside2);
    rDummy.MakeInterface3D(pPatch1, side1, pPatch2, side2, uv_or_vu, direction1, direction2);
}

template<typename TPatchPointerType>
void BSplinesPatchUtility_Reverse(BSplinesPatchUtility& rDummy,
                                  TPatchPointerType pPatch, std::size_t idir)
{
    rDummy.Reverse(pPatch, idir);
}

template<class TPatchType>
void BSplinesPatchUtility_Transpose2(BSplinesPatchUtility& rDummy,
                                     typename TPatchType::Pointer pPatch)
{
    rDummy.Transpose2D(pPatch);
}

template<class TPatchType>
void BSplinesPatchUtility_Transpose3(BSplinesPatchUtility& rDummy,
                                     typename TPatchType::Pointer pPatch, std::size_t idir, std::size_t jdir)
{
    rDummy.Transpose3D(pPatch, idir, jdir);
}

template<class TMultiPatchType>
void BSplinesPatchUtility_CheckRepeatedKnot(BSplinesPatchUtility& rDummy,
                                            typename TMultiPatchType::Pointer pMultiPatch)
{
    BSplinesPatchUtility::CheckRepeatedKnot<TMultiPatchType>(pMultiPatch);
}

//////////////////////////////////////////////////

template<int TDim>
typename Patch<TDim>::Pointer BendingStripUtility_CreateBendingStripNURBSPatch1(
    BendingStripUtility& rDummy,
    std::size_t Id,
    typename Patch<TDim>::Pointer pPatch1, const BoundarySide side1,
    typename Patch<TDim>::Pointer pPatch2, const BoundarySide side2,
    int Order)
{
    return rDummy.CreateBendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Order);
}

template<int TDim>
typename Patch<TDim>::Pointer BendingStripUtility_CreateBendingStripNURBSPatch2(
    BendingStripUtility& rDummy,
    std::size_t Id,
    typename Patch<TDim>::Pointer pPatch1, const BoundarySide side1,
    typename Patch<TDim>::Pointer pPatch2, const BoundarySide side2,
    const boost::python::list& order_list)
{
    std::vector<int> Orders(TDim);
    PythonUtils::Unpack<int, int>(order_list, Orders, TDim);

    return rDummy.CreateBendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Orders);
}

//////////////////////////////////////////////////

boost::python::list IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Two_Curves(IsogeometricIntersectionUtility& rDummy,
        double starting_point_1,
        double starting_point_2,
        Patch<1>::Pointer pPatch1,
        Patch<1>::Pointer pPatch2,
        int max_iters,
        double TOL,
        int option_space)
{
    double intersection_point_1, intersection_point_2;

    int stat = rDummy.ComputeIntersectionByNewtonRaphson(starting_point_1, starting_point_2,
               intersection_point_1, intersection_point_2,
               pPatch1, pPatch2, max_iters, TOL, option_space);

    boost::python::list point;
    point.append(intersection_point_1);
    point.append(intersection_point_2);

    boost::python::list output;
    output.append(stat);
    output.append(point);
    return output;
}

boost::python::list IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Curve_Plane(IsogeometricIntersectionUtility& rDummy,
        double starting_point,
        Patch<1>::Pointer pPatch,
        double A, double B, double C, double D,
        int max_iters,
        double TOL)
{
    // std::cout << "invoking " << __FUNCTION__ << std::endl;

    double intersection_point = starting_point;

    int stat = rDummy.ComputeIntersectionByNewtonRaphson(intersection_point, pPatch, A, B, C, D, max_iters, TOL);

    boost::python::list point;
    point.append(intersection_point);

    boost::python::list output;
    output.append(stat);
    output.append(point);
    return output;
}

boost::python::list IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Patch2_Plane(IsogeometricIntersectionUtility& rDummy,
        const boost::python::list& list_starting_points,
        Patch<2>::Pointer pPatch,
        double A, double B, double C, double D,
        int max_iters,
        double TOL)
{
    // std::cout << "invoking " << __FUNCTION__ << std::endl;

    std::vector<double> starting_points;
    PythonUtils::Unpack<double, double>(list_starting_points, starting_points);

    std::vector<std::vector<double> > intersection_points;

    std::vector<int> stat = rDummy.ComputeIntersectionByNewtonRaphson(starting_points, intersection_points,
                            pPatch, A, B, C, D, max_iters, TOL);

    boost::python::list list_points;
    for (std::size_t i = 0; i < intersection_points.size(); ++i)
    {
        boost::python::list point;
        point.append(intersection_points[i][0]);
        point.append(intersection_points[i][1]);
        list_points.append(point);
    }

    boost::python::list list_stat;
    for (std::size_t i = 0; i < stat.size(); ++i)
    {
        list_stat.append(stat[i]);
    }

    boost::python::list output;
    output.append(list_stat);
    output.append(list_points);
    return output;
}

boost::python::list IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Patch3_Plane(IsogeometricIntersectionUtility& rDummy,
        const boost::python::list& list_starting_points,
        Patch<3>::Pointer pPatch,
        double A, double B, double C, double D,
        int max_iters,
        double TOL)
{
    // std::cout << "invoking " << __FUNCTION__ << std::endl;

    std::vector<double> starting_points;
    PythonUtils::Unpack<double, double>(list_starting_points, starting_points);

    std::vector<std::vector<double> > intersection_points;

    std::vector<int> stat = rDummy.ComputeIntersectionByNewtonRaphson(starting_points, intersection_points,
                            pPatch, A, B, C, D, max_iters, TOL);

    boost::python::list list_points;
    for (std::size_t i = 0; i < intersection_points.size(); ++i)
    {
        boost::python::list point;
        point.append(intersection_points[i][0]);
        point.append(intersection_points[i][1]);
        point.append(intersection_points[i][2]);
        list_points.append(point);
    }

    boost::python::list list_stat;
    for (std::size_t i = 0; i < stat.size(); ++i)
    {
        list_stat.append(stat[i]);
    }

    boost::python::list output;
    output.append(list_stat);
    output.append(list_points);
    return output;
}

boost::python::list IsogeometricIntersectionUtility_ComputeIntersectionByBisection_Patch3_Plane(IsogeometricIntersectionUtility& rDummy,
        Patch<3>::Pointer pPatch,
        double A, double B, double C, double D,
        int max_iters,
        double TOL)
{
    // std::cout << "invoking " << __FUNCTION__ << std::endl;

    std::vector<array_1d<double, 3> > intersection_points;

    std::vector<int> stat = rDummy.ComputeIntersectionByBisection(intersection_points, pPatch, A, B, C, D, max_iters, TOL);

    boost::python::list list_points;
    for (std::size_t i = 0; i < intersection_points.size(); ++i)
    {
        list_points.append(intersection_points[i]);
    }

    boost::python::list list_stat;
    for (std::size_t i = 0; i < stat.size(); ++i)
    {
        list_stat.append(stat[i]);
    }

    boost::python::list output;
    output.append(list_stat);
    output.append(list_points);
    return output;
}

boost::python::list IsogeometricIntersectionUtility_ComputeIntersection_Curve_Surface(IsogeometricIntersectionUtility& rDummy,
        double starting_point_1,
        double starting_point_2_1,
        double starting_point_2_2,
        Patch<1>::Pointer pPatch1,
        Patch<2>::Pointer pPatch2,
        int max_iters,
        double TOL)
{
    // std::cout << "invoking " << __FUNCTION__ << std::endl;

    double intersection_point_1;
    std::vector<double> starting_point_2(2), intersection_point_2(2);

    starting_point_2[0] = starting_point_2_1;
    starting_point_2[1] = starting_point_2_2;

    int stat = rDummy.ComputeIntersectionByNewtonRaphson(starting_point_1, starting_point_2,
               intersection_point_1, intersection_point_2,
               pPatch1, pPatch2, max_iters, TOL);

    boost::python::list point;
    point.append(intersection_point_2[0]);
    point.append(intersection_point_2[1]);

    boost::python::list output;
    output.append(stat);
    output.append(intersection_point_1);
    output.append(point);
    return output;
}

template<int TDim>
boost::python::list IsogeometricIntersectionUtility_CheckIntersection(IsogeometricIntersectionUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch,
        double A, double B, double C, double D)
{
    std::pair<int, std::vector<int> > result = rDummy.CheckIntersection<TDim, 0>(pPatch, A, B, C, D);
    boost::python::list output;
    output.append(result.first);
    boost::python::list tmp;
    for (std::size_t i = 0; i < result.second.size(); ++i)
    {
        tmp.append(result.second[i]);
    }
    output.append(tmp);
    return output;
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////

void IsogeometricApplication_AddFrontendUtilitiesToPython()
{

    class_<MultiPatchUtility, MultiPatchUtility::Pointer, boost::noncopyable>
    ("MultiPatchUtility", init<>())
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<1>::RealPatch>)
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<2>::RealPatch>)
    .def("CreatePatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<3>::RealPatch>)
    .def("CreateEmptyPatch1DPointer", &MultiPatchUtility_CreateEmptyPatchPointer<PatchSelector<1>::RealPatch>)
    .def("CreateEmptyPatch2DPointer", &MultiPatchUtility_CreateEmptyPatchPointer<PatchSelector<2>::RealPatch>)
    .def("CreateEmptyPatch3DPointer", &MultiPatchUtility_CreateEmptyPatchPointer<PatchSelector<3>::RealPatch>)
    .def("CreateComplexPatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<1>::ComplexPatch>)
    .def("CreateComplexPatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<2>::ComplexPatch>)
    .def("CreateComplexPatchPointer", &MultiPatchUtility_CreatePatchPointer<PatchSelector<3>::ComplexPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<1>::RealPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<2>::RealPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<3>::RealPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<1>::ComplexPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<2>::ComplexPatch>)
    .def("MakeInterface", &MultiPatchUtility_MakeInterface<PatchSelector<3>::ComplexPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<1>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<2>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<3>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<1>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<2>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfaces<PatchSelector<3>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<1>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<2>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<3>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<1>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<2>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesWithDebug<PatchSelector<3>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<1>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<2>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<3>::RealMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<1>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<2>::ComplexMultiPatch>)
    .def("CheckInterfaces", &MultiPatchUtility_CheckInterfacesFull<PatchSelector<3>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<1>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<2>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<3>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<1>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<2>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates<PatchSelector<3>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<1>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<2>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<3>::RealMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<1>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<2>::ComplexMultiPatch>)
    .def("LocalCoordinates", &MultiPatchUtility_LocalCoordinates1<PatchSelector<3>::ComplexMultiPatch>)
    .def("ComputeSpatialDerivatives", &MultiPatchUtility_ComputeSpatialDerivatives<2>)
    .def("ComputeSpatialDerivatives", &MultiPatchUtility_ComputeSpatialDerivatives<3>)
    .def("GetLastNodeId", &MultiPatchUtility_GetLastNodeId)
    .def("GetLastElementId", &MultiPatchUtility_GetLastElementId)
    .def("GetLastConditionId", &MultiPatchUtility_GetLastConditionId)
    .def("GetLastConstraintId", &MultiPatchUtility_GetLastConstraintId)
    .def("CreateConditionFromElement", &MultiPatchUtility_CreateConditionFromElement)
    .def("ListModelPart", &MultiPatchUtility_ListModelPart)
    .def("GetEquationId", &MultiPatchUtility_GetEquationId<Variable<double> >)
    .def("GetEquationId", &MultiPatchUtility_GetEquationId<Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > >)
    .def("BoundaryFlag", &MultiPatchUtility_BoundaryFlag)
    .def("BoundaryFlag", &MultiPatchUtility_BoundaryFlag2D)
    .def("BoundaryFlag", &MultiPatchUtility_BoundaryFlag3D)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<PatchSelector<1>::RealPatch>)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<PatchSelector<2>::RealPatch>)
    .def("PrintAddress", &MultiPatchUtility_PrintAddress<PatchSelector<3>::RealPatch>)
    ;

    class_<MultiPatchRefinementUtility, MultiPatchRefinementUtility::Pointer, boost::noncopyable>
    ("MultiPatchRefinementUtility", init<>())
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<1>::RealPatch>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<2>::RealPatch>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<3>::RealPatch>)
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<1>::RealPatch>) // deprecated
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<2>::RealPatch>) // deprecated
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<3>::RealPatch>) // deprecated
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<1>::RealPatch>)
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<2>::RealPatch>)
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<3>::RealPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<1>::RealPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<2>::RealPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<3>::RealPatch>)
    // same as above but for ComplexPatch
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<1>::ComplexPatch>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<2>::ComplexPatch>)
    .def("InsertKnots", MultiPatchRefinementUtility_InsertKnots<PatchSelector<3>::ComplexPatch>)
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<1>::ComplexPatch>) // deprecated
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<2>::ComplexPatch>) // deprecated
    .def("InsertKnots2", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<3>::ComplexPatch>) // deprecated
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<1>::ComplexPatch>)
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<2>::ComplexPatch>)
    .def("InsertKnotsGetTrans", MultiPatchRefinementUtility_InsertKnots2<PatchSelector<3>::ComplexPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<1>::ComplexPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<2>::ComplexPatch>)
    .def("DegreeElevate", MultiPatchRefinementUtility_DegreeElevate<PatchSelector<3>::ComplexPatch>)
    // declare the static method. Same name is allowed for polymophism methods.
    .staticmethod("InsertKnots")
    .staticmethod("InsertKnots2")
    .staticmethod("InsertKnotsGetTrans")
    .staticmethod("DegreeElevate")
    ;

    class_<BSplinesPatchUtility, BSplinesPatchUtility::Pointer, boost::noncopyable>
    ("BSplinesPatchUtility", init<>())
    .def("CreatePatchFromGeo", &BSplinesPatchUtility_CreatePatchFromGeo)
    //
    .def("CreateLoftPatch", &BSplinesPatchUtility_CreateLoftPatch<PatchSelector<2>::RealPatch>)
    .def("CreateLoftPatch", &BSplinesPatchUtility_CreateLoftPatch<PatchSelector<3>::RealPatch>)
    .def("CreateLoftPatchFromList2D", &BSplinesPatchUtility_CreateLoftPatchFromList<PatchSelector<2>::RealPatch>)
    .def("CreateLoftPatchFromList3D", &BSplinesPatchUtility_CreateLoftPatchFromList<PatchSelector<3>::RealPatch>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface1D<PatchSelector<1>::RealPatch::Pointer>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface2D<PatchSelector<2>::RealPatch::Pointer>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface3D<PatchSelector<3>::RealPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<1>::RealPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<2>::RealPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<3>::RealPatch::Pointer>)
    .def("Transpose", &BSplinesPatchUtility_Transpose2<PatchSelector<2>::RealPatch>)
    .def("Transpose", &BSplinesPatchUtility_Transpose3<PatchSelector<3>::RealPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<1>::RealMultiPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<2>::RealMultiPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<3>::RealMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<1>::RealMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<2>::RealMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<3>::RealMultiPatch>)
    // same as above but for ComplexPatch
    .def("CreateLoftPatch", &BSplinesPatchUtility_CreateLoftPatch<PatchSelector<2>::ComplexPatch>)
    .def("CreateLoftPatch", &BSplinesPatchUtility_CreateLoftPatch<PatchSelector<3>::ComplexPatch>)
    .def("CreateComplexLoftPatchFromList2D", &BSplinesPatchUtility_CreateLoftPatchFromList<PatchSelector<2>::ComplexPatch>)
    .def("CreateComplexLoftPatchFromList3D", &BSplinesPatchUtility_CreateLoftPatchFromList<PatchSelector<3>::ComplexPatch>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface1D<PatchSelector<1>::ComplexPatch::Pointer>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface2D<PatchSelector<2>::ComplexPatch::Pointer>)
    .def("MakeInterface", &BSplinesPatchUtility_MakeInterface3D<PatchSelector<3>::ComplexPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<1>::ComplexPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<2>::ComplexPatch::Pointer>)
    .def("Reverse", &BSplinesPatchUtility_Reverse<PatchSelector<3>::ComplexPatch::Pointer>)
    .def("Transpose", &BSplinesPatchUtility_Transpose2<PatchSelector<2>::ComplexPatch>)
    .def("Transpose", &BSplinesPatchUtility_Transpose3<PatchSelector<3>::ComplexPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<1>::ComplexMultiPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<2>::ComplexMultiPatch>)
    .def("CreateInterfaces", &BSplinesPatchUtility::CreateInterfaces<PatchSelector<3>::ComplexMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<1>::ComplexMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<2>::ComplexMultiPatch>)
    .def("CheckRepeatedKnot", &BSplinesPatchUtility_CheckRepeatedKnot<PatchSelector<3>::ComplexMultiPatch>)
    ;

    class_<BendingStripUtility, BendingStripUtility::Pointer, boost::noncopyable>
    ("BendingStripUtility", init<>())
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch1<2>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch1<3>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch2<2>)
    .def("CreateBendingStripNURBSPatch", &BendingStripUtility_CreateBendingStripNURBSPatch2<3>)
    ;

    class_<IsogeometricIntersectionUtility, IsogeometricIntersectionUtility::Pointer, boost::noncopyable>
    ("IsogeometricIntersectionUtility", init<>())
    .def("ComputeIntersectionByNewtonRaphson", &IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Two_Curves)
    .def("ComputeIntersectionByNewtonRaphson", &IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Curve_Plane)
    .def("ComputeIntersectionByNewtonRaphson", &IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Patch2_Plane)
    .def("ComputeIntersectionByNewtonRaphson", &IsogeometricIntersectionUtility_ComputeIntersectionByNewtonRaphson_Patch3_Plane)
    .def("ComputeIntersectionByBisection", &IsogeometricIntersectionUtility_ComputeIntersectionByBisection_Patch3_Plane)
    .def("CheckIntersection", &IsogeometricIntersectionUtility_CheckIntersection<1>)
    .def("CheckIntersection", &IsogeometricIntersectionUtility_CheckIntersection<2>)
    .def("CheckIntersection", &IsogeometricIntersectionUtility_CheckIntersection<3>)
    ;

}

}  // namespace Python.

} // Namespace Kratos
