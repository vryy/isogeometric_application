/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 24, 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "python/containers_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/hbsplines/deprecated_hb_mesh.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_utilities/hbsplines/hbsplines_patch_utility.h"
#include "custom_utilities/hbsplines/hbsplines_refinement_utility.h"
#include "custom_utilities/import_export/multi_hbsplines_patch_matlab_exporter.h"
#include "custom_python/iga_define_python.h"
#include "custom_python3/add_hbsplines_to_python.h"
#include "custom_python3/add_point_based_control_grid_to_python.h"
#include "custom_python3/add_import_export_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

////////////////////////////////////////

// template<int TDim>
// std::size_t HBSplinesBasisFunction_GetId(HBSplinesBasisFunction<TDim>& rDummy)
// {
//     return rDummy.Id();
// }

template<int TDim>
void HBSplinesBasisFunction_SetId(HBSplinesBasisFunction<TDim>& rDummy, std::size_t Id)
{
    // DO NOTHING
}

template<int TDim>
std::size_t HBSplinesBasisFunction_GetLevel(HBSplinesBasisFunction<TDim>& rDummy)
{
    return rDummy.Level();
}

// template<int TDim>
// std::size_t HBSplinesBasisFunction_GetEquationId(HBSplinesBasisFunction<TDim>& rDummy)
// {
//     return rDummy.EquationId();
// }

// template<int TDim>
// void HBSplinesBasisFunction_SetEquationId(HBSplinesBasisFunction<TDim>& rDummy, std::size_t EquationId)
// {
//     rDummy.SetEquationId(EquationId);
// }

template<int TDim>
pybind11::list HBSplinesFESpace_ExtractBoundaryBfsByFlag(HBSplinesFESpace<TDim>& rDummy, std::size_t boundary_id)
{
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;

    std::vector<bf_t> bf_list = rDummy.ExtractBoundaryBfsByFlag(boundary_id);

    pybind11::list Output;
    for (std::size_t i = 0; i < bf_list.size(); ++i)
    {
        Output.append(bf_list[i]);
    }

    return Output;
}

template<int TDim>
std::size_t HBSplinesFESpace_MaxLevel(HBSplinesFESpace<TDim>& rDummy)
{
    return rDummy.MaxLevel();
}

template<int TDim>
typename HBSplinesBasisFunction<TDim>::Pointer HBSplinesFESpace_GetItem(HBSplinesFESpace<TDim>& rDummy, std::size_t i)
{
    return rDummy[i];
}

////////////////////////////////////////

template<int TDim>
void HBSplinesPatchUtility_ListBoundaryBfs(HBSplinesPatchUtility& rDummy,
        typename HBSplinesFESpace<TDim>::Pointer pFESpace, BoundarySide side)
{
    rDummy.ListBoundaryBfs<TDim>(std::cout, pFESpace, side);
}

template<int TDim>
typename Patch<TDim>::Pointer HBSplinesPatchUtility_CreatePatchFromBSplines(HBSplinesPatchUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch)
{
    return HBSplinesPatchUtility::CreatePatchFromBSplines<TDim>(pPatch);
}

template<int TDim>
typename HBSplinesFESpace<TDim>::bf_t HBSplinesPatchUtility_GetBfByEquationId(HBSplinesPatchUtility& rDummy,
        typename MultiPatch<TDim>::Pointer pMultiPatch, std::size_t EquationId)
{
    return rDummy.GetBfByEquationId<TDim>(pMultiPatch, EquationId);
}

template<int TDim>
void HBSplinesPatchUtility_ReportDuplicatedEquationId(HBSplinesPatchUtility& rDummy,
        typename MultiPatch<TDim>::Pointer pMultiPatch, const bool& throw_error)
{
    return rDummy.ReportDuplicatedEquationId<TDim>(pMultiPatch, throw_error);
}

////////////////////////////////////////

template<int TDim>
void HBSplinesRefinementUtility_Refine(HBSplinesRefinementUtility& rDummy,
                                       typename Patch<TDim>::Pointer pPatch, std::size_t Id, int EchoLevel)
{
    rDummy.Refine<TDim>(pPatch, Id, EchoLevel);
}

template<int TDim>
void HBSplinesRefinementUtility_RefineBf(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, typename HBSplinesFESpace<TDim>::bf_t p_bf, int EchoLevel)
{
    rDummy.Refine<TDim>(pPatch, p_bf, EchoLevel);
}

template<int TDim>
void HBSplinesRefinementUtility_RefineWindow(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, pybind11::list& window, int EchoLevel)
{
    std::vector<std::vector<double> > window_vector;
    std::size_t cnt1 = 0, cnt2 = 0;
    for (auto vect : window)
    {
        std::vector<double> win_vect;
        for (auto v : vect.cast<pybind11::list>())
        {
            win_vect.push_back(v.cast<double>());
        }
        window_vector.push_back(win_vect);
    }
    rDummy.RefineWindow<TDim>(pPatch, window_vector, EchoLevel);
}

template<int TDim>
void HBSplinesRefinementUtility_LinearDependencyRefine(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, std::size_t refine_cycle, int EchoLevel)
{
    rDummy.LinearDependencyRefine<TDim>(pPatch, refine_cycle, EchoLevel);
}

////////////////////////////////////////

// template<typename TDataType, class TFESpaceType>
// typename ControlGrid<TDataType>::Pointer ControlGridUtility_CreatePointBasedControlGrid(
//         ControlGridUtility& rDummy,
//         const Variable<TDataType>& rVariable, typename TFESpaceType::Pointer pFESpace)
// {
//     return rDummy.CreatePointBasedControlGrid<TDataType, TFESpaceType>(rVariable, pFESpace);
// }

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddHBSplinesSpaceToPython(pybind11::module& m)
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "HBSplinesBasisFunction" << TDim << "D";
    class_<HBSplinesBasisFunction<TDim>, typename HBSplinesBasisFunction<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<std::size_t, std::size_t>())
    // .add_property("Id", HBSplinesBasisFunction_GetId<TDim>, HBSplinesBasisFunction_SetId<TDim>)
    .def_property("Id", Isogeometric_GetId<HBSplinesBasisFunction<TDim> >, HBSplinesBasisFunction_SetId<TDim>)
    // .add_property("EquationId", HBSplinesBasisFunction_GetEquationId<TDim>, HBSplinesBasisFunction_SetEquationId<TDim>)
    .def_property("EquationId", Isogeometric_GetEquationId<HBSplinesBasisFunction<TDim>>, Isogeometric_SetEquationId<HBSplinesBasisFunction<TDim>>)
    .def_property_readonly("Weight", &HBSplinesBasisFunction<TDim>::Weight)
    .def_property_readonly("Level", &HBSplinesBasisFunction_GetLevel<TDim>)
    .def("__str__", &PrintObject<HBSplinesBasisFunction<TDim> >)
    ;

    ss.str(std::string());
    ss << "HBSplinesFESpace" << TDim << "D";
//    typename FESpace<TDim-1>::Pointer(HBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace1)(const BoundarySide& side) const = &HBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    // typename FESpace<TDim-1>::Pointer(HBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace2)(const BoundarySide& side, const BoundaryRotation& rotation) const = &HBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    class_<HBSplinesFESpace<TDim>, typename HBSplinesFESpace<TDim>::Pointer, FESpace<TDim> >
    (m, ss.str().c_str())
    .def(init<>())
    .def("__getitem__", &HBSplinesFESpace_GetItem<TDim>)
    .def("GetBoundaryBfs", &HBSplinesFESpace_ExtractBoundaryBfsByFlag<TDim>) // deprecated
    .def("ExtractBoundaryBfsByFlag", &HBSplinesFESpace_ExtractBoundaryBfsByFlag<TDim>)
//    .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace1)
    // .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace2)
    .def("UpdateCells", &HBSplinesFESpace<TDim>::UpdateCells)
    .def("MaxLevel", &HBSplinesFESpace_MaxLevel<TDim>)
    .def("SetMaxLevel", &HBSplinesFESpace<TDim>::SetMaxLevel)
    .def("GetBfByEquationId", &HBSplinesFESpace<TDim>::pGetBfByEquationId)
    .def("HasBfByEquationId", &HBSplinesFESpace<TDim>::HasBfByEquationId)
    .def("HasBfById", &HBSplinesFESpace<TDim>::HasBfById)
    .def("__str__", &PrintObject<HBSplinesFESpace<TDim> >)
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, HBSplinesFESpace<TDim> >::Execute(m);
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, HBSplinesFESpace<TDim> >::Execute(m);
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, HBSplinesFESpace<TDim> >::Execute(m);

    ////////////////////OLD H-SPLINES////////////////////////

    ss.str(std::string());
    ss << "DeprecatedHBMesh" << TDim << "D";
    // class_<DeprecatedHBMesh<TDim>, Patch<TDim> > // does not compile
    class_<DeprecatedHBMesh<TDim>, typename DeprecatedHBMesh<TDim>::Pointer, Patch<TDim> >
    (m, ss.str().c_str())
    .def(init<std::size_t, const std::string&>())
    .def("SetEchoLevel", &DeprecatedHBMesh<TDim>::SetEchoLevel)
    .def("ReadMesh", &DeprecatedHBMesh<TDim>::ReadMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("Refine", &DeprecatedHBMesh<TDim>::Refine) // use this for debugging only, use RefineNodes and LinearDependencyRefine instead
    .def("RefineNodes", &DeprecatedHBMesh<TDim>::RefineNodes)
    .def("LinearDependencyRefine", &DeprecatedHBMesh<TDim>::LinearDependencyRefine)
    .def("BuildMesh", &DeprecatedHBMesh<TDim>::BuildMesh)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportCellTopology", &DeprecatedHBMesh<TDim>::ExportCellTopology)
    .def("ExportCellGeology", &DeprecatedHBMesh<TDim>::ExportCellGeology)
    //    .def("ExportRefinedDomain", &DeprecatedHBMesh<TDim>::ExportRefinedDomain)
    .def("ExportSupportDomain", &DeprecatedHBMesh<TDim>::ExportSupportDomain)
    .def("ExportMatlab", &DeprecatedHBMesh<TDim>::ExportMatlab)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("ExportMDPA", &DeprecatedHBMesh<TDim>::ExportMDPA)
    .def("ExportMDPA2", &DeprecatedHBMesh<TDim>::ExportMDPA2)
    .def("ExportPostMDPA", &DeprecatedHBMesh<TDim>::ExportPostMDPA)
    .def("ExportCellGeologyAsPostMDPA", &DeprecatedHBMesh<TDim>::ExportCellGeologyAsPostMDPA)
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    .def("PrintKnotVectors", &DeprecatedHBMesh<TDim>::PrintKnotVectors)
    .def("PrintCells", &DeprecatedHBMesh<TDim>::PrintCells)
    .def("PrintBasisFuncs", &DeprecatedHBMesh<TDim>::PrintBasisFuncs)
    .def("PrintRefinementHistory", &DeprecatedHBMesh<TDim>::PrintRefinementHistory)
    .def("CheckNestedSpace", &DeprecatedHBMesh<TDim>::CheckNestedSpace)
    .def("__str__", &PrintObject<DeprecatedHBMesh<TDim> >)
    ;

    ss.str(std::string());
    ss << "DeprecatedHBMesh" << TDim << "DPointer";
    class_<typename DeprecatedHBMesh<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<typename DeprecatedHBMesh<TDim>::Pointer>())
    .def("GetReference", GetReference<DeprecatedHBMesh<TDim> >, return_value_policy::reference_internal)
    .def("__str__", &PrintObject<typename DeprecatedHBMesh<TDim>::Pointer>)
    ;

}

////////////////////////////////////////

void IsogeometricApplication_AddHBSplinesToPython(pybind11::module& m)
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////HIERARCHICAL BSplines/////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddHBSplinesSpaceToPython<1>(m);
    IsogeometricApplication_AddHBSplinesSpaceToPython<2>(m);
    IsogeometricApplication_AddHBSplinesSpaceToPython<3>(m);

    class_<HBSplinesPatchUtility, HBSplinesPatchUtility::Pointer>
    (m, "HBSplinesPatchUtility")
    .def(init<>())
    .def("CreatePatchFromBSplines", &HBSplinesPatchUtility_CreatePatchFromBSplines<2>)
    .def("CreatePatchFromBSplines", &HBSplinesPatchUtility_CreatePatchFromBSplines<3>)
    .def("ListBoundaryBfs", &HBSplinesPatchUtility_ListBoundaryBfs<2>)
    .def("ListBoundaryBfs", &HBSplinesPatchUtility_ListBoundaryBfs<3>)
    .def("GetBfByEquationId", &HBSplinesPatchUtility_GetBfByEquationId<2>)
    .def("GetBfByEquationId", &HBSplinesPatchUtility_GetBfByEquationId<3>)
    .def("ReportDuplicatedEquationId", &HBSplinesPatchUtility_ReportDuplicatedEquationId<2>)
    .def("ReportDuplicatedEquationId", &HBSplinesPatchUtility_ReportDuplicatedEquationId<3>)
    .def("__str__", &PrintObject<HBSplinesPatchUtility>)
    ;

    class_<HBSplinesRefinementUtility, typename HBSplinesRefinementUtility::Pointer>
    (m, "HBSplinesRefinementUtility")
    .def(init<>())
    .def("Refine", &HBSplinesRefinementUtility_Refine<2>)
    .def("Refine", &HBSplinesRefinementUtility_Refine<3>)
    .def("Refine", &HBSplinesRefinementUtility_RefineBf<2>)
    .def("Refine", &HBSplinesRefinementUtility_RefineBf<3>)
    .def("RefineWindow", &HBSplinesRefinementUtility_RefineWindow<2>)
    .def("RefineWindow", &HBSplinesRefinementUtility_RefineWindow<3>)
    .def("LinearDependencyRefine", &HBSplinesRefinementUtility_LinearDependencyRefine<2>)
    .def("LinearDependencyRefine", &HBSplinesRefinementUtility_LinearDependencyRefine<3>)
    .def("__str__", &PrintObject<HBSplinesRefinementUtility>)
    ;

    class_<MultiHBSplinesPatchMatlabExporter, MultiHBSplinesPatchMatlabExporter::Pointer>
    (m, "MultiHBSplinesPatchMatlabExporter")
    .def(init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiHBSplinesPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiHBSplinesPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, MultiPatch<3> >)
    .def("__str__", &PrintObject<MultiHBSplinesPatchMatlabExporter>)
    ;

}

}  // namespace Python.

} // Namespace Kratos
