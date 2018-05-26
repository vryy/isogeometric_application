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
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "python/pointer_vector_set_python_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/hbsplines/deprecated_hb_mesh.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_utilities/hbsplines/hbsplines_patch_utility.h"
#include "custom_utilities/hbsplines/hbsplines_refinement_utility.h"
#include "custom_utilities/import_export/multi_hbsplines_patch_matlab_exporter.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_point_based_control_grid_to_python.h"
#include "custom_python/add_import_export_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim>
std::size_t HBSplinesBasisFunction_GetId(HBSplinesBasisFunction<TDim>& rDummy)
{
    return rDummy.Id();
}

template<int TDim>
void HBSplinesBasisFunction_SetId(HBSplinesBasisFunction<TDim>& rDummy, std::size_t Id)
{
    // DO NOTHING
}

template<int TDim>
std::size_t HBSplinesBasisFunction_GetEquationId(HBSplinesBasisFunction<TDim>& rDummy)
{
    return rDummy.EquationId();
}

template<int TDim>
void HBSplinesBasisFunction_SetEquationId(HBSplinesBasisFunction<TDim>& rDummy, std::size_t EquationId)
{
    rDummy.SetEquationId(EquationId);
}

template<int TDim>
boost::python::list HBSplinesFESpace_GetBoundaryBfs(HBSplinesFESpace<TDim>& rDummy, std::size_t boundary_id)
{
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;

    std::vector<bf_t> bf_list = rDummy.GetBoundaryBfs(boundary_id);

    boost::python::list Output;
    for (std::size_t i = 0; i < bf_list.size(); ++i)
        Output.append(bf_list[i]);

    return Output;
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

////////////////////////////////////////

template<int TDim>
void HBSplinesRefinementUtility_Refine(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, const std::size_t& Id, const int& EchoLevel)
{
    rDummy.Refine<TDim>(pPatch, Id, EchoLevel);
}

template<int TDim>
void HBSplinesRefinementUtility_RefineWindow(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, boost::python::list& window, const int& EchoLevel)
{

    std::vector<std::vector<double> > window_vector;
    std::size_t cnt1 = 0, cnt2 = 0;
    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& vect, std::make_pair(iterator_value_type(window), iterator_value_type() ) )
    {
        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        std::vector<double> win_vect;
        BOOST_FOREACH(const iterator_value_type2::value_type& v, std::make_pair(iterator_value_type2(vect), iterator_value_type2() ) )
        {
            win_vect.push_back(v);
        }
        window_vector.push_back(win_vect);
    }
    rDummy.RefineWindow<TDim>(pPatch, window_vector, EchoLevel);
}

template<int TDim>
void HBSplinesRefinementUtility_LinearDependencyRefine(HBSplinesRefinementUtility& rDummy,
        typename Patch<TDim>::Pointer pPatch, const std::size_t& refine_cycle, const int& EchoLevel)
{
    rDummy.LinearDependencyRefine<TDim>(pPatch, refine_cycle, EchoLevel);
}

////////////////////////////////////////

template<typename TDataType, class TFESpaceType>
typename ControlGrid<TDataType>::Pointer ControlGridUtility_CreatePointBasedControlGrid(
        ControlGridUtility& rDummy,
        const Variable<TDataType>& rVariable, typename TFESpaceType::Pointer pFESpace)
{
    return rDummy.CreatePointBasedControlGrid<TDataType, TFESpaceType>(rVariable, pFESpace);
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddHBSplinesSpaceToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "HBSplinesBasisFunction" << TDim << "D";
    class_<HBSplinesBasisFunction<TDim>, typename HBSplinesBasisFunction<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<const std::size_t&, const std::size_t&>())
    .add_property("Id", HBSplinesBasisFunction_GetId<TDim>, HBSplinesBasisFunction_SetId<TDim>)
    .add_property("EquationId", HBSplinesBasisFunction_GetEquationId<TDim>, HBSplinesBasisFunction_SetEquationId<TDim>)
    .def("Weight", &HBSplinesBasisFunction<TDim>::Weight)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "HBSplinesFESpace" << TDim << "D";
    typename FESpace<TDim-1>::Pointer(HBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace1)(const BoundarySide& side) const = &HBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    // typename FESpace<TDim-1>::Pointer(HBSplinesFESpace<TDim>::*pointer_to_ConstructBoundaryFESpace2)(const BoundarySide& side, const BoundaryRotation& rotation) const = &HBSplinesFESpace<TDim>::ConstructBoundaryFESpace;
    class_<HBSplinesFESpace<TDim>, typename HBSplinesFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__getitem__", &HBSplinesFESpace_GetItem<TDim>)
    .def("GetBoundaryBfs", &HBSplinesFESpace_GetBoundaryBfs<TDim>)
    .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace1)
    // .def("ConstructBoundaryFESpace", pointer_to_ConstructBoundaryFESpace2)
    .def(self_ns::str(self))
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, HBSplinesFESpace<TDim> >::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, HBSplinesFESpace<TDim> >::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, HBSplinesFESpace<TDim> >::Execute();

    ////////////////////OLD H-SPLINES////////////////////////

    ss.str(std::string());
    ss << "DeprecatedHBMesh" << TDim << "D";
    class_<DeprecatedHBMesh<TDim>, bases<Patch<TDim> > >
    // class_<DeprecatedHBMesh<TDim>, typename DeprecatedHBMesh<TDim>::Pointer, bases<Patch<TDim> > >
    (ss.str().c_str(), init<const std::size_t&, const std::string&>())
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
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "DeprecatedHBMesh" << TDim << "DPointer";
    class_<typename DeprecatedHBMesh<TDim>::Pointer>
    (ss.str().c_str(), init<typename DeprecatedHBMesh<TDim>::Pointer>())
    .def("GetReference", GetReference<DeprecatedHBMesh<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

}

////////////////////////////////////////

void IsogeometricApplication_AddHBSplinesToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////HIERARCHICAL BSplines/////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddHBSplinesSpaceToPython<2>();
    IsogeometricApplication_AddHBSplinesSpaceToPython<3>();

    class_<ControlGridUtility, ControlGridUtility::Pointer, boost::noncopyable>
    ("PointBasedControlGridUtility", init<>())
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<double, HBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<3> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<1> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<2> >)
    .def("CreatePointBasedControlGrid", &ControlGridUtility_CreatePointBasedControlGrid<Vector, HBSplinesFESpace<3> >)
    ;

    class_<HBSplinesPatchUtility, HBSplinesPatchUtility::Pointer, boost::noncopyable>
    ("HBSplinesPatchUtility", init<>())
    .def("CreatePatchFromBSplines", &HBSplinesPatchUtility_CreatePatchFromBSplines<2>)
    .def("CreatePatchFromBSplines", &HBSplinesPatchUtility_CreatePatchFromBSplines<3>)
    .def("ListBoundaryBfs", &HBSplinesPatchUtility_ListBoundaryBfs<2>)
    .def("ListBoundaryBfs", &HBSplinesPatchUtility_ListBoundaryBfs<3>)
    .def(self_ns::str(self))
    ;

    class_<HBSplinesRefinementUtility, typename HBSplinesRefinementUtility::Pointer, boost::noncopyable>
    ("HBSplinesRefinementUtility", init<>())
    .def("Refine", &HBSplinesRefinementUtility_Refine<2>)
    .def("Refine", &HBSplinesRefinementUtility_Refine<3>)
    .def("RefineWindow", &HBSplinesRefinementUtility_RefineWindow<2>)
    .def("RefineWindow", &HBSplinesRefinementUtility_RefineWindow<3>)
    .def("LinearDependencyRefine", &HBSplinesRefinementUtility_LinearDependencyRefine<2>)
    .def("LinearDependencyRefine", &HBSplinesRefinementUtility_LinearDependencyRefine<3>)
    .def(self_ns::str(self))
    ;

    class_<MultiHBSplinesPatchMatlabExporter, MultiHBSplinesPatchMatlabExporter::Pointer, boost::noncopyable>
    ("MultiHBSplinesPatchMatlabExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

