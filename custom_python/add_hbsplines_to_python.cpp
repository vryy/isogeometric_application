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
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/control_grid_utility.h"
#include "custom_utilities/hbsplines/deprecated_hb_mesh.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_utilities/hbsplines/hbsplines_patch_utility.h"
#include "custom_utilities/hbsplines/hbsplines_refinement_utility.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_python/add_point_based_control_grid_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

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
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "HBSplinesFESpace" << TDim << "D";
    class_<HBSplinesFESpace<TDim>, typename HBSplinesFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<>())
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

void IsogeometricApplication_AddHBSplinesToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////TSPLINES//////////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<TsMesh2D, TsMesh2D::Pointer, boost::noncopyable>
    ("TsMesh2D", init<>())
    .def("BeginConstruct", &TsMesh2D::BeginConstruct)
    .def("EndConstruct", &TsMesh2D::EndConstruct)
    .def("ReadFromFile", &TsMesh2D::ReadFromFile)
    .def("ExportMatlab", &TsMesh2D::ExportMatlab)
    .def("BuildExtendedTmesh", &TsMesh2D::BuildExtendedTmesh)
    .def("IsAnalysisSuitable", &TsMesh2D::IsAnalysisSuitable)
    .def("BuildAnchors", &TsMesh2D::BuildAnchors)
    .def("BuildCells", &TsMesh2D::BuildCells)
    .def("ExportMDPA", &TsMesh2D::ExportMDPA)
    //    .def("FindKnots2", &TsMesh2D::FindKnots2)
    .def(self_ns::str(self))
    ;

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

}

}  // namespace Python.

} // Namespace Kratos

