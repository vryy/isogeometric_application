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
#include "custom_utilities/tsplines/tspline_utils.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_python/add_import_export_to_python.h"


namespace Kratos
{

void TSplineUtils_CreateFromBSplines(TSplineUtils& rDummy, TsMesh2D& tmesh,
    const BSplinesFESpace<2>& rFESpace)
{
    rDummy.CreateFromBSplines(tmesh, rFESpace);
}

void TSplineUtils_ReadFromFile(TSplineUtils& rDummy, TsMesh2D& tmesh,
    const std::string& fn)
{
    rDummy.ReadFromFile(tmesh, fn);
}

void TSplineUtils_ExportMatlab(TSplineUtils& rDummy, TsMesh2D& tmesh,
    const std::string& fn, const std::string& mesh_type)
{
    rDummy.ExportMatlab(tmesh, fn, mesh_type);
}

void TSplineUtils_ExportMDPA(TSplineUtils& rDummy, const TsMesh2D& tmesh,
    const std::string& fn, const int& division1, const int& division2)
{
    rDummy.ExportMDPA(tmesh, fn, std::vector<int>{division1, division2});
}

void TSplineUtils_ExportMDPA2(TSplineUtils& rDummy, const TsMesh2D& tmesh,
    const std::string& fn)
{
    rDummy.ExportMDPA(tmesh, fn, std::vector<int>{});
}

namespace Python
{

using namespace boost::python;

void IsogeometricApplication_AddTSplinesToPython()
{

    /////////////////////////////////////////////////////////////////
    /////////////////////////TMESH///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<TsMesh2D, TsMesh2D::Pointer, boost::noncopyable>
    ("TsMesh2D", init<>())
    .def("BeginConstruct", &TsMesh2D::BeginConstruct)
    .def("EndConstruct", &TsMesh2D::EndConstruct)
    .def("BuildExtendedTmesh", &TsMesh2D::BuildExtendedTmesh)
    .def("IsAnalysisSuitable", &TsMesh2D::IsAnalysisSuitable)
    .def("BuildAnchors", &TsMesh2D::BuildAnchors)
    .def("BuildCells", &TsMesh2D::BuildCells)
    //    .def("FindKnots2", &TsMesh2D::FindKnots2)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////UTILITIES///////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<TSplineUtils, TSplineUtils::Pointer, boost::noncopyable>
    ("TSplineUtils", init<>())
    .def("CreateFromBSplines", &TSplineUtils_CreateFromBSplines)
    .def("ReadFromFile", &TSplineUtils_ReadFromFile)
    .def("ExportMatlab", &TSplineUtils_ExportMatlab)
    .def("ExportMDPA", &TSplineUtils_ExportMDPA)
    .def("ExportMDPA", &TSplineUtils_ExportMDPA2)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

