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
#include "custom_utilities/multipatch.h"
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/tsplines/tcell.h"
#include "custom_utilities/tsplines/tsplines_utils.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_utilities/tsplines/tsplines_fespace.h"
#include "custom_utilities/tsplines/nonconforming_tsplines_multipatch_lagrange_mesh.h"
#include "custom_python/add_tsplines_to_python.h"
#include "custom_python/add_import_export_to_python.h"
#ifdef ISOGEOMETRIC_USE_TSPLINE
#include "custom_utilities/import_export/tsplines_patch_tsm_importer.h"
#endif


namespace Kratos
{

void TSplinesUtils_CreateFromBSplines(TSplinesUtils& rDummy, TsMesh2D& tmesh,
    const BSplinesFESpace<2>& rFESpace)
{
    rDummy.CreateFromBSplines(tmesh, rFESpace);
}

void TSplinesUtils_ReadFromFile(TSplinesUtils& rDummy, TsMesh2D& tmesh,
    const std::string& fn)
{
    rDummy.ReadFromFile(tmesh, fn);
}

void TSplinesUtils_ExportMatlab(TSplinesUtils& rDummy, TsMesh2D& tmesh,
    const std::string& fn, const std::string& mesh_type)
{
    rDummy.ExportMatlab(tmesh, fn, mesh_type);
}

void TSplinesUtils_ExportMDPA(TSplinesUtils& rDummy, const TsMesh2D& tmesh,
    const std::string& fn, int division1, int division2)
{
    rDummy.ExportMDPA(tmesh, fn, std::vector<int>{division1, division2});
}

void TSplinesUtils_ExportMDPA2(TSplinesUtils& rDummy, const TsMesh2D& tmesh,
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
    /////////////////////////T-SPLINES FESPACE///////////////////////
    /////////////////////////////////////////////////////////////////

    typedef PBBSplinesBasisFunction<2, TCell> PBBSplinesBasisFunctionType;
    typedef PBBSplinesFESpace<2, PBBSplinesBasisFunctionType, BCellManager<2, typename PBBSplinesBasisFunctionType::CellType> > PBBSplinesFESpaceType;
    typedef TSplinesFESpace<2, PBBSplinesBasisFunctionType, BCellManager<2, typename PBBSplinesBasisFunctionType::CellType> > TSplinesFESpaceType;
    class_<TSplinesFESpaceType, TSplinesFESpaceType::Pointer, bases<PBBSplinesFESpaceType>, boost::noncopyable>
    ("TSplinesFESpace2D", init<>())
    .def(self_ns::str(self))
    ;

    typedef NonConformingTSplinesMultipatchLagrangeMesh<TSplinesFESpaceType> NonConformingTSplinesMultipatchLagrangeMeshType;
    class_<NonConformingTSplinesMultipatchLagrangeMeshType, typename NonConformingTSplinesMultipatchLagrangeMeshType::Pointer, boost::noncopyable>
    ("NonConformingTSplinesMultipatchLagrangeMesh", init<typename MultiPatch<2>::Pointer>())
    .def("SetBaseElementName", &NonConformingTSplinesMultipatchLagrangeMeshType::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastNodeId)
    .def("SetLastElemId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastElemId)
    .def("SetLastPropId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastPropId)
    .def("SetDivision", &NonConformingTSplinesMultipatchLagrangeMeshType::SetDivision)
    .def("SetUniformDivision", &NonConformingTSplinesMultipatchLagrangeMeshType::SetUniformDivision)
    .def("WriteModelPart", &NonConformingTSplinesMultipatchLagrangeMeshType::WriteModelPart)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////UTILITIES///////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<TSplinesUtils, TSplinesUtils::Pointer, boost::noncopyable>
    ("TSplinesUtils", init<>())
    .def("CreateFromBSplines", &TSplinesUtils_CreateFromBSplines)
    .def("ReadFromFile", &TSplinesUtils_ReadFromFile)
    .def("ExportMatlab", &TSplinesUtils_ExportMatlab)
    .def("ExportMDPA", &TSplinesUtils_ExportMDPA)
    .def("ExportMDPA", &TSplinesUtils_ExportMDPA2)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////IMPORT/EXPORT///////////////////////////
    /////////////////////////////////////////////////////////////////

    #ifdef ISOGEOMETRIC_USE_TSPLINE
    class_<TSplinesPatchTSMImporter, TSplinesPatchTSMImporter::Pointer, boost::noncopyable>
    ("TSplinesPatchTSMImporter", init<>())
    .def("ImportSingle", &TSplinesPatchTSMImporter::ImportSingle)
    .def(self_ns::str(self))
    ;
    #endif

}

}  // namespace Python.

} // Namespace Kratos

