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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
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
    const std::string& fn, const int& division1, const int& division2)
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

void IsogeometricApplication_AddTSplinesToPython(pybind11::module& m)
{

    /////////////////////////////////////////////////////////////////
    /////////////////////////TMESH///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    pybind11::class_<TsMesh2D, TsMesh2D::Pointer>(m, "TsMesh2D")
    .def(pybind11::init<>())
    .def("BeginConstruct", &TsMesh2D::BeginConstruct)
    .def("EndConstruct", &TsMesh2D::EndConstruct)
    .def("BuildExtendedTmesh", &TsMesh2D::BuildExtendedTmesh)
    .def("IsAnalysisSuitable", &TsMesh2D::IsAnalysisSuitable)
    .def("BuildAnchors", &TsMesh2D::BuildAnchors)
    .def("BuildCells", &TsMesh2D::BuildCells)
    //    .def("FindKnots2", &TsMesh2D::FindKnots2)
    .def("__str__", &PrintObject<TsMesh2D>)
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////T-SPLINES FESPACE///////////////////////
    /////////////////////////////////////////////////////////////////

    typedef PBBSplinesBasisFunction<2, TCell> PBBSplinesBasisFunctionType;
    typedef PBBSplinesFESpace<2, PBBSplinesBasisFunctionType, BCellManager<2, typename PBBSplinesBasisFunctionType::CellType> > PBBSplinesFESpaceType;
    typedef TSplinesFESpace<2, PBBSplinesBasisFunctionType, BCellManager<2, typename PBBSplinesBasisFunctionType::CellType> > TSplinesFESpaceType;
    pybind11::class_<TSplinesFESpaceType, TSplinesFESpaceType::Pointer, PBBSplinesFESpaceType>
    (m, "TSplinesFESpace2D")
    .def(pybind11::init<>())
    .def("__str__", &PrintObject<TSplinesFESpaceType>)
    ;

    typedef NonConformingTSplinesMultipatchLagrangeMesh<TSplinesFESpaceType> NonConformingTSplinesMultipatchLagrangeMeshType;
    pybind11::class_<NonConformingTSplinesMultipatchLagrangeMeshType, typename NonConformingTSplinesMultipatchLagrangeMeshType::Pointer>
    (m, "NonConformingTSplinesMultipatchLagrangeMesh")
    .def(pybind11::init<typename MultiPatch<2>::Pointer>())
    .def("SetBaseElementName", &NonConformingTSplinesMultipatchLagrangeMeshType::SetBaseElementName)
    .def("SetLastNodeId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastNodeId)
    .def("SetLastElemId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastElemId)
    .def("SetLastPropId", &NonConformingTSplinesMultipatchLagrangeMeshType::SetLastPropId)
    .def("SetDivision", &NonConformingTSplinesMultipatchLagrangeMeshType::SetDivision)
    .def("SetUniformDivision", &NonConformingTSplinesMultipatchLagrangeMeshType::SetUniformDivision)
    .def("WriteModelPart", &NonConformingTSplinesMultipatchLagrangeMeshType::WriteModelPart)
    .def("__str__", &PrintObject<NonConformingTSplinesMultipatchLagrangeMeshType>)
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////UTILITIES///////////////////////////////
    /////////////////////////////////////////////////////////////////

    pybind11::class_<TSplinesUtils, TSplinesUtils::Pointer>(m, "TSplinesUtils")
    .def(pybind11::init<>())
    .def("CreateFromBSplines", &TSplinesUtils_CreateFromBSplines)
    .def("ReadFromFile", &TSplinesUtils_ReadFromFile)
    .def("ExportMatlab", &TSplinesUtils_ExportMatlab)
    .def("ExportMDPA", &TSplinesUtils_ExportMDPA)
    .def("ExportMDPA", &TSplinesUtils_ExportMDPA2)
    .def("__str__", &PrintObject<TSplinesUtils>)
    ;

    /////////////////////////////////////////////////////////////////
    /////////////////////////IMPORT/EXPORT///////////////////////////
    /////////////////////////////////////////////////////////////////

    #ifdef ISOGEOMETRIC_USE_TSPLINE
    pybind11::class_<TSplinesPatchTSMImporter, TSplinesPatchTSMImporter::Pointer>
    (m, "TSplinesPatchTSMImporter")
    .def(pybind11::init<>())
    .def("ImportSingle", &TSplinesPatchTSMImporter::ImportSingle)
    .def("__str__", &PrintObject<TSplinesPatchTSMImporter>)
    ;
    #endif

}

}  // namespace Python.

} // Namespace Kratos

