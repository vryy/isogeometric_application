/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Aug 6, 2026 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"
#include "custom_utilities/import_export/multi_pbbsplines_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_hbsplines_patch_matlab_exporter.h"
#ifdef ISOGEOMETRIC_USE_TSPLINE
#include "custom_utilities/import_export/tsplines_patch_tsm_importer.h"
#endif
#include "custom_python/iga_define_python.h"
#include "custom_python/add_import_export_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddImportToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "MultiNURBSPatchGeoImporter" << TDim << "D";
    class_<MultiNURBSPatchGeoImporter<TDim>, typename MultiNURBSPatchGeoImporter<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("ImportSingle", &MultiNURBSPatchGeoImporter<TDim>::ImportSingle)
    .def("Import", &MultiNURBSPatchGeoImporter<TDim>::Import)
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddExportToPython()
{
    class_<MultiNURBSPatchGeoExporter, MultiNURBSPatchGeoExporter::Pointer, boost::noncopyable>
    ("MultiNURBSPatchGeoExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

    class_<MultiNURBSPatchMatlabExporter, MultiNURBSPatchMatlabExporter::Pointer, boost::noncopyable>
    ("MultiNURBSPatchMatlabExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

    class_<MultiNURBSPatchGLVisExporter, MultiNURBSPatchGLVisExporter::Pointer, boost::noncopyable>
    ("MultiNURBSPatchGLVisExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGLVisExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export_Variable<1, MultiNURBSPatchGLVisExporter, Patch<1>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<1, MultiNURBSPatchGLVisExporter, Patch<1>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<1, MultiNURBSPatchGLVisExporter, Patch<1>, Variable<Vector> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export_Variable<2, MultiNURBSPatchGLVisExporter, Patch<2>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<2, MultiNURBSPatchGLVisExporter, Patch<2>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<2, MultiNURBSPatchGLVisExporter, Patch<2>, Variable<Vector> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export_Variable<3, MultiNURBSPatchGLVisExporter, Patch<3>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<3, MultiNURBSPatchGLVisExporter, Patch<3>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<3, MultiNURBSPatchGLVisExporter, Patch<3>, Variable<Vector> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGLVisExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export_Variable<1, MultiNURBSPatchGLVisExporter, MultiPatch<1>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<1, MultiNURBSPatchGLVisExporter, MultiPatch<1>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<1, MultiNURBSPatchGLVisExporter, MultiPatch<1>, Variable<Vector> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export_Variable<2, MultiNURBSPatchGLVisExporter, MultiPatch<2>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<2, MultiNURBSPatchGLVisExporter, MultiPatch<2>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<2, MultiNURBSPatchGLVisExporter, MultiPatch<2>, Variable<Vector> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, MultiPatch<3> >)
    .def("Export", &MultiPatchExporter_Export_Variable<3, MultiNURBSPatchGLVisExporter, MultiPatch<3>, Variable<double> >)
    .def("Export", &MultiPatchExporter_Export_Variable<3, MultiNURBSPatchGLVisExporter, MultiPatch<3>, Variable<array_1d<double, 3> > >)
    .def("Export", &MultiPatchExporter_Export_Variable_WithComponents<3, MultiNURBSPatchGLVisExporter, MultiPatch<3>, Variable<Vector> >)
    .def(self_ns::str(self))
    ;

    class_<MultiPBBSplinesPatchMatlabExporter, MultiPBBSplinesPatchMatlabExporter::Pointer, boost::noncopyable>
    ("MultiPBBSplinesPatchMatlabExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiPBBSplinesPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiPBBSplinesPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiPBBSplinesPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiPBBSplinesPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiPBBSplinesPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiPBBSplinesPatchMatlabExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

    class_<MultiHBSplinesPatchMatlabExporter, MultiHBSplinesPatchMatlabExporter::Pointer, boost::noncopyable>
    ("MultiHBSplinesPatchMatlabExporter", init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiHBSplinesPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiHBSplinesPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiHBSplinesPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiHBSplinesPatchMatlabExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

#ifdef ISOGEOMETRIC_USE_TSPLINE
    class_<TSplinesPatchTSMImporter, TSplinesPatchTSMImporter::Pointer, boost::noncopyable>
    ("TSplinesPatchTSMImporter", init<>())
    .def("ImportSingle", &TSplinesPatchTSMImporter::ImportSingle)
    .def(self_ns::str(self))
    ;
#endif
}

/// template instantiation
template void IsogeometricApplication_AddImportToPython<1>();
template void IsogeometricApplication_AddImportToPython<2>();
template void IsogeometricApplication_AddImportToPython<3>();

}  // namespace Python.

} // Namespace Kratos
