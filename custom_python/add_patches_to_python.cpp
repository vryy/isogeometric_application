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
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"
#include "custom_python/add_import_export_to_python.h"
#include "custom_python/add_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<class TPatchType>
std::size_t Patch_GetId(TPatchType& rDummy)
{
    return rDummy.Id();
}

template<class TPatchType>
void Patch_SetId(TPatchType& rDummy, const std::size_t& Id)
{
    rDummy.SetId(Id);
}

template<class TPatchType>
std::string Patch_GetPrefix(TPatchType& rDummy)
{
    return rDummy.Prefix();
}

template<class TPatchType>
void Patch_SetPrefix(TPatchType& rDummy, const std::string& Prefix)
{
    rDummy.SetPrefix(Prefix);
}

template<class TPatchType>
typename TPatchType::FESpaceType::Pointer Patch_pFESpace(TPatchType& rDummy)
{
    return rDummy.pFESpace();
}

template<int TDim, typename TDataType>
typename GridFunction<TDim, TDataType>::Pointer Patch_CreateGridFunction(Patch<TDim>& rDummy, typename ControlGrid<TDataType>::Pointer pControlGrid)
{
    return rDummy.template CreateGridFunction<TDataType>(pControlGrid);
}

template<int TDim, class TVariableType>
typename GridFunction<TDim, typename TVariableType::Type>::Pointer Patch_GridFunction(Patch<TDim>& rDummy, const TVariableType& rVariable)
{
    return rDummy.template pGetGridFunction<TVariableType>(rVariable);
}

template<class TMultiPatchType>
typename TMultiPatchType::PatchContainerType MultiPatch_GetPatches(TMultiPatchType& rDummy)
{
    return rDummy.Patches();
}

template<class TPatchType, class TMultiPatchType>
typename TPatchType::Pointer MultiPatch_GetItem(TMultiPatchType& rDummy, std::size_t index)
{
    return rDummy.pGetPatch(index);
}

template<class TMultiPatchType>
std::size_t MultiPatch_Len(TMultiPatchType& rDummy)
{
    return rDummy.size();
}

template<int TDim>
std::size_t MultiPatch_Enumerate1(MultiPatch<TDim>& rDummy)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate();
//    KRATOS_WATCH(system_size)

    return system_size;
}

template<int TDim>
std::size_t MultiPatch_Enumerate2(MultiPatch<TDim>& rDummy, const std::size_t& start)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate(start);
//    KRATOS_WATCH(system_size)

    return system_size;
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddPatchesToPython_Helper()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim> >
    // class_<Patch<TDim>, typename Patch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&, typename FESpace<TDim>::Pointer>())
    .add_property("Id", &Patch_GetId<Patch<TDim> >, &Patch_SetId<Patch<TDim> >)
    .add_property("Prefix", &Patch_GetPrefix<Patch<TDim> >, &Patch_SetPrefix<Patch<TDim> >)
    .def("WorkingSpaceDimension", &Patch<TDim>::WorkingSpaceDimension)
    .def("Name", &Patch<TDim>::Name)
    .def("CreateControlPointGridFunction", &Patch<TDim>::CreateControlPointGridFunction)
    .def("CreateGridFunction", &Patch_CreateGridFunction<TDim, double>)
    .def("CreateGridFunction", &Patch_CreateGridFunction<TDim, array_1d<double, 3> >)
    .def("CreateGridFunction", &Patch_CreateGridFunction<TDim, Vector>)
    .def("GridFunction", &Patch_GridFunction<TDim, Variable<ControlPoint<double> > >)
    .def("GridFunction", &Patch_GridFunction<TDim, Variable<double> >)
    .def("GridFunction", &Patch_GridFunction<TDim, Variable<array_1d<double, 3> > >)
    .def("GridFunction", &Patch_GridFunction<TDim, Variable<Vector> >)
    .def("ApplyTransformation", &Patch<TDim>::ApplyTransformation)
    .def("Order", &Patch<TDim>::Order)
    .def("TotalNumber", &Patch<TDim>::TotalNumber)
    .def("FESpace", &Patch_pFESpace<Patch<TDim> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Patch" << TDim << "DPointer";
    class_<typename Patch<TDim>::Pointer>
    (ss.str().c_str(), init<typename Patch<TDim>::Pointer>())
    .def("GetReference", GetReference<Patch<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Patch" << TDim << "DContainer";
    PointerVectorSetPythonInterface<typename MultiPatch<TDim>::PatchContainerType>::CreateInterface(ss.str());

    ss.str(std::string());
    ss << "PatchInterface" << TDim << "D";
    class_<PatchInterface<TDim> >
    // class_<PatchInterface<TDim>, typename PatchInterface<TDim>::Pointer >
    (ss.str().c_str(), init<>())
    .def(init<typename Patch<TDim>::Pointer, const BoundarySide&, typename Patch<TDim>::Pointer, const BoundarySide&>())
    // .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiPatch" << TDim << "D";
    class_<MultiPatch<TDim>, typename MultiPatch<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    // .def("ResetId", &MultiPatch<TDim>::ResetId) // this function is not really useful. One shall keep control over the id of the patch.
    .def("AddPatch", &MultiPatch<TDim>::AddPatch)
    .def("Patches", &MultiPatch_GetPatches<MultiPatch<TDim> >)
    .def("__getitem__", &MultiPatch_GetItem<Patch<TDim>, MultiPatch<TDim> >)
    .def("__len__", &MultiPatch_Len<MultiPatch<TDim> >)
    .def("EquationSystemSize", &MultiPatch<TDim>::EquationSystemSize)
    .def("ResetFunctionIndices", &MultiPatch<TDim>::ResetFunctionIndices)
    .def("Enumerate", &MultiPatch_Enumerate1<TDim>)
    .def("Enumerate", &MultiPatch_Enumerate2<TDim>)
    .def("IsEnumerated", &MultiPatch<TDim>::IsEnumerated)
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
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGLVisExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGLVisExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGLVisExporter, MultiPatch<3> >)
    .def(self_ns::str(self))
    ;

}

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

void IsogeometricApplication_AddPatchesToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////PATCHES///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    enum_<BoundarySide>("BoundarySide")
    .value("Left", _BLEFT_)
    .value("Right", _BRIGHT_)
    .value("Top", _BTOP_)
    .value("Bottom", _BBOTTOM_)
    .value("Front", _BFRONT_)
    .value("Back", _BBACK_)
    ;

    enum_<BoundarySide2D>("BoundarySide2D")
    .value("U0", _B2LEFT_)
    .value("U1", _B2RIGHT_)
    .value("V0", _B2BOTTOM_)
    .value("V1", _B2TOP_)
    ;

    enum_<BoundarySide3D>("BoundarySide3D")
    .value("U0", _B3LEFT_)
    .value("U1", _B3RIGHT_)
    .value("V0", _B3FRONT_)
    .value("V1", _B3BACK_)
    .value("W0", _B3BOTTOM_)
    .value("W1", _B3TOP_)
    ;

    enum_<BoundaryDirection>("BoundaryDirection")
    .value("Forward", _FORWARD_)
    .value("Reversed", _REVERSED_)
    ;

    enum_<BoundaryFlag>("BoundaryFlag")
    .value("Left", _FLEFT_)
    .value("Right", _FRIGHT_)
    .value("Top", _FTOP_)
    .value("Bottom", _FBOTTOM_)
    .value("Front", _FFRONT_)
    .value("Back", _FBACK_)
    ;

    enum_<IsogeometricEchoFlags>("IsogeometricEchoFlags")
    .value("ECHO_REFINEMENT", ECHO_REFINEMENT)
    ;

    IsogeometricApplication_AddPatchesToPython_Helper<1>();
    IsogeometricApplication_AddPatchesToPython_Helper<2>();
    IsogeometricApplication_AddPatchesToPython_Helper<3>();

    /////////////////////////////////////////////////////////////////
    ///////////////////////IMPORT/EXPORT/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddExportToPython();
    IsogeometricApplication_AddImportToPython<1>();
    IsogeometricApplication_AddImportToPython<2>();
    IsogeometricApplication_AddImportToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos

