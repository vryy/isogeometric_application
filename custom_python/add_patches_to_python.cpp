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
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "python/pointer_vector_set_python_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"
#include "custom_python/iga_define_python.h"
#include "custom_python/iga_python_utils.h"
#include "custom_python/add_import_export_to_python.h"
#include "custom_python/add_patches_to_python.h"

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
void Patch_SetId(TPatchType& rDummy, std::size_t Id)
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
int Patch_GetLayerIndex(TPatchType& rDummy)
{
    return rDummy.LayerIndex();
}

template<class TPatchType>
void Patch_SetLayerIndex(TPatchType& rDummy, int Index)
{
    rDummy.SetLayerIndex(Index);
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

template<class TPatchType>
boost::python::list Patch_Predict(TPatchType& rDummy, const boost::python::list& P, const boost::python::list& list_nsampling,
                                  const boost::python::list& list_xi_min, const boost::python::list& list_xi_max)
{
    std::vector<double> xi_min_vec;
    IsogeometricPythonUtils::Unpack<double, double>(list_xi_min, xi_min_vec);

    std::vector<double> xi_max_vec;
    IsogeometricPythonUtils::Unpack<double, double>(list_xi_max, xi_max_vec);

    std::vector<double> P_vec;
    IsogeometricPythonUtils::Unpack<double, double>(P, P_vec);

    array_1d<double, 3> point, xi, xi_min, xi_max;
    noalias(point) = ZeroVector(3);
    noalias(xi) = ZeroVector(3);
    noalias(xi_min) = ZeroVector(3);
    noalias(xi_max) = ZeroVector(3);

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), xi_min_vec.size()); ++i)
    {
        xi_min[i] = xi_min_vec[i];
    }

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), xi_max_vec.size()); ++i)
    {
        xi_max[i] = xi_max_vec[i];
    }

    std::vector<int> nsampling;
    IsogeometricPythonUtils::Unpack<int, int>(list_nsampling, nsampling);

    rDummy.Predict(point, xi, nsampling, xi_min, xi_max);

    boost::python::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    return out_point;
}

template<class TPatchType>
boost::python::list Patch_Predict1(TPatchType& rDummy, const boost::python::list& P, const boost::python::list& list_nsampling)
{
    std::vector<double> P_vec;
    IsogeometricPythonUtils::Unpack<double, double>(P, P_vec);

    array_1d<double, 3> point, xi;
    noalias(point) = ZeroVector(3);
    noalias(xi) = ZeroVector(3);

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    std::vector<int> nsampling;
    IsogeometricPythonUtils::Unpack<int, int>(list_nsampling, nsampling);

    rDummy.Predict(point, xi, nsampling);

    boost::python::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    return out_point;
}

template<class TPatchType>
boost::python::list Patch_LocalCoordinates(TPatchType& rDummy, const boost::python::list& P, const boost::python::list& xi0)
{
    std::vector<double> xi0_vec;
    IsogeometricPythonUtils::Unpack<double, double>(xi0, xi0_vec);

    std::vector<double> P_vec;
    IsogeometricPythonUtils::Unpack<double, double>(P, P_vec);

    array_1d<double, 3> point, xi;
    noalias(point) = ZeroVector(3);
    noalias(xi) = ZeroVector(3);

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), xi0_vec.size()); ++i)
    {
        xi[i] = xi0_vec[i];
    }

    int stat = rDummy.LocalCoordinates(point, xi);

    boost::python::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    boost::python::list output;
    output.append(stat);
    output.append(out_point);
    return output;
}

template<class TPatchType>
bool Patch_IsInside(TPatchType& rDummy, const boost::python::list& P, const boost::python::list& xi0)
{
    std::vector<double> xi0_vec;
    IsogeometricPythonUtils::Unpack<double, double>(xi0, xi0_vec);

    std::vector<double> P_vec;
    IsogeometricPythonUtils::Unpack<double, double>(P, P_vec);

    array_1d<double, 3> point, xi;
    noalias(point) = ZeroVector(3);
    noalias(xi) = ZeroVector(3);

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), xi0_vec.size()); ++i)
    {
        xi[i] = xi0_vec[i];
    }

    return rDummy.IsInside(point, xi);
}

template<int TDim>
typename PatchInterface<TDim>::Pointer Patch_GetInterface(Patch<TDim>& rDummy, std::size_t i)
{
    return rDummy.pInterface(i);
}

template<int TDim>
typename Patch < TDim - 1 >::Pointer Patch_ConstructBoundaryPatch(Patch<TDim>& rDummy, std::size_t iside)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.ConstructBoundaryPatch(side);
}

template<int TDim>
static typename MultiPatch<TDim>::Pointer MultiPatch_init(const boost::python::list& patch_list)
{
    typename MultiPatch<TDim>::Pointer pMultiPatch = typename MultiPatch<TDim>::Pointer(new MultiPatch<TDim>());

    std::vector<typename Patch<TDim>::Pointer> patches;
    IsogeometricPythonUtils::Unpack<typename Patch<TDim>::Pointer>(patch_list, patches);

    for (std::size_t i = 0; i < patches.size(); ++i)
        pMultiPatch->AddPatch(patches[i]);

    return pMultiPatch;
}

template<class TMultiPatchType>
typename TMultiPatchType::PatchContainerType MultiPatch_GetPatches(TMultiPatchType& rDummy)
{
    return rDummy.Patches();
}

template<class TMultiPatchType>
boost::python::list MultiPatch_GetPatchIndices(TMultiPatchType& rDummy)
{
    boost::python::list ids;
    for (typename TMultiPatchType::patch_const_iterator it = rDummy.Patches().begin();
            it != rDummy.Patches().end(); ++it)
    {
        ids.append(it->Id());
    }
    return ids;
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
std::size_t MultiPatch_Enumerate2(MultiPatch<TDim>& rDummy, std::size_t start)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate(start);
//    KRATOS_WATCH(system_size)

    return system_size;
}

template<int TDim>
typename Patch<TDim>::Pointer PatchInterface_pPatch1(PatchInterface<TDim>& rDummy)
{
    return rDummy.pPatch1();
}

template<class TPatchInterfaceType>
BoundarySide PatchInterface_Side1(TPatchInterfaceType& rDummy)
{
    return rDummy.Side1();
}

template<int TDim>
typename Patch<TDim>::Pointer PatchInterface_pPatch2(PatchInterface<TDim>& rDummy)
{
    return rDummy.pPatch2();
}

template<class TPatchInterfaceType>
BoundarySide PatchInterface_Side2(TPatchInterfaceType& rDummy)
{
    return rDummy.Side2();
}

////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddPatchesToPython_Helper()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim>, bases<Flags> >
    // class_<Patch<TDim>, typename Patch<TDim>::Pointer > // do not use this to export Patch pointer
    (ss.str().c_str(), init<std::size_t, typename FESpace<TDim>::Pointer>())
    .add_property("Id", &Patch_GetId<Patch<TDim> >, &Patch_SetId<Patch<TDim> >)
    .add_property("Prefix", &Patch_GetPrefix<Patch<TDim> >, &Patch_SetPrefix<Patch<TDim> >)
    .add_property("LayerIndex", &Patch_GetLayerIndex<Patch<TDim> >, &Patch_SetLayerIndex<Patch<TDim> >)
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
    .def("Predict", &Patch_Predict<Patch<TDim> >)
    .def("Predict", &Patch_Predict1<Patch<TDim> >)
    .def("LocalCoordinates", &Patch_LocalCoordinates<Patch<TDim> >)
    .def("IsInside", &Patch_IsInside<Patch<TDim> >)
    .def("NumberOfInterfaces", &Patch<TDim>::NumberOfInterfaces)
    .def("AddInterface", &Patch<TDim>::AddInterface)
    .def("RemoveInterface", &Patch<TDim>::RemoveInterface)
    .def("GetInterface", &Patch_GetInterface<TDim>)
    .def("ConstructBoundaryPatch", &Patch_ConstructBoundaryPatch<TDim>)
    .def("ConstructSlicedPatch", &Patch<TDim>::ConstructSlicedPatch)
    .def("FindBoundarySide", &Patch<TDim>::FindBoundarySide)
    .def("SetLocalSearchTolerance", &Patch<TDim>::SetLocalSearchTolerance)
    .def("SetLocalSearchMaxIters", &Patch<TDim>::SetLocalSearchMaxIters)
    .def("Validate", &Patch<TDim>::Validate)
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
    // class_<PatchInterface<TDim> >
    class_<PatchInterface<TDim>, typename PatchInterface<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def(init<typename Patch<TDim>::Pointer, const BoundarySide&, typename Patch<TDim>::Pointer, const BoundarySide&>())
    .def("Patch1", &PatchInterface_pPatch1<TDim>)
    .def("Patch2", &PatchInterface_pPatch2<TDim>)
    .def("Side1", &PatchInterface_Side1<PatchInterface<TDim> >)
    .def("Side2", &PatchInterface_Side2<PatchInterface<TDim> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "MultiPatch" << TDim << "D";
    class_<MultiPatch<TDim>, typename MultiPatch<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__init__", make_constructor(&MultiPatch_init<TDim>))
    // .def("ResetId", &MultiPatch<TDim>::ResetId) // this function is not really useful. One shall keep control over the id of the patch.
    .def("AddPatch", &MultiPatch<TDim>::AddPatch)
    .def("RemovePatch", &MultiPatch<TDim>::RemovePatch)
    .def("Patches", &MultiPatch_GetPatches<MultiPatch<TDim> >)
    .def("PatchIndices", &MultiPatch_GetPatchIndices<MultiPatch<TDim> >)
    .def("__getitem__", &MultiPatch_GetItem<Patch<TDim>, MultiPatch<TDim> >)
    .def("__len__", &MultiPatch_Len<MultiPatch<TDim> >)
    .def("EquationSystemSize", &MultiPatch<TDim>::EquationSystemSize)
    .def("ResetFunctionIndices", &MultiPatch<TDim>::ResetFunctionIndices)
    .def("Enumerate", &MultiPatch_Enumerate1<TDim>)
    .def("Enumerate", &MultiPatch_Enumerate2<TDim>)
    .def("IsEnumerated", &MultiPatch<TDim>::IsEnumerated)
    .def("LocalCoordinates", &Patch_LocalCoordinates<MultiPatch<TDim> >)
    .def("Validate", &MultiPatch<TDim>::Validate)
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

    enum_<ParametricAxis>("ParametricAxis")
    .value("U", _PA_U_)
    .value("V", _PA_V_)
    .value("W", _PA_W_)
    ;

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
    .value("Forward", BoundaryDirection::_FORWARD_)
    .value("Reversed", BoundaryDirection::_REVERSED_)
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
    .value("ECHO_REFINEMENT_DETAIL", ECHO_REFINEMENT_DETAIL)
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
