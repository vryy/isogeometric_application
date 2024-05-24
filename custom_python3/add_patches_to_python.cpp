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

// Project includes
#include "includes/define_python.h"
#include "python/containers_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"
#include "custom_python/iga_define_python.h"
#include "custom_python/add_import_export_to_python.h"
#include "custom_python/add_patches_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

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
pybind11::list Patch_Predict(TPatchType& rDummy, pybind11::list& P, pybind11::list& list_nsampling,
                             pybind11::list& list_xi_min, pybind11::list& list_xi_max)
{
    std::vector<double> xi_min_vec;
    for (auto v : list_xi_min)
    {
        xi_min_vec.push_back(v.cast<double>());
    }

    std::vector<double> xi_max_vec;
    for (auto v : list_xi_max)
    {
        xi_max_vec.push_back(v.cast<double>());
    }

    std::vector<double> P_vec;
    for (auto v : P)
    {
        P_vec.push_back(v.cast<double>());
    }

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
    for (auto v : list_nsampling)
    {
        nsampling.push_back(v.cast<int>());
    }

    rDummy.Predict(point, xi, nsampling, xi_min, xi_max);

    pybind11::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    return out_point;
}

template<class TPatchType>
pybind11::list Patch_LocalCoordinates(TPatchType& rDummy, pybind11::list& P, pybind11::list& xi0)
{
    std::vector<double> xi0_vec;
    for (auto v : xi0)
    {
        xi0_vec.push_back(v.cast<double>());
    }

    std::vector<double> P_vec;
    for (auto v : P)
    {
        P_vec.push_back(v.cast<double>());
    }

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

    pybind11::list out_point;
    out_point.append(xi[0]);
    out_point.append(xi[1]);
    out_point.append(xi[2]);

    pybind11::list output;
    output.append(stat);
    output.append(out_point);
    return output;
}

template<class TPatchType>
bool Patch_IsInside(TPatchType& rDummy, pybind11::list& P, pybind11::list& xi0)
{
    std::vector<double> xi0_vec;
    for (auto v : xi0)
    {
        xi0_vec.push_back(v.cast<double>());
    }

    std::vector<double> P_vec;
    for (auto v : P)
    {
        P_vec.push_back(v.cast<double>());
    }

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

#ifdef SD_APP_FORWARD_COMPATIBILITY
template<int TDim>
void MultiPatch_AddPatch(MultiPatch<TDim>& rDummy, iga::Wrapper<Patch<TDim>, typename Patch<TDim>::Pointer>& patch_wrapper)
{
    rDummy.AddPatch(patch_wrapper.GetPointer());
}

template<int TDim>
void MultiPatch_RemovePatch(MultiPatch<TDim>& rDummy, iga::Wrapper<Patch<TDim>, typename Patch<TDim>::Pointer>& patch_wrapper)
{
    rDummy.RemovePatch(patch_wrapper.GetPointer());
}
#endif

template<class TMultiPatchType>
typename TMultiPatchType::PatchContainerType MultiPatch_GetPatches(TMultiPatchType& rDummy)
{
    return rDummy.Patches();
}

template<class TMultiPatchType>
pybind11::list MultiPatch_GetPatchIndices(TMultiPatchType& rDummy)
{
    pybind11::list ids;
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
void IsogeometricApplication_AddPatchesToPython_Helper(pybind11::module& m)
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim>, typename Patch<TDim>::Pointer, Flags>
    (m, ss.str().c_str())
    .def(init<std::size_t, typename FESpace<TDim>::Pointer>())
    .def_property("Id", &Patch_GetId<Patch<TDim> >, &Patch_SetId<Patch<TDim> >)
    .def_property("Prefix", &Patch_GetPrefix<Patch<TDim> >, &Patch_SetPrefix<Patch<TDim> >)
    .def_property("LayerIndex", &Patch_GetLayerIndex<Patch<TDim> >, &Patch_SetLayerIndex<Patch<TDim> >)
    .def_property_readonly("WorkingSpaceDimension", &Patch<TDim>::WorkingSpaceDimension)
    .def_property_readonly("Name", &Patch<TDim>::Name)
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
    .def("LocalCoordinates", &Patch_LocalCoordinates<Patch<TDim> >)
    .def("IsInside", &Patch_IsInside<Patch<TDim> >)
    .def("NumberOfInterfaces", &Patch<TDim>::NumberOfInterfaces)
    .def("AddInterface", &Patch<TDim>::AddInterface)
    .def("RemoveInterface", &Patch<TDim>::RemoveInterface)
    .def("GetInterface", &Patch_GetInterface<TDim>)
    .def("ConstructBoundaryPatch", &Patch_ConstructBoundaryPatch<TDim>)
    .def("FindBoundarySide", &Patch<TDim>::FindBoundarySide)
    .def("__str__", &PrintObject<Patch<TDim> >)
    ;

    ss.str(std::string());
    ss << "Patch" << TDim << "DPointer";
    typedef iga::Wrapper<Patch<TDim>, typename Patch<TDim>::Pointer> PatchWrapperType;
    class_<PatchWrapperType>
    (m, ss.str().c_str())
    .def("GetReference", &PatchWrapperType::GetReference, return_value_policy::reference_internal)
    //.def("__str__", &PrintObject<PatchWrapperType>)
    ;

    ss.str(std::string());
    ss << "Patch" << TDim << "DContainer";
    PointerVectorSetPythonInterface<typename MultiPatch<TDim>::PatchContainerType> patch_container_python_interface;
    patch_container_python_interface.CreateInterface(m, ss.str());

    ss.str(std::string());
    ss << "PatchInterface" << TDim << "D";
    // class_<PatchInterface<TDim> >
    class_<PatchInterface<TDim>, typename PatchInterface<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<>())
    .def(init<typename Patch<TDim>::Pointer, const BoundarySide&, typename Patch<TDim>::Pointer, const BoundarySide&>())
    .def("Patch1", &PatchInterface_pPatch1<TDim>)
    .def("Patch2", &PatchInterface_pPatch2<TDim>)
    .def("Side1", &PatchInterface_Side1<PatchInterface<TDim> >)
    .def("Side2", &PatchInterface_Side2<PatchInterface<TDim> >)
    .def("__str__", &PrintObject<PatchInterface<TDim> >)
    ;

    ss.str(std::string());
    ss << "MultiPatch" << TDim << "D";
    class_<MultiPatch<TDim>, typename MultiPatch<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<>())
    // .def("ResetId", &MultiPatch<TDim>::ResetId) // this function is not really useful. One shall keep control over the id of the patch.
#ifdef SD_APP_FORWARD_COMPATIBILITY
    .def("AddPatch", &MultiPatch_AddPatch<TDim>)
    .def("RemovePatch", &MultiPatch_RemovePatch<TDim>)
#else
    .def("AddPatch", &MultiPatch<TDim>::AddPatch)
    .def("RemovePatch", &MultiPatch<TDim>::RemovePatch)
#endif
    .def("Patches", &MultiPatch_GetPatches<MultiPatch<TDim> >)
    .def("PatchIndices", &MultiPatch_GetPatchIndices<MultiPatch<TDim> >)
    .def("__getitem__", &MultiPatch_GetItem<Patch<TDim>, MultiPatch<TDim> >)
    .def("__len__", &MultiPatch_Len<MultiPatch<TDim> >)
    .def("EquationSystemSize", &MultiPatch<TDim>::EquationSystemSize)
    .def("ResetFunctionIndices", &MultiPatch<TDim>::ResetFunctionIndices)
    .def("Enumerate", &MultiPatch_Enumerate1<TDim>)
    .def("Enumerate", &MultiPatch_Enumerate2<TDim>)
    .def("IsEnumerated", &MultiPatch<TDim>::IsEnumerated)
    .def("__str__", &PrintObject<MultiPatch<TDim> >)
    ;
}

void IsogeometricApplication_AddExportToPython(pybind11::module& m)
{
    class_<MultiNURBSPatchGeoExporter, MultiNURBSPatchGeoExporter::Pointer>
    (m, "MultiNURBSPatchGeoExporter")
    .def(init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchGeoExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchGeoExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchGeoExporter, MultiPatch<3> >)
    .def("__str__", &PrintObject<MultiNURBSPatchGeoExporter>)
    ;

    class_<MultiNURBSPatchMatlabExporter, MultiNURBSPatchMatlabExporter::Pointer>
    (m, "MultiNURBSPatchMatlabExporter")
    .def(init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiNURBSPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiNURBSPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiNURBSPatchMatlabExporter, MultiPatch<3> >)
    .def("__str__", &PrintObject<MultiNURBSPatchMatlabExporter>)
    ;

    class_<MultiNURBSPatchGLVisExporter, MultiNURBSPatchGLVisExporter::Pointer>
    (m, "MultiNURBSPatchGLVisExporter")
    .def(init<>())
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
    .def("__str__", &PrintObject<MultiNURBSPatchGLVisExporter>)
    ;
}

template<int TDim>
void IsogeometricApplication_AddImportToPython(pybind11::module& m)
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "MultiNURBSPatchGeoImporter" << TDim << "D";
    class_<MultiNURBSPatchGeoImporter<TDim>, typename MultiNURBSPatchGeoImporter<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<>())
    .def("ImportSingle", &MultiNURBSPatchGeoImporter<TDim>::ImportSingle)
    .def("Import", &MultiNURBSPatchGeoImporter<TDim>::Import)
    .def("__str__", &PrintObject<MultiNURBSPatchGeoImporter<TDim> >)
    ;

}

void IsogeometricApplication_AddPatchesToPython(pybind11::module& m)
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////PATCHES///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    enum_<ParametricAxis>(m, "ParametricAxis")
    .value("U", _PA_U_)
    .value("V", _PA_V_)
    .value("W", _PA_W_)
    ;

    enum_<BoundarySide>(m, "BoundarySide")
    .value("Left", _BLEFT_)
    .value("Right", _BRIGHT_)
    .value("Top", _BTOP_)
    .value("Bottom", _BBOTTOM_)
    .value("Front", _BFRONT_)
    .value("Back", _BBACK_)
    ;

    enum_<BoundarySide2D>(m, "BoundarySide2D")
    .value("U0", _B2LEFT_)
    .value("U1", _B2RIGHT_)
    .value("V0", _B2BOTTOM_)
    .value("V1", _B2TOP_)
    ;

    enum_<BoundarySide3D>(m, "BoundarySide3D")
    .value("U0", _B3LEFT_)
    .value("U1", _B3RIGHT_)
    .value("V0", _B3FRONT_)
    .value("V1", _B3BACK_)
    .value("W0", _B3BOTTOM_)
    .value("W1", _B3TOP_)
    ;

    enum_<BoundaryDirection>(m, "BoundaryDirection")
    .value("Forward", _FORWARD_)
    .value("Reversed", _REVERSED_)
    ;

    enum_<BoundaryFlag>(m, "BoundaryFlag")
    .value("Left", _FLEFT_)
    .value("Right", _FRIGHT_)
    .value("Top", _FTOP_)
    .value("Bottom", _FBOTTOM_)
    .value("Front", _FFRONT_)
    .value("Back", _FBACK_)
    ;

    enum_<IsogeometricEchoFlags>(m, "IsogeometricEchoFlags")
    .value("ECHO_REFINEMENT", ECHO_REFINEMENT)
    .value("ECHO_REFINEMENT_DETAIL", ECHO_REFINEMENT_DETAIL)
    ;

    IsogeometricApplication_AddPatchesToPython_Helper<1>(m);
    IsogeometricApplication_AddPatchesToPython_Helper<2>(m);
    IsogeometricApplication_AddPatchesToPython_Helper<3>(m);

    /////////////////////////////////////////////////////////////////
    ///////////////////////IMPORT/EXPORT/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddExportToPython(m);
    IsogeometricApplication_AddImportToPython<1>(m);
    IsogeometricApplication_AddImportToPython<2>(m);
    IsogeometricApplication_AddImportToPython<3>(m);

}

}  // namespace Python.

} // Namespace Kratos
