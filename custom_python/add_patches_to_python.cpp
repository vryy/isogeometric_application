/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//
//

// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/python_utils.h"
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

template<class TPatchType>
typename TPatchType::Pointer Patch_Create(TPatchType& rDummy, std::size_t Id, typename TPatchType::FESpaceType::Pointer pFESpace)
{
    return rDummy.Create(Id, pFESpace);
}

template<class TPatchType, typename TDataType>
typename GridFunction<TPatchType::Dim, typename TPatchType::LocalCoordinateType, TDataType>::Pointer Patch_CreateGridFunction(TPatchType& rDummy, typename ControlGrid<TDataType>::Pointer pControlGrid)
{
    return rDummy.template CreateGridFunction<TDataType>(pControlGrid);
}

template<class TPatchType, class TVariableType>
typename GridFunction<TPatchType::Dim, typename TPatchType::LocalCoordinateType, typename TVariableType::Type>::Pointer Patch_GridFunction(TPatchType& rDummy, const TVariableType& rVariable)
{
    return rDummy.template pGetGridFunction<TVariableType>(rVariable);
}

template<class TPatchType>
boost::python::list Patch_Predict(TPatchType& rDummy, const boost::python::list& P, const boost::python::list& list_nsampling,
                                  const boost::python::list& list_xi_min, const boost::python::list& list_xi_max)
{
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;

    std::vector<LocalCoordinateType> xi_min_vec;
    PythonUtils::Unpack<double, LocalCoordinateType>(list_xi_min, xi_min_vec);

    std::vector<LocalCoordinateType> xi_max_vec;
    PythonUtils::Unpack<double, LocalCoordinateType>(list_xi_max, xi_max_vec);

    std::vector<CoordinateType> P_vec;
    PythonUtils::Unpack<CoordinateType, CoordinateType>(P, P_vec);

    array_1d<CoordinateType, 3> point;
    array_1d<LocalCoordinateType, 3> xi, xi_min, xi_max;

    point.clear();
    xi.clear();
    xi_min .clear();
    xi_max.clear();

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
    PythonUtils::Unpack<int, int>(list_nsampling, nsampling);

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
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;

    std::vector<CoordinateType> P_vec;
    PythonUtils::Unpack<CoordinateType, CoordinateType>(P, P_vec);

    array_1d<CoordinateType, 3> point;
    array_1d<LocalCoordinateType, 3> xi;

    point.clear();
    xi.clear();

    for (std::size_t i = 0; i < std::min(static_cast<std::size_t>(3), P_vec.size()); ++i)
    {
        point[i] = P_vec[i];
    }

    std::vector<int> nsampling;
    PythonUtils::Unpack<int, int>(list_nsampling, nsampling);

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
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;

    std::vector<LocalCoordinateType> xi0_vec;
    PythonUtils::Unpack<double, LocalCoordinateType>(xi0, xi0_vec);

    std::vector<CoordinateType> P_vec;
    PythonUtils::Unpack<CoordinateType, CoordinateType>(P, P_vec);

    array_1d<CoordinateType, 3> point;
    array_1d<LocalCoordinateType, 3> xi;

    point.clear();
    xi.clear();

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
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;

    std::vector<LocalCoordinateType> xi0_vec;
    PythonUtils::Unpack<double, LocalCoordinateType>(xi0, xi0_vec);

    std::vector<CoordinateType> P_vec;
    PythonUtils::Unpack<CoordinateType, CoordinateType>(P, P_vec);

    array_1d<CoordinateType, 3> point;
    array_1d<LocalCoordinateType, 3> xi;

    point.clear();
    xi.clear();

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

template<class TPatchType>
typename TPatchType::PatchInterfaceType::Pointer Patch_GetInterface(TPatchType& rDummy, std::size_t i)
{
    return rDummy.pInterface(i);
}

template<class TPatchType>
typename TPatchType::BoundaryPatchType::Pointer Patch_ConstructBoundaryPatch(TPatchType& rDummy, std::size_t iside)
{
    BoundarySide side = static_cast<BoundarySide>(iside);
    return rDummy.ConstructBoundaryPatch(side);
}

template<class TPatchType>
void Patch_save(TPatchType& rObject, Serializer& rSerializer, const std::string& rName)
{
    rSerializer.save(rName, rObject);
}

template<class TPatchType>
void Patch_load(TPatchType& rObject, Serializer& rSerializer, const std::string& rName)
{
    rSerializer.load(rName, rObject);
}

template<class TMultiPatchType>
static typename TMultiPatchType::Pointer MultiPatch_init(const boost::python::list& patch_list)
{
    typedef typename TMultiPatchType::PatchType PatchType;

    typename TMultiPatchType::Pointer pMultiPatch = typename TMultiPatchType::Pointer(new TMultiPatchType());

    std::vector<typename PatchType::Pointer> patches;
    PythonUtils::Unpack<typename PatchType::Pointer>(patch_list, patches);

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

template<class TMultiPatchType>
std::size_t MultiPatch_Enumerate1(TMultiPatchType& rDummy)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate();
//    KRATOS_WATCH(system_size)

    return system_size;
}

template<class TMultiPatchType>
std::size_t MultiPatch_Enumerate2(TMultiPatchType& rDummy, std::size_t start)
{
    std::size_t system_size;

    system_size = rDummy.Enumerate(start);
//    KRATOS_WATCH(system_size)

    return system_size;
}

template<class TPatchInterfaceType>
typename TPatchInterfaceType::PatchType::Pointer PatchInterface_pPatch1(TPatchInterfaceType& rDummy)
{
    return rDummy.pPatch1();
}

template<class TPatchInterfaceType>
BoundarySide PatchInterface_Side1(TPatchInterfaceType& rDummy)
{
    return rDummy.Side1();
}

template<class TPatchInterfaceType>
typename TPatchInterfaceType::PatchType::Pointer PatchInterface_pPatch2(TPatchInterfaceType& rDummy)
{
    return rDummy.pPatch2();
}

template<class TPatchInterfaceType>
BoundarySide PatchInterface_Side2(TPatchInterfaceType& rDummy)
{
    return rDummy.Side2();
}

////////////////////////////////////////

template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
void IsogeometricApplication_AddPatchesToPython_Helper(const std::string& Prefix)
{
    std::stringstream ss;

    typedef Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;
    typedef MultiPatch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> MultiPatchType;
    typedef PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchInterfaceType;

    typedef typename MatrixVectorTypeSelector<TDataType>::VectorType VectorType;

    ss.str(std::string());
    ss << Prefix << "Patch" << TDim << "D";
    auto patch_class = class_<PatchType, bases<Flags> >
    // class_<PatchType, typename PatchType::Pointer > // do not use this to export Patch pointer
    (ss.str().c_str(), init<std::size_t, typename FESpace<TDim>::Pointer>())
    .add_property("Id", &Patch_GetId<PatchType>, &Patch_SetId<PatchType>)
    .add_property("Prefix", &Patch_GetPrefix<PatchType>, &Patch_SetPrefix<PatchType>)
    .add_property("LayerIndex", &Patch_GetLayerIndex<PatchType>, &Patch_SetLayerIndex<PatchType>)
    .def("WorkingSpaceDimension", &PatchType::WorkingSpaceDimension)
    .def("Name", &PatchType::Name)
    .def("Create", &Patch_Create<PatchType>)
    .def("CreateControlPointGridFunction", &PatchType::CreateControlPointGridFunction)
    .def("CreateGridFunction", &Patch_CreateGridFunction<PatchType, TDataType>)
    .def("CreateGridFunction", &Patch_CreateGridFunction<PatchType, array_1d<TDataType, 3> >)
    .def("CreateGridFunction", &Patch_CreateGridFunction<PatchType, VectorType>)
    .def("GridFunction", &Patch_GridFunction<PatchType, Variable<ControlPoint<TCoordinateType> > >)
    .def("GridFunction", &Patch_GridFunction<PatchType, Variable<TDataType> >)
    .def("GridFunction", &Patch_GridFunction<PatchType, Variable<array_1d<TDataType, 3> > >)
    .def("GridFunction", &Patch_GridFunction<PatchType, Variable<VectorType> >)
    .def("ApplyTransformation", &PatchType::ApplyTransformation)
    .def("Order", &PatchType::Order)
    .def("TotalNumber", &PatchType::TotalNumber)
    .def("FESpace", &Patch_pFESpace<PatchType>)
    .def("Predict", &Patch_Predict<PatchType>)
    .def("Predict", &Patch_Predict1<PatchType>)
    .def("LocalCoordinates", &Patch_LocalCoordinates<PatchType>)
    .def("IsInside", &Patch_IsInside<PatchType>)
    .def("NumberOfInterfaces", &PatchType::NumberOfInterfaces)
    .def("AddInterface", &PatchType::AddInterface)
    .def("RemoveInterface", &PatchType::RemoveInterface)
    .def("GetInterface", &Patch_GetInterface<PatchType>)
    .def("ConstructBoundaryPatch", &Patch_ConstructBoundaryPatch<PatchType>)
    .def("ConstructSlicedPatch", &PatchType::ConstructSlicedPatch)
    .def("FindBoundarySide", &PatchType::FindBoundarySide)
    .def("SetLocalSearchTolerance", &PatchType::SetLocalSearchTolerance)
    .def("SetLocalSearchMaxIters", &PatchType::SetLocalSearchMaxIters)
    .def("Validate", &PatchType::Validate)
    .def("Save", &Patch_save<PatchType>)
    .def("Load", &Patch_load<PatchType>)
    .def(self_ns::str(self))
    ;

    if constexpr (!std::is_same<TCoordinateType, TDataType>::value)
    {
        patch_class.def("GridFunction", &Patch_GridFunction<PatchType, Variable<array_1d<TCoordinateType, 3> > >);
    }

    ss.str(std::string());
    ss << Prefix << "Patch" << TDim << "DPointer";
    class_<typename PatchType::Pointer>
    (ss.str().c_str(), init<typename PatchType::Pointer>())
    .def("GetReference", GetReference<PatchType>, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << Prefix << "Patch" << TDim << "DContainer";
    PointerVectorSetPythonInterface<typename MultiPatchType::PatchContainerType>::CreateInterface(ss.str());

    ss.str(std::string());
    ss << Prefix << "PatchInterface" << TDim << "D";
    // class_<PatchInterfaceType>
    class_<PatchInterfaceType, typename PatchInterfaceType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def(init<typename PatchType::Pointer, const BoundarySide&, typename PatchType::Pointer, const BoundarySide&>())
    .def("Patch1", &PatchInterface_pPatch1<PatchInterfaceType>)
    .def("Patch2", &PatchInterface_pPatch2<PatchInterfaceType>)
    .def("Side1", &PatchInterface_Side1<PatchInterfaceType>)
    .def("Side2", &PatchInterface_Side2<PatchInterfaceType>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << Prefix << "MultiPatch" << TDim << "D";
    class_<MultiPatchType, typename MultiPatchType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__init__", make_constructor(&MultiPatch_init<MultiPatchType>))
    // .def("ResetId", &MultiPatchType::ResetId) // this function is not really useful. One shall keep control over the id of the patch.
    .def("AddPatch", &MultiPatchType::AddPatch)
    .def("RemovePatch", &MultiPatchType::RemovePatch)
    .def("Patches", &MultiPatch_GetPatches<MultiPatchType>)
    .def("PatchIndices", &MultiPatch_GetPatchIndices<MultiPatchType>)
    .def("__getitem__", &MultiPatch_GetItem<PatchType, MultiPatchType>)
    .def("__len__", &MultiPatch_Len<MultiPatchType>)
    .def("EquationSystemSize", &MultiPatchType::EquationSystemSize)
    .def("ResetFunctionIndices", &MultiPatchType::ResetFunctionIndices)
    .def("Enumerate", &MultiPatch_Enumerate1<MultiPatchType>)
    .def("Enumerate", &MultiPatch_Enumerate2<MultiPatchType>)
    .def("IsEnumerated", &MultiPatchType::IsEnumerated)
    .def("LocalCoordinates", &Patch_LocalCoordinates<MultiPatchType>)
    .def("Validate", &MultiPatchType::Validate)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddPatchSelectorToPython_Helper()
{
    std::stringstream ss;
    ss << "Patch" << TDim << "DSelectorClass";
    class_<PatchSelector<TDim>, boost::noncopyable>
    (ss.str().c_str(), no_init)
    .add_property("RealPatch", make_function(&PatchSelector<TDim>::GetRealPatch, return_value_policy<reference_existing_object>()))
    .add_property("ComplexPatch", make_function(&PatchSelector<TDim>::GetComplexPatch, return_value_policy<reference_existing_object>()))
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

    // export the patch selector to select patch at runtime
    IsogeometricApplication_AddPatchSelectorToPython_Helper<1>();
    IsogeometricApplication_AddPatchSelectorToPython_Helper<2>();
    IsogeometricApplication_AddPatchSelectorToPython_Helper<3>();

    scope().attr("Patch1DSelector") = boost::ref(PatchSelector1DInstance);
    scope().attr("Patch2DSelector") = boost::ref(PatchSelector2DInstance);
    scope().attr("Patch3DSelector") = boost::ref(PatchSelector3DInstance);
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

    IsogeometricApplication_AddPatchesToPython_Helper<1, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE>("");
    IsogeometricApplication_AddPatchesToPython_Helper<2, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE>("");
    IsogeometricApplication_AddPatchesToPython_Helper<3, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE>("");
    IsogeometricApplication_AddPatchesToPython_Helper<1, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE>("Complex");
    IsogeometricApplication_AddPatchesToPython_Helper<2, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE>("Complex");
    IsogeometricApplication_AddPatchesToPython_Helper<3, KRATOS_DOUBLE_TYPE, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE>("Complex");
    IsogeometricApplication_AddPatchesToPython_Helper<1, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE>("GComplex");
    IsogeometricApplication_AddPatchesToPython_Helper<2, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE>("GComplex");
    IsogeometricApplication_AddPatchesToPython_Helper<3, KRATOS_DOUBLE_TYPE, KRATOS_COMPLEX_TYPE, KRATOS_COMPLEX_TYPE>("GComplex");

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
