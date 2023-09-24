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

// Project includes
#include "includes/define.h"
#include "python/containers_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/tsplines/tcell.h"
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/nurbs/pbbsplines_fespace.h"
#include "custom_utilities/import_export/multi_pbbsplines_patch_matlab_exporter.h"
#include "custom_python/iga_define_python.h"
#include "custom_python3/add_pbbsplines_to_python.h"
#include "custom_python3/add_point_based_control_grid_to_python.h"
#include "custom_python3/add_import_export_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

template<int TDim>
void IsogeometricApplication_AddPBBSplinesSpaceToPython(pybind11::module& m)
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "PBBSplinesBasisFunction" << TDim << "D";
    typedef PBBSplinesBasisFunction<TDim, TCell> PBBSplinesBasisFunctionType;
    class_<PBBSplinesBasisFunctionType, typename PBBSplinesBasisFunctionType::Pointer>
    (m, ss.str().c_str())
    .def(init<std::size_t>())
    .def_property("Id", Isogeometric_GetId<PBBSplinesBasisFunctionType>, Isogeometric_DoNotSetId<PBBSplinesBasisFunctionType>)
    .def_property("EquationId", Isogeometric_GetEquationId<PBBSplinesBasisFunctionType>, Isogeometric_SetEquationId<PBBSplinesBasisFunctionType>)
    .def_property_readonly("Weight", &PBBSplinesBasisFunctionType::Weight)
    .def("__str__", &PrintObject<PBBSplinesBasisFunctionType>)
    ;

    ss.str(std::string());
    ss << "PBBSplinesFESpace" << TDim << "D";
    typedef FESpace<TDim> FESpaceType;
    typedef PBBSplinesFESpace<TDim, PBBSplinesBasisFunctionType, BCellManager<TDim, typename PBBSplinesBasisFunctionType::CellType> > PBBSplinesFESpaceType;
    class_<PBBSplinesFESpaceType, typename PBBSplinesFESpaceType::Pointer, FESpaceType>
    (m, ss.str().c_str())
    .def(init<>())
    .def("__getitem__", &FESpace_GetItem<PBBSplinesFESpaceType>)
    .def("UpdateCells", &PBBSplinesFESpaceType::UpdateCells)
    .def("__str__", &PrintObject<PBBSplinesFESpaceType>)
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, PBBSplinesFESpaceType>::Execute(m);
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, PBBSplinesFESpaceType>::Execute(m);
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, PBBSplinesFESpaceType>::Execute(m);

}

////////////////////////////////////////

void IsogeometricApplication_AddPBBSplinesToPython(pybind11::module& m)
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////Point-based BSplines//////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddPBBSplinesSpaceToPython<1>(m);
    IsogeometricApplication_AddPBBSplinesSpaceToPython<2>(m);
    IsogeometricApplication_AddPBBSplinesSpaceToPython<3>(m);

    class_<MultiPBBSplinesPatchMatlabExporter, MultiPBBSplinesPatchMatlabExporter::Pointer>
    (m, "MultiPBBSplinesPatchMatlabExporter")
    .def(init<>())
    .def("Export", &MultiPatchExporter_Export<1, MultiPBBSplinesPatchMatlabExporter, Patch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiPBBSplinesPatchMatlabExporter, Patch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiPBBSplinesPatchMatlabExporter, Patch<3> >)
    .def("Export", &MultiPatchExporter_Export<1, MultiPBBSplinesPatchMatlabExporter, MultiPatch<1> >)
    .def("Export", &MultiPatchExporter_Export<2, MultiPBBSplinesPatchMatlabExporter, MultiPatch<2> >)
    .def("Export", &MultiPatchExporter_Export<3, MultiPBBSplinesPatchMatlabExporter, MultiPatch<3> >)
    .def("__str__", &PrintObject<MultiPBBSplinesPatchMatlabExporter>)
    ;

}

}  // namespace Python.

} // Namespace Kratos

