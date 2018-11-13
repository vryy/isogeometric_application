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
#include "python/pointer_vector_set_python_interface.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/tsplines/tcell.h"
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/nurbs/pbbsplines_fespace.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_python/add_point_based_control_grid_to_python.h"
#include "custom_python/add_import_export_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<int TDim>
void IsogeometricApplication_AddPBBSplinesSpaceToPython()
{

    std::stringstream ss;

    ss.str(std::string());
    ss << "PBBSplinesBasisFunction" << TDim << "D";
    typedef PBBSplinesBasisFunction<TDim, TCell> PBBSplinesBasisFunctionType;
    class_<PBBSplinesBasisFunctionType, typename PBBSplinesBasisFunctionType::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<const std::size_t&>())
    .add_property("Id", Isogeometric_GetId<PBBSplinesBasisFunctionType>, Isogeometric_DoNotSetId<PBBSplinesBasisFunctionType>)
    .add_property("EquationId", Isogeometric_GetEquationId<PBBSplinesBasisFunctionType>, Isogeometric_SetEquationId<PBBSplinesBasisFunctionType>)
    .def("Weight", &PBBSplinesBasisFunctionType::Weight)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "PBBSplinesFESpace" << TDim << "D";
    typedef FESpace<TDim> FESpaceType;
    typedef PBBSplinesFESpace<TDim, PBBSplinesBasisFunctionType, BCellManager<TDim, typename PBBSplinesBasisFunctionType::CellType> > PBBSplinesFESpaceType;
    class_<PBBSplinesFESpaceType, typename PBBSplinesFESpaceType::Pointer, bases<FESpaceType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("__getitem__", &FESpace_GetItem<PBBSplinesFESpaceType>)
    .def("UpdateCells", &PBBSplinesFESpaceType::UpdateCells)
    .def(self_ns::str(self))
    ;

    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<double>, PBBSplinesFESpaceType>::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<array_1d<double, 3> >, PBBSplinesFESpaceType>::Execute();
    IsogeometricApplication_AddPointBasedControlGrid_Helper<Variable<Vector>, PBBSplinesFESpaceType>::Execute();

}

////////////////////////////////////////

void IsogeometricApplication_AddPBBSplinesToPython()
{

    /////////////////////////////////////////////////////////////////
    ///////////////////////Point-based BSplines//////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddPBBSplinesSpaceToPython<1>();
    IsogeometricApplication_AddPBBSplinesSpaceToPython<2>();
    IsogeometricApplication_AddPBBSplinesSpaceToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos

