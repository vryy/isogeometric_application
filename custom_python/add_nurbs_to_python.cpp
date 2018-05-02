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
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/trans/translation.h"
#include "custom_utilities/trans/rotation.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/unstructured_control_grid.h"
#include "custom_utilities/nurbs/structured_control_grid.h"
#include "custom_utilities/control_grid_library.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/weighted_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/bending_strip_patch.h"
#include "custom_utilities/nurbs/bending_strip_nurbs_patch.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/import_export/multi_nurbs_patch_matlab_exporter.h"
#include "custom_utilities/import_export/multi_nurbs_patch_glvis_exporter.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

double ControlPoint_GetWX(ControlPoint<double>& rDummy)
{
    return rDummy.WX();
}

void ControlPoint_SetWX(ControlPoint<double>& rDummy, const double& newWX)
{
    rDummy.WX() = newWX;
}

double ControlPoint_GetWY(ControlPoint<double>& rDummy)
{
    return rDummy.WY();
}

void ControlPoint_SetWY(ControlPoint<double>& rDummy, const double& newWY)
{
    rDummy.WY() = newWY;
}

double ControlPoint_GetWZ(ControlPoint<double>& rDummy)
{
    return rDummy.WZ();
}

void ControlPoint_SetWZ(ControlPoint<double>& rDummy, const double& newWZ)
{
    rDummy.WZ() = newWZ;
}

double ControlPoint_GetW(ControlPoint<double>& rDummy)
{
    return rDummy.W();
}

void ControlPoint_SetW(ControlPoint<double>& rDummy, const double& newW)
{
    rDummy.W() = newW;
}

void ControlPoint_ApplyTransformation(ControlPoint<double>& rDummy, const Transformation<double>& trans)
{
    rDummy.ApplyTransformation(trans);
}

////////////////////////////////////////

template<typename TDataType>
TDataType ControlGrid_GetItem(ControlGrid<TDataType>& rDummy, int index)
{
    return rDummy.GetData(index);
}

template<typename TDataType>
void ControlGrid_SetItem(ControlGrid<TDataType>& rDummy, int index, const TDataType& value)
{
    rDummy.SetData(index, value);
}

template<typename TDataType>
struct ControlValue_Helper
{
    static boost::python::list GetValue(const TDataType& rValue)
    {
    }

    static TDataType GetValue(boost::python::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<double>
{
    static boost::python::list GetValue(const double& rValue)
    {
        boost::python::list v;
        v.append(rValue);
        return v;
    }

    static double GetValue(boost::python::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<ControlPoint<double> >
{
    static boost::python::list GetValue(const ControlPoint<double>& rValue)
    {
        boost::python::list v;
        v.append(rValue.WX());
        v.append(rValue.WY());
        v.append(rValue.WZ());
        v.append(rValue.W());
        return v;
    }

    static double GetValue(boost::python::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<array_1d<double, 3> >
{
    static boost::python::list GetValue(const array_1d<double, 3>& rValue)
    {
        boost::python::list v;
        v.append(rValue[0]);
        v.append(rValue[1]);
        v.append(rValue[2]);
        return v;
    }

    static double GetValue(boost::python::list rValue)
    {
    }
};

template<int TDim, typename TDataType>
struct StructuredControlGrid_Helper
{
    static boost::python::list GetValue(StructuredControlGrid<TDim, TDataType>& rDummy)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not implemented")
    }

    static void SetValue(StructuredControlGrid<TDim, TDataType>& rDummy, boost::python::list values)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not implemented")
    }
};

template<typename TDataType>
struct StructuredControlGrid_Helper<1, TDataType>
{
    static boost::python::list GetValue(StructuredControlGrid<1, TDataType>& rDummy)
    {
        boost::python::list output;

        for (std::size_t i = 0; i < rDummy.size(); ++i)
        {
            boost::python::list v = ControlValue_Helper<TDataType>::GetValue(rDummy.GetValue(i));
            output.append(v);
        }

        return output;
    }

    static void SetValue(StructuredControlGrid<1, TDataType>& rDummy, boost::python::list values)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not implemented")
    }
};

template<typename TDataType>
struct StructuredControlGrid_Helper<2, TDataType>
{
    static boost::python::list GetValue(StructuredControlGrid<2, TDataType>& rDummy)
    {
        boost::python::list output;

        for (std::size_t j = 0; j < rDummy.Size(1); ++j)
        {
            boost::python::list row;
            for (std::size_t i = 0; i < rDummy.Size(0); ++i)
            {
                boost::python::list v = ControlValue_Helper<TDataType>::GetValue(rDummy.GetValue(i, j));
                row.append(v);
            }
            output.append(row);
        }

        return output;
    }

    static void SetValue(StructuredControlGrid<2, TDataType>& rDummy, boost::python::list values)
    {
    }
};

template<typename TDataType>
struct StructuredControlGrid_Helper<3, TDataType>
{
    static boost::python::list GetValue(StructuredControlGrid<3, TDataType>& rDummy)
    {
        boost::python::list output;

        for (std::size_t k = 0; k < rDummy.Size(2); ++k)
        {
            boost::python::list row;
            for (std::size_t j = 0; j < rDummy.Size(1); ++j)
            {
                boost::python::list col;
                for (std::size_t i = 0; i < rDummy.Size(0); ++i)
                {
                    boost::python::list v = ControlValue_Helper<TDataType>::GetValue(rDummy.GetValue(i, j, k));
                    col.append(v);
                }
                row.append(col);
            }
            output.append(row);
        }

        return output;
    }

    static void SetValue(StructuredControlGrid<3, TDataType>& rDummy, boost::python::list values)
    {
    }
};

////////////////////////////////////////

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateLinearControlPointGrid(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u,
        const double& end_x, const double& end_y, const double& end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateStructuredControlPointGrid<1>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid1(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& end_x, const double& end_y)
{
    std::vector<double> start(2);
    start[0] = start_x;
    start[1] = start_y;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> end(2);
    end[0] = end_x;
    end[1] = end_y;

    return rDummy.CreateStructuredControlPointGrid<2>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid2(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& space1_x, const double& space1_y, const double& space1_z,
        const double& space2_x, const double& space2_y, const double& space2_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;

    std::vector<double> space1(3);
    space1[0] = space1_x;
    space1[1] = space1_y;
    space1[2] = space1_z;

    std::vector<double> space2(3);
    space2[0] = space2_x;
    space2[1] = space2_y;
    space2[2] = space2_z;

    std::vector<std::vector<double> > spacing_vectors(2);
    spacing_vectors[0] = space1;
    spacing_vectors[1] = space2;

    return rDummy.CreateStructuredControlPointGrid<2>(start, ngrid, spacing_vectors);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid1(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        const double& end_x, const double& end_y, const double& end_z)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<double> end(3);
    end[0] = end_x;
    end[1] = end_y;
    end[2] = end_z;

    return rDummy.CreateStructuredControlPointGrid<3>(start, ngrid, end);
}

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid2(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        boost::python::list spacing_vectors_data)
{
    std::vector<double> start(3);
    start[0] = start_x;
    start[1] = start_y;
    start[2] = start_z;

    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;

    std::vector<std::vector<double> > spacing_vectors;
    std::size_t cnt1 = 0, cnt2 = 0;
    typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& vect, std::make_pair(iterator_value_type(spacing_vectors_data), iterator_value_type() ) )
    {
        typedef boost::python::stl_input_iterator<double> iterator_value_type2;
        std::vector<double> space_vect;
        BOOST_FOREACH(const iterator_value_type2::value_type& v, std::make_pair(iterator_value_type2(vect), iterator_value_type2() ) )
        {
            space_vect.push_back(v);
        }
        spacing_vectors.push_back(space_vect);
    }

    return rDummy.CreateStructuredControlPointGrid<3>(start, ngrid, spacing_vectors);
}

template<class TVariableType>
typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateLinearZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u)
{
    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;
    return rDummy.CreateStructuredZeroControlGrid<1, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateRectangularZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u, const std::size_t& n_points_v)
{
    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    return rDummy.CreateStructuredZeroControlGrid<2, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateCubicZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w)
{
    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;
    return rDummy.CreateStructuredZeroControlGrid<3, TVariableType>(rVariable, ngrid);
}

////////////////////////////////////////

template<int TDim>
std::size_t FESpace_Enumerate(FESpace<TDim>& rDummy)
{
    std::size_t eq_size = 0;
    eq_size = rDummy.Enumerate(eq_size);
    return eq_size;
}

template<int TDim>
boost::python::list FESpace_BoundaryFunctionIndices(FESpace<TDim>& rDummy, const BoundarySide& side)
{
    boost::python::list indices;
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i]));
    return indices;
}

template<int TDim>
boost::python::list FESpace_BoundaryShiftedFunctionIndices(FESpace<TDim>& rDummy, const BoundarySide& side)
{
    boost::python::list indices;
    std::vector<std::size_t> boundary_indices = rDummy.ExtractBoundaryFunctionIndices(side);
    for (std::size_t i = 0; i < boundary_indices.size(); ++i)
        indices.append<int>(static_cast<int>(boundary_indices[i] + 1));
    return indices;
}

////////////////////////////////////////

template<int TDim, int TWhichDim>
boost::python::list BSplinesFESpace_GetKnotVector(BSplinesFESpace<TDim>& rDummy)
{
    boost::python::list knot_list;

    if (TWhichDim < TDim)
    {
        const typename BSplinesFESpace<TDim>::knot_container_t& knot_vector = rDummy.KnotVector(TWhichDim);

        for (std::size_t i = 0; i < knot_vector.size(); ++i)
            knot_list.append(knot_vector[i]);
    }

    return knot_list;
}

template<int TDim, int TWhichDim>
void BSplinesFESpace_SetKnotVector(BSplinesFESpace<TDim>& rDummy, boost::python::list knot_list)
{
    if (TWhichDim < TDim)
    {
        std::vector<double> knot_vec;
        typedef boost::python::stl_input_iterator<double> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& v, std::make_pair(iterator_value_type(knot_list), iterator_value_type() ) )
        {
            knot_vec.push_back(v);
        }

        rDummy.SetKnotVector(TWhichDim, knot_vec);
    }
}

////////////////////////////////////////

BSplinesFESpace<1>::Pointer BSplinesFESpaceLibrary_CreateLinearFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u)
{
    std::vector<std::size_t> orders(1);
    orders[0] = order_u;
    return rDummy.CreateRegularFESpace<1>(orders);
}

BSplinesFESpace<2>::Pointer BSplinesFESpaceLibrary_CreateRectangularFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v)
{
    std::vector<std::size_t> orders(2);
    orders[0] = order_u;
    orders[1] = order_v;
    return rDummy.CreateRegularFESpace<2>(orders);
}

BSplinesFESpace<3>::Pointer BSplinesFESpaceLibrary_CreateCubicFESpace(BSplinesFESpaceLibrary& rDummy, const std::size_t& order_u, const std::size_t& order_v, const std::size_t& order_w)
{
    std::vector<std::size_t> orders(3);
    orders[0] = order_u;
    orders[1] = order_v;
    orders[2] = order_w;
    return rDummy.CreateRegularFESpace<3>(orders);
}

////////////////////////////////////////

template<int TDim, typename TDataType>
typename FESpace<TDim>::Pointer GridFunction_GetFESpace(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pFESpace();
}

template<int TDim, typename TDataType>
void GridFunction_SetFESpace(GridFunction<TDim, TDataType>& rDummy, typename FESpace<TDim>::Pointer pNewFESpace)
{
    rDummy.SetFESpace(pNewFESpace);
}

template<int TDim, typename TDataType>
typename ControlGrid<TDataType>::Pointer GridFunction_GetControlGrid(GridFunction<TDim, TDataType>& rDummy)
{
    return rDummy.pControlGrid();
}

template<int TDim, typename TDataType>
void GridFunction_SetControlGrid(GridFunction<TDim, TDataType>& rDummy, typename ControlGrid<TDataType>::Pointer pNewControlGrid)
{
    rDummy.SetControlGrid(pNewControlGrid);
}

template<int TDim, typename TDataType>
TDataType GridFunction_GetValue(GridFunction<TDim, TDataType>& rDummy, const boost::python::list& xi)
{
    std::vector<double> xi_vec;
    typedef boost::python::stl_input_iterator<double> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& v, std::make_pair(iterator_value_type(xi), iterator_value_type() ) )
    {
        xi_vec.push_back(v);
    }

    return rDummy.GetValue(xi_vec);
}

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
typename TPatchType::Pointer Patch_pGetNeighbor(TPatchType& rDummy, const BoundarySide& side)
{
    return rDummy.pNeighbor(side);
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
void MultiPatch_MakeNeighbor(MultiPatch<TDim>& rDummy, typename Patch<TDim>::Pointer pPatch1, BoundarySide side1,
           typename Patch<TDim>::Pointer pPatch2, BoundarySide side2)
{
   rDummy.MakeNeighbor(pPatch1, side1, pPatch2, side2);
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

//////////////////////////////////////////////////

template<int TDim, class TExporter, class TPatchType>
void MultiPatchExporter_Export(TExporter& rDummy,
        typename TPatchType::Pointer pPatch, const std::string& filename)
{
    rDummy.template Export<TDim>(pPatch, filename);
}

//////////////////////////////////////////////////

void DomainManager2D_AddCell(DomainManager2D& rDummy, const double& x1, const double& x2, const double& y1, const double& y2)
{
    std::vector<double> box(4);
    box[0] = x1;
    box[1] = x2;
    box[2] = y1;
    box[3] = y2;

    rDummy.AddCell(box);
}

bool DomainManager2D_IsInside(DomainManager2D& rDummy, const double& x1, const double& x2, const double& y1, const double& y2)
{
    std::vector<double> box(4);
    box[0] = x1;
    box[1] = x2;
    box[2] = y1;
    box[3] = y2;

    return rDummy.IsInside(box);
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void IsogeometricApplication_AddTransformation()
{
    class_<Transformation<double>, Transformation<double>::Pointer, boost::noncopyable>
    ("Transformation", init<>())
    .def("AppendTransformation", &Transformation<double>::AppendTransformation)
    // .def(boost::python::operators<boost::python::op_mul>());
    .def(self_ns::str(self))
    ;

    class_<Translation<double>, Translation<double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("Translation", init<const double&, const double&, const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<0, double>, Rotation<0, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationX", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<1, double>, Rotation<1, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationY", init<const double&>())
    .def(self_ns::str(self))
    ;

    class_<Rotation<2, double>, Rotation<2, double>::Pointer, bases<Transformation<double> >, boost::noncopyable>
    ("RotationZ", init<const double&>())
    .def(self_ns::str(self))
    ;
}

void IsogeometricApplication_AddControlPoint()
{
    class_<ControlPoint<double>, ControlPoint<double>::Pointer>
    ("ControlPoint", init<>())
    .add_property("WX", ControlPoint_GetWX, ControlPoint_SetWX)
    .add_property("WY", ControlPoint_GetWY, ControlPoint_SetWY)
    .add_property("WZ", ControlPoint_GetWZ, ControlPoint_SetWZ)
    .add_property("W", ControlPoint_GetW, ControlPoint_SetW)
    .def("ApplyTransformation", &ControlPoint_ApplyTransformation)
    .def(self_ns::str(self))
    ;

    class_<Variable<ControlPoint<double> >, bases<VariableData>, boost::noncopyable >( "ControlPointVariable", no_init )
    .def( self_ns::str( self ) )
    ;
}

void IsogeometricApplication_AddControlGrids()
{
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<ControlGrid<ControlPoint<double> >, ControlGrid<ControlPoint<double> >::Pointer, boost::noncopyable>
    ("ControlPointControlGrid", init<>())
    .def("Size", &ControlGrid<ControlPoint<double> >::Size)
    .def("size", &ControlGrid<ControlPoint<double> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<ControlPoint<double> >)
    .def("__getitem__", &ControlGrid_GetItem<ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<double>, ControlGrid<double>::Pointer, boost::noncopyable>
    ("DoubleControlGrid", init<>())
    .def("Size", &ControlGrid<double>::Size)
    .def("size", &ControlGrid<double>::Size)
    .def("__setitem__", &ControlGrid_SetItem<double>)
    .def("__getitem__", &ControlGrid_GetItem<double>)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<array_1d<double, 3> >, ControlGrid<array_1d<double, 3> >::Pointer, boost::noncopyable>
    ("Array1DControlGrid", init<>())
    .def("Size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("size", &ControlGrid<array_1d<double, 3> >::Size)
    .def("__setitem__", &ControlGrid_SetItem<array_1d<double, 3> >)
    .def("__getitem__", &ControlGrid_GetItem<array_1d<double, 3> >)
    .def(self_ns::str(self))
    ;

    class_<ControlGrid<Vector>, ControlGrid<Vector>::Pointer, boost::noncopyable>
    ("VectorControlGrid", init<>())
    .def("Size", &ControlGrid<Vector>::Size)
    .def("size", &ControlGrid<Vector>::Size)
    .def("__setitem__", &ControlGrid_SetItem<Vector>)
    .def("__getitem__", &ControlGrid_GetItem<Vector>)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<BaseStructuredControlGrid<ControlPoint<double> >, BaseStructuredControlGrid<ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("BaseStructuredControlPointGrid", init<>())
    .def(self_ns::str(self))
    ;

    class_<BaseStructuredControlGrid<double>, BaseStructuredControlGrid<double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("BaseStructuredDoubleControlGrid", init<>())
    .def(self_ns::str(self))
    ;

    class_<BaseStructuredControlGrid<array_1d<double, 3> >, BaseStructuredControlGrid<array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("BaseStructuredArray1DControlGrid", init<>())
    .def(self_ns::str(self))
    ;

    class_<BaseStructuredControlGrid<Vector>, BaseStructuredControlGrid<Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("BaseStructuredVectorControlGrid", init<>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<StructuredControlGrid<1, ControlPoint<double> >, StructuredControlGrid<1, ControlPoint<double> >::Pointer, bases<BaseStructuredControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("StructuredControlPointGrid1D", init<const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<1, ControlPoint<double> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, double>, StructuredControlGrid<1, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid1D", init<const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, double>::GetValue, &StructuredControlGrid_Helper<1, double>::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, array_1d<double, 3> >, StructuredControlGrid<1, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid1D", init<const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<1, array_1d<double, 3> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, Vector>, StructuredControlGrid<1, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid1D", init<const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, Vector>::GetValue, &StructuredControlGrid_Helper<1, Vector>::SetValue)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<StructuredControlGrid<2, ControlPoint<double> >, StructuredControlGrid<2, ControlPoint<double> >::Pointer, bases<BaseStructuredControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("StructuredControlPointGrid2D", init<const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<2, ControlPoint<double> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, double>, StructuredControlGrid<2, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, double>::GetValue, &StructuredControlGrid_Helper<2, double>::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, array_1d<double, 3> >, StructuredControlGrid<2, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<2, array_1d<double, 3> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, Vector>, StructuredControlGrid<2, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid2D", init<const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, Vector>::GetValue, &StructuredControlGrid_Helper<2, Vector>::SetValue)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<StructuredControlGrid<3, ControlPoint<double> >, StructuredControlGrid<3, ControlPoint<double> >::Pointer, bases<BaseStructuredControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("StructuredControlPointGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<3, ControlPoint<double> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, double>, StructuredControlGrid<3, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, double>::GetValue, &StructuredControlGrid_Helper<3, double>::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, array_1d<double, 3> >, StructuredControlGrid<3, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<1, array_1d<double, 3> >::SetValue)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, Vector>, StructuredControlGrid<3, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid3D", init<const std::size_t&, const std::size_t&, const std::size_t&>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, Vector>::GetValue, &StructuredControlGrid_Helper<3, Vector>::SetValue)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<UnstructuredControlGrid<ControlPoint<double> >, UnstructuredControlGrid<ControlPoint<double> >::Pointer, bases<ControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("UnstructuredControlPointGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<double>, UnstructuredControlGrid<double>::Pointer, bases<ControlGrid<double> >, boost::noncopyable>
    ("UnstructuredDoubleControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<array_1d<double, 3> >, UnstructuredControlGrid<array_1d<double, 3> >::Pointer, bases<ControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("UnstructuredArray1DControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    class_<UnstructuredControlGrid<Vector>, UnstructuredControlGrid<Vector>::Pointer, bases<ControlGrid<Vector> >, boost::noncopyable>
    ("UnstructuredVectorControlGrid", init<const std::size_t&>())
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////
}

template<int TDim>
void IsogeometricApplication_AddFESpacesToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "FESpace" << TDim << "D";
    class_<FESpace<TDim>, typename FESpace<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Order", &FESpace<TDim>::Order)
    .def("TotalNumber", &FESpace<TDim>::TotalNumber)
    .def("Enumerate", &FESpace_Enumerate<TDim>)
    .def("BoundaryFunctionIndices", &FESpace_BoundaryFunctionIndices<TDim>)
    .def("BoundaryShiftedFunctionIndices", &FESpace_BoundaryShiftedFunctionIndices<TDim>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "WeightedFESpace" << TDim << "D";
    class_<WeightedFESpace<TDim>, typename WeightedFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, const std::vector<double>&>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "BSplinesFESpace" << TDim << "D";
    class_<BSplinesFESpace<TDim>, typename BSplinesFESpace<TDim>::Pointer, bases<FESpace<TDim> >, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("Number", &BSplinesFESpace<TDim>::Number)
    .add_property("KnotU", BSplinesFESpace_GetKnotVector<TDim, 0>, BSplinesFESpace_SetKnotVector<TDim, 0>)
    .add_property("KnotV", BSplinesFESpace_GetKnotVector<TDim, 1>, BSplinesFESpace_SetKnotVector<TDim, 1>)
    .add_property("KnotW", BSplinesFESpace_GetKnotVector<TDim, 2>, BSplinesFESpace_SetKnotVector<TDim, 2>)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddGridFunctionsToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "ControlPointGridFunction" << TDim << "D";
    class_<GridFunction<TDim, ControlPoint<double> >, typename GridFunction<TDim, ControlPoint<double> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<ControlPoint<double> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, ControlPoint<double> >, GridFunction_SetFESpace<TDim, ControlPoint<double> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, ControlPoint<double> >, GridFunction_SetControlGrid<TDim, ControlPoint<double> >)
    .def("GetValue", &GridFunction_GetValue<TDim, ControlPoint<double> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "DoubleGridFunction" << TDim << "D";
    class_<GridFunction<TDim, double>, typename GridFunction<TDim, double>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<double>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, double>, GridFunction_SetFESpace<TDim, double>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, double>, GridFunction_SetControlGrid<TDim, double>)
    .def("GetValue", &GridFunction_GetValue<TDim, double>)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "Array1DGridFunction" << TDim << "D";
    class_<GridFunction<TDim, array_1d<double, 3> >, typename GridFunction<TDim, array_1d<double, 3> >::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<array_1d<double, 3> >::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, array_1d<double, 3> >, GridFunction_SetFESpace<TDim, array_1d<double, 3> >)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, array_1d<double, 3> >, GridFunction_SetControlGrid<TDim, array_1d<double, 3> >)
    .def("GetValue", &GridFunction_GetValue<TDim, array_1d<double, 3> >)
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "VectorGridFunction" << TDim << "D";
    class_<GridFunction<TDim, Vector>, typename GridFunction<TDim, Vector>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<typename FESpace<TDim>::Pointer, typename ControlGrid<Vector>::Pointer>())
    .add_property("FESpace", GridFunction_GetFESpace<TDim, Vector>, GridFunction_SetFESpace<TDim, Vector>)
    .add_property("ControlGrid", GridFunction_GetControlGrid<TDim, Vector>, GridFunction_SetControlGrid<TDim, Vector>)
    .def("GetValue", &GridFunction_GetValue<TDim, Vector>)
    .def(self_ns::str(self))
    ;
}

template<int TDim>
void IsogeometricApplication_AddPatchesToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "Patch" << TDim << "D";
    class_<Patch<TDim> >
    // class_<Patch<TDim>, typename Patch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&, typename FESpace<TDim>::Pointer>())
    .add_property("Id", &Patch_GetId<Patch<TDim> >, &Patch_SetId<Patch<TDim> >)
    .def("WorkingSpaceDimension", &Patch<TDim>::WorkingSpaceDimension)
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
    .def("Neighbor", &Patch_pGetNeighbor<Patch<TDim> >)
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
    PointerVectorSetPythonInterface<typename MultiPatch<TDim>::PatchContainerType>::CreateInterface(ss.str())
    ;

    ss.str(std::string());
    ss << "BendingStripPatch" << TDim << "D";
    class_<BendingStripPatch<TDim>, bases<Patch<TDim> > >
    // class_<BendingStripPatch<TDim>, typename BendingStripPatch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&, const int&>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "BendingStripNURBSPatch" << TDim << "D";
    class_<BendingStripNURBSPatch<TDim>, bases<BendingStripPatch<TDim> > >
    // class_<BendingStripNURBSPatch<TDim>, typename BendingStripNURBSPatch<TDim>::Pointer >
    (ss.str().c_str(), init<const std::size_t&, const int&>())
    .def(init<const std::size_t&, typename Patch<TDim>::Pointer, const BoundarySide&, typename Patch<TDim>::Pointer, const BoundarySide&, const int&>())
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "BendingStripNURBSPatch" << TDim << "DPointer";
    class_<typename BendingStripNURBSPatch<TDim>::Pointer>
    (ss.str().c_str(), init<typename BendingStripNURBSPatch<TDim>::Pointer>())
    .def("GetReference", GetReference<BendingStripNURBSPatch<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
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
    .def("MakeNeighbor", &MultiPatch_MakeNeighbor<TDim>)
    .def("EquationSystemSize", &MultiPatch<TDim>::EquationSystemSize)
    .def("Enumerate", &MultiPatch_Enumerate1<TDim>)
    .def("Enumerate", &MultiPatch_Enumerate2<TDim>)
    .def("PrintAddress", &MultiPatch<TDim>::PrintAddress)
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

void IsogeometricApplication_AddNURBSToPython()
{
    /////////////////////////////////////////////////////////////////
    ///////////////////////SUPPORT DOMAIN////////////////////////////
    /////////////////////////////////////////////////////////////////

    class_<DomainManager, DomainManager::Pointer, boost::noncopyable>
    ("DomainManager", init<std::size_t>())
    ;

    class_<DomainManager2D, DomainManager2D::Pointer, boost::noncopyable>
    ("DomainManager2D", init<std::size_t>())
    .def("AddXcoord", &DomainManager2D::AddXcoord)
    .def("AddYcoord", &DomainManager2D::AddYcoord)
    .def("AddCell", &DomainManager2D_AddCell)
    .def("IsInside", &DomainManager2D_IsInside)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////TRANSFORMATION////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddTransformation();

    /////////////////////////////////////////////////////////////////
    ///////////////////CONTROL POINT/////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlPoint();

    /////////////////////////////////////////////////////////////////
    ///////////////////////CONTROL GRIDS/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddControlGrids();

    class_<ControlGridLibrary, ControlGridLibrary::Pointer, boost::noncopyable>
    ("ControlGridLibrary", init<>())
    .def("CreateLinearControlPointGrid", &ControlGridLibrary_CreateLinearControlPointGrid)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid1)
    .def("CreateRectangularControlPointGrid", &ControlGridLibrary_CreateRectangularControlPointGrid2)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid1)
    .def("CreateCubicControlPointGrid", &ControlGridLibrary_CreateCubicControlPointGrid2)
    .def("CreateLinearZeroDoubleControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<double> >)
    .def("CreateRectangularZeroDoubleControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<double> >)
    .def("CreateCubicZeroDoubleControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<double> >)
    .def("CreateLinearZeroArray1DControlGrid", &ControlGridLibrary_CreateLinearZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateRectangularZeroArray1DControlGrid", &ControlGridLibrary_CreateRectangularZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    .def("CreateCubicZeroArray1DControlGrid", &ControlGridLibrary_CreateCubicZeroControlGridWithVariable<Variable<array_1d<double, 3> > >)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////FESpace///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddFESpacesToPython<1>();
    IsogeometricApplication_AddFESpacesToPython<2>();
    IsogeometricApplication_AddFESpacesToPython<3>();

    class_<BSplinesFESpaceLibrary, BSplinesFESpaceLibrary::Pointer, boost::noncopyable>
    ("BSplinesFESpaceLibrary", init<>())
    .def("CreateLinearFESpace", &BSplinesFESpaceLibrary_CreateLinearFESpace)
    .def("CreateRectangularFESpace", &BSplinesFESpaceLibrary_CreateRectangularFESpace)
    .def("CreateCubicFESpace", &BSplinesFESpaceLibrary_CreateCubicFESpace)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////Grid Functions////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddGridFunctionsToPython<1>();
    IsogeometricApplication_AddGridFunctionsToPython<2>();
    IsogeometricApplication_AddGridFunctionsToPython<3>();

    /////////////////////////////////////////////////////////////////
    ///////////////////////PATCHES///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    enum_<BoundarySide>("BoundarySide")
    .value("Left", _LEFT_)
    .value("Right", _RIGHT_)
    .value("Top", _TOP_)
    .value("Bottom", _BOTTOM_)
    .value("Front", _FRONT_)
    .value("Back", _BACK_)
    ;

    enum_<IsogeometricEchoFlags>("IsogeometricEchoFlags")
    .value("ECHO_REFIMENT", ECHO_REFIMENT)
    ;

    IsogeometricApplication_AddPatchesToPython<1>();
    IsogeometricApplication_AddPatchesToPython<2>();
    IsogeometricApplication_AddPatchesToPython<3>();

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

