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
#include "python/python_utils.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/isogeometric_utility.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/nurbs/structured_control_grid.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/nurbs/bending_strip_nurbs_patch.h"
#include "custom_utilities/nurbs/nurbs_test_utils.h"
#include "custom_python/iga_define_python.h"
#include "custom_python/add_nurbs_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

template<int TDim, int TWhichDim>
boost::python::list BSplinesFESpace_GetKnotVector(BSplinesFESpace<TDim>& rDummy)
{
    boost::python::list knot_list;

    if (TWhichDim < TDim)
    {
        const typename BSplinesFESpace<TDim>::knot_container_t& knot_vector = rDummy.KnotVector(TWhichDim);

        for (std::size_t i = 0; i < knot_vector.size(); ++i)
        {
            knot_list.append(knot_vector[i]);
        }
    }

    return knot_list;
}

template<int TDim, int TWhichDim>
void BSplinesFESpace_SetKnotVector(BSplinesFESpace<TDim>& rDummy, const boost::python::list& knot_list)
{
    if (TWhichDim < TDim)
    {
        std::vector<double> knot_vec;
        PythonUtils::Unpack<double, double>(knot_list, knot_vec);

        rDummy.SetKnotVector(TWhichDim, knot_vec);
    }
}

////////////////////////////////////////

BSplinesFESpace<1>::Pointer BSplinesFESpaceLibrary_CreatePrimitiveFESpace1(BSplinesFESpaceLibrary& rDummy, std::size_t order_u)
{
    std::vector<std::size_t> orders(1);
    orders[0] = order_u;
    return rDummy.CreatePrimitiveFESpace<1>(orders);
}

BSplinesFESpace<2>::Pointer BSplinesFESpaceLibrary_CreatePrimitiveFESpace2(BSplinesFESpaceLibrary& rDummy, std::size_t order_u, std::size_t order_v)
{
    std::vector<std::size_t> orders(2);
    orders[0] = order_u;
    orders[1] = order_v;
    return rDummy.CreatePrimitiveFESpace<2>(orders);
}

BSplinesFESpace<3>::Pointer BSplinesFESpaceLibrary_CreatePrimitiveFESpace3(BSplinesFESpaceLibrary& rDummy, std::size_t order_u, std::size_t order_v, std::size_t order_w)
{
    std::vector<std::size_t> orders(3);
    orders[0] = order_u;
    orders[1] = order_v;
    orders[2] = order_w;
    return rDummy.CreatePrimitiveFESpace<3>(orders);
}

BSplinesFESpace<1>::Pointer BSplinesFESpaceLibrary_CreateUniformFESpace1(BSplinesFESpaceLibrary& rDummy,
        std::size_t number_u, std::size_t order_u)
{
    std::vector<std::size_t> numbers(1);
    numbers[0] = number_u;
    std::vector<std::size_t> orders(1);
    orders[0] = order_u;
    return rDummy.CreateUniformFESpace<1>(numbers, orders);
}

BSplinesFESpace<2>::Pointer BSplinesFESpaceLibrary_CreateUniformFESpace2(BSplinesFESpaceLibrary& rDummy,
        std::size_t number_u, std::size_t order_u,
        std::size_t number_v, std::size_t order_v)
{
    std::vector<std::size_t> numbers(2);
    numbers[0] = number_u;
    numbers[1] = number_v;
    std::vector<std::size_t> orders(2);
    orders[0] = order_u;
    orders[1] = order_v;
    return rDummy.CreateUniformFESpace<2>(numbers, orders);
}

BSplinesFESpace<3>::Pointer BSplinesFESpaceLibrary_CreateUniformFESpace3(BSplinesFESpaceLibrary& rDummy,
        std::size_t number_u, std::size_t order_u,
        std::size_t number_v, std::size_t order_v,
        std::size_t number_w, std::size_t order_w)
{
    std::vector<std::size_t> numbers(3);
    numbers[0] = number_u;
    numbers[1] = number_v;
    numbers[2] = number_w;
    std::vector<std::size_t> orders(3);
    orders[0] = order_u;
    orders[1] = order_v;
    orders[2] = order_w;
    return rDummy.CreateUniformFESpace<3>(numbers, orders);
}

//////////////////////////////////////////////////

void DomainManager2D_AddCell(DomainManager2D& rDummy, double x1, double x2, double y1, double y2)
{
    std::vector<double> box(4);
    box[0] = x1;
    box[1] = x2;
    box[2] = y1;
    box[3] = y2;

    rDummy.AddCell(box);
}

bool DomainManager2D_IsInside(DomainManager2D& rDummy, double x1, double x2, double y1, double y2)
{
    std::vector<double> box(4);
    box[0] = x1;
    box[1] = x2;
    box[2] = y1;
    box[3] = y2;

    return rDummy.IsInside(box);
}

//////////////////////////////////////////////////

template<int TDim>
void NURBSTestUtils_ProbeAndTestValuesOnPatch(NURBSTestUtils<TDim>& dummy, typename Patch<TDim>::Pointer pPatch,
        ModelPart::ConditionsContainerType& rConditions, int integration_order, double tol)
{
    GeometryData::IntegrationMethod integration_method = IsogeometricUtility::GetIntegrationMethod(integration_order);
    dummy.ProbeAndTestValuesOnPatch(pPatch, rConditions, integration_method, tol);
}

//////////////////////////////////////////////////

template<int TDim, typename TDataType>
struct StructuredControlGrid_Helper
{
    static boost::python::list GetValue(StructuredControlGrid<TDim, TDataType>& rDummy)
    {
        KRATOS_ERROR << "Not yet implemented";
    }

    static void SetValue(StructuredControlGrid<TDim, TDataType>& rDummy, const boost::python::list& values)
    {
        KRATOS_ERROR << "Not yet implemented";
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

    static void SetValue(StructuredControlGrid<1, TDataType>& rDummy, const boost::python::list& values)
    {
        KRATOS_ERROR << "Not yet implemented";
    }

    static void SetValue1D(StructuredControlGrid<1, TDataType>& rDummy, std::size_t i, const TDataType& value)
    {
        rDummy.SetValue(i, value);
    }

    static TDataType GetValue1D(StructuredControlGrid<1, TDataType>& rDummy, std::size_t i)
    {
        return rDummy.GetValue(i);
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

    static void SetValue(StructuredControlGrid<2, TDataType>& rDummy, const boost::python::list& values)
    {
        KRATOS_ERROR << "Not yet implemented";
    }

    static void SetValue2D(StructuredControlGrid<2, TDataType>& rDummy, std::size_t i, std::size_t j, const TDataType& value)
    {
        rDummy.SetValue(i, j, value);
    }

    static TDataType GetValue2D(StructuredControlGrid<2, TDataType>& rDummy, std::size_t i, std::size_t j)
    {
        return rDummy.GetValue(i, j);
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

    static void SetValue(StructuredControlGrid<3, TDataType>& rDummy, const boost::python::list& values)
    {
        KRATOS_ERROR << "Not yet implemented";
    }

    static void SetValue3D(StructuredControlGrid<3, TDataType>& rDummy, std::size_t i, std::size_t j, std::size_t k, const TDataType& value)
    {
        rDummy.SetValue(i, j, k, value);
    }

    static TDataType GetValue3D(StructuredControlGrid<3, TDataType>& rDummy, std::size_t i, std::size_t j, std::size_t k)
    {
        return rDummy.GetValue(i, j, k);
    }
};

void IsogeometricApplication_AddStructuredControlGridsToPython()
{
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
    ("StructuredControlPointGrid1D", init<std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<1, ControlPoint<double> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<1, ControlPoint<double> >::SetValue1D)
    .def("GetValue", &StructuredControlGrid_Helper<1, ControlPoint<double> >::GetValue1D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, double>, StructuredControlGrid<1, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid1D", init<std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, double>::GetValue, &StructuredControlGrid_Helper<1, double>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<1, double>::SetValue1D)
    .def("GetValue", &StructuredControlGrid_Helper<1, double>::GetValue1D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, array_1d<double, 3> >, StructuredControlGrid<1, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid1D", init<std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<1, array_1d<double, 3> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::SetValue1D)
    .def("GetValue", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::GetValue1D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<1, Vector>, StructuredControlGrid<1, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid1D", init<std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, Vector>::GetValue, &StructuredControlGrid_Helper<1, Vector>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<1, Vector>::SetValue1D)
    .def("GetValue", &StructuredControlGrid_Helper<1, Vector>::GetValue1D)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<StructuredControlGrid<2, ControlPoint<double> >, StructuredControlGrid<2, ControlPoint<double> >::Pointer, bases<BaseStructuredControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("StructuredControlPointGrid2D", init<std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<2, ControlPoint<double> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<2, ControlPoint<double> >::SetValue2D)
    .def("GetValue", &StructuredControlGrid_Helper<2, ControlPoint<double> >::GetValue2D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, double>, StructuredControlGrid<2, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid2D", init<std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, double>::GetValue, &StructuredControlGrid_Helper<2, double>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<2, double>::SetValue2D)
    .def("GetValue", &StructuredControlGrid_Helper<2, double>::GetValue2D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, array_1d<double, 3> >, StructuredControlGrid<2, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid2D", init<std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<2, array_1d<double, 3> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<2, array_1d<double, 3> >::SetValue2D)
    .def("GetValue", &StructuredControlGrid_Helper<2, array_1d<double, 3> >::GetValue2D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<2, Vector>, StructuredControlGrid<2, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid2D", init<std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<2, Vector>::GetValue, &StructuredControlGrid_Helper<2, Vector>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<2, Vector>::SetValue2D)
    .def("GetValue", &StructuredControlGrid_Helper<2, Vector>::GetValue2D)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////

    class_<StructuredControlGrid<3, ControlPoint<double> >, StructuredControlGrid<3, ControlPoint<double> >::Pointer, bases<BaseStructuredControlGrid<ControlPoint<double> > >, boost::noncopyable>
    ("StructuredControlPointGrid3D", init<std::size_t, std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, ControlPoint<double> >::GetValue, &StructuredControlGrid_Helper<3, ControlPoint<double> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<3, ControlPoint<double> >::SetValue3D)
    .def("GetValue", &StructuredControlGrid_Helper<3, ControlPoint<double> >::GetValue3D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, double>, StructuredControlGrid<3, double>::Pointer, bases<BaseStructuredControlGrid<double> >, boost::noncopyable>
    ("StructuredDoubleControlGrid3D", init<std::size_t, std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, double>::GetValue, &StructuredControlGrid_Helper<3, double>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<3, double>::SetValue3D)
    .def("GetValue", &StructuredControlGrid_Helper<3, double>::GetValue3D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, array_1d<double, 3> >, StructuredControlGrid<3, array_1d<double, 3> >::Pointer, bases<BaseStructuredControlGrid<array_1d<double, 3> > >, boost::noncopyable>
    ("StructuredArray1DControlGrid3D", init<std::size_t, std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<1, array_1d<double, 3> >::GetValue, &StructuredControlGrid_Helper<1, array_1d<double, 3> >::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<3, array_1d<double, 3> >::SetValue3D)
    .def("GetValue", &StructuredControlGrid_Helper<3, array_1d<double, 3> >::GetValue3D)
    .def(self_ns::str(self))
    ;

    class_<StructuredControlGrid<3, Vector>, StructuredControlGrid<3, Vector>::Pointer, bases<BaseStructuredControlGrid<Vector> >, boost::noncopyable>
    ("StructuredVectorControlGrid3D", init<std::size_t, std::size_t, std::size_t>())
    .add_property("ControlValues", &StructuredControlGrid_Helper<3, Vector>::GetValue, &StructuredControlGrid_Helper<3, Vector>::SetValue)
    .def("SetValue", &StructuredControlGrid_Helper<3, Vector>::SetValue3D)
    .def("GetValue", &StructuredControlGrid_Helper<3, Vector>::GetValue3D)
    .def(self_ns::str(self))
    ;

    /////////////////////////////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddBSplinesFESpaceToPython()
{
    std::stringstream ss;
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

//////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddBendingStripNURBSToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "BendingStripNURBSPatch" << TDim << "D";
    class_<BendingStripNURBSPatch<TDim>, bases<PatchInterface<TDim>, Patch<TDim> > >
    // class_<BendingStripNURBSPatch<TDim>, typename BendingStripNURBSPatch<TDim>::Pointer >
    (ss.str().c_str(), init<std::size_t, int>())
    .def(init<std::size_t, typename Patch<TDim>::Pointer, const BoundarySide&, typename Patch<TDim>::Pointer, const BoundarySide&, int>())
    // .def(self_ns::str(self_ns::self))
    .def(self_ns::str(self))
    ;

    ss.str(std::string());
    ss << "BendingStripNURBSPatch" << TDim << "DPointer";
    class_<typename BendingStripNURBSPatch<TDim>::Pointer>
    (ss.str().c_str(), init<typename BendingStripNURBSPatch<TDim>::Pointer>())
    .def("GetReference", GetReference<BendingStripNURBSPatch<TDim> >, return_value_policy<reference_existing_object>())
    .def(self_ns::str(self))
    ;
}

//////////////////////////////////////////////////

template<int TDim>
void IsogeometricApplication_AddNURBSTestUtilsToPython()
{
    std::stringstream ss;

    ss.str(std::string());
    ss << "NURBSTestUtils" << TDim << "D";
    class_<NURBSTestUtils<TDim>, typename NURBSTestUtils<TDim>::Pointer, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("ProbeAndTestValuesOnPatch", &NURBSTestUtils_ProbeAndTestValuesOnPatch<TDim>)
    .def(self_ns::str(self))
    ;
}

//////////////////////////////////////////////////

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
    ///////////////////////CONTROL GRIDS/////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddStructuredControlGridsToPython();

    /////////////////////////////////////////////////////////////////
    ///////////////////////FESpace///////////////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddBSplinesFESpaceToPython<1>();
    IsogeometricApplication_AddBSplinesFESpaceToPython<2>();
    IsogeometricApplication_AddBSplinesFESpaceToPython<3>();

    class_<BSplinesFESpaceLibrary, BSplinesFESpaceLibrary::Pointer, boost::noncopyable>
    ("BSplinesFESpaceLibrary", init<>())
    .def("CreateLinearFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace1) // backward compatibility
    .def("CreateRectangularFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace2) // backward compatibility
    .def("CreateCubicFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace3) // backward compatibility
    .def("CreatePrimitiveFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace1)
    .def("CreatePrimitiveFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace2)
    .def("CreatePrimitiveFESpace", &BSplinesFESpaceLibrary_CreatePrimitiveFESpace3)
    .def("CreateUniformFESpace", &BSplinesFESpaceLibrary_CreateUniformFESpace1)
    .def("CreateUniformFESpace", &BSplinesFESpaceLibrary_CreateUniformFESpace2)
    .def("CreateUniformFESpace", &BSplinesFESpaceLibrary_CreateUniformFESpace3)
    ;

    /////////////////////////////////////////////////////////////////
    ///////////////////////Bending Strip NURBS Patch/////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddBendingStripNURBSToPython<2>();
    IsogeometricApplication_AddBendingStripNURBSToPython<3>();

    /////////////////////////////////////////////////////////////////
    ///////////////////////NURBS Tets Utils//////////////////////////
    /////////////////////////////////////////////////////////////////

    IsogeometricApplication_AddNURBSTestUtilsToPython<1>();
    // IsogeometricApplication_AddNURBSTestUtilsToPython<2>();
    // IsogeometricApplication_AddNURBSTestUtilsToPython<3>();

}

}  // namespace Python.

} // Namespace Kratos
