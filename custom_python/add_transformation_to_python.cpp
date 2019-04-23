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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/trans/translation.h"
#include "custom_utilities/trans/rotation.h"
#include "custom_utilities/trans/mirror.h"
#include "custom_utilities/trans/transformation_utility.h"
#include "custom_python/add_transformation_to_python.h"


namespace Kratos
{

namespace Python
{

//////////////////////////////////////////////////

template<typename TDataType>
Transformation<TDataType> TransformationUtility_CreateAlignTransformation(TransformationUtility<TDataType>& rDummy,
    const typename Transformation<TDataType>::VectorType& a, const typename Transformation<TDataType>::VectorType& b)
{
    Transformation<TDataType> T;
    T = rDummy.CreateAlignTransformation(a, b);
    return T;
}

template<typename TDataType>
void Transformation_SetValue(Transformation<TDataType>& rDummy, const int& i, const int& j, const TDataType& v)
{
    rDummy(i, j) = v;
}

template<typename TDataType>
TDataType Transformation_GetValue(Transformation<TDataType>& rDummy, const int& i, const int& j)
{
    return rDummy(i, j);
}

template<typename TDataType, typename TVectorType>
TVectorType Transformation_Apply(Transformation<TDataType>& rDummy, const TVectorType& v)
{
    TVectorType newv = v;
    rDummy.template ApplyTransformation<TVectorType>(newv);
    return newv;
}

template<typename TDataType>
pybind11::list Transformation_Apply2(Transformation<TDataType>& rDummy, pybind11::list py_trans_list)
{
    std::vector<TDataType> newv;
    for (auto v : py_trans_list)
        newv.push_back(v.cast<TDataType>());

    rDummy.template ApplyTransformation<std::vector<TDataType> >(newv);

    pybind11::list res;
    for (std::size_t i = 0; i < newv.size(); ++i)
        res.append(newv[i]);

    return res;
}

template<typename TDataType>
array_1d<TDataType, 3> Transformation_P(Transformation<TDataType>& rDummy)
{
    return rDummy.P();
}

template<typename TDataType>
array_1d<TDataType, 3> Transformation_V1(Transformation<TDataType>& rDummy)
{
    return rDummy.V1();
}

template<typename TDataType>
array_1d<TDataType, 3> Transformation_V2(Transformation<TDataType>& rDummy)
{
    return rDummy.V2();
}

template<typename TDataType>
array_1d<TDataType, 3> Transformation_V3(Transformation<TDataType>& rDummy)
{
    return rDummy.V3();
}

//////////////////////////////////////////////////

void IsogeometricApplication_AddTransformationToPython(pybind11::module& m)
{

    typedef Transformation<double>::VectorType VectorType;

    pybind11::class_<Transformation<double>, Transformation<double>::Pointer>(m, "Transformation")
    .def(pybind11::init<>())
    .def(pybind11::init<const VectorType&, const VectorType&, const VectorType&>())
    .def(pybind11::init<const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&>())
    .def(pybind11::init<const VectorType&, const VectorType&, const VectorType&, const VectorType&>())
    .def(pybind11::init<const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&>())
    .def("AppendTransformation", &Transformation<double>::AppendTransformation)
    .def("PrependTransformation", &Transformation<double>::PrependTransformation)
    .def("Inverse", &Transformation<double>::Inverse)
    // .def(pybind11::operators<pybind11::op_mul>());
    .def("P", &Transformation_P<double>)
    .def("V1", &Transformation_V1<double>)
    .def("V2", &Transformation_V2<double>)
    .def("V3", &Transformation_V3<double>)
    .def("SetValue", &Transformation_SetValue<double>)
    .def("GetValue", &Transformation_GetValue<double>)
    .def("Apply", &Transformation_Apply<double, Vector>)
    .def("Apply", &Transformation_Apply<double, array_1d<double, 3> >)
    .def("Apply", &Transformation_Apply2<double>)
    .def("__str__", &PrintObject<Transformation<double> >)
    ;

    pybind11::class_<Translation<double>, Translation<double>::Pointer, Transformation<double> >
    (m, "Translation")
    .def(pybind11::init<const double&, const double&, const double&>())
    .def("__str__", &PrintObject<Translation<double> >)
    ;

    pybind11::class_<Rotation<0, double>, Rotation<0, double>::Pointer, Transformation<double> >
    (m, "RotationX")
    .def(pybind11::init<const double&>())
    .def("__str__", &PrintObject<Rotation<0, double> >)
    ;

    pybind11::class_<Rotation<1, double>, Rotation<1, double>::Pointer, Transformation<double> >
    (m, "RotationY")
    .def(pybind11::init<const double&>())
    .def("__str__", &PrintObject<Rotation<1, double> >)
    ;

    pybind11::class_<Rotation<2, double>, Rotation<2, double>::Pointer, Transformation<double> >
    (m, "RotationZ")
    .def(pybind11::init<const double&>())
    .def("__str__", &PrintObject<Rotation<2, double> >)
    ;

    pybind11::class_<Mirror<0, double>, Mirror<0, double>::Pointer, Transformation<double> >
    (m, "MirrorX")
    .def(pybind11::init<>())
    .def("__str__", &PrintObject<Mirror<0, double> >)
    ;

    pybind11::class_<Mirror<1, double>, Mirror<1, double>::Pointer, Transformation<double> >
    (m, "MirrorY")
    .def(pybind11::init<>())
    .def("__str__", &PrintObject<Mirror<1, double> >)
    ;

    pybind11::class_<Mirror<2, double>, Mirror<2, double>::Pointer, Transformation<double> >
    (m, "MirrorZ")
    .def(pybind11::init<>())
    .def("__str__", &PrintObject<Mirror<2, double> >)
    ;

    pybind11::class_<TransformationUtility<double>, TransformationUtility<double>::Pointer>
    (m, "TransformationUtility")
    .def(pybind11::init<>())
    .def("CreateAlignTransformation", &TransformationUtility_CreateAlignTransformation<double>)
    ;

}

}  // namespace Python.

} // Namespace Kratos

