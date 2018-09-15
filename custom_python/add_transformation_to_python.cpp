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
#include "includes/variables.h"
#include "custom_python/add_utilities_to_python.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/trans/translation.h"
#include "custom_utilities/trans/rotation.h"
#include "custom_utilities/trans/transformation_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

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
boost::python::list Transformation_Apply2(Transformation<TDataType>& rDummy, boost::python::list v)
{
    std::vector<TDataType> newv;
    typedef boost::python::stl_input_iterator<TDataType> iterator_value_type;
    BOOST_FOREACH(const typename iterator_value_type::value_type& d, std::make_pair(iterator_value_type(v), iterator_value_type() ) )
        newv.push_back(d);

    rDummy.template ApplyTransformation<std::vector<TDataType> >(newv);

    boost::python::list res;
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

void IsogeometricApplication_AddTransformationToPython()
{
    typedef Transformation<double>::VectorType VectorType;

    class_<Transformation<double>, Transformation<double>::Pointer>
    ("Transformation", init<>())
    .def(init<const VectorType&, const VectorType&, const VectorType&>())
    .def(init<const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&>())
    .def(init<const VectorType&, const VectorType&, const VectorType&, const VectorType&>())
    .def(init<const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&, const array_1d<double, 3>&>())
    .def("AppendTransformation", &Transformation<double>::AppendTransformation)
    .def("PrependTransformation", &Transformation<double>::PrependTransformation)
    .def("Inverse", &Transformation<double>::Inverse)
    // .def(boost::python::operators<boost::python::op_mul>());
    .def("P", &Transformation_P<double>)
    .def("V1", &Transformation_V1<double>)
    .def("V2", &Transformation_V2<double>)
    .def("V3", &Transformation_V3<double>)
    .def("SetValue", &Transformation_SetValue<double>)
    .def("GetValue", &Transformation_GetValue<double>)
    .def("Apply", &Transformation_Apply<double, Vector>)
    .def("Apply", &Transformation_Apply<double, array_1d<double, 3> >)
    .def("Apply", &Transformation_Apply2<double>)
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

    class_<TransformationUtility<double>, TransformationUtility<double>::Pointer, boost::noncopyable>
    ("TransformationUtility", init<>())
    .def("CreateAlignTransformation", &TransformationUtility_CreateAlignTransformation<double>)
    ;

}

}  // namespace Python.

} // Namespace Kratos

