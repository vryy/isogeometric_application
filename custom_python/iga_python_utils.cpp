//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 17 Sep 2023 $
//   Revision:            $Revision: 1.0 $
//
//

#include "custom_python/iga_python_utils.h"

namespace Kratos
{

namespace Python
{

template<typename TDataType>
void IsogeometricPythonUtils::Unpack(const boost::python::dict& rPatchNodalValues,
        std::map<std::size_t, std::map<std::size_t, TDataType> >& patch_nodal_values)
{
    boost::python::list keys = rPatchNodalValues.keys();

    typedef boost::python::stl_input_iterator<int> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& id,
                  std::make_pair(iterator_type(keys), // begin
                    iterator_type() ) ) // end
    {
        boost::python::object o = rPatchNodalValues.get(id);

        // here assumed that the patch_data given as a dict
        boost::python::dict values = boost::python::extract<boost::python::dict>(o);

        boost::python::list keys2 = values.keys();

        std::map<std::size_t, TDataType> nodal_values;
        BOOST_FOREACH(const typename iterator_type::value_type& id2,
                      std::make_pair(iterator_type(keys2), // begin
                        iterator_type() ) ) // end
        {
            boost::python::object o2 = values.get(id2);
            TDataType v = boost::python::extract<TDataType>(o2);
            nodal_values[static_cast<std::size_t>(id2)] = v;
        }

        patch_nodal_values[static_cast<std::size_t>(id)] = nodal_values;
    }
}

//// template function instantiation

template
void IsogeometricPythonUtils::Unpack<double>(const boost::python::dict& rPatchNodalValues,
        std::map<std::size_t, std::map<std::size_t, double> >& patch_nodal_values);
template
void IsogeometricPythonUtils::Unpack<array_1d<double, 3> >(const boost::python::dict& rPatchNodalValues,
        std::map<std::size_t, std::map<std::size_t, array_1d<double, 3> > >& patch_nodal_values);
template
void IsogeometricPythonUtils::Unpack<Vector>(const boost::python::dict& rPatchNodalValues,
        std::map<std::size_t, std::map<std::size_t, Vector> >& patch_nodal_values);

} // namespace Python

} // namespace Kratos.
