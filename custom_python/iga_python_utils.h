//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 17 Sep 2023 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_IGA_PYTHON_UTILS_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_IGA_PYTHON_UTILS_H_INCLUDED

// System includes
#include <vector>
#include <array>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/multipatch_model_part.h"

namespace Kratos
{

namespace Python
{

struct IsogeometricPythonUtils
{

    template<typename TInputValueType, typename TOutputValueType>
    static void Unpack(const boost::python::list& rValues, std::vector<TOutputValueType>& values)
    {
        typedef boost::python::stl_input_iterator<TInputValueType> iterator_value_type;
        BOOST_FOREACH(const typename iterator_value_type::value_type & v,
                      std::make_pair(iterator_value_type(rValues), iterator_value_type() ) )
        {
            values.push_back(static_cast<TOutputValueType>(v));
        }
    }

    template<typename TInputValueType, typename TOutputValueType>
    static std::size_t Unpack(const boost::python::list& rValues, std::vector<TOutputValueType>& values, std::size_t Dim)
    {
        if (values.size() != Dim)
        {
            values.resize(Dim);
        }

        std::size_t dim = 0;
        typedef boost::python::stl_input_iterator<TInputValueType> iterator_value_type;
        BOOST_FOREACH(const typename iterator_value_type::value_type & v,
                      std::make_pair(iterator_value_type(rValues), iterator_value_type() ) )
        {
            values[dim++] = static_cast<TOutputValueType>(v);
            if (dim == Dim)
            {
                break;
            }
        }

        return dim;
    }

    template<typename TInputValueType, typename TOutputValueType>
    static void Unpack(const boost::python::list& rValues, std::vector<std::vector<TOutputValueType> >& values)
    {
        typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type & rSubValues,
                      std::make_pair(iterator_value_type(rValues), // begin
                                     iterator_value_type() ) ) // end
        {
            std::vector<TOutputValueType> sub_values;

            typedef boost::python::stl_input_iterator<TInputValueType> iterator_value_type2;
            BOOST_FOREACH(const typename iterator_value_type2::value_type & x,
                          std::make_pair(iterator_value_type2(rSubValues), // begin
                                         iterator_value_type2() ) ) // end
            {
                sub_values.push_back(static_cast<TOutputValueType>(x));
            }

            values.push_back(sub_values);
        }
    }

    template<typename TInputValueType, typename TOutputValueType>
    static std::size_t Unpack(const boost::python::list& rValues, std::vector<std::vector<TOutputValueType> >& values, std::size_t Dim)
    {
        if (values.size() != Dim)
        {
            values.resize(Dim);
        }

        std::size_t dim = 0;
        typedef boost::python::stl_input_iterator<boost::python::list> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type & rSubValues,
                      std::make_pair(iterator_value_type(rValues), // begin
                                     iterator_value_type() ) ) // end
        {
            std::vector<TOutputValueType> sub_values;

            typedef boost::python::stl_input_iterator<TInputValueType> iterator_value_type2;
            BOOST_FOREACH(const typename iterator_value_type2::value_type & x,
                          std::make_pair(iterator_value_type2(rSubValues), // begin
                                         iterator_value_type2() ) ) // end
            {
                sub_values.push_back(static_cast<TOutputValueType>(x));
            }

            values[dim++] = sub_values;
            if (dim == Dim)
            {
                break;
            }
        }

        return dim;
    }

    template<class T>
    static void Unpack(const boost::python::list& rPatchList, std::vector<typename T::PatchType::Pointer>& pPatches)
    {
        typedef boost::python::stl_input_iterator<typename T::PatchType::Pointer> iterator_value_type;
        BOOST_FOREACH(const typename iterator_value_type::value_type & v,
                      std::make_pair(iterator_value_type(rPatchList), iterator_value_type() ) )
        {
            pPatches.push_back(v);
        }
    }

    template<typename TDataType>
    static void Unpack(const boost::python::dict& rPatchNodalValues,
                       std::map<std::size_t, std::map<std::size_t, TDataType> >& patch_nodal_values);

    template<typename TInputValueType, typename TOutputValueType, std::size_t TSize>
    static void Unpack(const boost::python::list& rValues, std::array<TOutputValueType, TSize>& values)
    {
        if (boost::python::len(rValues) < TSize)
            KRATOS_ERROR << "Input list does not have sufficient dimension";
        std::size_t cnt = 0;
        typedef boost::python::stl_input_iterator<TInputValueType> iterator_value_type;
        BOOST_FOREACH(const typename iterator_value_type::value_type & v,
                      std::make_pair(iterator_value_type(rValues), iterator_value_type() ) )
        {
            values[cnt] = (static_cast<TOutputValueType>(v));
            ++cnt;
        }
    }

};

} // namespace Python

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_IGA_PYTHON_UTILS_H_INCLUDED defined
