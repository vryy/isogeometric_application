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

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_NURBS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_NURBS_TO_PYTHON_H_INCLUDED

// System includes

// External includes
#include<pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/control_point.h"


namespace Kratos
{

namespace Python
{

/////////////////////////////////////////////////

template<typename TDataType>
struct ControlValue_Helper
{
    static pybind11::list GetValue(const TDataType& rValue)
    {
        pybind11::list dummy;
        return dummy;
    }

    static TDataType GetValue(pybind11::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<double>
{
    static pybind11::list GetValue(const double& rValue)
    {
        pybind11::list v;
        v.append(rValue);
        return v;
    }

    static double GetValue(pybind11::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<ControlPoint<double> >
{
    static pybind11::list GetValue(const ControlPoint<double>& rValue)
    {
        pybind11::list v;
        v.append(rValue.WX());
        v.append(rValue.WY());
        v.append(rValue.WZ());
        v.append(rValue.W());
        return v;
    }

    static double GetValue(pybind11::list rValue)
    {
    }
};

template<>
struct ControlValue_Helper<array_1d<double, 3> >
{
    static pybind11::list GetValue(const array_1d<double, 3>& rValue)
    {
        pybind11::list v;
        v.append(rValue[0]);
        v.append(rValue[1]);
        v.append(rValue[2]);
        return v;
    }

    static double GetValue(pybind11::list rValue)
    {
    }
};

/////////////////////////////////////////////////

void  IsogeometricApplication_AddNURBSToPython(pybind11::module& m);

}  // namespace Python.

} // Namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_NURBS_TO_PYTHON_H_INCLUDED
