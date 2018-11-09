/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 18, 2013 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"


namespace Kratos
{

namespace Python
{

template<class TPatchType>
TPatchType& GetReference(typename TPatchType::Pointer& dummy)
{
    return *dummy;
}

template<typename TDataType>
struct ControlValue_Helper
{
    static boost::python::list GetValue(const TDataType& rValue)
    {
        boost::python::list dummy;
        return dummy;
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

template<typename TType>
std::size_t Isogeometric_GetId(TType& rDummy)
{
    return rDummy.Id();
}

template<typename TType>
void Isogeometric_SetId(TType& rDummy, std::size_t Id)
{
    rDummy.SetId(Id);
}

template<typename TType>
void Isogeometric_DoNotSetId(TType& rDummy, std::size_t Id)
{
    // DO NOTHING
}

template<typename TType>
std::size_t Isogeometric_GetEquationId(TType& rDummy)
{
    return rDummy.EquationId();
}

template<typename TType>
void Isogeometric_SetEquationId(TType& rDummy, std::size_t EquationId)
{
    rDummy.SetEquationId(EquationId);
}


template<typename TFESpaceType>
boost::python::list FESpace_ExtractBoundaryBfsByFlag(TFESpaceType& rDummy, std::size_t boundary_id)
{
    typedef typename TFESpaceType::bf_t bf_t;

    std::vector<bf_t> bf_list = rDummy.ExtractBoundaryBfsByFlag(boundary_id);

    boost::python::list Output;
    for (std::size_t i = 0; i < bf_list.size(); ++i)
        Output.append(bf_list[i]);

    return Output;
}

template<typename TFESpaceType>
typename TFESpaceType::bf_t FESpace_GetItem(TFESpaceType& rDummy, std::size_t i)
{
    return rDummy[i];
}

void  IsogeometricApplication_AddBackendUtilitiesToPython();
void  IsogeometricApplication_AddFrontendUtilitiesToPython();
void  IsogeometricApplication_AddTransformationToPython();
void  IsogeometricApplication_AddControlGridsToPython();
void  IsogeometricApplication_AddPatchesToPython();
void  IsogeometricApplication_AddNURBSToPython();
void  IsogeometricApplication_AddPBSplinesToPython();
void  IsogeometricApplication_AddHBSplinesToPython();
void  IsogeometricApplication_AddTSplinesToPython();
void  IsogeometricApplication_AddMeshAndModelPartToPython();

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined

