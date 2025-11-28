//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2025 $
//
//

#include "custom_utilities/fespace_utility.h"
#include "custom_utilities/weighted_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"

namespace Kratos
{

template<int TDim, typename TLocalCoordinateType>
typename FESpace<TDim, TLocalCoordinateType>::Pointer FESpaceUtility::CreateEmptyFESpace(const std::string& type)
{
    const std::string wof_name = ("WeightedFESpace" + std::to_string(TDim) + "D_over_BSplinesFESpace" + std::to_string(TDim) + "D");

    if (type == BSplinesFESpace<TDim, TLocalCoordinateType>::StaticType())
    {
        return BSplinesFESpace<TDim, TLocalCoordinateType>::Create();
    }
    else if (type == WeightedFESpace<TDim, TLocalCoordinateType, double>::StaticType())
    {
        return WeightedFESpace<TDim, TLocalCoordinateType, double>::Create(FESpace<TDim, TLocalCoordinateType>::Create(), {});
    }
    else if (type == wof_name) // hack
    {
        return WeightedFESpace<TDim, TLocalCoordinateType, double>::Create(BSplinesFESpace<TDim, TLocalCoordinateType>::Create(), {});
    }
    else
        return FESpace<TDim, TLocalCoordinateType>::Create();
}

/* template class instantiation */
template KRATOS_API(ISOGEOMETRIC_APPLICATION) FESpace<0, double>::Pointer FESpaceUtility::CreateEmptyFESpace<0, double>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) FESpace<1, double>::Pointer FESpaceUtility::CreateEmptyFESpace<1, double>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) FESpace<2, double>::Pointer FESpaceUtility::CreateEmptyFESpace<2, double>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) FESpace<3, double>::Pointer FESpaceUtility::CreateEmptyFESpace<3, double>(const std::string&);

} // end namespace Kratos
