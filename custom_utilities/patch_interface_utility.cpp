//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2025 $
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/nurbs/bsplines_patch_interface.h"
#include "custom_utilities/patch_interface_utility.h"

namespace Kratos
{

template<typename TPatchInterfaceType>
typename TPatchInterfaceType::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface(const std::string& type)
{
    typedef PatchInterface<TPatchInterfaceType::Dim,
                typename TPatchInterfaceType::LocalCoordinateType,
                typename TPatchInterfaceType::CoordinateType,
                typename TPatchInterfaceType::DataType> PatchInterfaceType;

    typedef BSplinesPatchInterface<TPatchInterfaceType::Dim,
                typename TPatchInterfaceType::LocalCoordinateType,
                typename TPatchInterfaceType::CoordinateType,
                typename TPatchInterfaceType::DataType> BSplinesPatchInterfaceType;
    if (type == BSplinesPatchInterfaceType::StaticType())
    {
        return BSplinesPatchInterfaceType::Create();
    }
    else
        return PatchInterfaceType::Create();
}

/// template function instantiation

template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<1>::RealPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<1>::RealPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<2>::RealPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<2>::RealPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<3>::RealPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<3>::RealPatchInterface>(const std::string&);

template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<1>::ComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<1>::ComplexPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<2>::ComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<2>::ComplexPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<3>::ComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<3>::ComplexPatchInterface>(const std::string&);

template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<1>::GComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<1>::GComplexPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<2>::GComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<2>::GComplexPatchInterface>(const std::string&);
template KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchSelector<3>::GComplexPatchInterface::Pointer PatchInterfaceUtility::CreateEmptyPatchInterface<PatchSelector<3>::GComplexPatchInterface>(const std::string&);

} // namespace Kratos.
