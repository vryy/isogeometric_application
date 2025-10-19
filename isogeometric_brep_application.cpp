//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 26 Jan 2021 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "isogeometric_brep_application.h"

namespace Kratos
{

KratosIsogeometricBRepApplication::KratosIsogeometricBRepApplication()
{}

void KratosIsogeometricBRepApplication::Register()
{
    // calling base class to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosIsogeometricBRepApplication... " << std::endl;
}

}  // namespace Kratos.
