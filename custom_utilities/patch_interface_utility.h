//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2025 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
 * Utility class for patch interface handling and creation
 */
class KRATOS_API(ISOGEOMETRIC_APPLICATION) PatchInterfaceUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PatchInterfaceUtility);

    /// Default constructor
    PatchInterfaceUtility() {}

    /// Destructor
    virtual ~PatchInterfaceUtility() {}

    /// Helper function to create an empty patch interface based on name
    template<typename TPatchInterfaceType>
    static typename TPatchInterfaceType::Pointer CreateEmptyPatchInterface(const std::string& type);

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PatchInterfaceUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const PatchInterfaceUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_INTERFACE_UTILITY_H_INCLUDED defined
