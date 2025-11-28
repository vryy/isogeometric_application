//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2025 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/fespace.h"

namespace Kratos
{

/**
Utility class for FESpace.
 */
class FESpaceUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpaceUtility);

    /// Default constructor
    FESpaceUtility() {}

    /// Destructor
    virtual ~FESpaceUtility() {}

    /// Create the empty FESpace pointer from name
    template<int TDim, typename TLocalCoordinateType>
    static typename FESpace<TDim, TLocalCoordinateType>::Pointer CreateEmptyFESpace(const std::string& type);

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpaceUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const FESpaceUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_UTILITY_H_INCLUDED defined
