//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TRANSLATION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TRANSLATION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "custom_utilities/trans/transformation.h"

namespace Kratos
{

/**
Represent a Translation in homogeneous coordinates
 */
template<typename TDataType>
class Translation : public Transformation<TDataType>
{
public:
    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(Translation);

    /// Type definitions
    typedef Transformation<TDataType> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

    /// Default constructor
    Translation(const TDataType& tX, const TDataType& tY, const TDataType& tZ) : BaseType()
    {
        BaseType::mTransMat(0, 3) = tX;
        BaseType::mTransMat(1, 3) = tY;
        BaseType::mTransMat(2, 3) = tZ;
    }

    /// Destructor
    virtual ~Translation() {}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Homogeneous Translation";
    }

};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const Translation<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TRANSLATION_H_INCLUDED

