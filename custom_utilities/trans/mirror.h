//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2019 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MIRROR_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MIRROR_H_INCLUDED

// System includes
#include <math.h>
#include <cmath>
#include <string>

// External includes

// Project includes
#include "custom_utilities/trans/transformation.h"

namespace Kratos
{

template<int TAxis, typename TMatrixType, typename TDataType>
struct ComputeMirrorTransformationMatrix_Helper
{
    static void Execute(TMatrixType& trans_mat)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeMirrorTransformationMatrix_Helper<0, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat)
    {
        trans_mat(0, 0) = -1;
        trans_mat(1, 1) = 1;
        trans_mat(2, 2) = 1;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeMirrorTransformationMatrix_Helper<1, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat)
    {
        trans_mat(0, 0) = 1;
        trans_mat(1, 1) = -1;
        trans_mat(2, 2) = 1;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeMirrorTransformationMatrix_Helper<2, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat)
    {
        trans_mat(0, 0) = 1;
        trans_mat(1, 1) = 1;
        trans_mat(2, 2) = -1;
    }
};

/**
 * Represent a Mirror in homogeneous coordinates
 * TAxis represent the axis of Mirror
 */
template<int TAxis, typename TDataType>
class Mirror : public Transformation<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Mirror);

    /// Type definitions
    typedef Transformation<TDataType> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

    /// Default constructor
    Mirror() : BaseType()
    {
        ComputeMirrorTransformationMatrix_Helper<TAxis, MatrixType, TDataType>::Execute(BaseType::mTransMat);
    }

    /// Destructor
    virtual ~Mirror() {}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Homogeneous Mirror";
        if (TAxis == 0)
            rOStream << "_X";
        else if (TAxis == 1)
            rOStream << "_Y";
        else if (TAxis == 2)
            rOStream << "_Z";
    }

};

/// output stream function
template<int TAxis, class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const Mirror<TAxis, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MIRROR_H_INCLUDED

