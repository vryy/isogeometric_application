//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ROTATION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ROTATION_H_INCLUDED

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
struct ComputeRotationalTransformationMatrix_Helper
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<0, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle/180.0*M_PI);
        TDataType s = std::sin(angle/180.0*M_PI);
        trans_mat(1, 1) = c;
        trans_mat(1, 2) = -s;
        trans_mat(2, 1) = s;
        trans_mat(2, 2) = c;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<1, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle/180.0*M_PI);
        TDataType s = std::sin(angle/180.0*M_PI);
        trans_mat(0, 0) = c;
        trans_mat(0, 2) = s;
        trans_mat(2, 0) = -s;
        trans_mat(2, 2) = c;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<2, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle/180.0*M_PI);
        TDataType s = std::sin(angle/180.0*M_PI);
        trans_mat(0, 0) = c;
        trans_mat(0, 1) = -s;
        trans_mat(1, 0) = s;
        trans_mat(1, 1) = c;
    }
};

/**
 * Represent a Rotation in homogeneous coordinates
 * TAxis represent the axis of rotation
 */
template<int TAxis, typename TDataType>
class Rotation : public Transformation<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Rotation);

    /// Type definitions
    typedef Transformation<TDataType> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

    /// Default constructor
    Rotation(const TDataType& angle) : BaseType()
    {
        ComputeRotationalTransformationMatrix_Helper<TAxis, MatrixType, TDataType>::Execute(BaseType::mTransMat, angle);
    }

    /// Destructor
    virtual ~Rotation() {}

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Homogeneous Rotation";
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
inline std::ostream& operator <<(std::ostream& rOStream, const Rotation<TAxis, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ROTATION_H_INCLUDED

