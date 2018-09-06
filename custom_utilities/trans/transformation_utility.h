//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Aug 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/trans/transformation.h"

namespace Kratos
{

/**
Utility class to create useful transformation matrix.
 */
template<typename TDataType>
class TransformationUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TransformationUtility);

    /// Type definition
    typedef Transformation<TDataType> TransformationType;
    typedef typename TransformationType::VectorType VectorType;
    typedef typename TransformationType::MatrixType MatrixType;
    typedef typename TransformationType::IdentityMatrixType IdentityMatrixType;

    /// Default constructor
    TransformationUtility() {}

    /// Destructor
    virtual ~TransformationUtility() {}

    static MatrixType SkewSymmetric(const VectorType& v)
    {
        MatrixType M(3, 3);

        M(0, 0) = 0.0;
        M(0, 1) = -v(2);
        M(0, 2) = v(1);

        M(1, 0) = v(2);
        M(1, 1) = 0.0;
        M(1, 2) = -v(1);

        M(2, 0) = -v(1);
        M(2, 1) = v(0);
        M(2, 2) = 0.0;

        return M;
    }

    /// Create transformation matrix to transform vector a to vector b
    /// REF: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    static TransformationType CreateAlignTransformation(const VectorType& a, const VectorType& b)
    {
        MatrixType M = IdentityMatrixType(3, 3);
        TransformationType T;

        VectorType ua = a / norm_2(a);
        VectorType ub = b / norm_2(b);

        VectorType c = MathUtils<TDataType>::CrossProduct(ua, ub);
        double normc = norm_2(c);
        if (normc < 1.0e-10) return T;
        MatrixType ssc = SkewSymmetric(c);

        noalias(M) += ssc + prod(ssc, ssc) * (1.0 - inner_prod(ua, ub)) / pow(normc, 2);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                T(i, j) = M(i, j);

        return T;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TransformationUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
template<typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const TransformationUtility<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_UTILITY_H_INCLUDED defined

