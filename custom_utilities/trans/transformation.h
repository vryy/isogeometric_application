//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/iga_define.h"


namespace Kratos
{

/**
Represent a transformation in homogeneous coordinates
 */
template<typename TDataType>
class Transformation
{
public:
    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(Transformation);

    /// Type definitions
    typedef boost::numeric::ublas::matrix<TDataType> MatrixType;
    typedef boost::numeric::ublas::vector<TDataType> VectorType;
    typedef boost::numeric::ublas::zero_matrix<TDataType> ZeroMatrixType;
    typedef boost::numeric::ublas::identity_matrix<TDataType> IdentityMatrixType;

    /// Default constructor
    Transformation()
    {
        mTransMat.resize(4, 4, false);
        noalias(mTransMat) = IdentityMatrixType(4, 4);
    }

    /// Constructor from location and three directional vectors
    /// B & N are the normal vectors, T is tangent vector
    Transformation(const VectorType& B, const VectorType& N, const VectorType& T, const VectorType& P)
    {
        mTransMat.resize(4, 4, false);
        noalias(column(mTransMat, 0)) = B;
        noalias(column(mTransMat, 1)) = N;
        noalias(column(mTransMat, 2)) = T;
        noalias(column(mTransMat, 3)) = P;
        mTransMat(3, 0) = 0.0;
        mTransMat(3, 1) = 0.0;
        mTransMat(3, 2) = 0.0;
        mTransMat(3, 3) = 1.0;
    }

    /// Constructor from location and three directional vectors
    /// B & N are the normal vectors, T is tangent vector
    Transformation(const array_1d<TDataType, 3>& B, const array_1d<TDataType, 3>& N, const array_1d<TDataType, 3>& T, const array_1d<TDataType, 3>& P)
    {
        mTransMat.resize(4, 4, false);
        noalias(column(mTransMat, 0)) = B;
        noalias(column(mTransMat, 1)) = N;
        noalias(column(mTransMat, 2)) = T;
        noalias(column(mTransMat, 3)) = P;
        mTransMat(3, 0) = 0.0;
        mTransMat(3, 1) = 0.0;
        mTransMat(3, 2) = 0.0;
        mTransMat(3, 3) = 1.0;
    }

    /// Constructor from location and two directional vectors. The remaining one is computed by cross product.
    Transformation(const VectorType& B, const VectorType& T, const VectorType& P)
    {
        VectorType N = MathUtils<TDataType>::CrossProduct(T, B);
        N *= 1.0/norm_2(N);

        mTransMat.resize(4, 4, false);
        noalias(column(mTransMat, 0)) = B;
        noalias(column(mTransMat, 1)) = N;
        noalias(column(mTransMat, 2)) = T;
        noalias(column(mTransMat, 3)) = P;
        mTransMat(3, 0) = 0.0;
        mTransMat(3, 1) = 0.0;
        mTransMat(3, 2) = 0.0;
        mTransMat(3, 3) = 1.0;
    }

    /// Constructor from location and two directional vectors. The remaining one is computed by cross product.
    Transformation(const array_1d<TDataType, 3>& B, const array_1d<TDataType, 3>& T, const array_1d<TDataType, 3>& P)
    {
        VectorType N = MathUtils<TDataType>::CrossProduct(T, B);
        N *= 1.0/norm_2(N);

        mTransMat.resize(4, 4);
        noalias(column(mTransMat, 0)) = B;
        noalias(column(mTransMat, 1)) = N;
        noalias(column(mTransMat, 2)) = T;
        noalias(column(mTransMat, 3)) = P;
        mTransMat(3, 0) = 0.0;
        mTransMat(3, 1) = 0.0;
        mTransMat(3, 2) = 0.0;
        mTransMat(3, 3) = 1.0;
    }

    /// Copy constructor
    Transformation(const Transformation& rOther)
    : mTransMat(rOther.mTransMat)
    {}

    /// Destructor
    virtual ~Transformation() {}

    /// Get the homogeneous transformation matrix
    const MatrixType& Mat() const {return mTransMat;}

    /// Prepend the transformation
    void PrependTransformation(const Transformation& rOther)
    {
        Matrix tmp = prod(this->mTransMat, rOther.Mat());
        noalias(this->mTransMat) = tmp;
    }

    /// Append the transformation
    void AppendTransformation(const Transformation& rOther)
    {
        Matrix tmp = prod(rOther.Mat(), this->mTransMat);
        noalias(this->mTransMat) = tmp;
    }

    /// Apply the transformation for a point in 3D
    template<typename TVectorType>
    void ApplyTransformation(TVectorType& value) const
    {
        TDataType new_value[3];
        for (std::size_t i = 0; i < 3; ++i)
        {
            new_value[i] = 0.0;
            for (std::size_t j = 0; j < 3; ++j)
                new_value[i] += mTransMat(i, j) * value[j];
            new_value[i] += mTransMat(i, 3);
        }
        for (std::size_t i = 0; i < 3; ++i)
            value[i] = new_value[i];
    }

    /// Compute the inverse of the transformation matrix
    /// REF: https://www.springer.com/cda/content/document/cda_downloaddocument/9783319197876-c2.pdf?SGWID=0-0-45-1514078-p177405814
    Transformation Inverse() const
    {
        Transformation trans;

        array_1d<TDataType, 3> s, n, a, r;
        this->V1(s);
        this->V2(n);
        this->V3(a);
        this->P(r);

        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                trans(i, j) = mTransMat(j, i);
        trans(0, 3) = -inner_prod(s, r);
        trans(1, 3) = -inner_prod(n, r);
        trans(2, 3) = -inner_prod(a, r);

        return trans;
    }

    /// Get the origin point
    void P(array_1d<TDataType, 3>& V) const { noalias(V) = subrange(column(mTransMat, 3), 0, 3); }
    array_1d<TDataType, 3> P() const
    {
        array_1d<TDataType, 3> V;
        this->P(V);
        return V;
    }

    /// Get the first vector
    void V1(array_1d<TDataType, 3>& V) const { noalias(V) = subrange(column(mTransMat, 0), 0, 3); }
    array_1d<TDataType, 3> V1() const
    {
        array_1d<TDataType, 3> V;
        this->V1(V);
        return V;
    }

    /// Get the second vector
    void V2(array_1d<TDataType, 3>& V) const { noalias(V) = subrange(column(mTransMat, 1), 0, 3); }
    array_1d<TDataType, 3> V2() const
    {
        array_1d<TDataType, 3> V;
        this->V2(V);
        return V;
    }

    /// Get the third vector
    void V3(array_1d<TDataType, 3>& V) const { noalias(V) = subrange(column(mTransMat, 2), 0, 3); }
    array_1d<TDataType, 3> V3() const
    {
        array_1d<TDataType, 3> V;
        this->V3(V);
        return V;
    }

    /// overload operator ()
    TDataType& operator() (const int& i, const int& j)
    {
        return mTransMat(i, j);
    }

    /// overload operator ()
    const TDataType& operator() (const int& i, const int& j) const
    {
        return mTransMat(i, j);
    }

    /// Assignment operator
    Transformation& operator=(const Transformation& rOther)
    {
        this->mTransMat = rOther.mTransMat;
        return *this;
    }

    /// Multiplication operator
    Transformation& operator*=(const Transformation& rOther)
    {
        Matrix tmp = prod(rOther.mTransMat, this->mTransMat);
        noalias(this->mTransMat) = tmp;
        return *this;
    }

    /// Multiplication operator
    friend Transformation operator*(const Transformation& t1, const Transformation& t2)
    {
        Transformation t = t1;
        t *= t2;
        return t;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Homogeneous Transformation";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mTransMat;
    }

protected:

    MatrixType mTransMat;

};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const Transformation<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TRANSFORMATION_H_INCLUDED

