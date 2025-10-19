//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Aug 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_MATH_UTILS_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_MATH_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <ctime>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
 * Math operations for IGA
 */
template<class TDataType> class IsogeometricMathUtils
{
public:

    /// Pointer definition of IsogeometricMathUtils
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricMathUtils);

    typedef boost::numeric::ublas::matrix<TDataType> MatrixType;

    typedef boost::numeric::ublas::vector<TDataType> VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    /// Default constructor.
    IsogeometricMathUtils()
    {}

    /// Destructor.
    virtual ~IsogeometricMathUtils()
    {}

    static std::string current_time()
    {
        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        std::stringstream ss;
        ss << now->tm_mday << "/" << now->tm_mon + 1 << "/" << (now->tm_year + 1900)
           << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec;
        return ss.str();
    }

    static void timestamp(std::ostream& rOStream)
    {
        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        rOStream << "//This file is created on " << current_time() << std::endl << std::endl;
    }

    // void IsogeometricMathUtils::compute_extended_knot_vector(
    //    VectorType& Ubar,       // extended knot vector (OUTPUT)
    //    int& nt,            // relative location of the basis function w.r.t extended knot vector (OUTPUT)
    //    const std::vector<TDataType>& Xi,   // local knot vector (INPUT)
    //    const int p)        // degree of the basis function (INPUT)
    // {
    //     // count the multiplicity of the first knot
    //     int n = Xi.size();
    //     int a = 0;
    //     for(std::size_t i = 0; i < n; ++i)
    //     {
    //        if(Xi[i] == Xi[0])
    //            ++a;
    //        else
    //            break;
    //     }

    //     // compute the index of the basis function w.r.t the extended knot vector
    //     nt = p - a + 1;

    //     // count the multiplicity of the last knot
    //     int b = 0;
    //     for(std::size_t i = n - 1; i >= 0; --i)
    //     {
    //        if(Xi[i] == Xi[n-1])
    //            ++b;
    //        else
    //            break;
    //     }

    //     // compute the extended knot vector
    //     Ubar.resize(nt + n + (p-b+1));

    //     for(std::size_t i = 0; i < nt; ++i)
    //        Ubar[i] = Xi[0];

    //     for(std::size_t i = nt; i < nt + n; ++i)
    //        Ubar[i] = Xi[i - nt];

    //     for(std::size_t i = nt + n; i < nt + n + (p-b+1); ++i)
    //        Ubar[i] = Xi[n-1];
    // }

    /**
     Compute the extended knot vector for a local knot vector
     */
    template<class VectorType, class ValuesContainerType>
    static void compute_extended_knot_vector(
        VectorType& Ubar,               // extended knot vector (OUTPUT)
        int& nt,                        // relative location of the basis function w.r.t extended knot vector (OUTPUT)
        const ValuesContainerType& Xi,  // local knot vector (INPUT)
        const int p)                    // degree of the basis function (INPUT)
    {
        // count the multiplicity of the first knot
        int n = Xi.size();
        int a = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (Xi[i] == Xi[0])
            {
                ++a;
            }
            else
            {
                break;
            }
        }

        // compute the index of the basis function w.r.t the extended knot vector
        nt = p - a + 1;

        // count the multiplicity of the last knot
        int b = 0;
        for (std::size_t i = n - 1; i >= 0; --i)
        {
            if (Xi[i] == Xi[n - 1])
            {
                ++b;
            }
            else
            {
                break;
            }
        }

        // compute the extended knot vector
        std::size_t len = nt + n + (p - b + 1);
        if (Ubar.size() != len)
        {
            Ubar.resize(len);
        }

        for (std::size_t i = 0; i < nt; ++i)
        {
            Ubar[i] = Xi[0];
        }

        for (std::size_t i = nt; i < nt + n; ++i)
        {
            Ubar[i] = Xi[i - nt];
        }

        for (std::size_t i = nt + n; i < nt + n + (p - b + 1); ++i)
        {
            Ubar[i] = Xi[n - 1];
        }
    }

    /**
        Compute outer product of 2 matrices
     */
    template<class MatrixType>
    static void outer_prod_mat(MatrixType& C, const MatrixType& A, const MatrixType& B)
    {
        std::size_t dimA1 = A.size1();
        std::size_t dimA2 = A.size2();
        std::size_t dimB1 = B.size1();
        std::size_t dimB2 = B.size2();
        if (C.size1() != dimA1 * dimB1 || C.size2() != dimA2 * dimB2)
        {
            C.resize(dimA1 * dimB1, dimA2 * dimB2);
        }
        for (std::size_t i = 0; i < dimA1; ++i)
            for (std::size_t j = 0; j < dimA2; ++j)
                for (std::size_t k = 0; k < dimB1; ++k)
                    for (std::size_t l = 0; l < dimB2; ++l)
                    {
                        C(i * dimB1 + k, j * dimB2 + l) = A(i, j) * B(k, l);
                    }
    }

    /**
        Compute outer product of 2 vectors
     */
    template<class VectorType>
    static void outer_prod_vec(VectorType& C, const VectorType& A, const VectorType& B)
    {
        std::size_t dimA = A.size();
        std::size_t dimB = B.size();
        if (C.size() != dimA * dimB)
        {
            C.resize(dimA * dimB);
        }
        for (std::size_t i = 0; i < dimA; ++i)
            for (std::size_t j = 0; j < dimB; ++j)
            {
                C(i * dimB + j) = A(i) * B(j);
            }
    }

    /**
     * Convert a modified compressed sparse row matrix to compressed sparse row matrix B <- A
     * TODO check if compressed_matrix is returned
     */
    static MatrixType MCSR2CSR(const MatrixType& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);

        CompressedMatrix B(n, n);
        noalias( B ) = ZeroMatrix(n, n);

        for (unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            // traverse left for off-diagonal entries
            for (unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if (col < i)
                {
                    B.push_back(i, col, A(1, A(0, i) + j));
                }
                else
                {
                    break;
                }
            }

            // insert the diagonal
            B.push_back(i, i, A(1, i));

            // traverse left for off-diagonal entries
            for (unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if (col > i)
                {
                    B.push_back(i, col, A(1, A(0, i) + j));
                }
            }
        }

        B.complete_index1_data();

        return B;
    }

    /**
     * Convert a modified compressed sparse row matrix to trivial matrix
     */
    static MatrixType MCSR2MAT(const MatrixType& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);

        MatrixType B = ZeroMatrix(n, n);

        for (unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            B(i, i) = A(1, i);

            // traverse left for off-diagonal entries
            for (unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                B(i, col) = A(1, A(0, i) + j);
            }
        }

        return B;
    }

    /**
     * Convert a trivial matrix to modified compressed sparse row matrix
     */
    static MatrixType MAT2MCSR(const MatrixType& A)
    {
        int n = A.size1();

        if (A.size2() != n)
            KRATOS_ERROR << "The matrix is not square";

        std::vector<int> idx;
        std::vector<TDataType> val;

        // firstly write the diagonal part of the matrix
        for (int i = 0; i < n; ++i)
        {
            val.push_back(A(i, i));
            idx.push_back(0);
        }

        // unused value
        val.push_back(0.0);

        // secondly traverse to off-diagonal element and write to idx and val
        int cnt = n + 1;
        for (int i = 0; i < n; ++i)
        {
            idx[i] = cnt;
            for (int j = 0; j < n; ++j)
            {
                if (j != i && A(i, j) != 0)
                {
                    val.push_back(A(i, j));
                    idx.push_back(j);
                    ++cnt;
                }
            }
        }

        idx[n] = cnt;

        MatrixType B(2, cnt - 1);
        for (int i = 0; i < cnt - 1; ++i)
        {
            B(0, i) = static_cast<TDataType>(idx[i] - 1);
            B(1, i) = val[i];
        }

        return B;
    }

    /**
     * Convert a triplet to compressed sparse row matrix
     */
    static MatrixType Triplet2CSR(const VectorType& rowPtr, const VectorType& colInd, const VectorType& values)
    {
        int m = rowPtr.size() - 1; // number of rows
        int n = *(std::max_element(colInd.begin(), colInd.end())) + 1; // number of columns
        return Triplet2CSR(m, n, rowPtr, colInd, values);
    }

    /**
     * Convert a triplet to compressed sparse row matrix
     */
    static MatrixType Triplet2CSR(int m, int n, const VectorType& rowPtr, const VectorType& colInd, const VectorType& values)
    {
        CompressedMatrix M(m, n);
        noalias(M) = ZeroMatrix(m, n);

        int i, j, nz, rowptr_this, rowptr_next, col;
        TDataType val;
        for (i = 0; i < m; ++i)
        {
            rowptr_this = static_cast<int>(rowPtr[i]);
            rowptr_next = static_cast<int>(rowPtr[i + 1]);
            nz = rowptr_next - rowptr_this;
            for (j = 0; j < nz; ++j)
            {
                col = static_cast<int>(colInd[rowptr_this + j]);
                val = values[rowptr_this + j];
                M.push_back(i, col, val);
//                M(row_this, col) = val;
            }
        }

        M.complete_index1_data();
        return M;
    }

    /**
     * Compute the intersection between two lines in 2D
     */
    static int ComputeIntersection2D(const TDataType a1, const TDataType b1, const TDataType c1,
                                     const TDataType a2, const TDataType b2, const TDataType c2,
                                     TDataType& x, TDataType& y, const TDataType TOL)
    {
        const TDataType det = a1*b2 - a2*b1;
        if (std::abs(det) < TOL)
            return -1; // the two lines are parallel

        x = (b1*c2 - b2*c1) / det;
        y = (c1*a2 - c2*a1) / det;

        return 0;
    }

    /**
     * Compute the intersection between two lines in 2D
     */
    template<typename TPointType>
    static int ComputeIntersection2D(const TPointType& P0, const TPointType& P1,
            const TPointType& P2, const TPointType& P3, TPointType& P, const TDataType TOL)
    {
        const TDataType a1 = P1[1] - P0[1];
        const TDataType b1 = -P1[0] + P0[0];
        const TDataType c1 = -a1*P0[0] - b1*P0[1];

        const TDataType a2 = P3[1] - P2[1];
        const TDataType b2 = -P3[0] + P2[0];
        const TDataType c2 = -a2*P2[0] - b2*P2[1];

        TDataType x, y;
        int error_code = ComputeIntersection2D(a1, b1, c1, a2, b2, c2, x, y, TOL);
        if (error_code != 0)
            return error_code;

        P[0] = x;
        P[1] = y;
        P[2] = 0.0;

        TDataType d1, d2, d;
        int i1 = 0, i2 = 0;
        d1 = sqrt(pow(x - P0[0], 2) + pow(y - P0[1], 2));
        d2 = sqrt(pow(x - P1[0], 2) + pow(y - P1[1], 2));
        d = sqrt(pow(P1[0] - P0[0], 2) + pow(P1[1] - P0[1], 2));
        if (abs(d1 + d2 - d) < TOL) i1 = 1;
        d1 = sqrt(pow(x - P2[0], 2) + pow(y - P2[1], 2));
        d2 = sqrt(pow(x - P3[0], 2) + pow(y - P3[1], 2));
        d = sqrt(pow(P3[0] - P2[0], 2) + pow(P3[1] - P2[1], 2));
        if (abs(d1 + d2 - d) < TOL) i2 = 1;
        if (i1 == 1 && i2 == 1) return 0; // the intersection is in the middle of both P0-P1 and P2-P3
        if (i1 == 1 && i2 == 0) return 1; // the intersection is in the middle of P0-P1 (line 1) and not P2-P3
        if (i1 == 0 && i2 == 1) return 2; // the intersection is in the middle of P2-P3 (line 2) and not P0-P1
        KRATOS_ERROR << "Can't go here. Something's wrong";
        return -1; // can't go here
    }

    /**
     * Compute the projection of a point on the line in 2D
     */
    static int ComputeProjection2D(const TDataType x0, const TDataType y0,
                                     const TDataType a, const TDataType b, const TDataType c,
                                     TDataType& x, TDataType& y, const TDataType TOL)
    {
        const double norm = a*a + b*b;
        if (norm < TOL)
            return -1; // the line is singular

        x = (-b*(a*y0 - b*x0) - c*a) / norm;
        y = (a*(a*y0 - b*x0) - c*b) / norm;
        return 0;
    }

    /**
     * Compute the projection of a point on the line in 2D
     */
    template<typename TPointType>
    static int ComputeProjection2D(const TPointType& P0,
            const TPointType& P1, const TPointType& P2,
            TPointType& P, const TDataType TOL)
    {
        const TDataType a = P2[1] - P1[1];
        const TDataType b = -P2[0] + P1[0];
        const TDataType c = -a*P1[0] - b*P1[1];

        TDataType x, y;
        int error_code = ComputeProjection2D(P0[0], P0[1], a, b, c, x, y, TOL);
        if (error_code != 0)
            return error_code;

        P[0] = x;
        P[1] = y;
        P[2] = 0.0;

        TDataType d1, d2, d;
        int i = 0;
        d1 = sqrt(pow(x - P1[0], 2) + pow(y - P1[1], 2));
        d2 = sqrt(pow(x - P2[0], 2) + pow(y - P2[1], 2));
        d = sqrt(pow(P2[0] - P1[0], 2) + pow(P2[1] - P1[1], 2));
        if (abs(d1 + d2 - d) < TOL)
            return 0; // the projection is in the middle of both P1-P2
        else
            return 1; // the projection is outside P1-P2
        KRATOS_ERROR << "Can't go here. Something's wrong";
        return -1; // can't go here
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "IsogeometricMathUtils";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

private:

    /// Assignment operator.
    IsogeometricMathUtils& operator=(IsogeometricMathUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricMathUtils(IsogeometricMathUtils const& rOther)
    {
    }

}; // Class IsogeometricMathUtils

/// input stream function
template<class TDataType>
inline std::istream& operator >>(std::istream& rIStream, IsogeometricMathUtils<TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const IsogeometricMathUtils<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_MATH_UTILS_H_INCLUDED  defined
