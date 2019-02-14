//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Aug 2015 $
//   Revision:            $Revision: 1.0 $
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

/// Short class definition.
/** Detail class definition.
 */
class IsogeometricMathUtils
{
public:

    /// Pointer definition of IsogeometricMathUtils
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricMathUtils);

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
        ss << now->tm_mday << "/" << now->tm_mon+1 << "/" << (now->tm_year + 1900)
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
    //    Vector& Ubar,       // extended knot vector (OUTPUT)
    //    int& nt,            // relative location of the basis function w.r.t extended knot vector (OUTPUT)
    //    const std::vector<double>& Xi,   // local knot vector (INPUT)
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
        for(std::size_t i = 0; i < n; ++i)
        {
            if(Xi[i] == Xi[0])
                ++a;
            else
                break;
        }

        // compute the index of the basis function w.r.t the extended knot vector
        nt = p - a + 1;

        // count the multiplicity of the last knot
        int b = 0;
        for(std::size_t i = n - 1; i >= 0; --i)
        {
            if(Xi[i] == Xi[n-1])
                ++b;
            else
                break;
        }

        // compute the extended knot vector
        std::size_t len = nt + n + (p-b+1);
        if(Ubar.size() != len)
            Ubar.resize(len);

        for(std::size_t i = 0; i < nt; ++i)
            Ubar[i] = Xi[0];

        for(std::size_t i = nt; i < nt + n; ++i)
            Ubar[i] = Xi[i - nt];

        for(std::size_t i = nt + n; i < nt + n + (p-b+1); ++i)
            Ubar[i] = Xi[n-1];
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
        if(C.size1() != dimA1*dimB1 || C.size2() != dimA2*dimB2)
            C.resize(dimA1*dimB1, dimA2*dimB2);
        for(std::size_t i = 0; i < dimA1; ++i)
            for(std::size_t j = 0; j < dimA2; ++j)
                for(std::size_t k = 0; k < dimB1; ++k)
                    for(std::size_t l = 0; l < dimB2; ++l)
                        C(i*dimB1 + k, j*dimB2 + l) = A(i, j) * B(k, l);
    }


    /**
        Compute outer product of 2 vectors
     */
    template<class VectorType>
    static void outer_prod_vec(VectorType& C, const VectorType& A, const VectorType& B)
    {
        std::size_t dimA = A.size();
        std::size_t dimB = B.size();
        if(C.size() != dimA * dimB)
            C.resize(dimA * dimB);
        for(std::size_t i = 0; i < dimA; ++i)
            for(std::size_t j = 0; j < dimB; ++j)
                C(i * dimB + j) = A(i) * B(j);
    }


    /**
     * Convert a modified compressed sparse row matrix to compressed sparse row matrix B <- A
     * TODO check if compressed_matrix is returned
     */
    static Matrix MCSR2CSR(const Matrix& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);

        CompressedMatrix B(n, n);
        noalias( B ) = ZeroMatrix(n, n);

        for(unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if(col < i)
                {
                    B.push_back(i, col, A(1, A(0, i) + j));
                }
                else
                    break;
            }

            // insert the diagonal
            B.push_back(i, i, A(1, i));

            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if(col > i)
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
    static Matrix MCSR2MAT(const Matrix& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);

        Matrix B = ZeroMatrix(n, n);

        for(unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            B(i, i) = A(1, i);

            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
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
    static Matrix MAT2MCSR(const Matrix& A)
    {
        int n = A.size1();

        if(A.size2() != n)
            KRATOS_THROW_ERROR(std::logic_error, "The matrix needs to be square", "")

        std::vector<int> idx;
        std::vector<double> val;

        // firstly write the diagonal part of the matrix
        for(int i = 0; i < n; ++i)
        {
            val.push_back(A(i, i));
            idx.push_back(0);
        }

        // unused value
        val.push_back(0.0);

        // secondly traverse to off-diagonal element and write to idx and val
        int cnt = n + 1;
        for(int i = 0; i < n; ++i)
        {
            idx[i] = cnt;
            for(int j = 0; j < n; ++j)
            {
                if(j != i && A(i, j) != 0)
                {
                    val.push_back(A(i, j));
                    idx.push_back(j);
                    ++cnt;
                }
            }
        }

        idx[n] = cnt;

        Matrix B(2, cnt - 1);
        for(int i = 0; i < cnt - 1; ++i)
        {
            B(0, i) = static_cast<double>(idx[i] - 1);
            B(1, i) = val[i];
        }

        return B;
    }

    /**
     * Convert a triplet to compressed sparse row matrix
     */
    static Matrix Triplet2CSR(const Vector& rowPtr, const Vector& colInd, const Vector& values)
    {
        int m = rowPtr.size() - 1; // number of rows
        int n = *(std::max_element(colInd.begin(), colInd.end())) + 1; // number of columns
        return Triplet2CSR(m, n, rowPtr, colInd, values);
    }

    /**
     * Convert a triplet to compressed sparse row matrix
     */
    static Matrix Triplet2CSR(int m, int n, const Vector& rowPtr, const Vector& colInd, const Vector& values)
    {
        CompressedMatrix M(m, n);
        noalias(M) = ZeroMatrix(m, n);

        int i, j, nz, rowptr_this, rowptr_next, col;
        double val;
        for(i = 0; i < m; ++i)
        {
            rowptr_this = static_cast<int>(rowPtr[i]);
            rowptr_next = static_cast<int>(rowPtr[i + 1]);
            nz = rowptr_next - rowptr_this;
            for(j = 0; j < nz; ++j)
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


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "IsogeometricMathUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricMathUtils";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

protected:

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
inline std::istream& operator >>(std::istream& rIStream, IsogeometricMathUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const IsogeometricMathUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_MATH_UTILS_H_INCLUDED  defined
