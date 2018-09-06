//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-08-18 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BSPLINE_UTILS_H_INCLUDED )
#define  KRATOS_BSPLINE_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/isogeometric_math_utils.h"

#define FUNCTIONALITY_CHECK // this macro is used to enable the compulsory functionality check of the program. If we are sure if it works, then we can disable it.

namespace Kratos
{
///@addtogroup IsogeometricApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * Utility for various operations with Bsplines
 */
class BSplineUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BSplineUtils
    KRATOS_CLASS_POINTER_DEFINITION(BSplineUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BSplineUtils()
    {}

    /// Destructor.
    virtual ~BSplineUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // original implementation by hbui
//    template<class ValuesContainerType>
//      static int FindSpan(
//        const int rN,
//        const int rP,
//        const double rXi,
//        const ValueContainerType& rU
//        )
//      {
//        int size = rU.size();
//
//        if(size != rN + rP + 1)
//        {
//            KRATOS_THROW_ERROR(std::logic_error, "Size of rU must be n + p + 1", __FUNCTION__);
//        }
//
//        double tol = 1e-12;
//
//        if(fabs(rXi - 1.0) < tol)
//            return rN - 1;
//
//        int span_a = 0;
//        int span_b = rN + rP;
//        int span_mid;
//
//        while( abs(span_a - span_b) != 1 )
//        {
//            span_mid = (span_a + span_b) / 2;
//
//            if(rXi < rU[span_mid])
//            {
//                span_b = span_mid;
//            }
//            else
//            {
//                span_a = span_mid;
//            }
//        }
//
//        return span_a;
//      }

    //implementation from the NURBS book
    // Note:
    // Algorithm A2.1 from 'The NURBS BOOK' pg68
    // as that algorithm only works for nonperiodic
    // knot vectors, nonetheless the results should
    // be EXACTLY the same if U is nonperiodic
    // sample: rN = 4, rP = 2, rU = [0 0 0 0.5 1 1 1], rXi = 0.1 (rU.size() == rN + rP + 1)
    // note: this only works with open knot vector (i.e. 0  and 1 is repeated p + 1 times)
    template<class ValuesContainerType>
    static int FindSpan(
            const int& rN,
            const int& rP,
            const double& rXi,
            const ValuesContainerType& rU
    )
    {
        if(rXi == rU[rN]) return rN - 1;

        int low = rP;
        int high = rN;
        int mid = (low + high) / 2;

        while( rXi < rU[mid] || rXi >= rU[mid+1] )
        {
            if(rXi < rU[mid])
            {
                high = mid;
            }
            else
            {
                low = mid;
            }
            mid = (low + high) / 2;
        }

        return mid;
    }

    // implementation in GeoPde, low_level_functions.cc
    // Note: this implementation has linear, rather than log complexity
    template<class ValuesContainerType>
    static int FindSpan2(
            const int rN,
            const int rP,
            const double rXi,
            const ValuesContainerType& rU
    )
    {
        int ret = 0;
//        KRATOS_WATCH(rN)
        while ((ret++ < rN) && (rU[ret] <= rXi))
        {
//            KRATOS_WATCH(ret)
        }
        if(ret == rN+1) --ret;
//        std::cout << "--------" << std::endl;
        return (ret - 1);
    }

    /// Find the span of knot in the local knot vector
    /// Remarks: it will give the based-1 index
    ///          it only works if U is non-repeated
    template<typename TDataType, class TValuesContainerType>
    static std::size_t FindSpanLocal(const TDataType& Xi, const TValuesContainerType& U)
    {
        //        int low = 0;
//        int high = U.size() - 1;
//        int mid = (low + high) / 2;
//
//        while( Xi < U[mid] || Xi >= U[mid+1] )
//        {
//            if(Xi < U[mid])
//            {
//                high = mid;
//            }
//            else
//            {
//                low = mid;
//            }
//            mid = (low + high) / 2;
//        }
//
//        return mid + 1;

        if(!U.empty())
        {
            if(Xi < U[0])
                return 0;

            if(Xi > U[U.size()-1])
                return U.size();

            for(std::size_t i = 0; i < U.size()-1; ++i)
                if(Xi >= U[i] && Xi < U[i + 1])
                    return i + 1;
        }

        return 0;
    }

    // Remark: this function only works correctly when rXi is in the knot span rI.
//    static void BasisFuns(
//            ValueContainerType& rS,
//            const int rI,
//            const double rXi,
//            const int rP,
//            const ValueContainerType& rU
//    )
//    {
//        rS[0] = 1.0;
//
//        for(unsigned int j = 1; j < rP + 1; ++j)
//        {
//            rS[j] = (rXi - rU[rI]) / (rU[rI + j] - rU[rI]) * rS[j-1];
//
//            for(unsigned int i = j - 1; i > 0; --i)
//            {
//                int t = j - i;
//                rS[i] = (rXi - rU[rI - t]) / (rU[rI - t + j] - rU[rI - t]) * rS[i - 1]
//                + (rU[rI - t + j + 1] - rXi) / (rU[rI - t + j + 1] - rU[rI - t + 1]) * rS[i];
//            }
//
//            rS[0] *= (rU[rI + 1] - rXi) / (rU[rI + 1] - rU[rI - j + 1]);
//        }
//    }
    //time:
    //N = 100: 0.000123978
    //N = 1000: 0.00120997
    //N = 10000: 0.013212
    //N = 100000: 0.119535
    //N = 1000000: 0.96382
    //N = 10000000: 8.45895
    //N = 100000000: 81.9643, 81.9494

    // This function works slightly faster than the above implementation
    template<class ValuesContainerType1, class ValuesContainerType2>
    static void BasisFuns(ValuesContainerType1& rS,
                          const int& rI,
                          const double& rXi,
                          const int& rP,
                          const ValuesContainerType2& rU)
    {
        unsigned int j, r;
        double saved, temp;

        double* left = new double[rP + 1];
        double* right = new double[rP + 1];

        std::fill(left, left + rP + 1, 0.0);
        std::fill(right, right + rP + 1, 0.0);

        rS[0] = 1.0;
        for (j = 1; j <= rP; ++j)
        {
            left[j] = rXi - rU[rI + 1 - j];
            right[j] = rU[rI + j] - rXi;
            saved = 0.0;

            for (r = 0; r < j; ++r)
            {
                temp = rS[r] / (right[r + 1] + left[j - r]);
                rS[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }

            rS[j] = saved;
        }

        delete left, right;
    }
    //N = 10000000: 8.18893
    //N = 100000000: 76.7523, 81.9167

    struct MatrixOp
    {
        void InitZero(Matrix& S, const int& m, const int& n) const
        {
            if (S.size1() != m || S.size2() != n)
                S.resize(m, n, false);
            noalias(S) = ZeroMatrix(m, n);
        }

        double& Get(Matrix& S, const int& i, const int& j) const
        {
            return S(i, j);
        }
    };

    template<typename TDataType>
    struct StdVector2DOp
    {
        typedef std::vector<std::vector<TDataType> > Matrix;

        void InitZero(Matrix& S, const int& m, const int& n) const
        {
            if (S.size() != m) S.resize(m);
            for (int i = 0; i < m; ++i)
            {
                if (S[i].size() != n) S[i].resize(n);
                for (int j = 0; j < n; ++j) S[i][j] = 0.0;
            }
        }

        double& Get(Matrix& S, const int& i, const int& j) const
        {
            return S[i][j];
        }
    };

    /**
     * Computes b-spline function derivatives
     */
    template<class ValuesContainerType1, class ValuesContainerType1Operator, class ValuesContainerType2>
    static void BasisFunsDer(ValuesContainerType1& rS,
                             int rI,
                             const double& rXi,
                             const int& rP,
                             const ValuesContainerType2& rU,
                             const int& rD,
                             const ValuesContainerType1Operator& Op) // use MatrixOp() or StdVector2DOp()
    {
        int j, r, k, rk, pk, j1, j2, s1, s2;
        double d, temp, saved;

//        ders = zeros(nders+1,pl+1);
//        noalias(rS) = ZeroMatrix(rD + 1, rP + 1);
        Op.InitZero(rS, rD + 1, rP + 1);

//        ndu = zeros(pl+1,pl+1);
        Matrix ndu = ZeroMatrix(rP + 1, rP + 1);
//        left = zeros(pl+1);
        Vector left = ZeroVector(rP + 1);
//        right = zeros(pl+1);
        Vector right = ZeroVector(rP + 1);

//        a = zeros(2,pl+1);
        Matrix a = ZeroMatrix(2, rP + 1);

//        ndu(1,1) = 1;
        ndu(0, 0) = 1.0;
        rI += 1;

//        for j = 1:pl
        for (j = 1; j <= rP; ++j)
        {
//            left(j+1) = u - u_knotl(i+1-j);
            left(j) = rXi - rU[rI - j];
//            right(j+1) = u_knotl(i+j) - u;
            right(j) = rU[rI + j - 1] - rXi;
//            saved = 0;
            saved = 0.0;
//            for r = 0:j-1
            for (r = 0; r <= j - 1; ++r)
            {
//                ndu(j+1,r+1) = right(r+2) + left(j-r+1);
                ndu(j, r) = right(r + 1) + left(j - r);
//                temp = ndu(r+1,j)/ndu(j+1,r+1);
                temp = ndu(r, j - 1) / ndu(j, r);
//                ndu(r+1,j+1) = saved + right(r+2)*temp;
                ndu(r, j) = saved + right(r + 1) * temp;
//                saved = left(j-r+1)*temp;
                saved = left(j - r) * temp;
//            end
            }
//            ndu(j+1,j+1) = saved;
            ndu(j, j) = saved;
//        end
        }

//        for j = 0:pl
        for (j = 0; j <= rP; ++j)
        {
//            ders(1,j+1) = ndu(j+1,pl+1);
//            rS(0, j) = ndu(j, rP);
            Op.Get(rS, 0, j) = ndu(j, rP);
//        end
        }

//        for r = 0:pl
        for (r = 0; r <= rP; ++r)
        {
            s1 = 0;
            s2 = 1;
//            a(1,1) = 1;
            a(0, 0) = 1.0;
//            for k = 1:nders //compute kth derivative
            for (k = 1; k <= rD; ++k)
            {
                d = 0.0;
                rk = r - k;
                pk = rP - k;
                if(r >= k)
                {
//                    a(s2+1 , 1) = a(s1+1,1)/ndu(pk+2,rk+1);
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
//                    d = a(s2+1,1)*ndu(rk+1,pk+1);
                    d = a(s2, 0) * ndu(rk, pk);
                }
                if(rk >= -1)
                {
                    j1 = 1;
                }
                else
                {
                    j1 = -rk;
                }
                if((r-1) <= pk)
                {
                    j2 = k - 1;
                }
                else
                {
                    j2 = rP - r;
                }
//                for j = j1:j2
                for (j = j1; j <= j2; ++j)
                {
//                    a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
//                    d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
                    d = d + a(s2, j) * ndu(rk + j, pk);
//                end
                }
                if (r <= pk)
                {
//                    a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                    a(s2, k) = -a(s1, k - 1) / ndu(pk+1, r);
//                    d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
                    d = d + a(s2, k) * ndu(r, pk);
//                end
                }
//                ders(k+1,r+1) = d;
//                rS(k, r) = d;
                Op.Get(rS, k, r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
//            end
            }
//        end
        }

        r = rP;
//        for k = 1:nders
        for (k = 1; k <= rD; ++k)
        {
//            for j = 0:pl
            for (j = 0; j <= rP; ++j)
            {
//                ders(k+1,j+1) = ders(k+1,j+1)*r;
//                rS(k, j) = rS(k, j) * r;
                Op.Get(rS, k, j) *= r;
//            end
            }
            r = r * (rP - k);
//        end
        }
    }

    /// Compute the B-spline basis function based on Cox-de-Boor algorithm
    //    % Input:
    //    %   u       knot to be compute the function value
    //    %   i       index of the basis function; for local knot vector, it must be 0
    //    %   p       B-spline degree
    //    %   knots   local knot vector, must be ascending
    //    % Output: function value
    template<class ValuesContainerType>
    static double CoxDeBoor(const double& u, const int& i, const int& p, const ValuesContainerType& knots)
    {
        double y = 0.0;

        if (p == 0)
        {
            if (u >= knots[i] && u < knots[i+1])
                y = 1.0;
            else
                y = 0.0;

            if (i+2 <= knots.size()-1)
            {
                if (u >= knots[i+1] && knots[i+2] == knots[i+1])
                    y = 1.0;
            }

            return y;
        }

        if (knots[i+p] > knots[i])
            y = y + (u - knots[i]) / (knots[i+p] - knots[i]) * CoxDeBoor(u, i, p-1, knots);

        if (knots[i+p+1] > knots[i+1])
            y = y + (knots[i+p+1] - u) / (knots[i+p+1] - knots[i+1]) * CoxDeBoor(u, i+1, p-1, knots);

        return y;
    }

    /// Compute the refinement coefficients for one knot insertion B-Splines refinement in 1D
    /// REF: Eq (5.10) the NURBS books
    template<class MatrixType, class ValuesContainerType>
    static void ComputeBsplinesKnotInsertionCoefficients1D(MatrixType& D,
                                                           ValuesContainerType& new_knots,
                                                           const int& p,
                                                           const ValuesContainerType& knots,
                                                           const double& k)
    {
        // compute the number of basis function
        int n = knots.size() - p - 1;

        // find the span of the inserted knot
        int s = FindSpan(n, p, k, knots);
//        std::cout << "knots:";
//        for(int i = 0; i < knots.size(); ++i)
//            std::cout << " " << knots[i];
//        std::cout << std::endl;
//        std::cout << "k=" << k << ",s=" << s << std::endl;

        // form the new knot vector
        new_knots.resize(knots.size() + 1);
        for(int i = 0; i < s+1; ++i)
            new_knots[i] = knots[i];
        new_knots[s+1] = k;
        for(int i = s+2; i < knots.size() + 1; ++i)
            new_knots[i] = knots[i-1];

//        std::cout << "new_knots:";
//        for(int i = 0; i < new_knots.size(); ++i)
//            std::cout << " " << new_knots[i];
//        std::cout << std::endl;

        // initialize and compute coefficient matrix
        D.resize(n, n+1);
        noalias(D) = ZeroMatrix(n, n+1);
        for(int i = 0; i < s-p; ++i)
            D(i, i) = 1.0;
        for(int i = s-p; i < s+1; ++i)
        {
            D(i, i) = (k - new_knots[i]) / (new_knots[i+p+1] - new_knots[i]);
            D(i, i+1) = (new_knots[i+p+2] - k) / (new_knots[i+p+2] - new_knots[i+1]);
        }
        for(int i = s+1; i < n; ++i)
            D(i, i+1) = 1.0;
//        KRATOS_WATCH(D)
//        std::cout << "-----------------" << std::endl;
    }

    /// Compute the refinement coefficients for multiple knots insertion B-Splines refinement in 1D
    /// The algorithm used here is not efficient, since multiple matrix multiplication is performed. I shall use the algorithm in the NURBS book.
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients1D(MatrixType& D,
                                                           ValuesContainerType& new_knots,
                                                           const int& p,
                                                           const ValuesContainerType2& knots,
                                                           const ValuesContainerType3& ins_knots)
    {
        // compute the number of basis function
        int n = knots.size() - p - 1;

        // initialize coefficient matrix
        D.resize(n, n);
        noalias(D) = IdentityMatrix(n, n);

        // copy the knot vector
        new_knots.resize(knots.size());
        // std::copy(knots.begin(), knots.end(), new_knots.begin());
        for (std::size_t i = 0; i < knots.size(); ++i) new_knots[i] = knots[i];

        ValuesContainerType tmp_knots;
        tmp_knots.resize(new_knots.size());
        std::copy(new_knots.begin(), new_knots.end(), tmp_knots.begin());

        // insert individual knots and concatenate coefficient matrix
        MatrixType Di, Dt;
        for(std::size_t i = 0; i < ins_knots.size(); ++i)
        {
            ComputeBsplinesKnotInsertionCoefficients1D(Di, new_knots, p, tmp_knots, ins_knots[i]);

            Dt.resize(D.size1(), Di.size2());
            noalias(Dt) = prod(D, Di);

            D.resize(Dt.size1(), Dt.size2());
            noalias(D) = Dt;

            tmp_knots.resize(new_knots.size());
            std::copy(new_knots.begin(), new_knots.end(), tmp_knots.begin());
        }
//        KRATOS_WATCH(D)
    }

    /// Compute the refinement coefficients for multiple knots insertion B-Splines refinement in 2D
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients2D(MatrixType& D,
                                                           ValuesContainerType& new_knots1,
                                                           ValuesContainerType& new_knots2,
                                                           const int& p1,
                                                           const int& p2,
                                                           const ValuesContainerType2& knots1,
                                                           const ValuesContainerType2& knots2,
                                                           const ValuesContainerType3& ins_knots1,
                                                           const ValuesContainerType3& ins_knots2)
    {
        MatrixType D1, D2;
        ComputeBsplinesKnotInsertionCoefficients1D(D1, new_knots1, p1, knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1D(D2, new_knots2, p2, knots2, ins_knots2);
//        KRATOS_WATCH(D1)
//        KRATOS_WATCH(D2)
        IsogeometricMathUtils::outer_prod_mat(D, D2, D1);
//        KRATOS_WATCH(D)
    }

    /// Compute the refinement coefficients for multiple knots insertion B-Splines refinement in 3D
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients3D(MatrixType& D,
                                                           ValuesContainerType& new_knots1,
                                                           ValuesContainerType& new_knots2,
                                                           ValuesContainerType& new_knots3,
                                                           const int& p1,
                                                           const int& p2,
                                                           const int& p3,
                                                           const ValuesContainerType2& knots1,
                                                           const ValuesContainerType2& knots2,
                                                           const ValuesContainerType2& knots3,
                                                           const ValuesContainerType3& ins_knots1,
                                                           const ValuesContainerType3& ins_knots2,
                                                           const ValuesContainerType3& ins_knots3)
    {
        MatrixType D1, D2, D3;
        ComputeBsplinesKnotInsertionCoefficients1D(D1, new_knots1, p1, knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1D(D2, new_knots2, p2, knots2, ins_knots2);
        ComputeBsplinesKnotInsertionCoefficients1D(D3, new_knots3, p3, knots3, ins_knots3);

        MatrixType tmp;
        IsogeometricMathUtils::outer_prod_mat(tmp, D2, D1);
        IsogeometricMathUtils::outer_prod_mat(D, D3, tmp);
    }

    /// Compute the refinement coefficients for one knot insertion NURBS refinement in 1D
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3, class ValuesContainerType4>
    static void ComputeNURBSKnotInsertionCoefficients1D(MatrixType& D,
                                                        ValuesContainerType& new_knots,
                                                        ValuesContainerType& new_weights,
                                                        const int& p,
                                                        const ValuesContainerType2& knots,
                                                        const ValuesContainerType3& ins_knots,
                                                        const ValuesContainerType4& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients1D(D, new_knots, p, knots, ins_knots);

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
    }

    /// Compute the refinement coefficients for multiple knots insertion NURBS refinement in 2D
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3, class ValuesContainerType4>
    static void ComputeNURBSKnotInsertionCoefficients2D(MatrixType& D,
                                                        ValuesContainerType& new_knots1,
                                                        ValuesContainerType& new_knots2,
                                                        ValuesContainerType& new_weights,
                                                        const int& p1,
                                                        const int& p2,
                                                        const ValuesContainerType2& knots1,
                                                        const ValuesContainerType2& knots2,
                                                        const ValuesContainerType3& ins_knots1,
                                                        const ValuesContainerType3& ins_knots2,
                                                        const ValuesContainerType4& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients2D(D, new_knots1, new_knots2, p1, p2, knots1, knots2, ins_knots1, ins_knots2);
//        std::cout << "B-splines D:" << D << std::endl;

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }
//        std::cout << "new_weights:";
//        for(int i = 0; i < new_weights.size(); ++i)
//            std::cout << " " << new_weights[i];
//        std::cout << std::endl;

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
//        KRATOS_WATCH(D)
    }

    /// Compute the refinement coefficients for multiple knots insertion NURBS refinement in 3D
    template<class MatrixType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3, class ValuesContainerType4>
    static void ComputeNURBSKnotInsertionCoefficients3D(MatrixType& D,
                                                        ValuesContainerType& new_knots1,
                                                        ValuesContainerType& new_knots2,
                                                        ValuesContainerType& new_knots3,
                                                        ValuesContainerType& new_weights,
                                                        const int& p1,
                                                        const int& p2,
                                                        const int& p3,
                                                        const ValuesContainerType2& knots1,
                                                        const ValuesContainerType2& knots2,
                                                        const ValuesContainerType2& knots3,
                                                        const ValuesContainerType3& ins_knots1,
                                                        const ValuesContainerType3& ins_knots2,
                                                        const ValuesContainerType3& ins_knots3,
                                                        const ValuesContainerType4& weights)
    {
        // firstly compute refinement matrix for B-splines refinement
        ComputeBsplinesKnotInsertionCoefficients3D(D, new_knots1, new_knots2, new_knots3, p1, p2, p3, knots1, knots2, knots3, ins_knots1, ins_knots2, ins_knots3);

        // secondly compute new weights
        new_weights.resize(D.size2());
        for(int j = 0; j < D.size2(); ++j)
        {
            new_weights[j] = 0.0;
            for(int i = 0; i < D.size1(); ++i)
                new_weights[j] += D(i, j) * weights[i];
        }

        // thirdly update the coefficient matrix
        for(int i = 0; i < D.size1(); ++i)
            for(int j = 0; j < D.size2(); ++j)
                D(i, j) *= weights[i]/new_weights[j];
    }

    /**
     * Compute the refinement coefficients for multiple knot insertion local refinement in 1D
     * REF: M. Scott T-splines paper
     */
    template<class VectorType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients1DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots,
                                                                const int& p,
                                                                const ValuesContainerType2& local_knots,
                                                                const ValuesContainerType3& ins_knots)
    {
        // compute the extended knot vector
        ValuesContainerType Ubar;
        int nt;
        IsogeometricMathUtils::compute_extended_knot_vector(Ubar, nt, local_knots, p);
//        KRATOS_WATCH(nt)

        // compute the refinement matrix for extended knot vector
        Matrix M;
        ComputeBsplinesKnotInsertionCoefficients1D(M, new_knots, p, Ubar, ins_knots);
//        KRATOS_WATCH(M)

        // extract the correct row from the refinement matrix
        int num_basis = ins_knots.size() + 1;
        if(D.size() != num_basis)
            D.resize(num_basis);
        for(unsigned int i = 0; i < num_basis; ++i)
            D[i] = M(nt, nt + i);

        #ifdef FUNCTIONALITY_CHECK
        for(unsigned int i = 0; i < nt; ++i)
            if(M(nt, i) != 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "M(nt,..) is nonzero at", i)
        for(unsigned int i = 0; i < num_basis; ++i)
            if(D(i) == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "D(..) is zero at", i)
        for(unsigned int i = nt + num_basis; i < M.size2(); ++i)
            if(M(nt, i) != 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "M(nt,..) is nonzero at", i)
        #endif
    }

    /// Compute the refinement coefficients for multiple knot insertion local refinement in 2D
    template<class VectorType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients2DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots1,
                                                                ValuesContainerType& new_knots2,
                                                                const int& p1,
                                                                const int& p2,
                                                                const ValuesContainerType2& local_knots1,
                                                                const ValuesContainerType2& local_knots2,
                                                                const ValuesContainerType3& ins_knots1,
                                                                const ValuesContainerType3& ins_knots2)
    {
        VectorType D1, D2;
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D1, new_knots1, p1, local_knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D2, new_knots2, p2, local_knots2, ins_knots2);
        IsogeometricMathUtils::outer_prod_vec(D, D2, D1);
    }

    /// Compute the refinement coefficients for multiple knot insertion local refinement in 3D
    template<class VectorType, class ValuesContainerType, class ValuesContainerType2, class ValuesContainerType3>
    static void ComputeBsplinesKnotInsertionCoefficients3DLocal(VectorType& D,
                                                                ValuesContainerType& new_knots1,
                                                                ValuesContainerType& new_knots2,
                                                                ValuesContainerType& new_knots3,
                                                                const int& p1,
                                                                const int& p2,
                                                                const int& p3,
                                                                const ValuesContainerType2& local_knots1,
                                                                const ValuesContainerType2& local_knots2,
                                                                const ValuesContainerType2& local_knots3,
                                                                const ValuesContainerType3& ins_knots1,
                                                                const ValuesContainerType3& ins_knots2,
                                                                const ValuesContainerType3& ins_knots3)
    {
        VectorType D1, D2, D3;
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D1, new_knots1, p1, local_knots1, ins_knots1);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D2, new_knots2, p2, local_knots2, ins_knots2);
        ComputeBsplinesKnotInsertionCoefficients1DLocal(D3, new_knots3, p3, local_knots3, ins_knots3);

        VectorType aux;
        IsogeometricMathUtils::outer_prod_vec(aux, D2, D1);
        IsogeometricMathUtils::outer_prod_vec(D, D3, aux);
    }

    /* "Adapted" modified version of Algorithm A5.9 from 'The NURBS BOOK' pg206. */
    /// @param d    degree of the B-Splines
    /// @param c    value container of the control points. Control points are the matrix [mc x nc]
    /// @param mc   dimension of the control point
    /// @param nc   number of control points
    /// @param k    value container of the knot vector
    /// @param nk   number of knots
    /// @param t    the order increment
    /// @param nh   number of resulted control points (maybe)
    /// @param ic   new control points
    /// @param ik   new knots
    /// REF: int bspdegelev(int d, double *c, int mc, int nc, double *k, int nk,
    ///                          int t, int *nh, double *ic, double *ik)
    /// REMARKS: This function can also be used to elevate the degree of NURBS
    template<typename TDataType, class ValuesContainerType, class ValuesContainerType1, class ValuesContainerType2>
    static int ComputeBsplinesDegreeElevation1D(const int& d, // order of B-Splines
            const ValuesContainerType& ctrl, // control values
            const ValuesContainerType1& k, // knot vector
            const int& t, // order increment
            ValuesContainerType& ictrl, // new control values
            ValuesContainerType2& ik, // new knot vector
            const TDataType& zero) // a sample zero control value to avoid explicit declaring TDataType
    {
      int nc = ctrl.size();
      int nk = k.size();
      int row, col;

      int ierr = 0;
      int i, j, q, s, m, ph, ph2, mpi, mh, nh, r, a, b, cind, oldr, mul;
      int n, lbz, rbz, save, tr, kj, first, kind, last, bet, ii;
      int nic, nik;
      double inv, ua, ub, numer, den, alf, gam;
      double **bezalfs, *alfs;

      /* allocate work space t times larger than original number */
      /* of control points and knots */
      ik.resize(nk*(t+1));
      ictrl.resize(nc*(t+1));

      n = nc - 1;

      // bezalfs = matrix(d+1,d+t+1);
      // bezalfs = (double**) calloc(d+1, sizeof(double*));
      // for (int i = 0; i < d+1; ++i)
      //   bezalfs[i] = (double*) calloc(d+t+1, sizeof(double));
      bezalfs = (double**) calloc(d+t+1, sizeof(double*));
      for (int i = 0; i < d+t+1; ++i)
        bezalfs[i] = (double*) calloc(d+1, sizeof(double));

      ValuesContainerType bpts(d+1);
      ValuesContainerType ebpts(d+t+1);
      ValuesContainerType Nextbpts(d+1);
      alfs = (double *) calloc(d, sizeof(double));

      m = n + d + 1;
      ph = d + t;
      ph2 = ph / 2;

      /* compute bezier degree elevation coefficients   */
      bezalfs[0][0] = bezalfs[ph][d] = 1.0;

      for (i = 1; i <= ph2; i++)
      {
        inv = 1.0 / bincoeff(ph,i);
        mpi = std::min(d,i);

        for (j = std::max(0,i-t); j <= mpi; j++)
          bezalfs[i][j] = inv * bincoeff(d,j) * bincoeff(t,i-j);
      }

      for (i = ph2+1; i <= ph-1; i++)
      {
        mpi = std::min(d, i);
        for (j = std::max(0,i-t); j <= mpi; j++)
          bezalfs[i][j] = bezalfs[ph-i][d-j];
      }

      mh = ph;
      kind = ph+1;
      r = -1;
      a = d;
      b = d+1;
      cind = 1;
      ua = k[0];

      ictrl[0] = ctrl[0];

      for (i = 0; i <= ph; i++)
        ik[i] = ua;

      /* initialise first bezier seg */
      for (i = 0; i <= d; i++)
          bpts[i] = ctrl[i];

      /* big loop thru knot vector */
      while (b < m)
      {
        i = b;
        while (b < m && k[b] == k[b+1])
          b++;

        mul = b - i + 1;
        mh += mul + t;
        ub = k[b];
        oldr = r;
        r = d - mul;

        /* insert knot u(b) r times */
        if (oldr > 0)
          lbz = (oldr+2) / 2;
        else
          lbz = 1;

        if (r > 0)
          rbz = ph - (r+1)/2;
        else
          rbz = ph;

        if (r > 0)
        {
          /* insert knot to get bezier segment */
          numer = ub - ua;
          for (q = d; q > mul; q--)
            alfs[q-mul-1] = numer / (k[a+q]-ua);
          for (j = 1; j <= r; j++)
          {
            save = r - j;
            s = mul + j;

            for (q = d; q >= s; q--)
                bpts[q] = alfs[q-s]*bpts[q]+(1.0-alfs[q-s])*bpts[q-1];

            Nextbpts[save] = bpts[d];
          }
        }
        /* end of insert knot */

        /* degree elevate bezier */
        for (i = lbz; i <= ph; i++)
        {
          ebpts[i] = zero;
          mpi = std::min(d, i);
          for (j = std::max(0,i-t); j <= mpi; j++)
            ebpts[i] = ebpts[i] + bezalfs[i][j]*bpts[j];
        }
        /* end of degree elevating bezier */

        if (oldr > 1)
        {
          /* must remove knot u=k[a] oldr times */
          first = kind - 2;
          last = kind;
          den = ub - ua;
          bet = (ub-ik[kind-1]) / den;

          /* knot removal loop */
          for (tr = 1; tr < oldr; tr++)
          {
            i = first;
            j = last;
            kj = j - kind + 1;
            while (j - i > tr)
            {
              /* loop and compute the new control points */
             /* for one removal step    */
              if (i < cind)
              {
                alf = (ub-ik[i])/(ua-ik[i]);
                ictrl[i] = alf * ictrl[i] + (1.0-alf) * ictrl[i-1];
              }
              if (j >= lbz)
              {
                if (j-tr <= kind-ph+oldr)
                {
                  gam = (ub-ik[j-tr]) / den;
                  ebpts[kj] = gam*ebpts[kj] + (1.0-gam)*ebpts[kj+1];
                }
                else
                {
                  ebpts[kj] = bet*ebpts[kj] + (1.0-bet)*ebpts[kj+1];
                }
              }
              i++;
              j--;
              kj--;
            }

            first--;
            last++;
          }
        }
        /* end of removing knot n=k[a] */

        /* load the knot ua  */
        if (a != d)
          for (i = 0; i < ph-oldr; i++)
          {
            ik[kind] = ua;
            kind++;
          }

        /* load ctrl pts into ic  */
        for (j = lbz; j <= rbz; j++)
        {
          ictrl[cind] = ebpts[j];
          cind++;
        }

        if (b < m)
        {
          /* setup for next pass thru loop  */
          for (j = 0; j < r; j++)
            bpts[j] = Nextbpts[j];
          for (j = r; j <= d; j++)
            bpts[j] = ctrl[b-d+j];
          a = b;
          b++;
          ua = ub;
        }
        else
          /* end knot  */
          for (i = 0; i <= ph; i++)
            ik[kind+i] = ub;
      }
      /* end while loop  */

      nh = mh - ph - 1;
      nic = nh + 1;
      nik = nic + d + t + 1;

      // resize to the new size
      ictrl.resize(nic);
      ik.resize(nik);

      free(alfs);
      // for (int i = 0; i < d+1; ++i)
      for (int i = 0; i < d+t+1; ++i)
        free(bezalfs[i]);
      free(bezalfs);

      return(ierr);
    }

    /// Degree elevation for B-Splines surface
    /// REMARKS: This function can also be used to elevate the degree of NURBS
    template<typename TDataType, class ValuesContainerType, class ValuesContainerType1, class ValuesContainerType2>
    static int ComputeBsplinesDegreeElevation2D(const int& d1, const int& d2, // order of B-Splines
            const ValuesContainerType& ctrl, // control values
            const ValuesContainerType1& k1, // knot vector
            const ValuesContainerType1& k2, // knot vector
            const int& t1, // order increment
            const int& t2, // order increment
            ValuesContainerType& ictrl, // new control values
            ValuesContainerType2& ik1, // new knot vector
            ValuesContainerType2& ik2, // new knot vector
            const TDataType& zero) // a sample zero control value to avoid explicit declaring TDataType
    {
        // firstly elevate the degree along the v-direction
        // this can be done by collecting the v-line of control values and then elevate the degree

        int n1 = k1.size() - d1 - 1;
        int n2 = k2.size() - d2 - 1;

        int new_n1, new_n2;

        ValuesContainerType temp_ictrl(n1, n2); // n1, n2 is just a temporary size

        for (int i = 0; i < n1; ++i)
        {
            std::vector<typename ValuesContainerType::DataType> ctrl_v(n2);
            std::vector<typename ValuesContainerType::DataType> ictrl_v(n2); // n2 is just temporary size

            for (int j = 0; j < n2; ++j)
                ctrl_v[j] = ctrl(i, j);

            ComputeBsplinesDegreeElevation1D(d2, ctrl_v, k2, t2, ictrl_v, ik2, zero);
            // here it's assumed that newly ik2 is the same for each iteration

            if (i == 0)
            {
                new_n2 = ik2.size() - (d2 + t2) - 1;
                temp_ictrl.resize(n1, new_n2);
            }

            for (int j = 0; j < new_n2; ++j)
                temp_ictrl(i, j) = ictrl_v[j];
        }

        // the degree in the second dimension can be elevated in the same way
        for (int j = 0; j < new_n2; ++j)
        {
            std::vector<typename ValuesContainerType::DataType> ctrl_u(n1);
            std::vector<typename ValuesContainerType::DataType> ictrl_u(n1); // n1 is just temporary size

            for (int i = 0; i < n1; ++i)
                ctrl_u[i] = temp_ictrl(i, j);

            ComputeBsplinesDegreeElevation1D(d1, ctrl_u, k1, t1, ictrl_u, ik1, zero);
            // here it's assumed that newly ik1 is the same for each iteration

            if (j == 0)
            {
                new_n1 = ik1.size() - (d1 + t1) - 1;
                ictrl.resize(new_n1, new_n2);
            }

            for (int i = 0; i < new_n1; ++i)
                ictrl(i, j) = ictrl_u[i];
        }

        return 0;
    }

    /// Degree elevation for B-Splines volume
    /// REMARKS: This function can also be used to elevate the degree of NURBS
    template<typename TDataType, class ValuesContainerType, class ValuesContainerType1, class ValuesContainerType2>
    static int ComputeBsplinesDegreeElevation3D(const int& d1, const int& d2, const int& d3, // order of B-Splines
            const ValuesContainerType& ctrl, // control values
            const ValuesContainerType1& k1, // knot vector
            const ValuesContainerType1& k2, // knot vector
            const ValuesContainerType1& k3, // knot vector
            const int& t1, // order increment
            const int& t2, // order increment
            const int& t3, // order increment
            ValuesContainerType& ictrl, // new control values
            ValuesContainerType2& ik1, // new knot vector
            ValuesContainerType2& ik2, // new knot vector
            ValuesContainerType2& ik3, // new knot vector
            const TDataType& zero) // a sample zero control value to avoid explicit declaring TDataType
    {
        // firstly elevate the degree along the w-direction
        // this can be done by collecting the w-line of control values and then elevate the degree

        int n1 = k1.size() - d1 - 1;
        int n2 = k2.size() - d2 - 1;
        int n3 = k3.size() - d3 - 1;

        int new_n1, new_n2, new_n3;

        ValuesContainerType temp_ictrl(n1, n2, n3); // n1, n2, n3 is just a temporary size

        for (int i = 0; i < n1; ++i)
        {
            for (int j = 0; j < n2; ++j)
            {
                std::vector<typename ValuesContainerType::DataType> ctrl_w(n3);
                std::vector<typename ValuesContainerType::DataType> ictrl_w(n3); // n3 is just temporary size

                for (int k = 0; k < n3; ++k)
                    ctrl_w[k] = ctrl(i, j, k);

                ComputeBsplinesDegreeElevation1D(d3, ctrl_w, k3, t3, ictrl_w, ik3, zero);
                // here it's assumed that newly ik3 is the same for each iteration

                if (i == 0 && j == 0)
                {
                    new_n3 = ik3.size() - (d3 + t3) - 1;
                    temp_ictrl.resize(n1, n2, new_n3);
                }

                for (int k = 0; k < new_n3; ++k)
                    temp_ictrl(i, j, k) = ictrl_w[k];
            }
        }

        ValuesContainerType temp_ictrl2(n1, n2, new_n3); // n1, n2, new_n3 is just a temporary size

        // the degree in the second dimension can be elevated in the same way
        for (int i = 0; i < n1; ++i)
        {
            for (int k = 0; k < new_n3; ++k)
            {
                std::vector<typename ValuesContainerType::DataType> ctrl_v(n2);
                std::vector<typename ValuesContainerType::DataType> ictrl_v(n2); // n2 is just temporary size

                for (int j = 0; j < n2; ++j)
                    ctrl_v[j] = temp_ictrl(i, j, k);

                ComputeBsplinesDegreeElevation1D(d2, ctrl_v, k2, t2, ictrl_v, ik2, zero);
                // here it's assumed that newly ik2 is the same for each iteration

                if (i == 0 && k == 0)
                {
                    new_n2 = ik2.size() - (d2 + t2) - 1;
                    temp_ictrl2.resize(n1, new_n2, new_n3);
                }

                for (int j = 0; j < new_n2; ++j)
                    temp_ictrl2(i, j, k) = ictrl_v[j];
            }
        }

        // lastly elevate the degree in the first dimension
        for (int j = 0; j < new_n2; ++j)
        {
            for (int k = 0; k < new_n3; ++k)
            {
                std::vector<typename ValuesContainerType::DataType> ctrl_u(n1);
                std::vector<typename ValuesContainerType::DataType> ictrl_u(n1); // n1 is just temporary size

                for (int i = 0; i < n1; ++i)
                    ctrl_u[i] = temp_ictrl2(i, j, k);

                ComputeBsplinesDegreeElevation1D(d1, ctrl_u, k1, t1, ictrl_u, ik1, zero);
                // here it's assumed that newly ik1 is the same for each iteration

                if (j == 0 && k == 0)
                {
                    new_n1 = ik1.size() - (d1 + t1) - 1;
                    ictrl.resize(new_n1, new_n2, new_n3);
                }

                for (int i = 0; i < new_n1; ++i)
                    ictrl(i, j, k) = ictrl_u[i];
            }
        }

        return 0;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Testing
    ///@{

    void test_ComputeBsplinesKnotInsertionCoefficients1DLocal();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BSplineUtils";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215. */
    static double bincoeff(const int& n, const int& k)
    {
        return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
    }

    /* computes ln(n!) */
    /* Numerical Recipes in C */
    /* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215. */
    static double factln(const int& n)
    {
        static int ntop = 0;
        static double a[101];

        if (n <= 1) return 0.0;
        while (n > ntop)
        {
            ++ntop;
            a[ntop] = gammaln(ntop+1.0);
        }
        return a[n];
    }

    /* Compute logarithm of the gamma function */
    /* Algorithm from 'Numerical Recipes in C, 2nd Edition' pg214. */
    static double gammaln(const double& xx)
    {
        double x,y,tmp,ser;
        static double cof[6] = {76.18009172947146,-86.50532032291677,
                                24.01409824083091,-1.231739572450155,
                                0.12086650973866179e-2, -0.5395239384953e-5};
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x+0.5) * log(tmp);
        ser = 1.000000000190015;
        for (j=0; j<=5; j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BSplineUtils& operator=(BSplineUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BSplineUtils(BSplineUtils const& rOther)
    {
    }

    ///@}

}; // Class BSplineUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BSplineUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const BSplineUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef FUNCTIONALITY_CHECK

#endif // KRATOS_BSPLINE_UTILS_H_INCLUDED  defined
