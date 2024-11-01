//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014-01-28 $
//   Revision:            $Revision: 1.1 $
//
//

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/bspline_utils.h"

namespace Kratos
{

const int BezierUtils::msBernsteinCoefs[] =
{
    2                 // p = 2
    , 3                 // p = 3
    , 4, 6              // p = 4
    , 5, 10
    , 6, 15, 20
    , 7, 21, 35
    , 8, 28, 56, 70
    , 9, 36, 84, 126
    , 10, 45, 120, 210, 252
    , 11, 55, 165, 330, 462
    , 12, 66, 220, 495, 792, 924
    , 13, 78, 186, 715, 1287, 1716
    , 14, 91, 364, 1001, 2002, 3003, 3432
    , 15, 105, 455, 1365, 3003, 5005, 6435
    , 16, 120, 560, 1820, 4368, 8008, 11440, 12870
};

BezierUtils::MapType BezierUtils::mIntegrationMethods;

void BezierUtils::bezier_extraction_tsplines_1d(
    std::vector<Vector>& Crows,     // bezier extraction operator, each row of the operator contain the Bezier decomposition coefficients on each knot span (OUTPUT). The operator is of size nb x (p+1)
    int& nb,                        // number of knot spans (number of rows of the extraction operator) of the filled extended knot vector (OUTPUT)
    Vector& Ubar,                   // filled extended knot vector (OUTPUT)
    const std::vector<double>& Xi,  // local knot vector (INPUT)
    const std::vector<double>& U,   // vector of inserted knots (INPUT)
    const std::vector<int>& spans,  // vector of knot span of the inserted knots (INPUT) (1-based index)
    const int p)                    // degree of the basis function (INPUT)
{
    // bound checking
    if (Xi.size() != p + 2)
        KRATOS_THROW_ERROR(std::logic_error, "local knot vector must be of length p + 2", "")

        if (U.size() != spans.size())
            KRATOS_THROW_ERROR(std::logic_error, "spans array must be the same length as the number of interior knots", "")

            // compute the extended knot vector and identify the relative position of the T-splines basis function
            int nt;
    IsogeometricMathUtils<double>::compute_extended_knot_vector(Ubar, nt, Xi, p);

    // algorithm: modify from the paper: Isogeometric finite element data structure based on Bezier extraction of T-splines, Scott et al
    int a = p + 1;
    int b = a + 1;
    nb = 1;

    if (!Crows.empty())
    {
        Crows.clear();
    }
    Crows.push_back(ZeroVector(p + 1));
    Crows[0][nt] = 1.0;

    int m = U.size();
    int mbar = p + 2 + nt + m;
    int ki = 1;
    int si = 1;
    int i, j, k, loc, r, save, s, mult, add, total_add = 0;
    double numer, alpha;//, alphas[p];
    std::vector<double> alphas(p + 1);
    while (b < mbar)
    {
        // count multiplicity of knots at location b
        add = 0;
        if ((si <= m) && (spans[si - 1] == ki))
        {
            mult = 0;
            add = 1;
            // add the new knot to the knot vector
            Ubar.resize(Ubar.size() + 1);
//            for(i = mbar+p-m+si-1; i >= b; --i)
            for (i = Ubar.size() - 1; i >= b; --i)
            {
                Ubar[i] = Ubar[i - 1];
            }
            Ubar[b - 1] = U[si - 1];
            ++si;
        }
        else
        {
            ++ki;
            i = b;
            while (b < mbar)
            {
                if (Ubar[b] == Ubar[b - 1])
                {
                    ++b;
                }
                else
                {
                    break;
                }
            }
            mult = b - i + 1;
        }
        total_add += add; // count the total number of additional knots

        if (mult <= p)
        {
            // initialize the next extraction operator row
            Crows.push_back(ZeroVector(p + 1));
            loc = nt + 1 - nb + total_add; // identify the next location to be 1.0
            if ((loc >= 1) && (loc <= p + 1))
            {
                Crows[nb][loc - 1] = 1.0;
            }

            // compute the interpolation coefficients
            numer = Ubar[b - 1] - Ubar[a - 1];
            for (j = p; j >= mult + 1; --j)
            {
                alphas[j - mult] = numer / (Ubar[a + j + add - 1] - Ubar[a - 1]);
            }
            r = p - mult;

            // update the matrix coefficients for r new knots
            for (j = 1; j <= r; ++j)
            {
                save = r - j + 1;
                s = mult + j;
                for (k = p + 1; k >= s + 1; --k)
                {
                    alpha = alphas[k - s];
                    Crows[nb - 1][k - 1] = alpha * Crows[nb - 1][k - 1] + (1.0 - alpha) * Crows[nb - 1][k - 2];
                }
                if (b <= mbar)
                    // update overlapping coefficients of the next operator row
                {
                    Crows[nb][save - 1] = Crows[nb - 1][p];
                }
            }
            ++nb;
            if (b <= mbar)
            {
                // update indices for the next operator
                a = b;
                ++b;
            }
        }
//        else if(mult == p)
//        {
//            if(b <= mbar)
//            {
//                ++nb;
//                a = b;
//                ++b;
//            }
//        }
    }
}

void BezierUtils::bezier_extraction_tsplines_2d(std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const std::vector<int>& spans_xi,
        const std::vector<int>& spans_eta,
        const int p,
        const int q)
{
    std::vector<Vector> Cxi;
    std::vector<Vector> Ceta;
    bezier_extraction_tsplines_1d(Cxi, nb_xi, Ubar_xi, Xi, Uxi, spans_xi, p);
    bezier_extraction_tsplines_1d(Ceta, nb_eta, Ubar_eta, Eta, Ueta, spans_eta, q);
//    KRATOS_WATCH(nb_xi)
//    KRATOS_WATCH(nb_eta)
//    for(std::size_t i = 0; i < Cxi.size(); ++i)
//        std::cout << "Cxi[" << i << "]:" << Cxi[i] << std::endl;
//    for(std::size_t i = 0; i < Ceta.size(); ++i)
//        std::cout << "Ceta[" << i << "]:" << Ceta[i] << std::endl;

    if (Crows.size() != nb_xi * nb_eta)
    {
        Crows.resize(nb_xi * nb_eta);
    }
    for (std::size_t i = 0; i < Crows.size(); ++i)
        if (Crows[i].size() != (p + 1) * (q + 1))
        {
            Crows[i].resize((p + 1) * (q + 1));
        }

    int row_start, col_start;
    double aux;
    for (std::size_t i = 0; i < nb_xi; ++i)
        for (std::size_t j = 0; j < p + 1; ++j)
        {
            aux = Cxi[i][j];
            for (std::size_t k = 0; k < nb_eta; ++k)
                for (std::size_t l = 0; l < q + 1; ++l)
                {
                    row_start = i * nb_eta;
                    col_start = j * (q + 1);
                    Crows[row_start + k][col_start + l] = aux * Ceta[k][l];
                }
        }
}

void BezierUtils::bezier_extraction_local_1d(std::vector<Vector>& Crows,
        int& nb,
        Vector& Ubar,
        const std::vector<double>& Xi,
        const std::vector<double>& U,
        const int p)
{
    // bound checking
    if (Xi.size() != p + 2)
        KRATOS_THROW_ERROR(std::logic_error, "local knot vector must be of length p + 2", "")

        // compute the extended knot vector and identify the relative position of the T-splines basis function
        int nt;
    Vector Uextended;
    IsogeometricMathUtils<double>::compute_extended_knot_vector(Uextended, nt, Xi, p);

    // count the multiplicity of inner knots
    int i = p + 1;
    std::vector<double> Ud;
    std::vector<int> Um;
    int num_inner_knots = 0;
    while (i <= Uextended.size() - p - 2)
    {
        int mult = 0;
        while (Uextended(i + mult) == Uextended(i))
        {
            ++mult;
        }
        ++num_inner_knots;
        Ud.push_back(Uextended(i));
        Um.push_back(mult);
        i += mult;
    }

    // compute the inserted knots
    std::vector<double> ins_knots;
    for (int i = 0; i < num_inner_knots; ++i)
        if (Um[i] < p)
            for (int j = 0; j < p - Um[i]; ++j)
            {
                ins_knots.push_back(Ud[i]);
            }

    // append the external knots
    for (int i = 0; i < U.size(); ++i)
    {
        bool found = false;
        for (int j = 0; j < Ud.size(); ++j)
            if (Ud[j] == U[i])
            {
                found = true;
                break;
            }
        if (!found)
        {
            ++num_inner_knots;
            for (int j = 0; j < p; ++j)
            {
                ins_knots.push_back(U[i]);
            }
        }
    }

    // compute the full extraction operator
    Matrix D;
    BSplineUtils::ComputeBsplinesKnotInsertionCoefficients1D(D, Ubar, p, Uextended, ins_knots);

    // extract the local extraction operator
    Vector C(D.size2());
    noalias(C) = row(D, nt);
    nb = num_inner_knots + 1;
    int j = 0;
    Crows.clear();
    for (int i = 0; i < nb; ++i)
    {
        Vector Crow(p + 1);
        for (int k = 0; k < p + 1; ++k)
        {
            Crow(k) = C(j + k);
        }
        j += p;
        Crows.push_back(Crow);
    }
}

void BezierUtils::bezier_extraction_local_2d(std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const int p,
        const int q)
{
    std::vector<Vector> Cxi;
    std::vector<Vector> Ceta;

    bezier_extraction_local_1d(Cxi, nb_xi, Ubar_xi, Xi, Uxi, p);
    bezier_extraction_local_1d(Ceta, nb_eta, Ubar_eta, Eta, Ueta, q);

    if (Crows.size() != nb_xi * nb_eta)
    {
        Crows.resize(nb_xi * nb_eta);
    }
    for (std::size_t i = 0; i < Crows.size(); ++i)
        if (Crows[i].size() != (p + 1) * (q + 1))
        {
            Crows[i].resize((p + 1) * (q + 1));
        }

    unsigned int row, col;
    for (std::size_t i = 0; i < nb_xi; ++i)
        for (std::size_t j = 0; j < p + 1; ++j)
        {
            for (std::size_t k = 0; k < nb_eta; ++k)
                for (std::size_t l = 0; l < q + 1; ++l)
                {
                    row = i * nb_eta + k;
                    col = j * (q + 1) + l;
                    Crows[row][col] = Cxi[i][j] * Ceta[k][l];
                }
        }
}

void BezierUtils::bezier_extraction_local_3d(std::vector<Vector>& Crows,
        int& nb_xi,
        int& nb_eta,
        int& nb_zeta,
        Vector& Ubar_xi,
        Vector& Ubar_eta,
        Vector& Ubar_zeta,
        const std::vector<double>& Xi,
        const std::vector<double>& Eta,
        const std::vector<double>& Zeta,
        const std::vector<double>& Uxi,
        const std::vector<double>& Ueta,
        const std::vector<double>& Uzeta,
        const int p,
        const int q,
        const int r)
{
    std::vector<Vector> Cxi;
    std::vector<Vector> Ceta;
    std::vector<Vector> Czeta;

    bezier_extraction_local_1d(Cxi, nb_xi, Ubar_xi, Xi, Uxi, p);
    bezier_extraction_local_1d(Ceta, nb_eta, Ubar_eta, Eta, Ueta, q);
    bezier_extraction_local_1d(Czeta, nb_zeta, Ubar_zeta, Zeta, Uzeta, r);

    if (Crows.size() != nb_xi * nb_eta * nb_zeta)
    {
        Crows.resize(nb_xi * nb_eta * nb_zeta);
    }
    for (std::size_t i = 0; i < Crows.size(); ++i)
        if (Crows[i].size() != (p + 1) * (q + 1) * (r + 1))
        {
            Crows[i].resize((p + 1) * (q + 1) * (r + 1));
        }

    unsigned int row, col;
    for (std::size_t i = 0; i < nb_xi; ++i)
    {
        for (std::size_t j = 0; j < p + 1; ++j)
        {
            for (std::size_t k = 0; k < nb_eta; ++k)
                for (std::size_t l = 0; l < q + 1; ++l)
                {
                    for (std::size_t m = 0; m < nb_zeta; ++m)
                        for (std::size_t n = 0; n < r + 1; ++n)
                        {
                            row = (i * nb_eta + k) * nb_zeta + m;
                            col = (j * (q + 1) + l) * (r + 1) + n;
                            Crows[row][col] = Cxi[i][j] * Ceta[k][l] * Czeta[m][n];
                        }
                }
        }
    }
}

}// namespace Kratos.
