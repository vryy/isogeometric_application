//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <vector>

// External includes

// Project includes
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"

namespace Kratos
{

template<>
void BSplinesFESpace_Helper<1>::GetValues(const BSplinesFESpace<1>& rFESpace,
        std::vector<double>& values, const std::vector<double>& xi)
{
    // initialize the shape functions
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    // locate the knot span
    int Span;
    Span = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));

    if ((Span >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span == 0))
    {
        return;
    }

    // compute the non-zero shape function values
    std::vector<double> ShapeFunctionValues(rFESpace.Order(0) + 1);

    BSplineUtils::BasisFuns(ShapeFunctionValues, Span, xi[0], rFESpace.Order(0), rFESpace.KnotVector(0));

    int Start;
    Start = Span - rFESpace.Order(0);

    double N;

    unsigned int i, Index;
    for (i = Start; i <= Span; ++i)
    {
        Index = BSplinesIndexingUtility_Helper::Index1D(i + 1, rFESpace.Number(0));

        N = ShapeFunctionValues[i - Start];

        values[Index] = N;
    }
}

template<>
void BSplinesFESpace_Helper<1>::GetValuesAndDerivatives(const BSplinesFESpace<1>& rFESpace,
        std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<double>& xi)
{
    // initialize the shape functions and derivatives
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    if (derivatives.size() != rFESpace.TotalNumber())
    {
        derivatives.resize(rFESpace.TotalNumber());
    }
    for (unsigned int i = 0; i < derivatives.size(); ++i)
    {
        if (derivatives[i].size() != 1)
        {
            derivatives[i].resize(1);
            derivatives[i][0] = 0.0;
        }
    }

    // locate the knot span
    int Span;
    Span = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));

    if ((Span >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span == 0))
    {
        return;
    }

    // compute the non-zero shape function values and derivatives
    const int NumberOfDerivatives = 1;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives;

    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives, Span, xi[0], rFESpace.Order(0), rFESpace.KnotVector(0), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());

    // for (int i = 0; i < ShapeFunctionsValuesAndDerivatives.size(); ++i)
    // {
    //     std::cout << "ShapeFunctionsValuesAndDerivatives[" << i << "]:";
    //     for (int j = 0; j < ShapeFunctionsValuesAndDerivatives[i].size(); ++j)
    //         std::cout << " " << ShapeFunctionsValuesAndDerivatives[i][j];
    //     std::cout << std::endl;
    // }

    // distribute it into the array

    int Start;
    Start = Span - rFESpace.Order(0);

    double N, dN;

    unsigned int i, Index;
    for (i = Start; i <= Span; ++i)
    {
        Index = BSplinesIndexingUtility_Helper::Index1D(i + 1, rFESpace.Number(0));

        N = ShapeFunctionsValuesAndDerivatives[0][i - Start];
        dN = ShapeFunctionsValuesAndDerivatives[1][i - Start];

        values[Index] = N;
        derivatives[Index][0] = dN;
    }
}

template<>
void BSplinesFESpace_Helper<2>::GetValues(const BSplinesFESpace<2>& rFESpace,
        std::vector<double>& values, const std::vector<double>& xi)
{
    // inititialize the shape functions
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    // locate the knot span
    int Span[2];
    Span[0] = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(rFESpace.Number(1), rFESpace.Order(1), xi[1], rFESpace.KnotVector(1));

    if ((Span[0] >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span[0] == 0)
            || (Span[1] >= rFESpace.Number(1) + rFESpace.Order(1)) || (Span[1] == 0))
    {
        return;
    }

    // compute the non-zero shape function values
    std::vector<double> ShapeFunctionValues1(rFESpace.Order(0) + 1);
    std::vector<double> ShapeFunctionValues2(rFESpace.Order(1) + 1);

    BSplineUtils::BasisFuns(ShapeFunctionValues1, Span[0], xi[0], rFESpace.Order(0), rFESpace.KnotVector(0));
    BSplineUtils::BasisFuns(ShapeFunctionValues2, Span[1], xi[1], rFESpace.Order(1), rFESpace.KnotVector(1));

    int Start[2];
    Start[0] = Span[0] - rFESpace.Order(0);
    Start[1] = Span[1] - rFESpace.Order(1);

    double N1, N2;

    unsigned int i, j, Index;
    for (i = Start[0]; i <= Span[0]; ++i)
    {
        for (j = Start[1]; j <= Span[1]; ++j)
        {
            Index = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, rFESpace.Number(0), rFESpace.Number(1));

            N1 = ShapeFunctionValues1[i - Start[0]];
            N2 = ShapeFunctionValues2[j - Start[1]];

            values[Index] = N1 * N2;
        }
    }
}

template<>
void BSplinesFESpace_Helper<2>::GetValuesAndDerivatives(const BSplinesFESpace<2>& rFESpace,
        std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<double>& xi)
{
    // initialize the shape functions and derivatives
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    if (derivatives.size() != rFESpace.TotalNumber())
    {
        derivatives.resize(rFESpace.TotalNumber());
    }
    for (unsigned int i = 0; i < rFESpace.TotalNumber(); ++i)
    {
        if (derivatives[i].size() != 2)
        {
            derivatives[i].resize(2);
            derivatives[i][0] = 0.0;
            derivatives[i][1] = 0.0;
        }
    }

    // locate the knot span
    int Span[2];
    Span[0] = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(rFESpace.Number(1), rFESpace.Order(1), xi[1], rFESpace.KnotVector(1));

    if ((Span[0] >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span[0] == 0)
            || (Span[1] >= rFESpace.Number(1) + rFESpace.Order(1)) || (Span[1] == 0))
    {
        return;
    }

    // compute the non-zero shape function values and derivatives
    const int NumberOfDerivatives = 1;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives1;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives2;

    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives1, Span[0], xi[0], rFESpace.Order(0), rFESpace.KnotVector(0), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());
    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives2, Span[1], xi[1], rFESpace.Order(1), rFESpace.KnotVector(1), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());

    // distribute the values to arrays

    int Start[2];
    Start[0] = Span[0] - rFESpace.Order(0);
    Start[1] = Span[1] - rFESpace.Order(1);

    double N1, N2, dN1, dN2;

    unsigned int i, j, Index;
    for (i = Start[0]; i <= Span[0]; ++i)
    {
        for (j = Start[1]; j <= Span[1]; ++j)
        {
            Index = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, rFESpace.Number(0), rFESpace.Number(1));

            N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
            dN1 = ShapeFunctionsValuesAndDerivatives1[1][i - Start[0]];
            N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];
            dN2 = ShapeFunctionsValuesAndDerivatives2[1][j - Start[1]];

            values[Index] = N1 * N2;
            derivatives[Index][0] = dN1 * N2;
            derivatives[Index][1] = N1 * dN2;
        }
    }
}

template<>
void BSplinesFESpace_Helper<3>::GetValues(const BSplinesFESpace<3>& rFESpace,
        std::vector<double>& values, const std::vector<double>& xi)
{
    // initialize the shape functions
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    // locate the knot span
    int Span[3];
    Span[0] = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(rFESpace.Number(1), rFESpace.Order(1), xi[1], rFESpace.KnotVector(1));
    Span[2] = BSplineUtils::FindSpan(rFESpace.Number(2), rFESpace.Order(2), xi[2], rFESpace.KnotVector(2));

    if ((Span[0] >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span[0] == 0)
            || (Span[1] >= rFESpace.Number(1) + rFESpace.Order(1)) || (Span[1] == 0)
            || (Span[2] >= rFESpace.Number(2) + rFESpace.Order(2)) || (Span[2] == 0))
    {
        return;
    }

    // compute the non-zero shape function values
    std::vector<double> ShapeFunctionValues1(rFESpace.Order(0) + 1);
    std::vector<double> ShapeFunctionValues2(rFESpace.Order(1) + 1);
    std::vector<double> ShapeFunctionValues3(rFESpace.Order(2) + 1);

    BSplineUtils::BasisFuns(ShapeFunctionValues1, Span[0], xi[0], rFESpace.Order(0), rFESpace.KnotVector(0));
    BSplineUtils::BasisFuns(ShapeFunctionValues2, Span[1], xi[1], rFESpace.Order(1), rFESpace.KnotVector(1));
    BSplineUtils::BasisFuns(ShapeFunctionValues3, Span[2], xi[2], rFESpace.Order(2), rFESpace.KnotVector(2));

    int Start[3];
    Start[0] = Span[0] - rFESpace.Order(0);
    Start[1] = Span[1] - rFESpace.Order(1);
    Start[2] = Span[2] - rFESpace.Order(2);

    double N1, N2, N3;

    unsigned int i, j, k, Index;
    for (i = Start[0]; i <= Span[0]; ++i)
    {
        for (j = Start[1]; j <= Span[1]; ++j)
        {
            for (k = Start[2]; k <= Span[2]; ++k)
            {
                Index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, rFESpace.Number(0), rFESpace.Number(1), rFESpace.Number(2));

                N1 = ShapeFunctionValues1[i - Start[0]];
                N2 = ShapeFunctionValues2[j - Start[1]];
                N3 = ShapeFunctionValues3[k - Start[2]];

                values[Index] = N1 * N2 * N3;
            }
        }
    }
}

template<>
void BSplinesFESpace_Helper<3>::GetValuesAndDerivatives(const BSplinesFESpace<3>& rFESpace,
        std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<double>& xi)
{
    // initialize the shape functions and derivatives
    if (values.size() != rFESpace.TotalNumber())
    {
        values.resize(rFESpace.TotalNumber());
    }
    std::fill(values.begin(), values.end(), 0.0);

    if (derivatives.size() != rFESpace.TotalNumber())
    {
        derivatives.resize(rFESpace.TotalNumber());
    }
    for (unsigned int i = 0; i < derivatives.size(); ++i)
    {
        if (derivatives[i].size() != 3)
        {
            derivatives[i].resize(3);
            derivatives[i][0] = 0.0;
            derivatives[i][1] = 0.0;
            derivatives[i][2] = 0.0;
        }
    }

    // locate the knot span
    int Span[3];
    Span[0] = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));
    Span[1] = BSplineUtils::FindSpan(rFESpace.Number(1), rFESpace.Order(1), xi[1], rFESpace.KnotVector(1));
    Span[2] = BSplineUtils::FindSpan(rFESpace.Number(2), rFESpace.Order(2), xi[2], rFESpace.KnotVector(2));

    if ((Span[0] >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span[0] == 0)
            || (Span[1] >= rFESpace.Number(1) + rFESpace.Order(1)) || (Span[1] == 0)
            || (Span[2] >= rFESpace.Number(2) + rFESpace.Order(2)) || (Span[2] == 0))
    {
        return;
    }

    // compute the non-zero shape function values and derivatives
    const int NumberOfDerivatives = 1;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives1;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives2;
    std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives3;

    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives1, Span[0], xi[0], rFESpace.Order(0), rFESpace.KnotVector(0), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());
    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives2, Span[1], xi[1], rFESpace.Order(1), rFESpace.KnotVector(1), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());
    BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives3, Span[2], xi[2], rFESpace.Order(2), rFESpace.KnotVector(2), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());

    // distribute the values to arrays

    int Start[3];
    Start[0] = Span[0] - rFESpace.Order(0);
    Start[1] = Span[1] - rFESpace.Order(1);
    Start[2] = Span[2] - rFESpace.Order(2);

    double N1, N2, N3, dN1, dN2, dN3;

    unsigned int i, j, k, Index;
    for (i = Start[0]; i <= Span[0]; ++i)
    {
        for (j = Start[1]; j <= Span[1]; ++j)
        {
            for (k = Start[2]; k <= Span[2]; ++k)
            {
                Index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, rFESpace.Number(0), rFESpace.Number(1), rFESpace.Number(2));

                N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
                dN1 = ShapeFunctionsValuesAndDerivatives1[1][i - Start[0]];
                N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];
                dN2 = ShapeFunctionsValuesAndDerivatives2[1][j - Start[1]];
                N3 = ShapeFunctionsValuesAndDerivatives3[0][k - Start[2]];
                dN3 = ShapeFunctionsValuesAndDerivatives3[1][k - Start[2]];

                values[Index] = N1 * N2 * N3;
                derivatives[Index][0] = dN1 * N2 * N3;
                derivatives[Index][1] = N1 * dN2 * N3;
                derivatives[Index][2] = N1 * N2 * dN3;
            }
        }
    }
}

} // namespace Kratos.
