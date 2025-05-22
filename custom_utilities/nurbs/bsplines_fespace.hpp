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

template<int TDim>
template<typename TLocalCoordinateType>
void BSplinesFESpace_Helper<TDim>::GetValues(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
        std::vector<double>& values, const std::vector<TLocalCoordinateType>& xi)
{
    if constexpr (TDim == 1)
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
    else if constexpr (TDim == 2)
    {
        // inititialize the shape functions
        if (values.size() != rFESpace.TotalNumber())
        {
            values.resize(rFESpace.TotalNumber());
        }
        std::fill(values.begin(), values.end(), 0.0);

        // locate the knot span
        std::vector<int> Span(2);
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

        std::vector<int> Start(2);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                const unsigned int Index = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, rFESpace.Number(0), rFESpace.Number(1));

                const double N1 = ShapeFunctionValues1[i - Start[0]];
                const double N2 = ShapeFunctionValues2[j - Start[1]];

                values[Index] = N1 * N2;
            }
        }
    }
    else if constexpr (TDim == 3)
    {
        // initialize the shape functions
        if (values.size() != rFESpace.TotalNumber())
        {
            values.resize(rFESpace.TotalNumber());
        }
        std::fill(values.begin(), values.end(), 0.0);

        // locate the knot span
        std::vector<int> Span(3);
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

        std::vector<int> Start(3);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);
        Start[2] = Span[2] - rFESpace.Order(2);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                for (unsigned int k = Start[2]; k <= Span[2]; ++k)
                {
                    const unsigned int Index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, rFESpace.Number(0), rFESpace.Number(1), rFESpace.Number(2));

                    const double N1 = ShapeFunctionValues1[i - Start[0]];
                    const double N2 = ShapeFunctionValues2[j - Start[1]];
                    const double N3 = ShapeFunctionValues3[k - Start[2]];

                    values[Index] = N1 * N2 * N3;
                }
            }
        }
    }
}

template<int TDim>
template<typename TLocalCoordinateType>
void BSplinesFESpace_Helper<TDim>::GetValuesAndDerivatives(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
        std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<TLocalCoordinateType>& xi)
{
    if constexpr (TDim == 1)
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

        // distribute it into the array

        int Start;
        Start = Span - rFESpace.Order(0);

        for (unsigned int i = Start; i <= Span; ++i)
        {
            const unsigned int Index = BSplinesIndexingUtility_Helper::Index1D(i + 1, rFESpace.Number(0));

            const double N = ShapeFunctionsValuesAndDerivatives[0][i - Start];
            const double dN = ShapeFunctionsValuesAndDerivatives[1][i - Start];

            values[Index] = N;
            derivatives[Index][0] = dN;
        }
    }
    else if constexpr (TDim == 2)
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
        std::vector<int> Span(2);
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

        std::vector<int> Start(2);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                const unsigned int Index = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, rFESpace.Number(0), rFESpace.Number(1));

                const double N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
                const double dN1 = ShapeFunctionsValuesAndDerivatives1[1][i - Start[0]];
                const double N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];
                const double dN2 = ShapeFunctionsValuesAndDerivatives2[1][j - Start[1]];

                values[Index] = N1 * N2;
                derivatives[Index][0] = dN1 * N2;
                derivatives[Index][1] = N1 * dN2;
            }
        }
    }
    else if constexpr (TDim == 3)
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
        std::vector<int> Span(3);
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

        std::vector<int> Start(3);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);
        Start[2] = Span[2] - rFESpace.Order(2);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                for (unsigned int k = Start[2]; k <= Span[2]; ++k)
                {
                    const unsigned int Index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, rFESpace.Number(0), rFESpace.Number(1), rFESpace.Number(2));

                    const double N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
                    const double dN1 = ShapeFunctionsValuesAndDerivatives1[1][i - Start[0]];
                    const double N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];
                    const double dN2 = ShapeFunctionsValuesAndDerivatives2[1][j - Start[1]];
                    const double N3 = ShapeFunctionsValuesAndDerivatives3[0][k - Start[2]];
                    const double dN3 = ShapeFunctionsValuesAndDerivatives3[1][k - Start[2]];

                    values[Index] = N1 * N2 * N3;
                    derivatives[Index][0] = dN1 * N2 * N3;
                    derivatives[Index][1] = N1 * dN2 * N3;
                    derivatives[Index][2] = N1 * N2 * dN3;
                }
            }
        }
    }
}

template<int TDim>
template<typename TLocalCoordinateType>
void BSplinesFESpace_Helper<TDim>::GetValuesAndDerivatives(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
        const unsigned int nd,
        std::vector<double>& values,
        std::vector<std::vector<std::vector<double> > >& derivatives,
        const std::vector<TLocalCoordinateType>& xi)
{
    if constexpr (TDim == 1)
    {
        // initialize the shape functions and derivatives
        if (values.size() != rFESpace.TotalNumber())
        {
            values.resize(rFESpace.TotalNumber());
        }
        std::fill(values.begin(), values.end(), 0.0);

        if (derivatives.size() != nd)
            derivatives.resize(nd);

        for (unsigned int i = 0; i < nd; ++i)
        {
            if (derivatives[i].size() != rFESpace.TotalNumber())
                derivatives[i].resize(rFESpace.TotalNumber());

            for (unsigned int j = 0; j < derivatives[i].size(); ++j)
            {
                if (derivatives[i][j].size() != 1)
                {
                    derivatives[i][j].resize(1);
                    derivatives[i][j][0] = 0.0;
                }
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
        const int NumberOfDerivatives = static_cast<int>(nd);
        std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives;

        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives, Span, xi[0], rFESpace.Order(0), rFESpace.KnotVector(0), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());

        // distribute it into the array

        int Start;
        Start = Span - rFESpace.Order(0);

        for (unsigned int i = Start; i <= Span; ++i)
        {
            const unsigned int Index = BSplinesIndexingUtility_Helper::Index1D(i + 1, rFESpace.Number(0));

            const double N = ShapeFunctionsValuesAndDerivatives[0][i - Start];
            values[Index] = N;

            for (unsigned int j = 0; j < nd; ++j)
            {
                const double d = ShapeFunctionsValuesAndDerivatives[j+1][i - Start];
                derivatives[j][Index][0] = d;
            }
        }
    }
    else if constexpr (TDim == 2)
    {
        // initialize the shape functions and derivatives
        if (values.size() != rFESpace.TotalNumber())
        {
            values.resize(rFESpace.TotalNumber());
        }
        std::fill(values.begin(), values.end(), 0.0);

        if (derivatives.size() != nd)
            derivatives.resize(nd);

        for (unsigned int i = 0; i < nd; ++i)
        {
            if (derivatives[i].size() != rFESpace.TotalNumber())
                derivatives[i].resize(rFESpace.TotalNumber());

            unsigned int nders;
            if (i == 0) nders = 2;
            else if (i == 1) nders = 3;
            else KRATOS_ERROR << "Third order derivatives is to be implemented";

            for (unsigned int j = 0; j < derivatives[i].size(); ++j)
            {
                if (derivatives[i][j].size() != nders)
                {
                    derivatives[i][j].resize(nders);
                    for (unsigned int n = 0; n < nders; ++n)
                        derivatives[i][j][n] = 0.0;
                }
            }
        }

        // locate the knot span
        std::vector<int> Span(2);
        Span[0] = BSplineUtils::FindSpan(rFESpace.Number(0), rFESpace.Order(0), xi[0], rFESpace.KnotVector(0));
        Span[1] = BSplineUtils::FindSpan(rFESpace.Number(1), rFESpace.Order(1), xi[1], rFESpace.KnotVector(1));

        if ((Span[0] >= rFESpace.Number(0) + rFESpace.Order(0)) || (Span[0] == 0)
                || (Span[1] >= rFESpace.Number(1) + rFESpace.Order(1)) || (Span[1] == 0))
        {
            return;
        }

        // compute the non-zero shape function values and derivatives
        const int NumberOfDerivatives = static_cast<int>(nd);
        std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives1;
        std::vector<std::vector<double> > ShapeFunctionsValuesAndDerivatives2;

        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives1, Span[0], xi[0], rFESpace.Order(0), rFESpace.KnotVector(0), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());
        BSplineUtils::BasisFunsDer(ShapeFunctionsValuesAndDerivatives2, Span[1], xi[1], rFESpace.Order(1), rFESpace.KnotVector(1), NumberOfDerivatives, BSplineUtils::StdVector2DOp<double>());

        // distribute the values to arrays

        std::vector<int> Start(2);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                const unsigned int Index = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, rFESpace.Number(0), rFESpace.Number(1));

                const double N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
                const double N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];

                values[Index] = N1 * N2;

                // first derivatives
                unsigned int k = 0;
                {
                    const double d1 = ShapeFunctionsValuesAndDerivatives1[k+1][i - Start[0]];
                    const double d2 = ShapeFunctionsValuesAndDerivatives2[k+1][j - Start[1]];

                    derivatives[k][Index][0] = d1 * N2;
                    derivatives[k][Index][1] = N1 * d2;
                }

                // second derivatives
                k = 1;
                {
                    const double d2_1 = ShapeFunctionsValuesAndDerivatives1[k+1][i - Start[0]];
                    const double d2_2 = ShapeFunctionsValuesAndDerivatives2[k+1][j - Start[1]];

                    const double d1 = ShapeFunctionsValuesAndDerivatives1[k][i - Start[0]];
                    const double d2 = ShapeFunctionsValuesAndDerivatives2[k][j - Start[1]];

                    derivatives[k][Index][0] = d2_1 * N2;
                    derivatives[k][Index][1] = N1 * d2_2;
                    derivatives[k][Index][2] = d1 * d2;
                }

                // TODO higher derivatives
            }
        }
    }
    else if constexpr (TDim == 3)
    {
        // initialize the shape functions and derivatives
        if (values.size() != rFESpace.TotalNumber())
        {
            values.resize(rFESpace.TotalNumber());
        }
        std::fill(values.begin(), values.end(), 0.0);

        if (derivatives.size() != nd)
            derivatives.resize(nd);

        for (unsigned int i = 0; i < nd; ++i)
        {
            if (derivatives[i].size() != rFESpace.TotalNumber())
                derivatives[i].resize(rFESpace.TotalNumber());

            unsigned int nders;
            if (i == 0) nders = 3;
            else if (i == 1) nders = 6;
            else KRATOS_ERROR << "Higher order derivatives are to be implemented";

            for (unsigned int j = 0; j < derivatives[i].size(); ++j)
            {
                if (derivatives[i][j].size() != nders)
                {
                    derivatives[i][j].resize(nders);
                    for (unsigned int n = 0; n < nders; ++n)
                        derivatives[i][j][n] = 0.0;
                }
            }
        }

        // locate the knot span
        std::vector<int> Span(3);
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

        std::vector<int> Start(3);
        Start[0] = Span[0] - rFESpace.Order(0);
        Start[1] = Span[1] - rFESpace.Order(1);
        Start[2] = Span[2] - rFESpace.Order(2);

        for (unsigned int i = Start[0]; i <= Span[0]; ++i)
        {
            for (unsigned int j = Start[1]; j <= Span[1]; ++j)
            {
                for (unsigned int k = Start[2]; k <= Span[2]; ++k)
                {
                    const unsigned int Index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, rFESpace.Number(0), rFESpace.Number(1), rFESpace.Number(2));

                    const double N1 = ShapeFunctionsValuesAndDerivatives1[0][i - Start[0]];
                    const double N2 = ShapeFunctionsValuesAndDerivatives2[0][j - Start[1]];
                    const double N3 = ShapeFunctionsValuesAndDerivatives3[0][k - Start[2]];

                    values[Index] = N1 * N2 * N3;

                    // first derivatives
                    unsigned int l = 0;
                    {
                        const double d1 = ShapeFunctionsValuesAndDerivatives1[l+1][i - Start[0]];
                        const double d2 = ShapeFunctionsValuesAndDerivatives2[l+1][j - Start[1]];
                        const double d3 = ShapeFunctionsValuesAndDerivatives2[l+1][k - Start[2]];

                        derivatives[l][Index][0] = d1 * N2 * N3;
                        derivatives[l][Index][1] = N1 * d2 * N3;
                        derivatives[l][Index][2] = N1 * N2 * d3;
                    }

                    // second derivatives
                    l = 1;
                    {
                        const double d2_1 = ShapeFunctionsValuesAndDerivatives1[l+1][i - Start[0]];
                        const double d2_2 = ShapeFunctionsValuesAndDerivatives2[l+1][j - Start[1]];
                        const double d2_3 = ShapeFunctionsValuesAndDerivatives2[l+1][k - Start[2]];

                        const double d1 = ShapeFunctionsValuesAndDerivatives1[l][i - Start[0]];
                        const double d2 = ShapeFunctionsValuesAndDerivatives2[l][j - Start[1]];
                        const double d3 = ShapeFunctionsValuesAndDerivatives2[l][k - Start[2]];

                        derivatives[l][Index][0] = d2_1 * N2 * N3;
                        derivatives[l][Index][1] = N1 * d2_2 * N3;
                        derivatives[l][Index][2] = N1 * N2 * d2_3;

                        derivatives[l][Index][3] = d1 * d2 * N3;
                        derivatives[l][Index][4] = N1 * d2 * d3;
                        derivatives[l][Index][5] = d1 * N2 * d3;
                    }

                    // higher derivatives TODO
                }
            }
        }
    }
}

} // namespace Kratos.
