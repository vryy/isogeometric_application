//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

namespace Kratos
{

template<class TFESpaceType, typename TDataType, typename TLocalCoordinatesType>
void GridFunction_Helper<TFESpaceType, TDataType, TLocalCoordinatesType>::GetValue(TDataType& v, const TLocalCoordinatesType& xi,
                        const TFESpaceType& rFESpace, const ControlGridType& r_control_grid)
{
    // firstly get the values of all the basis functions
    std::vector<double> f_values;
    std::vector<typename TFESpaceType::LocalCoordinateType> xin(xi.size());
    std::copy(xi.begin(), xi.end(), xin.begin());
    rFESpace.GetValues(f_values, xin);

    // then interpolate the value at local coordinates using the control values
    v = f_values[0] * r_control_grid.GetData(0);
    for (std::size_t i = 1; i < r_control_grid.size(); ++i)
    {
        v += f_values[i] * r_control_grid.GetData(i);
    }
}

template<class TFESpaceType, typename TDataType, typename TLocalCoordinatesType>
void GridFunction_Helper<TFESpaceType, TDataType, TLocalCoordinatesType>::GetDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi,
                        const TFESpaceType& rFESpace, const ControlGridType& r_control_grid)
{
    constexpr int Dim = TFESpaceType::Dim();

    std::vector<typename TFESpaceType::LocalCoordinateType> xin(xi.size());
    std::copy(xi.begin(), xi.end(), xin.begin());

    // firstly get the values and derivatives of all the basis functions
    std::vector<std::vector<double> > f_derivatives;
    rFESpace.GetDerivatives(f_derivatives, xin);

    // then interpolate the derivative at local coordinates using the control values
    if (dv.size() != Dim)
    {
        dv.resize(Dim);
    }

    for (unsigned int dim = 0; dim < Dim; ++dim)
    {
        dv[dim] = f_derivatives[0][dim] * r_control_grid.GetData(0);
    }

    for (std::size_t i = 1; i < r_control_grid.size(); ++i)
    {
        for (unsigned int dim = 0; dim < Dim; ++dim)
        {
            dv[dim] += f_derivatives[i][dim] * r_control_grid.GetData(i);
        }
    }
}

template<class TFESpaceType, typename TDataType, typename TLocalCoordinatesType>
void GridFunction_Helper<TFESpaceType, TDataType, TLocalCoordinatesType>::GetSecondDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi,
                        const TFESpaceType& rFESpace, const ControlGridType& r_control_grid)
{
    constexpr int Dim = TFESpaceType::Dim();

    std::vector<typename TFESpaceType::LocalCoordinateType> xin(xi.size());
    std::copy(xi.begin(), xi.end(), xin.begin());

    // firstly get the values and derivatives of all the basis functions
    std::vector<std::vector<std::vector<double> > > f_derivatives;
    rFESpace.GetDerivatives(2, f_derivatives, xin);

    // then interpolate the second derivative at local coordinates using the control values
    const unsigned int size = Dim*(Dim+1)/2;
    if (dv.size() != size)
    {
        dv.resize(size);
    }

    for (unsigned int dim = 0; dim < size; ++dim)
    {
        dv[dim] = f_derivatives[1][0][dim] * r_control_grid.GetData(0);
    }

    for (std::size_t i = 1; i < r_control_grid.size(); ++i)
    {
        for (unsigned int dim = 0; dim < size; ++dim)
        {
            dv[dim] += f_derivatives[1][i][dim] * r_control_grid.GetData(i);
        }
    }
}

template<typename TDataType, typename TLocalCoordinateType>
struct GridFunction_Array1D_Predict_Helper<1, TDataType, TLocalCoordinateType>
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<TDataType, 3>& v, array_1d<TLocalCoordinateType, 3>& xi, const std::vector<int>& nsampling)
    {
        if (nsampling.size() < 1)
            KRATOS_ERROR << "sampling array must have dimension of at least 1";

        array_1d<TLocalCoordinateType, 3> xi0;
        array_1d<TDataType, 3> p;
        xi0[1] = 0;
        xi0[2] = 0;
        double dist, min_dist = 1.0e99;
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi0[0] = ((TLocalCoordinateType) i) / nsampling[0];
            noalias(p) = rGridFunc.GetValue(xi0);
            dist = std::abs(norm_2(p - v));
            if (dist < min_dist)
            {
                noalias(xi) = xi0;
                min_dist = dist;
            }
        }
    }
};

template<typename TDataType, typename TLocalCoordinateType>
struct GridFunction_Array1D_Predict_Helper<2, TDataType, TLocalCoordinateType>
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<TDataType, 3>& v, array_1d<TLocalCoordinateType, 3>& xi, const std::vector<int>& nsampling)
    {
        if (nsampling.size() < 2)
            KRATOS_ERROR << "sampling array must have dimension of at least 2";

        array_1d<TLocalCoordinateType, 3> xi0;
        array_1d<TDataType, 3> p;
        xi0[2] = 0.0;
        double dist, min_dist = 1.0e99;
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi0[0] = ((TLocalCoordinateType) i) / nsampling[0];
            for (int j = 0; j < nsampling[1] + 1; ++j)
            {
                xi0[1] = ((TLocalCoordinateType) j) / nsampling[1];
                noalias(p) = rGridFunc.GetValue(xi0);
                dist = std::abs(norm_2(p - v));
                if (dist < min_dist)
                {
                    noalias(xi) = xi0;
                    min_dist = dist;
                }
            }
        }
    }
};

template<typename TDataType, typename TLocalCoordinateType>
struct GridFunction_Array1D_Predict_Helper<3, TDataType, TLocalCoordinateType>
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<TDataType, 3>& v, array_1d<TLocalCoordinateType, 3>& xi, const std::vector<int>& nsampling)
    {
        if (nsampling.size() < 3)
            KRATOS_ERROR << "sampling array must have dimension of at least 3";

        array_1d<TLocalCoordinateType, 3> xi0;
        array_1d<TDataType, 3> p;
        double dist, min_dist = 1.0e99;
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi0[0] = ((TLocalCoordinateType) i) / nsampling[0];
            for (int j = 0; j < nsampling[1] + 1; ++j)
            {
                xi0[1] = ((TLocalCoordinateType) j) / nsampling[1];
                for (int k = 0; k < nsampling[2] + 1; ++k)
                {
                    xi0[2] = ((TLocalCoordinateType) k) / nsampling[2];
                    noalias(p) = rGridFunc.GetValue(xi0);
                    dist = std::abs(norm_2(p - v));
                    if (dist < min_dist)
                    {
                        noalias(xi) = xi0;
                        min_dist = dist;
                    }
                }
            }
        }
    }
};

template<>
struct GridFunction_Predict_Helper<1, array_1d<double, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<double, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<1, double, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

template<>
struct GridFunction_Predict_Helper<2, array_1d<double, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<double, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<2, double, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

template<>
struct GridFunction_Predict_Helper<3, array_1d<double, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<double, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<3, double, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

//

template<>
struct GridFunction_Predict_Helper<1, array_1d<std::complex<double>, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<std::complex<double>, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<1, std::complex<double>, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

template<>
struct GridFunction_Predict_Helper<2, array_1d<std::complex<double>, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<std::complex<double>, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<2, std::complex<double>, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

template<>
struct GridFunction_Predict_Helper<3, array_1d<std::complex<double>, 3>, array_1d<double, 3> >
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<std::complex<double>, 3>& v, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        GridFunction_Array1D_Predict_Helper<3, std::complex<double>, double>::Execute(rGridFunc, v, xi, nsampling);
    }
};

template<class TGridFunctionType, typename TDataType, typename TLocalCoordinatesType>
int GridFunction_Scalar_LocalCoordinates_Helper<TGridFunctionType, TDataType, TLocalCoordinatesType>::Execute(
        const TGridFunctionType& rGridFunc, const TDataType& v,
        TLocalCoordinatesType& xi, const int max_iters, const double TOL)
{
    KRATOS_ERROR << "To be implemented";

    // TODO

    // typedef Vector VectorType;
    // typedef Matrix MatrixType;

    // TDataType val;
    // std::vector<TDataType> ders(TDim);

    // VectorType res(TDim), dxi(TDim);
    // MatrixType J(TDim, TDim), InvJ(TDim, TDim);
    // double DetJ;

    // int it = 0;
    // bool converged = false;

    // do
    // {
    //     this->GetValue(val, xi);
    //     noalias(res) = v - val;
    //     if (norm_2(res) < TOL)
    //     {
    //         break;
    //     }

    //     this->GetDerivative(ders, xi);

    //     for (std::size_t i = 0; i < TDim; ++i)
    //     {
    //         for (std::size_t j = 0; j < TDim; ++j)
    //         {
    //             J(i, j) = ders[j][i];
    //         }
    //     }

    //     MathUtils<double>::InvertMatrix(J, InvJ, DetJ);
    //     noalias(dxi) = prod(InvJ, res);
    //     for (std::size_t i = 0; i < TDim; ++i)
    //         xi[i] += dxi[i];
    // }
    // while (++it < max_iters);

    // if ((it >= max_iters) && !converged)
    // {
    //     return 1;
    // }

    // return 0;
}

template<class TGridFunctionType, typename TVectorType, typename TLocalCoordinatesType>
int GridFunction_Vector_LocalCoordinates_Helper<TGridFunctionType, TVectorType, TLocalCoordinatesType>::Execute(
        const TGridFunctionType& rGridFunc, const TVectorType& v,
        TLocalCoordinatesType& xi, const int max_iters, const double TOL)
{
    typedef typename TVectorType::value_type DataType;

    typedef typename MatrixVectorTypeSelector<DataType>::VectorType VectorType;
    typedef typename MatrixVectorTypeSelector<DataType>::MatrixType MatrixType;

    constexpr int Dim = TGridFunctionType::Dim;

    TVectorType val;
    std::vector<TVectorType> ders(Dim);

    VectorType res(Dim), dxi(Dim);
    MatrixType J(Dim, Dim), InvJ(Dim, Dim);
    DataType DetJ;

    int it = 0;
    bool converged = false;

    do
    {
        rGridFunc.GetValue(val, xi);
        noalias(res) = v - val;
        if (std::abs(norm_2(res)) < TOL)
        {
            break;
        }

        rGridFunc.GetDerivative(ders, xi);

        for (std::size_t i = 0; i < Dim; ++i)
        {
            for (std::size_t j = 0; j < Dim; ++j)
            {
                J(i, j) = ders[j][i];
            }
        }

        MathUtils<DataType>::InvertMatrix(J, InvJ, DetJ);
        noalias(dxi) = prod(InvJ, res);
        for (std::size_t i = 0; i < Dim; ++i)
        {
            if constexpr (std::is_arithmetic<DataType>::value)
            {
                xi[i] += dxi[i];
            }
            else if constexpr (Internals::is_extended_arithmetic<DataType>::value)
            {
                if (std::abs(dxi[i].imag() > 1e-13))
                    KRATOS_ERROR << "The imaginary part is not zero in complex number Newton-Raphson";
                xi[i] += dxi[i].real(); // for complex number we use the real part, assuming the imaginary part is zero
            }
        }
    }
    while (++it < max_iters);

    if ((it >= max_iters) && !converged)
    {
        return 1;
    }

    return 0;
}

template<class TGridFunctionType, typename TDataType, typename TLocalCoordinatesType>
int GridFunction_LocalCoordinates_Helper<TGridFunctionType, TDataType, TLocalCoordinatesType>::Execute(
        const TGridFunctionType& rGridFunc, const TDataType& v,
        TLocalCoordinatesType& xi, const int max_iters, const double TOL)
{
    if constexpr (Internals::is_extended_arithmetic<TDataType>::value)
    {
        return GridFunction_Scalar_LocalCoordinates_Helper<TGridFunctionType, TDataType, TLocalCoordinatesType>::Execute(rGridFunc, v, xi, max_iters, TOL);
    }
    else
    {
        return GridFunction_Vector_LocalCoordinates_Helper<TGridFunctionType, TDataType, TLocalCoordinatesType>::Execute(rGridFunc, v, xi, max_iters, TOL);
    }
}

} // namespace Kratos.
