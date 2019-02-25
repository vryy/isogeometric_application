//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/control_grid.h"

namespace Kratos
{

// forward declaration
template<int TDim, typename TDataType> class GridFunction;

template<int TDim, typename TDataType, typename TCoordinatesType>
struct GridFunction_Helper
{
    typedef FESpace<TDim> FESpaceType;
    typedef ControlGrid<TDataType> ControlGridType;

    static void GetValue(TDataType& v, const TCoordinatesType& xi,
        const FESpaceType& rFESpace, const ControlGridType& r_control_grid)
    {
        // firstly get the values of all the basis functions
        std::vector<double> f_values;
        std::vector<double> xin(xi.size());
        std::copy(xi.begin(), xi.end(), xin.begin());
        rFESpace.GetValue(f_values, xin);

        // then interpolate the value at local coordinates using the control values
        v = f_values[0] * r_control_grid.GetData(0);
        for (std::size_t i = 1; i < r_control_grid.size(); ++i)
            v += f_values[i] * r_control_grid.GetData(i);
    }

    static void GetDerivative(std::vector<TDataType>& dv, const TCoordinatesType& xi,
        const FESpaceType& rFESpace, const ControlGridType& r_control_grid)
    {
        std::vector<double> xin(xi.size());
        std::copy(xi.begin(), xi.end(), xin.begin());

        // firstly get the values and derivatives of all the basis functions
        std::vector<std::vector<double> > f_derivatives;
        rFESpace.GetDerivative(f_derivatives, xin);

        // rFESpace.PrintInfo(std::cout); std::cout << std::endl;

        // for (int i = 0; i < f_derivatives.size(); ++i)
        // {
        //     std::cout << "f_derivatives[" << i << "]:";
        //     for (int j = 0; j < f_derivatives[i].size(); ++j)
        //         std::cout << " " << f_derivatives[i][j];
        //     std::cout << std::endl;
        // }

        // then interpolate the derivative at local coordinates using the control values
        if (dv.size() != TDim)
            dv.resize(TDim);

        for (int dim = 0; dim < TDim; ++dim)
            dv[dim] = f_derivatives[0][dim] * r_control_grid.GetData(0);

        for (std::size_t i = 1; i < r_control_grid.size(); ++i)
        {
            for (int dim = 0; dim < TDim; ++dim)
                dv[dim] += f_derivatives[i][dim] * r_control_grid.GetData(i);
        }
    }
};

template<int TDim, typename TDataType>
struct GridFunction_Helper<TDim, TDataType, std::vector<double> >
{
    typedef FESpace<TDim> FESpaceType;
    typedef ControlGrid<TDataType> ControlGridType;

    static void GetValue(TDataType& v, const std::vector<double>& xi,
        const FESpaceType& rFESpace, const ControlGridType& r_control_grid)
    {
        // firstly get the values of all the basis functions
        std::vector<double> f_values;
        rFESpace.GetValue(f_values, xi);

        // then interpolate the value at local coordinates using the control values
        v = f_values[0] * r_control_grid.GetData(0);
        for (std::size_t i = 1; i < r_control_grid.size(); ++i)
            v += f_values[i] * r_control_grid.GetData(i);
    }

    static void GetDerivative(std::vector<TDataType>& dv, const std::vector<double>& xi,
        const FESpaceType& rFESpace, const ControlGridType& r_control_grid)
    {
        // firstly get the values and derivatives of all the basis functions
        std::vector<std::vector<double> > f_derivatives;
        rFESpace.GetDerivative(f_derivatives, xi);

        // rFESpace.PrintInfo(std::cout); std::cout << std::endl;

        // for (int i = 0; i < f_derivatives.size(); ++i)
        // {
        //     std::cout << "f_derivatives[" << i << "]:";
        //     for (int j = 0; j < f_derivatives[i].size(); ++j)
        //         std::cout << " " << f_derivatives[i][j];
        //     std::cout << std::endl;
        // }

        // then interpolate the derivative at local coordinates using the control values
        if (dv.size() != TDim)
            dv.resize(TDim);

        for (int dim = 0; dim < TDim; ++dim)
            dv[dim] = f_derivatives[0][dim] * r_control_grid.GetData(0);

        for (std::size_t i = 1; i < r_control_grid.size(); ++i)
        {
            for (int dim = 0; dim < TDim; ++dim)
                dv[dim] += f_derivatives[i][dim] * r_control_grid.GetData(i);
        }
    }
};

/**
 * A grid function is a function defined over the parametric domain. It takes the control values at grid point and interpolate the corresponding physical terms.
 */
template<int TDim, typename TDataType>
class GridFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);

    /// Type definition
    typedef TDataType DataType;
    typedef std::vector<TDataType> DataContainerType;
    typedef FESpace<TDim> FESpaceType;
    typedef ControlGrid<TDataType> ControlGridType;

    /// Default constructor
    GridFunction(typename FESpaceType::Pointer pFESpace, typename ControlGridType::Pointer pControlGrid)
    : mpFESpace(pFESpace), mpControlGrid(pControlGrid) {}

    /// Destructor
    virtual ~GridFunction() {}

    /// Helper to create the new instance of grid function
    static typename GridFunction<TDim, TDataType>::Pointer Create(typename FESpaceType::Pointer pFESpace, typename ControlGridType::Pointer pControlGrid)
    {
        return typename GridFunction<TDim, TDataType>::Pointer(new GridFunction<TDim, TDataType>(pFESpace, pControlGrid));
    }

    /// Set the FESpace
    void SetFESpace(typename FESpaceType::Pointer pNewFESpace) {mpFESpace = pNewFESpace;} // use this with care

    /// Get the FESpace pointer
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Get the FESpace pointer
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Set the control grid
    /// Remarks: this function will effectively replace the underlying control grid, so use this with care
    void SetControlGrid(typename ControlGridType::Pointer pNewControlGrid) {mpControlGrid = pNewControlGrid;}

    /// Get the control grid pointer
    typename ControlGridType::Pointer pControlGrid() {return mpControlGrid;}

    /// Get the control grid pointer
    typename ControlGridType::ConstPointer pControlGrid() const {return mpControlGrid;}

    /// Get the value of the grid at specific local coordinates
    /// The function values to interpolate the grid value are provided by FESpace. In the case of NURBS,
    /// either FESpace or grid value must contain weight information. For examples, if TDataType is CONTROL_POINT,
    /// the FESpace must be an unweighted one.
    template<typename TCoordinatesType>
    void GetValue(TDataType& v, const TCoordinatesType& xi) const
    {
        GridFunction_Helper<TDim, TDataType, TCoordinatesType>::GetValue(v, xi, *pFESpace(), *pControlGrid());
    }

    /// Get the value of the grid at specific local coordinates
    template<typename TCoordinatesType>
    TDataType GetValue(const TCoordinatesType& xi) const
    {
        TDataType v;
        this->GetValue(v, xi);
        return v;
    }

    /// Get the derivatives of the grid at specific local coordinates
    /// The output values has the form: [d(values(xi)) / d(xi_0), d(values(xi)) / d(xi_1), ...]
    /// The function derivatives to interpolate the grid value are provided by FESpace. Hence, the TDataType must
    /// be unweighted type, in order to have correct derivative values. Because, homogeous transformation
    /// is not applied for derivatives. If TDataType is a weighted type, e.g. CONTROL_POINT, one must change to
    /// use the unweighted one, e.g. CONTROL_POINT_COORDINATES.
    template<typename TCoordinatesType>
    void GetDerivative(std::vector<TDataType>& dv, const TCoordinatesType& xi) const
    {
        GridFunction_Helper<TDim, TDataType, TCoordinatesType>::GetDerivative(dv, xi, *pFESpace(), *pControlGrid());
    }

    /// Get the derivatives of the grid at specific local coordinates
    /// The return values has the form: [d_values(xi) / d_xi_0, d_values(xi) / d_xi_1, ...]
    template<typename TCoordinatesType>
    std::vector<TDataType> GetDerivative(const TCoordinatesType& xi) const
    {
        std::vector<TDataType> dv(TDim);
        this->GetDerivative(dv, xi);
        return dv;
    }

    /// Compute the local coordinate of point that has a specific interpolated values
    /// It is noted that this function only works with TDataType and TCoordinatesType as a Vector-compatible type
    /// On the output:
    ///     0: the local NR converged successfully, rLocalCoordinates is the found point
    ///     1: the NR does not converge
    template<typename TCoordinatesType>
    int LocalCoordinates(const TDataType& v, TCoordinatesType& xi) const
    {
        typedef Vector VectorType;
        typedef Matrix MatrixType;

        TDataType val;
        std::vector<TDataType> ders(TDim);

        VectorType res(TDim), dxi(TDim);
        MatrixType J(TDim, TDim), InvJ(TDim, TDim);
        double DetJ;

        int it = 0;
        const int max_iters = 30;
        const double TOL = 1.0e-8;
        bool converged = false;

        do
        {
            this->GetValue(val, xi);
            noalias(res) = v - val;
            if (norm_2(res) < TOL)
                break;

            this->GetDerivative(ders, xi);

            for (std::size_t i = 0; i < TDim; ++i)
                for (std::size_t j = 0; j < TDim; ++j)
                    J(i, j) = ders[j][i];

            MathUtils<double>::InvertMatrix(J, InvJ, DetJ);
            noalias(dxi) = prod(InvJ, res);
            xi += dxi;
        } while (++it < max_iters);

        if ((it >= max_iters) && !converged)
            return 1;

        return 0;
    }

    /// Check the compatibility between the underlying control grid and fe space.
    bool Validate() const
    {
        if (mpFESpace == NULL)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The FESpace is not defined for ", Info())
            return false;
        }

        if (mpControlGrid == NULL)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The control grid is not defined for ", Info())
            return false;
        }

        if (mpFESpace->TotalNumber() != mpControlGrid->Size())
        {
            KRATOS_THROW_ERROR(std::logic_error, "The control grid and the FESpace does not have the same size", "")
            return false;
        }

        return true;
    }

    /// Information
    std::string Info() const
    {
        std::stringstream ss;
        ss << "GridFunction" << TDim << "D_" << mpControlGrid->Name();
        return ss.str();
    }

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "<<<Listing of grid function " << mpControlGrid->Name() << ":" << std::endl;
        rOStream << "-----FESPace:" << std::endl;
        rOStream << *mpFESpace << std::endl;
        rOStream << "-----Control point grid:" << std::endl;
        rOStream << *mpControlGrid << std::endl;
        rOStream << ">>>End Listing of grid function " << mpControlGrid->Name() << std::endl;
    }

private:

    typename FESpaceType::Pointer mpFESpace;
    typename ControlGridType::Pointer mpControlGrid;

};


/// output stream function
template<int TDim, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const GridFunction<TDim, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED defined
