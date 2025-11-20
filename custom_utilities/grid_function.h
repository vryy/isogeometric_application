//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
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
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/fespace_utility.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/control_grid_utility.h"

namespace Kratos
{

template<class TFESpaceType, typename TDataType, typename TLocalCoordinatesType>
struct GridFunction_Helper
{
    typedef ControlGrid<TDataType> ControlGridType;

    static void GetValue(TDataType& v, const TLocalCoordinatesType& xi,
                    const TFESpaceType& rFESpace, const ControlGridType& r_control_grid);

    static void GetDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi,
                    const TFESpaceType& rFESpace, const ControlGridType& r_control_grid);

    static void GetSecondDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi,
                    const TFESpaceType& rFESpace, const ControlGridType& r_control_grid);
};

template<int TDim, typename TDataType, typename TLocalCoordinatesType>
struct GridFunction_Predict_Helper
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const TDataType& v, TLocalCoordinatesType& xi, const std::vector<int>& nsampling)
    {
        KRATOS_ERROR << "Not yet implemented";
    }
};

template<int TDim, typename TDataType, typename TLocalCoordinateType>
struct GridFunction_Array1D_Predict_Helper
{
    template<class TGridFunctionType>
    static void Execute(TGridFunctionType& rGridFunc,
                        const array_1d<TDataType, 3>& v, array_1d<TLocalCoordinateType, 3>& xi, const std::vector<int>& nsampling)
    {
        KRATOS_ERROR << "Not yet implemented";
    }
};

template<class TGridFunctionType, typename TDataType, typename TLocalCoordinatesType>
struct GridFunction_LocalCoordinates_Helper
{
    static int Execute(const TGridFunctionType& rGridFunc, const TDataType& v, TLocalCoordinatesType& xi, const int max_iters = 30, const double TOL = 1.0e-8);
};

template<class TGridFunctionType, typename TDataType, typename TLocalCoordinatesType>
struct GridFunction_Scalar_LocalCoordinates_Helper
{
    static int Execute(const TGridFunctionType& rGridFunc, const TDataType& v, TLocalCoordinatesType& xi, const int max_iters = 30, const double TOL = 1.0e-8);
};

template<class TGridFunctionType, typename TVectorType, typename TLocalCoordinatesType>
struct GridFunction_Vector_LocalCoordinates_Helper
{
    static int Execute(const TGridFunctionType& rGridFunc, const TVectorType& v, TLocalCoordinatesType& xi, const int max_iters = 30, const double TOL = 1.0e-8);
};

/**
 * A grid function is a function defined over the parametric domain. It takes the control values at grid point and interpolate the corresponding physical terms.
 */
template<int TDim, typename TLocalCoordinateType, typename TDataType>
class GridFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GridFunction);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const GridFunction> ConstPointer;
#endif

    /// Type definition
    typedef TDataType DataType;
    typedef std::vector<TDataType> DataContainerType;
    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;
    typedef ControlGrid<TDataType> ControlGridType;

    typedef GridFunction<TDim, TLocalCoordinateType, TDataType> ThisType;

    /// Constant
    static constexpr int Dim = TDim;

    /// Default constructor
    GridFunction(typename FESpaceType::Pointer pFESpace, typename ControlGridType::Pointer pControlGrid)
        : mpFESpace(pFESpace), mpControlGrid(pControlGrid) {}

    /// Destructor
    virtual ~GridFunction() {}

    /// Helper to create the new instance of grid function
    static typename ThisType::Pointer Create(typename FESpaceType::Pointer pFESpace, typename ControlGridType::Pointer pControlGrid)
    {
        return typename ThisType::Pointer(new ThisType(pFESpace, pControlGrid));
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

    /// Get the value of the grid function at specific local coordinates
    /// The function values to interpolate the grid function value are provided by FESpace. In the case of NURBS,
    /// either FESpace or grid function value must contain weight information. For examples, if TDataType is ControlPoint,
    /// the FESpace must be an unweighted one.
    template<typename TLocalCoordinatesType>
    void GetValue(TDataType& v, const TLocalCoordinatesType& xi) const
    {
        GridFunction_Helper<FESpaceType, TDataType, TLocalCoordinatesType>::GetValue(v, xi, *pFESpace(), *pControlGrid());
    }

    /// Get the value of the grid function at specific local coordinates
    template<typename TLocalCoordinatesType>
    TDataType GetValue(const TLocalCoordinatesType& xi) const
    {
        TDataType v;
        this->GetValue(v, xi);
        return v;
    }

    /// Get the derivatives of the grid function at specific local coordinates
    /// The output values has the form: [d(values(xi)) / d(xi_0), d(values(xi)) / d(xi_1), ...]
    /// The function derivatives to interpolate the grid function value are provided by FESpace. Hence, the TDataType must
    /// be unweighted type, in order to have correct derivative values. Because, homogeous transformation
    /// is not applied for derivatives. If TDataType is a weighted type, e.g. CONTROL_POINT, one must change to
    /// use the unweighted one, e.g. CONTROL_POINT_COORDINATES.
    template<typename TLocalCoordinatesType>
    void GetDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi) const
    {
        GridFunction_Helper<FESpaceType, TDataType, TLocalCoordinatesType>::GetDerivative(dv, xi, *pFESpace(), *pControlGrid());
    }

    /// Get the derivatives of the grid function at specific local coordinates
    /// The return values has the form: [d_values(xi) / d_xi_0, d_values(xi) / d_xi_1, ...]
    template<typename TLocalCoordinatesType>
    std::vector<TDataType> GetDerivative(const TLocalCoordinatesType& xi) const
    {
        std::vector<TDataType> dv(TDim);
        this->GetDerivative(dv, xi);
        return dv;
    }

    /// Get the second derivatives of the grid function at specific local coordinates
    /// The output values has the form:
    /// +   in 1D: [d2(values(xi)) / d(xi_0) d(xi_0)]
    /// +   in 2D: [d2(values(xi)) / d(xi_0) d(xi_0), d2(values(xi)) / d(xi_1) d(xi_1), d2(values(xi)) / d(xi_0) d(xi_1)]
    /// +   in 3D: [d2(values(xi)) / d(xi_0) d(xi_0), d2(values(xi)) / d(xi_1) d(xi_1), d2(values(xi)) / d(xi_2) d(xi_2), d2(values(xi)) / d(xi_0) d(xi_1), d2(values(xi)) / d(xi_0) d(xi_2), , d2(values(xi)) / d(xi_1) d(xi_2)]
    /// The function derivatives to interpolate the grid function value are provided by FESpace. Hence, the TDataType must
    /// be unweighted type, in order to have correct derivative values. Because, homogeous transformation
    /// is not applied for derivatives. If TDataType is a weighted type, e.g. CONTROL_POINT, one must change to
    /// use the unweighted one, e.g. CONTROL_POINT_COORDINATES.
    template<typename TLocalCoordinatesType>
    void GetSecondDerivative(std::vector<TDataType>& dv, const TLocalCoordinatesType& xi) const
    {
        GridFunction_Helper<FESpaceType, TDataType, TLocalCoordinatesType>::GetSecondDerivative(dv, xi, *pFESpace(), *pControlGrid());
    }

    /// Get the second derivatives of the grid function at specific local coordinates
    /// The output values has the form:
    /// +   in 1D: [d2(values(xi)) / d(xi_0) d(xi_0)]
    /// +   in 2D: [d2(values(xi)) / d(xi_0) d(xi_0), d2(values(xi)) / d(xi_1) d(xi_1), d2(values(xi)) / d(xi_0) d(xi_1)]
    /// +   in 3D: [d2(values(xi)) / d(xi_0) d(xi_0), d2(values(xi)) / d(xi_1) d(xi_1), d2(values(xi)) / d(xi_2) d(xi_2), d2(values(xi)) / d(xi_0) d(xi_1), d2(values(xi)) / d(xi_0) d(xi_2), , d2(values(xi)) / d(xi_1) d(xi_2)]
    template<typename TLocalCoordinatesType>
    std::vector<TDataType> GetSecondDerivative(const TLocalCoordinatesType& xi) const
    {
        std::vector<TDataType> dv;
        GridFunction_Helper<FESpaceType, TDataType, TLocalCoordinatesType>::GetSecondDerivative(dv, xi, *pFESpace(), *pControlGrid());
        return dv;
    }

    /// Compute a prediction for LocalCoordinates algorithm. Because LocalCoordinates uses Newton-Raphson algorithm to compute
    /// the inversion, it requires a good initial starting point
    template<typename TLocalCoordinatesType>
    void Predict(const TDataType& v, TLocalCoordinatesType& xi, const std::vector<int>& nsampling,
                 const TLocalCoordinatesType& xi_min, const TLocalCoordinatesType& xi_max) const
    {
        std::cout << "WARNING!!!Predict on range {xi_min, xi_max} is not yet implemented. {0, 1} is used for now." << std::endl;
        GridFunction_Predict_Helper<TDim, TDataType, TLocalCoordinatesType>::Execute(*this, v, xi, nsampling);
    }

    /// Compute a prediction for LocalCoordinates algorithm. Because LocalCoordinates uses Newton-Raphson algorithm to compute
    /// the inversion, it requires a good initial starting point
    template<typename TLocalCoordinatesType>
    void Predict(const TDataType& v, TLocalCoordinatesType& xi, const std::vector<int>& nsampling) const
    {
        GridFunction_Predict_Helper<TDim, TDataType, TLocalCoordinatesType>::Execute(*this, v, xi, nsampling);
    }

    /// Compute the local coordinate of point that has a specific interpolated values
    /// It is noted that this function only works with TDataType and TLocalCoordinatesType as a Vector-compatible type
    /// Assumming the size of TLocalCoordinatesType is TDim
    /// On the output:
    ///     0: the local NR converged successfully, rLocalCoordinates is the found point
    ///     1: the NR does not converge
    template<typename TLocalCoordinatesType>
    int LocalCoordinates(const TDataType& v, TLocalCoordinatesType& xi, const int max_iters = 30, const double TOL = 1.0e-8) const
    {
        return GridFunction_LocalCoordinates_Helper<ThisType, TDataType, TLocalCoordinatesType>::Execute(*this, v, xi, max_iters, TOL);
    }

    /// Check the compatibility between the underlying control grid and fe space.
    bool Validate() const
    {
        if (mpFESpace == NULL)
        {
            KRATOS_ERROR << "The FESpace is not defined for " << Info();
            return false;
        }

        if (mpControlGrid == NULL)
        {
            KRATOS_ERROR << "The control grid is not defined for " << Info();
            return false;
        }

        if (mpFESpace->TotalNumber() != mpControlGrid->Size())
        {
            KRATOS_ERROR << "The control grid and the FESpace does not have the same size";
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

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("FESpaceType", mpFESpace->Type());
        rSerializer.save("FESpace", *mpFESpace);
        rSerializer.save("ControlGridType", mpControlGrid->Type());
        rSerializer.save("ControlGrid", *mpControlGrid);
    }

    virtual void load(Serializer& rSerializer)
    {
        std::string fespace_type;
        rSerializer.load("FESpaceType", fespace_type);
        mpFESpace = FESpaceUtility<TDim, TLocalCoordinateType>::CreateEmptyFESpace(fespace_type);
        rSerializer.load("FESpace", *mpFESpace);

        std::string control_grid_type;
        rSerializer.load("ControlGridType", control_grid_type);
        mpControlGrid = ControlGridUtility::CreateEmptyControlGrid<TDim, TDataType>(control_grid_type);
        rSerializer.load("ControlGrid", *mpControlGrid);
    }
    ///@}

};

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const GridFunction<TDim, TLocalCoordinateType, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#include "grid_function.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_GRID_FUNCTION_H_INCLUDED defined
