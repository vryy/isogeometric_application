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
#include "custom_utilities/fespace.h"
#include "custom_utilities/control_grid.h"

namespace Kratos
{

/**
A grid function is a function defined over the parametric domain. It takes the control values at grid point and interpolate the corresponding physical terms.
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
    template<typename TCoordinatesType>
    TDataType GetValue(const TCoordinatesType& xi) const
    {
        // firstly get the values of all the basis functions
        std::vector<double> f_values;
        pFESpace()->GetValue(f_values, xi);

        // then interpolate the value at local coordinates using the control values
        const ControlGridType& r_control_grid = *pControlGrid();

        TDataType v = f_values[0] * r_control_grid.GetData(0);
        for (std::size_t i = 1; i < r_control_grid.size(); ++i)
            v += f_values[i] * r_control_grid.GetData(i);

        return v;
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
        ss << "GridFunction" << TDim << "D";
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
