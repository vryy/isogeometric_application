//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/nurbs/structured_control_grid.h"


namespace Kratos
{

template<int TDim, typename TDataType>
struct ControlGridLibrary_Helper
{
    /// Generate regular control grid with a specific data type
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid);
};

/**
Helper library to generate control grid for isogeometric analysis
 */
class ControlGridLibrary
{
public:
    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(ControlGridLibrary);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;

    /// Default constructor
    ControlGridLibrary() {}

    /// Destructor
    virtual ~ControlGridLibrary() {}



    /// Generate the regular equidistant control point grid. All the point has unit weight.
    template<int TDim>
    static ControlGrid<ControlPointType>::Pointer CreateStructuredControlPointGrid(
            const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<double>& end)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }



    /// Generate the regular control point grid based on starting point and the director vector in each direction. All the point has unit weight.
    template<int TDim>
    static ControlGrid<ControlPointType>::Pointer CreateStructuredControlPointGrid(
            const std::vector<double>& start, const std::vector<std::size_t>& ngrid, const std::vector<std::vector<double> >& spacing_vectors)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }



    /// Generate regular control grid with a specific data type
    template<int TDim, typename TDataType>
    static typename ControlGrid<TDataType>::Pointer CreateStructuredZeroControlGrid(const std::string& Name, const std::vector<std::size_t>& ngrid)
    {
        return ControlGridLibrary_Helper<TDim, TDataType>::CreateStructuredZeroControlGrid(Name, ngrid);
    }



    /// Generate regular control grid with variable
    template<int TDim, class TVariableType>
    static typename ControlGrid<typename TVariableType::Type>::Pointer CreateStructuredZeroControlGrid(const TVariableType& rVariable, const std::vector<std::size_t>& ngrid)
    {
        return CreateStructuredZeroControlGrid<TDim, typename TVariableType::Type>(rVariable.Name(), ngrid);
    }



    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlGridLibrary";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const ControlGridLibrary& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#include "control_grid_library.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_LIBRARY_H_INCLUDED defined

