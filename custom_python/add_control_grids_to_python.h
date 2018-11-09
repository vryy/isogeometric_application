/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 11, 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED

// System includes
#include <string>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/unstructured_control_grid.h"
#include "custom_utilities/nurbs/structured_control_grid.h"
#include "custom_utilities/control_grid_library.h"
#include "custom_utilities/control_grid_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

////////////////////////////////////////

inline double ControlPoint_GetWX(ControlPoint<double>& rDummy)
{
    return rDummy.WX();
}

inline void ControlPoint_SetWX(ControlPoint<double>& rDummy, const double& newWX)
{
    rDummy.WX() = newWX;
}

inline double ControlPoint_GetWY(ControlPoint<double>& rDummy)
{
    return rDummy.WY();
}

inline void ControlPoint_SetWY(ControlPoint<double>& rDummy, const double& newWY)
{
    rDummy.WY() = newWY;
}

inline double ControlPoint_GetWZ(ControlPoint<double>& rDummy)
{
    return rDummy.WZ();
}

inline void ControlPoint_SetWZ(ControlPoint<double>& rDummy, const double& newWZ)
{
    rDummy.WZ() = newWZ;
}

inline double ControlPoint_GetW(ControlPoint<double>& rDummy)
{
    return rDummy.W();
}

inline void ControlPoint_SetW(ControlPoint<double>& rDummy, const double& newW)
{
    rDummy.W() = newW;
}

inline void ControlPoint_ApplyTransformation(ControlPoint<double>& rDummy, const Transformation<double>& trans)
{
    rDummy.ApplyTransformation(trans);
}

////////////////////////////////////////

template<typename TDataType>
TDataType ControlGrid_GetItem(ControlGrid<TDataType>& rDummy, int index)
{
    return rDummy.GetData(index);
}

template<typename TDataType>
void ControlGrid_SetItem(ControlGrid<TDataType>& rDummy, int index, const TDataType& value)
{
    rDummy.SetData(index, value);
}

////////////////////////////////////////

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateLinearControlPointGrid(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u,
        const double& end_x, const double& end_y, const double& end_z);

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid1(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& end_x, const double& end_y);

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateRectangularControlPointGrid2(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v,
        const double& space1_x, const double& space1_y, const double& space1_z,
        const double& space2_x, const double& space2_y, const double& space2_z);

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid1(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        const double& end_x, const double& end_y, const double& end_z);

ControlGrid<ControlPoint<double> >::Pointer ControlGridLibrary_CreateCubicControlPointGrid2(
        ControlGridLibrary& rDummy,
        const double& start_x, const double& start_y, const double& start_z,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w,
        boost::python::list spacing_vectors_data);

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateLinearZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u)
{
    std::vector<std::size_t> ngrid(1);
    ngrid[0] = n_points_u;
    return rDummy.CreateStructuredZeroControlGrid<1, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateRectangularZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u, const std::size_t& n_points_v)
{
    std::vector<std::size_t> ngrid(2);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    return rDummy.CreateStructuredZeroControlGrid<2, TVariableType>(rVariable, ngrid);
}

template<class TVariableType>
inline typename ControlGrid<typename TVariableType::Type>::Pointer ControlGridLibrary_CreateCubicZeroControlGridWithVariable(
        ControlGridLibrary& rDummy,
        const TVariableType& rVariable,
        const std::size_t& n_points_u, const std::size_t& n_points_v, const std::size_t& n_points_w)
{
    std::vector<std::size_t> ngrid(3);
    ngrid[0] = n_points_u;
    ngrid[1] = n_points_v;
    ngrid[2] = n_points_w;
    return rDummy.CreateStructuredZeroControlGrid<3, TVariableType>(rVariable, ngrid);
}

template<typename TDataType, class TFESpaceType>
inline typename ControlGrid<TDataType>::Pointer ControlGridUtility_CreatePointBasedControlGrid(
        ControlGridUtility& rDummy,
        const Variable<TDataType>& rVariable, typename TFESpaceType::Pointer pFESpace)
{
    return rDummy.CreatePointBasedControlGrid<TDataType, TFESpaceType>(rVariable, pFESpace);
}

}  // namespace Python.

} // Namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_CONTROL_GRIDS_TO_PYTHON_H_INCLUDED
