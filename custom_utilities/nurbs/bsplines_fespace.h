//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED

// System includes
#include <vector>
#include <algorithm>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/bcell.h"
#include "custom_utilities/nurbs/bcell_manager.h"

// #define DEBUG_GEN_CELL

namespace Kratos
{

template<int TDim, typename TLocalCoordinateType>
class BSplinesFESpace;

template<int TDim>
struct BSplinesFESpace_Helper
{
    /// Get the values of the basis functions at point xi
    template<typename TLocalCoordinateType>
    static void GetValues(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
                          std::vector<double>& values, const std::vector<TLocalCoordinateType>& xi);

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    template<typename TLocalCoordinateType>
    static void GetValuesAndDerivatives(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
                                        std::vector<double>& values,
                                        std::vector<std::vector<double> >& derivatives,
                                        const std::vector<TLocalCoordinateType>& xi);

    /// the output derivatives has the form of values[der_index][func_index][param_index]
    /// For the first derivatives, param_index is 0,1,..,dim
    /// For the second derivatives, param_index is 00,11,..,01,02,..,11,12,..
    /// For the higher derivatives, param_index is [d1][d2][d3] (d1<=d2<=d3) (TODO)
    template<typename TLocalCoordinateType>
    static void GetValuesAndDerivatives(const BSplinesFESpace<TDim, TLocalCoordinateType>& rFESpace,
                                        const unsigned int nd,
                                        std::vector<double>& values,
                                        std::vector<std::vector<std::vector<double> > >& derivatives,
                                        const std::vector<TLocalCoordinateType>& xi);
};

/**
This class represents the FESpace for a single BSplines patch defined over parametric domain.
 */
template<int TDim, typename TLocalCoordinateType = double>
class BSplinesFESpace : public FESpace<TDim, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesFESpace);

    /// Type definition
    typedef FESpace<TDim, TLocalCoordinateType> BaseType;
    typedef BSplinesFESpace<TDim, TLocalCoordinateType> ThisType;
    typedef KnotArray1D<TLocalCoordinateType> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef BCellManager<TDim, BCell> cell_container_t;

    /// Default constructor
    BSplinesFESpace() : BaseType() {}

    /// Destructor
    ~BSplinesFESpace() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << this->Type() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper to create new BSplinesFESpace pointer
    static typename ThisType::Pointer Create()
    {
        return typename ThisType::Pointer(new ThisType());
    }

    /// Get the order of the BSplines patch in specific direction
    std::size_t Order(std::size_t idir) const final
    {
        if (idir >= TDim) { return 0; }
        else { return mOrders[idir]; }
    }

    /// Get the number of control points of the BSplines in all direction
    std::vector<std::size_t> Numbers() const
    {
        std::vector<std::size_t> numbers(TDim);
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            numbers[dim] = mNumbers[dim];
        }
        return numbers;
    }

    /// Get the number of control points of the BSplines in all direction except the specified one
    /// This function is useful to extract the number of shape functions in each direction of a sliced patch
    std::vector<std::size_t> Numbers(std::size_t idir) const
    {
        std::vector<std::size_t> numbers(TDim - 1);
        std::size_t cnt = 0;
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            if (dim != idir)
            {
                numbers[cnt++] = mNumbers[dim];
            }
        }
        return numbers;
    }

    /// Get the number of control points of the BSplines in specific direction
    std::size_t Number(std::size_t idir) const {return mNumbers[idir];}

    /// Get the number of basis functions defined over the BSplines
    std::size_t TotalNumber() const final
    {
        std::size_t Number = 1;
        for (std::size_t i = 0; i < TDim; ++i)
        {
            Number *= mNumbers[i];
        }
        return Number;
    }

    /// Get the lower and upper bound of the parametric space in a specific direction
    std::vector<TLocalCoordinateType> ParametricBounds(std::size_t idir) const final
    {
        std::vector<TLocalCoordinateType> bound(2);
        bound[0] = (*(mKnotVectors[idir].begin()))->Value();
        bound[1] = (*(mKnotVectors[idir].end() - 1))->Value();
        return bound;
    }

    /// Get the string representing the type of the patch
    std::string Type() const final
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "BSplinesFESpace" << TDim << "D";
        return ss.str();
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(std::size_t idir, const knot_container_t& p_knot_vector)
    {
        mKnotVectors[idir] = p_knot_vector;
    }

    /// Create and set the knot vector in direction i.
    void SetKnotVector(std::size_t idir, const std::vector<TLocalCoordinateType>& values)
    {
        if (idir >= TDim)
        {
            KRATOS_ERROR << "Invalid direction " << idir;
        }
        else
        {
            mKnotVectors[idir].clear();
            for (std::size_t j = 0; j < values.size(); ++j)
            {
                mKnotVectors[idir].pCreateKnot(values[j]);
            }
        }
    }

    /// Get the knot vector in i-direction
    const knot_container_t& KnotVector(std::size_t idir) const {return mKnotVectors[idir];}

    /// Reverse the evaluation in i-direction
    void Reverse(std::size_t idir) override
    {
        // reverse the knot vector
        mKnotVectors[idir].Reverse();

        // also change the function indices
        BSplinesIndexingUtility::Reverse<TDim, std::vector<std::size_t>, std::vector<std::size_t> >(mFunctionsIds, this->Numbers(), idir);

        // and the global to local map
        BaseType::mGlobalToLocal.clear();
        for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
        {
            BaseType::mGlobalToLocal[mFunctionsIds[i]] = i;
        }
    }

    /// Transpose (swap) the i- and j- direction
    void Transpose(std::size_t idir, std::size_t jdir) override
    {
        // reverse the knot vector
        const auto kvi = mKnotVectors[idir];
        const auto kvj = mKnotVectors[jdir];

        std::size_t order_i = mOrders[idir];
        std::size_t order_j = mOrders[jdir];

        std::size_t number_i = mNumbers[idir];
        std::size_t number_j = mNumbers[jdir];

        mKnotVectors[idir] = kvj;
        mKnotVectors[jdir] = kvi;

        mOrders[idir] = order_j;
        mOrders[jdir] = order_i;

        mNumbers[idir] = number_j;
        mNumbers[jdir] = number_i;

        // also change the function indices
        BSplinesIndexingUtility::Transpose<TDim, std::vector<std::size_t>, std::vector<std::size_t> >(mFunctionsIds, this->Numbers(), idir, jdir);

        // and the global to local map
        BaseType::mGlobalToLocal.clear();
        for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
        {
            BaseType::mGlobalToLocal[mFunctionsIds[i]] = i;
        }
    }

    /// Set the BSplines information in the direction i
    void SetInfo(std::size_t idir, std::size_t Number, std::size_t Order)
    {
        mOrders[idir] = Order;
        mNumbers[idir] = Number;
    }

    /// Validate the BSplinesFESpace
    bool Validate() const final
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (mKnotVectors[i].size() != mNumbers[i] + mOrders[i] + 1)
            {
                KRATOS_ERROR << "The knot vector is incompatible at dimension " << i;
            }
        }

        return BaseType::Validate();
    }

    /// Get the values of the basis function i at point xi
    void GetValue(double& v, std::size_t i, const std::vector<TLocalCoordinateType>& xi) const final
    {
        // TODO the current approach is expensive (all is computed). We have to find the way to optimize it.
        std::vector<double> values;
        this->GetValues(values, xi);
        v = values[i];
    }

    /// Get the values of the basis functions at point xi
    void GetValues(std::vector<double>& values, const std::vector<TLocalCoordinateType>& xi) const final
    {
        BSplinesFESpace_Helper<TDim>::template GetValues(*this, values, xi);
    }

    /// Get the derivatives of the basis function i at point xi
    void GetDerivative(std::vector<double>& values, std::size_t i, const std::vector<TLocalCoordinateType>& xi) const final
    {
        // TODO the current approach is expensive (all is computed). Find the way to optimize it.
        std::vector<std::vector<double> > tmp;
        this->GetDerivatives(tmp, xi);
        values = tmp[i];
    }

    /// Get the derivatives of the basis functions at point xi
    void GetDerivatives(std::vector<std::vector<double> >& values, const std::vector<TLocalCoordinateType>& xi) const final
    {
        std::vector<double> dummy;
        this->GetValuesAndDerivatives(dummy, values, xi);
    }

    /// [derived]
    void GetDerivatives(const unsigned int nd, std::vector<std::vector<std::vector<double> > >& values, const std::vector<TLocalCoordinateType>& xi) const final
    {
        std::vector<double> dummy;
        this->GetValuesAndDerivatives(nd, dummy, values, xi);
    }

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    void GetValuesAndDerivatives(std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<TLocalCoordinateType>& xi) const final
    {
        BSplinesFESpace_Helper<TDim>::template GetValuesAndDerivatives(*this, values, derivatives, xi);
    }

    /// [derived]
    void GetValuesAndDerivatives(const unsigned int nd, std::vector<double>& values, std::vector<std::vector<std::vector<double> > >& derivatives, const std::vector<TLocalCoordinateType>& xi) const final
    {
        BSplinesFESpace_Helper<TDim>::template GetValuesAndDerivatives(*this, nd, values, derivatives, xi);
    }

    /// Check if a point lies inside the parametric domain of the BSplinesFESpace
    bool IsInside(const std::vector<TLocalCoordinateType>& xi) const final
    {
        bool is_inside = true;
        for (std::size_t i = 0; i < TDim; ++i)
        {
            is_inside = is_inside && mKnotVectors[i].IsInside(xi[i]);
            if (!is_inside) { break; }
        }
        return is_inside;
    }

    /// Compare between two BSplines patches in terms of parametric information
    bool IsCompatible(const BaseType& rOtherFESpace) const final
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const ThisType& rOtherBSplinesFESpace = dynamic_cast<const ThisType&>(rOtherFESpace);

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Number(i) == rOtherBSplinesFESpace.Number(i)))
            {
                return false;
            }
            if (!(this->Order(i) == rOtherBSplinesFESpace.Order(i)))
            {
                return false;
            }
            if (!(this->KnotVector(i) == rOtherBSplinesFESpace.KnotVector(i)))
            {
                return false;
            }
        }

        return true;
    }

    /// Reset all the dof numbers for each grid function to -1.
    void ResetFunctionIndices() final
    {
        BaseType::mGlobalToLocal.clear();
        if (mFunctionsIds.size() != this->TotalNumber())
        {
            mFunctionsIds.resize(this->TotalNumber());
        }
        std::fill(mFunctionsIds.begin(), mFunctionsIds.end(), -1);
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    void ResetFunctionIndices(const std::vector<std::size_t>& func_indices) final
    {
        if (func_indices.size() != this->TotalNumber())
        {
            KRATOS_WATCH(this->TotalNumber())
            std::cout << "func_indices:";
            for (std::size_t i = 0; i < func_indices.size(); ++i)
            {
                std::cout << " " << func_indices[i];
            }
            std::cout << std::endl;
            KRATOS_ERROR << "The func_indices vector does not have the same size as total number of basis functions";
        }
        if (mFunctionsIds.size() != this->TotalNumber())
        {
            mFunctionsIds.resize(this->TotalNumber());
        }
        // std::copy(func_indices.begin(), func_indices.end(), mFunctionsIds.begin());
        for (std::size_t i = 0; i < func_indices.size(); ++i)
        {
            mFunctionsIds[i] = func_indices[i];
            BaseType::mGlobalToLocal[mFunctionsIds[i]] = i;
        }
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    std::size_t& Enumerate(std::size_t& start) final
    {
        BaseType::mGlobalToLocal.clear();
        for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
        {
            if (mFunctionsIds[i] == -1) { mFunctionsIds[i] = start++; }
            BaseType::mGlobalToLocal[mFunctionsIds[i]] = i;
        }

        return start;
    }

    /// Access the function indices (aka global ids)
    std::vector<std::size_t> FunctionIndices() const final
    {
        return mFunctionsIds;
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map) final
    {
        for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
        {
            std::map<std::size_t, std::size_t>::const_iterator it = indices_map.find(mFunctionsIds[i]);

            if (it == indices_map.end())
            {
                std::cout << "WARNING!!! the indices_map does not contain " << mFunctionsIds[i] << std::endl;
                continue;
            }

            mFunctionsIds[i] = it->second;
        }

        BaseType::mGlobalToLocal.clear();
        for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
        {
            BaseType::mGlobalToLocal[mFunctionsIds[i]] = i;
        }
    }

    /// Get the first equation_id in this space
    std::size_t GetFirstEquationId() const final
    {
        std::size_t first_id, size = mFunctionsIds.size();

        for (std::size_t i = 0; i < size; ++i)
        {
            if (mFunctionsIds[i] == -1)
            {
                return -1;
            }
            else
            {
                if (i == 0)
                {
                    first_id = mFunctionsIds[i];
                }
                else if (mFunctionsIds[i] < first_id)
                {
                    first_id = mFunctionsIds[i];
                }
            }
        }

        return first_id;
    }

    /// Get the last equation_id in this space
    std::size_t GetLastEquationId() const final
    {
        std::size_t last_id = -1, size = mFunctionsIds.size();
        bool hit = false;

        for (std::size_t i = 0; i < size; ++i)
        {
            if (mFunctionsIds[i] != -1)
            {
                if (!hit)
                {
                    hit = true;
                    last_id = mFunctionsIds[i];
                }
                else
                {
                    if (mFunctionsIds[i] > last_id)
                    {
                        last_id = mFunctionsIds[i];
                    }
                }
            }
        }

        return last_id;
    }

    /// Extract the index of the functions on the boundaries
    std::vector<std::size_t> ExtractBoundaryFunctionIndicesByFlag(int boundary_id) const final
    {
        std::set<std::size_t> bf_id_set;
        bool first = false;

        for (int iside = _BLEFT_; iside < _NUMBER_OF_BOUNDARY_SIDE; ++iside)
        {
            BoundarySide side = static_cast<BoundarySide>(iside);
            if ((BOUNDARY_FLAG(side) & boundary_id) == BOUNDARY_FLAG(side))
            {
                std::vector<std::size_t> func_indices = this->ExtractBoundaryFunctionIndices(side);
                if (func_indices.size() != 0)
                {
                    if (!first)
                    {
                        // initialize the new set
                        bf_id_set.insert(func_indices.begin(), func_indices.end());
                        first = true;
                    }
                    else
                    {
                        // take the intersection
                        std::vector<std::size_t> temp;
                        std::set_intersection(bf_id_set.begin(), bf_id_set.end(), func_indices.begin(), func_indices.end(), std::back_inserter(temp));
                        bf_id_set.clear();
                        bf_id_set.insert(temp.begin(), temp.end());
                    }
                }
            }
        }

        return std::vector<std::size_t>(bf_id_set.begin(), bf_id_set.end());
    }

    /// Extract the index of the functions on the boundary
    std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const final
    {
        std::vector<std::size_t> size_info;
        return this->ExtractBoundaryFunctionIndices(size_info, side);
    }

    std::vector<std::size_t> ExtractBoundaryFunctionIndices(std::vector<std::size_t>& size_info, const BoundarySide& side) const final
    {
        std::vector<std::size_t> func_indices;

        if constexpr (TDim == 1)
        {
            size_info.resize(0);
            func_indices.resize(1);
            if (side == _BLEFT_)
            {
                func_indices[0] = mFunctionsIds[BSplinesIndexingUtility_Helper::Index1D(1, this->Number(0))];
            }
            else if (side  == _BRIGHT_)
            {
                func_indices[0] = mFunctionsIds[BSplinesIndexingUtility_Helper::Index1D(this->Number(0), this->Number(0))];
            }
        }
        else if constexpr (TDim == 2)
        {
            size_info.resize(1);
            if (side == _BLEFT_)
            {
                size_info[0] = this->Number(1);
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(1, j + 1, this->Number(0), this->Number(1))];
            }
            else if (side == _BRIGHT_)
            {
                size_info[0] = this->Number(1);
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(this->Number(0), j + 1, this->Number(0), this->Number(1))];
            }
            else if (side == _BBOTTOM_)
            {
                size_info[0] = this->Number(0);
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(i + 1, 1, this->Number(0), this->Number(1))];
            }
            else if (side == _BTOP_)
            {
                size_info[0] = this->Number(0);
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(i + 1, this->Number(1), this->Number(0), this->Number(1))];
            }
        }
        else if constexpr (TDim == 3)
        {
            size_info.resize(2);
            if (side == _BLEFT_)
            {
                size_info[0] = this->Number(1);
                size_info[1] = this->Number(2);
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(1, j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
            else if (side == _BRIGHT_)
            {
                size_info[0] = this->Number(1);
                size_info[1] = this->Number(2);
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(this->Number(0), j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
            else if (side == _BBOTTOM_)
            {
                size_info[0] = this->Number(0);
                size_info[1] = this->Number(1);
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, 1, this->Number(0), this->Number(1), this->Number(2))];
            }
            else if (side == _BTOP_)
            {
                size_info[0] = this->Number(0);
                size_info[1] = this->Number(1);
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, this->Number(2), this->Number(0), this->Number(1), this->Number(2))];
            }
            else if (side == _BFRONT_)
            {
                size_info[0] = this->Number(0);
                size_info[1] = this->Number(2);
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
            else if (side == _BBACK_)
            {
                size_info[0] = this->Number(0);
                size_info[1] = this->Number(2);
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, this->Number(1), k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }

        return func_indices;
    }

    /// Extract the index of the functions on the boundary down to some level
    std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side, std::size_t level) const final
    {
        std::vector<std::size_t> func_indices;

        if (side == _BLEFT_)
        {
            if constexpr (TDim == 1)
            {
                func_indices.resize(1);
                func_indices[0] = mFunctionsIds[BSplinesIndexingUtility_Helper::Index1D(1 + level, this->Number(0))];
            }
            else if constexpr (TDim == 2)
            {
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(1 + level, j + 1, this->Number(0), this->Number(1))];
            }
            else if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(1 + level, j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BRIGHT_)
        {
            if constexpr (TDim == 1)
            {
                func_indices.resize(1);
                func_indices[0] = mFunctionsIds[BSplinesIndexingUtility_Helper::Index1D(this->Number(0) - level, this->Number(0))];
            }
            else if constexpr (TDim == 2)
            {
                func_indices.resize(this->Number(1));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(this->Number(0) - level, j + 1, this->Number(0), this->Number(1))];
            }
            else if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(1)*this->Number(2));
                for (std::size_t j = 0; j < this->Number(1); ++j)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(this->Number(0) - level, j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BBOTTOM_)
        {
            if constexpr (TDim == 2)
            {
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(i + 1, 1 + level, this->Number(0), this->Number(1))];
            }
            else if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, 1 + level, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BTOP_)
        {
            if constexpr (TDim == 2)
            {
                func_indices.resize(this->Number(0));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))]
                        = mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(i + 1, this->Number(1) - level, this->Number(0), this->Number(1))];
            }
            else if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(1));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, this->Number(2) - level, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BFRONT_)
        {
            if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, 1 + level, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }
        else if (side == _BBACK_)
        {
            if constexpr (TDim == 3)
            {
                func_indices.resize(this->Number(0)*this->Number(2));
                for (std::size_t i = 0; i < this->Number(0); ++i)
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                        func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))]
                            = mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, this->Number(1) - level, k + 1, this->Number(0), this->Number(1), this->Number(2))];
            }
        }

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices, const bool override) final
    {
        if (side == _BLEFT_)
        {
            if constexpr (TDim == 1)
            {
                if (func_indices[0] != -1)
                {
                    const auto index = BSplinesIndexingUtility_Helper::Index1D(1, this->Number(0));
                    if (override)
                    {
                        const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                        mFunctionsIds[index] = func_indices[0];
                        BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                    }
                    else
                    {
                        if (mFunctionsIds[index] == -1)
                            mFunctionsIds[index] = func_indices[0];
                    }
                }
            }
            else if constexpr (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))];
                    if (aux != -1)
                    {
                        const auto index = BSplinesIndexingUtility_Helper::Index2D(1, j + 1, this->Number(0), this->Number(1));
                        if (override)
                        {
                            const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                            mFunctionsIds[index] = aux;
                            BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                        }
                        else
                        {
                            if (mFunctionsIds[index] == -1)
                                mFunctionsIds[index] = aux;
                        }
                    }
                }
            }
            else if constexpr (TDim == 3)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(1, j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
        else if (side == _BRIGHT_)
        {
            if constexpr (TDim == 1)
            {
                if (func_indices[0] != -1)
                {
                    const auto index = BSplinesIndexingUtility_Helper::Index1D(this->Number(0), this->Number(0));
                    if (override)
                    {
                        const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                        mFunctionsIds[index] = func_indices[0];
                        BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                    }
                    else
                    {
                        if (mFunctionsIds[index] == -1)
                            mFunctionsIds[index] = func_indices[0];
                    }
                }
            }
            else if constexpr (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index1D(j + 1, this->Number(1))];
                    if (aux != -1)
                    {
                        const auto index = BSplinesIndexingUtility_Helper::Index2D(this->Number(0), j + 1, this->Number(0), this->Number(1));
                        if (override)
                        {
                            const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                            mFunctionsIds[index] = aux;
                            BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                        }
                        else
                        {
                            if (mFunctionsIds[index] == -1)
                                mFunctionsIds[index] = aux;
                        }
                    }
                }
            }
            else if constexpr (TDim == 3)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(j + 1, k + 1, this->Number(1), this->Number(2))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(this->Number(0), j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
        else if (side == _BBOTTOM_)
        {
            if constexpr (TDim == 2)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))];
                    if (aux != -1)
                    {
                        const auto index = BSplinesIndexingUtility_Helper::Index2D(i + 1, 1, this->Number(0), this->Number(1));
                        if (override)
                        {
                            const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                            mFunctionsIds[index] = aux;
                            BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                        }
                        else
                        {
                            if (mFunctionsIds[index] == -1)
                                mFunctionsIds[index] = aux;
                        }
                    }
                }
            }
            else if constexpr (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, 1, this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
        else if (side == _BTOP_)
        {
            if constexpr (TDim == 2)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index1D(i + 1, this->Number(0))];
                    if (aux != -1)
                    {
                        const auto index = BSplinesIndexingUtility_Helper::Index2D(i + 1, this->Number(1), this->Number(0), this->Number(1));
                        if (override)
                        {
                            const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                            mFunctionsIds[index] = aux;
                            BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                        }
                        else
                        {
                            if (mFunctionsIds[index] == -1)
                                mFunctionsIds[index] = aux;
                        }
                    }
                }
            }
            else if constexpr (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, this->Number(2), this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
        else if (side == _BFRONT_)
        {
            if constexpr (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(i + 1, 1, k + 1, this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
        else if (side == _BBACK_)
        {
            if constexpr (TDim == 3)
            {
                for (std::size_t i = 0; i < this->Number(0); ++i)
                {
                    for (std::size_t k = 0; k < this->Number(2); ++k)
                    {
                        std::size_t aux = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, k + 1, this->Number(0), this->Number(2))];
                        if (aux != -1)
                        {
                            const auto index = BSplinesIndexingUtility_Helper::Index3D(i + 1, this->Number(1), k + 1, this->Number(0), this->Number(1), this->Number(2));
                            if (override)
                            {
                                const auto local_id = this->LocalId(mFunctionsIds[index]); // record the local id
                                mFunctionsIds[index] = aux;
                                BaseType::mGlobalToLocal[mFunctionsIds[index]] = local_id; // reassign the local id
                            }
                            else
                            {
                                if (mFunctionsIds[index] == -1)
                                    mFunctionsIds[index] = aux;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Construct the boundary patch based on side
    typename FESpace<TDim-1, TLocalCoordinateType>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const final
    {
        typename BSplinesFESpace<TDim-1, TLocalCoordinateType>::Pointer pBFESpace = typename BSplinesFESpace<TDim-1, TLocalCoordinateType>::Pointer(new BSplinesFESpace<TDim-1, TLocalCoordinateType> ());

        // assign the knot vectors
        if constexpr (TDim == 2)
        {
            if (side == _BLEFT_ || side == _BRIGHT_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(1).Clone());
                pBFESpace->SetInfo(0, Number(1), Order(1));
            }
            else if (side == _BTOP_ || side == _BBOTTOM_)
            {
                pBFESpace->SetKnotVector(0, KnotVector(0).Clone());
                pBFESpace->SetInfo(0, Number(0), Order(0));
            }
        }
        else if constexpr (TDim == 3)
        {
            if (side == _BLEFT_ || side == _BRIGHT_)
            {
                // on the left and right side, the boundary patch is always (v,w) -> (u,v)
                pBFESpace->SetKnotVector(0, KnotVector(1).Clone());
                pBFESpace->SetKnotVector(1, KnotVector(2).Clone());
                pBFESpace->SetInfo(0, Number(1), Order(1));
                pBFESpace->SetInfo(1, Number(2), Order(2));
            }
            else if (side == _BTOP_ || side == _BBOTTOM_)
            {
                // on the top and bottom side, the boundary patch is always (u,v) -> (u,v)
                pBFESpace->SetKnotVector(0, KnotVector(0).Clone());
                pBFESpace->SetKnotVector(1, KnotVector(1).Clone());
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(1), Order(1));
            }
            else if (side == _BFRONT_ || side == _BBACK_)
            {
                // on the front and back side, the boundary patch is always (u,w) -> (u,v)
                pBFESpace->SetKnotVector(0, KnotVector(0).Clone());
                pBFESpace->SetKnotVector(1, KnotVector(2).Clone());
                pBFESpace->SetInfo(0, Number(0), Order(0));
                pBFESpace->SetInfo(1, Number(2), Order(2));
            }
        }

        // transfer the function indices
        std::vector<std::size_t> b_func_indices = this->ExtractBoundaryFunctionIndices(side);
        pBFESpace->ResetFunctionIndices(b_func_indices);

        return pBFESpace;
    }

    /// Construct the boundary patch based on side and direction
    typename FESpace<TDim-1, TLocalCoordinateType>::Pointer ConstructBoundaryFESpace(const BoundarySide& side,
            const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const final
    {
        typename BSplinesFESpace<TDim-1, TLocalCoordinateType>::Pointer pBFESpace = typename BSplinesFESpace<TDim-1, TLocalCoordinateType>::Pointer(new BSplinesFESpace<TDim-1, TLocalCoordinateType> ());
        std::vector<int> param_dirs = ParameterDirection<TDim>::Get(side);

        // assign the knot vectors
        if constexpr (TDim == 2)
        {
            pBFESpace->SetKnotVector(0, KnotVector(param_dirs[0]).Clone(directions[0]));
            pBFESpace->SetInfo(0, Number(param_dirs[0]), Order(param_dirs[0]));
        }
        else if constexpr (TDim == 3)
        {
            std::vector<int> map_param_dirs(2);

            map_param_dirs[0] = param_dirs[local_parameter_map.at(0)];
            map_param_dirs[1] = param_dirs[local_parameter_map.at(1)];

            pBFESpace->SetKnotVector(0, KnotVector(map_param_dirs[0]).Clone(directions[0]));
            pBFESpace->SetKnotVector(1, KnotVector(map_param_dirs[1]).Clone(directions[1]));
            pBFESpace->SetInfo(0, Number(map_param_dirs[0]), Order(map_param_dirs[0]));
            pBFESpace->SetInfo(1, Number(map_param_dirs[1]), Order(map_param_dirs[1]));
        }

        // transfer the function indices
        std::vector<std::size_t> size_info;
        std::vector<std::size_t> b_func_indices = this->ExtractBoundaryFunctionIndices(size_info, side);

        // here transform the function indices
        if constexpr (TDim == 2)
        {
            if (directions[0] == BoundaryDirection::_REVERSED_)
            {
                // reverse the function indices
                BSplinesIndexingUtility::Reverse<1, std::vector<std::size_t>, std::vector<std::size_t> >(b_func_indices, size_info, 0);
            }
        }
        else if constexpr (TDim == 3)
        {
            const bool uv_or_vu = (local_parameter_map.at(0) == 0);
            BSplinesIndexingUtility::Transform(b_func_indices, size_info, uv_or_vu, directions[local_parameter_map.at(0)], directions[local_parameter_map.at(1)]);
        }

        pBFESpace->ResetFunctionIndices(b_func_indices);

        return pBFESpace;
    }

    /// Construct the sliced FESpace
    typename FESpace<TDim-1>::Pointer ConstructSlicedFESpace(int idir, TLocalCoordinateType xi) const final
    {
        return this->ConstructSlicedFESpace(idir);
    }

    /// Construct the sliced FESpace
    typename FESpace<TDim-1>::Pointer ConstructSlicedFESpace(int idir) const
    {
        typename FESpace<TDim-1>::Pointer pSFESpace;

        if (idir == 0)
        {
            pSFESpace = this->ConstructBoundaryFESpace(_BLEFT_);
        }
        else if (idir == 1)
        {
            pSFESpace = this->ConstructBoundaryFESpace(_BBOTTOM_);
        }
        else if (idir == 2)
        {
            pSFESpace = this->ConstructBoundaryFESpace(_BBACK_);
        }
        else
        {
            KRATOS_ERROR << "Invalid direction " << idir;
        }

        pSFESpace->ResetFunctionIndices(); // the sliced FESpace should contain no function indices
        return pSFESpace;
    }

    /// Construct the FESpace in single direction
    typename FESpace<1, TLocalCoordinateType>::Pointer ConstructUniaxialFESpace(int idir) const
    {
        typename BSplinesFESpace<1, TLocalCoordinateType>::Pointer pBFESpace = typename BSplinesFESpace<1, TLocalCoordinateType>::Pointer(new BSplinesFESpace<1, TLocalCoordinateType>());
        pBFESpace->SetKnotVector(0, KnotVector(idir).Clone());
        pBFESpace->SetInfo(0, Number(idir), Order(idir));
        pBFESpace->ResetFunctionIndices();
        return pBFESpace;
    }

    /// Create the cell manager for all the cells in the support domain of the BSplinesFESpace
    typename BaseType::cell_container_t::Pointer ConstructCellManager() const final
    {
        typename cell_container_t::Pointer pCellManager;

        std::vector<std::size_t> func_indices = this->FunctionIndices();

        pCellManager = typename cell_container_t::Pointer(new BCellManager<TDim, BCell>());

        if constexpr (TDim == 1)
        {
            // firstly compute the Bezier extraction operator on the patch
            std::vector<Matrix> C;
            int ne1;
            // BezierUtils::bezier_extraction_2d(C, ne1, ne2, this->KnotVector(0), this->KnotVector(1), this->Order(0), this->Order(1));
            BezierUtils::bezier_extraction_1d(C, ne1, this->KnotVector(0), this->Order(0));

#ifdef DEBUG_GEN_CELL
            KRATOS_WATCH(ne1)
            KRATOS_WATCH(C.size())
            KRATOS_WATCH(C[0].size1())
            KRATOS_WATCH(C[0].size2())
#endif

            // construct cells and add to the manager
            std::size_t n1 = this->Number(0);
            std::size_t p1 = this->Order(0);
            std::size_t i, j, k, l, b1 = p1 + 1, tmp, mul1, sum_mul1 = 0, id1, id;
            std::size_t cnt = 0; // cell counter

            for (i = 0; i < ne1; ++i)
            {
                // check the multiplicity
                tmp = b1;
                while (b1 <= (n1 + p1 + 1) && this->KnotVector(0)[b1] == this->KnotVector(0)[b1 - 1]) { ++b1; }
                mul1 = b1 - tmp + 1;
                b1 = b1 + 1;
                sum_mul1 = sum_mul1 + (mul1 - 1);

                std::vector<std::size_t> anchors;
                anchors.reserve(p1 + 1);
                for (k = 0; k < p1 + 1; ++k)
                {
                    id1 = i + k + sum_mul1;
                    id = id1; // this is the local id
                    anchors.push_back(id);
                }

                // add the cell
                std::tuple<knot_t, knot_t> span1 = this->KnotVector(0).span(i + 1);
                BCell::Pointer p_cell = BCell::Pointer(new BCell(cnt, std::get<0>(span1), std::get<1>(span1)));
                double W = 1.0; // here we set to one because B-Splines space does not have weight
                for (std::size_t r = 0; r < (p1 + 1); ++r)
                {
                    p_cell->AddAnchor(func_indices[anchors[r]], W, row(C[cnt], r));
                }
                pCellManager->insert(p_cell);
                ++cnt;
            }
        }
        else if constexpr (TDim == 2)
        {
            // firstly compute the Bezier extraction operator on the patch
            std::vector<Matrix> C;
            int ne1, ne2;
            // BezierUtils::bezier_extraction_2d(C, ne1, ne2, this->KnotVector(0), this->KnotVector(1), this->Order(0), this->Order(1));
            BezierUtils::bezier_extraction_2d(C, ne2, ne1, this->KnotVector(1), this->KnotVector(0), this->Order(1), this->Order(0)); // we rotate the order of input

#ifdef DEBUG_GEN_CELL
            KRATOS_WATCH(ne1)
            KRATOS_WATCH(ne2)
            KRATOS_WATCH(C.size())
            KRATOS_WATCH(C[0].size1())
            KRATOS_WATCH(C[0].size2())
#endif

            // construct cells and add to the manager
            std::size_t n1 = this->Number(0);
            std::size_t n2 = this->Number(1);
            std::size_t p1 = this->Order(0);
            std::size_t p2 = this->Order(1);
            std::size_t i, j, k, l, b1 = p1 + 1, b2, tmp, mul1, mul2, sum_mul1 = 0, sum_mul2 = 0, id1, id2, id;
            std::size_t cnt = 0; // cell counter

            for (i = 0; i < ne1; ++i)
            {
                // check the multiplicity
                tmp = b1;
                while (b1 <= (n1 + p1 + 1) && this->KnotVector(0)[b1] == this->KnotVector(0)[b1 - 1]) { ++b1; }
                mul1 = b1 - tmp + 1;
                b1 = b1 + 1;
                sum_mul1 = sum_mul1 + (mul1 - 1);

                b2 = p2 + 1;
                sum_mul2 = 0;
                for (j = 0; j < ne2; ++j)
                {
                    // check the multiplicity
                    tmp = b2;
                    while (b2 <= (n2 + p2 + 1) && this->KnotVector(1)[b2] == this->KnotVector(1)[b2 - 1]) { ++b2; }
                    mul2 = b2 - tmp + 1;
                    b2 = b2 + 1;
                    sum_mul2 = sum_mul2 + (mul2 - 1);

                    std::vector<std::size_t> anchors;
                    anchors.reserve((p1 + 1) * (p2 + 1));
                    for (k = 0; k < p1 + 1; ++k)
                    {
                        for (l = 0; l < p2 + 1; ++l)
                        {
                            id1 = i + k + sum_mul1;
                            id2 = j + l + sum_mul2;
                            id = id1 + id2 * n1; // this is the local id
                            anchors.push_back(id);
                        }
                    }

                    // add the cell
                    std::tuple<knot_t, knot_t> span1 = this->KnotVector(0).span(i + 1);
                    std::tuple<knot_t, knot_t> span2 = this->KnotVector(1).span(j + 1);
                    BCell::Pointer p_cell = BCell::Pointer(new BCell(cnt, std::get<0>(span1), std::get<1>(span1), std::get<0>(span2), std::get<1>(span2)));
                    double W = 1.0; // here we set to one because B-Splines space does not have weight
                    for (std::size_t r = 0; r < (p1 + 1) * (p2 + 1); ++r)
                    {
                        p_cell->AddAnchor(func_indices[anchors[r]], W, row(C[cnt], r));
                    }
                    pCellManager->insert(p_cell);
                    ++cnt;
                }
            }
        }
        else if constexpr (TDim == 3)
        {
            // firstly compute the Bezier extraction operator on the patch
            std::vector<Matrix> C;
            int ne1, ne2, ne3;
            // BezierUtils::bezier_extraction_3d(C, ne1, ne2, ne3,
            //     this->KnotVector(0), this->KnotVector(1), this->KnotVector(2),
            //     this->Order(0), this->Order(1), this->Order(2));
            BezierUtils::bezier_extraction_3d(C, ne3, ne2, ne1,
                                              this->KnotVector(2), this->KnotVector(1), this->KnotVector(0),
                                              this->Order(2), this->Order(1), this->Order(0)); // we rotate the order of input

            // construct cells and add to the manager
            std::size_t n1 = this->Number(0);
            std::size_t n2 = this->Number(1);
            std::size_t n3 = this->Number(2);
            std::size_t p1 = this->Order(0);
            std::size_t p2 = this->Order(1);
            std::size_t p3 = this->Order(2);
            std::size_t i, j, k, l, u, v, w, b1 = p1 + 1, b2, b3, tmp, mul1, mul2, mul3;
            std::size_t sum_mul1 = 0, sum_mul2 = 0, sum_mul3 = 0, id1, id2, id3, id;
            std::size_t cnt = 0; // cell counter

            for (i = 0; i < ne1; ++i)
            {
                // check the multiplicity
                tmp = b1;
                while (b1 <= (n1 + p1 + 1) && this->KnotVector(0)[b1] == this->KnotVector(0)[b1 - 1]) { ++b1; }
                mul1 = b1 - tmp + 1;
                b1 = b1 + 1;
                sum_mul1 = sum_mul1 + (mul1 - 1);

                b2 = p2 + 1;
                sum_mul2 = 0;
                for (j = 0; j < ne2; ++j)
                {
                    // check the multiplicity
                    tmp = b2;
                    while (b2 <= (n2 + p2 + 1) && this->KnotVector(1)[b2] == this->KnotVector(1)[b2 - 1]) { ++b2; }
                    mul2 = b2 - tmp + 1;
                    b2 = b2 + 1;
                    sum_mul2 = sum_mul2 + (mul2 - 1);

                    b3 = p3 + 1;
                    sum_mul3 = 0;
                    for (k = 0; k < ne3; ++k)
                    {
                        // check the multiplicity
                        tmp = b3;
                        while (b3 <= (n3 + p3 + 1) && this->KnotVector(2)[b3] == this->KnotVector(2)[b3 - 1]) { ++b3; }
                        mul3 = b3 - tmp + 1;
                        b3 = b3 + 1;
                        sum_mul3 = sum_mul3 + (mul3 - 1);

                        std::vector<std::size_t> anchors;
                        anchors.reserve((p1 + 1) * (p2 + 1) * (p3 + 1));
                        for (u = 0; u < p1 + 1; ++u)
                        {
                            for (v = 0; v < p2 + 1; ++v)
                            {
                                for (v = 0; v < p2 + 1; ++v)
                                {
                                    for (w = 0; w < p3 + 1; ++w)
                                    {
                                        id1 = i + u + sum_mul1;
                                        id2 = j + v + sum_mul2;
                                        id3 = k + w + sum_mul3;
                                        id = id1 + (id2 + id3 * n2) * n1; // this is the local id
                                        anchors.push_back(id);
                                    }
                                }
                            }
                        }

                        // add the cell
                        std::tuple<knot_t, knot_t> span1 = this->KnotVector(0).span(i + 1);
                        std::tuple<knot_t, knot_t> span2 = this->KnotVector(1).span(j + 1);
                        std::tuple<knot_t, knot_t> span3 = this->KnotVector(2).span(k + 1);
                        BCell::Pointer p_cell = BCell::Pointer(new BCell(cnt, std::get<0>(span1), std::get<1>(span1),
                                                               std::get<0>(span2), std::get<1>(span2), std::get<0>(span3), std::get<1>(span3)));
                        double W = 1.0; // here we set to one because B-Splines space does not have weight
                        for (std::size_t r = 0; r < (p1 + 1) * (p2 + 1) * (p3 + 1); ++r)
                        {
                            p_cell->AddAnchor(func_indices[anchors[r]], W, row(C[cnt], r));
                        }
                        pCellManager->insert(p_cell);
                        ++cnt;
                    }
                }
            }
        }

        return pCellManager;
    }

    /// Overload assignment operator
    ThisType& operator=(const ThisType& rOther)
    {
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            this->SetKnotVector(dim, rOther.KnotVector(dim));
            this->SetInfo(dim, rOther.Number(dim), rOther.Order(dim));
        }
        this->mFunctionsIds = rOther.mFunctionsIds;
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    typename BaseType::Pointer Clone() const final
    {
        typename ThisType::Pointer pNewFESpace = typename ThisType::Pointer(new ThisType());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Type() << ", Addr = " << this << ", n = (";
        for (std::size_t i = 0; i < TDim; ++i)
        {
            rOStream << " " << this->Number(i);
        }
        rOStream << "), p = (";
        for (std::size_t i = 0; i < TDim; ++i)
        {
            rOStream << " " << this->Order(i);
        }
        rOStream << ")";
    }

    void PrintData(std::ostream& rOStream) const final
    {
        for (std::size_t i = 0; i < TDim; ++i)
        {
            rOStream << " knot vector " << i << ":";
            for (std::size_t j = 0; j < mKnotVectors[i].size(); ++j)
            {
                rOStream << " " << mKnotVectors[i].pKnotAt(j)->Value();
            }
            rOStream << std::endl;
        }
        if (mFunctionsIds.size() == this->TotalNumber())
        {
            rOStream << " Function Indices:";
            if constexpr (TDim == 1)
            {
                for (std::size_t i = 0; i < mFunctionsIds.size(); ++i)
                {
                    rOStream << " " << mFunctionsIds[i];
                }
            }
            else if constexpr (TDim == 2)
            {
                for (std::size_t j = 0; j < this->Number(1); ++j)
                {
                    for (std::size_t i = 0; i < this->Number(0); ++i)
                    {
                        rOStream << " " << mFunctionsIds[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, this->Number(0), this->Number(1))];
                    }
                    rOStream << std::endl;
                }
            }
            else if constexpr (TDim == 3)
            {
                for (std::size_t k = 0; k < this->Number(2); ++k)
                {
                    for (std::size_t j = 0; j < this->Number(1); ++j)
                    {
                        for (std::size_t i = 0; i < this->Number(0); ++i)
                        {
                            rOStream << " " << mFunctionsIds[BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, this->Number(0), this->Number(1), this->Number(2))];
                        }
                        rOStream << std::endl;
                    }
                    rOStream << std::endl;
                }
            }
        }
    }

private:

    /**
     * internal data to construct the shape functions on the BSplines
     */
    boost::array<std::size_t, TDim> mOrders;
    boost::array<std::size_t, TDim> mNumbers;
    boost::array<knot_container_t, TDim> mKnotVectors;

    /**
     * data for grid function interpolation
     */
    std::vector<std::size_t> mFunctionsIds; // this is to store a unique number of the shape function over the forest of FESpace(s).
};

/**
 * Template specific instantiation for null-D BSplines patch to terminate the compilation
 */
template<typename TLocalCoordinateType>
class BSplinesFESpace<0, TLocalCoordinateType> : public FESpace<0, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesFESpace);

    /// Type definition
    typedef FESpace<0, TLocalCoordinateType> BaseType;
    typedef KnotArray1D<TLocalCoordinateType> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    BSplinesFESpace() : BaseType() {}

    /// Destructor
    ~BSplinesFESpace() override {}

    /// Get the order of the BSplines patch in specific direction
    std::size_t Order(std::size_t i) const final {return 0;}

    /// Get the number of basis functions defined over the BSplinesFESpace on one direction
    std::size_t Number(std::size_t i) const {return 1;}

    /// Get the number of basis functions defined over the BSplinesFESpace
    std::size_t Number() const {return 1;}

    /// Get the string describing the type of the patch
    std::string Type() const final
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "BSplinesFESpace0D";
    }

    /// Set the knot vector in direction i.
    void SetKnotVector(std::size_t i, const knot_container_t& p_knot_vector)
    {}

    /// Get the knot vector in i-direction
    const knot_container_t KnotVector(std::size_t i) const
    {
        return knot_container_t();
    }

    std::vector<std::size_t> Numbers() const
    {
        return std::vector<std::size_t>();
    }

    /// Set the BSplines information in the direction i
    void SetInfo(std::size_t i, std::size_t Number, std::size_t Order)
    {}

    /// Validate the BSplinesFESpace
    bool Validate() const final
    {
        return BaseType::Validate();
    }

    /// Compare between two BSplines patches in terms of parametric information
    bool IsCompatible(const BaseType& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other FESpace type is not " << Type() << std::endl;
            return false;
        }

        return true;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BSplinesFESpace<0>";
    }
};

/// output stream function
template<int TDim, typename TLocalCoordinateType>
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesFESpace<TDim, TLocalCoordinateType>& rThis)
{
    rOStream << "-------------Begin BSplinesFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End BSplinesFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_GEN_CELL

#include "bsplines_fespace.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_H_INCLUDED defined
