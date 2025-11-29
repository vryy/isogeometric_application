//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/bsplines_patch_interface.h"
// #include "custom_utilities/hierarchical_bsplines/hb_mesh.h"

namespace Kratos
{

template<int TDim, typename TDataType>
struct ComputeBsplinesDegreeElevation_Helper
{
    static void Compute(const StructuredControlGrid<TDim, TDataType>& ControlValues,
                        const BSplinesFESpace<TDim>& rFESpace,
                        const std::vector<std::size_t>& order_increment,
                        StructuredControlGrid<TDim, TDataType>& NewControlValues,
                        std::vector<std::vector<double> >& new_knots)
    {
        KRATOS_ERROR << "Not implemented for dimension " << TDim;
    }
};

template<int TDim>
struct ComputeBsplinesKnotInsertionCoefficients_Helper
{
    static void Compute(Matrix& T,
                        std::vector<std::vector<double> >& new_knots,
                        typename BSplinesFESpace<TDim>::Pointer& pFESpace,
                        const std::vector<std::vector<double> >& ins_knots)
    {
        KRATOS_ERROR << "Not implemented for dimension " << TDim;
    }
};

/**
 * Utility to perform the refinement on multipatch structure
 */
class MultiPatchRefinementUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchRefinementUtility);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /*************************************************************************
                                    B-SPLINES
    *************************************************************************/

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    template<class TPatchType>
    static void InsertKnots(typename TPatchType::Pointer& pPatch,
                     const std::vector<std::vector<typename TPatchType::LocalCoordinateType> >& ins_knots)
    {
        std::map<std::size_t, std::vector<int> > refined_patches;
        std::map<std::size_t, Matrix> trans_mats;
        bool record_trans_mat = false;
        InsertKnots<TPatchType>(pPatch, refined_patches, ins_knots, trans_mats, record_trans_mat);
    }

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    /// The transformation matrix will be stored. It will be useful for geometric multigrid.
    template<class TPatchType>
    static void InsertKnots(typename TPatchType::Pointer& pPatch,
                     const std::vector<std::vector<typename TPatchType::LocalCoordinateType> >& ins_knots,
                     std::map<std::size_t, Matrix>& trans_mats)
    {
        std::map<std::size_t, std::vector<int> > refined_patches;
        bool record_trans_mat = true;
        InsertKnots<TPatchType>(pPatch, refined_patches, ins_knots, trans_mats, record_trans_mat);
    }

    /// Insert the knots to the NURBS patch and make it compatible across neighbors
    /// if record_trans_mat is true, the transformation matrix for each patch will be stored in trans_mats
    template<class TPatchType>
    static void InsertKnots(typename TPatchType::Pointer& pPatch,
                     std::map<std::size_t, std::vector<int> >& refined_patches,
                     const std::vector<std::vector<typename TPatchType::LocalCoordinateType> >& ins_knots,
                     std::map<std::size_t, Matrix>& trans_mats,
                     bool record_trans_mat = false);

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<class TPatchType>
    static void DegreeElevate(typename TPatchType::Pointer& pPatch,
                       const std::vector<std::size_t>& order_increment)
    {
        std::map<std::size_t, std::vector<int> > refined_patches;
        DegreeElevate<TPatchType>(pPatch, refined_patches, order_increment);
    }

    /// Degree elevation for the NURBS patch and make it compatible across neighbors
    template<class TPatchType>
    static void DegreeElevate(typename TPatchType::Pointer& pPatch,
                       std::map<std::size_t, std::vector<int> >& refined_patches,
                       const std::vector<std::size_t>& order_increment);

    /*************************************************************************
                              HIERARCHICAL B-SPLINES
                                    (TODO)
    *************************************************************************/

private:

    /// Compute the transformation matrix for knot insertion (NURBS version)
    template<int TDim>
    static void ComputeBsplinesKnotInsertionCoefficients(
        Matrix& T,
        std::vector<std::vector<double> >& new_knots,
        typename BSplinesFESpace<TDim>::Pointer& pFESpace,
        const std::vector<std::vector<double> >& ins_knots)
    {
        ComputeBsplinesKnotInsertionCoefficients_Helper<TDim>::Compute(T, new_knots, pFESpace, ins_knots);
    }

    template<int TDim, typename TDataType>
    static void ComputeBsplinesDegreeElevation(
        const StructuredControlGrid<TDim, TDataType>& ControlValues,
        const BSplinesFESpace<TDim>& rFESpace,
        const std::vector<std::size_t>& order_increment,
        StructuredControlGrid<TDim, TDataType>& NewControlValues,
        std::vector<std::vector<double> >& new_knots)
    {
        ComputeBsplinesDegreeElevation_Helper<TDim, TDataType>::Compute(ControlValues,
                rFESpace, order_increment, NewControlValues, new_knots);
    }

};

} // namespace Kratos.

#include "multipatch_refinement_utility.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_H_INCLUDED defined
