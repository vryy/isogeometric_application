//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/isogeometric_utility.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/multipatch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class MultiPatchUtility : public IsogeometricUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchUtility);

    /// Type definition

    /// Default constructor
    MultiPatchUtility() {}

    /// Destructor
    ~MultiPatchUtility() override {}

    /// Create new patch from a FESpace and wrap it with pointer
#ifdef SD_APP_FORWARD_COMPATIBILITY
    template<class TPatchType>
    static iga::Wrapper<TPatchType, typename TPatchType::Pointer> CreatePatchPointer(std::size_t Id, typename TPatchType::FESpaceType::Pointer pFESpace)
    {
        return iga::Wrapper<TPatchType, typename TPatchType::Pointer>(typename TPatchType::Pointer(new TPatchType(Id, pFESpace)));
    }
#else
    template<class TPatchType>
    static typename TPatchType::Pointer CreatePatchPointer(std::size_t Id, typename TPatchType::FESpaceType::Pointer pFESpace)
    {
        return typename TPatchType::Pointer(new TPatchType(Id, pFESpace));
    }
#endif

    /// Make a simple interface between two patches
    /// For BSplines patch, one shall use BSplinesPatchUtility::MakeInterfacexD instead
    template<int TDim>
    static void MakeInterface(typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
                              typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2)
    {
        typename PatchInterface<TDim>::Pointer pInterface12;
        typename PatchInterface<TDim>::Pointer pInterface21;

        pInterface12 = iga::make_shared<PatchInterface<TDim> >(pPatch1, side1, pPatch2, side2);
        pInterface21 = iga::make_shared<PatchInterface<TDim> >(pPatch2, side2, pPatch1, side1);

        pInterface12->SetOtherInterface(pInterface21);
        pInterface21->SetOtherInterface(pInterface12);

        pPatch1->AddInterface(pInterface12);
        pPatch2->AddInterface(pInterface21);
    }

    /// Check all the interfaces of the multipatch for compatibility issue
    template<int TDim>
    static void CheckInterfaces(const MultiPatch<TDim>& rMultiPatch, const bool debug = false, const double dist_tol = 0.0)
    {
        // loop through all the patches
        for (typename MultiPatch<TDim>::patch_const_iterator it = rMultiPatch.begin(); it != rMultiPatch.end(); ++it)
        {
            // loop through all interfaces of the patch
            for (auto it_interface = it->InterfaceBegin(); it_interface != it->InterfaceEnd(); ++it_interface)
            {
                // validate the interface
                const bool is_valid = (*it_interface)->Validate(debug, dist_tol);

                if (!is_valid)
                {
                    KRATOS_ERROR << "Interface between patch " << (*it_interface)->pPatch1()->Id()
                                 << " and patch " << (*it_interface)->pPatch2()->Id()
                                 << " is not valid.";
                }
            }
        }

        std::cout << "MultiPatch interfaces are all valid." << std::endl;
    }

    /// Construct the 12 edge patches of a 3D patch
    static std::vector<Patch<1>::Pointer> ConstructEdgePatches(Patch<3>::ConstPointer pPatch)
    {
        Patch<2>::Pointer pTopPatch = pPatch->ConstructBoundaryPatch(_BTOP_);
        Patch<2>::Pointer pBottomPatch = pPatch->ConstructBoundaryPatch(_BBOTTOM_);
        Patch<2>::Pointer pLeftPatch = pPatch->ConstructBoundaryPatch(_BLEFT_);
        Patch<2>::Pointer pRightPatch = pPatch->ConstructBoundaryPatch(_BRIGHT_);

        std::vector<Patch<1>::Pointer> pEdgePatches;

        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BBOTTOM_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BTOP_));

        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BBOTTOM_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BTOP_));

        pEdgePatches.push_back(pLeftPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pLeftPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pRightPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pRightPatch->ConstructBoundaryPatch(_BRIGHT_));

        return pEdgePatches;
    }

    /// Find the local coordinates of a point on a patch. The sampling is performed
    /// on the patch to determine the best initial point.
    /// On output return the error code (0: successful)
    template<int TDim>
    static int LocalCoordinates(const Patch<TDim>& rPatch,
                                const array_1d<double, 3>& point, array_1d<double, 3>& xi, const std::vector<int>& nsampling)
    {
        /// find the best initial point by sampling
        rPatch.Predict(point, xi, nsampling);

        /// compute the local coordinates
        return rPatch.LocalCoordinates(point, xi);
    }

    /// Find the local coordinates of a point on a multipatch. The sampling is performed
    /// on each patch to determine the best initial point.
    /// On output return the patch_id if sucessful; -1 otherwise
    template<int TDim>
    static int LocalCoordinates(const MultiPatch<TDim>& rMultiPatch,
                                const array_1d<double, 3>& point, array_1d<double, 3>& xi,
                                const std::vector<int>& nsampling,
                                const int echo_level = 0)
    {
        int error_code;

        if (echo_level > 0) std::cout << "searching local coordinates for point " << point << std::endl;

        for (typename MultiPatch<TDim>::patch_const_iterator it = rMultiPatch.begin(); it != rMultiPatch.end(); ++it)
        {
            error_code = LocalCoordinates( *it, point, xi, nsampling);
            if (echo_level > 1) std::cout << "error code for patch " << it->Id() << ": " << error_code << std::endl;
            if (error_code == 0)
            {
                if (echo_level > 0) std::cout << "searching local coordinates for point " << point << " successfully" << std::endl;
                return it->Id();
            }
        }

        return -1;
    }

    /// Compute the spatial derivatives of a scalar grid function. The values are
    /// a matrix that is organized such as
    ///     [dv/dx dv/dy dv/dz] in 3D
    /// and
    ///     [dv/dx dv/dy] in 2D
    template<int TDim>
    static Vector ComputeSpatialDerivatives(const GridFunction<TDim, array_1d<double, 3> >& rControlPointGridFunction,
                                            const GridFunction<TDim, double>& rControlValueGridFunction,
                                            const std::vector<double>& rCoordinates)
    {
        // compute the Jacobian
        std::vector<array_1d<double, 3> > tmp;
        rControlPointGridFunction.GetDerivative(tmp, rCoordinates);

        Matrix Jac(TDim, TDim);

        for (int i = 0; i < TDim; ++i)
        {
            for (int j = 0; j < TDim; ++j)
            {
                Jac(i, j) = tmp[j][i];
            }
        }

        // compute the spatial derivatives
        std::vector<double> tmp1;
        rControlValueGridFunction.GetDerivative(tmp1, rCoordinates);

        Vector dv(TDim);

        for (int i = 0; i < TDim; ++i)
        {
            dv(i) = tmp1[i];
        }

        Matrix InvJ(TDim, TDim);
        double DetJ;
        MathUtils<double>::InvertMatrix(Jac, InvJ, DetJ);

        Vector Result(TDim);
        noalias(Result) = prod(InvJ, dv);

        return Result;
    }

    /// Compute the spatial derivatives of a vector 3d grid function. The values are
    /// a matrix that is organized such as
    ///     [dvx/dx dvx/dy dvx/dz
    ///      dvy/dx dvy/dy dvy/dz
    ///      dvz/dx dvz/dy dvz/dz] in 3D
    /// and
    ///     [dvx/dx dvx/dy
    ///      dvy/dx dvy/dy
    ///      dvz/dx dvz/dy] in 2D
    template<int TDim>
    static Matrix ComputeSpatialDerivatives(const GridFunction<TDim, array_1d<double, 3> >& rControlPointGridFunction,
                                            const GridFunction<TDim, array_1d<double, 3> >& rControlValueGridFunction,
                                            const std::vector<double>& rCoordinates)
    {
        // compute the Jacobian
        std::vector<array_1d<double, 3> > tmp;
        rControlPointGridFunction.GetDerivative(tmp, rCoordinates);

        Matrix Jac(TDim, TDim);

        for (int i = 0; i < TDim; ++i)
        {
            for (int j = 0; j < TDim; ++j)
            {
                Jac(i, j) = tmp[j][i];
            }
        }

        // compute the spatial derivatives
        std::vector<array_1d<double, 3> > tmp1;
        rControlValueGridFunction.GetDerivative(tmp1, rCoordinates);

        Matrix Dv(3, TDim);

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < TDim; ++j)
            {
                Dv(i, j) = tmp1[j][i];
            }
        }

        Matrix InvJ(TDim, TDim);
        double DetJ;
        MathUtils<double>::InvertMatrix(Jac, InvJ, DetJ);

        Matrix Result(3, TDim);
        noalias(Result) = prod(Dv, InvJ);

        return Result;
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MultiPatchUtility";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

}; // end class MultiPatchUtility

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED defined
