//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Oct 2024 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <cstddef>
#include <iostream>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "isogeometric_projection_utility.h"

// #define CHECK_DERIVATIVES

namespace Kratos
{

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::PredictVerticalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, typename Patch<TDim>::Pointer pPatch,
        const std::array<unsigned int, TDim>& nsampling
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    double dist = 1.0e99;
    std::vector<double> xi(TDim);
    array_1d<double, 3> P;

    rLocalPoint.clear();

    if constexpr (TDim == 1)
    {
        for (unsigned int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double)i) / nsampling[0];

            pControlGridFunc->GetValue(P, xi);

            double d = sqrt(pow(rPoint[0] - P[0], 2) + pow(rPoint[1] - P[1], 2));

            if (d < dist)
            {
                dist = d;
                rLocalPoint[0] = xi[0];
            }
        }
    }
    else if constexpr (TDim == 2)
    {
        for (unsigned int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double)i) / nsampling[0];
            for (unsigned int j = 0; j < nsampling[1] + 1; ++j)
            {
                xi[1] = ((double)j) / nsampling[1];

                pControlGridFunc->GetValue(P, xi);

                double d = sqrt(pow(rPoint[0] - P[0], 2) + pow(rPoint[1] - P[1], 2));

                if (d < dist)
                {
                    dist = d;
                    rLocalPoint[0] = xi[0];
                    rLocalPoint[1] = xi[1];
                }
            }
        }
    }

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeVerticalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
        typename Patch<TDim>::Pointer pPatch,
        double TOL, int max_iters,
        int echo_level
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    std::vector<array_1d<double, 3> > dP;

    if (echo_level > 0)
    {
        std::cout << "Global point to be projected: " << rPoint << std::endl;
    }

    int it = 0;
    bool converged = false;

    if constexpr (TDim == 1)
    {
        // TODO
        KRATOS_ERROR << "To be implemented";
    }
    else if constexpr (TDim == 2)
    {
        double a00, a01, a10, a11;

        while (!converged && (it < max_iters))
        {
            pControlGridFunc->GetValue(rGlobalPoint, rLocalPoint);

            double dist = sqrt(pow(rPoint[0] - rGlobalPoint[0], 2) + pow(rPoint[1] - rGlobalPoint[1], 2));

            if (echo_level > 1)
            {
                std::cout << "At iteration " << (it + 1) << ", local point = " << rLocalPoint[0] << ", " << rLocalPoint[1]
                          << ", projected global point = " << rGlobalPoint
                          << ", distance = " << dist
                          << std::endl;
            }

            if (dist < TOL)
            {
                converged = true;
            }
            else
            {
                // pPatch->pFESpace()->PrintInfo(std::cout); std::cout << std::endl;
                // std::vector<double> f_values;
                // std::vector<std::vector<double> > f_derivatives;
                // pPatch->pFESpace()->GetValues(f_values, rLocalPoint);
                // pPatch->pFESpace()->GetDerivatives(f_derivatives, rLocalPoint);
                // // pPatch->pFESpace()->GetValuesAndDerivatives(f_values, f_derivatives, rLocalPoint);
                // std::cout << "f_values:";
                // for (std::size_t i = 0; i < f_values.size(); ++i)
                //     std::cout << " " << f_values[i];
                // std::cout << std::endl;
                // for (std::size_t i = 0; i < f_derivatives.size(); ++i)
                // {
                //     std::cout << "f_derivatives[" << i << "]:";
                //     for (std::size_t j = 0; j < f_derivatives[i].size(); ++j)
                //         std::cout << " " << f_derivatives[i][j];
                //     std::cout << std::endl;
                // }
                // std::cout << "control points:" << std::endl;
                // for (std::size_t i = 0; i < f_values.size(); ++i)
                //     std::cout << i << ": " << pControlGridFunc->pControlGrid()->GetData(i) << std::endl;

                pControlGridFunc->GetDerivative(dP, rLocalPoint);
                a00 = dP[0][0];
                a01 = dP[1][0];
                a10 = dP[0][1];
                a11 = dP[1][1];

                double det = a00 * a11 - a01 * a10;
                rLocalPoint[0] += ( a11 * (rPoint[0] - rGlobalPoint[0]) - a01 * (rPoint[1] - rGlobalPoint[1])) / det;
                rLocalPoint[1] += (-a10 * (rPoint[0] - rGlobalPoint[0]) + a00 * (rPoint[1] - rGlobalPoint[1])) / det;

                if (rLocalPoint[0] < 0.0) { rLocalPoint[0] = 0.0; }
                if (rLocalPoint[0] > 1.0) { rLocalPoint[0] = 1.0; }
                if (rLocalPoint[1] < 0.0) { rLocalPoint[1] = 0.0; }
                if (rLocalPoint[1] > 1.0) { rLocalPoint[1] = 1.0; }

                if (echo_level > 2)
                {
                    std::cout << "tangent matrix: [[" << a00 << ", " << a01 << "], [" << a10 << ", " << a11 << "]" << std::endl;
                }
            }

            ++it;
        }
    }

    if (it >= max_iters && !converged)
    {
        return 1;
    }

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeVerticalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<TDim>::Pointer pMultiPatch,
        double TOL, int max_iters,
        int echo_level
)
{
    typedef MultiPatch<TDim> MultiPatchType;
    typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

    std::vector<double> InitialLocalPoint = rLocalPoint;

    for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
    {
        // make a check to make sure the point can project into the bounding box
        std::vector<double> bounding_box;
        (*it)->GetBoundingBox(bounding_box);

        bool is_in = true;
        if constexpr (TDim == 1)
        {
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL)
                          && (rPoint[0] < bounding_box[1] + TOL);
        }
        else if constexpr (TDim == 2)
        {
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL)
                          && (rPoint[0] < bounding_box[1] + TOL)
                          && (rPoint[1] > bounding_box[2] - TOL)
                          && (rPoint[1] < bounding_box[3] + TOL);
        }

        if (!is_in)
        {
            continue;
        }

        // compute the vertical projection
        rLocalPoint = InitialLocalPoint; // reset the initial point
        int error_code = ComputeVerticalProjection(rPoint, rLocalPoint, rGlobalPoint, *it, TOL, max_iters, echo_level);

        if (error_code == 0)
        {
            patch_id = (*it)->Id();
            return 0;
        }
    }

    patch_id = -1;
    return 1; // can't find the projection point
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeVerticalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<TDim>::Pointer pMultiPatch,
        double TOL, int max_iters,
        const std::array<unsigned int, TDim>& nsampling,
        int echo_level
)
{
    typedef MultiPatch<TDim> MultiPatchType;
    typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

    for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
    {
        // make a check to make sure the point can project into the bounding box
        std::vector<double> bounding_box;
        (*it)->GetBoundingBox(bounding_box);

        bool is_in = true;
        if constexpr (TDim == 1)
        {
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL)
                          && (rPoint[0] < bounding_box[1] + TOL);
        }
        else if constexpr (TDim == 2)
        {
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL)
                          && (rPoint[0] < bounding_box[1] + TOL)
                          && (rPoint[1] > bounding_box[2] - TOL)
                          && (rPoint[1] < bounding_box[3] + TOL);
        }

        if (!is_in)
        {
            continue;
        }

        // compute a closest prediction
        int error_code = PredictVerticalProjection(rPoint, rLocalPoint, *it, nsampling);

        if (echo_level > 0)
        {
            std::cout << "Prediction local point for patch " << (*it)->Id() << ": " << rLocalPoint[0] << ", " << rLocalPoint[1] << std::endl;
        }

        // compute the vertical projection
        error_code = ComputeVerticalProjection(rPoint, rLocalPoint, rGlobalPoint, *it, TOL, max_iters, echo_level);

        if (error_code == 0)
        {
            patch_id = (*it)->Id();
            return 0;
        }
    }

    patch_id = -1;
    return 1; // can't find the projection point
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::PredictRayProjection(
        const TPointType& rPoint, const TPointType& rDirection,
        std::vector<double>& rLocalPoint,
        typename Patch<TDim>::Pointer pPatch,
        double TOL,
        const std::array<unsigned int, TDim>& nsampling
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    double cmin = 1.0e99;
    std::vector<double> xi(TDim);
    array_1d<double, 3> P;

    const double norm_dir = norm_2(rDirection);

    if constexpr (TDim == 1)
    {
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double) i) / nsampling[0];

            pControlGridFunc->GetValue(P, xi);

            const double d = norm_2(rPoint - P);

            if (d < TOL)
            {
                rLocalPoint[0] = xi[0];
                return 0; // the point is on the curve
            }

            // compute the angle
            double c = inner_prod(rPoint - P, rDirection) / (d * norm_dir);

            // select the smallest angle
            if (c < cmin)
            {
                cmin = c;
                rLocalPoint[0] = xi[0];
            }
        }
    }
    else if constexpr (TDim == 2)
    {
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double) i) / nsampling[0];

            for (int j = 0; j < nsampling[1] + 1; ++j)
            {
                xi[1] = ((double) j) / nsampling[1];

                pControlGridFunc->GetValue(P, xi);

                const double d = norm_2(rPoint - P);

                if (d < TOL)
                {
                    rLocalPoint[0] = xi[0];
                    rLocalPoint[1] = xi[1];
                    return 0; // the point is on the surface
                }

                // compute the angle
                double c = std::abs(inner_prod(rPoint - P, rDirection)) / (d * norm_dir);

                // select the smallest angle
                if (c < cmin)
                {
                    cmin = c;
                    rLocalPoint[0] = xi[0];
                    rLocalPoint[1] = xi[1];
                }
            }
        }
    }

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeRayProjection(
        const TPointType& rPoint, const TPointType& rDirection,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
        typename Patch<TDim>::Pointer pPatch,
        double TOL, int max_iters,
        int echo_level
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    std::vector<array_1d<double, 3> > dP;

    int it = 0;
    bool converged = false;

    if (echo_level > 0)
    {
        std::cout << "Global point to be projected: " << rPoint << std::endl;
    }

    if constexpr (TDim == 1)
    {
        const TPointType& v1 = rDirection;
        TPointType v2, v12, dv12dxi;
        while (!converged && (it < max_iters))
        {
            pControlGridFunc->GetValue(rGlobalPoint, rLocalPoint);

            noalias(v2) = rGlobalPoint - rPoint;

            MathUtils<double>::CrossProduct(v12, v1, v2);
            double r = norm_2(v12);

            if (echo_level > 1)
            {
                std::cout << "At iteration " << (it + 1) << ", local point = " << rLocalPoint[0] << ", " << rLocalPoint[1]
                          << ", projected global point = " << rGlobalPoint
                          << ", residual = " << r
                          << std::endl;
            }

            if (r < TOL)
            {
                converged = true;
            }
            else
            {
                pControlGridFunc->GetDerivative(dP, rLocalPoint);
                dv12dxi[0] = v1[1]*dP[0][2] - v1[2]*dP[0][1];
                dv12dxi[1] = v1[2]*dP[0][0] - v1[0]*dP[0][2];
                dv12dxi[2] = v1[0]*dP[0][1] - v1[1]*dP[0][0];

                double drdxi = inner_prod(v12, dv12dxi) / r;

                rLocalPoint[0] -= r / drdxi;

                if (echo_level > 2)
                {
                    std::cout << "tangent: " << drdxi << std::endl;
                }
            }

            ++it;
        }
    }
    else if constexpr (TDim == 2)
    {
        // TODO
        KRATOS_ERROR << "To be implemented";
    }

    if (it >= max_iters && !converged)
    {
        return 1;
    }

    if constexpr (TDim == 1)
    {
        if (rLocalPoint[0] < 0.0 || rLocalPoint[0] > 1.0)
            return 2;
    }
    else if constexpr (TDim == 2)
    {
        for (int i = 0; i < 2; ++i)
            if (rLocalPoint[i] < 0.0 || rLocalPoint[i] > 1.0)
                return 2;
    }

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeRayProjection(
        const TPointType& rPoint, const TPointType& rDirection,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<TDim>::Pointer pMultiPatch,
        double TOL, int max_iters,
        const std::array<unsigned int, TDim>& nsampling,
        int echo_level
)
{
    typedef MultiPatch<TDim> MultiPatchType;
    typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

    bool found = false;

    std::vector<std::vector<double> > LocalPoints;
    std::vector<TPointType> GlobalPoints;
    std::vector<int> patch_ids;

    TPointType GlobalPoint;
    std::vector<double> LocalPoint = rLocalPoint;
    for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
    {
        // compute a closest prediction
        int error_code = PredictRayProjection(rPoint, rDirection, LocalPoint, *it, TOL, nsampling);

        if (echo_level > 0)
        {
            std::cout << "Predicted local point for patch " << (*it)->Id() << ":";
            for (std::size_t i = 0; i < LocalPoint.size(); ++i)
                std::cout << " " << LocalPoint[i];
            std::cout << std::endl;
        }

        // compute the ray projection
        error_code = ComputeRayProjection(rPoint, rDirection, LocalPoint, GlobalPoint, *it, TOL, max_iters, echo_level);

        // store the one with acute angle
        if (error_code == 0)
        {
            found = true;
            double c = inner_prod(GlobalPoint - rPoint, rDirection);
            if (c > 0.0)
            {
                LocalPoints.push_back(LocalPoint);
                GlobalPoints.push_back(GlobalPoint);
                patch_ids.push_back((*it)->Id());
            }
        }
    }

    // look for the closest points in valid projection points
    double dist = 1e99;
    for (std::size_t i = 0; i < LocalPoints.size(); ++i)
    {
        double d = norm_2(GlobalPoints[i] - rPoint);

        if (d < dist)
        {
            patch_id = patch_ids[i];
            rLocalPoint = LocalPoints[i];
            noalias(rGlobalPoint) = GlobalPoints[i];
        }
    }

    if (found)
        return 0;

    patch_id = -1;
    return 1; // can't find the projection point
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::PredictNormalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint,
        typename Patch<TDim>::Pointer pPatch,
        const std::array<unsigned int, TDim>& nsampling
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    double dist = 1.0e99;
    std::vector<double> xi(TDim);
    array_1d<double, 3> P;

    if constexpr (TDim == 1)
    {
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double) i) / nsampling[0];

            pControlGridFunc->GetValue(P, xi);

            double d = norm_2(rPoint - P);

            if (d < dist)
            {
                dist = d;
                rLocalPoint[0] = xi[0];
            }
        }
    }
    else if constexpr (TDim == 2)
    {
        for (int i = 0; i < nsampling[0] + 1; ++i)
        {
            xi[0] = ((double) i) / nsampling[0];
            for (int j = 0; j < nsampling[1] + 1; ++j)
            {
                xi[1] = ((double) j) / nsampling[1];

                pControlGridFunc->GetValue(P, xi);

                double d = norm_2(rPoint - P);

                if (d < dist)
                {
                    dist = d;
                    rLocalPoint[0] = xi[0];
                    rLocalPoint[1] = xi[1];
                }
            }
        }
    }

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeNormalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
        typename Patch<TDim>::Pointer pPatch,
        double TOL, int max_iters,
        int echo_level
)
{
    typedef Patch<TDim> PatchType;
    typedef typename PatchType::ControlPointType ControlPointType;

    auto pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

    std::vector<array_1d<double, 3> > dP;

    int it = 0;
    bool converged = false;
    TPointType v;

    if (echo_level > 0)
    {
        std::cout << "Global point to be projected: " << rPoint << std::endl;
    }

    if constexpr (TDim == 1)
    {
        std::vector<array_1d<double, 3> > d2P;

        while (!converged && (it < max_iters))
        {
            pControlGridFunc->GetValue(rGlobalPoint, rLocalPoint);
            pControlGridFunc->GetDerivative(dP, rLocalPoint); // tangential vector

            noalias(v) = rPoint - rGlobalPoint;
            const double r = inner_prod(v, dP[0]);

            if (echo_level > 1)
            {
                std::cout << "At iteration " << (it + 1) << ", local point = " << rLocalPoint[0] << ", " << rLocalPoint[1]
                          << ", projected global point = " << rGlobalPoint
                          << ", residual = " << r
                          << std::endl;
            }

            if (fabs(r) < TOL)
            {
                converged = true;
            }
            else
            {
                pControlGridFunc->GetSecondDerivative(d2P, rLocalPoint);

                const double drdxi = -inner_prod(dP[0], dP[0]) + inner_prod(v, d2P[0]);

                rLocalPoint[0] -= r / drdxi;

                if (echo_level > 2)
                {
                    std::cout << "tangent:" << drdxi << std::endl;
                }
            }

            ++it;
        }
    }
    else if constexpr (TDim == 2)
    {
        std::vector<array_1d<double, 3> > d2P;

        Vector r(2), dr(2);
        Matrix drdxi(2, 2), inv_drdxi(2, 2);
        double det;
        while (!converged && (it < max_iters))
        {
            pControlGridFunc->GetValue(rGlobalPoint, rLocalPoint);
            pControlGridFunc->GetDerivative(dP, rLocalPoint); // tangential vector

            #ifdef CHECK_DERIVATIVES
            TPointType NewGlobalPoint;
            std::vector<double> RefLocalPoint = rLocalPoint;

            const double epsilon = 1e-6;
            rLocalPoint = RefLocalPoint;
            rLocalPoint[0] += epsilon;
            pControlGridFunc->GetValue(NewGlobalPoint, rLocalPoint);
            KRATOS_WATCH(dP[0])
            KRATOS_WATCH((NewGlobalPoint - rGlobalPoint) / epsilon)

            rLocalPoint = RefLocalPoint;
            rLocalPoint[1] += epsilon;
            pControlGridFunc->GetValue(NewGlobalPoint, rLocalPoint);
            KRATOS_WATCH(dP[1])
            KRATOS_WATCH((NewGlobalPoint - rGlobalPoint) / epsilon)

            rLocalPoint = RefLocalPoint;
            #endif

            noalias(v) = rPoint - rGlobalPoint;
            for (unsigned int i = 0; i < 2; ++i)
                r[i] = inner_prod(v, dP[i]);

            if (echo_level > 1)
            {
                std::cout << "At iteration " << (it + 1) << ", local point = " << rLocalPoint[0] << ", " << rLocalPoint[1]
                          << ", projected global point = " << rGlobalPoint
                          << ", residual = " << r
                          << std::endl;
            }

            if (norm_2(r) < TOL)
            {
                converged = true;
            }
            else
            {
                pControlGridFunc->GetSecondDerivative(d2P, rLocalPoint); // second derivatives

                #ifdef CHECK_DERIVATIVES
                std::vector<array_1d<double, 3> > NewdP(TDim);
                std::vector<double> RefLocalPoint = rLocalPoint;

                rLocalPoint = RefLocalPoint;
                rLocalPoint[0] += epsilon;
                pControlGridFunc->GetDerivative(NewdP, rLocalPoint);
                KRATOS_WATCH(d2P[0])
                KRATOS_WATCH((NewdP[0] - dP[0]) / epsilon)
                KRATOS_WATCH(d2P[2])
                KRATOS_WATCH((NewdP[1] - dP[1]) / epsilon)

                rLocalPoint = RefLocalPoint;
                rLocalPoint[1] += epsilon;
                pControlGridFunc->GetDerivative(NewdP, rLocalPoint);
                KRATOS_WATCH(d2P[2])
                KRATOS_WATCH((NewdP[0] - dP[0]) / epsilon)
                KRATOS_WATCH(d2P[1])
                KRATOS_WATCH((NewdP[1] - dP[1]) / epsilon)

                rLocalPoint = RefLocalPoint;
                #endif

                drdxi(0, 0) = -inner_prod(dP[0], dP[0]) + inner_prod(v, d2P[0]);
                drdxi(0, 1) = -inner_prod(dP[0], dP[1]) + inner_prod(v, d2P[2]);
                drdxi(1, 0) = -inner_prod(dP[1], dP[0]) + inner_prod(v, d2P[2]);
                drdxi(1, 1) = -inner_prod(dP[1], dP[1]) + inner_prod(v, d2P[1]);

                MathUtils<double>::InvertMatrix2(drdxi, inv_drdxi, det);

                noalias(dr) = -prod(inv_drdxi, r);

                for (int i = 0; i < 2; ++i)
                    rLocalPoint[i] += dr(i);

                if (echo_level > 2)
                {
                    std::cout << "tangent:" << drdxi << std::endl;
                }
            }

            ++it;
        }
    }

    #ifdef CHECK_DERIVATIVES
    KRATOS_WATCH("------------------------------")
    #endif

    if (it >= max_iters && !converged)
    {
        return 1;
    }

    if constexpr (TDim == 1)
    {
        if (rLocalPoint[0] < 0.0 || rLocalPoint[0] > 1.0)
            return 2;
    }
    else if constexpr (TDim == 2)
    {
        for (int i = 0; i < 2; ++i)
            if (rLocalPoint[i] < 0.0 || rLocalPoint[i] > 1.0)
                return 2;
    }

    if (echo_level > 0)
    {
        std::cout << "Valid projection point found on patch " << pPatch->Id() << " after " << it << " iterations" << std::endl;
    }

    #ifdef CHECK_DERIVATIVES
    KRATOS_WATCH("------------------------------")
    #endif

    return 0;
}

template<typename TPointType, int TDim>
int IsogeometricProjectionUtility<TPointType, TDim>::ComputeNormalProjection(
        const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<TDim>::Pointer pMultiPatch,
        double TOL, int max_iters,
        const std::array<unsigned int, TDim>& nsampling,
        int echo_level
)
{
    typedef MultiPatch<TDim> MultiPatchType;
    typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

    double dist = 1e99;
    bool found = false;

    TPointType GlobalPoint;
    std::vector<double> LocalPoint = rLocalPoint;
    for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
    {
        // compute a closest prediction
        int error_code = PredictNormalProjection(rPoint, LocalPoint, *it, nsampling);

        if (echo_level > 0)
        {
            std::cout << "Predicted local point for patch " << (*it)->Id() << ":";
            for (std::size_t i = 0; i < LocalPoint.size(); ++i)
                std::cout << " " << LocalPoint[i];
            std::cout << std::endl;
        }

        // compute the normal projection
        error_code = ComputeNormalProjection(rPoint, LocalPoint, GlobalPoint, *it, TOL, max_iters, echo_level);

        // update the one with shortest distance
        if (error_code == 0)
        {
            found = true;
            double d = norm_2(rPoint - GlobalPoint);
            if (d < dist)
            {
                dist = d;
                patch_id = (*it)->Id();
                rLocalPoint = LocalPoint;
                noalias(rGlobalPoint) = GlobalPoint;
            }
        }
    }

    if (found)
        return 0;

    patch_id = -1;
    return 1; // can't find the projection point
}


//// template class instantiation

template class IsogeometricProjectionUtility<Element::GeometryType::PointType::PointType, 1>;
template class IsogeometricProjectionUtility<Element::GeometryType::PointType::PointType, 2>;

template class IsogeometricProjectionUtility<array_1d<double, 3>, 1>;
template class IsogeometricProjectionUtility<array_1d<double, 3>, 2>;

} // end namespace Kratos

#undef CHECK_DERIVATIVES
