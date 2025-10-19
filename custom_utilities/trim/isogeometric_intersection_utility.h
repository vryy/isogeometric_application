//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Feb 2019 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_INTERSECTION_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_INTERSECTION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch_utility.h"

// #define DEBUG_INTERSECT_CURVE_PLANE

namespace Kratos
{
///@addtogroup IsogeometricApplication

///@{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// An intersection utility for Isogeometric Analysis
/**
 * Compute class to compute the intersection between patches
 */
class IsogeometricIntersectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsogeometricIntersectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricIntersectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricIntersectionUtility()
    {
    }

    /// Destructor.
    virtual ~IsogeometricIntersectionUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Check if a plane intersects with a patch. A plane does not intersect with the patch if all control points of the patch are on one side of the plane
     * Return 0 if the patch intersects with the plane; 1 if the patch is on the positive side; -1 if on the negative side
     * In the case TOption == 0, the patch intersects the plane if there is a point lie on the plane
     */
    template<int TDim, int TOption = 0>
    static std::pair<int, std::vector<int> > CheckIntersection(typename Patch<TDim>::Pointer pPatch, double A, double B, double C, double D)
    {
        typedef GridFunction<TDim, double, array_1d<double, 3> > GridFunctionType;

        typename GridFunctionType::ConstPointer pGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

        typename GridFunctionType::ControlGridType::ConstPointer pControlGrid = pGridFunc->pControlGrid();

        std::size_t nneg = 0, npos = 0, nzero = 0, npoints = pControlGrid->size();
        for (std::size_t i = 0; i < npoints; ++i)
        {
            const array_1d<double, 3>& point = (*pControlGrid)[i];
            double f = A * point[0] + B * point[1] + C * point[2] + D;
            if (f > 0.0) { ++npos; }
            else if (f < 0.0) { ++nneg; }
            else { ++nzero; }
        }

        if (TOption != 0)
        {
            if (nzero != 0)
            {
                std::pair<int, std::vector<int> > output;
                output.first = 0;
                output.second.push_back(nzero);
                output.second.push_back(nneg);
                output.second.push_back(npos);
                return output;
            }
        }

        if (nneg + nzero == npoints)
        {
            std::pair<int, std::vector<int> > output;
            output.first = -1;
            output.second.push_back(nzero);
            output.second.push_back(nneg);
            return output;
        }
        else if (npos + nzero == npoints)
        {
            std::pair<int, std::vector<int> > output;
            output.first = 1;
            output.second.push_back(nzero);
            output.second.push_back(npos);
            return output;
        }
        else
        {
            std::pair<int, std::vector<int> > output;
            output.first = 0;
            output.second.push_back(nzero);
            output.second.push_back(nneg);
            output.second.push_back(npos);
            return output;
        }
    }

    /**
     * Compute the intersection between two 1D patch
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in each curve must be given.
     * If the option_space = 0, the intersection point will be search in XY plane, and then check again on the Z direction.
     * If the option_space = 1, the intersection point will be search in YZ plane, and then check again on the X direction.
     * If the option_space = 2, the intersection point will be search in XZ plane, and then check again on the Y direction.
     * Return code definition:
     *   0: the intersection point is found and will be returned
     *   1: the intersection point is found but intersection point on curve 1 does not lie in the parametric boundary of the patch (usually [0, 1])
     *   2: the intersection point is found but intersection point on curve 2 does not lie in the parametric boundary of the patch (usually [0, 1])
     *   3: the intersection point is found but intersection points on curve 1 and 2 do not lie in the parametric boundary of the patch (usually [0, 1])
     *   4: the intersection point is not found. It is because the curves intersect in the specified plane, but does not cut in the remaining direction.
     *   5: the Newton-Raphson iteration does not converge within maximum iterations. The existence of intersection point is unknown.
     */
    static int ComputeIntersectionByNewtonRaphson(double starting_point_1,
            double starting_point_2,
            double& intersection_point_1,
            double& intersection_point_2,
            Patch<1>::Pointer pPatch1,
            Patch<1>::Pointer pPatch2,
            int max_iters,
            double TOL,
            int option_space = 1)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, double, array_1d<double, 3> >::Pointer pGridFunc1 = pPatch1->pGetGridFunction(CONTROL_POINT_COORDINATES);
        typename GridFunction<1, double, array_1d<double, 3> >::Pointer pGridFunc2 = pPatch2->pGetGridFunction(CONTROL_POINT_COORDINATES);
        // KRATOS_WATCH(typeid(*pGridFunc1->pFESpace()).name())
        // KRATOS_WATCH(typeid(*pGridFunc2->pFESpace()).name())

        std::vector<double> xi1(1), xi2(1);
        xi1[0] = starting_point_1;
        xi2[0] = starting_point_2;
        Vector res(2), dxi(2);
        Matrix J(2, 2), InvJ(2, 2);
        double DetJ;
        array_1d<double, 3> p1, p2;
        std::vector<array_1d<double, 3> > dp1, dp2;

        // determine the intersection point in flat plane
        bool converged = false;
        int it = 0;
        // KRATOS_WATCH(xi1[0])
        // KRATOS_WATCH(xi2[0])
        while (!converged && (it++ < max_iters))
        {
            pGridFunc1->GetValue(p1, xi1);
            // KRATOS_WATCH(p1)
            pGridFunc1->GetDerivative(dp1, xi1);
            // KRATOS_WATCH(dp1[0])
            pGridFunc2->GetValue(p2, xi2);
            // KRATOS_WATCH(p2)
            pGridFunc2->GetDerivative(dp2, xi2);
            // KRATOS_WATCH(dp2[0])

            if (option_space == 0)
            {
                res(0) = p1[0] - p2[0];
                res(1) = p1[1] - p2[1];

                J(0, 0) = dp1[0][0];
                J(0, 1) = -dp2[0][0];
                J(1, 0) = dp1[0][1];
                J(1, 1) = -dp2[0][1];
            }
            else if (option_space == 1)
            {
                res(0) = p1[1] - p2[1];
                res(1) = p1[2] - p2[2];

                J(0, 0) = dp1[0][1];
                J(0, 1) = -dp2[0][1];
                J(1, 0) = dp1[0][2];
                J(1, 1) = -dp2[0][2];
            }
            else if (option_space == 2)
            {
                res(0) = p1[0] - p2[0];
                res(1) = p1[2] - p2[2];

                J(0, 0) = dp1[0][0];
                J(0, 1) = -dp2[0][0];
                J(1, 0) = dp1[0][2];
                J(1, 1) = -dp2[0][2];
            }

            if (norm_2(res) < TOL)
            {
                converged = true;
                break;
            }

            // KRATOS_WATCH(norm_2(res))
            // KRATOS_WATCH(J)
            MathUtils<double>::InvertMatrix(J, InvJ, DetJ);

            noalias(dxi) = prod(InvJ, res);
            // KRATOS_WATCH(dxi)
            xi1[0] -= dxi(0);
            xi2[0] -= dxi(1);
        }

        if (!converged && (it >= max_iters))
        {
            return 5;
        }

        // KRATOS_WATCH(it)
        // KRATOS_WATCH(max_iters)

        intersection_point_1 = xi1[0];
        intersection_point_2 = xi2[0];

        // check the remaining direction
        pGridFunc1->GetValue(p1, xi1);
        pGridFunc2->GetValue(p2, xi2);

        // KRATOS_WATCH(p1)
        // KRATOS_WATCH(p2)

        if (option_space == 0)
        {
            if (std::fabs(p1[2] - p2[2]) > TOL)
            {
                return 4;    // the intersection point is found but the coordinate in Z direction is not the same
            }
        }
        else if (option_space == 1)
        {
            if (std::fabs(p1[0] - p2[0]) > TOL)
            {
                return 4;    // the intersection point is found but the coordinate in X direction is not the same
            }
        }
        else if (option_space == 2)
        {
            if (std::fabs(p1[1] - p2[1]) > TOL)
            {
                return 4;    // the intersection point is found but the coordinate in Y direction is not the same
            }
        }

        // KRATOS_WATCH(typeid(*pPatch1->pFESpace()).name())
        // KRATOS_WATCH(typeid(*pPatch2->pFESpace()).name())
        std::vector<double> pb1 = pPatch1->pFESpace()->ParametricBounds(0);
        std::vector<double> pb2 = pPatch2->pFESpace()->ParametricBounds(0);
        // KRATOS_WATCH(pb1[0])
        // KRATOS_WATCH(pb1[1])
        // KRATOS_WATCH(pb2[0])
        // KRATOS_WATCH(pb2[1])

        if ( ( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) )
                && ( (intersection_point_2 >= pb2[0]) && (intersection_point_2 <= pb2[1]) ) )
        {
            return 0;
        }

        if (  ( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) )
                && !( (intersection_point_2 >= pb2[0]) && (intersection_point_2 <= pb2[1]) ) )
        {
            return 2;
        }

        if ( !( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) )
                &&  ( (intersection_point_2 >= pb2[0]) && (intersection_point_2 <= pb2[1]) ) )
        {
            return 1;
        }

        if ( !( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) )
                && !( (intersection_point_2 >= pb2[0]) && (intersection_point_2 <= pb2[1]) ) )
        {
            return 3;
        }

        return -1; // should not come here, just to make the compiler happy
    }

    /**
     * Compute the intersection between 1D patch and a plane.
     * The plane is defined by equation Ax + By + Cz + D = 0.
     * It is assumed that the curve does not lie on the plane. Otherwise the starting point will be returned.
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve must be given.
     * Return code definition:
     *   0: the intersection point is found and will be returned
     *   1: the intersection point is found but does not lie in the parametric boundary of the patch (usually [0, 1])
     *   2: the Newton-Raphson iteration does not converge within maximum iterations. The existence of intersection point is unknown.
     */
    static int ComputeIntersectionByNewtonRaphson(double& intersection_point,
            Patch<1>::Pointer pPatch,
            double A, double B, double C, double D,
            int max_iters,
            double TOL)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, double, array_1d<double, 3> >::Pointer pGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);
#ifdef DEBUG_INTERSECT_CURVE_PLANE
        KRATOS_WATCH(typeid(*pGridFunc->pFESpace()).name())
        KRATOS_WATCH(*pPatch)
#endif

#ifdef DEBUG_INTERSECT_CURVE_PLANE
        std::vector<double> weights = pPatch->GetControlWeights();
        std::cout << "weights:";
        for (std::size_t i = 0; i < weights.size(); ++i)
        {
            std::cout << " " << weights[i];
        }
        std::cout << std::endl;
#endif

        std::vector<double> xi(1);
        xi[0] = intersection_point;
        double res, J;
        array_1d<double, 3> p;
        std::vector<array_1d<double, 3> > dp;

        // determine the intersection point in flat plane
        bool converged = false;
        int it = 0;
        while (!converged && (it++ < max_iters))
        {
#ifdef DEBUG_INTERSECT_CURVE_PLANE
            KRATOS_WATCH(xi[0])
#endif
            pGridFunc->GetValue(p, xi);
#ifdef DEBUG_INTERSECT_CURVE_PLANE
            KRATOS_WATCH(p)
#endif
            pGridFunc->GetDerivative(dp, xi);
#ifdef DEBUG_INTERSECT_CURVE_PLANE
            KRATOS_WATCH(dp[0])
#endif

            res = A * p[0] + B * p[1] + C * p[2] + D;
            J = A * dp[0][0] + B * dp[0][1] + C * dp[0][2];

            if (std::fabs(res) < TOL)
            {
                converged = true;
                break;
            }

#ifdef DEBUG_INTERSECT_CURVE_PLANE
            KRATOS_WATCH(res)
            KRATOS_WATCH(J)
#endif

            if (std::fabs(J) < TOL)
            {
                return 2;
            }

            xi[0] -= res / J;
        }

        if (!converged && (it >= max_iters))
        {
            return 2;
        }

#ifdef DEBUG_INTERSECT_CURVE_PLANE
        KRATOS_WATCH(it)
        KRATOS_WATCH(max_iters)
#endif

        intersection_point = xi[0];

        // KRATOS_WATCH(typeid(*pPatch->pFESpace()).name())
        std::vector<double> pb = pPatch->pFESpace()->ParametricBounds(0);
        // KRATOS_WATCH(pb[0])
        // KRATOS_WATCH(pb[1])

        if ( ( (intersection_point >= pb[0]) && (intersection_point <= pb[1]) ) )
        {
            return 0;
        }
        else
        {
            return 1;
        }

        return -1; // should not come here, just to make the compiler happy
    }

    /**
     * Compute the intersection between 1D patch and a plane.
     * The plane is defined by equation Ax + By + Cz + D = 0.
     * It is assumed that the curve does not lie on the plane. Otherwise the starting point will be returned.
     * This subroutine uses the bisection procedure to compute the intersection point. The two end points of the curve must lie on other side of the surface in order to return a value.
     * Return code definition:
     *   0: the intersection point is found and will be returned
     *   1: the intersection point is found but does not lie in the parametric boundary of the patch (usually [0, 1])
     *   2: the Newton-Raphson iteration does not converge within maximum iterations. The existence of intersection point is unknown.
     *   3: the two end point is on the same side of the plane. The intersection point can't be found.
     */
    static int ComputeIntersectionByBisection(double& intersection_point,
            Patch<1>::Pointer pPatch,
            double A, double B, double C, double D,
            int max_iters,
            double TOL)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, double, array_1d<double, 3> >::Pointer pGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);
#ifdef DEBUG_INTERSECT_CURVE_PLANE
        KRATOS_WATCH(typeid(*pGridFunc->pFESpace()).name())
        KRATOS_WATCH(*pPatch)
#endif

#ifdef DEBUG_INTERSECT_CURVE_PLANE
        std::vector<double> weights = pPatch->GetControlWeights();
        std::cout << "weights:";
        for (std::size_t i = 0; i < weights.size(); ++i)
        {
            std::cout << " " << weights[i];
        }
        std::cout << std::endl;
#endif

        std::vector<double> xi(1);
        double left = 0.0, right = 1.0, mid, fleft, fright, fmid;
        array_1d<double, 3> p;

        xi[0] = left;
        pGridFunc->GetValue(p, xi);
        fleft = A * p[0] + B * p[1] + C * p[2] + D;

        xi[0] = right;
        pGridFunc->GetValue(p, xi);
        fright = A * p[0] + B * p[1] + C * p[2] + D;

#ifdef DEBUG_INTERSECT_CURVE_PLANE
        KRATOS_WATCH(fleft)
        KRATOS_WATCH(fright)
#endif

        if (std::fabs(fleft) < TOL)
        {
            xi[0] = left;
            return 0;
        }

        if (std::fabs(fright) < TOL)
        {
            xi[0] = right;
            return 0;
        }

        if (fleft * fright > 0.0)
        {
            return 3;
        }

        // determine the intersection point in flat plane
        bool converged = false;
        int it = 0;
        while (!converged && (it++ < max_iters))
        {
            mid = 0.5 * (left + right);
            xi[0] = mid;
            pGridFunc->GetValue(p, xi);
            fmid = A * p[0] + B * p[1] + C * p[2] + D;

            if (std::fabs(fmid) < TOL)
            {
                converged = true;
                break;
            }

            if (fleft * fmid > 0.0)
            {
                left = mid;
            }
            else
            {
                right = mid;
            }
        }

        if (!converged && (it >= max_iters))
        {
            return 2;
        }

#ifdef DEBUG_INTERSECT_CURVE_PLANE
        KRATOS_WATCH(it)
        KRATOS_WATCH(max_iters)
#endif

        intersection_point = xi[0];

        // KRATOS_WATCH(typeid(*pPatch->pFESpace()).name())
        std::vector<double> pb = pPatch->pFESpace()->ParametricBounds(0);
        // KRATOS_WATCH(pb[0])
        // KRATOS_WATCH(pb[1])

        if ( ( (intersection_point >= pb[0]) && (intersection_point <= pb[1]) ) )
        {
            return 0;
        }
        else
        {
            return 1;
        }

        return -1; // should not come here, just to make the compiler happy
    }

    /**
     * Compute the intersection between boundaries of a 2D patch and a plane.
     * The plane is defined by equation Ax + By + Cz + D = 0.
     * It is assumed that the curve does not lie on the plane. Otherwise the starting point will be returned.
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve must be given.
     * On the return, the status code of each intersection in each boundary will be return.
     */
    static std::vector<int> ComputeIntersectionByNewtonRaphson(const std::vector<double>& starting_points,
            std::vector<std::vector<double> >& intersection_points,
            Patch<2>::Pointer pPatch,
            double A, double B, double C, double D,
            int max_iters,
            double TOL)
    {
        if (intersection_points.size() != 4)
        {
            intersection_points.resize(4);
        }

        for (std::size_t i = 0; i < 4; ++i)
            if (intersection_points[i].size() != 2)
            {
                intersection_points[i].resize(2);
            }

        std::vector<int> status(4);

        Patch<1>::Pointer pBoundaryPatch;

        pBoundaryPatch = pPatch->ConstructBoundaryPatch(_BLEFT_);
        intersection_points[0][0] = 0.0;
        intersection_points[0][1] = starting_points[0];
        status[0] = ComputeIntersectionByNewtonRaphson(intersection_points[0][1], pBoundaryPatch, A, B, C, D, max_iters, TOL);

        pBoundaryPatch = pPatch->ConstructBoundaryPatch(_BRIGHT_);
        intersection_points[1][0] = 1.0;
        intersection_points[1][1] = starting_points[1];
        status[1] = ComputeIntersectionByNewtonRaphson(intersection_points[1][1], pBoundaryPatch, A, B, C, D, max_iters, TOL);

        pBoundaryPatch = pPatch->ConstructBoundaryPatch(_BBOTTOM_);
        intersection_points[2][1] = 0.0;
        intersection_points[2][0] = starting_points[2];
        status[2] = ComputeIntersectionByNewtonRaphson(intersection_points[2][0], pBoundaryPatch, A, B, C, D, max_iters, TOL);

        pBoundaryPatch = pPatch->ConstructBoundaryPatch(_BTOP_);
        intersection_points[3][1] = 1.0;
        intersection_points[3][0] = starting_points[3];
        status[3] = ComputeIntersectionByNewtonRaphson(intersection_points[3][0], pBoundaryPatch, A, B, C, D, max_iters, TOL);

        return status;
    }

    /**
     * Compute the intersection between boundaries of a 3D patch and a plane.
     * The plane is defined by equation Ax + By + Cz + D = 0.
     * It is assumed that the curve does not lie on the plane. Otherwise the starting point will be returned.
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve must be given.
     * On the return, the status code of each intersection in each boundary will be return.
     */
    template<typename TCoordinatesType>
    static std::vector<int> ComputeIntersectionByNewtonRaphson(const std::vector<double>& starting_points,
            std::vector<TCoordinatesType>& intersection_points,
            Patch<3>::Pointer pPatch,
            double A, double B, double C, double D,
            int max_iters,
            double TOL)
    {
        std::vector<Patch<1>::Pointer> pEdgePatches = MultiPatchUtility::ConstructEdgePatches(pPatch);

        std::size_t nedges = pEdgePatches.size();
        if (intersection_points.size() != nedges)
        {
            intersection_points.resize(nedges);
        }

        for (std::size_t i = 0; i < nedges; ++i)
            if (intersection_points[i].size() != 3)
            {
                intersection_points[i].resize(3);
            }

        std::vector<int> status(nedges);

        for (std::size_t i = 0; i < nedges; ++i)
        {
            double tmp_point = starting_points[i];
            status[i] = ComputeIntersectionByNewtonRaphson(tmp_point, pEdgePatches[i], A, B, C, D, max_iters, TOL);

            if (i == 0) {intersection_points[i][0] = 0.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 0.0;}
            else if (i == 1) {intersection_points[i][0] = 1.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 0.0;}
            else if (i == 2) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 0.0; intersection_points[i][2] = 0.0;}
            else if (i == 3) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 1.0; intersection_points[i][2] = 0.0;}
            else if (i == 4) {intersection_points[i][0] = 0.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 1.0;}
            else if (i == 5) {intersection_points[i][0] = 1.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 1.0;}
            else if (i == 6) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 0.0; intersection_points[i][2] = 1.0;}
            else if (i == 7) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 1.0; intersection_points[i][2] = 1.0;}
            else if (i == 8) {intersection_points[i][0] = 0.0; intersection_points[i][1] = 0.0; intersection_points[i][2] = tmp_point;}
            else if (i == 9) {intersection_points[i][0] = 0.0; intersection_points[i][1] = 1.0; intersection_points[i][2] = tmp_point;}
            else if (i == 10) {intersection_points[i][0] = 1.0; intersection_points[i][1] = 0.0; intersection_points[i][2] = tmp_point;}
            else if (i == 11) {intersection_points[i][0] = 1.0; intersection_points[i][1] = 1.0; intersection_points[i][2] = tmp_point;}
        }

        return status;
    }

    /**
     * Compute the intersection between boundaries of a 3D patch and a plane.
     * The plane is defined by equation Ax + By + Cz + D = 0.
     * It is assumed that the curve does not lie on the plane. Otherwise the starting point will be returned.
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve must be given.
     * On the return, the status code of each intersection in each boundary will be return.
     */
    template<typename TCoordinatesType>
    static std::vector<int> ComputeIntersectionByBisection(std::vector<TCoordinatesType>& intersection_points,
            Patch<3>::Pointer pPatch,
            double A, double B, double C, double D,
            int max_iters,
            double TOL)
    {
        std::vector<Patch<1>::Pointer> pEdgePatches = MultiPatchUtility::ConstructEdgePatches(pPatch);

        std::size_t nedges = pEdgePatches.size();
        if (intersection_points.size() != nedges)
        {
            intersection_points.resize(nedges);
        }

        for (std::size_t i = 0; i < nedges; ++i)
            if (intersection_points[i].size() != 3)
            {
                intersection_points[i].resize(3);
            }

        std::vector<int> status(nedges);

        for (std::size_t i = 0; i < nedges; ++i)
        {
            double tmp_point;
            status[i] = ComputeIntersectionByBisection(tmp_point, pEdgePatches[i], A, B, C, D, max_iters, TOL);

            if (i == 0) {intersection_points[i][0] = 0.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 0.0;}
            else if (i == 1) {intersection_points[i][0] = 1.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 0.0;}
            else if (i == 2) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 0.0; intersection_points[i][2] = 0.0;}
            else if (i == 3) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 1.0; intersection_points[i][2] = 0.0;}
            else if (i == 4) {intersection_points[i][0] = 0.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 1.0;}
            else if (i == 5) {intersection_points[i][0] = 1.0; intersection_points[i][1] = tmp_point; intersection_points[i][2] = 1.0;}
            else if (i == 6) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 0.0; intersection_points[i][2] = 1.0;}
            else if (i == 7) {intersection_points[i][0] = tmp_point; intersection_points[i][1] = 1.0; intersection_points[i][2] = 1.0;}
            else if (i == 8) {intersection_points[i][0] = 0.0; intersection_points[i][1] = 0.0; intersection_points[i][2] = tmp_point;}
            else if (i == 9) {intersection_points[i][0] = 0.0; intersection_points[i][1] = 1.0; intersection_points[i][2] = tmp_point;}
            else if (i == 10) {intersection_points[i][0] = 1.0; intersection_points[i][1] = 0.0; intersection_points[i][2] = tmp_point;}
            else if (i == 11) {intersection_points[i][0] = 1.0; intersection_points[i][1] = 1.0; intersection_points[i][2] = tmp_point;}
        }

        return status;
    }

    /**
     * Compute the intersection between a 1D patch and a 2D patch
     * A 1D patch representes a 3D curve and a 2D patch represents a 3D surface. In case that the curve and surface is co-planar, the result may be unpredictable.
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve and surface must be given.
     * Return code definition:
     *   0: the intersection point is found and will be returned
     *   1: the intersection point is found but intersection point on the curve does not lie in the parametric boundary of the patch (usually [0, 1])
     *   2: the intersection point is found but intersection point on the surface does not lie in the parametric boundary of the patch (usually [0, 1] x [0, 1])
     *   3: the Newton-Raphson iteration does not converge within maximum iterations. The existence of intersection point is unknown.
     */
    static int ComputeIntersectionByNewtonRaphson(double starting_point_1,
            const std::vector<double>& starting_point_2,
            double& intersection_point_1,
            std::vector<double>& intersection_point_2,
            Patch<1>::Pointer pPatch1,
            Patch<2>::Pointer pPatch2,
            int max_iters,
            double TOL)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, double, array_1d<double, 3> >::Pointer pGridFunc1 = pPatch1->pGetGridFunction(CONTROL_POINT_COORDINATES);
        typename GridFunction<2, double, array_1d<double, 3> >::Pointer pGridFunc2 = pPatch2->pGetGridFunction(CONTROL_POINT_COORDINATES);
        // KRATOS_WATCH(typeid(*pGridFunc1->pFESpace()).name())
        // KRATOS_WATCH(typeid(*pGridFunc2->pFESpace()).name())

        std::vector<double> xi1(1), xi2(2);
        xi1[0] = starting_point_1;
        xi2[0] = starting_point_2[0];
        xi2[1] = starting_point_2[1];
        Vector res(3), dxi(3);
        Matrix J(3, 3), InvJ(3, 3);
        double DetJ;
        array_1d<double, 3> p1, p2;
        std::vector<array_1d<double, 3> > dp1, dp2;

        // determine the intersection point in flat plane
        bool converged = false;
        int it = 0;
        // KRATOS_WATCH(xi1[0])
        // KRATOS_WATCH(xi2[0])
        while (!converged && (it++ < max_iters))
        {
            pGridFunc1->GetValue(p1, xi1);
            // KRATOS_WATCH(p1)
            pGridFunc1->GetDerivative(dp1, xi1);
            // KRATOS_WATCH(dp1[0])
            pGridFunc2->GetValue(p2, xi2);
            // KRATOS_WATCH(p2)
            pGridFunc2->GetDerivative(dp2, xi2);
            // KRATOS_WATCH(dp2[0])

            res(0) = p1[0] - p2[0];
            res(1) = p1[1] - p2[1];
            res(2) = p1[2] - p2[2];

            J(0, 0) = dp1[0][0];
            J(0, 1) = -dp2[0][0];
            J(0, 2) = -dp2[1][0];

            J(1, 0) = dp1[0][1];
            J(1, 1) = -dp2[0][1];
            J(1, 2) = -dp2[1][1];

            J(2, 0) = dp1[0][2];
            J(2, 1) = -dp2[0][2];
            J(2, 2) = -dp2[1][2];

            if (norm_2(res) < TOL)
            {
                converged = true;
                break;
            }

            // KRATOS_WATCH(norm_2(res))
            // KRATOS_WATCH(J)
            MathUtils<double>::InvertMatrix(J, InvJ, DetJ);

            if (std::fabs(DetJ) < TOL)
            {
                return 3;
            }

            noalias(dxi) = prod(InvJ, res);
            // KRATOS_WATCH(dxi)
            xi1[0] -= dxi(0);
            xi2[0] -= dxi(1);
            xi2[1] -= dxi(2);
        }

        if (!converged && (it >= max_iters))
        {
            return 3;
        }

        // KRATOS_WATCH(it)
        // KRATOS_WATCH(max_iters)

        intersection_point_1 = xi1[0];

        if (intersection_point_2.size() != 2)
        {
            intersection_point_2.resize(2);
        }
        intersection_point_2[0] = xi2[0];
        intersection_point_2[1] = xi2[1];

        // KRATOS_WATCH(p1)
        // KRATOS_WATCH(p2)

        // KRATOS_WATCH(typeid(*pPatch1->pFESpace()).name())
        // KRATOS_WATCH(typeid(*pPatch2->pFESpace()).name())
        std::vector<double> pb1 = pPatch1->pFESpace()->ParametricBounds(0);
        std::vector<double> pb2 = pPatch2->pFESpace()->ParametricBounds(0);
        // KRATOS_WATCH(pb1[0])
        // KRATOS_WATCH(pb1[1])
        // KRATOS_WATCH(pb2[0])
        // KRATOS_WATCH(pb2[1])

        if ( ( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) )
                && ( (intersection_point_2[0] >= pb2[0]) && (intersection_point_2[0] <= pb2[1]) )
                && ( (intersection_point_2[1] >= pb2[2]) && (intersection_point_2[1] <= pb2[3]) ) )
        {
            return 0;
        }

        if ( !( (intersection_point_1 >= pb1[0]) && (intersection_point_1 <= pb1[1]) ) )
        {
            return 1;
        }

        if ( !( ( (intersection_point_2[0] >= pb2[0]) && (intersection_point_2[0] <= pb2[1]) )
                && ( (intersection_point_2[1] >= pb2[2]) && (intersection_point_2[1] <= pb2[3]) ) ) )
        {
            return 2;
        }

        return -1; // should not come here, just to make the compiler happy
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "IsogeometricIntersectionUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricIntersectionUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IsogeometricIntersectionUtility& operator=(IsogeometricIntersectionUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricIntersectionUtility(IsogeometricIntersectionUtility const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricIntersectionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, IsogeometricIntersectionUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const IsogeometricIntersectionUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef DEBUG_INTERSECT_CURVE_PLANE

#endif // KRATOS_ISOGEOMETRIC_INTERSECTION_UTILITY_H_INCLUDED
