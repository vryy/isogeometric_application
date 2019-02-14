//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Feb 2019 $
//   Revision:            $Revision: 1.0 $
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

/// Short class definition.
/**
 * Abstract class to compute the intersection between patches
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
    static int ComputeIntersection(const double& starting_point_1,
            const double& starting_point_2,
            double& intersection_point_1,
            double& intersection_point_2,
            Patch<1>::Pointer pPatch1,
            Patch<1>::Pointer pPatch2,
            const int& max_iters,
            const double& TOL,
            const int& option_space = 1)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, array_1d<double, 3> >::Pointer pGridFunc1 = pPatch1->pGetGridFunction(CONTROL_POINT_COORDINATES);
        typename GridFunction<1, array_1d<double, 3> >::Pointer pGridFunc2 = pPatch2->pGetGridFunction(CONTROL_POINT_COORDINATES);
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
                return 4; // the intersection point is found but the coordinate in Z direction is not the same
        }
        else if (option_space == 1)
        {
            if (std::fabs(p1[0] - p2[0]) > TOL)
                return 4; // the intersection point is found but the coordinate in X direction is not the same
        }
        else if (option_space == 2)
        {
            if (std::fabs(p1[1] - p2[1]) > TOL)
                return 4; // the intersection point is found but the coordinate in Y direction is not the same
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
    static int ComputeIntersection(const double& starting_point,
            double& intersection_point,
            Patch<1>::Pointer pPatch,
            const double& A, const double& B, const double& C, const double& D,
            const int& max_iters,
            const double& TOL)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, array_1d<double, 3> >::Pointer pGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);
        // KRATOS_WATCH(typeid(*pGridFunc->pFESpace()).name())

        std::vector<double> xi(1);
        xi[0] = starting_point;
        double res, J;
        array_1d<double, 3> p;
        std::vector<array_1d<double, 3> > dp;

        // determine the intersection point in flat plane
        bool converged = false;
        int it = 0;
        // KRATOS_WATCH(xi[0])
        while (!converged && (it++ < max_iters))
        {
            pGridFunc->GetValue(p, xi);
            // KRATOS_WATCH(p)
            pGridFunc->GetDerivative(dp, xi);
            // KRATOS_WATCH(dp[0])

            res = A*p[0] + B*p[1] + C*p[2] + D;
            J = A*dp[0][0] + B*dp[0][1] + C*dp[0][2];

            if (std::fabs(res) < TOL)
            {
                converged = true;
                break;
            }

            KRATOS_WATCH(res)
            // KRATOS_WATCH(J)

            xi[0] -= res/J;
        }

        if (!converged && (it >= max_iters))
        {
            return 2;
        }

        // KRATOS_WATCH(it)
        // KRATOS_WATCH(max_iters)

        intersection_point = xi[0];

        // KRATOS_WATCH(typeid(*pPatch->pFESpace()).name())
        std::vector<double> pb = pPatch->pFESpace()->ParametricBounds(0);
        // KRATOS_WATCH(pb[0])
        // KRATOS_WATCH(pb[1])

        if ( ( (intersection_point >= pb[0]) && (intersection_point <= pb[1]) ) )
            return 0;
        else
            return 1;

        return -1; // should not come here, just to make the compiler happy
    }

    /**
     * Compute the intersection between a 1D patch and a 2D patch
     * This subroutine uses the Newton-Raphson procedure to compute the intersection point. A starting value in the curve and surface must be given.
     * Return code definition:
     *   0: the intersection point is found and will be returned
     *   1: the intersection point is found but intersection point on the curve does not lie in the parametric boundary of the patch (usually [0, 1])
     *   2: the intersection point is found but intersection point on the surface does not lie in the parametric boundary of the patch (usually [0, 1] x [0, 1])
     *   3: the Newton-Raphson iteration does not converge within maximum iterations. The existence of intersection point is unknown.
     */
    static int ComputeIntersection(const double& starting_point_1,
            const std::vector<double>& starting_point_2,
            double& intersection_point_1,
            std::vector<double>& intersection_point_2,
            Patch<1>::Pointer pPatch1,
            Patch<2>::Pointer pPatch2,
            const int& max_iters,
            const double& TOL)
    {
        // extract the control point grid function of both curves
        typename GridFunction<1, array_1d<double, 3> >::Pointer pGridFunc1 = pPatch1->pGetGridFunction(CONTROL_POINT_COORDINATES);
        typename GridFunction<2, array_1d<double, 3> >::Pointer pGridFunc2 = pPatch2->pGetGridFunction(CONTROL_POINT_COORDINATES);
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
            intersection_point_2.resize(2);
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
            return 1;

        if ( !( ( (intersection_point_2[0] >= pb2[0]) && (intersection_point_2[0] <= pb2[1]) )
             && ( (intersection_point_2[1] >= pb2[2]) && (intersection_point_2[1] <= pb2[3]) ) ) )
        {
            return 2;
        }

        return -1; // should not come here, just to make the compiler happy
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

#endif // KRATOS_ISOGEOMETRIC_INTERSECTION_UTILITY_H_INCLUDED

