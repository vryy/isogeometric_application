//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jan 2021 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch_utility.h"


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

/// Utility to compute the projection for Isogeometric Analysis
/**
 * Compute class to compute the projection on patch and multipatch
 */
class IsogeometricProjectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsogeometricProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricProjectionUtility);

    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef NodeType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricProjectionUtility()
    {
    }

    /// Destructor.
    virtual ~IsogeometricProjectionUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Compute a prediction of vertical projection of a point on patch
    static int PredictVerticalProjection(const PointType& rPoint,
        std::vector<double>& rLocalPoint,
        typename Patch<2>::Pointer pPatch,
        std::size_t nsampling1, std::size_t nsampling2)
    {
        typedef Patch<2> PatchType;
        typedef typename PatchType::ControlPointType ControlPointType;

        typename GridFunction<2, array_1d<double, 3> >::ConstPointer pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

        double dist = 1.0e99;
        std::vector<double> xi(2);
        array_1d<double, 3> P;

        // KRATOS_WATCH(nsampling1)
        // KRATOS_WATCH(nsampling2)

        for (std::size_t i = 0; i < nsampling1+1; ++i)
        {
            xi[0] = ((double) i) / nsampling1;
            for (std::size_t j = 0; j < nsampling2+1; ++j)
            {
                xi[1] = ((double) j) / nsampling2;

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

        return 0;
    }

    /// Compute the vertical projection of a point on patch
    static int ComputeVerticalProjection(const PointType& rPoint,
        std::vector<double>& rLocalPoint, PointType& rGlobalPoint,
        typename Patch<2>::Pointer pPatch,
        double TOL, int max_iters,
        int echo_level)
    {
        typedef Patch<2> PatchType;
        typedef typename PatchType::ControlPointType ControlPointType;

        typename GridFunction<2, array_1d<double, 3> >::ConstPointer pControlGridFunc = pPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

        std::vector<array_1d<double, 3> > dP;

        int it = 0;
        bool converged = false;
        double a00, a01, a10, a11;

        if (echo_level > 0)
        {
            std::cout << "Global point to be projected: " << rPoint << std::endl;
        }

        while(!converged && (it < max_iters))
        {
            pControlGridFunc->GetValue(rGlobalPoint, rLocalPoint);

            double dist = sqrt(pow(rPoint[0] - rGlobalPoint[0], 2) + pow(rPoint[1] - rGlobalPoint[1], 2));

            if (echo_level > 1)
            {
                std::cout << "At iteration " << (it+1) << ", local point = " << rLocalPoint[0] << ", " << rLocalPoint[1]
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

                double det = a00*a11 - a01*a10;
                rLocalPoint[0] += ( a11*(rPoint[0] - rGlobalPoint[0]) - a01*(rPoint[1] - rGlobalPoint[1])) / det;
                rLocalPoint[1] += (-a10*(rPoint[0] - rGlobalPoint[0]) + a00*(rPoint[1] - rGlobalPoint[1])) / det;

                if (rLocalPoint[0] < 0.0) rLocalPoint[0] = 0.0;
                if (rLocalPoint[0] > 1.0) rLocalPoint[0] = 1.0;
                if (rLocalPoint[1] < 0.0) rLocalPoint[1] = 0.0;
                if (rLocalPoint[1] > 1.0) rLocalPoint[1] = 1.0;

                if (echo_level > 2)
                {
                    std::cout << "tangent matrix: [[" << a00 << ", " << a01 << "], [" << a10 << ", " << a11 << "]" << std::endl;
                }
            }

            ++it;
        }

        if (it >= max_iters && !converged)
            return 1;

        return 0;
    }

    /// Compute the vertical projection of a point on multipatch
    /// The rLocalPoint shall be initialized to a good value to find out the vertical projection
    static int ComputeVerticalProjection(const PointType& rPoint,
        std::vector<double>& rLocalPoint, PointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<2>::Pointer pMultiPatch,
        double TOL, int max_iters,
        int echo_level)
    {
        typedef MultiPatch<2> MultiPatchType;
        typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

        std::vector<double> InitialLocalPoint = rLocalPoint;

        std::vector<double> bounding_box;
        for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
        {
            // make a check to make sure the point can project into the bounding box
            (*it)->GetBoundingBox(bounding_box);

            bool is_in = true;
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL);
            is_in = is_in && (rPoint[0] < bounding_box[1] + TOL);
            is_in = is_in && (rPoint[1] > bounding_box[2] - TOL);
            is_in = is_in && (rPoint[1] < bounding_box[3] + TOL);

            if (!is_in)
                continue;

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

    /// Compute the vertical projection of a point on multipatch
    static int ComputeVerticalProjection(const PointType& rPoint,
        std::vector<double>& rLocalPoint, PointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<2>::Pointer pMultiPatch,
        double TOL, int max_iters,
        int nsampling1, int nsampling2,
        int echo_level)
    {
        typedef MultiPatch<2> MultiPatchType;
        typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

        std::vector<double> bounding_box;
        for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
        {
            // make a check to make sure the point can project into the bounding box
            (*it)->GetBoundingBox(bounding_box);

            bool is_in = true;
            is_in = is_in && (rPoint[0] > bounding_box[0] - TOL);
            is_in = is_in && (rPoint[0] < bounding_box[1] + TOL);
            is_in = is_in && (rPoint[1] > bounding_box[2] - TOL);
            is_in = is_in && (rPoint[1] < bounding_box[3] + TOL);

            if (!is_in)
                continue;

            // compute a closest prediction
            int error_code = PredictVerticalProjection(rPoint, rLocalPoint, *it, nsampling1, nsampling2);

            if (echo_level > 0)
                std::cout << "Prediction local point for patch " << (*it)->Id() << ": " << rLocalPoint[0] << ", " << rLocalPoint[1] << std::endl;

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
        buffer << "IsogeometricProjectionUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricProjectionUtility";
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
    IsogeometricProjectionUtility& operator=(IsogeometricProjectionUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricProjectionUtility(IsogeometricProjectionUtility const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricProjectionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, IsogeometricProjectionUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const IsogeometricProjectionUtility& rThis)
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

#endif // KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED

