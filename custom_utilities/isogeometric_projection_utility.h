//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jan 2021 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"

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
template<typename TPointType, int TDim>
class KRATOS_API(ISOGEOMETRIC_APPLICATION) IsogeometricProjectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsogeometricProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricProjectionUtility);

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

    /******************************************************************/
    /************* VERTICAL PROJECTION (ALONG Z AXIS) *****************/
    /******************************************************************/

    /// Compute the prediction of vertical projection (along z axis) of a point on a surface patch in 3D
    /// and along y axis on a line patch in 2D
    static int PredictVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint,
                                         typename Patch<TDim>::Pointer pPatch,
                                         const std::array<unsigned int, TDim>& nsampling);

    /// Compute the vertical projection (along z axis) of a point on a surface patch in 3D
    /// and along y axis on a line patch in 2D
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                         typename Patch<TDim>::Pointer pPatch,
                                         double TOL, int max_iters,
                                         int echo_level);

    /// Compute the vertical projection (along z axis) of a point on a surface multipatch in 3D
    /// and along y axis on a line patch in 2D
    /// The rLocalPoint shall be initialized to a good value to find out the vertical projection
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                         typename MultiPatch<TDim>::Pointer pMultiPatch,
                                         double TOL, int max_iters,
                                         int echo_level);

    /// Compute the vertical projection (along z axis) of a point on a surface multipatch in 3D
    /// and along y axis on a line patch in 2D
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                         typename MultiPatch<TDim>::Pointer pMultiPatch,
                                         double TOL, int max_iters,
                                         const std::array<unsigned int, TDim>& nsampling,
                                         int echo_level);

    /******************************************************************/
    /************************** RAY TRACING ***************************/
    /******************************************************************/

    /// Compute the prediction of the intersection of a ray on a line/surface patch
    static int PredictRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint,
                                    typename Patch<TDim>::Pointer pPatch,
                                    double TOL,
                                    const std::array<unsigned int, TDim>& nsampling);

    /// Compute the intersection of a ray on a line/surface patch
    static int ComputeRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                    typename Patch<TDim>::Pointer pPatch,
                                    double TOL, int max_iters,
                                    int echo_level);

    /// Compute the intersection of a ray on a line/surface multipatch
    static int ComputeRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                    typename MultiPatch<TDim>::Pointer pMultiPatch,
                                    double TOL, int max_iters,
                                    const std::array<unsigned int, TDim>& nsampling,
                                    int echo_level);

    /******************************************************************/
    /*********************** NORMAL PROJECTION ************************/
    /******************************************************************/

    /// Compute the prediction of normal projection of a point on a line patch
    static int PredictNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint,
                                       typename Patch<TDim>::Pointer pPatch,
                                       const std::array<unsigned int, TDim>& nsampling);

    /// Compute the normal projection of a point on a line patch
    static int ComputeNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                       typename Patch<TDim>::Pointer pPatch,
                                       double TOL, int max_iters,
                                       int echo_level);

    /// Compute the normal projection of a point on a line/surface multipatch
    static int ComputeNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                       typename MultiPatch<TDim>::Pointer pMultiPatch,
                                       double TOL, int max_iters,
                                       const std::array<unsigned int, TDim>& nsampling,
                                       int echo_level);

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
template<typename TPointType, int TDim>
inline std::istream& operator >>(std::istream& rIStream, IsogeometricProjectionUtility<TPointType, TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<typename TPointType, int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const IsogeometricProjectionUtility<TPointType, TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_PROJECTION_UTILITY_H_INCLUDED
