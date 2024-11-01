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
class IsogeometricProjectionUtility
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

    /// Compute the prediction of vertical projection (along z axis) of a point on a surface patch
    template<typename TPointType>
    static int PredictVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint,
                                         typename Patch<2>::Pointer pPatch,
                                         std::size_t nsampling1, std::size_t nsampling2);

    /// Compute the vertical projection (along z axis) of a point on a surface patch
    template<typename TPointType>
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                         typename Patch<2>::Pointer pPatch,
                                         double TOL, int max_iters,
                                         int echo_level);

    /// Compute the vertical projection (along z axis) of a point on a surface multipatch
    /// The rLocalPoint shall be initialized to a good value to find out the vertical projection
    template<typename TPointType>
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                         typename MultiPatch<2>::Pointer pMultiPatch,
                                         double TOL, int max_iters,
                                         int echo_level);

    /// Compute the vertical projection (along z axis) of a point on a surface multipatch
    template<typename TPointType>
    static int ComputeVerticalProjection(const TPointType& rPoint,
                                         std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                         typename MultiPatch<2>::Pointer pMultiPatch,
                                         double TOL, int max_iters,
                                         std::size_t nsampling1, std::size_t nsampling2,
                                         int echo_level);

    /******************************************************************/
    /************************** RAY TRACING ***************************/
    /******************************************************************/

    /// Compute the prediction of the intersection of a ray on a line/surface patch
    template<typename TPointType, int TDim>
    static int PredictRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint,
                                    typename Patch<TDim>::Pointer pPatch,
                                    double TOL,
                                    const std::array<int, TDim>& nsampling);

    /// Compute the intersection of a ray on a line/surface patch
    template<typename TPointType, int TDim>
    static int ComputeRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                    typename Patch<TDim>::Pointer pPatch,
                                    double TOL, int max_iters,
                                    int echo_level);

    /// Compute the intersection of a ray on a line/surface multipatch
    template<typename TPointType, int TDim>
    static int ComputeRayProjection(const TPointType& rPoint, const TPointType& rDirection,
                                    std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                    typename MultiPatch<TDim>::Pointer pMultiPatch,
                                    double TOL, int max_iters,
                                    const std::array<int, TDim>& nsampling,
                                    int echo_level);

    /******************************************************************/
    /************** NORMAL PROJECTION (ALONG Z AXIS) ******************/
    /******************************************************************/

    /// Compute the prediction of normal projection of a point on a line patch
    template<typename TPointType, int TDim>
    static int PredictNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint,
                                       typename Patch<TDim>::Pointer pPatch,
                                       const std::array<int, TDim>& nsampling);

    /// Compute the normal projection of a point on a line patch
    template<typename TPointType, int TDim>
    static int ComputeNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint, TPointType& rGlobalPoint,
                                       typename Patch<TDim>::Pointer pPatch,
                                       double TOL, int max_iters,
                                       int echo_level);

    /// Compute the normal projection of a point on a line/surface multipatch
    template<typename TPointType, int TDim>
    static int ComputeNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                       typename MultiPatch<TDim>::Pointer pMultiPatch,
                                       double TOL, int max_iters,
                                       const std::array<int, TDim>& nsampling,
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
inline std::istream& operator >>(std::istream& rIStream, IsogeometricProjectionUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const IsogeometricProjectionUtility& rThis)
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
