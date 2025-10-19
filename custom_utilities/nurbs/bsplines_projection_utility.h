//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Oct 2024 $
//
//

#if !defined(KRATOS_BSPLINES_PROJECTION_UTILITY_H_INCLUDED )
#define  KRATOS_BSPLINES_PROJECTION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/isogeometric_projection_utility.h"

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

/**
 * Compute class to compute the projection on B-Splines patch and multipatch
 */
template<typename TPointType, int TDim>
class KRATOS_API(ISOGEOMETRIC_APPLICATION) BSplinesProjectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BSplinesProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesProjectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BSplinesProjectionUtility()
    {
    }

    /// Destructor.
    virtual ~BSplinesProjectionUtility()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Compute the normal projection of a point on a line multipatch
    static int ComputeNormalProjection(const TPointType& rPoint,
                                       std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
                                       typename MultiPatch<TDim>::Pointer pMultiPatch, const int nsampling,
                                       double TOL, int max_iters,
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
        return "BSplinesProjectionUtility";
    }

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
    BSplinesProjectionUtility& operator=(BSplinesProjectionUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BSplinesProjectionUtility(BSplinesProjectionUtility const& rOther)
    {
    }

    ///@}

}; // Class BSplinesProjectionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_BSPLINES_PROJECTION_UTILITY_H_INCLUDED
