//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014 Jan 28 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BEZIER_TEST_UTILS_H_INCLUDED )
#define  KRATOS_BEZIER_TEST_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream> 

// External includes 
#include "boost/numeric/ublas/vector.hpp"
// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/bezier_utils.h"

namespace Kratos
{
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
/** Detail class definition.
 */
class BezierTestUtils
{
public:
    ///@name Type Definitions
    ///@{

//      typedef std::vector<double>       ValueContainerType;
    typedef boost::numeric::ublas::vector<double> ValueContainerType;

    typedef ModelPart::NodesContainerType NodesArrayType;
    
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename Element::GeometryType::PointType::CoordinatesArrayType CoordinatesArrayType;

    /// Pointer definition of BezierTestUtils
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(BezierTestUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BezierTestUtils()
    {}

    /// Destructor.
    virtual ~BezierTestUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void TestBernsteinBasis()
    {
        int n = 3;
        Vector v(n + 1);
        BezierUtils::bernstein(v, n, 0.5);
        
        KRATOS_WATCH(v)
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
        buffer << "BezierTestUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BezierTestUtils";
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BezierTestUtils& operator=(BezierTestUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BezierTestUtils(BezierTestUtils const& rOther)
    {
    }

    ///@}

}; // Class BezierTestUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BezierTestUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const BezierTestUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_BEZIER_UTILS_H_INCLUDED  defined 
