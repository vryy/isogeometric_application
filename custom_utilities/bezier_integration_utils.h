//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2014-01-28 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_BEZIER_UTILS_H_INCLUDED )
#define  KRATOS_BEZIER_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream> 

// External includes 
#include "boost/numeric/ublas/vector.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

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
class BezierUtils
{
public:
    
    static int const msBernsteinCoefs[];
    
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;
    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;

    /// Pointer definition of BezierUtils
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(BezierUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BezierUtils()
    {}

    /// Destructor.
    virtual ~BezierUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        
    /**
     * Computes Bernstein basis function B(i, p)(x) on [-1, 1]
     */
    static inline double bernstein(int i, int p, double x)
    {
        //bound checking
        if(i < 0 || i > p)
            return 0.0;
            
        double a = x;
        double b = 1 - x;
        double coeff;
        
        if(p == 0)
        {
            return 1.0;
        }
        else if(p == 1)
        {
            coeff = 1.0;
        }
        else
        {
            int cnt;
            for(int j = 0; j < p - 2; ++j)
            {
                cnt += (j + 2) / 2;
            }
            coeff = msBernsteinCoefs[cnt - 1];
        }
        
        return coeff * pow(a, i) * pow(b, p - i);
    }
    
    /**
     * Computes Bernstein basis function B(i, p)(x) on [-1, 1] using recursive iteration
     */
    static double bernstein2(int i, int p, double x)
    {
        if(i < 0 || i > p)
            return 0.0;

        if(p == 0)
        {
            return 1.0;
        }
        
        double a = x;
        double b = 1 - x;
        
        double tmp1 = bernstein2(i, p - 1, x);
        double tmp2 = bernstein2(i - 1, p - 1, x);
        return b * tmp1 + a * tmp2;
    }
    
    /**
     * Computes Bernstein basis function value & derivative B(i, p)(x) on [-1, 1]
     */
    static inline void bernstein(double& v, double& d, int i, int p, double x)
    {
        if(i < 0 || i > p)
        {
            v = 0.0;
            d = 0.0;
            return;
        }
        
        if(p == 0)
        {
            v = 1.0;
            d = 0.0;
            return;
        }
        
        double a = x;
        double b = 1 - x;
        
        double tmp1 = bernstein2(i, p - 1, x);
        double tmp2 = bernstein2(i - 1, p - 1, x);
        
        v = b * tmp1 + a * tmp2;
        d = p * (tmp2 - tmp1);
    }
    
    template<class ValuesContainerType>
    static inline void bernstein(ValuesContainerType& rS, int p, double x)
    {
        double a = x;
        double b = 1 - x;
        
        if(p == 0)
        {
            rS[0] = 1.0;
            return;
        }
        else if(p == 1)
        {
            rS[0] = b;
            rS[1] = a;
            return;
        }
        else
        {
            rS[0] = 1.0;
            rS[p] = 1.0;
            int cnt = 0;
            
            for(int j = 0; j < p - 2; ++j)
            {
                cnt += (j + 2) / 2;
            }
            
            for(int j = 0; j < p / 2; ++j)
            {
                rS[j + 1] = msBernsteinCoefs[cnt + j];
            }
            
            for(int j = 0; j < (p - 1) / 2; ++j)
            {
                rS[p - j - 1] = rS[j + 1];
            }
            
            for(int j = 0; j < p + 1; ++j)
            {
                rS[j] *= (pow(a, j) * pow(b, p - j));
            }

            return;
        }
    }
    
    template<class ValuesContainerType>
    static inline void bernstein(ValuesContainerType& rS, ValuesContainerType& rD, int p, double x)
    {
        for(int i = 0; i < p + 1; ++i)
        {
            bernstein(rS[i], rD[i], i, p, x);
        }
    }

    /**
     * Convert a modified compressed sparse row matrix to compressed sparse row matrix B <- A
     */
    static Matrix MCSR2CSR(const Matrix& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);
        
        CompressedMatrix B(n, n);
        noalias( B ) = ZeroMatrix(n, n);
        
        for(unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if(col < i)
                {
                    B.push_back(i, col, A(1, A(0, i) + j));
                }
                else
                    break;
            }
            
            // insert the diagonal
            B.push_back(i, i, A(1, i));
            
            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                if(col > i)
                {
                    B.push_back(i, col, A(1, A(0, i) + j));
                }
            }
        }
        
        B.complete_index1_data();
        
        return B;
    }
    
    /**
     * Convert a modified compressed sparse row matrix to trivial matrix
     */
    static Matrix MCSR2MAT(const Matrix& A)
    {
        unsigned int n = (unsigned int)(A(0, 0) - 1);
        
        Matrix B = ZeroMatrix(n, n);
        
        for(unsigned int i = 0; i < n; ++i)
        {
            // compute number of nonzeros for this row
            int nnz = (unsigned int)(A(0, i + 1) - A(0, i));

            B(i, i) = A(1, i);
            
            // traverse left for off-diagonal entries
            for(unsigned int j = 0; j < nnz; ++j)
            {
                unsigned int col = (unsigned int)( A(0, A(0, i) + j) );
                B(i, col) = A(1, A(0, i) + j);
            }
        }
        
        return B;
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
        buffer << "BezierUtils";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BezierUtils";
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
    BezierUtils& operator=(BezierUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BezierUtils(BezierUtils const& rOther)
    {
    }

    ///@}

}; // Class BezierUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BezierUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const BezierUtils& rThis)
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
