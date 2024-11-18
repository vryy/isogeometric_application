//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PBSPLINES_BASIS_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PBSPLINES_BASIS_FUNCTION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "containers/data_value_container.h"
#include "custom_utilities/control_point.h"
#include "isogeometric_application_variables.h"

namespace Kratos
{

template<class TBasisFunctionType, class TVariableType>
struct PBSplinesBasisFunction_InitializeValue_Helper
{
    static void Initialize(TBasisFunctionType& r_bf, const TVariableType& rVariable)
    {
        KRATOS_ERROR << "Not yet implemented";
    }

    static void Initialize(TBasisFunctionType& r_bf, const TVariableType& rVariable, typename TBasisFunctionType::Pointer p_ref_bf)
    {
        KRATOS_ERROR << "Not yet implemented";
    }
};

/**
 * Abstract class for point-based Splines basis function.
 * Each basis function associates with a control point via CONTROL_POINT variable.
 */
template<typename TCellType>
class PBSplinesBasisFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PBSplinesBasisFunction);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;

    typedef TCellType CellType;
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::ConstPointer const_cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename cell_container_t::const_iterator cell_const_iterator;

    /// Empty constructor for serialization
    PBSplinesBasisFunction()
        : mId(0), mEquationId(-1)
    {}

    /// Constructor with Id
    PBSplinesBasisFunction(std::size_t Id)
        : mId(Id), mEquationId(-1)
    {}

    /// Destructor
    ~PBSplinesBasisFunction()
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "PBSplinesBasisFunction" << TDim << "D " << this->Id() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /**************************************************************************
                            MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Add a cell support this basis function to the list
    cell_iterator AddCell(cell_t p_cell)
    {
        for (cell_iterator it = cell_begin(); it != cell_end(); ++it)
            if (*it == p_cell)
            {
                return it;
            }
        return mpCells.insert(p_cell).first;
    }

    /// Remove the cell from the list
    void RemoveCell(cell_t p_cell)
    {
        for (cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if (*it == p_cell)
            {
                mpCells.erase(it);
                break;
            }
        }
    }

    /// Remove the cell from the list
    void RemoveCell(CellType& r_cell)
    {
        for (cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if (&(*(*it)) == &r_cell)
            {
                mpCells.erase(it);
                break;
            }
        }
    }

    /**************************************************************************
                            ACCESS SUBROUTINES
    **************************************************************************/

    /// Iterators to the supporting cell of this basis function
    cell_iterator cell_begin() {return mpCells.begin();}
    cell_const_iterator cell_begin() const {return mpCells.begin();}
    cell_iterator cell_end() {return mpCells.end();}
    cell_const_iterator cell_end() const {return mpCells.end();}

    /// Get the Id of this basis function. Each basis function should have unique Id within a patch
    std::size_t Id() const {return mId;}

    /// Get the equation Id of this basis function. Each basis function should have unique equation Id accross patches
    std::size_t EquationId() const {return mEquationId;}

    /// Set the equation Id for this basis function. One shall use this function only in the enumeration process
    void SetEquationId(std::size_t EquationId) {mEquationId = EquationId;}

    /// Get the weight associated with the control point
    double Weight() const
    {
        return this->GetValue(CONTROL_POINT).W();
    }

    /// Get the bounding box (=support domain) of this basis function
    virtual std::vector<double> GetBoundingBox() const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the value of point-based B-splines basis function
    virtual double GetValueAt(const std::vector<double>& xi) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the value of point-based B-splines basis function
    virtual void GetValueAt(double& res, const std::vector<double>& xi) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the derivative of point-based B-splines basis function
    virtual std::vector<double> GetDerivativeAt(const std::vector<double>& xi) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the derivative of point-based B-splines basis function
    virtual void GetDerivativeAt(std::vector<double>& res, const std::vector<double>& xi) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Access the internal data container, be very careful with this function
    DataValueContainer& Data() {return mData;}
    const DataValueContainer& Data() const {return mData;}

    /**************************************************************************
                            CONTROL VALUES
      IT IS IMPORTANT TO NOTE THAT THE CONTROL VALUE MUST BE THE WEIGHTED ONE
    **************************************************************************/

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType>
    const typename TVariableType::Type& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType>
    void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    template<class TDataType>
    bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

#ifndef SD_APP_FORWARD_COMPATIBILITY
    template<class TAdaptorType>
    bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }
#endif

    /**************************************************************************
                            COMPARISON SUBROUTINES
    **************************************************************************/

    /// Implement relational operator for automatic arrangement in container
    inline bool operator==(const PBSplinesBasisFunction& rA) const
    {
        return (this->Id() == rA.Id()) && (this->EquationId() == rA.EquationId());
    }

    inline bool operator<(const PBSplinesBasisFunction& rA) const
    {
        if (this->Id() != rA.Id())
        {
            return this->Id() < rA.Id();
        }
        else
        {
            return this->EquationId() < rA.EquationId();
        }
    }

    /// Compare the two basis functions, in terms of the equation_id and parameter space
    inline bool IsSame(const PBSplinesBasisFunction& rA) const
    {
        if (this->EquationId() != rA.EquationId())
        {
            return false;
        }
        else
        {
            // TODO
        }

        return true;
    }

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Print information of this basis function
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PBSplinesBasisFunction (id: " << this->Id() << "), eq_id: " << this->EquationId();
    }

    /// Print data of this basis function
    virtual void PrintData(std::ostream& rOStream) const
    {
        // Print the cells
        rOStream << " Supporting cells:";
        std::size_t cnt = 0;
        for (cell_const_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            rOStream << std::endl << "  " << ++cnt << ": " << *(*it);
        }
        if (cell_end() == cell_begin())
        {
            rOStream << " none";
        }
        rOStream << std::endl;
    }

protected:

    std::size_t mId;
    std::size_t mEquationId; // this variable stores the equation id of this basis function accross patches for the case of scalar PDE.
    cell_container_t mpCells; // list of cells support this basis function

    /** A pointer to data related to this basis function. */
    DataValueContainer mData;

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Data", mData);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Data", mData);
    }

};

template<class TBasisFunctionType>
struct PBSplinesBasisFunction_InitializeValue_Helper<TBasisFunctionType, Variable<double> >
{
    static void Initialize(TBasisFunctionType& r_bf, const Variable<double>& rVariable)
    {
        if (!r_bf.Has(rVariable))
        {
            r_bf.SetValue(rVariable, 0.0);
        }
    }

    static void Initialize(TBasisFunctionType& r_bf, const Variable<double>& rVariable, typename TBasisFunctionType::Pointer p_ref_bf)
    {
        Initialize(r_bf, rVariable);
    }
};

template<class TBasisFunctionType>
struct PBSplinesBasisFunction_InitializeValue_Helper<TBasisFunctionType, Variable<array_1d<double, 3> > >
{
    static void Initialize(TBasisFunctionType& r_bf, const Variable<array_1d<double, 3> >& rVariable)
    {
        if (!r_bf.Has(rVariable))
        {
            array_1d<double, 3> zero_v;
            zero_v[0] = 0.0; zero_v[1] = 0.0; zero_v[2] = 0.0;
            r_bf.SetValue(rVariable, zero_v);
        }
    }

    static void Initialize(TBasisFunctionType& r_bf, const Variable<array_1d<double, 3> >& rVariable, typename TBasisFunctionType::Pointer p_ref_bf)
    {
        Initialize(r_bf, rVariable);
    }
};

template<class TBasisFunctionType>
struct PBSplinesBasisFunction_InitializeValue_Helper<TBasisFunctionType, Variable<Vector> >
{
    static void Initialize(TBasisFunctionType& r_bf, const Variable<Vector>& rVariable)
    {
        KRATOS_ERROR << "Not supported";
    }

    static void Initialize(TBasisFunctionType& r_bf, const Variable<Vector>& rVariable, typename TBasisFunctionType::Pointer p_ref_bf)
    {
        if (!r_bf.Has(rVariable))
        {
            Vector zero_v = p_ref_bf->GetValue(rVariable);
            for (std::size_t i = 0; i < zero_v.size(); ++i) { zero_v[i] = 0.0; }
            r_bf.SetValue(rVariable, zero_v);
        }
    }
};

/// output stream function
template<typename TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const PBSplinesBasisFunction<TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PBSPLINES_BASIS_FUNCTION_H_INCLUDED
