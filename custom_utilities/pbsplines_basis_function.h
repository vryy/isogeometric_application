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
#include <ctime>
#include <cmath>
#include <climits>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

// External includes
#include <omp.h>
#include "boost/progress.hpp"
#include "boost/algorithm/string.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "containers/data_value_container.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/hbsplines/hb_cell.h"
#include "isogeometric_application/isogeometric_application.h"


namespace Kratos
{

template<int TDim>
struct PBSplinesBasisFunction_Helper
{
    template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
    static void ComputeExtractionOperator(TVectorType& Crow, const TIArrayType& orders,
        const TKnotContainerType& local_knots, const TCellType& r_cell);

    static bool CheckBoundingBox(const std::vector<double>& bounding_box, const std::vector<std::vector<double> >& window);
};

template<class TBasisFunctionType, class TVariableType>
struct PBSplinesBasisFunction_InitializeValue_Helper
{
    static void Initialize(TBasisFunctionType& r_bf, const TVariableType& rVariable)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }

    static void Initialize(TBasisFunctionType& r_bf, const TVariableType& rVariable, typename TBasisFunctionType::Pointer p_ref_bf)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }
};

/**
Class represents a point-based basis function.
Each basis function associates with a control point via CONTROL_POINT variable.
*/
template<int TDim, typename TCellType>
class PBSplinesBasisFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PBSplinesBasisFunction);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef Knot<double>::Pointer knot_t;

    typedef TCellType CellType;
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::ConstPointer const_cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename cell_container_t::const_iterator cell_const_iterator;

    /// Empty constructor for serialization
    PBSplinesBasisFunction()
    : mId(0), mEquationId(-1), mBoundaryId(0)
    {}

    /// Constructor with Id
    PBSplinesBasisFunction(const std::size_t& Id)
    : mId(Id), mEquationId(-1), mBoundaryId(0)
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
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
            if(*it == p_cell)
                return it;
        return mpCells.insert(p_cell).first;
    }

    /// Remove the cell from the list
    void RemoveCell(cell_t p_cell)
    {
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if(*it == p_cell)
            {
                mpCells.erase(it);
                break;
            }
        }
    }

    /// Remove the cell from the list
    void RemoveCell(CellType& r_cell)
    {
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if(&(*(*it)) == &r_cell)
            {
                mpCells.erase(it);
                break;
            }
        }
    }

    /// Set the local knot vectors to this basis function
    void SetLocalKnotVectors(const int& dim, const std::vector<knot_t>& rpKnots)
    {
        mpLocalKnots[dim].clear();
        for(std::size_t i = 0; i < rpKnots.size(); ++i)
            mpLocalKnots[dim].push_back(rpKnots[i]);
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
    const std::size_t& Id() const {return mId;}

    /// Get the equation Id of this basis function. Each basis function should have unique equation Id accross patches
    const std::size_t& EquationId() const {return mEquationId;}

    /// Set the equation Id for this basis function. One shall use this function only in the enumeration process
    void SetEquationId(const std::size_t& EquationId) {mEquationId = EquationId;}

    /// Get the boundary information of this basis function
    const std::size_t& BoundaryId() const {return mBoundaryId;}

    /// Add the boundary information to this basis function
    void AddBoundary(const std::size_t& BoundaryInfo) {mBoundaryId |= BoundaryInfo;}

    /// Remove the boundary information from this basis function
    void RemoveBoundary(const std::size_t& BoundaryInfo) {mBoundaryId &= (0xff - BoundaryInfo);}

    /// Function to identify the location of this basis function
    bool IsOnSide(const std::size_t& BoundaryInfo) const {return (this->BoundaryId() & BoundaryInfo) == BoundaryInfo;}

    /// Set the information in each direction
    void SetInfo(const int& dim, const std::size_t& Order) {mOrders[dim] = Order;}

    /// Get the order in specific direction
    const std::size_t& Order(const int& dim) const {return mOrders[dim];}

    /// return the internal reference of the knot vectors; use it with care
    const std::vector<knot_t>& LocalKnots(const int& dim) const {return mpLocalKnots[dim];}

    /// Get the local knot vectors
    template<class ValuesContainerType>
    void LocalKnots(const int& dim, ValuesContainerType& rKnots) const
    {
        if(rKnots.size() != mpLocalKnots[dim].size())
            rKnots.resize(mpLocalKnots[dim].size());
        for(std::size_t i = 0; i < mpLocalKnots[dim].size(); ++i)
            rKnots[i] = mpLocalKnots[dim][i]->Value();
    }

    /// Get the weight associated with the control point
    double Weight() const
    {
        return this->GetValue(CONTROL_POINT).W();
    }

    /// Get the bounding box (=support domain) of this basis function
    std::vector<double> GetBoundingBox() const
    {
        double Xmin = static_cast<double>(INT_MAX);
        double Xmax = -Xmin;
        double Ymin = Xmin;
        double Ymax = -Ymin;
        double Zmin = Xmin;
        double Zmax = -Zmin;
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if((*it)->LeftValue() < Xmin)
                Xmin = (*it)->LeftValue();
            if((*it)->RightValue() > Xmax)
                Xmax = (*it)->RightValue();
            if((*it)->DownValue() < Ymin)
                Ymin = (*it)->DownValue();
            if((*it)->UpValue() > Ymax)
                Ymax = (*it)->UpValue();
            if((*it)->BelowValue() < Zmin)
                Zmin = (*it)->BelowValue();
            if((*it)->AboveValue() > Zmax)
                Zmax = (*it)->AboveValue();
        }

        if (TDim == 1)
            return std::vector<double>{Xmin, Xmax};
        else if (TDim == 2)
            return std::vector<double>{Xmin, Xmax, Ymin, Ymax};
        else if (TDim == 3)
            return std::vector<double>{Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
    }

    /// check if this bf contain this local knot vectors. Two bfs are the same if they have exactly the same local knot vector. The order of the basis function is implied.
    // TODO change function name
    bool Contain(const std::vector<std::vector<knot_t> >& rpKnots) const
    {
        for (int dim = 0; dim < TDim; ++dim)
        {
            if (mpLocalKnots[dim].size() != rpKnots[dim].size())
                return false;

            for(std::size_t i = 0; i < mpLocalKnots[dim].size(); ++i)
                if(mpLocalKnots[dim][i] != rpKnots[dim][i])
                    return false;
        }

        return true;
    }

    /// Get the value of point-based B-splines basis function
    double GetValue(const std::vector<double>& xi) const
    {
        double res;
        this->GetValue(res, xi);
        return res;
    }

    /// Get the value of point-based B-splines basis function
    void GetValue(double& res, const std::vector<double>& xi) const
    {
        res = 1.0;
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            std::vector<double> value(1);
            std::vector<double> local_knots;
            this->LocalKnots(dim, local_knots);
            int order = this->Order(dim);
            res *= BSplineUtils::CoxDeBoor(xi[dim], 0, order, local_knots);
        }
    }

    /// Get the derivative of point-based B-splines basis function
    std::vector<double> GetDerivative(const std::vector<double>& xi) const
    {
        std::vector<double> res;
        this->GetDerivative(res, xi);
        return res;
    }

    /// Get the derivative of point-based B-splines basis function
    void GetDerivative(std::vector<double>& res, const std::vector<double>& xi) const
    {
        // TODO
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }

    /// Access the internal data container, be very careful with this function
    DataValueContainer& Data() {return mData;}
    const DataValueContainer& Data() const {return mData;}

    /**************************************************************************
                            CONTROL VALUES
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

    template<class TAdaptorType>
    bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    /**************************************************************************
                            COMPUTATION SUBROUTINES
    **************************************************************************/

    /// Compute the Bezier extraction operator of this basis function on the cell
    void ComputeExtractionOperator(Vector& Crow, const_cell_t p_cell)
    {
        std::vector<std::vector<double> > LocalKnots(TDim);
        std::vector<std::size_t> orders(TDim);
        for (int dim = 0; dim < TDim; ++dim)
        {
            orders[dim] = this->Order(dim);
            this->LocalKnots(dim, LocalKnots[dim]);
        }

        PBSplinesBasisFunction_Helper<TDim>::ComputeExtractionOperator(Crow, orders, LocalKnots, *p_cell);
    }

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
            return this->Id() < rA.Id();
        else
            return this->EquationId() < rA.EquationId();
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
        rOStream << "PBSplinesBasisFunction" << TDim << "D (id: " << this->Id() << "), eq_id: " << this->EquationId() << ", p = (";
        for (int dim = 0; dim < TDim; ++dim)
            rOStream << " " << this->Order(dim);
        rOStream << ")";
        rOStream << ", boundary info:";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BLEFT_) ) ) rOStream << " left";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BRIGHT_) ) ) rOStream << " right";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BFRONT_) ) ) rOStream << " front";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BBACK_) ) ) rOStream << " back";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BTOP_) ) ) rOStream << " top";
        if ( this->IsOnSide( BOUNDARY_FLAG(_BBOTTOM_) ) ) rOStream << " bottom";
    }

    /// Print data of this basis function
    virtual void PrintData(std::ostream& rOStream) const
    {
        // Print the local knot vectors
        rOStream << " Local knot vectors:\n";
        for (int dim = 0; dim < TDim; ++dim)
        {
            rOStream << "  " << dim+1 << ":";
            for(std::size_t i = 0; i < mpLocalKnots[dim].size(); ++i)
                rOStream << " " << mpLocalKnots[dim][i]->Value();
            rOStream << std::endl;
        }

        // Print the cells
        rOStream << " Supporting cells:";
        std::size_t cnt = 0;
        for(cell_const_iterator it = cell_begin(); it != cell_end(); ++it)
            rOStream << std::endl << "  " << ++cnt << ": " << *(*it);
        if(cell_end() == cell_begin())
            rOStream << " none";
        rOStream << std::endl;
    }

protected:

    std::size_t mId;
    std::size_t mEquationId; // this variable stores the equation id of this basis function accross patches for the case of scalar PDE.
    std::size_t mBoundaryId; // this variable stores the boundary information associated with this basis function.
                // By default, the new basis function is considerred inside of the patch.
                // User can add/remove the boundary information by using AddBoundary/RemoveBoundary
    boost::array<std::size_t, TDim> mOrders;
    cell_container_t mpCells; // list of cells support this basis function
    boost::array<std::vector<knot_t>, TDim> mpLocalKnots;

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

/// output stream function
template<int TDim, typename TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const PBSplinesBasisFunction<TDim, TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#include "pbsplines_basis_function.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PBSPLINES_BASIS_FUNCTION_H_INCLUDED

