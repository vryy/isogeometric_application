//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_H_INCLUDED

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

namespace Kratos
{

template<int TDim>
struct HBSplinesBasisFunction_Helper
{
    template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
    static void ComputeExtractionOperator(TVectorType& Crow, const TIArrayType& orders,
        const TKnotContainerType& local_knots, const TCellType& r_cell);
};

/**
Class represents a basis function in hierarchical B-Splines mesh.
Each basis function associates with a control point via CONTROL_POINT variable.
*/
template<int TDim>
class HBSplinesBasisFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesBasisFunction);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef Knot<double>::Pointer knot_t;

    typedef typename HBSplinesBasisFunction<TDim>::Pointer bf_t;
    typedef std::vector<bf_t> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef HBCell<HBSplinesBasisFunction<TDim> > CellType;
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::ConstPointer const_cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename cell_container_t::const_iterator cell_const_iterator;

    /// Empty constructor for serialization
    HBSplinesBasisFunction()
    : mId(0), mEquationId(-1), mBoundaryId(0), mLevel(0)
    {}

    /// Empty constructor for serialization
    HBSplinesBasisFunction(const std::size_t& Id)
    : mId(Id), mEquationId(-1), mBoundaryId(0), mLevel(0)
    {}

    /// Default constructor
    HBSplinesBasisFunction(const std::size_t& Id, const std::size_t& Level)
    : mId(Id), mEquationId(-1), mBoundaryId(0), mLevel(Level)
    {}

    /// Destructor
    ~HBSplinesBasisFunction()
    {}

    /**************************************************************************
                            MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Get the number of children of this basis function
    std::size_t NumberOfChildren() const {return mpChilds.size();}

    /// Add a child which support this basis function
    void AddChild(bf_t p_bf, const double& RefinedCoefficient)
    {
        mpChilds.push_back(p_bf);
        mRefinedCoefficients[p_bf->Id()] = RefinedCoefficient;
    }

    /// Remove the cell from the list
    void RemoveChild(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if(*it == p_bf)
            {
                mpChilds.erase(it);
                mRefinedCoefficients.erase(p_bf->Id());
                break;
            }
        }
    }

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

    /// Iterators to the child of this basis function
    bf_iterator bf_begin() {return mpChilds.begin();}
    bf_const_iterator bf_begin() const {return mpChilds.begin();}
    bf_iterator bf_end() {return mpChilds.end();}
    bf_const_iterator bf_end() const {return mpChilds.end();}

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

    /// Get the level of this basis function
    const std::size_t& Level() const {return mLevel;}

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

    /// Get the refined coefficient of a child
    double GetRefinedCoefficient(const int& child_id) const
    {
        std::map<int, double>::const_iterator it = mRefinedCoefficients.find(child_id);
        if(it != mRefinedCoefficients.end())
            return it->second;
        else
        {
            std::stringstream ss;
            ss << "The basis function " << child_id << " is not the child of basis function " << Id();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
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

    /// Get the value of B-splines basis function
    double GetValue(const std::vector<double>& xi) const
    {
        double res;
        this->GetValue(res, xi);
        return res;
    }

    /// Get the value of B-splines basis function
    void GetValue(double& res, const std::vector<double>& xi) const
    {
        res = 1.0;

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            std::vector<double> value(1);
            std::vector<double> local_knots;
            this->LocalKnots(dim, local_knots);
            int order = this->Order(dim);
            double v = BSplineUtils::CoxDeBoor(xi[dim], 0, order, local_knots);
            res *= v;
        }
    }

    /// Get the derivative of B-splines basis function
    std::vector<double> GetDerivative(const std::vector<double>& xi) const
    {
        std::vector<double> res;
        this->GetDerivative(res, xi);
        return res;
    }

    /// Get the derivative of B-splines basis function
    void GetDerivative(std::vector<double>& res, const std::vector<double>& xi) const
    {
        // TODO
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

        HBSplinesBasisFunction_Helper<TDim>::ComputeExtractionOperator(Crow, orders, LocalKnots, *p_cell);
    }

    /// Construct the hierarchical B-Splines basis function in subspace
    typename HBSplinesBasisFunction<TDim-1>::Pointer Project(const std::size_t& dim) const
    {
        typename HBSplinesBasisFunction<TDim-1>::Pointer pNewSubBf;

        if (dim >= TDim)
            KRATOS_THROW_ERROR(std::logic_error, "The dimension is invalid", "")

        pNewSubBf = typename HBSplinesBasisFunction<TDim-1>::Pointer(new HBSplinesBasisFunction<TDim-1>(this->Id(), this->Level()));
        pNewSubBf->SetEquationId(this->EquationId());

        if (TDim == 2)
        {
            pNewSubBf->SetLocalKnotVectors(0, this->mpLocalKnots[dim]);
            pNewSubBf->SetInfo(0, this->mOrders[dim]);

            if (dim == 0)
            {
                if (this->IsOnSide(BOUNDARY_FLAG(_BLEFT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BRIGHT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));
            }
        }
        else if (TDim == 3)
        {
            if (dim == 0)
            {
                pNewSubBf->SetLocalKnotVectors(0, this->mpLocalKnots[1]);
                pNewSubBf->SetLocalKnotVectors(1, this->mpLocalKnots[2]);
                pNewSubBf->SetInfo(0, this->mOrders[1]);
                pNewSubBf->SetInfo(1, this->mOrders[2]);
                if (this->IsOnSide(BOUNDARY_FLAG(_BBOTTOM_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BTOP_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BTOP_));
                if (this->IsOnSide(BOUNDARY_FLAG(_BFRONT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BBACK_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));
            }
            else if (dim == 1)
            {
                pNewSubBf->SetLocalKnotVectors(0, this->mpLocalKnots[2]);
                pNewSubBf->SetLocalKnotVectors(1, this->mpLocalKnots[0]);
                pNewSubBf->SetInfo(0, this->mOrders[2]);
                pNewSubBf->SetInfo(1, this->mOrders[0]);
                if (this->IsOnSide(BOUNDARY_FLAG(_BLEFT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BRIGHT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BTOP_));
                if (this->IsOnSide(BOUNDARY_FLAG(_BBOTTOM_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BTOP_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));
            }
            else if (dim == 2)
            {
                pNewSubBf->SetLocalKnotVectors(0, this->mpLocalKnots[0]);
                pNewSubBf->SetLocalKnotVectors(1, this->mpLocalKnots[1]);
                pNewSubBf->SetInfo(0, this->mOrders[0]);
                pNewSubBf->SetInfo(1, this->mOrders[1]);
                if (this->IsOnSide(BOUNDARY_FLAG(_BFRONT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BBACK_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BTOP_));
                if (this->IsOnSide(BOUNDARY_FLAG(_BLEFT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
                else if (this->IsOnSide(BOUNDARY_FLAG(_BRIGHT_))) pNewSubBf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));
            }
        }

        // transfer the values
        pNewSubBf->Data() = this->Data();

        return pNewSubBf;
    }

    /**************************************************************************
                            COMPARISON SUBROUTINES
    **************************************************************************/

    /// Implement relational operator for automatic arrangement in container
    inline bool operator==(const HBSplinesBasisFunction& rA) const
    {
        return (this->Id() == rA.Id()) && (this->EquationId() == rA.EquationId());
    }

    inline bool operator<(const HBSplinesBasisFunction& rA) const
    {
        if (this->Id() != rA.Id())
            return this->Id() < rA.Id();
        else
            return this->EquationId() < rA.EquationId();
    }

    /// Compare the two basis functions, in terms of the equation_id and parameter space
    inline bool IsSame(const HBSplinesBasisFunction& rA) const
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
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HBSplinesBasisFunction" << TDim << "D (id: " << this->Id() << "), eq_id: " << this->EquationId() << ", p = (";
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
    void PrintData(std::ostream& rOStream) const
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
        rOStream << "List of children:";
        cnt = 0;
        for(bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            std::map<int, double>::const_iterator it_coeff = mRefinedCoefficients.find((*it)->Id());
            rOStream << "  " << ++cnt << ": (" << (*it)->Id() << "," << it_coeff->second << ")";
        }
        if(bf_end() == bf_begin())
            rOStream << " none";
        rOStream << std::endl;
    }

private:

    std::size_t mId;
    std::size_t mEquationId; // this variable stores the equation id of this basis function accross patches for the case of scalar PDE.
    std::size_t mLevel;
    std::size_t mBoundaryId; // this variable embeddeds the boundary information associated with this basis function. By default, the new basis function is considerred inside of the patch.
    boost::array<std::size_t, TDim> mOrders;
    bf_container_t mpChilds; // list of refined basis functions that constitute this basis function
    std::map<int, double> mRefinedCoefficients; // store the coefficient of refined basis functions
    cell_container_t mpCells; // list of cells support this basis function at the level of this basis function
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


/// Template Specialization to terminate the compilation
template<>
class HBSplinesBasisFunction<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesBasisFunction);

    /// Type definitions
    typedef Knot<double>::Pointer knot_t;
    typedef HBCell<HBSplinesBasisFunction<0> > CellType;
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::ConstPointer const_cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename cell_container_t::const_iterator cell_const_iterator;

    typedef typename HBSplinesBasisFunction<0>::Pointer bf_t;
    typedef std::vector<bf_t> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    /// Default constructor
    HBSplinesBasisFunction(const std::size_t& Id) {}

    /// Default constructor
    HBSplinesBasisFunction(const std::size_t& Id, const std::size_t& Level) {}

    /// Get the Id of the basis function
    std::size_t Id() const {return 0;}

    /// Set the equation Id for this basis function. One shall use this function only in the enumeration process
    void SetEquationId(const std::size_t& EquationId) {}

    /// Set the information in each direction
    void SetInfo(const int& dim, const std::size_t& Order) {}

    /// Set the local knot vectors to this basis function
    void SetLocalKnotVectors(const int& dim, const std::vector<knot_t>& rpKnots) {}

    /// return the internal reference of the knot vectors; use it with care
    std::vector<knot_t> LocalKnots(const int& dim) const {return std::vector<knot_t>{};}

    /// Add a cell support this basis function to the list
    cell_iterator AddCell(cell_t p_cell) {return mpCells.end();}

    /// Dummy function to return an empty data value container
    DataValueContainer Data() {return DataValueContainer();}

    /// Add the boundary information to this basis function
    void AddBoundary(const std::size_t& BoundaryInfo) {}

private:

    bf_container_t mpChilds; // list of refined basis functions that constitute this basis function
    cell_container_t mpCells; // list of cells support this basis function at the level of this basis function
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesBasisFunction<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#include "hbsplines_basis_function.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_H_INCLUDED

