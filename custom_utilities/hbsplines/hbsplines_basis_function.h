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
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hb_cell.h"
#include "isogeometric_application/isogeometric_application.h"


namespace Kratos
{

/**
Class represents a basis function in hierarchical B-Splines mesh.
Each basis function associates with a control point via CONTROL_POINT variable.
*/
template<int TDim>
class HBSplinesBasisFunction : public PBBSplinesBasisFunction<TDim, HBCell<HBSplinesBasisFunction<TDim> > >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesBasisFunction);

    /// Type definition
    typedef PBBSplinesBasisFunction<TDim, HBCell<HBSplinesBasisFunction<TDim> > > BaseType;
    typedef typename BaseType::ControlPointType ControlPointType;
    typedef typename BaseType::knot_t knot_t;

    typedef typename HBSplinesBasisFunction<TDim>::Pointer bf_t;
    typedef std::vector<bf_t> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef typename BaseType::CellType CellType; // HBCell<HBSplinesBasisFunction<TDim> >
    typedef typename CellType::Pointer cell_t;
    typedef typename CellType::ConstPointer const_cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef typename cell_container_t::iterator cell_iterator;
    typedef typename cell_container_t::const_iterator cell_const_iterator;

    /// Empty constructor for serialization
    HBSplinesBasisFunction() : BaseType(), mLevel(0)
    {}

    /// Constructor with Id
    HBSplinesBasisFunction(const std::size_t& Id)
    : BaseType(Id), mLevel(0)
    {}

    /// Constructor with Id and Level
    HBSplinesBasisFunction(const std::size_t& Id, const std::size_t& Level)
    : BaseType(Id), mLevel(Level)
    {}

    /// Destructor
    ~HBSplinesBasisFunction()
    {
        #ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "HBSplinesBasisFunction" << TDim << "D " << this->Id() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    static typename HBSplinesBasisFunction::Pointer Create(const std::size_t& Id)
    {
        return typename HBSplinesBasisFunction::Pointer(new HBSplinesBasisFunction(Id));
    }

    static typename HBSplinesBasisFunction::Pointer Create(const std::size_t& Id, const std::size_t& Level)
    {
        return typename HBSplinesBasisFunction::Pointer(new HBSplinesBasisFunction(Id, Level));
    }

    /**************************************************************************
                            MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Get the number of parent that this basis function supports
    std::size_t NumberOfParent() const {return mpParents.size();}

    /// Add a parent that this basis function supports
    void AddParent(bf_t p_bf)
    {
        mpParents.push_back(p_bf);
    }

    /// Remove the parent from the list
    void RemoveParent(bf_t p_bf)
    {
        for(bf_iterator it = bf_parent_begin(); it != bf_parent_end(); ++it)
        {
            if(*it == p_bf)
            {
                mpParents.erase(it);
                break;
            }
        }
    }

    /// Get the number of children of this basis function
    std::size_t NumberOfChildren() const {return mpChilds.size();}

    /// Add a child which support this basis function
    void AddChild(bf_t p_bf, const double& RefinedCoefficient)
    {
        mpChilds.push_back(p_bf);
        mRefinedCoefficients[p_bf->Id()] = RefinedCoefficient;
    }

    /// Remove the child from the list
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

    /**************************************************************************
                            ACCESS SUBROUTINES
    **************************************************************************/

    /// Iterators to the parent of this basis function
    bf_iterator bf_parent_begin() {return mpParents.begin();}
    bf_const_iterator bf_parent_begin() const {return mpParents.begin();}
    bf_iterator bf_parent_end() {return mpParents.end();}
    bf_const_iterator bf_parent_end() const {return mpParents.end();}

    /// Iterators to the children of this basis function
    bf_iterator bf_begin() {return mpChilds.begin();}
    bf_const_iterator bf_begin() const {return mpChilds.begin();}
    bf_iterator bf_end() {return mpChilds.end();}
    bf_const_iterator bf_end() const {return mpChilds.end();}

    /// Get and set the level of this basis function
    void SetLevel(const std::size_t& Level) {mLevel = Level;}
    const std::size_t& Level() const {return mLevel;}

    /// Get the refined coefficient of a child
    double GetRefinedCoefficient(const int& child_id) const
    {
        std::map<int, double>::const_iterator it = mRefinedCoefficients.find(child_id);
        if(it != mRefinedCoefficients.end())
            return it->second;
        else
        {
            std::stringstream ss;
            ss << "The basis function " << child_id << " is not the child of basis function " << this->Id();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /**************************************************************************
                            COMPUTATION SUBROUTINES
    **************************************************************************/

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

//    /**************************************************************************
//                            COMPARISON SUBROUTINES
//    **************************************************************************/

//    /// Implement relational operator for automatic arrangement in container
//    inline bool operator==(const HBSplinesBasisFunction& rA) const
//    {
//        return (this->Id() == rA.Id()) && (this->EquationId() == rA.EquationId());
//    }

//    inline bool operator<(const HBSplinesBasisFunction& rA) const
//    {
//        if (this->Id() != rA.Id())
//            return this->Id() < rA.Id();
//        else
//            return this->EquationId() < rA.EquationId();
//    }

//    /// Compare the two basis functions, in terms of the equation_id and parameter space
//    inline bool IsSame(const HBSplinesBasisFunction& rA) const
//    {
//        if (this->EquationId() != rA.EquationId())
//        {
//            return false;
//        }
//        else
//        {
//            // TODO
//        }

//        return true;
//    }

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
        BaseType::PrintData(rOStream);

        rOStream << " List of children:";
        std::size_t cnt = 0;
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

    std::size_t mLevel;
    bf_container_t mpParents; // list of refined basis functions that this basis function is composed from
    bf_container_t mpChilds; // list of refined basis functions that composes this basis function
    std::map<int, double> mRefinedCoefficients; // store the coefficient of refined basis functions

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

    /// Get the equation Id of this basis function. Each basis function should have unique equation Id accross patches
    std::size_t EquationId() const {return -1;}

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

    /// Remove the cell from the list
    void RemoveCell(CellType& r_cell) {}

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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_H_INCLUDED

