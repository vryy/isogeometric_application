//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 May 2015 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_HB_BASIS_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_HB_BASIS_FUNCTION_H_INCLUDED

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
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/hbsplines/hb_cell.h"

namespace Kratos
{

/**
    Class represents a basis function in hierarchical B-Splines mesh
*/
class DeprecatedHBBasisFunction
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DeprecatedHBBasisFunction);

    /// Type definition
    typedef Knot<double> KnotType;
    typedef ControlPoint<double> ControlPointType;
    typedef KnotType::Pointer knot_t;

    typedef DeprecatedHBBasisFunction::Pointer bf_t;
    typedef std::vector<bf_t> bf_container_t;
    typedef bf_container_t::iterator bf_iterator;
    typedef bf_container_t::const_iterator bf_const_iterator;

    typedef HBCell<DeprecatedHBBasisFunction> CellType;
    typedef typename CellType::Pointer cell_t;
    typedef std::set<cell_t> cell_container_t;
    typedef cell_container_t::iterator cell_iterator;
    typedef cell_container_t::const_iterator cell_const_iterator;

    /// Default constructor
    DeprecatedHBBasisFunction(std::size_t Id, unsigned int Level) : mId(Id), mLevel(Level)
    {}

    /// Destructor
    ~DeprecatedHBBasisFunction()
    {}

    /**************************************************************************
                            MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Get the number of children of this basis function
    std::size_t NumberOfChildren() const {return mpChilds.size();}

    /// Add a child which support this basis function
    void AddChild(bf_t p_bf, double RefinedCoefficient)
    {
        mpChilds.push_back(p_bf);
        mRefinedCoefficients[p_bf->Id()] = RefinedCoefficient;
    }

    /// Remove the cell from the list
    void RemoveChild(bf_t p_bf)
    {
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            if(*it == p_bf)
            {
                mpChilds.erase(it);
                mRefinedCoefficients.erase(p_bf->Id());
                break;
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
    void SetLocalKnotVectors(int Dim, const std::vector<knot_t>& rpKnots)
    {
        if(Dim == 1)
        {
            mpLocalKnots1.clear();
            for(std::size_t i = 0; i < rpKnots.size(); ++i)
                mpLocalKnots1.push_back(rpKnots[i]);
        }
        else if(Dim == 2)
        {
            mpLocalKnots2.clear();
            for(std::size_t i = 0; i < rpKnots.size(); ++i)
                mpLocalKnots2.push_back(rpKnots[i]);
        }
        else if(Dim == 3)
        {
            mpLocalKnots3.clear();
            for(std::size_t i = 0; i < rpKnots.size(); ++i)
                mpLocalKnots3.push_back(rpKnots[i]);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
    }

    /**************************************************************************
                            ACCESS SUBROUTINES
    **************************************************************************/

    /// Iterators to the child of this basis function
    bf_iterator bf_begin() {return mpChilds.begin();}
    bf_iterator bf_end() {return mpChilds.end();}
    bf_const_iterator bf_const_begin() const {return mpChilds.begin();}
    bf_const_iterator bf_const_end() const {return mpChilds.end();}

    /// Iterators to the supporting cell of this basis function
    cell_iterator cell_begin() {return mpCells.begin();}
    cell_iterator cell_end() {return mpCells.end();}
    cell_const_iterator cell_const_begin() const {return mpCells.begin();}
    cell_const_iterator cell_const_end() const {return mpCells.end();}

    /// Get the Id of this basis function. Each basis function should have unique Id
    std::size_t Id() const {return mId;}

    /// Get the level id of this basis function
    unsigned int Level() const {return mLevel;}

    /// Get the underlying control point
    /// TODO move control point to the grid
    ControlPointType& GetControlPoint() {return mControlPoint;}
    const ControlPointType& GetControlPoint() const {return mControlPoint;}

    /// Get the local knot vectors
    template<class ValuesContainerType>
    void GetLocalKnots(int dim, ValuesContainerType& rKnots) const
    {
        if(dim == 1)
        {
            if(rKnots.size() != mpLocalKnots1.size())
                rKnots.resize(mpLocalKnots1.size());
            for(std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
                rKnots[i] = mpLocalKnots1[i]->Value();
        }
        else if(dim == 2)
        {
            if(rKnots.size() != mpLocalKnots2.size())
                rKnots.resize(mpLocalKnots2.size());
            for(std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
                rKnots[i] = mpLocalKnots2[i]->Value();
        }
        else if(dim == 3)
        {
            if(rKnots.size() != mpLocalKnots3.size())
                rKnots.resize(mpLocalKnots3.size());
            for(std::size_t i = 0; i < mpLocalKnots3.size(); ++i)
                rKnots[i] = mpLocalKnots3[i]->Value();
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension parameters", "")
    }

    /// Get the bounding box (=support domain) of this basis function
    void GetBoundingBox(double& Xmin, double& Xmax, double& Ymin, double& Ymax, double& Zmin, double& Zmax)
    {
        Xmin = static_cast<double>(INT_MAX);
        Xmax = -Xmin;
        Ymin = Xmin;
        Ymax = -Ymin;
        Zmin = Xmin;
        Zmax = -Zmin;
        for(cell_iterator it = cell_begin(); it != cell_end(); ++it)
        {
            if((*it)->XiMinValue() < Xmin)
                Xmin = (*it)->XiMinValue();
            if((*it)->XiMaxValue() > Xmax)
                Xmax = (*it)->XiMaxValue();
            if((*it)->EtaMinValue() < Ymin)
                Ymin = (*it)->EtaMinValue();
            if((*it)->EtaMaxValue() > Ymax)
                Ymax = (*it)->EtaMaxValue();
            if((*it)->ZetaMinValue() < Zmin)
                Zmin = (*it)->ZetaMinValue();
            if((*it)->ZetaMaxValue() > Zmax)
                Zmax = (*it)->ZetaMaxValue();
        }
    }

    /// Get the refined coefficient of a child
    double GetRefinedCoefficient(int ChildId) const
    {
        std::map<int, double>::const_iterator it = mRefinedCoefficients.find(ChildId);
        if(it != mRefinedCoefficients.end())
            return it->second;
        else
        {
            std::stringstream ss;
            ss << "The basis function " << ChildId << " is not the child of basis function " << Id();
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /// check if this bf contain this local knot vectors. Two bfs are the same if they have exactly the same local knot vector. The order of the basis function is implied.
    bool Contain(const std::vector<knot_t>& rpKnots1,
                 const std::vector<knot_t>& rpKnots2,
                 const std::vector<knot_t>& rpKnots3) const
    {
        if( mpLocalKnots1.size() != rpKnots1.size()
         || mpLocalKnots2.size() != rpKnots2.size()
         || mpLocalKnots3.size() != rpKnots3.size() )
            return false;

        for(std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
            if(mpLocalKnots1[i] != rpKnots1[i])
                return false;

        for(std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
            if(mpLocalKnots2[i] != rpKnots2[i])
                return false;

        for(std::size_t i = 0; i < mpLocalKnots3.size(); ++i)
            if(mpLocalKnots3[i] != rpKnots3[i])
                return false;

        return true;
    }

    /// return the internal reference of the knot vectors; use it with care
    const std::vector<knot_t>& LocalKnots(unsigned int i) const
    {
        if(i == 1)
            return mpLocalKnots1;
        else if(i == 2)
            return mpLocalKnots2;
        else if(i == 3)
            return mpLocalKnots3;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", "")
    }

    /**************************************************************************
                            COMPUTATION SUBROUTINES
    **************************************************************************/

    /// Compute the Bezier extraction operator of this basis function on the cell
    void ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2);
    void ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2, int p3);

    /**************************************************************************
                            COMPARISON SUBROUTINES
    **************************************************************************/

    /// Implement relational operator for automatic arrangement in container
    inline bool operator==(const DeprecatedHBBasisFunction& rA) const
    {
        return this->Id() == rA.Id();
    }

    inline bool operator<(const DeprecatedHBBasisFunction& rA) const
    {
        return this->Id() < rA.Id();
    }

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Print information of this basis function
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Bf(id:" << Id()
                 << ",(x:" << GetControlPoint().X()
                 << ",y:"  << GetControlPoint().Y()
                 << ",z:"  << GetControlPoint().Z()
                 << ",w:"  << GetControlPoint().W()
                 << "))";
    }

    /// Print data of this basis function
    void PrintData(std::ostream& rOStream) const
    {
        // Print the local knot vectors
        rOStream << "Local knot vectors:\n";
        std::cout << " 1:";
        for(std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
            std::cout << " " << mpLocalKnots1[i]->Value();
        std::cout << std::endl;

        std::cout << " 2:";
        for(std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
            std::cout << " " << mpLocalKnots2[i]->Value();
        std::cout << std::endl;

        std::cout << " 3:";
        for(std::size_t i = 0; i < mpLocalKnots3.size(); ++i)
            std::cout << " " << mpLocalKnots3[i]->Value();
        std::cout << std::endl;

        // Print the cells
        rOStream << "Supporting cells:";
        for(cell_const_iterator it = cell_const_begin(); it != cell_const_end(); ++it)
            rOStream << std::endl << *(*it);
        if(cell_const_end() == cell_const_begin())
            rOStream << " none";
        rOStream << std::endl;
        rOStream << "List of children:";
        for(bf_const_iterator it = bf_const_begin(); it != bf_const_end(); ++it)
        {
            std::map<int, double>::const_iterator it_coeff = mRefinedCoefficients.find((*it)->Id());
            rOStream << " (" << (*it)->Id() << "," << it_coeff->second << ")";
        }
        if(bf_const_end() == bf_const_begin())
            rOStream << " none";
        rOStream << std::endl;
    }

private:
    std::size_t mId;
    unsigned int mLevel;
    ControlPointType mControlPoint;
    bf_container_t mpChilds; // list of refined basis functions that constitute this basis function
    std::map<int, double> mRefinedCoefficients; // coefficient of refined basis functions
    cell_container_t mpCells; // list of cells support this basis function at the level of this basis function
    std::vector<knot_t> mpLocalKnots1;
    std::vector<knot_t> mpLocalKnots2;
    std::vector<knot_t> mpLocalKnots3;

    /// Find the span of knot in the local knot vector
    /// Remarks: it will give the based-1 index
    ///          it only works if U is non-repeated
    template<class ValuesContainerType>
    int FindSpanLocal(double Xi, const ValuesContainerType& U)
    {
        if(!U.empty())
        {
            if(Xi < U[0])
                return 0;

            if(Xi > U[U.size()-1])
                return U.size();

            for(std::size_t i = 0; i < U.size()-1; ++i)
                if(Xi >= U[i] && Xi < U[i + 1])
                    return i + 1;
        }

        return 0;
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const DeprecatedHBBasisFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_HB_BASIS_FUNCTION_H_INCLUDED

