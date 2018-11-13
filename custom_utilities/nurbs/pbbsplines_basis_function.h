//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_BASIS_FUNCTION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_BASIS_FUNCTION_H_INCLUDED

// System includes
#include <cmath>

// External includes

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
#include "custom_utilities/pbsplines_basis_function.h"


namespace Kratos
{

template<int TDim>
struct PBBSplinesBasisFunction_Helper
{
    template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
    static void ComputeExtractionOperator(TVectorType& Crow, const TIArrayType& orders,
        const TKnotContainerType& local_knots, const TCellType& r_cell);

    static bool CheckBoundingBox(const std::vector<double>& bounding_box, const std::vector<std::vector<double> >& window);
};

/**
 * Abstract class for point-based B-Splines basis function.
 * Each basis function associates with a control point via CONTROL_POINT variable.
 */
template<int TDim, typename TCellType>
class PBBSplinesBasisFunction : public PBSplinesBasisFunction<TCellType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PBBSplinesBasisFunction);

    /// Type definition
    typedef PBSplinesBasisFunction<TCellType> BaseType;
    typedef typename BaseType::ControlPointType ControlPointType;
    typedef typename BaseType::CellType CellType;
    typedef typename CellType::knot_t knot_t;
    typedef typename BaseType::cell_t cell_t;
    typedef typename BaseType::const_cell_t const_cell_t;
    typedef typename BaseType::cell_container_t cell_container_t;
    typedef typename BaseType::cell_iterator cell_iterator;
    typedef typename BaseType::cell_const_iterator cell_const_iterator;

    /// Empty constructor for serialization
    PBBSplinesBasisFunction() : BaseType(), mBoundaryId(0)
    {}

    /// Constructor with Id
    PBBSplinesBasisFunction(const std::size_t& Id) : BaseType(Id), mBoundaryId(0)
    {}

    /// Destructor
    ~PBBSplinesBasisFunction()
    {
        #ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << "PBBSplinesBasisFunction" << TDim << "D " << this->Id() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    static typename PBBSplinesBasisFunction::Pointer Create(const std::size_t& Id)
    {
        return typename PBBSplinesBasisFunction::Pointer(new PBBSplinesBasisFunction(Id));
    }

    /**************************************************************************
                            ACCESS SUBROUTINES
    **************************************************************************/

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
            rKnots[i] = CellType::GetValue(mpLocalKnots[dim][i]);
    }

    /// Set the local knot vectors to this basis function
    void SetLocalKnotVectors(const int& dim, const std::vector<knot_t>& rpKnots)
    {
        mpLocalKnots[dim].clear();
        for(std::size_t i = 0; i < rpKnots.size(); ++i)
            mpLocalKnots[dim].push_back(rpKnots[i]);
    }

    /// Get the bounding box (=support domain) of this basis function
    virtual std::vector<double> GetBoundingBox() const
    {
        double Xmin = static_cast<double>(INT_MAX);
        double Xmax = -Xmin;
        double Ymin = Xmin;
        double Ymax = -Ymin;
        double Zmin = Xmin;
        double Zmax = -Zmin;
        for(cell_iterator it = BaseType::cell_begin(); it != BaseType::cell_end(); ++it)
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
    virtual double GetValueAt(const std::vector<double>& xi) const
    {
        double res;
        this->GetValueAt(res, xi);
        return res;
    }

    /// Get the value of point-based B-splines basis function
    virtual void GetValueAt(double& res, const std::vector<double>& xi) const
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
    virtual std::vector<double> GetDerivativeAt(const std::vector<double>& xi) const
    {
        std::vector<double> res;
        this->GetDerivativeAt(res, xi);
        return res;
    }

    /// Get the derivative of point-based B-splines basis function
    virtual void GetDerivativeAt(std::vector<double>& res, const std::vector<double>& xi) const
    {
        // TODO
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
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

        PBBSplinesBasisFunction_Helper<TDim>::ComputeExtractionOperator(Crow, orders, LocalKnots, *p_cell);
    }

    /**************************************************************************
                            COMPARISON SUBROUTINES
    **************************************************************************/

    /// Implement relational operator for automatic arrangement in container
    inline bool operator==(const PBBSplinesBasisFunction& rA) const
    {
        return BaseType::operator==(rA);
    }

    inline bool operator<(const PBBSplinesBasisFunction& rA) const
    {
        return BaseType::operator<(rA);
    }

    /// Compare the two basis functions, in terms of the equation_id and parameter space
    inline bool IsSame(const PBBSplinesBasisFunction& rA) const
    {
        return BaseType::IsSame(rA);
    }

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Print information of this basis function
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PBBSplinesBasisFunction" << TDim << "D (id: " << this->Id() << "), eq_id: " << this->EquationId() << ", p = (";
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
                // rOStream << " " << mpLocalKnots[dim][i]->Value();
                rOStream << " " << CellType::GetValue(mpLocalKnots[dim][i]);
            rOStream << std::endl;
        }

        // Print the cells
        rOStream << " Supporting cells:";
        std::size_t cnt = 0;
        for(cell_const_iterator it = BaseType::cell_begin(); it != BaseType::cell_end(); ++it)
            rOStream << std::endl << "  " << ++cnt << ": " << *(*it);
        if(BaseType::cell_end() == BaseType::cell_begin())
            rOStream << " none";
        rOStream << std::endl;
    }

protected:

    std::size_t mBoundaryId; // this variable stores the boundary information associated with this basis function.
                // By default, the new basis function is considerred inside of the patch.
                // User can add/remove the boundary information by using AddBoundary/RemoveBoundary
    boost::array<std::size_t, TDim> mOrders;
    cell_container_t mpCells; // list of cells support this basis function
    boost::array<std::vector<knot_t>, TDim> mpLocalKnots;

    /** A pointer to data related to this basis function. */

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        BaseType::save(rSerializer);
    }

    virtual void load(Serializer& rSerializer)
    {
        BaseType::load(rSerializer);
    }

};

/// output stream function
template<int TDim, typename TCellType>
inline std::ostream& operator <<(std::ostream& rOStream, const PBBSplinesBasisFunction<TDim, TCellType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#include "pbbsplines_basis_function.hpp"

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_BASIS_FUNCTION_H_INCLUDED

