//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGES_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGES_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "custom_utilities/tsplines/tsvertex.h"


namespace Kratos
{

/**
    Represent an edge in Tsplines mesh topology
    The edge is defined as (XiL, EtaL) -> (XiU, EtaU)
 */
class TsEdge
{
public:
    /// constant definition
    static const int UNDEFINED_EDGE = 0;
    static const int VERTICAL_EDGE = 1;
    static const int HORIZONTAL_EDGE = 2;
    static const int VIRTUAL_VERTICAL_EDGE = 3;
    static const int VIRTUAL_HORIZONTAL_EDGE = 4;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsEdge);

    /// Default constructor
    TsEdge(unsigned int Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2) : mId(Id), mpV1(pV1), mpV2(pV2)
    {}

    /// Get and Set for edge identification
    unsigned int Id() const {return mId;}
    void SetId(unsigned int Id) {mId = Id;}

    /// This function is derived by TsHEdge and TsVEdge
    virtual bool IsCut(double v) const {return false;}

    /// Get associated vertices
    TsVertex::Pointer pV1() const {return mpV1;}
    TsVertex::Pointer pV2() const {return mpV2;}

    /// Get the edge type of this specific edge
    /// 0: a general edge
    /// 1: a vertical edge
    /// 2: a horizontal edge
    virtual int EdgeType() const
    {
        return UNDEFINED_EDGE;
    }

    /// Check the activeness of the edge. The edge is active if two vertices are active.
    bool IsActive() const
    {
        return pV1()->IsActive() && pV2()->IsActive();
    }

    /// Return the corresponding index of the edge
    virtual int Index() const {return -1;}

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Edge (" << pV1()->Id() << ", " << pV2()->Id() << "), type = " << EdgeType();
    }

private:
    unsigned int mId;
    TsVertex::Pointer mpV1, mpV2;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsEdge& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}



/**
    Represent a horizontal edge in Tsplines mesh
 */
class TsHEdge : public TsEdge
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsHEdge);

    /// Type definitions
    typedef TsEdge BaseType;

    /// Constructor
    TsHEdge(unsigned int Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2) : BaseType(Id, pV1, pV2)
    {
        // check if the index is valid
        if(BaseType::pV1()->Index2() != BaseType::pV2()->Index2())
            KRATOS_THROW_ERROR(std::logic_error, "The edge is not a horizonal edge", "")
    }

    /// check if the edge was cut by a vertical ray
    virtual bool IsCut(double anchor_xi_index) const
    {
        int edge_xi_index1 = BaseType::pV1()->Index1();
        int edge_xi_index2 = BaseType::pV2()->Index1();
        int edge_xi_index_min = std::min(edge_xi_index1, edge_xi_index2);
        int edge_xi_index_max = std::max(edge_xi_index1, edge_xi_index2);
        return (anchor_xi_index >= edge_xi_index_min) && (anchor_xi_index <= edge_xi_index_max);
    }

    /// Get the edge type of this specific edge
    /// 0: a general edge
    /// 1: a vertical edge
    /// 2: a horizontal edge
    virtual int EdgeType() const
    {
        return BaseType::HORIZONTAL_EDGE;
    }

    /// Return the vertical index of this edge
    virtual int Index() const
    {
        return pV1()->Index2();
    }
};


/**
    Represent a virtual horizontal edge in Tsplines mesh
 */
class TsVirtualHEdge : public TsHEdge
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsVirtualHEdge);

    /// Type definitions
    typedef TsHEdge BaseType;

    /// Constructor
    TsVirtualHEdge(unsigned int Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2) : BaseType(Id, pV1, pV2)
    {}

    /// Get the edge type of this specific edge
    virtual int EdgeType() const
    {
        return TsEdge::VIRTUAL_HORIZONTAL_EDGE;
    }
};


/**
    Represent a vertical edge in Tsplines mesh
 */
class TsVEdge : public TsEdge
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsVEdge);

    /// Type definitions
    typedef TsEdge BaseType;

    /// Constructor
    TsVEdge(unsigned int Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2) : BaseType(Id, pV1, pV2)
    {
        // check if the index is valid
        if(BaseType::pV1()->Index1() != BaseType::pV2()->Index1())
            KRATOS_THROW_ERROR(std::logic_error, "The edge is not a vertical edge", "")
    }

    /// check if the edge was cut by a horizontal ray
    virtual bool IsCut(double anchor_eta_index) const
    {
        int edge_eta_index1 = BaseType::pV1()->Index2();
        int edge_eta_index2 = BaseType::pV2()->Index2();
        int edge_eta_index_min = std::min(edge_eta_index1, edge_eta_index2);
        int edge_eta_index_max = std::max(edge_eta_index1, edge_eta_index2);
        return (anchor_eta_index >= edge_eta_index_min) && (anchor_eta_index <= edge_eta_index_max);
    }

    /// Get the edge type of this specific edge
    /// 0: a general edge
    /// 1: a vertical edge
    /// 2: a horizontal edge
    virtual int EdgeType() const
    {
        return BaseType::VERTICAL_EDGE;
    }

    /// Return the horizontal index of this edge
    virtual int Index() const
    {
        return pV1()->Index1();
    }
};


/**
    Represent a virtual vertical edge in Tsplines mesh
 */
class TsVirtualVEdge : public TsVEdge
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsVirtualVEdge);

    /// Type definitions
    typedef TsVEdge BaseType;

    /// Constructor
    TsVirtualVEdge(unsigned int Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2) : BaseType(Id, pV1, pV2)
    {}

    /// Get the edge type of this specific edge
    virtual int EdgeType() const
    {
        return TsEdge::VIRTUAL_VERTICAL_EDGE;
    }
};


}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGES_H_INCLUDED

