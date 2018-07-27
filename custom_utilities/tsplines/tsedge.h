//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGE_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/tsplines/tsvertex.h"


namespace Kratos
{

/**
    Abstract class to represent an edge in Tsplines mesh topology
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
    TsEdge(const std::size_t& Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2)
    : mId(Id), mpV1(pV1), mpV2(pV2)
    {}

    /// Get and Set for edge identification
    const std::size_t& Id() const {return mId;}
    void SetId(const std::size_t& Id) {mId = Id;}

    /// check if the edge was cut by a straight ray
    /// It is noted that the value to be checked against is the index value, not the knot value (which is real)
    virtual bool IsCut(const double& index) const {return false;}

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
    virtual std::size_t Index() const {return -1;}

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Edge (" << pV1()->Id() << ", " << pV2()->Id() << "), type = " << EdgeType();
    }

private:
    std::size_t mId;
    TsVertex::Pointer mpV1, mpV2;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsEdge& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_EDGE_H_INCLUDED

