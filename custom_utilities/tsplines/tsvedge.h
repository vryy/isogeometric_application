//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_VEDGE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_VEDGE_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/tsplines/tsedge.h"


namespace Kratos
{

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
    TsVEdge(const std::size_t& Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2)
    : BaseType(Id, pV1, pV2)
    {
        // check if the index is valid
        if(BaseType::pV1()->Index1() != BaseType::pV2()->Index1())
            KRATOS_THROW_ERROR(std::logic_error, "The edge is not a vertical edge", "")
    }

    /// check if the edge was cut by a horizontal ray
    virtual bool IsCut(const double& anchor_eta_index) const
    {
        std::size_t edge_eta_index1 = BaseType::pV1()->Index2();
        std::size_t edge_eta_index2 = BaseType::pV2()->Index2();
        std::size_t edge_eta_index_min = std::min(edge_eta_index1, edge_eta_index2);
        std::size_t edge_eta_index_max = std::max(edge_eta_index1, edge_eta_index2);
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
    virtual std::size_t Index() const
    {
        return pV1()->Index1();
    }
};

}

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_VEDGE_H_INCLUDED
