//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_VIRTUAL_VEDGE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_VIRTUAL_VEDGE_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/tsplines/tsvedge.h"


namespace Kratos
{

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
    TsVirtualVEdge(const std::size_t& Id, TsVertex::Pointer pV1, TsVertex::Pointer pV2)
    : BaseType(Id, pV1, pV2)
    {}

    /// Get the edge type of this specific edge
    virtual int EdgeType() const
    {
        return TsEdge::VIRTUAL_VERTICAL_EDGE;
    }
};

}

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_VIRTUAL_VEDGE_H_INCLUDED
