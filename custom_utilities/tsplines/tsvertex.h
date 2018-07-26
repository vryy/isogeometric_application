//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_VERTEX_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_VERTEX_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"

namespace Kratos
{

/**
    Represent a vertex in Tsplines mesh topology
    A vertex is determined by its topology coordinates (its indexing in the know vectors in each direction)
    This vertex definition allows for repeated knot values to be detected
    Optional values are knot values
 */
class TsVertex
{
public:
    /// Constant declaration
    static const int UNDEFINED_JOINT = -2;
    static const int BORDER_JOINT = -1;
    static const int NORMAL_JOINT = 0;
    static const int T_JOINT_LEFT = 1;
    static const int T_JOINT_RIGHT = 2;
    static const int T_JOINT_UP = 3;
    static const int T_JOINT_DOWN = 4;
    static const int L_JOINT_LEFT_DOWN = 5;
    static const int L_JOINT_LEFT_UP = 6;
    static const int L_JOINT_RIGHT_DOWN = 7;
    static const int L_JOINT_RIGHT_UP = 8;
    static const int I_JOINT_LEFT_RIGHT = 9;
    static const int I_JOINT_UP_DOWN = 10;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsVertex);

    /// Type definitions
    typedef Knot<double>::Pointer knot_t;

    /// Default constructor
    TsVertex(const std::size_t& Id, knot_t pXi, knot_t pEta)
    : mId(Id), mpXi(pXi), mpEta(pEta), mpZeta(knot_t(new Knot<double>(0.0))), mType(UNDEFINED_JOINT)
    {}

    TsVertex(const std::size_t& Id, knot_t pXi, knot_t pEta, knot_t pZeta)
    : mId(Id), mpXi(pXi), mpEta(pEta), mpZeta(pZeta), mType(UNDEFINED_JOINT)
    {}

    /// Get and Set for edge identification
    std::size_t Id() const {return mId;}
    void SetId(const std::size_t& Id) {mId = Id;}

    /// Get the index
    std::size_t Index1() const {return mpXi->Index();}
    std::size_t Index2() const {return mpEta->Index();}
    std::size_t Index3() const {return mpZeta->Index();}

    /// Get the respective knots
    knot_t pXi() const {return mpXi;}
    knot_t pEta() const {return mpEta;}
    knot_t pZeta() const {return mpZeta;}

    /// Check the activeness of the vertex
    bool IsActive() const
    {
        return (pXi()->IsActive()) && (pEta()->IsActive()) && (pZeta()->IsActive());
    }

    /// Get and Set the type of this joint
    int Type() const {return mType;}
    void SetType(const int& Type) {mType = Type;}

    /// Check if this vertex is a T-joint
    bool IsTJoint() const
    {
        return (mType == T_JOINT_LEFT)
            || (mType == T_JOINT_RIGHT)
            || (mType == T_JOINT_UP)
            || (mType == T_JOINT_DOWN);
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Vertex Id " << mId;
        rOStream << ", Index = {" << Index1() << ", " << Index2() << ", " << Index3() << "}";
    }

private:
    std::size_t mId;
    knot_t mpXi; // topology coordinates
    knot_t mpEta;
    knot_t mpZeta;
    int mType;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsVertex& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_VERTEX_H_INCLUDED

