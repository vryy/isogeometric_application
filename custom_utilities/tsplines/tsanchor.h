//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Apr 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_ANCHOR_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_ANCHOR_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
    Represent an anchor in the T-splines mesh topology.
    In the T-splines mesh, an anchor represents a basis function. Control values and coordinates are defined at the anchor.
 */
class TsAnchor
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsAnchor);

    /// Default constructor
    TsAnchor(std::size_t Id,
             double Xi, double Eta,
             double X, double Y,
             double W)
        : mId(Id), mXi(Xi), mEta(Eta), mZeta(0.0), mX(X), mY(Y), mZ(0.0), mW(W)
    {}

    TsAnchor(std::size_t Id,
             double Xi, double Eta, double Zeta,
             double X, double Y, double Z,
             double W)
        : mId(Id), mXi(Xi), mEta(Eta), mZeta(Zeta), mX(X), mY(Y), mZ(Z), mW(W)
    {}

    /// Get the topology coordinates of the anchor
    double Xi()   const {return mXi;}
    double Eta()  const {return mEta;}
    double Zeta() const {return mZeta;}

    /// Get the coordinates
    double X() const {return mX;}
    double Y() const {return mY;}
    double Z() const {return mZ;}
    double W() const {return mW;}

    /// Get the Id of the anchor
    std::size_t Id() const {return mId;}

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "(" << Id() << ": " << Xi() << ", " << Eta() << ", " << Zeta() << ", " << W() << ")";
    }

private:

    std::size_t mId;

    double mXi;
    double mEta;
    double mZeta;

    double mX;
    double mY;
    double mZ;
    double mW; // weight of the shape function at the anchor
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsAnchor& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_ANCHOR_H_INCLUDED
