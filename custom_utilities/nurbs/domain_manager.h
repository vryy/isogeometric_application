//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 June 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED

// System includes
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
This class represents a union of domains (rectangle/cuboid) in the parameter space of NURBS in 2D. THis is useful to manage the refinement domain in the hierarchical NURBS mesh
 */
class DomainManager
{

public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DomainManager);

    /// Type definition
    typedef std::set<double> coords_container_t;

    /// Default constructor
    DomainManager(std::size_t Id) : mId(Id) {}

    /// Destructor
    virtual ~DomainManager() {}

    /// Get the Id of this domain manager
    std::size_t Id() const {return mId;}

    /// Set the id for this domain manager
    void SetId(std::size_t Id) {mId = Id;}

    /// Fill the internal array of X & Y-coordinates. It must be done before cells are added to this container.
    virtual void AddXcoord(double X) {mXcoords.insert(X);}
    virtual void AddYcoord(double Y) {mYcoords.insert(Y);}
    virtual void AddZcoord(double Z) {mZcoords.insert(Z);}

    /// Add the cell to the set
    virtual void AddCell(const std::vector<double>& box)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Check if a cell if it is contained in the set.
    virtual bool IsInside(const std::vector<double>& bounding_box) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Export the domain to Matlab for visualization
    virtual void ExportDomain(const std::string& fn, const std::string& color, double distance) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const {}
    virtual void PrintData(std::ostream& rOStream) const {}

protected:

    coords_container_t mXcoords;
    coords_container_t mYcoords;
    coords_container_t mZcoords;

private:

    std::size_t mId;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const DomainManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_H_INCLUDED defined

