//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 June 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_2D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_2D_H_INCLUDED

// System includes
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/domain_manager.h"

namespace Kratos
{

/**
    This class represents a union of domains (rectangle) in the parameter space of NURBS in 2D. THis is useful to manage the refinement domain in the hierarchical NURBS mesh
 */
class DomainManager2D : public DomainManager
{

public:
    /// Type definition
    typedef DomainManager BaseType;
    typedef std::size_t key_t;
    typedef std::set<std::size_t> index_container_t;
    typedef std::map<key_t, index_container_t> map_t;

    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(DomainManager2D);

    /// Default constructor
    DomainManager2D(const std::size_t& Id) : BaseType(Id) {}

    /// Destructor
    virtual ~DomainManager2D() {}

    /// Add the Cuboid to the Cuboid set
    virtual void AddCell(const std::vector<double>& box);

    /// Check if a Cuboid if it is contained in the Cuboid set.
    virtual bool IsInside(const std::vector<double>& bounding_box) const;

    /// Export the domain to Matlab for visualization
    virtual void ExportDomain(const std::string& fn, const std::string& color, const double& distance) const;

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const;
    virtual void PrintData(std::ostream& rOStream) const;

private:

    map_t mActiveCells;
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const DomainManager2D& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_DOMAIN_MANAGER_2D_H_INCLUDED defined

