//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_IMPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_IMPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/multipatch.h"

namespace Kratos
{

/**
Abstract class to export the MultiPatch Geometry to various visualization framework
 */
template<int TDim>
class MultiPatchImporter
{
public:
    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(MultiPatchImporter);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiPatchImporter() {}

    /// Destructor
    virtual ~MultiPatchImporter() {}

    /// Import a single patch from file
    virtual typename Patch<TDim>::Pointer ImportSingle(const std::string& filename) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Import a multipatch from file
    virtual typename MultiPatch<TDim>::Pointer Import(const std::string& filename) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchImporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    std::size_t mAccuracy; // the number of digits after comma

}; // end class MultiPatchImporter

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchImporter<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_IMPORTER_H_INCLUDED defined

