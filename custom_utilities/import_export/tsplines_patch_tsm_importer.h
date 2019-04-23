//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Nov 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_PATCH_TSM_IMPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_PATCH_TSM_IMPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/import_export/multipatch_importer.h"


namespace Kratos
{

class TSplinesPatchTSMImporter : public MultiPatchImporter<2>
{
public:
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(TSplinesPatchTSMImporter);

    virtual Patch<2>::Pointer ImportSingle(const std::string& filename) const;

    virtual MultiPatch<2>::Pointer Import(const std::string& filename) const;

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TSplinesPatchTSMImporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TSplinesPatchTSMImporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TSPLINES_PATCH_TSM_IMPORTER_H_INCLUDED defined

