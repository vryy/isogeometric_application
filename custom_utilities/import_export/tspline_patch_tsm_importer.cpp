//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Nov 2018 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "custom_utilities/import_export/tspline_patch_tsm_importer.h"
// #include "custom_utilities/pbbsplines_basis_function.h"
// #include "custom_utilities/pbbsplines_fespace.h"

// External includes
#include "rhbuilder.h"
#include "extractor.h"

namespace Kratos
{

Patch<2>::Pointer TSplinePatchTSMImporter::ImportSingle(const std::string& filename) const
{
    RhBuilderPtr reader = TSPLINE::makePtr<RhBuilder>(filename);
    TSplinePtr tspline = reader->findTSpline();

    // create the PBSplinesFESpaceType
    // typedef PBBSplinesBasisFunction<2, Cell> PBBSplinesBasisFunctionType;
    // typedef FESpace<2> FESpaceType;
    // typedef PBBSplinesFESpace<2, PBBSplinesBasisFunctionType> PBBSplinesFESpaceType;

    // typename PBBSplinesFESpaceType::Pointer pNewFESpace = PBBSplinesFESpaceType::Create();




}

MultiPatch<2>::Pointer TSplinePatchTSMImporter::Import(const std::string& filename) const
{
    KRATOS_THROW_ERROR(std::logic_error, "TSM reader does not support multipatch", "")
}

} // namespace Kratos.

