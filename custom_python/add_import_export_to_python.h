/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Nov 30, 2017 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <sstream>

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

template<int TDim, class TExporter, class TPatchType>
void MultiPatchExporter_Export(TExporter& rDummy,
        typename TPatchType::Pointer pPatch, const std::string& filename)
{
    rDummy.template Export<TDim>(pPatch, filename);
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED  defined

