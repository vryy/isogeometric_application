/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Nov 30, 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED

// System includes
#include <string>
#include <sstream>

// External includes

// Project includes

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

template<int TDim, class TExporter, class TPatchType, typename TVariableType>
void MultiPatchExporter_Export_Variable(TExporter& rDummy,
                                        typename TPatchType::Pointer pPatch, const TVariableType& rVariable, const std::string& filename)
{
    rDummy.template Export<TDim, TVariableType>(pPatch, rVariable, filename);
}

template<int TDim, class TExporter, class TPatchType, typename TVariableType>
void MultiPatchExporter_Export_Variable_WithComponents(TExporter& rDummy,
        typename TPatchType::Pointer pPatch, const TVariableType& rVariable, const std::string& filename, std::size_t ncomponents)
{
    rDummy.template Export<TDim, TVariableType>(pPatch, rVariable, filename, ncomponents);
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_ADD_IMPORT_EXPORT_TO_PYTHON_H_INCLUDED  defined
