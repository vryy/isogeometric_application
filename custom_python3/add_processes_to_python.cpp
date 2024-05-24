/*
LICENSE: see isogeometric_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 21 Jan 2015 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"
#include "processes/process.h"
#include "custom_python3/add_processes_to_python.h"

#ifdef ISOGEOMETRIC_USE_PARMETIS
#include "custom_processes/isogeometric_partitioning_process.h"
#endif

namespace Kratos
{

namespace Python
{

void IsogeometricApplication_AddProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

#ifdef ISOGEOMETRIC_USE_PARMETIS
    class_<IsogeometricPartitioningProcess, Process>
    (m, "IsogeometricPartitioningProcess")
    .def(init<ModelPart&, IO&, unsigned int>())
    ;
#endif
}

} // namespace Python.

} // Namespace Kratos
