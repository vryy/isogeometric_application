//
//   Project Name:        KratosIsogeometricApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 4 Jul 2015 $
//
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#ifdef _residualbased_elimination_builder_and_solver_deactivation_existed_
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_deactivation.h"
#endif
#include "custom_strategies/builder_and_solvers/row_constraint_builder_and_solver.h"
#include "custom_python/add_strategies_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void IsogeometricApplication_AddStrategiesToPython()
{

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType, ModelPart> LinearSolverType;

#ifdef _residualbased_elimination_builder_and_solver_deactivation_existed_
    typedef ResidualBasedEliminationBuilderAndSolverDeactivation<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> ResidualBasedEliminationBuilderAndSolverDeactivationType;

    typedef RowConstraintBuilderAndSolver<ResidualBasedEliminationBuilderAndSolverDeactivationType> RowConstraintResidualBasedEliminationBuilderAndSolverDeactivationType;
    class_<RowConstraintResidualBasedEliminationBuilderAndSolverDeactivationType, bases<ResidualBasedEliminationBuilderAndSolverDeactivationType>, boost::noncopyable>
    ("RowConstraintResidualBasedEliminationBuilderAndSolverDeactivation", init<typename LinearSolverType::Pointer>())
    .def("AddConstraint", &RowConstraintResidualBasedEliminationBuilderAndSolverDeactivationType::AddConstraint)
    ;
#endif

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
