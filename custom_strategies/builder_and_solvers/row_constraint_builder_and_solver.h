/*
see isogeometric_application/LICENSE.txt
*/

/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 21 May 2018 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_ROW_CONSTRAINT_BUILDER_AND_SOLVER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_ROW_CONSTRAINT_BUILDER_AND_SOLVER_H_INCLUDED


/* System includes */
#include <set>

/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * A special BuilderAndSolver that maps a specified row of the stiffness matrix to other row.
 */
template<class TBuilderAndSolverType>
class RowConstraintBuilderAndSolver : public TBuilderAndSolverType
{
public:
    /**@name Const Definitions */
    /*@{ */

    /*@} */

    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( RowConstraintBuilderAndSolver );

    typedef TBuilderAndSolverType BaseType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TLinearSolverType TLinearSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef typename TSparseSpaceType::IndexType IndexType;

    typedef typename TSparseSpaceType::SizeType SizeType;

    typedef std::size_t KeyType; // For Dof->GetVariable().Key()

/*    typedef typename TSparseSpaceType::CommType CommType;*/

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /// Default Constructor.
    RowConstraintBuilderAndSolver(typename TLinearSolverType::Pointer pNewLinearSystemSolver)
    : BaseType(pNewLinearSystemSolver)
    {
    }

    /// Destructor
    virtual ~RowConstraintBuilderAndSolver()
    {
    }


    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    //**************************************************************************
    //**************************************************************************
    void AddConstraint(const std::size_t& row1, const std::size_t& row2)
    {
        mRowConstraints[row1] = row2;
    }

    //**************************************************************************
    //**************************************************************************
    virtual void Build(typename TSchemeType::Pointer pScheme,
                       ModelPart& r_model_part,
                       TSystemMatrixType& A,
                       TSystemVectorType& b)
    {
        KRATOS_TRY

        BaseType::Build(pScheme, r_model_part, A, b);

        for (std::map<std::size_t, std::size_t>::iterator it = mRowConstraints.begin(); it != mRowConstraints.end(); ++it)
        {
            const std::size_t& row1 = it->first;
            const std::size_t& row2 = it->second;

            noalias( row(A, row1) ) = row(A, row2);
            noalias( column(A, row1) ) = column(A, row2);
        }

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Access */
    /*@{ */

    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:

    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    // The variable to construct the row constraint
    std::map<std::size_t, std::size_t> mRowConstraints;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class RowConstraintBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_ISOGEOMETRIC_APPLICATION_ROW_CONSTRAINT_BUILDER_AND_SOLVER_H_INCLUDED defined */

