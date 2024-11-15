//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-10-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_BEZIER_POST_UTILITY_H_INCLUDED )
#define  KRATOS_BEZIER_POST_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>

#ifdef ISOGEOMETRIC_USE_MPI
#include "mpi.h"
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "includes/legacy_structural_app_vars.h"
#endif
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_utilities/spatial_containers/auto_collapse_spatial_binning.h"
#else
#include "utilities/auto_collapse_spatial_binning.h"
#include "utilities/progress.h"
#endif
#include "custom_geometries/isogeometric_geometry.h"
#include "custom_utilities/isogeometric_post_utility.h"
#include "isogeometric_application_variables.h"

// #define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_MULTISOLVE
//#define DEBUG_GENERATE_MESH
// #define ENABLE_PROFILING

namespace Kratos
{
///@addtogroup IsogeometricApplication

///@{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template<typename TDataType>
struct BezierPostUtility_Helper
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Interpolation on element
    static TDataType& CalculateOnPoint(const Variable<TDataType>& rVariable,
                                       TDataType& rResult, Element::Pointer& pElement, const CoordinatesArrayType& rCoordinates)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling unimplemented function", __FUNCTION__)
    }
};

/// Short class definition.
/**
An advanced utility to export directly the FEM mesh out from isogeometric Bezier mesh. Each Bezier element will generate its own set of FEM elements. Therefore a large amount of nodes and elements may be generated.
One shall carefully use this utility for large problem. Previously, this class is named IsogeometricPostUtility.
 */
class BezierPostUtility : public IsogeometricPostUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;

    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;

    typedef typename ModelPart::NodesContainerType NodesArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;

    typedef typename ModelPart::ConditionsContainerType ConditionsArrayType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef IsogeometricGeometry<NodeType> IsogeometricGeometryType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename NodeType::DofsContainerType DofsContainerType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SerialSparseSpaceType;

    typedef UblasSpace<double, Matrix, Vector> SerialDenseSpaceType;

    typedef LinearSolver<SerialSparseSpaceType, SerialDenseSpaceType> LinearSolverType;

    typedef std::size_t IndexType;

    /// Pointer definition of BezierPostUtility
    KRATOS_CLASS_POINTER_DEFINITION(BezierPostUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BezierPostUtility()
    {
    }

    /// Destructor.
    virtual ~BezierPostUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // Synchronize post model_part with the reference model_part
    template<class TVariableType>
    void TransferNodalResults(
        const TVariableType& rThisVariable,
        ModelPart& r_model_part,
        ModelPart& r_model_part_post) const
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
#endif

        NodesArrayType& pTargetNodes = r_model_part_post.Nodes();

        ElementsContainerType& pElements = r_model_part.Elements();

        typename TVariableType::Type Results;
        CoordinatesArrayType LocalPos;
        int ElementId;

#ifdef SD_APP_FORWARD_COMPATIBILITY
        const Variable<int>& PARENT_ELEMENT_ID_var = KratosComponents<const Variable<int> >::Get("PARENT_ELEMENT_ID");
#else
        const Variable<int>& PARENT_ELEMENT_ID_var = PARENT_ELEMENT_ID;
#endif

//        #pragma omp parallel for
        //TODO: to be parallelized.
        for (typename NodesArrayType::ptr_iterator it = pTargetNodes.ptr_begin(); it != pTargetNodes.ptr_end(); ++it)
        {
            ElementId = (*it)->GetSolutionStepValue(PARENT_ELEMENT_ID_var);
            noalias(LocalPos) = (*it)->GetSolutionStepValue(LOCAL_COORDINATES);
            Results = BezierPostUtility_Helper<typename TVariableType::Type>::CalculateOnPoint(rThisVariable, Results, pElements(ElementId), LocalPos);
            (*it)->GetSolutionStepValue(rThisVariable) = Results;
        }

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer nodal point results for " << rThisVariable.Name() << " completed: " << end_compute - start_compute << " s" << std::endl;
#else
        std::cout << "Transfer nodal point results for " << rThisVariable.Name() << " completed" << std::endl;
#endif
    }

    // Synchronize post model_part with the reference model_part
    template<class TVariableType>
    void TransferIntegrationPointResults(
        const TVariableType& rThisVariable,
        ModelPart& r_model_part,
        ModelPart& r_model_part_post,
        LinearSolverType::Pointer pSolver) const
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "########################################" << std::endl;
        std::cout << "Transfer integration point results for "
                  << rThisVariable.Name() << " starts" << std::endl;
#endif

        // firstly transfer rThisVariable from integration points of reference model_part to its nodes
        TransferVariablesToNodes(pSolver, r_model_part, r_model_part.Elements(), rThisVariable);

        // secondly transfer new nodal variables results to the post model_part
        TransferNodalResults(rThisVariable, r_model_part, r_model_part_post);

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer integration point results for "
                  << rThisVariable.Name() << " completed: "
                  << end_compute - start_compute << "s" << std::endl;
        std::cout << "########################################" << std::endl;
#endif
    }

    /// Transfer the variable to nodes for model_part
    template<class TVariableType>
    void TransferVariablesToNodes(
        const TVariableType& rThisVariable,
        ModelPart& r_model_part,
        LinearSolverType::Pointer pSolver) const
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "########################################" << std::endl;
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " starts" << std::endl;
#endif

        TransferVariablesToNodes(pSolver, r_model_part, r_model_part.Elements(), rThisVariable);

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " completed: "
                  << end_compute - start_compute << "s" << std::endl;
        std::cout << "########################################" << std::endl;
#endif
    }

    /// Transfer the variable to nodes for model_part
    template<class TVariableType>
    void TransferVariablesToNodes(
        const TVariableType& rThisVariable,
        ModelPart& r_model_part,
        const ElementsContainerType& ElementsArray,
        LinearSolverType::Pointer pSolver) const
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "########################################" << std::endl;
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " starts" << std::endl;
#endif

        TransferVariablesToNodes(pSolver, r_model_part, ElementsArray, rThisVariable);

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " completed: "
                  << end_compute - start_compute << "s" << std::endl;
        std::cout << "########################################" << std::endl;
#endif
    }

    /// Compute the nodal values from the integration values
    void TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
                                       std::map<std::size_t, std::size_t>& node_row_id,
                                       SerialSparseSpaceType::VectorType& rValues, LinearSolverType::Pointer pSolver,
                                       const ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                       const Variable<double>& rThisVariable, bool check_active) const;

    /// Compute the nodal values from the integration values
    void TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
                                       std::map<std::size_t, std::size_t>& node_row_id,
                                       SerialDenseSpaceType::MatrixType& rValues, LinearSolverType::Pointer pSolver,
                                       const ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                       const Variable<array_1d<double, 3> >& rThisVariable, bool check_active) const;

    /// Compute the nodal values from the integration values
    void TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
                                       std::map<std::size_t, std::size_t>& node_row_id,
                                       SerialDenseSpaceType::MatrixType& rValues, LinearSolverType::Pointer pSolver,
                                       const ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                       const Variable<Vector>& rThisVariable, std::size_t ncomponents, bool check_active) const;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "BezierPostUtility";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BezierPostUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Transfer of rThisVariable defined on integration points to corresponding
     * nodal values. The transformation is done in a form that ensures a minimization
     * of L_2-norm error (/sum{rThisVariable- f(x)) whereas
     * f(x)= /sum{shape_func_i*rThisVariable_i}
     * @param model_part model_part on which the transfer should be done
     * @param rThisVariable Vector-Variable which should be transferred
     * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
     * Journal for numer. meth. in eng. 61 (2004) 2402--2427
     * WARNING: this may cause segmentation faults as the respective variables
     * will be created on nodal level while they are originally intended to be
     * stored on integration points!
     * @param pSolver       the solver used for solving the local system matrix
     * @param pModelPart    pointer to model_part that we wish to transfer the result from its integration points to its nodes
     * @param rThisVariable the variable need to transfer the respected values
     * @param check_active  if false the activeness of the elements will not be checked; true otherwise
     * REMARKS: this subroutine will only transfer the variables to nodes connecting with the mesh defined by ElementsArray
     */
    void TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
                                  ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                  const Variable<double>& rThisVariable,
                                  bool check_active = false) const;

    /**
     * Transfer of rThisVariable defined on integration points to corresponding
     * nodal values. The transformation is done in a form that ensures a minimization
     * of L_2-norm error (/sum{rThisVariable- f(x)) whereas
     * f(x)= /sum{shape_func_i*rThisVariable_i}
     * @param model_part model_part on which the transfer should be done
     * @param rThisVariable Vector-Variable which should be transferred
     * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
     * Journal for numer. meth. in eng. 61 (2004) 2402--2427
     * WARNING: this may cause segmentation faults as the respective variables
     * will be created on nodal level while they are originally intended to be
     * stored on integration points!
     * @param pSolver       the solver used for solving the local system matrix
     * @param pModelPart    pointer to model_part that we wish to transfer the result from its integration points to its nodes
     * @param rThisVariable the variable need to transfer the respected values
     * @param check_active  if false the activeness of the elements will not be checked; true otherwise
     * REMARKS: this subroutine will only transfer the variables to nodes connecting with the mesh defined by ElementsArray
     */
    void TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
                                  ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                  const Variable<array_1d<double, 3> >& rThisVariable,
                                  bool check_active = false) const;

    /**
     * Transfer of rThisVariable defined on integration points to corresponding
     * nodal values. The transformation is done in a form that ensures a minimization
     * of L_2-norm error (/sum{rThisVariable- f(x)) whereas
     * f(x)= /sum{shape_func_i*rThisVariable_i}
     * @param model_part model_part on which the transfer should be done
     * @param rThisVariable Vector-Variable which should be transferred
     * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
     * Journal for numer. meth. in eng. 61 (2004) 2402--2427
     * WARNING: this may cause segmentation faults as the respective variables
     * will be created on nodal level while they are originally intended to be
     * stored on integration points!
     * @param pSolver       the solver used for solving the local system matrix
     * @param pModelPart    pointer to model_part that we wish to transfer the result from its integration points to its nodes
     * @param rThisVariable the variable need to transfer the respected values
     * @param ncomponents   number of components of the nodal vector
     * @param check_active  if false the activeness of the elements will not be checked; true otherwise
     * REMARKS: this subroutine will only transfer the variables to nodes connecting with the mesh defined by ElementsArray
     */
    void TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
                                  ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
                                  const Variable<Vector>& rThisVariable,
                                  std::size_t ncomponents = 6,
                                  bool check_active = false) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BezierPostUtility& operator=(BezierPostUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    BezierPostUtility(BezierPostUtility const& rOther)
    {
    }

    ///@}

}; // Class BezierPostUtility

///@}

template<>
struct BezierPostUtility_Helper<double>
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Interpolation on element
    static double& CalculateOnPoint(const Variable<double>& rVariable,
                                    double& rResult, Element::Pointer& pElement, const CoordinatesArrayType& rCoordinates);
};

template<>
struct BezierPostUtility_Helper<Vector>
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Interpolation on element
    static Vector& CalculateOnPoint(const Variable<Vector>& rVariable,
                                    Vector& rResult, Element::Pointer& pElement, const CoordinatesArrayType& rCoordinates);
};

template<>
struct BezierPostUtility_Helper<array_1d<double, 3> >
{
    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Interpolation on element
    static array_1d<double, 3>& CalculateOnPoint(const Variable<array_1d<double, 3> >& rVariable,
            array_1d<double, 3>& rResult, Element::Pointer& pElement, const CoordinatesArrayType& rCoordinates);
};

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BezierPostUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const BezierPostUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
#undef DEBUG_MULTISOLVE
#undef DEBUG_GENERATE_MESH
#undef ENABLE_PROFILING

#endif
