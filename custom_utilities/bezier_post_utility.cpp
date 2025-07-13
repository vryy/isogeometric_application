#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
#include "custom_utilities/isogeometric_math_utils.h"
#include "custom_utilities/bezier_post_utility.h"

// #define DEBUG_MULTISOLVE
// #define ENABLE_DEBUG
// #define ENABLE_PROFILING

namespace Kratos
{
double& BezierPostUtility_Helper<double>::CalculateOnPoint(
    const Variable<double>& rVariable,
    double& rResult,
    Element::Pointer& pElement,
    const Element::GeometryType::CoordinatesArrayType& rCoordinates)
{
    Vector N;
    pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);

    rResult = 0.0;
    for (unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
    {
        double NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
        rResult += N( i ) * NodalValues;
    }
    return rResult;
}

Vector& BezierPostUtility_Helper<Vector>::CalculateOnPoint(
    const Variable<Vector>& rVariable,
    Vector& rResult,
    Element::Pointer& pElement,
    const Element::GeometryType::CoordinatesArrayType& rCoordinates)
{
    Vector N;
    pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);

    for (unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
    {
        Vector& NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);

        if (i == 0)
        {
            rResult = N( i ) * NodalValues;
        }
        else
        {
            noalias(rResult) += N( i ) * NodalValues;
        }
    }
    return rResult;
}

array_1d<double, 3>& BezierPostUtility_Helper<array_1d<double, 3> >::CalculateOnPoint(
    const Variable<array_1d<double, 3> >& rVariable,
    array_1d<double, 3>& rResult,
    Element::Pointer& pElement,
    const Element::GeometryType::CoordinatesArrayType& rCoordinates)
{
    Vector N;
    pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);

    rResult[0] = 0.0;
    rResult[1] = 0.0;
    rResult[2] = 0.0;
    for (unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
    {
        array_1d<double, 3> NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
        rResult += N( i ) * NodalValues;
    }

    return rResult;
}

void BezierPostUtility::TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
        std::map<std::size_t, std::size_t>& node_row_id, SerialSparseSpaceType::VectorType& rValues,
        LinearSolverType::Pointer pSolver, const ModelPart& r_model_part,
        const ElementsContainerType& ElementsArray,
        const Variable<double>& rThisVariable) const
{
    active_nodes.clear();
    for ( ElementsContainerType::const_iterator it = ElementsArray.begin(); it != ElementsArray.end(); ++it )
    {
        bool is_inactive = false;
        if (mCheckActive)
        {
            if (it->IsDefined ( ACTIVE ))
            {
                is_inactive = !(it->Is( ACTIVE ));
            }
#ifndef SD_APP_FORWARD_COMPATIBILITY
            if (it->Has ( IS_INACTIVE ))
            {
                is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
            }
#endif
        }

        if ( !is_inactive )
        {
            for ( std::size_t i = 0; i < it->GetGeometry().size(); ++i )
            {
                active_nodes.insert( it->GetGeometry()[i].Id() );
            }
        }
    }

#ifdef ENABLE_DEBUG
    KRATOS_WATCH(active_nodes.size())
#endif

    // assign each node an id. That id is the row of this node in the global L2 projection matrix
    std::size_t cnt = 0;
    node_row_id.clear();
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        node_row_id[*it] = cnt++;
    }

    // create and initialize matrix and vectors
    std::size_t NumberOfNodes = active_nodes.size();
    SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
    noalias(M) = ZeroMatrix(NumberOfNodes, NumberOfNodes);

    if (rValues.size() != NumberOfNodes)
    {
        rValues.resize(NumberOfNodes, false);
    }
    noalias(rValues) = ZeroVector(NumberOfNodes);

    SerialSparseSpaceType::VectorType b(NumberOfNodes);
    noalias(b) = ZeroVector(NumberOfNodes);

    // create the structure for M a priori
    ConstructL2MatrixStructure<Element>(M, ElementsArray, node_row_id);

    // Transfer of GaussianVariables to Nodal Variables via L_2-Minimization
    // see Jiao + Heath "Common-refinement-based data tranfer ..."
    // International Journal for numerical methods in engineering 61 (2004) 2402--2427
    // for general description of L_2-Minimization

    // set up the system of equations
    //create a partition of the element array
    #ifdef ENABLE_MULTITHREAD
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
    std::vector<unsigned int> element_partition;
    OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

#ifdef ENABLE_DEBUG
    KRATOS_WATCH( number_of_threads )
    std::cout << "element_partition:";
    for (std::size_t i = 0; i < element_partition.size(); ++i)
    {
        std::cout << " " << element_partition[i];
    }
    std::cout << std::endl;
#endif

    //create the array of lock for matrix/vector assembly
    std::vector< omp_lock_t > lock_array(NumberOfNodes);
    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_init_lock(&lock_array[i]);
    }

    #pragma omp parallel for
    for (int k = 0; k < number_of_threads; ++k)
    {
        Matrix InvJ(3, 3);
        double DetJ;

        typename ElementsContainerType::const_iterator it_begin = ElementsArray.begin() + element_partition[k];
        typename ElementsContainerType::const_iterator it_end = ElementsArray.begin() + element_partition[k + 1];

        for ( ElementsContainerType::const_iterator it = it_begin; it != it_end; ++it )
        {
            bool is_inactive = false;
            if (mCheckActive)
            {
                if (it->IsDefined ( ACTIVE ))
                {
                    is_inactive = !(it->Is( ACTIVE ));
                }
#ifndef SD_APP_FORWARD_COMPATIBILITY
                if (it->Has ( IS_INACTIVE ))
                {
                    is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
                }
#endif
            }

            const auto& rElementGeometry = it->GetGeometry();

            if (!is_inactive)
            {
                const IntegrationPointsArrayType& integration_points
                    = rElementGeometry.IntegrationPoints(it->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                Matrix Ncontainer;

                const IsogeometricGeometryType* pIsogeometricGeometry = dynamic_cast<const IsogeometricGeometryType*>(&rElementGeometry);
                if (pIsogeometricGeometry != nullptr)
                {
                    J = pIsogeometricGeometry->Jacobian0(J, it->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    pIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        it->GetIntegrationMethod()
                    );
                }
                else
                {
                    typename GeometryType::MatrixType DeltaPosition(rElementGeometry.size(), 3);

                    for ( unsigned int node = 0; node < rElementGeometry.size(); ++node )
                        noalias( row( DeltaPosition, node ) ) = rElementGeometry[node].Coordinates()
                                                              - rElementGeometry[node].GetInitialPosition();

                    J = rElementGeometry.Jacobian( J, it->GetIntegrationMethod(), DeltaPosition );

                    Ncontainer = rElementGeometry.ShapeFunctionsValues( it->GetIntegrationMethod() );
                }

                // get the values at the integration_points
                std::vector<double> ValuesOnIntPoint(integration_points.size());
                it->CalculateOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo());

                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                    double dV = DetJ * integration_points[point].Weight();
                    for (unsigned int prim = 0 ; prim < rElementGeometry.size(); ++prim)
                    {
                        unsigned int row = node_row_id[rElementGeometry[prim].Id()];
                        omp_set_lock(&lock_array[row]);
                        b(row) += (ValuesOnIntPoint[point]) * Ncontainer(point, prim) * dV;
                        for (unsigned int sec = 0 ; sec < rElementGeometry.size(); ++sec)
                        {
                            unsigned int col = node_row_id[rElementGeometry[sec].Id()];
                            M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                        }
                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
            else
            {
                // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                for (unsigned int prim = 0 ; prim < rElementGeometry.size(); ++prim)
                {
                    unsigned int row = node_row_id[rElementGeometry[prim].Id()];

                    omp_set_lock(&lock_array[row]);

                    for (unsigned int sec = 0 ; sec < rElementGeometry.size(); ++sec)
                    {
                        unsigned int col = node_row_id[rElementGeometry[sec].Id()];

                        if (col == row)
                        {
                            M(row, col) += 1.0;
                        }
                    }

                    omp_unset_lock(&lock_array[row]);
                }
            }
        }
    }

    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_destroy_lock(&lock_array[i]);
    }

    // solve the system
    pSolver->Solve(M, rValues, b);
}

void BezierPostUtility::TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
        std::map<std::size_t, std::size_t>& node_row_id,
        SerialDenseSpaceType::MatrixType& rValues, LinearSolverType::Pointer pSolver,
        const ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
        const Variable<array_1d<double, 3> >& rThisVariable) const
{
#ifdef ENABLE_PROFILING
    //profiling variables
    double start_compute, end_compute;
    start_compute = OpenMPUtils::GetCurrentTime();
#endif

    active_nodes.clear();
    for ( ElementsContainerType::const_iterator it = ElementsArray.begin(); it != ElementsArray.end(); ++it )
    {
        bool is_inactive = false;
        if (mCheckActive)
        {
            if (it->IsDefined ( ACTIVE ))
            {
                is_inactive = !(it->Is( ACTIVE ));
            }
#ifndef SD_APP_FORWARD_COMPATIBILITY
            if (it->Has ( IS_INACTIVE ))
            {
                is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
            }
#endif
        }

        const auto& rElementGeometry = it->GetGeometry();

        if ( !is_inactive )
        {
            for ( std::size_t i = 0; i < rElementGeometry.size(); ++i )
            {
                active_nodes.insert( rElementGeometry[i].Id() );
            }
        }
    }

#ifdef ENABLE_DEBUG
    KRATOS_WATCH(active_nodes.size())
    // std::cout << "active_nodes:";
    // for( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    //     std::cout << " " << *it;
    // std::cout << std::endl;
#endif

    // assign each node an id. That id is the row of this node in the global L2 projection matrix
    std::size_t cnt = 0;
    node_row_id.clear();
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        node_row_id[*it] = cnt++;
    }

    // create and initialize matrix
    std::size_t NumberOfNodes = active_nodes.size();
    SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
    noalias(M) = ZeroMatrix(NumberOfNodes, NumberOfNodes);
    ConstructL2MatrixStructure<Element>(M, ElementsArray, node_row_id);

#ifdef ENABLE_PROFILING
    end_compute = OpenMPUtils::GetCurrentTime();
    std::cout << "ConstructMatrixStructure completed: " << end_compute - start_compute << " s" << std::endl;
    start_compute = end_compute;
#endif

    // create and initialize vectors
    if (rValues.size1() != NumberOfNodes || rValues.size2() != 3 )
    {
        rValues.resize(NumberOfNodes, 3);
    }
    noalias(rValues) = ZeroMatrix(NumberOfNodes, 3);
    SerialDenseSpaceType::MatrixType b(NumberOfNodes, 3);
    noalias(b) = ZeroMatrix(NumberOfNodes, 3);

    //create a partition of the elements
    #ifdef ENABLE_MULTITHREAD
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
    std::vector<unsigned int> element_partition;
    OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

#ifdef ENABLE_DEBUG
    KRATOS_WATCH( number_of_threads )
    std::cout << "element_partition:";
    for (std::size_t i = 0; i < element_partition.size(); ++i)
    {
        std::cout << " " << element_partition[i];
    }
    std::cout << std::endl;
#endif

    // create a lock array for parallel matrix fill
    std::vector< omp_lock_t > lock_array(NumberOfNodes);
    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_init_lock(&lock_array[i]);
    }

    const unsigned int Dim = (*(ElementsArray.ptr_begin()))->GetGeometry().WorkingSpaceDimension();

    #pragma omp parallel for
    for (int k = 0; k < number_of_threads; ++k)
    {
        Matrix InvJ(Dim, Dim);
        double DetJ;

        typename ElementsContainerType::const_iterator it_begin = ElementsArray.begin() + element_partition[k];
        typename ElementsContainerType::const_iterator it_end = ElementsArray.begin() + element_partition[k + 1];

        for ( ElementsContainerType::const_iterator it = it_begin; it != it_end; ++it )
        {
            bool is_inactive = false;
            if (mCheckActive)
            {
                if (it->IsDefined ( ACTIVE ))
                {
                    is_inactive = !(it->Is( ACTIVE ));
                }
#ifndef SD_APP_FORWARD_COMPATIBILITY
                if (it->Has ( IS_INACTIVE ))
                {
                    is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
                }
#endif
            }

            const auto& rElementGeometry = it->GetGeometry();

            if (!is_inactive)
            {
                const IntegrationPointsArrayType& integration_points
                    = rElementGeometry.IntegrationPoints(it->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                Matrix Ncontainer;

                const IsogeometricGeometryType* pIsogeometricGeometry = dynamic_cast<const IsogeometricGeometryType*>(&rElementGeometry);
                if (pIsogeometricGeometry != nullptr)
                {
                    J = pIsogeometricGeometry->Jacobian0(J, it->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    pIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        it->GetIntegrationMethod()
                    );
                }
                else
                {
                    typename GeometryType::MatrixType DeltaPosition(rElementGeometry.size(), 3);

                    for ( unsigned int node = 0; node < rElementGeometry.size(); ++node )
                        noalias( row( DeltaPosition, node ) ) = rElementGeometry[node].Coordinates()
                                                              - rElementGeometry[node].GetInitialPosition();

                    J = rElementGeometry.Jacobian( J, it->GetIntegrationMethod(), DeltaPosition );

                    Ncontainer = rElementGeometry.ShapeFunctionsValues( it->GetIntegrationMethod() );
                }

                // get the values at the integration_points
                std::vector<array_1d<double, 3> > ValuesOnIntPoint(integration_points.size());
                it->CalculateOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo());

                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                    double dV = DetJ * integration_points[point].Weight();

                    for (unsigned int prim = 0; prim < rElementGeometry.size(); ++prim)
                    {
                        unsigned int row = node_row_id[rElementGeometry[prim].Id()];

                        omp_set_lock(&lock_array[row]);

                        for (unsigned int i = 0; i < 3; ++i)
                        {
                            b(row, i) += ValuesOnIntPoint[point][i] * Ncontainer(point, prim) * dV;
                        }

                        for (unsigned int sec = 0; sec < rElementGeometry.size(); ++sec)
                        {
                            unsigned int col = node_row_id[rElementGeometry[sec].Id()];
                            M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                        }

                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
            else
            {
                // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                for (unsigned int prim = 0; prim < rElementGeometry.size(); ++prim)
                {
                    unsigned int row = node_row_id[rElementGeometry[prim].Id()];

                    omp_set_lock(&lock_array[row]);

//                        for(unsigned int i = 0; i < ncomponents; ++i)
//                            b(row, i) += 0.0;

                    for (unsigned int sec = 0; sec < rElementGeometry.size(); ++sec)
                    {
                        unsigned int col = node_row_id[rElementGeometry[sec].Id()];
                        if (col == row)
                        {
                            M(row, col) += 1.0;
                        }
//                            else
//                                M(row, col) += 0.0;
                    }

                    omp_unset_lock(&lock_array[row]);
                }
            }
        }
    }

    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_destroy_lock(&lock_array[i]);
    }

#ifdef ENABLE_PROFILING
    end_compute = OpenMPUtils::GetCurrentTime();
    std::cout << "Assemble the matrix completed: " << end_compute - start_compute << " s" << std::endl;
    start_compute = end_compute;
#endif

#ifdef DEBUG_MULTISOLVE
    KRATOS_WATCH(M)
    KRATOS_WATCH(b)
    KRATOS_WATCH(*pSolver)
#endif

    // solve the system
    // solver must support the multisove method
    pSolver->Solve(M, rValues, b);

#ifdef DEBUG_MULTISOLVE
    KRATOS_WATCH(g)
#endif
}

void BezierPostUtility::TransferVariablesToNodalArray(std::set<std::size_t>& active_nodes,
        std::map<std::size_t, std::size_t>& node_row_id,
        SerialDenseSpaceType::MatrixType& rValues, LinearSolverType::Pointer pSolver,
        const ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
        const Variable<Vector>& rThisVariable, std::size_t ncomponents) const
{
#ifdef ENABLE_PROFILING
    //profiling variables
    double start_compute, end_compute;
    start_compute = OpenMPUtils::GetCurrentTime();
#endif

    active_nodes.clear();
    for ( ElementsContainerType::const_iterator it = ElementsArray.begin(); it != ElementsArray.end(); ++it )
    {
        bool is_inactive = false;
        if (mCheckActive)
        {
            if (it->IsDefined ( ACTIVE ))
            {
                is_inactive = !(it->Is( ACTIVE ));
            }
#ifndef SD_APP_FORWARD_COMPATIBILITY
            if (it->Has ( IS_INACTIVE ))
            {
                is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
            }
#endif
        }

        if ( !is_inactive )
        {
            for ( std::size_t i = 0; i < it->GetGeometry().size(); ++i )
            {
                active_nodes.insert( it->GetGeometry()[i].Id() );
            }
        }
    }

#ifdef ENABLE_DEBUG
    KRATOS_WATCH(active_nodes.size())
    // std::cout << "active_nodes:";
    // for( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    //     std::cout << " " << *it;
    // std::cout << std::endl;
#endif

    // assign each node an id. That id is the row of this node in the global L2 projection matrix
    std::size_t cnt = 0;
    node_row_id.clear();
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        node_row_id[*it] = cnt++;
    }

    // create and initialize matrix
    std::size_t NumberOfNodes = active_nodes.size();
    SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
    noalias(M) = ZeroMatrix(NumberOfNodes, NumberOfNodes);
    ConstructL2MatrixStructure<Element>(M, ElementsArray, node_row_id);

#ifdef ENABLE_PROFILING
    end_compute = OpenMPUtils::GetCurrentTime();
    std::cout << "ConstructMatrixStructure completed: " << end_compute - start_compute << " s" << std::endl;
    start_compute = end_compute;
#endif

    // create and initialize vectors
    if (rValues.size1() != NumberOfNodes || rValues.size2() != ncomponents)
    {
        rValues.resize(NumberOfNodes, ncomponents, false);
    }
    noalias(rValues) = ZeroMatrix(NumberOfNodes, ncomponents);
    SerialDenseSpaceType::MatrixType b(NumberOfNodes, ncomponents);
    noalias(b) = ZeroMatrix(NumberOfNodes, ncomponents);

    //create a partition of the elements
    #ifdef ENABLE_MULTITHREAD
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
    std::vector<unsigned int> element_partition;
    OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

#ifdef ENABLE_DEBUG
    KRATOS_WATCH( number_of_threads )
    KRATOS_WATCH_STD_CON( element_partition )
#endif

    // create a lock array for parallel matrix fill
    std::vector< omp_lock_t > lock_array(NumberOfNodes);
    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_init_lock(&lock_array[i]);
    }

    const unsigned int Dim = (*(ElementsArray.ptr_begin()))->GetGeometry().WorkingSpaceDimension();

    #pragma omp parallel for
    for (int k = 0; k < number_of_threads; ++k)
    {
        Matrix InvJ(Dim, Dim);
        double DetJ;

        typename ElementsContainerType::const_iterator it_begin = ElementsArray.begin() + element_partition[k];
        typename ElementsContainerType::const_iterator it_end = ElementsArray.begin() + element_partition[k + 1];

        for ( ElementsContainerType::const_iterator it = it_begin; it != it_end; ++it )
        {
            bool is_inactive = false;
            if (mCheckActive)
            {
                if (it->IsDefined ( ACTIVE ))
                {
                    is_inactive = !(it->Is( ACTIVE ));
                }
#ifndef SD_APP_FORWARD_COMPATIBILITY
                if (it->Has ( IS_INACTIVE ))
                {
                    is_inactive = is_inactive && it->GetValue( IS_INACTIVE );
                }
#endif
            }

            const auto& rElementGeometry =  it->GetGeometry();

            if (!is_inactive)
            {
                const IntegrationPointsArrayType& integration_points
                    = rElementGeometry.IntegrationPoints(it->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                Matrix Ncontainer;

                const IsogeometricGeometryType* pIsogeometricGeometry = dynamic_cast<const IsogeometricGeometryType*>(&rElementGeometry);
                if (pIsogeometricGeometry != nullptr)
                {
                    J = pIsogeometricGeometry->Jacobian0(J, it->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    pIsogeometricGeometry->CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        it->GetIntegrationMethod()
                    );
                }
                else
                {
                    typename GeometryType::MatrixType DeltaPosition(rElementGeometry.size(), 3);

                    for ( unsigned int node = 0; node < rElementGeometry.size(); ++node )
                        noalias( row( DeltaPosition, node ) ) = rElementGeometry[node].Coordinates()
                                                              - rElementGeometry[node].GetInitialPosition();

                    J = rElementGeometry.Jacobian( J, it->GetIntegrationMethod(), DeltaPosition );

                    Ncontainer = rElementGeometry.ShapeFunctionsValues( it->GetIntegrationMethod() );
                }

                // get the values at the integration_points
                std::vector<Vector> ValuesOnIntPoint(integration_points.size());
                it->CalculateOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo());

                for (unsigned int point = 0; point < integration_points.size(); ++point)
                {
                    MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                    double dV = DetJ * integration_points[point].Weight();

                    for (unsigned int prim = 0; prim < rElementGeometry.size(); ++prim)
                    {
                        unsigned int row = node_row_id[rElementGeometry[prim].Id()];

                        omp_set_lock(&lock_array[row]);

                        for (unsigned int i = 0; i < ncomponents; ++i)
                        {
                            b(row, i) += ValuesOnIntPoint[point][i] * Ncontainer(point, prim) * dV;
                        }

                        for (unsigned int sec = 0; sec < rElementGeometry.size(); ++sec)
                        {
                            unsigned int col = node_row_id[rElementGeometry[sec].Id()];
                            M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                        }

                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
            else
            {
                // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                for (unsigned int prim = 0; prim < rElementGeometry.size(); ++prim)
                {
                    unsigned int row = node_row_id[rElementGeometry[prim].Id()];

                    omp_set_lock(&lock_array[row]);

                    for (unsigned int sec = 0; sec < rElementGeometry.size(); ++sec)
                    {
                        unsigned int col = node_row_id[rElementGeometry[sec].Id()];
                        if (col == row)
                        {
                            M(row, col) += 1.0;
                        }
                    }

                    omp_unset_lock(&lock_array[row]);
                }
            }
        }
    }

    for (unsigned int i = 0; i < NumberOfNodes; ++i)
    {
        omp_destroy_lock(&lock_array[i]);
    }

#ifdef ENABLE_PROFILING
    end_compute = OpenMPUtils::GetCurrentTime();
    std::cout << "Assemble the matrix completed: " << end_compute - start_compute << " s" << std::endl;
    start_compute = end_compute;
#endif

#ifdef DEBUG_MULTISOLVE
    KRATOS_WATCH(M)
    KRATOS_WATCH(rValues)
    KRATOS_WATCH(*pSolver)
#endif

    // solve the system
    // solver must support the multisove method
    pSolver->Solve(M, rValues, b);

#ifdef DEBUG_MULTISOLVE
    KRATOS_WATCH(rValues)
#endif
}

void BezierPostUtility::TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
        ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
        const Variable<double>& rThisVariable) const
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    SerialSparseSpaceType::VectorType g;

    this->TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part, ElementsArray, rThisVariable);

    // transfer the solution to the nodal variables
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        std::size_t row = node_row_id[*it];
        NodeType& r_node = r_model_part.GetMesh().GetNode(*it);
        r_node.GetSolutionStepValue(rThisVariable) = g(row);
    }

#ifdef ENABLE_DEBUG
    std::cout << "Transfer variable to node for " << rThisVariable.Name() << " completed" << std::endl;
#endif
}

void BezierPostUtility::TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
        ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
        const Variable<Vector>& rThisVariable, std::size_t ncomponents) const
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    SerialDenseSpaceType::MatrixType g;

    this->TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part,
                                        ElementsArray, rThisVariable, ncomponents);

    // transfer the solution to the nodal variables
    Vector tmp(ncomponents);
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        std::size_t this_row = node_row_id[*it];
        NodeType& r_node = r_model_part.GetMesh().GetNode(*it);
        noalias(tmp) = row(g, this_row);
        r_node.GetSolutionStepValue(rThisVariable) = tmp;
    }

#ifdef ENABLE_DEBUG
    std::cout << "Transfer variable to node for " << rThisVariable.Name() << " completed" << std::endl;
#endif
}

void BezierPostUtility::TransferVariablesToNodes(LinearSolverType::Pointer pSolver,
        ModelPart& r_model_part, const ElementsContainerType& ElementsArray,
        const Variable<array_1d<double, 3> >& rThisVariable) const
{
    // compute the nodal values
    std::set<std::size_t> active_nodes;
    std::map<std::size_t, std::size_t> node_row_id;
    SerialDenseSpaceType::MatrixType g;

    this->TransferVariablesToNodalArray(active_nodes, node_row_id, g, pSolver, r_model_part, ElementsArray, rThisVariable);

    // transfer the solution to the nodal variables
    array_1d<double, 3> tmp;
    for ( std::set<std::size_t>::iterator it = active_nodes.begin(); it != active_nodes.end(); ++it )
    {
        std::size_t this_row = node_row_id[*it];
        NodeType& r_node = r_model_part.GetMesh().GetNode(*it);
        noalias(tmp) = row(g, this_row);
        noalias(r_node.GetSolutionStepValue(rThisVariable)) = tmp;
    }

#ifdef ENABLE_DEBUG
    std::cout << "Transfer variable to node for " << rThisVariable.Name() << " completed" << std::endl;
#endif
}

}

#undef ENABLE_DEBUG
#undef DEBUG_MULTISOLVE
#undef ENABLE_PROFILING
#undef ENABLE_MULTITHREAD
