#include "includes/deprecated_variables.h"
#include "custom_utilities/bezier_post_utility.h"

//#define DEBUG_MULTISOLVE

namespace Kratos
{
    double BezierPostUtility::CalculateOnPoint(
        const Variable<double>& rVariable,
        double& rResult,
        Element::Pointer& pElement,
        const CoordinatesArrayType& rCoordinates
    )
    {
        Vector N;
        pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);
        
        rResult = 0.0;
        for(unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
        {
            double NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
            rResult += N( i ) * NodalValues;
        }
        return rResult;
    }

    Vector& BezierPostUtility::CalculateOnPoint(
        const Variable<Vector>& rVariable,
        Vector& rResult,
        Element::Pointer& pElement,
        const CoordinatesArrayType& rCoordinates
    )
    {
        Vector N;
        pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);
        
        for(unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
        {
            Vector& NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
            
            if(i == 0)
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
    
    array_1d<double, 3>& BezierPostUtility::CalculateOnPoint(
        const Variable<array_1d<double, 3> >& rVariable,
        array_1d<double, 3>& rResult,
        Element::Pointer& pElement,
        const CoordinatesArrayType& rCoordinates
    )
    {
        Vector N;
        pElement->GetGeometry().ShapeFunctionsValues(N, rCoordinates);
        
        rResult[0] = 0.0;
        rResult[1] = 0.0;
        rResult[2] = 0.0;
        for(unsigned int i = 0; i < pElement->GetGeometry().size(); ++i)
        {
            array_1d<double, 3> NodalValues = pElement->GetGeometry()[i].GetSolutionStepValue(rVariable);
            rResult += N( i ) * NodalValues;
        }
        
        return rResult;
    }
    
    void BezierPostUtility::TransferVariablesToNodes(LinearSolverType::Pointer& pSolver,
                                                           ModelPart& r_model_part,
                                                           const Variable<double>& rThisVariable)
    {
        ElementsArrayType& ElementsArray= r_model_part.Elements();

        // create and initialize matrix and vectors
        int NumberOfNodes = r_model_part.NumberOfNodes();
        SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
        noalias(M)= ZeroMatrix(NumberOfNodes, NumberOfNodes);
        
        SerialSparseSpaceType::VectorType g(NumberOfNodes);
        noalias(g)= ZeroVector(NumberOfNodes);
        
        SerialSparseSpaceType::VectorType b(NumberOfNodes);
        noalias(b)= ZeroVector(NumberOfNodes);

        // create a map from node Id to matrix/vector row
        std::map<unsigned int, unsigned int> MapNodeIdToVec;
        unsigned int cnt = 0;
        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
            MapNodeIdToVec[it->Id()] = cnt++;

        // create the structure for M a priori
        ConstructMatrixStructure(M, ElementsArray, MapNodeIdToVec, r_model_part.GetProcessInfo());
        
        // Transfer of GaussianVariables to Nodal Variables via L_2-Minimization
        // see Jiao + Heath "Common-refinement-based data tranfer ..."
        // International Journal for numerical methods in engineering 61 (2004) 2402--2427
        // for general description of L_2-Minimization
        
        // set up the system of equations
        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
        
        KRATOS_WATCH( number_of_threads )
        KRATOS_WATCH( element_partition )

        //create the array of lock for matrix/vector assembly
        std::vector< omp_lock_t > lock_array(NumberOfNodes);
        for(unsigned int i = 0; i < NumberOfNodes; ++i)
            omp_init_lock(&lock_array[i]);

        #pragma omp parallel for
        for(int k = 0; k < number_of_threads; ++k)
        {
            Matrix InvJ(3, 3);
            double DetJ;
            unsigned int row, col;

            typename ElementsArrayType::ptr_iterator it_begin = ElementsArray.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = ElementsArray.ptr_begin() + element_partition[k + 1];
            
            for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                if(!(*it)->GetValue(IS_INACTIVE))
                {
                    const IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    
    //                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
    //                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                    IsogeometricGeometryType& rIsogeometricGeometry = dynamic_cast<IsogeometricGeometryType&>((*it)->GetGeometry());
                    J = rIsogeometricGeometry.Jacobian0(J, (*it)->GetIntegrationMethod());
                    
                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    Matrix Ncontainer;
                    rIsogeometricGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        (*it)->GetIntegrationMethod()
                    );

                    // get the values at the integration_points
                    std::vector<double> ValuesOnIntPoint(integration_points.size());
                    (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo());

                    for(unsigned int point = 0; point< integration_points.size(); ++point)
                    {
                        MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                        double dV = DetJ * integration_points[point].Weight();
                        for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                        {
                            row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];
                            omp_set_lock(&lock_array[row]);
                            b(row) += (ValuesOnIntPoint[point]) * Ncontainer(point, prim) * dV;
                            for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                            {
                                col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                                M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                            }
                            omp_unset_lock(&lock_array[row]);
                        }
                    }
                }
                else
                {
                    // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                    for(unsigned int prim = 0 ; prim < (*it)->GetGeometry().size(); ++prim)
                    {
                        row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];
                        omp_set_lock(&lock_array[row]);
//                        b(row) += 0.0;
                        for(unsigned int sec = 0 ; sec < (*it)->GetGeometry().size(); ++sec)
                        {
                            col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                            if(col == row)
                                M(row, col) += 1.0;
//                            else
//                                M(row, col) += 0.0;
                        }
                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
        }
        
        for(unsigned int i = 0; i < NumberOfNodes; ++i)
            omp_destroy_lock(&lock_array[i]);

        // solver the system
        pSolver->Solve(M, g, b);
        
        // transfer the solution to the nodal variables
        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
        {
            unsigned int row = MapNodeIdToVec[it->Id()];
            it->GetSolutionStepValue(rThisVariable) = g(row);
        }
        std::cout << "Transfer variable to node for " << rThisVariable.Name() << " completed" << std::endl;
    }

    void BezierPostUtility::TransferVariablesToNodes(LinearSolverType::Pointer& pSolver,
                                                           ModelPart& r_model_part,
                                                           const Variable<Vector>& rThisVariable,
                                                           std::size_t ncomponents)
    {
        ElementsArrayType& ElementsArray = r_model_part.Elements();

        const unsigned int& Dim = (*(ElementsArray.ptr_begin()))->GetGeometry().WorkingSpaceDimension();

        #ifdef ENABLE_PROFILING
        //profiling variables
        double start_compute, end_compute;
        start_compute = OpenMPUtils::GetCurrentTime();
        #endif

        // create a map from node Id to matrix/vector row
        std::map<unsigned int, unsigned int> MapNodeIdToVec;
        unsigned int cnt = 0;
        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
            MapNodeIdToVec[it->Id()] = cnt++;

        // create and initialize matrix
        unsigned int NumberOfNodes = r_model_part.NumberOfNodes();
        SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
        noalias(M)= ZeroMatrix(NumberOfNodes, NumberOfNodes);
        ConstructMatrixStructure(M, ElementsArray, MapNodeIdToVec, r_model_part.GetProcessInfo());
        
        #ifdef ENABLE_PROFILING
        end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "ConstructMatrixStructure completed: " << end_compute - start_compute << " s" << std::endl;
        start_compute = end_compute;
        #endif

        // create and initialize vectors        
        SerialDenseSpaceType::MatrixType g(NumberOfNodes, ncomponents);
        noalias(g)= ZeroMatrix(NumberOfNodes, ncomponents);
        SerialDenseSpaceType::MatrixType b(NumberOfNodes, ncomponents);
        noalias(b)= ZeroMatrix(NumberOfNodes, ncomponents);

        //create a partition of the elements
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

        KRATOS_WATCH( number_of_threads )
        KRATOS_WATCH( element_partition )

        // create a lock array for parallel matrix fill
        std::vector< omp_lock_t > lock_array(NumberOfNodes);
        for(unsigned int i = 0; i < NumberOfNodes; ++i)
            omp_init_lock(&lock_array[i]);

        #pragma omp parallel for
        for(int k = 0; k < number_of_threads; ++k)
        {
            Matrix InvJ(Dim, Dim);
            double DetJ;
            unsigned int row, col;

            typename ElementsArrayType::ptr_iterator it_begin = ElementsArray.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = ElementsArray.ptr_begin() + element_partition[k + 1];

            for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                if(!(*it)->GetValue(IS_INACTIVE))
                {
                    const IntegrationPointsArrayType& integration_points
                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());

    //                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
    //                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                    IsogeometricGeometryType& rIsogeometricGeometry = dynamic_cast<IsogeometricGeometryType&>((*it)->GetGeometry());
                    J = rIsogeometricGeometry.Jacobian0(J, (*it)->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    Matrix Ncontainer;
                    rIsogeometricGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        (*it)->GetIntegrationMethod()
                    );

                    // get the values at the integration_points
                    std::vector<Vector> ValuesOnIntPoint(integration_points.size());
                    (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, r_model_part.GetProcessInfo());

                    for(unsigned int point = 0; point < integration_points.size(); ++point)
                    {
                        MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                        double dV = DetJ * integration_points[point].Weight();

                        for(unsigned int prim = 0; prim < (*it)->GetGeometry().size(); ++prim)
                        {
                            row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];

                            omp_set_lock(&lock_array[row]);

                            for(unsigned int i = 0; i < ncomponents; ++i)
                                b(row, i) += ValuesOnIntPoint[point][i] * Ncontainer(point, prim) * dV;

                            for(unsigned int sec = 0; sec < (*it)->GetGeometry().size(); ++sec)
                            {
                                col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                                M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                            }

                            omp_unset_lock(&lock_array[row]);
                        }
                    }
                }
                else
                {
                    // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                    for(unsigned int prim = 0; prim < (*it)->GetGeometry().size(); ++prim)
                    {
                        row = MapNodeIdToVec[(*it)->GetGeometry()[prim].Id()];
                            
                        omp_set_lock(&lock_array[row]);
                                
//                        for(unsigned int i = 0; i < ncomponents; ++i)
//                            b(row, i) += 0.0;
                                
                        for(unsigned int sec = 0; sec < (*it)->GetGeometry().size(); ++sec)
                        {
                            col = MapNodeIdToVec[(*it)->GetGeometry()[sec].Id()];
                            if(col == row)
                                M(row, col) += 1.0;
//                            else
//                                M(row, col) += 0.0;
                        }

                        omp_unset_lock(&lock_array[row]);
                    }
                }
            }
        }
        
        for(unsigned int i = 0; i < NumberOfNodes; ++i)
            omp_destroy_lock(&lock_array[i]);

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
        pSolver->Solve(M, g, b);

        #ifdef DEBUG_MULTISOLVE
        KRATOS_WATCH(g)
        #endif

        // transfer the solution to the nodal variables
        Vector tmp(ncomponents);
        for(ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
        {
            unsigned int r = MapNodeIdToVec[it->Id()];
            noalias(tmp) = row(g, r);
            it->GetSolutionStepValue(rThisVariable) = tmp;
        }
        std::cout << "Transfer variable to node for " << rThisVariable.Name() << " completed" << std::endl;
    }

}
