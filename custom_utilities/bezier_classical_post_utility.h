//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-10-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_BEZIER_CLASSICAL_POST_UTILITY_H_INCLUDED )
#define  KRATOS_BEZIER_CLASSICAL_POST_UTILITY_H_INCLUDED

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
#include "includes/deprecated_variables.h"
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "includes/legacy_structural_app_vars.h"
#endif
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "custom_utilities/spatial_containers/auto_collapse_spatial_binning.h"
#include "layer_application/custom_utilities/containers/vector_map.h"
#else
#include "utilities/auto_collapse_spatial_binning.h"
#include "containers/vector_map.h"
#endif
#include "custom_utilities/iga_define.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "custom_utilities/isogeometric_utility.h"
#include "custom_utilities/isogeometric_post_utility.h"
#include "isogeometric_application_variables.h"

//#define DEBUG_LEVEL1
//#define DEBUG_LEVEL2
//#define DEBUG_MULTISOLVE
//#define DEBUG_GENERATE_MESH
#define ENABLE_PROFILING

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

template<class T> void AddToModelPart(ModelPart& rModelPart, typename T::Pointer pE);

template<> void AddToModelPart<Element>(ModelPart& rModelPart, typename Element::Pointer pE)
{
    rModelPart.AddElement(pE);
}

template<> void AddToModelPart<Condition>(ModelPart& rModelPart, typename Condition::Pointer pC)
{
    rModelPart.AddCondition(pC);
}

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
A simple utility to export directly the FEM mesh out from isogeometric Bezier mesh. Each Bezier element will generate its own set of FEM elements. Therefore a large amount of nodes and elements may be generated.
One shall carefully use this utility for large problem. Previously, this class is named IsogeometricClassicalPostUtility.
 */
class BezierClassicalPostUtility : public IsogeometricPostUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;

    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;

    typedef typename ModelPart::NodesContainerType NodesArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsArrayType;

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

    /// Pointer definition of BezierClassicalPostUtility
    KRATOS_CLASS_POINTER_DEFINITION(BezierClassicalPostUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BezierClassicalPostUtility(ModelPart& r_model_part)
        : mr_model_part(r_model_part)
    {
    }

    /// Destructor.
    virtual ~BezierClassicalPostUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Generate the post model_part from reference model_part
    /// Deprecated
    void GenerateModelPart(ModelPart& rModelPartPost, PostElementType postElementType)
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
#endif

#ifdef DEBUG_LEVEL1
        std::cout << typeid(*this).name() << "::GenerateModelPart" << std::endl;
#endif

        ElementsArrayType& pElements = mr_model_part.Elements();

#ifdef DEBUG_LEVEL1
        std::cout << "Retrieved pElements" << std::endl;
#endif

        std::string NodeKey = std::string("Node");

        //select the correct post element type
        std::string element_name;
        if (postElementType == _TRIANGLE_)
        {
            element_name = std::string("KinematicLinear2D3N");
        }
        else if (postElementType == _QUADRILATERAL_)
        {
            element_name = std::string("KinematicLinear2D4N");
        }
        else if (postElementType == _TETRAHEDRA_)
        {
            element_name = std::string("KinematicLinear3D4N");
        }
        else if (postElementType == _HEXAHEDRA_)
        {
            element_name = std::string("KinematicLinear3D8N");
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "This element type is not supported for isogeometric post-processing", __FUNCTION__);
        }

        if (!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
            KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
        }

        Element const& rCloneElement = KratosComponents<Element>::Get(element_name);

        IndexType NodeCounter = 0;
        IndexType ElementCounter = 0;
        Kratos::progress_display show_progress( pElements.size() );
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            bool is_active = true;
            if (it->IsDefined ( ACTIVE ))
            {
                is_active = it->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
            }
#else
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
            }
            else
            {
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = !it->GetValue( IS_INACTIVE );
                }
            }
#endif

            if (!is_active)
            {
//                std::cout << "Element " << it->Id() << " is inactive" << std::endl;
                continue;
            }

            int Dim = it->GetGeometry().WorkingSpaceDimension();
            IndexType NodeCounter_old = NodeCounter;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(Dim)
#endif

            //get the properties
            Properties::Pointer pProperties = rModelPartPost.pGetProperties(it->pGetProperties()->Id());

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(*pProperties)
#endif

            // generate list of nodes
            if (Dim == 1)
            {
                // TODO
            }
            else if (Dim == 2)
            {
                IndexType NumDivision1 = static_cast<IndexType>( it->GetValue(NUM_DIVISION_1) );
                IndexType NumDivision2 = static_cast<IndexType>( it->GetValue(NUM_DIVISION_2) );
                IndexType i, j;
                CoordinatesArrayType p_ref;
                CoordinatesArrayType p;

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(NumDivision1)
                KRATOS_WATCH(NumDivision2)
                std::cout << "Generating Nodes..." << std::endl;
#endif

                // create and add nodes
                p_ref[2] = 0.0;
                for (i = 0; i <= NumDivision1; ++i)
                {
                    p_ref[0] = ((double) i) / NumDivision1;
                    for (j = 0; j <= NumDivision2; ++j)
                    {
                        p_ref[1] = ((double) j) / NumDivision2;

                        p = GlobalCoordinates(it->GetGeometry(), p, p_ref);

                        NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                        pNewNode->SetId(++NodeCounter);

#ifdef DEBUG_GENERATE_MESH
//                        if(NodeCounter == 585 || NodeCounter == 588 || NodeCounter == 589)
                        if (NodeCounter)
                        {
                            std::cout << "Node " << NodeCounter << " p_ref: " << p_ref << ", p: " << p << std::endl;
                        }
#endif

                        // Giving model part's variables list to the node
                        pNewNode->SetSolutionStepVariablesList(&rModelPartPost.GetNodalSolutionStepVariablesList());

                        //set buffer size
                        pNewNode->SetBufferSize(rModelPartPost.GetBufferSize());

                        rModelPartPost.AddNode(pNewNode);

                        mNodeToLocalCoordinates(pNewNode->Id()) = p_ref;
                        mNodeToElement(pNewNode->Id()) = it->Id();
                    }
                }

                //for correct mapping to element, the repetitive node is allowed.
//                rModelPartPost.Nodes().Unique();

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(rModelPartPost.Nodes().size())
                std::cout << "Generating Elements..." << std::endl;
#endif

                // create and add element
                // Element::NodesArrayType temp_element_nodes;
                // for(i = 0; i < NumDivision1; ++i)
                // {
                //     for(j = 0; j < NumDivision2; ++j)
                //     {
                //         int Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
                //         int Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
                //         int Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
                //         int Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

                //         if(postElementType == _TRIANGLE_)
                //         {
                //             // TODO: check if jacobian checking is necessary
                //             temp_element_nodes.clear();
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node1, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node2, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node4, NodeKey).base()));

                //             Element::Pointer NewElement1 = rCloneElement.Create(++ElementCounter, temp_element_nodes, pProperties);
                //             rModelPartPost.AddElement(NewElement1);
                //             mOldToNewElements[it->Id()].insert(ElementCounter);

                //             temp_element_nodes.clear();
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node1, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node4, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node3, NodeKey).base()));

                //             Element::Pointer NewElement2 = rCloneElement.Create(++ElementCounter, temp_element_nodes, pProperties);
                //             rModelPartPost.AddElement(NewElement2);
                //             mOldToNewElements[it->Id()].insert(ElementCounter);
                //         }
                //         else if(postElementType == _QUADRILATERAL_)
                //         {
                //             // TODO: check if jacobian checking is necessary
                //             temp_element_nodes.clear();
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node1, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node2, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node4, NodeKey).base()));
                //             temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node3, NodeKey).base()));

                //             Element::Pointer NewElement = rCloneElement.Create(++ElementCounter, temp_element_nodes, pProperties);
                //             rModelPartPost.AddElement(NewElement);
                //             mOldToNewElements[it->Id()].insert(ElementCounter);
                //         }
                //     }
                // }

                std::vector<std::vector<IndexType> > connectivities;

                for (i = 0; i < NumDivision1; ++i)
                {
                    for (j = 0; j < NumDivision2; ++j)
                    {
                        IndexType Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
                        IndexType Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
                        IndexType Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
                        IndexType Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

                        if (postElementType == _TRIANGLE_)
                        {
                            connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4});
                            connectivities.push_back(std::vector<IndexType> {Node1, Node4, Node3});
                        }
                        else if (postElementType == _QUADRILATERAL_)
                        {
                            connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4, Node3});
                        }
                    }
                }

                ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                                                     connectivities, rModelPartPost, rCloneElement, ElementCounter, pProperties, NodeKey);

                for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
                {
                    rModelPartPost.AddElement(*it2);
                    mOldToNewElements[it->Id()].insert((*it2)->Id());
                }

                rModelPartPost.Elements().Unique();

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(rModelPartPost.Elements().size())
#endif
            }
            else if (Dim == 3)
            {
                IndexType NumDivision1 = static_cast<IndexType>( it->GetValue(NUM_DIVISION_1) );
                IndexType NumDivision2 = static_cast<IndexType>( it->GetValue(NUM_DIVISION_2) );
                IndexType NumDivision3 = static_cast<IndexType>( it->GetValue(NUM_DIVISION_3) );
                IndexType i, j, k;
                CoordinatesArrayType p_ref;
                CoordinatesArrayType p;

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(it->Id())
                KRATOS_WATCH(NumDivision1)
                KRATOS_WATCH(NumDivision2)
                KRATOS_WATCH(NumDivision3)
                std::cout << "Generating Nodes..." << std::endl;
#endif

                // create and add nodes
                for (i = 0; i <= NumDivision1; ++i)
                {
                    p_ref[0] = ((double) i) / NumDivision1;
                    for (j = 0; j <= NumDivision2; ++j)
                    {
                        p_ref[1] = ((double) j) / NumDivision2;
                        for (k = 0; k <= NumDivision3; ++k)
                        {
                            p_ref[2] = ((double) k) / NumDivision3;

                            p = GlobalCoordinates(it->GetGeometry(), p, p_ref);

                            NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                            pNewNode->SetId(++NodeCounter);

#ifdef DEBUG_GENERATE_MESH
                            if (NodeCounter)
                            {
                                std::cout << "Node " << NodeCounter << " p_ref: " << p_ref << ", p: " << p << std::endl;
                            }
#endif

                            // Giving model part's variables list to the node
                            pNewNode->SetSolutionStepVariablesList(&rModelPartPost.GetNodalSolutionStepVariablesList());

                            //set buffer size
                            pNewNode->SetBufferSize(rModelPartPost.GetBufferSize());

                            rModelPartPost.AddNode(pNewNode);

                            mNodeToLocalCoordinates(pNewNode->Id()) = p_ref;
                            mNodeToElement(pNewNode->Id()) = it->Id();
                        }
                    }
                }

                //for correct mapping to element, the repetitive node is allowed.
//                rModelPartPost.Nodes().Unique();

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(rModelPartPost.Nodes().size())
                std::cout << "Generating Elements..." << std::endl;
#endif

                // create and add element
                // Element::NodesArrayType temp_element_nodes;
                // for(i = 0; i < NumDivision1; ++i)
                // {
                //     for(j = 0; j < NumDivision2; ++j)
                //     {
                //         for(k = 0; k < NumDivision3; ++k)
                //         {
                //             int Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                //             int Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                //             int Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                //             int Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                //             int Node5 = Node1 + 1;
                //             int Node6 = Node2 + 1;
                //             int Node7 = Node3 + 1;
                //             int Node8 = Node4 + 1;

                //             if(postElementType == _TETRAHEDRA_)
                //             {
                //                 // TODO: check if jacobian checking is necessary
                //             }
                //             else if(postElementType == _HEXAHEDRA_)
                //             {
                //                 // TODO: check if jacobian checking is necessary
                //                 temp_element_nodes.clear();
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node1, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node2, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node4, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node3, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node5, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node6, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node8, NodeKey).base()));
                //                 temp_element_nodes.push_back(*(FindKey(rModelPartPost.Nodes(), Node7, NodeKey).base()));

                //                 Element::Pointer NewElement = rCloneElement.Create(++ElementCounter, temp_element_nodes, pProperties);
                //                 rModelPartPost.AddElement(NewElement);
                //                 mOldToNewElements[it->Id()].insert(ElementCounter);
                //             }
                //         }
                //     }
                // }

                std::vector<std::vector<IndexType> > connectivities;

                for (i = 0; i < NumDivision1; ++i)
                {
                    for (j = 0; j < NumDivision2; ++j)
                    {
                        for (k = 0; k < NumDivision3; ++k)
                        {
                            IndexType Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                            IndexType Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                            IndexType Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                            IndexType Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                            IndexType Node5 = Node1 + 1;
                            IndexType Node6 = Node2 + 1;
                            IndexType Node7 = Node3 + 1;
                            IndexType Node8 = Node4 + 1;

                            if (postElementType == _TETRAHEDRA_)
                            {
                                // TODO: check if creating Tetrahedra is correct
                                connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4});
                                connectivities.push_back(std::vector<IndexType> {Node1, Node4, Node3});
                            }
                            else if (postElementType == _HEXAHEDRA_)
                            {
                                connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4, Node3, Node5, Node6, Node8, Node7});
                            }
                        }
                    }
                }

                ElementsArrayType pNewElements = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, Element, ElementsArrayType>(
                                                     connectivities, rModelPartPost, rCloneElement, ElementCounter, pProperties, NodeKey);

                for (typename ElementsArrayType::ptr_iterator it2 = pNewElements.ptr_begin(); it2 != pNewElements.ptr_end(); ++it2)
                {
                    rModelPartPost.AddElement(*it2);
                    mOldToNewElements[it->Id()].insert((*it2)->Id());
                }

                rModelPartPost.Elements().Unique();

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(rModelPartPost.Elements().size())
#endif
            }
            ++show_progress;
        }

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GeneratePostModelPart completed: " << (end_compute - start_compute) << " s" << std::endl;
#else
        std::cout << "GeneratePostModelPart completed" << std::endl;
#endif
        std::cout << NodeCounter << " nodes and " << ElementCounter << " elements are created" << std::endl;
    }

    /// Generate the post model_part from reference model_part
    /// this is the improved version of GenerateModelPart
    /// which uses template function to generate post Elements for both Element and Condition
    void GenerateModelPart2(ModelPart& rModelPartPost, const bool& generate_for_condition)
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
#endif

#ifdef DEBUG_LEVEL1
        std::cout << typeid(*this).name() << "::GenerateModelPart" << std::endl;
#endif

        ElementsArrayType& pElements = mr_model_part.Elements();
        ConditionsArrayType& pConditions = mr_model_part.Conditions();

        std::string NodeKey = std::string("Node");

        IndexType NodeCounter = 0;
        IndexType ElementCounter = 0;
        Kratos::progress_display show_progress( pElements.size() );
        std::vector<std::size_t> dummy_ids;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            // This is wrong, we will not skill the IS_INACTIVE elements
            // TODO: to be deleted
//            if(it->GetValue( IS_INACTIVE ))
//            {
////                std::cout << "Element " << it->Id() << " is inactive" << std::endl;
//                ++show_progress;
//                continue;
//            }
            if (it->pGetGeometry() == 0)
                KRATOS_THROW_ERROR(std::logic_error, "Error: geometry is NULL at element", it->Id())

                int Dim = it->GetGeometry().WorkingSpaceDimension(); // global dimension of the geometry that it works on
            int ReducedDim = it->GetGeometry().LocalSpaceDimension(); // reduced dimension of the geometry
            IndexType NodeCounter_old = NodeCounter;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(Dim)
            KRATOS_WATCH(ReducedDim)
#endif

            //select the correct post element type
            std::string element_name;
            if ((Dim == 2) && (ReducedDim == 2))
            {
                element_name = std::string("KinematicLinear2D4N");
            }
            else if ((Dim == 3) && (ReducedDim == 2))
            {
                element_name = std::string("KinematicLinear2D4N");
            }
            else if ((Dim == 3) && (ReducedDim == 3))
            {
                element_name = std::string("KinematicLinear3D8N");
            }
            else
            {
                std::stringstream ss;
                ss << "Invalid dimension of ";
                ss << typeid(*it).name();
                ss << ", Dim = " << Dim;
                ss << ", ReducedDim = " << ReducedDim;
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__);
            }

            if (!KratosComponents<Element>::Has(element_name))
            {
                std::stringstream buffer;
                buffer << "Element " << element_name << " is not registered in Kratos.";
                buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
                KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
            }

            Element const& rCloneElement = KratosComponents<Element>::Get(element_name);

            GenerateForOneEntity<Element, ElementsArrayType, 1>(rModelPartPost,
                    *it, rCloneElement, NodeCounter_old, NodeCounter, ElementCounter, NodeKey, false,
                    dummy_ids, dummy_ids, false);

            ++show_progress;
        }
        KRATOS_WATCH(ElementCounter)

#ifdef DEBUG_LEVEL1
        std::cout << "Done generating for elements" << std::endl;
#endif

        IndexType ConditionCounter = 0;
        if (generate_for_condition)
        {
            Kratos::progress_display show_progress2( pConditions.size() );
            for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
            {
                // This is wrong, we will not kill the IS_INACTIVE conditions
                // TODO: to be deleted
                //            if(it->GetValue( IS_INACTIVE ))
                //            {
                ////                std::cout << "Condition " << it->Id() << " is inactive" << std::endl;
                //                ++show_progress2;
                //                continue;
                //            }
                if (it->pGetGeometry() == 0)
                    KRATOS_THROW_ERROR(std::logic_error, "Error: geometry is NULL at condition", it->Id())

                    int Dim = it->GetGeometry().WorkingSpaceDimension(); // global dimension of the geometry that it works on
                int ReducedDim = it->GetGeometry().LocalSpaceDimension(); // reduced dimension of the geometry
                IndexType NodeCounter_old = NodeCounter;

#ifdef DEBUG_LEVEL1
                KRATOS_WATCH(typeid(it->GetGeometry()).name())
                KRATOS_WATCH(Dim)
                KRATOS_WATCH(ReducedDim)
#endif

                //select the correct post condition type
                std::string condition_name;
                if (Dim == 3 && ReducedDim == 1)
                {
                    condition_name = std::string("LineForce3D2N");
                }
                else if (Dim == 3 && ReducedDim == 2)
                {
                    condition_name = std::string("FaceForce3D4N");
                }
                else
                {
                    std::stringstream ss;
                    ss << "Invalid dimension of ";
                    ss << typeid(*it).name();
                    ss << ", Dim = " << Dim;
                    ss << ", ReducedDim = " << ReducedDim;
                    ss << ". Condition " << it->Id() << " will be skipped.";
                    //                KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__);
                    continue;
                }

                if (!KratosComponents<Condition>::Has(condition_name))
                {
                    std::stringstream buffer;
                    buffer << "Condition " << condition_name << " is not registered in Kratos.";
                    buffer << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
                    KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
                }

                Condition const& rCloneCondition = KratosComponents<Condition>::Get(condition_name);

                GenerateForOneEntity<Condition, ConditionsArrayType, 2>(rModelPartPost,
                        *it, rCloneCondition, NodeCounter_old, NodeCounter, ConditionCounter, NodeKey, false,
                        dummy_ids, dummy_ids, false);

                ++show_progress2;
            }
            KRATOS_WATCH(ConditionCounter)
        }

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "GeneratePostModelPart2 completed: " << (end_compute - start_compute) << " s" << std::endl;
#else
        std::cout << "GeneratePostModelPart2 completed" << std::endl;
#endif
        std::cout << NodeCounter << " nodes and " << ElementCounter << " elements";
        if (generate_for_condition)
        {
            std::cout << ", " << ConditionCounter << " conditions";
        }
        std::cout << " are created" << std::endl;
    }

    // Generate the post model_part from reference model_part
    // this is the improved version of GenerateModelPart
    // which uses template function to generate post Elements for both Element and Condition
    // this version used a collapsing utility to collapse nodes automatically
    void GenerateModelPart2AutoCollapse(ModelPart& rModelPartPost,
                                        double dx, double dy, double dz, double tol)
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
#endif

#ifdef DEBUG_LEVEL1
        std::cout << typeid(*this).name() << "::GenerateModelPart" << std::endl;
#endif

        AutoCollapseSpatialBinning collapse_util(0.0, 0.0, 0.0, dx, dy, dz, tol);

        ElementsArrayType& pElements = mr_model_part.Elements();
        ConditionsArrayType& pConditions = mr_model_part.Conditions();

        std::string NodeKey = std::string("Node");

        IndexType NodeCounter = 0;
        IndexType ElementCounter = 0;
        Kratos::progress_display show_progress( pElements.size() );
        VectorMap<IndexType, IndexType> MapToCollapseNode;
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            bool is_active = true;
            if (it->IsDefined ( ACTIVE ))
            {
                is_active = it->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
            }
#else
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
            }
            else
            {
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = !it->GetValue( IS_INACTIVE );
                }
            }
#endif

            if (!is_active)
            {
//                std::cout << "Element " << it->Id() << " is inactive" << std::endl;
                ++show_progress;
                continue;
            }

            int Dim = it->GetGeometry().WorkingSpaceDimension(); // global dimension of the geometry that it works on
            int ReducedDim = it->GetGeometry().LocalSpaceDimension(); // reduced dimension of the geometry
            IndexType NodeCounter_old = NodeCounter;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(Dim)
            KRATOS_WATCH(ReducedDim)
#endif

            //select the correct post element type
            std::string element_name;
            if (Dim == 2 && ReducedDim == 2)
            {
                element_name = std::string("KinematicLinear2D4N");
            }
            else if (Dim == 3 && ReducedDim == 3)
            {
                element_name = std::string("KinematicLinear3D8N");
            }
            else
            {
                std::stringstream ss;
                ss << "Invalid dimension of ";
                ss << typeid(*it).name();
                ss << ", Dim = " << Dim;
                ss << ", ReducedDim = " << ReducedDim;
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__);
            }

            if (!KratosComponents<Element>::Has(element_name))
            {
                std::stringstream buffer;
                buffer << "Element " << element_name << " is not registered in Kratos.";
                buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
                KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
            }

            Element const& rCloneElement = KratosComponents<Element>::Get(element_name);

            GenerateForOneEntityAutoCollapse<Element, ElementsArrayType, 1>(collapse_util,
                    rModelPartPost, *it, rCloneElement, MapToCollapseNode, NodeCounter_old,
                    NodeCounter, ElementCounter, NodeKey);

            ++show_progress;
        }

#ifdef DEBUG_LEVEL1
        std::cout << "Done generating for elements" << std::endl;
#endif

        IndexType ConditionCounter = 0;
        Kratos::progress_display show_progress2( pConditions.size() );
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            bool is_active = true;
            if (it->IsDefined ( ACTIVE ))
            {
                is_active = it->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
            }
#else
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                }
            }
            else
            {
                if (it->Has ( IS_INACTIVE ))
                {
                    is_active = !it->GetValue( IS_INACTIVE );
                }
            }
#endif

            int Dim = it->GetGeometry().WorkingSpaceDimension(); // global dimension of the geometry that it works on
            int ReducedDim = it->GetGeometry().LocalSpaceDimension(); // reduced dimension of the geometry
            IndexType NodeCounter_old = NodeCounter;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(typeid(it->GetGeometry()).name())
            KRATOS_WATCH(Dim)
            KRATOS_WATCH(ReducedDim)
#endif

            //select the correct post condition type
            std::string condition_name;
            if (Dim == 3 && ReducedDim == 1)
            {
                condition_name = std::string("LineForce3D2N");
            }
            else if (Dim == 3 && ReducedDim == 2)
            {
                condition_name = std::string("FaceForce3D4N");
            }
            else
            {
                std::stringstream ss;
                ss << "Invalid dimension of ";
                ss << typeid(*it).name();
                ss << ", Dim = " << Dim;
                ss << ", ReducedDim = " << ReducedDim;
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), __FUNCTION__);
            }

            if (!KratosComponents<Condition>::Has(condition_name))
            {
                std::stringstream buffer;
                buffer << "Condition " << condition_name << " is not registered in Kratos.";
                buffer << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
                KRATOS_THROW_ERROR(std::runtime_error, buffer.str(), "");
            }

            Condition const& rCloneCondition = KratosComponents<Condition>::Get(condition_name);

            GenerateForOneEntityAutoCollapse<Condition, ConditionsArrayType, 2>(collapse_util,
                    rModelPartPost, *it, rCloneCondition, MapToCollapseNode, NodeCounter_old,
                    NodeCounter, ConditionCounter, NodeKey);

            ++show_progress2;
        }

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Generate PostModelPart completed: " << (end_compute - start_compute) << " s" << std::endl;
#else
        std::cout << "Generate PostModelPart completed" << std::endl;
#endif
        std::cout << NodeCounter << " nodes and " << ElementCounter << " elements" << ", " << ConditionCounter << " conditions are created" << std::endl;
    }

    /**
     * Utility function to generate elements/conditions for element/condition
     * if TEntityType==Element, type must be 1; if T==Condition, type is 2
     */
    template<class TEntityType, class TEntityContainerType, std::size_t type>
    void GenerateForOneEntity(ModelPart& rModelPart,
                              const TEntityType& rE,
                              TEntityType const& rSample,
                              IndexType NodeCounter_old,
                              IndexType& NodeCounter,
                              IndexType& EntityCounter,
                              const std::string& NodeKey,
                              bool transfer_nodal_var,
                              std::vector<std::size_t>& node_ids,
                              std::vector<std::size_t>& element_ids,
                              bool get_indices)
    {
        int ReducedDim = rE.GetGeometry().LocalSpaceDimension();

        //get the properties
        Properties::Pointer pProperties = rModelPart.pGetProperties(rE.pGetProperties()->Id());

#ifdef DEBUG_LEVEL1
        std::cout << "Generating for " << rE.Info() << std::endl;
        KRATOS_WATCH(*pProperties)
        KRATOS_WATCH(EntityCounter)
#endif

        // generate list of nodes
        if (ReducedDim == 1)
        {
            // TODO
        }
        else if (ReducedDim == 2)
        {
            IndexType NumDivision1 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_1) );
            IndexType NumDivision2 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_2) );
            IndexType i, j;
            CoordinatesArrayType p_ref;
            CoordinatesArrayType p;
            Vector shape_values;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(NumDivision1)
            KRATOS_WATCH(NumDivision2)
            std::cout << "Generating Nodes..." << std::endl;
#endif

            // create and add nodes
            p_ref[2] = 0.0;
            for (i = 0; i <= NumDivision1; ++i)
            {
                p_ref[0] = ((double) i) / NumDivision1;
                for (j = 0; j <= NumDivision2; ++j)
                {
                    p_ref[1] = ((double) j) / NumDivision2;

                    p = GlobalCoordinates(rE.GetGeometry(), p, p_ref);

                    NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                    pNewNode->SetId(++NodeCounter);

                    // Giving model part's variables list to the node
                    pNewNode->SetSolutionStepVariablesList(&rModelPart.GetNodalSolutionStepVariablesList());

                    //set buffer size
                    pNewNode->SetBufferSize(rModelPart.GetBufferSize());

                    rModelPart.AddNode(pNewNode);

                    if (type == 1)
                    {
                        mNodeToLocalCoordinates(pNewNode->Id()) = p_ref;
                        mNodeToElement(pNewNode->Id()) = rE.Id();
                    }

                    if (transfer_nodal_var)
                    {
                        shape_values = rE.GetGeometry().ShapeFunctionsValues(shape_values, p_ref);

                        VariablesList& var_list = rModelPart.GetNodalSolutionStepVariablesList();

                        for (VariablesList::const_iterator it = var_list.begin(); it != var_list.end(); ++it)
                        {
                            if (typeid(*it) == typeid(Variable<double>))
                            {
                                const Variable<double>& my_variable = dynamic_cast<const Variable<double>&>(*it);
                                double value = 0.0;
                                for (std::size_t n = 0; n < rE.GetGeometry().size(); ++n)
                                {
                                    value += shape_values[n] * rE.GetGeometry()[n].GetSolutionStepValue(my_variable);
                                }
                                pNewNode->GetSolutionStepValue(my_variable) = value;
                            }
                            else if (typeid(*it) == typeid(Variable<array_1d<double, 3> >))
                            {
                                const Variable<array_1d<double, 3> >& my_variable = dynamic_cast<const Variable<array_1d<double, 3> >&>(*it);
                                array_1d<double, 3> value;
                                noalias(value) = ZeroVector(3);
                                for (std::size_t n = 0; n < rE.GetGeometry().size(); ++n)
                                {
                                    noalias(value) += shape_values[n] * rE.GetGeometry()[n].GetSolutionStepValue(my_variable);
                                }
                                pNewNode->GetSolutionStepValue(my_variable) = value;
                            }
                        }
                    }

                    if (get_indices)
                    {
                        node_ids.push_back(pNewNode->Id());
                    }
                }
            }

            //for correct mapping to element, the repetitive node is allowed.
//            rModelPart.Nodes().Unique();

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rModelPart.Nodes().size())
            if (type == 1)
            {
                std::cout << "Generating Elements..." << std::endl;
            }
            else
            {
                std::cout << "Generating Conditions..." << std::endl;
            }
#endif

            // create and add element
//             typename T::NodesArrayType temp_nodes;
//             for(i = 0; i < NumDivision1; ++i)
//             {
//                 for(j = 0; j < NumDivision2; ++j)
//                 {
//                     int Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
//                     int Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
//                     int Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
//                     int Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

//                     // TODO: check if jacobian checking is necessary
//                     temp_nodes.clear();
//                     temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node1, NodeKey).base()));
//                     temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node2, NodeKey).base()));
//                     temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node4, NodeKey).base()));
//                     temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node3, NodeKey).base()));

//                     int NewEntityId = ++EntityCounter;
// //                    int NewEntityId = rE.Id(); ++EntityCounter;
//                     typename TEntityType::Pointer NewEntity = rSample.Create(NewEntityId, temp_nodes, pProperties);
//                     AddToModelPart<TEntityType>(rModelPart, NewEntity);
//                     if(type == 1)
//                         mOldToNewElements[rE.Id()].insert(NewEntityId);
//                     else if(type == 2)
//                         mOldToNewConditions[rE.Id()].insert(NewEntityId);
//                 }
//             }

            std::vector<std::vector<IndexType> > connectivities;

            for (i = 0; i < NumDivision1; ++i)
            {
                for (j = 0; j < NumDivision2; ++j)
                {
                    IndexType Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
                    IndexType Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
                    IndexType Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
                    IndexType Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

                    connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4, Node3});
                }
            }

            TEntityContainerType pNewEntities = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, TEntityType, TEntityContainerType>(
                                                    connectivities, rModelPart, rSample, EntityCounter, pProperties, NodeKey);

            for (typename TEntityContainerType::ptr_iterator it2 = pNewEntities.ptr_begin(); it2 != pNewEntities.ptr_end(); ++it2)
            {
                AddToModelPart<TEntityType>(rModelPart, *it2);
                if (type == 1)
                {
                    mOldToNewElements[rE.Id()].insert((*it2)->Id());
                }
                else if (type == 2)
                {
                    mOldToNewConditions[rE.Id()].insert((*it2)->Id());
                }
            }

            if (type == 1)
            {
                rModelPart.Elements().Unique();
            }
            else if (type == 2)
            {
                rModelPart.Conditions().Unique();
            }

#ifdef DEBUG_LEVEL1
            if (type == 1)
                KRATOS_WATCH(rModelPart.Elements().size())
                else
                    KRATOS_WATCH(rModelPart.Conditions().size())
#endif

                    if (get_indices)
                    {
                        for (typename TEntityContainerType::iterator it2 = pNewEntities.begin(); it2 != pNewEntities.end(); ++it2)
                        {
                            element_ids.push_back(it2->Id());
                        }
                    }
        }
        else if (ReducedDim == 3)
        {
            IndexType NumDivision1 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_1) );
            IndexType NumDivision2 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_2) );
            IndexType NumDivision3 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_3) );
            IndexType i, j, k;
            CoordinatesArrayType p_ref;
            CoordinatesArrayType p;
            Vector shape_values;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rE.Id())
            KRATOS_WATCH(NumDivision1)
            KRATOS_WATCH(NumDivision2)
            KRATOS_WATCH(NumDivision3)
            std::cout << "Generating Nodes..." << std::endl;
#endif

            // create and add nodes
            for (i = 0; i <= NumDivision1; ++i)
            {
                p_ref[0] = ((double) i) / NumDivision1;
                for (j = 0; j <= NumDivision2; ++j)
                {
                    p_ref[1] = ((double) j) / NumDivision2;
                    for (k = 0; k <= NumDivision3; ++k)
                    {
                        p_ref[2] = ((double) k) / NumDivision3;

                        p = GlobalCoordinates(rE.GetGeometry(), p, p_ref);

                        NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                        pNewNode->SetId(++NodeCounter);

#ifdef DEBUG_GENERATE_MESH
                        if (NodeCounter)
                        {
                            std::cout << "Node " << NodeCounter << " p_ref: " << p_ref << ", p: " << p << std::endl;
                        }
#endif

                        // Giving model part's variables list to the node
                        pNewNode->SetSolutionStepVariablesList(&rModelPart.GetNodalSolutionStepVariablesList());

                        //set buffer size
                        pNewNode->SetBufferSize(rModelPart.GetBufferSize());

                        rModelPart.AddNode(pNewNode);

                        if (type == 1)
                        {
                            mNodeToLocalCoordinates(pNewNode->Id()) = p_ref;
                            mNodeToElement(pNewNode->Id()) = rE.Id();
                        }

                        if (transfer_nodal_var)
                        {
                            shape_values = rE.GetGeometry().ShapeFunctionsValues(shape_values, p_ref);

                            VariablesList& var_list = rModelPart.GetNodalSolutionStepVariablesList();

                            for (VariablesList::const_iterator it = var_list.begin(); it != var_list.end(); ++it)
                            {
                                if (typeid(*it) == typeid(Variable<double>))
                                {
                                    const Variable<double>& my_variable = dynamic_cast<const Variable<double>&>(*it);
                                    double value = 0.0;
                                    for (std::size_t n = 0; n < rE.GetGeometry().size(); ++n)
                                    {
                                        value += shape_values[n] * rE.GetGeometry()[n].GetSolutionStepValue(my_variable);
                                    }
                                    pNewNode->GetSolutionStepValue(my_variable) = value;
                                }
                                else if (typeid(*it) == typeid(Variable<array_1d<double, 3> >))
                                {
                                    const Variable<array_1d<double, 3> >& my_variable = dynamic_cast<const Variable<array_1d<double, 3> >&>(*it);
                                    array_1d<double, 3> value;
                                    noalias(value) = ZeroVector(3);
                                    for (std::size_t n = 0; n < rE.GetGeometry().size(); ++n)
                                    {
                                        noalias(value) += shape_values[n] * rE.GetGeometry()[n].GetSolutionStepValue(my_variable);
                                    }
                                    pNewNode->GetSolutionStepValue(my_variable) = value;
                                }
                            }
                        }

                        if (get_indices)
                        {
                            node_ids.push_back(pNewNode->Id());
                        }
                    }
                }
            }

            //for correct mapping to element, the repetitive node is allowed.
//           rModelPart.Nodes().Unique();

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rModelPart.Nodes().size())
            if (type == 1)
            {
                std::cout << "Generating Elements..." << std::endl;
            }
            else
            {
                std::cout << "Generating Conditions..." << std::endl;
            }
#endif

            // create and add element
            // typename T::NodesArrayType temp_nodes;
            // for(i = 0; i < NumDivision1; ++i)
            // {
            //     for(j = 0; j < NumDivision2; ++j)
            //     {
            //         for(k = 0; k < NumDivision3; ++k)
            //         {
            //             int Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
            //             int Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
            //             int Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
            //             int Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
            //             int Node5 = Node1 + 1;
            //             int Node6 = Node2 + 1;
            //             int Node7 = Node3 + 1;
            //             int Node8 = Node4 + 1;

            //             // TODO: check if jacobian checking is necessary
            //             temp_nodes.clear();
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node1, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node2, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node4, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node3, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node5, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node6, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node8, NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), Node7, NodeKey).base()));

            //             int NewEntityId = ++EntityCounter;
            //             typename TEntityType::Pointer NewEntity = rSample.Create(NewEntityId, temp_nodes, pProperties);
            //             AddToModelPart<TEntityType>(rModelPart, NewEntity);
            //             if(type == 1)
            //                 mOldToNewElements[rE.Id()].insert(NewEntityId);
            //             else if(type == 2)
            //                 mOldToNewConditions[rE.Id()].insert(NewEntityId);
            //         }
            //     }
            // }

            std::vector<std::vector<IndexType> > connectivities;

            for (i = 0; i < NumDivision1; ++i)
            {
                for (j = 0; j < NumDivision2; ++j)
                {
                    for (k = 0; k < NumDivision3; ++k)
                    {
                        IndexType Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        IndexType Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        IndexType Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        IndexType Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        IndexType Node5 = Node1 + 1;
                        IndexType Node6 = Node2 + 1;
                        IndexType Node7 = Node3 + 1;
                        IndexType Node8 = Node4 + 1;

                        connectivities.push_back(std::vector<IndexType> {Node1, Node2, Node4, Node3, Node5, Node6, Node8, Node7});
                    }
                }
            }

            TEntityContainerType pNewEntities = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, TEntityType, TEntityContainerType>(
                                                    connectivities, rModelPart, rSample, EntityCounter, pProperties, NodeKey);

            for (typename TEntityContainerType::ptr_iterator it2 = pNewEntities.ptr_begin(); it2 != pNewEntities.ptr_end(); ++it2)
            {
                AddToModelPart<TEntityType>(rModelPart, *it2);
                if (type == 1)
                {
                    mOldToNewElements[rE.Id()].insert((*it2)->Id());
                }
                else if (type == 2)
                {
                    mOldToNewConditions[rE.Id()].insert((*it2)->Id());
                }
            }

            if (type == 1)
            {
                rModelPart.Elements().Unique();
            }
            else if (type == 2)
            {
                rModelPart.Conditions().Unique();
            }

#ifdef DEBUG_LEVEL1
            if (type == 1)
                KRATOS_WATCH(rModelPart.Elements().size())
                else
                    KRATOS_WATCH(rModelPart.Conditions().size())
#endif

                    if (get_indices)
                    {
                        for (typename TEntityContainerType::iterator it2 = pNewEntities.begin(); it2 != pNewEntities.end(); ++it2)
                        {
                            element_ids.push_back(it2->Id());
                        }
                    }
        }
    }

    /**
     * Utility function to generate elements/conditions for element/condition.
     * This uses a collapse utility to automatically merge the coincident nodes
     * if T==Element, type must be 1; otherwise type=2
     */
    template<class TEntityType, class TEntityContainerType, std::size_t type>
    void GenerateForOneEntityAutoCollapse(AutoCollapseSpatialBinning& collapse_util,
                                          ModelPart& rModelPart,
                                          const TEntityType& rE,
                                          TEntityType const& rSample,
                                          VectorMap<IndexType, IndexType>& rMapToCollapseNode,
                                          IndexType NodeCounter_old,
                                          IndexType& NodeCounter,
                                          IndexType& EntityCounter,
                                          const std::string& NodeKey)
    {
        int ReducedDim = rE.GetGeometry().LocalSpaceDimension();

        //get the properties
        Properties::Pointer pProperties = rModelPart.pGetProperties(rE.pGetProperties()->Id());

#ifdef DEBUG_LEVEL1
        if (type == 1)
        {
            std::cout << "Generating for element " << rE.Id() << std::endl;
        }
        else
        {
            std::cout << "Generating for condition " << rE.Id() << std::endl;
        }
        KRATOS_WATCH(*pProperties)
#endif

        // generate list of nodes
        if (ReducedDim == 1)
        {
            // TODO
        }
        else if (ReducedDim == 2)
        {
            IndexType NumDivision1 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_1) );
            IndexType NumDivision2 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_2) );
            IndexType i, j;
            CoordinatesArrayType p_ref;
            CoordinatesArrayType p;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(NumDivision1)
            KRATOS_WATCH(NumDivision2)
            std::cout << "Generating Nodes..." << std::endl;
#endif

            // create and add nodes
            p_ref[2] = 0.0;
            for (i = 0; i <= NumDivision1; ++i)
            {
                p_ref[0] = ((double) i) / NumDivision1;
                for (j = 0; j <= NumDivision2; ++j)
                {
                    p_ref[1] = ((double) j) / NumDivision2;

                    p = GlobalCoordinates(rE.GetGeometry(), p, p_ref);

                    IndexType id = static_cast<IndexType>( collapse_util.AddNode(p[0], p[1], p[2]) );
                    ++NodeCounter;
                    rMapToCollapseNode[NodeCounter] = id;

                    if (rModelPart.Nodes().find(id) == rModelPart.Nodes().end())
                    {
                        // this is a new node
                        NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                        pNewNode->SetId(id);

                        // Giving model part's variables list to the node
                        pNewNode->SetSolutionStepVariablesList(&rModelPart.GetNodalSolutionStepVariablesList());

                        //set buffer size
                        pNewNode->SetBufferSize(rModelPart.GetBufferSize());

                        rModelPart.AddNode(pNewNode);
                    }
                    else
                    {
                        // this is an old node, not required to add to model_part
                        // so do nothing
                    }

                    // in this way, the node will always point to the last local coodinates and element
                    if (type == 1)
                    {
                        mNodeToLocalCoordinates(id) = p_ref;
                        mNodeToElement(id) = rE.Id();
                    }
                }
            }

            //for correct mapping to element, the repetitive node is allowed.
//            rModelPart.Nodes().Unique();

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rModelPart.Nodes().size())
            if (type == 1)
            {
                std::cout << "Generating Elements..." << std::endl;
            }
            else
            {
                std::cout << "Generating Conditions..." << std::endl;
            }
#endif

            // create and add element
            // typename T::NodesArrayType temp_nodes;
            // for(i = 0; i < NumDivision1; ++i)
            // {
            //     for(j = 0; j < NumDivision2; ++j)
            //     {
            //         int Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
            //         int Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
            //         int Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
            //         int Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

            //         // TODO: check if jacobian checking is necessary
            //         temp_nodes.clear();
            //         temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node1], NodeKey).base()));
            //         temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node2], NodeKey).base()));
            //         temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node4], NodeKey).base()));
            //         temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node3], NodeKey).base()));

            //         typename T::Pointer NewEntity = rSample.Create(++EntityCounter, temp_nodes, pProperties);
            //         AddToModelPart<T>(rModelPart, NewEntity);
            //         if(type == 1)
            //             mOldToNewElements[rE.Id()].insert(EntityCounter);
            //         else if(type == 2)
            //             mOldToNewConditions[rE.Id()].insert(EntityCounter);
            //     }
            // }

            std::vector<std::vector<IndexType> > connectivities;

            for (i = 0; i < NumDivision1; ++i)
            {
                for (j = 0; j < NumDivision2; ++j)
                {
                    IndexType Node1 = NodeCounter_old + i * (NumDivision2 + 1) + j + 1;
                    IndexType Node2 = NodeCounter_old + i * (NumDivision2 + 1) + j + 2;
                    IndexType Node3 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 1;
                    IndexType Node4 = NodeCounter_old + (i + 1) * (NumDivision2 + 1) + j + 2;

                    connectivities.push_back(std::vector<IndexType>
                    {
                        rMapToCollapseNode[Node1],
                        rMapToCollapseNode[Node2],
                        rMapToCollapseNode[Node4],
                        rMapToCollapseNode[Node3]
                    });
                }
            }

            TEntityContainerType pNewEntities = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, TEntityType, TEntityContainerType>(
                                                    connectivities, rModelPart, rSample, EntityCounter, pProperties, NodeKey);

            for (typename TEntityContainerType::ptr_iterator it2 = pNewEntities.ptr_begin(); it2 != pNewEntities.ptr_end(); ++it2)
            {
                AddToModelPart<TEntityType>(rModelPart, *it2);
                if (type == 1)
                {
                    mOldToNewElements[rE.Id()].insert((*it2)->Id());
                }
                else if (type == 2)
                {
                    mOldToNewConditions[rE.Id()].insert((*it2)->Id());
                }
            }

            if (type == 1)
            {
                rModelPart.Elements().Unique();
            }
            else if (type == 2)
            {
                rModelPart.Conditions().Unique();
            }

#ifdef DEBUG_LEVEL1
            if (type == 1)
                KRATOS_WATCH(rModelPart.Elements().size())
                else
                    KRATOS_WATCH(rModelPart.Conditions().size())
#endif
                }
        else if (ReducedDim == 3)
        {
            IndexType NumDivision1 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_1) );
            IndexType NumDivision2 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_2) );
            IndexType NumDivision3 = static_cast<IndexType>( rE.GetValue(NUM_DIVISION_3) );
            IndexType i, j, k;
            CoordinatesArrayType p_ref;
            CoordinatesArrayType p;

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rE.Id())
            KRATOS_WATCH(NumDivision1)
            KRATOS_WATCH(NumDivision2)
            KRATOS_WATCH(NumDivision3)
            std::cout << "Generating Nodes..." << std::endl;
#endif

            // create and add nodes
            for (i = 0; i <= NumDivision1; ++i)
            {
                p_ref[0] = ((double) i) / NumDivision1;
                for (j = 0; j <= NumDivision2; ++j)
                {
                    p_ref[1] = ((double) j) / NumDivision2;
                    for (k = 0; k <= NumDivision3; ++k)
                    {
                        p_ref[2] = ((double) k) / NumDivision3;

                        p = GlobalCoordinates(rE.GetGeometry(), p, p_ref);

                        IndexType id = static_cast<IndexType>( collapse_util.AddNode(p[0], p[1], p[2]) );
                        ++NodeCounter;
                        rMapToCollapseNode[NodeCounter] = id;

                        if (rModelPart.Nodes().find(id) == rModelPart.Nodes().end())
                        {
                            // this is a new node
                            NodeType::Pointer pNewNode( new NodeType( 0, p ) );
                            pNewNode->SetId(id);

#ifdef DEBUG_GENERATE_MESH
                            if (NodeCounter)
                            {
                                std::cout << "Node " << NodeCounter << " p_ref: " << p_ref << ", p: " << p << std::endl;
                            }
#endif

                            // Giving model part's variables list to the node
                            pNewNode->SetSolutionStepVariablesList(&rModelPart.GetNodalSolutionStepVariablesList());

                            //set buffer size
                            pNewNode->SetBufferSize(rModelPart.GetBufferSize());

                            rModelPart.AddNode(pNewNode);
                        }
                        else
                        {
                            // this is an old node, not required to add to model_part
                            // so do nothing
                        }

                        // in this way, the node will always point to the last local coodinates and element
                        if (type == 1)
                        {
                            mNodeToLocalCoordinates(id) = p_ref;
                            mNodeToElement(id) = rE.Id();
                        }
                    }
                }
            }

            //for correct mapping to element, the repetitive node is allowed.
//           rModelPart.Nodes().Unique();

#ifdef DEBUG_LEVEL1
            KRATOS_WATCH(rModelPart.Nodes().size())
            if (type == 1)
            {
                std::cout << "Generating Elements..." << std::endl;
            }
            else
            {
                std::cout << "Generating Conditions..." << std::endl;
            }
#endif

            // create and add element
            // typename T::NodesArrayType temp_nodes;
            // for(i = 0; i < NumDivision1; ++i)
            // {
            //     for(j = 0; j < NumDivision2; ++j)
            //     {
            //         for(k = 0; k < NumDivision3; ++k)
            //         {
            //             int Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
            //             int Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
            //             int Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
            //             int Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
            //             int Node5 = Node1 + 1;
            //             int Node6 = Node2 + 1;
            //             int Node7 = Node3 + 1;
            //             int Node8 = Node4 + 1;

            //             // TODO: check if jacobian checking is necessary
            //             temp_nodes.clear();
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node1], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node2], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node4], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node3], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node5], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node6], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node8], NodeKey).base()));
            //             temp_nodes.push_back(*(FindKey(rModelPart.Nodes(), rMapToCollapseNode[Node7], NodeKey).base()));

            //             typename T::Pointer NewEntity = rSample.Create(++EntityCounter, temp_nodes, pProperties);
            //             AddToModelPart<T>(rModelPart, NewEntity);
            //             if(type == 1)
            //                 mOldToNewElements[rE.Id()].insert(EntityCounter);
            //             else if(type == 2)
            //                 mOldToNewConditions[rE.Id()].insert(EntityCounter);
            //         }
            //     }
            // }

            std::vector<std::vector<IndexType> > connectivities;

            for (i = 0; i < NumDivision1; ++i)
            {
                for (j = 0; j < NumDivision2; ++j)
                {
                    for (k = 0; k < NumDivision3; ++k)
                    {
                        IndexType Node1 = NodeCounter_old + (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        IndexType Node2 = NodeCounter_old + (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        IndexType Node3 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        IndexType Node4 = NodeCounter_old + ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        IndexType Node5 = Node1 + 1;
                        IndexType Node6 = Node2 + 1;
                        IndexType Node7 = Node3 + 1;
                        IndexType Node8 = Node4 + 1;

                        connectivities.push_back(std::vector<IndexType>
                        {
                            rMapToCollapseNode[Node1],
                            rMapToCollapseNode[Node2],
                            rMapToCollapseNode[Node4],
                            rMapToCollapseNode[Node3],
                            rMapToCollapseNode[Node5],
                            rMapToCollapseNode[Node6],
                            rMapToCollapseNode[Node8],
                            rMapToCollapseNode[Node7]
                        });
                    }
                }
            }

            TEntityContainerType pNewEntities = IsogeometricPostUtility::CreateEntities<std::vector<std::vector<IndexType> >, TEntityType, TEntityContainerType>(
                                                    connectivities, rModelPart, rSample, EntityCounter, pProperties, NodeKey);

            for (typename TEntityContainerType::ptr_iterator it2 = pNewEntities.ptr_begin(); it2 != pNewEntities.ptr_end(); ++it2)
            {
                AddToModelPart<TEntityType>(rModelPart, *it2);
                if (type == 1)
                {
                    mOldToNewElements[rE.Id()].insert((*it2)->Id());
                }
                else if (type == 2)
                {
                    mOldToNewConditions[rE.Id()].insert((*it2)->Id());
                }
            }

            if (type == 1)
            {
                rModelPart.Elements().Unique();
            }
            else if (type == 2)
            {
                rModelPart.Conditions().Unique();
            }

#ifdef DEBUG_LEVEL1
            if (type == 1)
                KRATOS_WATCH(rModelPart.Elements().size())
                else
                    KRATOS_WATCH(rModelPart.Conditions().size())
#endif
                }
    }

    // Synchronize the activation between model_parts
    void SynchronizeActivation(ModelPart& rModelPartPost)
    {
        ElementsArrayType& pElements = mr_model_part.Elements();
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            std::set<IndexType> NewElements = mOldToNewElements[it->Id()];
            for (std::set<IndexType>::iterator it2 = NewElements.begin(); it2 != NewElements.end(); ++it2)
            {
                if (it->IsDefined(ACTIVE))
                {
                    rModelPartPost.GetElement(*it2).Set(ACTIVE, it->Is( ACTIVE ));
                }
#ifndef SD_APP_FORWARD_COMPATIBILITY
                if (it->Has(IS_INACTIVE))
                {
                    rModelPartPost.GetElement(*it2).GetValue(IS_INACTIVE) = it->GetValue( IS_INACTIVE );
                }
#endif
            }
        }
        ConditionsArrayType& pConditions = mr_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            std::set<IndexType> NewConditions = mOldToNewConditions[it->Id()];
            for (std::set<IndexType>::iterator it2 = NewConditions.begin(); it2 != NewConditions.end(); ++it2)
            {
                if (it->IsDefined(ACTIVE))
                {
                    rModelPartPost.GetElement(*it2).Set(ACTIVE, it->Is( ACTIVE ));
                }
#ifndef SD_APP_FORWARD_COMPATIBILITY
                if (it->Has(IS_INACTIVE))
                {
                    rModelPartPost.GetCondition(*it2).GetValue(IS_INACTIVE) = it->GetValue( IS_INACTIVE );
                }
#endif
            }
        }
    }

    // transfer the elemental data
    template<class TVariableType>
    void TransferElementalData(const TVariableType& rThisVariable, ModelPart& rModelPartPost)
    {
        ElementsArrayType& pElements = mr_model_part.Elements();
        for (typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            std::set<IndexType> NewElements = mOldToNewElements[it->Id()];
            for (std::set<IndexType>::iterator it2 = NewElements.begin(); it2 != NewElements.end(); ++it2)
            {
                rModelPartPost.GetElement(*it2).GetValue(rThisVariable) = it->GetValue(rThisVariable);
            }
        }
    }

    // transfer the conditional data
    template<class TVariableType>
    void TransferConditionalData(const TVariableType& rThisVariable, ModelPart& rModelPartPost)
    {
        ConditionsArrayType& pConditions = mr_model_part.Conditions();
        for (typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            std::set<IndexType> NewConditions = mOldToNewConditions[it->Id()];
            for (std::set<IndexType>::iterator it2 = NewConditions.begin(); it2 != NewConditions.end(); ++it2)
            {
                rModelPartPost.GetCondition(*it2).GetValue(rThisVariable) = it->GetValue(rThisVariable);
            }
        }
    }

    // Synchronize post model_part with the reference model_part
    template<class TVariableType>
    void TransferNodalResults(
        const TVariableType& rThisVariable,
        ModelPart& rModelPartPost
    )
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
#endif

        NodesArrayType& pTargetNodes = rModelPartPost.Nodes();

        ElementsArrayType& pElements = mr_model_part.Elements();

        typename TVariableType::Type Results;
        CoordinatesArrayType LocalPos;
        IndexType ElementId;

//        #pragma omp parallel for
        //TODO: check this. This is not parallelized.
        for (NodesArrayType::iterator it = pTargetNodes.begin(); it != pTargetNodes.end(); ++it)
        {
            IndexType key = it->Id();
            if (mNodeToElement.find(key) != mNodeToElement.end())
            {
                ElementId = mNodeToElement[key];
                Element::Pointer pElement = pElements(ElementId);

                bool is_active = true;
                if (pElement->IsDefined ( ACTIVE ))
                {
                    is_active = pElement->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
                }
#else
                    if (pElement->Has ( IS_INACTIVE ))
                    {
                        is_active = is_active && (!pElement->GetValue( IS_INACTIVE ));
                    }
                }
                else
                {
                    if (pElement->Has ( IS_INACTIVE ))
                    {
                        is_active = !pElement->GetValue( IS_INACTIVE );
                    }
                }
#endif

                if ( is_active ) // skip the inactive elements
                {
                    noalias(LocalPos) = mNodeToLocalCoordinates[key];
                    Results = CalculateOnPoint(rThisVariable, Results, pElements(ElementId), LocalPos);
                    it->GetSolutionStepValue(rThisVariable) = Results;
                }
            }
        }

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer nodal point results for " << rThisVariable.Name() << " completed: " << end_compute - start_compute << " s" << std::endl;
#endif
    }

    // Synchronize post model_part with the reference model_part
    template<class TVariableType>
    void TransferIntegrationPointResults(
        const TVariableType& rThisVariable,
        ModelPart& rModelPartPost,
        LinearSolverType::Pointer pSolver
    )
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "########################################" << std::endl;
        std::cout << "Transfer integration point results for "
                  << rThisVariable.Name() << " starts" << std::endl;
#endif

        // firstly transfer rThisVariable from integration points of reference model_part to its nodes
        TransferVariablesToNodes(rThisVariable, mr_model_part, pSolver);

        // secondly transfer new nodal variables results to the post model_part
        TransferNodalResults(rThisVariable, rModelPartPost);

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer integration point results for "
                  << rThisVariable.Name() << " completed: "
                  << end_compute - start_compute << "s" << std::endl;
        std::cout << "########################################" << std::endl;
#endif
    }

    // Transfer the variable to nodes for model_part
    template<class TVariableType>
    void TransferVariablesToNodes(const TVariableType& rThisVariable,
                                  ModelPart& rModelPart, LinearSolverType::Pointer pSolver) const
    {
#ifdef ENABLE_PROFILING
        double start_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "########################################" << std::endl;
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " starts" << std::endl;
#endif

        TransferVariablesToNodes(rThisVariable, rModelPart, pSolver);

#ifdef ENABLE_PROFILING
        double end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "Transfer integration point results to nodes for "
                  << rThisVariable.Name() << " completed: "
                  << end_compute - start_compute << "s" << std::endl;
        std::cout << "########################################" << std::endl;
#endif
    }

    /**
     * Utility function to renumber the nodes of the post model_part (for parallel merge)
     */
    void GlobalNodalRenumbering(ModelPart& rModelPartPost)
    {
#ifdef ISOGEOMETRIC_USE_MPI
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // gather the number of nodes on each process
        int NumberOfNodes[size];
        int MyNumberOfNodes = rModelPartPost.NumberOfNodes();
        MPI_Allgather(&MyNumberOfNodes, 1, MPI_INT, NumberOfNodes, 1, MPI_INT, MPI_COMM_WORLD);
//        std::cout << "NumberOfNodes:";
//        for(int i = 0; i < size; ++i)
//            std::cout << " " << NumberOfNodes[i];
//        std::cout << std::endl;

        // compute the numbering offset
        int offset = 0;
        for (int i = 0; i < rank; ++i)
        {
            offset += NumberOfNodes[i];
        }

        // renumber the nodes of the current process
        for (ModelPart::NodeIterator it = rModelPartPost.NodesBegin(); it != rModelPartPost.NodesEnd(); ++it)
        {
            it->SetId(++offset);
            it->GetSolutionStepValue(PARTITION_INDEX) = rank;
        }
        if (rank == 0)
        {
            std::cout << "Global renumbering completed" << std::endl;
        }
#endif
    }

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "BezierClassicalPostUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BezierClassicalPostUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    ModelPart& mr_model_part; // reference to a model_part

    VectorMap<IndexType, CoordinatesArrayType> mNodeToLocalCoordinates; // vector map to store local coordinates of node on a NURBS entity
    VectorMap<IndexType, IndexType> mNodeToElement; // vector map to store local coordinates of node on a NURBS entity
    std::map<IndexType, std::set<IndexType> > mOldToNewElements; // vector map to store id map from old element to new elements
    std::map<IndexType, std::set<IndexType> > mOldToNewConditions; // vector map to store id map from old condition to new conditions

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Calculate global coodinates w.r.t initial configuration
     */
    CoordinatesArrayType& GlobalCoordinates(
        const GeometryType& rGeometry,
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates
    ) const
    {
        noalias( rResult ) = ZeroVector( 3 );

        Vector ShapeFunctionsValues;

        rGeometry.ShapeFunctionsValues(ShapeFunctionsValues, LocalCoordinates);

        for ( IndexType i = 0 ; i < rGeometry.size() ; ++i )
        {
            noalias( rResult ) += ShapeFunctionsValues( i ) * rGeometry.GetPoint( i ).GetInitialPosition();
        }

        return rResult;
    }

    /**
     * Interpolation on element
     */
    double CalculateOnPoint(
        const Variable<double>& rVariable,
        double& rResult,
        Element::Pointer pElement,
        const CoordinatesArrayType& rCoordinates
    ) const
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

    /**
     * Interpolation on element
     */
    Vector& CalculateOnPoint(
        const Variable<Vector>& rVariable,
        Vector& rResult,
        Element::Pointer pElement,
        const CoordinatesArrayType& rCoordinates
    ) const
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

    /**
     * Interpolation on element
     */
    array_1d<double, 3>& CalculateOnPoint(
        const Variable<array_1d<double, 3> >& rVariable,
        array_1d<double, 3>& rResult,
        Element::Pointer pElement,
        const CoordinatesArrayType& rCoordinates
    ) const
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

    /**
     * Transfer variable at integration points to nodes
     *
     * @param pSolver       the solver used for solving the local system matrix
     * @param pModelPart    pointer to model_part that we wish to transfer the result from its integration points to its nodes
     * @param rThisVariable the variable need to transfer the respected values
     */
    void TransferVariablesToNodes(const Variable<double>& rThisVariable,
                                  ModelPart& rModelPart, LinearSolverType::Pointer pSolver) const
    {
        ElementsArrayType& ElementsArray = rModelPart.Elements();

        //Initialize system of equations
        int NumberOfNodes = rModelPart.NumberOfNodes();
        SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
        noalias(M) = ZeroMatrix(NumberOfNodes, NumberOfNodes);

        SerialSparseSpaceType::VectorType g(NumberOfNodes);
        noalias(g) = ZeroVector(NumberOfNodes);

        SerialSparseSpaceType::VectorType b(NumberOfNodes);
        noalias(b) = ZeroVector(NumberOfNodes);

        // create the structure for M a priori
        ConstructL2MatrixStructure<Element>(M, ElementsArray);

        // Transfer of GaussianVariables to Nodal Variables via L_2-Minimization
        // see Jiao + Heath "Common-refinement-based data tranfer ..."
        // International Journal for numerical methods in engineering 61 (2004) 2402--2427
        // for general description of L_2-Minimization
        // set up the system of equations

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

        KRATOS_WATCH( number_of_threads )
        std::cout << "element_partition:";
        for (std::size_t i = 0; i < element_partition.size(); ++i)
        {
            std::cout << " " << element_partition[i];
        }
        std::cout << std::endl;

        //create the array of lock
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
            unsigned int row, col;

            typename ElementsArrayType::iterator it_begin = ElementsArray.begin() + element_partition[k];
            typename ElementsArrayType::iterator it_end = ElementsArray.begin() + element_partition[k + 1];

            for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
            {
                bool is_active = true;
                if (it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
                }
#else
                    if (it->Has ( IS_INACTIVE ))
                    {
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                    }
                }
                else
                {
                    if (it->Has ( IS_INACTIVE ))
                    {
                        is_active = !it->GetValue( IS_INACTIVE );
                    }
                }
#endif

                if (is_active)
                {
                    const IntegrationPointsArrayType& integration_points
                        = it->GetGeometry().IntegrationPoints(it->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());

                    //                J = it->GetGeometry().Jacobian(J, it->GetIntegrationMethod());
                    //                const Matrix& Ncontainer = it->GetGeometry().ShapeFunctionsValues(it->GetIntegrationMethod());

                    IsogeometricGeometryType& rIsogeometricGeometry = dynamic_cast<IsogeometricGeometryType&>(it->GetGeometry());
                    J = rIsogeometricGeometry.Jacobian0(J, it->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    Matrix Ncontainer;
                    rIsogeometricGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        it->GetIntegrationMethod()
                    );

                    // get the values at the integration_points
                    std::vector<double> ValuesOnIntPoint(integration_points.size());
                    it->CalculateOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, rModelPart.GetProcessInfo());

                    for (unsigned int point = 0; point < integration_points.size(); ++point)
                    {
                        MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                        double dV = DetJ * integration_points[point].Weight();
                        for (unsigned int prim = 0 ; prim < it->GetGeometry().size(); ++prim)
                        {
                            row = it->GetGeometry()[prim].Id() - 1;
                            omp_set_lock(&lock_array[row]);
                            b(row) += (ValuesOnIntPoint[point]) * Ncontainer(point, prim) * dV;
                            for (unsigned int sec = 0 ; sec < it->GetGeometry().size(); ++sec)
                            {
                                col = it->GetGeometry()[sec].Id() - 1;
                                M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                            }
                            omp_unset_lock(&lock_array[row]);
                        }
                    }
                }
                else
                {
                    // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                    for (unsigned int prim = 0 ; prim < it->GetGeometry().size(); ++prim)
                    {
                        row = it->GetGeometry()[prim].Id() - 1;
                        omp_set_lock(&lock_array[row]);
//                        b(row) += 0.0;
                        for (unsigned int sec = 0 ; sec < it->GetGeometry().size(); ++sec)
                        {
                            col = it->GetGeometry()[sec].Id() - 1;
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

        // solver the system
        pSolver->Solve(M, g, b);

        // transfer the solution to the nodal variables
        for (ModelPart::NodeIterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            it->GetSolutionStepValue(rThisVariable) = g((it->Id() - 1));
        }
    }

    /**
     * Transfer of rThisVariable defined on integration points to corresponding
     * nodal values. The transformation is done in a form that ensures a minimization
     * of L_2-norm error (/sum{rThisVariable- f(x)) whereas
     * f(x)= /sum{shape_func_i*rThisVariable_i}
     * @param model_part model_part on which the transfer should be done
     * @param rThisVariable Matrix-Variable which should be transferred
     * @see TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
     * @see TransferVariablesToNodes(ModelPart& model_part, Variable<double>& rThisVariable)
     * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
     * Journal for numer. meth. in eng. 61 (2004) 2402--2427
     * WARNING: this may cause segmentation faults as the respective variables
     * will be created on nodal level while they are originally intended to be
     * stored on integration points!
     *
     * @param pSolver       the solver used for solving the local system matrix
     * @param pModelPart    pointer to model_part that we wish to transfer the result from its integration points to its nodes
     * @param rThisVariable the variable need to transfer the respected values
     */
    void TransferVariablesToNodes(const Variable<Vector>& rThisVariable,
                                  ModelPart& rModelPart, LinearSolverType::Pointer& pSolver) const
    {
        ElementsArrayType& ElementsArray = rModelPart.Elements();

        const unsigned int& Dim = (*(ElementsArray.ptr_begin()))->GetGeometry().WorkingSpaceDimension();
        unsigned int VariableSize;
        if (rThisVariable.Name() == std::string("STRESSES")
                || rThisVariable.Name() == std::string("PLASTIC_STRAIN_VECTOR")
                || rThisVariable.Name() == std::string("PRESTRESS")
                || rThisVariable.Name() == std::string("STRAIN")
                // TODO: extend for more variables
           )
        {
            VariableSize = Dim * (Dim + 1) / 2;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, rThisVariable.Name(), " is not a supported variable for TransferVariablesToNodes routine.")

#ifdef ENABLE_PROFILING
            //profiling variables
            double start_compute, end_compute;
        start_compute = OpenMPUtils::GetCurrentTime();
#endif

        //Initialize system of equations
        unsigned int NumberOfNodes = rModelPart.NumberOfNodes();
        SerialSparseSpaceType::MatrixType M(NumberOfNodes, NumberOfNodes);
        noalias(M) = ZeroMatrix(NumberOfNodes, NumberOfNodes);

        // create the structure for M a priori
        ConstructL2MatrixStructure<Element>(M, ElementsArray);

#ifdef ENABLE_PROFILING
        end_compute = OpenMPUtils::GetCurrentTime();
        std::cout << "ConstructMatrixStructure completed: " << end_compute - start_compute << " s" << std::endl;
        start_compute = end_compute;
#endif

        SerialDenseSpaceType::MatrixType g(NumberOfNodes, VariableSize);
        noalias(g) = ZeroMatrix(NumberOfNodes, VariableSize);
        SerialDenseSpaceType::MatrixType b(NumberOfNodes, VariableSize);
        noalias(b) = ZeroMatrix(NumberOfNodes, VariableSize);

        std::vector< omp_lock_t > lock_array(NumberOfNodes);

        for (unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            omp_init_lock(&lock_array[i]);
        }

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, ElementsArray.size(), element_partition);

        KRATOS_WATCH( number_of_threads )
        std::cout << "element_partition:";
        for (std::size_t i = 0; i < element_partition.size(); ++i)
        {
            std::cout << " " << element_partition[i];
        }
        std::cout << std::endl;

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; ++k)
        {
            Matrix InvJ(Dim, Dim);
            double DetJ;
            unsigned int row, col;

            typename ElementsArrayType::iterator it_begin = ElementsArray.begin() + element_partition[k];
            typename ElementsArrayType::iterator it_end = ElementsArray.begin() + element_partition[k + 1];

            for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
            {
                bool is_active = true;
                if (it->IsDefined ( ACTIVE ))
                {
                    is_active = it->Is( ACTIVE );
#ifdef SD_APP_FORWARD_COMPATIBILITY
                }
#else
                    if (it->Has ( IS_INACTIVE ))
                    {
                        is_active = is_active && (!it->GetValue( IS_INACTIVE ));
                    }
                }
                else
                {
                    if (it->Has ( IS_INACTIVE ))
                    {
                        is_active = !it->GetValue( IS_INACTIVE );
                    }
                }
#endif

                if (is_active)
                {
                    const IntegrationPointsArrayType& integration_points
                        = it->GetGeometry().IntegrationPoints(it->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());

                    //                J = it->GetGeometry().Jacobian(J, it->GetIntegrationMethod());
                    //                const Matrix& Ncontainer = it->GetGeometry().ShapeFunctionsValues(it->GetIntegrationMethod());

                    IsogeometricGeometryType& rIsogeometricGeometry = dynamic_cast<IsogeometricGeometryType&>(it->GetGeometry());
                    J = rIsogeometricGeometry.Jacobian0(J, it->GetIntegrationMethod());

                    GeometryType::ShapeFunctionsGradientsType DN_De;
                    Matrix Ncontainer;
                    rIsogeometricGeometry.CalculateShapeFunctionsIntegrationPointsValuesAndLocalGradients(
                        Ncontainer,
                        DN_De,
                        it->GetIntegrationMethod()
                    );

                    // get the values at the integration_points
                    std::vector<Vector> ValuesOnIntPoint(integration_points.size());
                    it->CalculateOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, rModelPart.GetProcessInfo());

                    for (unsigned int point = 0; point < integration_points.size(); ++point)
                    {
                        MathUtils<double>::InvertMatrix(J[point], InvJ, DetJ);

                        double dV = DetJ * integration_points[point].Weight();

                        for (unsigned int prim = 0; prim < it->GetGeometry().size(); ++prim)
                        {
                            row = it->GetGeometry()[prim].Id() - 1;

                            omp_set_lock(&lock_array[row]);

                            for (unsigned int i = 0; i < VariableSize; ++i)
                            {
                                b(row, i) += ValuesOnIntPoint[point][i] * Ncontainer(point, prim) * dV;
                            }

                            for (unsigned int sec = 0; sec < it->GetGeometry().size(); ++sec)
                            {
                                col = it->GetGeometry()[sec].Id() - 1;
                                M(row, col) += Ncontainer(point, prim) * Ncontainer(point, sec) * dV;
                            }

                            omp_unset_lock(&lock_array[row]);
                        }
                    }
                }
                else
                {
                    // for inactive elements the contribution to LHS is identity matrix and RHS is zero
                    for (unsigned int prim = 0; prim < it->GetGeometry().size(); ++prim)
                    {
                        row = it->GetGeometry()[prim].Id() - 1;

                        omp_set_lock(&lock_array[row]);

//                        for(unsigned int i = 0; i < VariableSize; ++i)
//                            b(row, i) += 0.0;

                        for (unsigned int sec = 0; sec < it->GetGeometry().size(); ++sec)
                        {
                            col = it->GetGeometry()[sec].Id() - 1;
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
        pSolver->Solve(M, g, b);

#ifdef DEBUG_MULTISOLVE
        KRATOS_WATCH(g)
#endif

        // transfer the solution to the nodal variables
        for (ModelPart::NodeIterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            Vector tmp(VariableSize);
            for (unsigned int i = 0; i < VariableSize; ++i)
            {
                tmp(i) = g((it->Id() - 1), i);
            }
            it->GetSolutionStepValue(rThisVariable) = tmp;
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // /// Assignment operator.
    // BezierClassicalPostUtility& operator=(BezierClassicalPostUtility const& rOther)
    // {
    //     this->mr_model_part = rOther.mr_model_part;
    //     return *this;
    // }

    /// Copy constructor.
    BezierClassicalPostUtility(BezierClassicalPostUtility const& rOther)
        : mr_model_part(rOther.mr_model_part)
    {
    }

    ///@}

}; // Class BezierClassicalPostUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, BezierClassicalPostUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const BezierClassicalPostUtility& rThis)
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
