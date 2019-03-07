//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-10-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_POST_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_POST_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>
#include "boost/progress.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/isogeometric_utility.h"

#define USE_TRIANGULATION_UTILS_FOR_TRIANGULATION

#if defined(USE_TRIANGULATION_UTILS_FOR_TRIANGULATION)
#include "custom_utilities/triangulation_utils.h"
#elif defined(USE_CGAL_FOR_TRIANGULATION) && defined(ISOGEOMETRIC_APPLICATION_USE_CGAL)
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#endif

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

/// Short class definition.
/**
 * Abstract class for all utility to export mesh from NURBS. Also to provide basic utility functions.
 */
class IsogeometricPostUtility : public IsogeometricUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename ModelPart::NodesContainerType NodesArrayType;

    typedef typename ModelPart::ElementsContainerType ElementsArrayType;

    typedef typename ModelPart::ConditionsContainerType ConditionsArrayType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename NodeType::DofsContainerType DofsContainerType;

    typedef std::size_t IndexType;

    /// Pointer definition of IsogeometricPostUtility
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricPostUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricPostUtility()
    {
    }

    /// Destructor.
    virtual ~IsogeometricPostUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Find the entity of the same type in the list of entities
    template<class TEntityType, class TEntitiesContainerType>
    static TEntitiesContainerType FindEntities(TEntitiesContainerType& pEntities, TEntityType const& r_sample_entity)
    {
        TEntitiesContainerType pFoundEntities;

        for (typename TEntitiesContainerType::ptr_iterator it = pEntities.ptr_begin(); it != pEntities.ptr_end(); ++it)
        {
            if (typeid(*(*it)) == typeid(r_sample_entity))
                if (typeid((*it)->GetGeometry()) == typeid(r_sample_entity.GetGeometry()))
                    pFoundEntities.push_back(*it);
        }

        return pFoundEntities;
    }

    /// Create the entities based on the connectivities
    /// It is noted that the newly created entities are not added to the other model_part. User must do it manually.
    template<typename TConnectivityType, typename TEntityType, typename TEntitiesContainerType>
    static TEntitiesContainerType CreateEntities(
        const TConnectivityType& r_connectivities,
        ModelPart& r_model_part,
        TEntityType const& r_sample_entity,
        std::size_t& last_entity_id,
        Properties::Pointer pProperties,
        const std::string& NodeKey)
    {
        TEntitiesContainerType pNewEntities;
        typename TEntityType::NodesArrayType temp_entity_nodes;

        for (typename TConnectivityType::const_iterator it = r_connectivities.begin(); it != r_connectivities.end(); ++it)
        {
            temp_entity_nodes.clear();

            for (typename TConnectivityType::value_type::const_iterator it2 = it->begin(); it2 != it->end(); ++it2)
                temp_entity_nodes.push_back(*(FindKey(r_model_part.Nodes(), *it2, NodeKey).base()));

            typename TEntityType::Pointer pNewEntity = r_sample_entity.Create(++last_entity_id, temp_entity_nodes, pProperties);
            pNewEntities.push_back(pNewEntity);
        }

        return pNewEntities;
    }

    /// Create a list of entities (element/condition) (from a model_part) to another model_part
    /// It is noted that the newly created entities are not added to the other model_part. User must do it manually.
    template<class TEntityType, class TEntitiesContainerType>
    static TEntitiesContainerType CreateEntities(
        TEntitiesContainerType& pEntities,
        ModelPart& r_other_model_part,
        TEntityType const& r_sample_entity,
        std::size_t& last_entity_id,
        Properties::Pointer pProperties)
    {
        // first collect all the nodes from the elements
        std::map<std::size_t, NodeType::Pointer> pNodes;

        for (typename TEntitiesContainerType::ptr_iterator it = pEntities.ptr_begin(); it != pEntities.ptr_end(); ++it)
        {
            for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
            {
                pNodes[(*it)->GetGeometry()[i].Id()] = (*it)->GetGeometry().pGetPoint(i);
            }
        }

        // create the new nodes in the other model_part
        std::size_t last_node_id = GetLastNodeId(r_other_model_part);
        std::map<std::size_t, std::size_t> MapOldToNew;
        for (std::map<std::size_t, NodeType::Pointer>::iterator it = pNodes.begin(); it != pNodes.end(); ++it)
        {
            const PointType& rPoint = it->second->GetInitialPosition();
            NodeType::Pointer pNewNode = r_other_model_part.CreateNewNode(++last_node_id, rPoint[0], rPoint[1], rPoint[2]);
            MapOldToNew[it->second->Id()] = last_node_id;
        }

        // create new elements in the other model_part
        const std::string NodeKey = std::string("Node");
        typename TEntityType::NodesArrayType temp_entity_nodes;
        TEntitiesContainerType pNewEntities;
        for (typename TEntitiesContainerType::ptr_iterator it = pEntities.ptr_begin(); it != pEntities.ptr_end(); ++it)
        {
            temp_entity_nodes.clear();
            for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
            {
                std::size_t node_id = MapOldToNew[(*it)->GetGeometry()[i].Id()];
                temp_entity_nodes.push_back(*(FindKey(r_other_model_part.Nodes(), node_id, NodeKey).base()));
            }
            pNewEntities.push_back(r_sample_entity.Create(++last_entity_id, temp_entity_nodes, pProperties));
        }

        return pNewEntities;
    }

    /// Create a list of conditions/elements from the list of points.
    /// The point list will be triangulated before the conditions are created.
    /// It is noted that the newly created entities are not added to the other model_part. User must do it manually. Nevertheless, the nodes are added to the model_part.
    template<typename TPointType, typename TVectorType, class TEntityType, class TNodesContainerType, class TEntitiesContainerType>
    static std::pair<TNodesContainerType, TEntitiesContainerType> CreateEntities(
        const std::vector<TPointType>& points,
        const TVectorType& rCenter,
        const TVectorType& rNormal,
        const TVectorType& rTangent1,
        const TVectorType& rTangent2,
        ModelPart& r_model_part,
        TEntityType const& r_sample_entity,
        std::size_t& last_entity_id,
        Properties::Pointer pProperties)
    {
        // create the 2D coordinates for points, in order to triangulate
        std::vector<double> XY;
        TPointType Projection;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            noalias(Projection) = points[i] - inner_prod(points[i] - rCenter, rNormal) * rNormal;
            XY.push_back(inner_prod(Projection - rCenter, rTangent1));
            XY.push_back(inner_prod(Projection - rCenter, rTangent2));
        }

        // std::cout << "XY:" << std::endl;
        // for (std::size_t i = 0; i < XY.size()/2; ++i)
        //     std::cout << " " << XY[2*i] << " " << XY[2*i+1] << std::endl;

        // compute the triangulation
        std::vector<std::vector<unsigned int> > Connectivities;
        #if defined(USE_CGAL_FOR_TRIANGULATION)
        typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
        typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>   Vb;
        typedef CGAL::Triangulation_data_structure_2<Vb>                            Tds;
        typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                         Delaunay;
        typedef Kernel::Point_2                                                     Point2;

        std::vector< std::pair<Point2, unsigned int> > clipped_points;
        for(std::size_t i = 0; i < XY.size() / 2; ++i)
        {
            clipped_points.push_back( std::make_pair( Point2(XY[2*i], XY[2*i+1]), i ) );
        }

        Delaunay triangulation;
        triangulation.insert(clipped_points.begin(), clipped_points.end());

        for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); ++fit)
        {
            Delaunay::Face_handle face = fit;
            std::vector<unsigned int> con(3);
            con[0] = face->vertex(0)->info();
            con[1] = face->vertex(1)->info();
            con[2] = face->vertex(2)->info();
            Connectivities.push_back(con);
        }
        #elif defined(USE_TRIANGULATION_UTILS_FOR_TRIANGULATION)
        TriangulationUtils tri_util;
        tri_util.ComputeDelaunayTriangulation(XY, Connectivities);
        #else
        // REMARK: a tool to perform triangulation is not defined. You must define it.
        KRATOS_THROW_ERROR(std::logic_error, "A triangulation method must be specialized", "")
        #endif

        // std::cout << "Connectivities:" << std::endl;
        // for (std::size_t i = 0; i < Connectivities.size(); ++i)
        // {
        //     std::cout << "  " << i << ":";
        //     for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
        //         std::cout << " " << Connectivities[i][j];
        //     std::cout << std::endl;
        // }

        // create the nodes
        std::vector<std::size_t> map_con_to_mp(points.size());
        std::size_t last_node_id = GetLastNodeId(r_model_part);
        TNodesContainerType pNewNodes;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id, points[i][0], points[i][1], points[i][2]);
            map_con_to_mp[i] = pNewNode->Id();
            pNewNodes.push_back(pNewNode);
        }

        // create the entities based on connectivity
        const std::string NodeKey = std::string("Node");
        typename TEntityType::NodesArrayType temp_entity_nodes;
        TEntitiesContainerType pNewEntities;
        for (std::size_t i = 0; i < Connectivities.size(); ++i)
        {
            temp_entity_nodes.clear();
            for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
            {
                std::size_t node_id = map_con_to_mp[Connectivities[i][j]];
                temp_entity_nodes.push_back(*(FindKey(r_model_part.Nodes(), node_id, NodeKey).base()));
            }
            pNewEntities.push_back(r_sample_entity.Create(++last_entity_id, temp_entity_nodes, pProperties));
        }

        return std::make_pair(pNewNodes, pNewEntities);
    }

    //**********AUXILIARY FUNCTION**************************************************************
    // Construct the matrix structure for high performance assembling
    // This subroutine shall only be used to construct the matrix structure for L2 projection
    // using in post-processing
    //******************************************************************************************
    template<typename TElementType, typename TCompressedMatrixType, typename TElementsArrayType>
    static void ConstructL2MatrixStructure (
        TCompressedMatrixType& A,
        TElementsArrayType& rElements,
        std::map<std::size_t, std::size_t> MapNodeIdToVec)
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        typename TElementType::EquationIdVectorType ids;
        for(typename TElementsArrayType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; ++i_element)
        {
            ids.resize((i_element)->GetGeometry().size());
            for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                ids[i] = MapNodeIdToVec[(i_element)->GetGeometry()[i].Id()];

            for(std::size_t i = 0 ; i < ids.size() ; ++i)
            {
                if(ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for(std::size_t j = 0 ; j < ids.size() ; ++j)
                    {
                        if(ids[j] < equation_size)
                            AddUnique(row_indices, ids[j]);
                    }
                }
            }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; ++i)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
#ifndef _OPENMP
        for(std::size_t i = 0 ; i < indices.size() ; ++i)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; ++it)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> matrix_partition;
        OpenMPUtils::CreatePartition(number_of_threads, indices.size(), matrix_partition);
        for( int k=0; k < number_of_threads; ++k )
        {
            #pragma omp parallel
            if( omp_get_thread_num() == k )
            {
                for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; i++ )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for(std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; ++it)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    // Construct the matrix structure for high performance assembling
    // This subroutine shall only be used to construct the matrix structure for L2 projection
    // using in post-processing
    //******************************************************************************************
    template<typename TElementType, typename TCompressedMatrixType, typename TElementsArrayType>
    static void ConstructL2MatrixStructure (
        TCompressedMatrixType& A,
        TElementsArrayType& rElements)
    {
        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);

        typename TElementType::EquationIdVectorType ids;
        for(typename TElementsArrayType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; ++i_element)
        {
            ids.resize((i_element)->GetGeometry().size());
            for(unsigned int i = 0; i < (i_element)->GetGeometry().size();  ++i)
                ids[i] = (i_element)->GetGeometry()[i].Id() - 1;

            for(std::size_t i = 0 ; i < ids.size() ; ++i)
            {
                if(ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for(std::size_t j = 0 ; j < ids.size() ; ++j)
                    {
                        if(ids[j] < equation_size)
                            AddUnique(row_indices, ids[j]);
                    }
                }
            }
        }

        //allocating the memory needed
        int data_size = 0;
        for(std::size_t i = 0 ; i < indices.size() ; ++i)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
#ifndef _OPENMP
        for(std::size_t i = 0 ; i < indices.size() ; i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        std::vector<unsigned int> matrix_partition;
        OpenMPUtils::CreatePartition(number_of_threads, indices.size(), matrix_partition);
        for( int k=0; k < number_of_threads; ++k )
        {
            #pragma omp parallel
            if( omp_get_thread_num() == k )
            {
                for( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; i++ )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
    }

    //**********AUXILIARY FUNCTION**************************************************************
    // Support function for ConstructMatrixStructure
    //******************************************************************************************
    static inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while ( i != endit && (*i) != candidate)
        {
            ++i;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }
    }

    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    static inline double CoordinateScaling(const double& x, const int& Type)
    {
        if(Type == _NURBS_)
        {
            return x;
        }
        else if(Type == _BEZIER_)
        {
            return 2 * x - 1;
        }
        else
            return 0.0;
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
        buffer << "IsogeometricPostUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricPostUtility";
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IsogeometricPostUtility& operator=(IsogeometricPostUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricPostUtility(IsogeometricPostUtility const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricPostUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, IsogeometricPostUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const IsogeometricPostUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_POST_UTILITY_H_INCLUDED

