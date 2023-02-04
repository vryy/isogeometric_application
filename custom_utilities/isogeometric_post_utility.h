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
#include <tuple>
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

    /// Create a node for a model_part with a specific Id and transfer the values
    template<class TPatchType, typename TCoordinatesType, typename TIndexType>
    static typename NodeType::Pointer CreateNodeAndTransferValues(const TCoordinatesType& p_ref, const TPatchType& rPatch,
        ModelPart& r_model_part, const TIndexType& NodeId)
    {
        typename NodeType::Pointer pNewNode = CreateNode(p_ref, rPatch, r_model_part, NodeId);
        TransferValuesToNodes(*pNewNode, p_ref, rPatch);
        return pNewNode;
    }

    /// Create a node for a model_part with a specific Id
    template<class TPatchType, typename TCoordinatesType, typename TIndexType>
    static typename NodeType::Pointer CreateNode(const TCoordinatesType& p_ref, const TPatchType& rPatch,
        ModelPart& r_model_part, const TIndexType& NodeId)
    {
        typename TPatchType::ControlPointType p = rPatch.pControlPointGridFunction()->GetValue(p_ref);
        typename NodeType::Pointer pNewNode = r_model_part.CreateNewNode(NodeId, p.X(), p.Y(), p.Z());
        return pNewNode;
    }

    /// Transfer the control values from patch to node
    /// The node has to be inside the patch
    template<class TPatchType>
    static void TransferValuesToNodes(NodeType& rNode, const TPatchType& rPatch)
    {
        typename TPatchType::Array1DGridFunctionType::ConstPointer pControlPointCoordinatesGridFunction = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES);

        typename TPatchType::Array1DGridFunctionType::DataType p_ref;
        pControlPointCoordinatesGridFunction->LocalCoordinates(rNode, p_ref);

        TransferValuesToNodes(rNode, p_ref, rPatch);
    }

    /// Transfer the control values from patch to node
    /// p_ref is the local coordinates of the node in patch
    template<class TPatchType, typename TCoordinatesType>
    static void TransferValuesToNodes(NodeType& rNode, const TCoordinatesType& p_ref, const TPatchType& rPatch)
    {
        typedef typename TPatchType::DoubleGridFunctionContainerType DoubleGridFunctionContainerType;
        typedef typename TPatchType::Array1DGridFunctionContainerType Array1DGridFunctionContainerType;
        typedef typename TPatchType::VectorGridFunctionContainerType VectorGridFunctionContainerType;

        // transfer the control values
        DoubleGridFunctionContainerType DoubleGridFunctions_ = rPatch.DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it_gf = DoubleGridFunctions_.begin();
                it_gf != DoubleGridFunctions_.end(); ++it_gf)
        {
            typedef double DataType;
            typedef Variable<DataType> VariableType;
            const std::string& var_name = (*it_gf)->pControlGrid()->Name();
            if (KratosComponents<VariableData>::Has(var_name))
            {
                VariableType* pVariable = dynamic_cast<VariableType*>(&KratosComponents<VariableData>::Get(var_name));
                DataType value = (*it_gf)->GetValue(p_ref);
                if (rNode.SolutionStepsDataHas(*pVariable))
                    rNode.GetSolutionStepValue(*pVariable) = value;
            }
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = rPatch.Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it_gf = Array1DGridFunctions_.begin();
                it_gf != Array1DGridFunctions_.end(); ++it_gf)
        {
            typedef array_1d<double, 3> DataType;
            typedef Variable<DataType> VariableType;
            const std::string& var_name = (*it_gf)->pControlGrid()->Name();
            if (var_name == "CONTROL_POINT_COORDINATES") continue;
            if (KratosComponents<VariableData>::Has(var_name))
            {
                VariableType* pVariable = dynamic_cast<VariableType*>(&KratosComponents<VariableData>::Get(var_name));
                DataType value = (*it_gf)->GetValue(p_ref);
                if (rNode.SolutionStepsDataHas(*pVariable))
                    rNode.GetSolutionStepValue(*pVariable) = value;
            }
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = rPatch.VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it_gf = VectorGridFunctions_.begin();
                it_gf != VectorGridFunctions_.end(); ++it_gf)
        {
            typedef Vector DataType;
            typedef Variable<DataType> VariableType;
            const std::string& var_name = (*it_gf)->pControlGrid()->Name();
            if (KratosComponents<VariableData>::Has(var_name))
            {
                VariableType* pVariable = dynamic_cast<VariableType*>(&KratosComponents<VariableData>::Get(var_name));
                DataType value = (*it_gf)->GetValue(p_ref);
                if (rNode.SolutionStepsDataHas(*pVariable))
                    rNode.GetSolutionStepValue(*pVariable) = value;
            }
        }
    }

    /// Transfer the control values from patch to Gauss points
    template<class TEntityType, typename TVariableType, class TPatchType>
    static void TransferValuesToGaussPoints(TEntityType& rElement, const TVariableType& rVariable,
        const TPatchType& rPatch, const ProcessInfo& rProcessInfo)
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod = rElement.GetIntegrationMethod();

        GeometryType& rGeometry = rElement.GetGeometry();

        typename TPatchType::Array1DGridFunctionType::ConstPointer pControlPointCoordinatesGridFunction
            = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES);

        typename GridFunction<TPatchType::FESpaceType::Dim(), typename TVariableType::Type>::ConstPointer pGridFunc
            = rPatch.pGetGridFunction(rVariable);

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        rGeometry.Initialize(ThisIntegrationMethod);
        #endif

        const IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints(ThisIntegrationMethod);

        std::vector<typename TVariableType::Type> ValuesOnIntPoint(integration_points.size());

        CoordinatesArrayType GlobalCoords;
        typename TPatchType::Array1DGridFunctionType::DataType p_ref;

        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
        {
            rGeometry.GlobalCoordinates(GlobalCoords, integration_points[PointNumber]);

            typename TPatchType::Array1DGridFunctionType::ConstPointer pControlPointCoordinatesGridFunction = rPatch.pGetGridFunction(CONTROL_POINT_COORDINATES);

            pControlPointCoordinatesGridFunction->LocalCoordinates(GlobalCoords, p_ref);

            ValuesOnIntPoint[PointNumber] = pGridFunc->GetValue(p_ref);
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry
        rGeometry.Clean();
        #endif

        rElement.SetValuesOnIntegrationPoints( rVariable, ValuesOnIntPoint, rProcessInfo);
    }

    /// Generate corner points for regular geometry
    template<int TDim, typename TCoordinatesType, typename TValueType>
    static void GenerateRegular(std::vector<TCoordinatesType>& points,
        const std::vector<TCoordinatesType>& cmin, const std::vector<TCoordinatesType>& cmax)
    {
        if (TDim == 2)
        {
            GenerateRectangle(points, cmin[0], cmax[0], cmin[1], cmax[1]);
        }
        else if (TDim == 3)
        {
            GenerateBox(points, cmin[0], cmax[0], cmin[1], cmax[1], cmin[2], cmax[2]);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid dimension", TDim)
    }

    /// Generate a single rectangle. The 4 corner points are denoted as
    ///  4---3
    ///  |   |
    ///  1---2
    template<typename TCoordinatesType, typename TValueType>
    static void GenerateRectangle(std::vector<TCoordinatesType>& points,
        const TValueType& xmin, const TValueType& xmax,
        const TValueType& ymin, const TValueType& ymax)
    {
        points[0][0] = xmin;
        points[0][1] = ymin;
        // points[0][2] = 0.0;

        points[1][0] = xmax;
        points[1][1] = ymin;
        // points[1][2] = 0.0;

        points[2][0] = xmax;
        points[2][1] = ymax;
        // points[2][2] = 0.0;

        points[3][0] = xmin;
        points[3][1] = ymax;
        // points[3][2] = 0.0;
    }

    /// Generate the triangulation for a list of points in 3D
    /// The triangulation will be performed on the physical points with the information {cemter, normal, t1, t2}
    template<typename TCoordinatesType, typename TVectorType, typename TIndexType>
    static std::vector<std::vector<TIndexType> >
    GenerateTriangleGrid(const std::vector<TCoordinatesType>& points,
        const TVectorType& rCenter,
        const TVectorType& rNormal,
        const TVectorType& rTangent1,
        const TVectorType& rTangent2)
    {
        // create the 2D coordinates for points, in order to triangulate
        std::vector<double> XY;
        TCoordinatesType Projection;
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
        typedef std::vector<std::vector<TIndexType> > connectivity_t;
        connectivity_t Connectivities;
        #if defined(USE_CGAL_FOR_TRIANGULATION) && defined(ISOGEOMETRIC_APPLICATION_USE_CGAL)
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

        return Connectivities;
    }

    /// Generate the triangulation for a list of points in 3D
    /// The triangulation will be performed on the physical points with the information {cemter, normal, t1, t2}
    /// The refinement is performed on the local points instead.
    template<typename TCoordinatesType, typename TVectorType, typename TIndexType>
    static std::pair<std::vector<TCoordinatesType>, std::vector<std::vector<TIndexType> > >
    GenerateTriangleGrid(const std::vector<TCoordinatesType>& physical_points,
        const TVectorType& rCenter,
        const TVectorType& rNormal,
        const TVectorType& rTangent1,
        const TVectorType& rTangent2,
        const std::vector<TCoordinatesType>& local_points,
        const TIndexType& offset,
        const std::size_t& nrefine)
    {
        // compute the triangulation
        typedef std::vector<std::vector<TIndexType> > connectivity_t;
        connectivity_t Connectivities = GenerateTriangleGrid<TCoordinatesType, TVectorType, TIndexType>(physical_points, rCenter, rNormal, rTangent1, rTangent2);

        // refine if needed
        std::vector<TCoordinatesType> new_points = local_points;
        for (std::size_t i = 0; i < nrefine; ++i)
            RefineTriangleGrid<TIndexType, TCoordinatesType>(new_points, Connectivities);

        // offset the connectivity
        for (std::size_t i = 0; i < Connectivities.size(); ++i)
            for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
                Connectivities[i][j] += offset;

        return std::make_pair(new_points, Connectivities);
    }

    /// Generate the quadrilateral grid. The 4 corner points are denoted as
    ///  4---3
    ///  |   |
    ///  1---2
    template<typename TCoordinatesType, typename TIndexType>
    static std::pair<std::vector<TCoordinatesType>, std::vector<std::vector<TIndexType> > >
    GenerateQuadGrid(const TCoordinatesType& p1, const TCoordinatesType& p2,
        const TCoordinatesType& p3, const TCoordinatesType& p4,
        const TIndexType& starting_node_id,
        const std::size_t& num_div_1, const std::size_t& num_div_2)
    {
        TCoordinatesType p, pm, pn;

        std::vector<TCoordinatesType> points;
        std::vector<std::vector<TIndexType> > connectivities;

        double xi, eta;
        std::size_t i, j;
        for (i = 0; i <= num_div_1; ++i)
        {
            xi = ((double) i) / num_div_1;
            pm = p1 + xi*(p2 - p1);
            pn = p4 + xi*(p3 - p4);

            for (j = 0; j <= num_div_2; ++j)
            {
                eta = ((double) j) / num_div_2;
                p = pm + eta*(pn - pm);

                points.push_back(p);
            }
        }

        TIndexType n1, n2, n3, n4;
        for (i = 0; i < num_div_1; ++i)
        {
            for(j = 0; j < num_div_2; ++j)
            {
                n1 = starting_node_id + i * (num_div_2 + 1) + j;
                n2 = starting_node_id + i * (num_div_2 + 1) + j + 1;
                n3 = starting_node_id + (i + 1) * (num_div_2 + 1) + j;
                n4 = starting_node_id + (i + 1) * (num_div_2 + 1) + j + 1;

                connectivities.push_back(std::vector<std::size_t>{n1, n2, n4, n3});
            }
        }

        return std::make_pair(points, connectivities);
    }

    /// Generate a single box. The 8 corner points are denoted as
    ///  4---3         8---7
    ///  |   |   -->   |   |
    ///  1---2         5---6
    template<typename TCoordinatesType, typename TValueType>
    static void GenerateBox(std::vector<TCoordinatesType>& points,
        const TValueType& xmin, const TValueType& xmax,
        const TValueType& ymin, const TValueType& ymax,
        const TValueType& zmin, const TValueType& zmax)
    {
        points[0][0] = xmin;
        points[0][1] = ymin;
        points[0][2] = zmin;

        points[1][0] = xmax;
        points[1][1] = ymin;
        points[1][2] = zmin;

        points[2][0] = xmax;
        points[2][1] = ymax;
        points[2][2] = zmin;

        points[3][0] = xmin;
        points[3][1] = ymax;
        points[3][2] = zmin;

        points[4][0] = xmin;
        points[4][1] = ymin;
        points[4][2] = zmax;

        points[5][0] = xmax;
        points[5][1] = ymin;
        points[5][2] = zmax;

        points[6][0] = xmax;
        points[6][1] = ymax;
        points[6][2] = zmax;

        points[7][0] = xmin;
        points[7][1] = ymax;
        points[7][2] = zmax;
    }

    /// Generate the hexahedral grid. The 8 corner points are denoted as
    ///  4---3         8---7
    ///  |   |   -->   |   |
    ///  1---2         5---6
    template<typename TCoordinatesType, typename TIndexType>
    static std::pair<std::vector<TCoordinatesType>, std::vector<std::vector<TIndexType> > >
    GenerateHexGrid(const TCoordinatesType& p1, const TCoordinatesType& p2,
        const TCoordinatesType& p3, const TCoordinatesType& p4,
        const TCoordinatesType& p5, const TCoordinatesType& p6,
        const TCoordinatesType& p7, const TCoordinatesType& p8,
        const TIndexType& starting_node_id,
        const std::size_t& num_div_1, const std::size_t& num_div_2, const std::size_t& num_div_3)
    {
        TCoordinatesType p, pm1, pn1, pm2, pn2, pq1, pq2;

        std::vector<TCoordinatesType> points;
        std::vector<std::vector<TIndexType> > connectivities;

        double xi, eta, zeta;
        std::size_t i, j, k;
        for (i = 0; i <= num_div_1; ++i)
        {
            xi = ((double) i) / num_div_1;
            pm1 = p1 + xi*(p2 - p1);
            pn1 = p4 + xi*(p3 - p4);
            pm2 = p5 + xi*(p6 - p5);
            pn2 = p8 + xi*(p7 - p8);

            for (j = 0; j <= num_div_2; ++j)
            {
                eta = ((double) j) / num_div_2;
                pq1 = pm1 + eta*(pn1 - pm1);
                pq2 = pm2 + eta*(pn2 - pm2);

                for (k = 0; k <= num_div_3; ++k)
                {
                    zeta = ((double) k) / num_div_3;
                    p = pq1 + zeta*(pq2-pq1);

                    points.push_back(p);
                }
            }
        }

        // std::cout << "points:" << std::endl;
        // for (std::size_t i = 0; i < points.size(); ++i)
        //     std::cout << "  " << points[i] << std::endl;
        // std::cout << std::endl;

        TIndexType n1, n2, n3, n4, n5, n6, n7, n8;
        for (i = 0; i < num_div_1; ++i)
        {
            for (j = 0; j < num_div_2; ++j)
            {
                for (k = 0; k < num_div_3; ++k)
                {
                    IndexType n1 = starting_node_id + (i * (num_div_2 + 1) + j) * (num_div_3 + 1) + k;
                    IndexType n2 = starting_node_id + (i * (num_div_2 + 1) + j + 1) * (num_div_3 + 1) + k;
                    IndexType n3 = starting_node_id + ((i + 1) * (num_div_2 + 1) + j) * (num_div_3 + 1) + k;
                    IndexType n4 = starting_node_id + ((i + 1) * (num_div_2 + 1) + j + 1) * (num_div_3 + 1) + k;
                    IndexType n5 = n1 + 1;
                    IndexType n6 = n2 + 1;
                    IndexType n7 = n3 + 1;
                    IndexType n8 = n4 + 1;

                    connectivities.push_back(std::vector<std::size_t>{n1, n2, n4, n3, n5, n6, n8, n7});
                }
            }
        }

        // std::cout << "connectivities:" << std::endl;
        // for (std::size_t i = 0; i < connectivities.size(); ++i)
        // {
        //     std::cout << " ";
        //     for (std::size_t j = 0; j < connectivities[i].size(); ++j)
        //         std::cout << " " << connectivities[i][j];
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;

        return std::make_pair(points, connectivities);
    }

    /// Refine a triangle grid by sub-divide a triangle into 4 sub-triangles.
    template<typename TIndexType = std::size_t,
        typename TCoordinatesType = std::vector<double>,
        typename TCoordinatesListType = std::vector<TCoordinatesType>,
        typename TConnectivityType = std::vector<TIndexType>,
        typename TConnectivityListType = std::vector<TConnectivityType> >
    static void RefineTriangleGrid(TCoordinatesListType& Points, TConnectivityListType& Connectivities)
    {
        std::size_t npoints = Points.size();
        TIndexType last_id = static_cast<TIndexType>(npoints-1);

        // generate the new middle points
        typedef std::pair<TIndexType, TIndexType> key_t;
        std::map<key_t, TIndexType> map_corner_to_middle;
        key_t key1, key2;
        TIndexType n1, n2, n3;
        for (typename TConnectivityListType::iterator it = Connectivities.begin(); it != Connectivities.end(); ++it)
        {
            n1 = (*it)[0];
            n2 = (*it)[1];
            n3 = (*it)[2];

            key1 = std::make_pair(n1, n2);
            key2 = std::make_pair(n2, n1);
            if (map_corner_to_middle.find(key1) == map_corner_to_middle.end())
            {
                Points.push_back(0.5*(Points[n1] + Points[n2]));
                map_corner_to_middle[key1] = ++last_id;
                map_corner_to_middle[key2] = last_id;
            }

            key1 = std::make_pair(n2, n3);
            key2 = std::make_pair(n3, n2);
            if (map_corner_to_middle.find(key1) == map_corner_to_middle.end())
            {
                Points.push_back(0.5*(Points[n2] + Points[n3]));
                map_corner_to_middle[key1] = ++last_id;
                map_corner_to_middle[key2] = last_id;
            }

            key1 = std::make_pair(n3, n1);
            key2 = std::make_pair(n1, n3);
            if (map_corner_to_middle.find(key1) == map_corner_to_middle.end())
            {
                Points.push_back(0.5*(Points[n3] + Points[n1]));
                map_corner_to_middle[key1] = ++last_id;
                map_corner_to_middle[key2] = last_id;
            }
        }

        // generate new triangles
        TIndexType m1, m2, m3;
        TConnectivityListType Connectivities_old = Connectivities;
        Connectivities.clear();
        for (typename TConnectivityListType::iterator it = Connectivities_old.begin(); it != Connectivities_old.end(); ++it)
        {
            n1 = (*it)[0];
            n2 = (*it)[1];
            n3 = (*it)[2];

            m1 = map_corner_to_middle[std::make_pair(n1, n2)];
            m2 = map_corner_to_middle[std::make_pair(n2, n3)];
            m3 = map_corner_to_middle[std::make_pair(n3, n1)];

            Connectivities.push_back(TConnectivityType{n1, m1, m3});
            Connectivities.push_back(TConnectivityType{m1, n2, m2});
            Connectivities.push_back(TConnectivityType{m1, m2, m3});
            Connectivities.push_back(TConnectivityType{m2, n3, m3});
        }
    }

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

    /// Create a list of entities (element/condition) from a model_part to another model_part
    /// It is noted that the newly created entities are not added to the other model_part. User must do it manually.
    template<class TEntityType, class TEntitiesContainerType>
    static TEntitiesContainerType CreateEntities(
        TEntitiesContainerType& pEntities,
        ModelPart& r_other_model_part,
        TEntityType const& r_sample_entity,
        std::size_t& last_entity_id,
        Properties::Pointer pProperties,
        const bool& retain_prop_id = false)
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
            if (!retain_prop_id)
            {
                pNewEntities.push_back(r_sample_entity.Create(++last_entity_id, temp_entity_nodes, pProperties));
            }
            else
            {
                Properties::Pointer pNewProperties = r_other_model_part.pGetProperties((*it)->GetProperties().Id());
                pNewEntities.push_back(r_sample_entity.Create(++last_entity_id, temp_entity_nodes, pNewProperties));
            }
        }

        return pNewEntities;
    }

    // /// Create conditions/elements from the list of points.
    // /// The triangulation will be performed on the physical points with the information {cemter, normal, t1, t2}
    // /// The refinement is performed on the local points instead.
    // /// The point list will be triangulated before the conditions are created.
    // /// It is noted that the newly created entities are not added to the other model_part. User must do it manually. Nevertheless, the nodes are added to the model_part.
    // /// The new node will be created from last_node_id+1
    // template<typename TPointType, typename TVectorType, class TEntityType, class TPointsContainerType, class TNodesContainerType, class TEntitiesContainerType>
    // static std::tuple<TPointsContainerType, TNodesContainerType, TEntitiesContainerType> CreateEntities(
    //     const TPointsContainerType& physical_points,
    //     const TVectorType& rCenter,
    //     const TVectorType& rNormal,
    //     const TVectorType& rTangent1,
    //     const TVectorType& rTangent2,
    //     const TPointsContainerType& local_points,
    //     const std::size_t& nrefine,
    //     ModelPart& r_model_part,
    //     TEntityType const& r_sample_entity,
    //     std::size_t& last_node_id,
    //     std::size_t& last_entity_id,
    //     Properties::Pointer pProperties)
    // {
    //     // compute the triangulation
    //     typedef unsigned int IndexType;
    //     typedef std::vector<std::vector<IndexType> > connectivity_t;
    //     connectivity_t Connectivities = GenerateTriangleGrid<TPointType, TVectorType, IndexType>(physical_points, rCenter, rNormal, rTangent1, rTangent2);

    //     // refine if needed
    //     TPointsContainerType new_points = local_points;
    //     for (std::size_t i = 0; i < nrefine; ++i)
    //         RefineTriangleGrid<unsigned int, TPointType>(new_points, Connectivities);

    //     // offset the connectivity
    //     for (std::size_t i = 0; i < Connectivities.size(); ++i)
    //         for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
    //             Connectivities[i][j] += last_node_id+1;

    //     // std::cout << "Connectivities:" << std::endl;
    //     // for (std::size_t i = 0; i < Connectivities.size(); ++i)
    //     // {
    //     //     std::cout << "  " << i << ":";
    //     //     for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
    //     //         std::cout << " " << Connectivities[i][j];
    //     //     std::cout << std::endl;
    //     // }

    //     // create the nodes
    //     std::vector<std::size_t> map_con_to_mp(new_points.size());
    //     TNodesContainerType pNewNodes;
    //     for (std::size_t i = 0; i < new_points.size(); ++i)
    //     {
    //         NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id, new_points[i][0], new_points[i][1], new_points[i][2]);
    //         map_con_to_mp[i] = pNewNode->Id();
    //         pNewNodes.push_back(pNewNode);
    //     }

    //     // create the entities based on connectivity
    //     const std::string NodeKey = std::string("Node");
    //     TEntitiesContainerType pNewEntities = CreateEntities<connectivity_t, TEntityType, TEntitiesContainerType>(Connectivities, r_model_part,
    //         r_sample_entity, last_entity_id, pProperties, NodeKey);

    //     return std::make_tuple(new_points, pNewNodes, pNewEntities);
    // }

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

#ifdef USE_TRIANGULATION_UTILS_FOR_TRIANGULATION
#undef USE_TRIANGULATION_UTILS_FOR_TRIANGULATION
#endif

#ifdef USE_CGAL_FOR_TRIANGULATION
#undef USE_CGAL_FOR_TRIANGULATION
#endif

#endif // KRATOS_ISOGEOMETRIC_POST_UTILITY_H_INCLUDED

