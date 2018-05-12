//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 4 Feb 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_TS_MESH_2D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_TS_MESH_2D_H_INCLUDED

// System includes
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

// External includes
#include <omp.h>
#include "boost/progress.hpp"
#include "boost/algorithm/string.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/nurbs/cell.h"
#include "custom_utilities/tsplines/tsedges.h"
#include "custom_utilities/tsplines/tsanchor.h"

namespace Kratos
{

/**
    Class represents a Tsplines mesh in 2D
    TODO:
    +   add subroutines to identify the list of anchors in the Tsplines mesh (hard)
        -   simple for odd order in u & v, since anchors and vertex are identical
        -   non-trivial if one of order in u or v is even, or both.
            According to Fig. 1, ANALYSIS-SUITABLE T-SPLINES OF ARBITRARY DEGREE: DEFINITION, LINEAR INDEPENDENCE AND APPROXIMATION PROPERTIES, Veiga et al. For p even, q odd, the anchors are the mids of horizontal edges. For p odd, q even, the anchors are the mids of the vertical edges. For p odd, q odd, the anchors are the mid of the cell => need to identify all atomic cells in the T-splines topology mesh.
            Note that anchors are always located in the active region (i.e excluding the repetitive boundary)
            Question on p even, q even if T-splines mesh contain L-junctions.
        -> need to add routines to identify all cells in the T-mesh


    +   add subroutines to check the validation of the Tsplines mesh (hard)
        -   need to identify all the cells within the topology mesh

    +   add subroutines to build the extended topology mesh (hard)
        -> need to add vertices/edges

    +   add subroutines to check for the analysis-suitable T-splines (medium)

*/
class TsMesh2D
{
public:
    /// Const definition
    static const int NO_READ      = 0;
    static const int READ_ORDER   = 1;
    static const int READ_KNOTS   = 2;
    static const int READ_H_EDGES = 3;
    static const int READ_V_EDGES = 4;
    static const int READ_ANCHORS = 5;

    /// Type definition
    typedef std::pair<std::pair<int, int>, std::pair<int, int> > cell_t;
    typedef std::pair<double, double>       anchor_t;
    typedef Knot<double>::Pointer           knot_t;

    typedef std::vector<knot_t>      knot_container_t;
    typedef std::list<Cell::Pointer>        cell_container_t;
    typedef std::list<TsAnchor::Pointer>    anchor_container_t;
    typedef std::list<TsVertex::Pointer>    vertex_container_t;
    typedef std::list<TsEdge::Pointer>      edge_container_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsMesh2D);

    /// Default constructor
    TsMesh2D() : mOrder1(1), mOrder2(1), mLastEdge(0), mLastVertex(0), mLockConstruct(true), mIsExtended(false) {}

    /// Destructor
    ~TsMesh2D() {}


    /// Subroutines to modify the T-splines mesh
    void BeginConstruct();
    void SetOrder(int Dim, int Order);
    knot_t InsertKnot(int Dim, double Knot);
    TsVertex::Pointer AddVertex(knot_t pXi, knot_t pEta);
    TsEdge::Pointer AddHEdge(TsVertex::Pointer pV1, TsVertex::Pointer pV2);
    TsEdge::Pointer AddVEdge(TsVertex::Pointer pV1, TsVertex::Pointer pV2);
    void ReadFromFile(std::string fn);
    void EndConstruct();

    /// Subroutines to query the T-splines mesh
    void FindCells(std::set<cell_t>& rCells, bool _extend = false) const;
    bool IsAnalysisSuitable();

    /// Subroutines to modify the T-splines mesh
    void RenumberMesh();

    /// Auxilliary subroutines
    void ClearExtendedTmesh();
    void BuildExtendedTmesh();
    void BuildAnchors(std::string fn);
    void BuildCells();
    void PrintInfo(std::ostream& rOStream) const;
    void ExportMatlab(std::string fn, std::string mesh_type) const;
    void ExportMDPA(std::string fn, int Division1, int Division2) const;

private:
    vertex_container_t mVertices; // list of vertices
    vertex_container_t mVirtualVertices; // list of virtual vertices
    edge_container_t mEdges; // list of edges
    cell_container_t mCells; // list of cells
    anchor_container_t mAnchors; // list of anchors

    int mOrder1, mOrder2; // order of the Tsplines mesh in horizontal and vertical direction
    int mLastVertex; // internal variable point to the last vertex identification in the T-splines mesh
    int mLastEdge; // internal variable point to the last edge identification in the T-splines mesh
    double mKnots1Min;
    double mKnots1Max;
    double mKnots2Min;
    double mKnots2Max;

    knot_container_t mKnots1; // knot vector in horizontal direction
    knot_container_t mKnots2; // knot vector in vertical direction

    bool mLockConstruct; // lock variable to control the build process
    bool mIsExtended; // variable to keep track with the construction of extended topology mesh

    void LockQuery()
    {
        if(mLockConstruct)
            KRATOS_THROW_ERROR(std::logic_error, "The T-splines mesh is currently locked. Please call BeginConstruct() to unlock", "")
    }

    /// Get the list of anchors associated with the T-splines topology mesh
    /// Remarks: this is the anchors in the topology coordinates, not the anchors in knot coordinates
    void FindAnchors(std::vector<anchor_t>& rAnchors) const
    {
        double anchor_xi;
        double anchor_eta;
        if((mOrder1 % 2 != 0) && (mOrder2 % 2 != 0))
        {
            // for odd order T-splines topology mesh, the vertex is also the anchor
            for(vertex_container_t::const_iterator it = mVertices.begin(); it != mVertices.end(); ++it)
                if((*it)->pXi()->IsActive() && (*it)->pEta()->IsActive())
                {
                    anchor_xi = static_cast<double>((*it)->Index1());
                    anchor_eta = static_cast<double>((*it)->Index2());
                    rAnchors.push_back(anchor_t(anchor_xi, anchor_eta));
                }
        }
        else if((mOrder1 % 2 == 0) && (mOrder2 % 2 != 0))
        {
            // the anchors are the middle of all horizontal edges
            for(edge_container_t::const_iterator it = mEdges.begin(); it != mEdges.end(); ++it)
                if((*it)->EdgeType() == TsEdge::HORIZONTAL_EDGE && (*it)->IsActive())
                {
                    anchor_xi = 0.5 * static_cast<double>((*it)->pV1()->Index1() + (*it)->pV2()->Index1());
                    anchor_eta = static_cast<double>((*it)->Index());
                    rAnchors.push_back(anchor_t(anchor_xi, anchor_eta));
                }
        }
        else if((mOrder1 % 2 != 0) && (mOrder2 % 2 == 0))
        {
            // the anchors are the middle of all vertical edges
            for(edge_container_t::const_iterator it = mEdges.begin(); it != mEdges.end(); ++it)
                if((*it)->EdgeType() == TsEdge::VERTICAL_EDGE && (*it)->IsActive())
                {
                    anchor_xi = static_cast<double>((*it)->Index());
                    anchor_eta = 0.5 * static_cast<double>((*it)->pV1()->Index2() + (*it)->pV2()->Index2());
                    rAnchors.push_back(anchor_t(anchor_xi, anchor_eta));
                }
        }
        else if((mOrder1 % 2 == 0) && (mOrder2 % 2 == 0))
        {
            // the anchors are the middle of the cells
            std::set<cell_t> cells;
            this->FindCells(cells, false);
            for(std::set<cell_t>::const_iterator it = cells.begin(); it != cells.end(); ++it)
            {
                anchor_xi = 0.5 * (mKnots1[it->first.first]->Value() + mKnots1[it->first.second]->Value());
                anchor_eta = 0.5 * (mKnots2[it->second.first]->Value() + mKnots2[it->second.second]->Value());
                rAnchors.push_back(anchor_t(anchor_xi, anchor_eta));
            }
        }
    }

    /// Find the local knot vectors for an arbitrary anchor
    /// Algorithm: ray marching, Isogeometric analysis using T-splines
    /// Id1, Id2: topology coordinates of the anchor (IN)
    ///     In the case that Order is odd, the topology coordinate is at the vertex, so it will be an integer
    ///     In the case that Order is even, the topology coordinate is in the middle of the edge, so it will be a double
    /// Knots1, Knots2: knot vectors in horizontal and vertical direction
    /// Usage:
    ///     Call FindKnots<1, double> if one wants to find the local knot vector associated with the anchor
    ///     Call FindKnots<2, int> if one wants to find the index in topology coordinates of local knot vector associated with the anchor
    template<int FuncType, class DataType>
    void FindKnots(double Anchor_xi_index, double Anchor_eta_index, std::vector<DataType>& Knots1, std::vector<DataType>& Knots2) const
    {
        std::set<int> tmp_knot_index_left;
        std::set<int> tmp_knot_index_right;
        std::set<int> tmp_knot_index_up;
        std::set<int> tmp_knot_index_down;

        // marching to the all directions and find the intersecting edges
        for(edge_container_t::const_iterator it = mEdges.begin(); it != mEdges.end(); ++it)
        {
            if((*it)->EdgeType() == TsEdge::VERTICAL_EDGE) //vertical edge
            {
                int edge_xi_index = (*it)->Index();
                if((*it)->IsCut(Anchor_eta_index) && edge_xi_index < Anchor_xi_index)
                    tmp_knot_index_left.insert(edge_xi_index);
                if((*it)->IsCut(Anchor_eta_index) && edge_xi_index > Anchor_xi_index)
                    tmp_knot_index_right.insert(edge_xi_index);
            }
            if((*it)->EdgeType() == TsEdge::HORIZONTAL_EDGE) //horizontal edge
            {
                int edge_eta_index = (*it)->Index();
                if((*it)->IsCut(Anchor_xi_index) && edge_eta_index < Anchor_eta_index)
                    tmp_knot_index_down.insert(edge_eta_index);
                if((*it)->IsCut(Anchor_xi_index) && edge_eta_index > Anchor_eta_index)
                    tmp_knot_index_up.insert(edge_eta_index);
            }
        }
//        std::cout << "marching completed" << std::endl;

//        std::cout << "tmp_knot_index_left:";
//        for(std::set<int>::iterator it = tmp_knot_index_left.begin(); it != tmp_knot_index_left.end(); ++it)
//            std::cout << " " << (*it);
//        std::cout << std::endl;
//
//        std::cout << "tmp_knot_index_right:";
//        for(std::set<int>::iterator it = tmp_knot_index_right.begin(); it != tmp_knot_index_right.end(); ++it)
//            std::cout << " " << (*it);
//        std::cout << std::endl;
//
//        std::cout << "tmp_knot_index_up:";
//        for(std::set<int>::iterator it = tmp_knot_index_up.begin(); it != tmp_knot_index_up.end(); ++it)
//            std::cout << " " << (*it);
//        std::cout << std::endl;
//
//        std::cout << "tmp_knot_index_down:";
//        for(std::set<int>::iterator it = tmp_knot_index_down.begin(); it != tmp_knot_index_down.end(); ++it)
//            std::cout << " " << (*it);
//        std::cout << std::endl;

        // fill in the knot vectors
        if(mOrder1 % 2 == 0)
        {
            int span = mOrder1/2 + 1;
            int k_index;
            std::vector<int> tmp_left(tmp_knot_index_left.begin(), tmp_knot_index_left.end());
            std::vector<int> tmp_right(tmp_knot_index_right.begin(), tmp_knot_index_right.end());

            if(Knots1.size() != 2*span)
                Knots1.resize(2*span);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_left.end() - span + i);
                if(FuncType == 1)
                    Knots1[i] = static_cast<DataType>(mKnots1[k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i] = static_cast<DataType>(mKnots1[k_index]->Index());
            }
            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_right.begin() + i);
                if(FuncType == 1)
                    Knots1[i + span] = static_cast<DataType>(mKnots1[k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i + span] = static_cast<DataType>(mKnots1[k_index]->Index());
            }
        }
        else
        {
            int span = (mOrder1 + 1)/2;
            int k_index;
            std::vector<int> tmp_left(tmp_knot_index_left.begin(), tmp_knot_index_left.end());
            std::vector<int> tmp_right(tmp_knot_index_right.begin(), tmp_knot_index_right.end());

            if(Knots1.size() != 2*span + 1)
                Knots1.resize(2*span + 1);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_left.end() - span + i);
                if(FuncType == 1)
                    Knots1[i] = static_cast<DataType>(mKnots1[k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i] = static_cast<DataType>(mKnots1[k_index]->Index());
            }

            if(FuncType == 1)
                Knots1[span] = static_cast<DataType>(mKnots1[Anchor_xi_index]->Value());
            else if(FuncType == 2)
                Knots1[span] = static_cast<DataType>(mKnots1[Anchor_xi_index]->Index());

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_right.begin() + i);
                if(FuncType == 1)
                    Knots1[i + span + 1] = static_cast<DataType>(mKnots1[k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i + span + 1] = static_cast<DataType>(mKnots1[k_index]->Index());
            }
        }

        if(mOrder2 % 2 == 0)
        {
            int span = mOrder2/2 + 1;
            int k_index;
            std::vector<int> tmp_down(tmp_knot_index_down.begin(), tmp_knot_index_down.end());
            std::vector<int> tmp_up(tmp_knot_index_up.begin(), tmp_knot_index_up.end());

            if(Knots2.size() != 2*span)
                Knots2.resize(2*span);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_down.end() - span + i);
                if(FuncType == 1)
                    Knots2[i] = static_cast<DataType>(mKnots2[k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i] = static_cast<DataType>(mKnots2[k_index]->Index());
            }
            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_up.begin() + i);
                if(FuncType == 1)
                    Knots2[i + span] = static_cast<DataType>(mKnots2[k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i + span] = static_cast<DataType>(mKnots2[k_index]->Index());
            }
        }
        else
        {
            int span = (mOrder2 + 1)/2;
            int k_index;
            std::vector<int> tmp_down(tmp_knot_index_down.begin(), tmp_knot_index_down.end());
            std::vector<int> tmp_up(tmp_knot_index_up.begin(), tmp_knot_index_up.end());

            if(Knots2.size() != 2*span + 1)
                Knots2.resize(2*span + 1);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_down.end() - span + i);
                if(FuncType == 1)
                    Knots2[i] = static_cast<DataType>(mKnots2[k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i] = static_cast<DataType>(mKnots2[k_index]->Index());
            }

            if(FuncType == 1)
                Knots2[span] = static_cast<DataType>(mKnots2[Anchor_eta_index]->Value());
            else if(FuncType == 2)
                Knots2[span] = static_cast<DataType>(mKnots2[Anchor_eta_index]->Index());

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_up.begin() + i);
                if(FuncType == 1)
                    Knots2[i + span + 1] = static_cast<DataType>(mKnots2[k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i + span + 1] = static_cast<DataType>(mKnots2[k_index]->Index());
            }
        }
    }

    /// For debugging only
    void FindKnots2(double Anchor_xi_index, double Anchor_eta_index, Vector& Knots1, Vector& Knots2) const
    {
        std::vector<double> tmpKnots1;
        std::vector<double> tmpKnots2;
        this->FindKnots<1, double>(Anchor_xi_index, Anchor_eta_index, tmpKnots1, tmpKnots2);

        if(Knots1.size() != tmpKnots1.size())
            Knots1.resize(tmpKnots1.size());
        std::copy(tmpKnots1.begin(), tmpKnots1.end(), Knots1.begin());

        if(Knots2.size() != tmpKnots2.size())
            Knots2.resize(tmpKnots2.size());
        std::copy(tmpKnots2.begin(), tmpKnots2.end(), Knots2.begin());
    }

    /// Find the span of knot in the local knot vector
    /// Remarks: it will give the based-1 index
    int FindSpanLocal(double Xi, const std::vector<double>& U)
    {
//        int low = 0;
//        int high = U.size() - 1;
//        int mid = (low + high) / 2;
//
//        while( Xi < U[mid] || Xi >= U[mid+1] )
//        {
//            if(Xi < U[mid])
//            {
//                high = mid;
//            }
//            else
//            {
//                low = mid;
//            }
//            mid = (low + high) / 2;
//        }
//
//        return mid + 1;

        if(!U.empty())
        {
            if(Xi < U.front())
                return 0;

            if(Xi > U.back())
                return U.size();

            for(std::size_t i = 0; i < U.size() - 1; ++i)
                if(Xi >= U[i] && Xi < U[i + 1])
                    return i + 1;
        }

        return 0;
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsMesh2D& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_MESH_2D_H_INCLUDED

