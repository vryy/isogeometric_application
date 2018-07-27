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
#include "custom_utilities/tsplines/tsedge.h"
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

    typedef std::vector<knot_t>             knot_container_t;
    typedef std::list<Cell::Pointer>        cell_container_t;
    typedef std::list<TsAnchor::Pointer>    anchor_container_t;
    typedef std::list<TsVertex::Pointer>    vertex_container_t;
    typedef std::list<TsEdge::Pointer>      edge_container_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(TsMesh2D);

    /// Default constructor
    TsMesh2D();

    /// Destructor
    ~TsMesh2D();

    /// Subroutines to modify the T-splines mesh
    void BeginConstruct();
    void SetOrder(const int& dim, const int& order);
    knot_t InsertKnot(const int& dim, const double& knot);
    TsVertex::Pointer AddVertex(knot_t pXi, knot_t pEta);
    TsEdge::Pointer AddHEdge(TsVertex::Pointer pV1, TsVertex::Pointer pV2);
    TsEdge::Pointer AddVEdge(TsVertex::Pointer pV1, TsVertex::Pointer pV2);
    void EndConstruct();

    /// Subroutines to query the T-splines mesh
    int Order(const int& dim) const;
    std::size_t NumberOfKnots(const int& dim) const;
    knot_t GetKnot(const int& dim, const std::size_t& index) const;
    const edge_container_t& Edges() const;
    const anchor_container_t& Anchors() const;
    const cell_container_t& Cells() const;
    void FindCells(std::set<cell_t>& rCells, bool _extend = false) const;
    void FindAnchors(std::vector<anchor_t>& rAnchors) const;
    bool IsAnalysisSuitable();

    /// Subroutines to modify the T-splines mesh
    void RenumberMesh();

    /// Auxilliary subroutines
    void ClearExtendedTmesh();
    void BuildExtendedTmesh();
    void BuildAnchors(std::string fn);
    void BuildCells();

    /// Print out
    void PrintInfo(std::ostream& rOStream) const;

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
    void FindKnots(const double& Anchor_xi_index, const double& Anchor_eta_index,
    	std::vector<DataType>& Knots1, std::vector<DataType>& Knots2) const
    {
        std::set<std::size_t> tmp_knot_index_left;
        std::set<std::size_t> tmp_knot_index_right;
        std::set<std::size_t> tmp_knot_index_up;
        std::set<std::size_t> tmp_knot_index_down;

        // marching to the all directions and find the intersecting edges
        for(edge_container_t::const_iterator it = mEdges.begin(); it != mEdges.end(); ++it)
        {
            if((*it)->EdgeType() == TsEdge::VERTICAL_EDGE) //vertical edge
            {
                std::size_t edge_xi_index = (*it)->Index();
                if((*it)->IsCut(Anchor_eta_index) && edge_xi_index < Anchor_xi_index)
                    tmp_knot_index_left.insert(edge_xi_index);
                if((*it)->IsCut(Anchor_eta_index) && edge_xi_index > Anchor_xi_index)
                    tmp_knot_index_right.insert(edge_xi_index);
            }
            if((*it)->EdgeType() == TsEdge::HORIZONTAL_EDGE) //horizontal edge
            {
                std::size_t edge_eta_index = (*it)->Index();
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
        if(this->mOrder[0] % 2 == 0)
        {
            std::size_t span = this->mOrder[0]/2 + 1;
            int k_index;
            std::vector<std::size_t> tmp_left(tmp_knot_index_left.begin(), tmp_knot_index_left.end());
            std::vector<std::size_t> tmp_right(tmp_knot_index_right.begin(), tmp_knot_index_right.end());

            if(Knots1.size() != 2*span)
                Knots1.resize(2*span);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_left.end() - span + i);
                if(FuncType == 1)
                    Knots1[i] = static_cast<DataType>(mKnots[0][k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i] = static_cast<DataType>(mKnots[0][k_index]->Index());
            }
            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_right.begin() + i);
                if(FuncType == 1)
                    Knots1[i + span] = static_cast<DataType>(mKnots[0][k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i + span] = static_cast<DataType>(mKnots[0][k_index]->Index());
            }
        }
        else
        {
            std::size_t span = (this->mOrder[0] + 1)/2;
            int k_index;
            std::vector<std::size_t> tmp_left(tmp_knot_index_left.begin(), tmp_knot_index_left.end());
            std::vector<std::size_t> tmp_right(tmp_knot_index_right.begin(), tmp_knot_index_right.end());

            if(Knots1.size() != 2*span + 1)
                Knots1.resize(2*span + 1);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_left.end() - span + i);
                if(FuncType == 1)
                    Knots1[i] = static_cast<DataType>(mKnots[0][k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i] = static_cast<DataType>(mKnots[0][k_index]->Index());
            }

            if(FuncType == 1)
                Knots1[span] = static_cast<DataType>(mKnots[0][Anchor_xi_index]->Value());
            else if(FuncType == 2)
                Knots1[span] = static_cast<DataType>(mKnots[0][Anchor_xi_index]->Index());

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_right.begin() + i);
                if(FuncType == 1)
                    Knots1[i + span + 1] = static_cast<DataType>(mKnots[0][k_index]->Value());
                else if(FuncType == 2)
                    Knots1[i + span + 1] = static_cast<DataType>(mKnots[0][k_index]->Index());
            }
        }

        if(this->mOrder[1] % 2 == 0)
        {
            std::size_t span = this->mOrder[1]/2 + 1;
            int k_index;
            std::vector<std::size_t> tmp_down(tmp_knot_index_down.begin(), tmp_knot_index_down.end());
            std::vector<std::size_t> tmp_up(tmp_knot_index_up.begin(), tmp_knot_index_up.end());

            if(Knots2.size() != 2*span)
                Knots2.resize(2*span);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_down.end() - span + i);
                if(FuncType == 1)
                    Knots2[i] = static_cast<DataType>(mKnots[1][k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i] = static_cast<DataType>(mKnots[1][k_index]->Index());
            }
            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_up.begin() + i);
                if(FuncType == 1)
                    Knots2[i + span] = static_cast<DataType>(mKnots[1][k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i + span] = static_cast<DataType>(mKnots[1][k_index]->Index());
            }
        }
        else
        {
            std::size_t span = (this->mOrder[1] + 1)/2;
            int k_index;
            std::vector<std::size_t> tmp_down(tmp_knot_index_down.begin(), tmp_knot_index_down.end());
            std::vector<std::size_t> tmp_up(tmp_knot_index_up.begin(), tmp_knot_index_up.end());

            if(Knots2.size() != 2*span + 1)
                Knots2.resize(2*span + 1);

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_down.end() - span + i);
                if(FuncType == 1)
                    Knots2[i] = static_cast<DataType>(mKnots[1][k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i] = static_cast<DataType>(mKnots[1][k_index]->Index());
            }

            if(FuncType == 1)
                Knots2[span] = static_cast<DataType>(mKnots[1][Anchor_eta_index]->Value());
            else if(FuncType == 2)
                Knots2[span] = static_cast<DataType>(mKnots[1][Anchor_eta_index]->Index());

            for(std::size_t i = 0; i < span; ++i)
            {
                k_index = *(tmp_up.begin() + i);
                if(FuncType == 1)
                    Knots2[i + span + 1] = static_cast<DataType>(mKnots[1][k_index]->Value());
                else if(FuncType == 2)
                    Knots2[i + span + 1] = static_cast<DataType>(mKnots[1][k_index]->Index());
            }
        }
    }

private:
    vertex_container_t mVertices; // list of vertices
    vertex_container_t mVirtualVertices; // list of virtual vertices
    edge_container_t mEdges; // list of edges
    cell_container_t mCells; // list of cells
    anchor_container_t mAnchors; // list of anchors

    boost::array<int, 2> mOrder; // order of the Tsplines mesh in horizontal and vertical direction

    std::size_t mLastVertex; // internal variable point to the last vertex identification in the T-splines mesh
    std::size_t mLastEdge; // internal variable point to the last edge identification in the T-splines mesh

    boost::array<double, 2> mKnotsMin;
    boost::array<double, 2> mKnotsMax;

    boost::array<knot_container_t, 2> mKnots; // 0: knot vector in horizontal direction
                                              // 1: knot vector in vertical direction

    bool mLockConstruct; // lock variable to control the build process
    bool mIsExtended; // variable to keep track with the construction of extended topology mesh

    void LockQuery()
    {
        if(mLockConstruct)
            KRATOS_THROW_ERROR(std::logic_error, "The T-splines mesh is currently locked. Please call BeginConstruct() to unlock", "")
    }

    /// For debugging only
    void FindKnots2(const double& Anchor_xi_index, const double& Anchor_eta_index,
    	Vector& Knots1, Vector& Knots2) const
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

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const TsMesh2D& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_TS_MESH_2D_H_INCLUDED

