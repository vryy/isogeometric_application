//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Jul 2018 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_TSPLINE_UTILS_H_INCLUDED )
#define  KRATOS_TSPLINE_UTILS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <ctime>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/tsplines/tsmesh_2d.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"


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
 * Utility for various operations on Tsplines
 */
class TSplineUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TSplineUtils
    KRATOS_CLASS_POINTER_DEFINITION(TSplineUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TSplineUtils()
    {}

    /// Destructor.
    virtual ~TSplineUtils()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create the T-mesh from BSplines
    static void CreateFromBSplines(TsMesh2D& tmesh, const BSplinesFESpace<2>& rFESpace)
    {
        // TODO
    }

    /// Read the T-mesh topology from .tmesh file
    static void ReadFromFile(TsMesh2D& tmesh, const std::string& fn)
    {
        tmesh.BeginConstruct();

        std::ifstream infile(fn.c_str());

        std::string line;
        std::vector<std::string> words;
        int ReadMode = TsMesh2D::NO_READ;
        int dim = 0;
        std::vector<std::pair<int, std::pair<int, int> > > HEdgeList;
        std::vector<std::pair<int, std::pair<int, int> > > VEdgeList;
        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

//            for(std::size_t i = 0; i < words.size(); ++i)
//                std::cout << " " << words[i];
//            std::cout << std::endl;

            if(words.size() != 0)
            {
                if(words[0] == std::string("Begin"))
                {
                    if(words.size() < 2)
                        KRATOS_THROW_ERROR(std::logic_error, "Missing statement for Begin", "")
                    if(words[1] == "Order")
                        ReadMode = TsMesh2D::READ_ORDER;
                    else if(words[1] == std::string("Knots"))
                        ReadMode = TsMesh2D::READ_KNOTS;
                    else if(words[1] == std::string("H-edges"))
                        ReadMode = TsMesh2D::READ_H_EDGES;
                    else if(words[1] == std::string("V-edges"))
                        ReadMode = TsMesh2D::READ_V_EDGES;
                    continue;
                }

                if(words[0] == std::string("End"))
                {
                    dim = 0;
                    ReadMode = TsMesh2D::NO_READ;
                    continue;
                }

                if(ReadMode == TsMesh2D::READ_ORDER)
                {
                    int order = atoi(words[0].c_str());
                    tmesh.SetOrder(dim, order);
                    ++dim;
                }
                else if(ReadMode == TsMesh2D::READ_KNOTS)
                {
                    for(std::size_t i = 0; i < words.size(); ++i)
                        tmesh.InsertKnot(dim, atof(words[i].c_str()));
                    ++dim;
                }
                else if(ReadMode == TsMesh2D::READ_H_EDGES)
                {
                    int v  = atoi(words[0].c_str());
                    int u1 = atoi(words[1].c_str());
                    int u2 = atoi(words[2].c_str());
                    HEdgeList.push_back(std::pair<int, std::pair<int, int> >(v, std::pair<int, int>(u1, u2)));
                }
                else if(ReadMode == TsMesh2D::READ_V_EDGES)
                {
                    int u  = atoi(words[0].c_str());
                    int v1 = atoi(words[1].c_str());
                    int v2 = atoi(words[2].c_str());
                    VEdgeList.push_back(std::pair<int, std::pair<int, int> >(u, std::pair<int, int>(v1, v2)));
                }
            }
        }

        infile.close();

        // assign the index to the knot vectors
        int cnt = -1;
        for(std::size_t i = 0; i < tmesh.NumberOfKnots(0); ++i)
            tmesh.GetKnot(0, i)->UpdateIndex(++cnt);
        cnt = -1;
        for(std::size_t i = 0; i < tmesh.NumberOfKnots(1); ++i)
            tmesh.GetKnot(1, i)->UpdateIndex(++cnt);

        // create the vertex list
        std::set<std::pair<int, int> > VertexList;
        for(std::size_t i = 0; i < HEdgeList.size(); ++i)
        {
            VertexList.insert(std::pair<int, int>(HEdgeList[i].second.first, HEdgeList[i].first));
            VertexList.insert(std::pair<int, int>(HEdgeList[i].second.second, HEdgeList[i].first));
        }
        for(std::size_t i = 0; i < VEdgeList.size(); ++i)
        {
            VertexList.insert(std::pair<int, int>(VEdgeList[i].first, VEdgeList[i].second.first));
            VertexList.insert(std::pair<int, int>(VEdgeList[i].first, VEdgeList[i].second.second));
        }

        // insert the vertices to the T-splines mesh
        std::map<std::pair<int, int>, TsVertex::Pointer> VertexMap;
        for(std::set<std::pair<int, int> >::iterator it = VertexList.begin(); it != VertexList.end(); ++it)
        {
            TsVertex::Pointer pV = tmesh.AddVertex(tmesh.GetKnot(0, it->first), tmesh.GetKnot(1, it->second));
            VertexMap[std::pair<int, int>(it->first, it->second)] = pV;
        }

        // insert the edges to the T-splines mesh
        for(std::size_t i = 0; i < HEdgeList.size(); ++i)
        {
            TsVertex::Pointer pV1 = VertexMap[std::pair<int, int>(HEdgeList[i].second.first, HEdgeList[i].first)];
            TsVertex::Pointer pV2 = VertexMap[std::pair<int, int>(HEdgeList[i].second.second, HEdgeList[i].first)];
            tmesh.AddHEdge(pV1, pV2);
        }
        for(std::size_t i = 0; i < VEdgeList.size(); ++i)
        {
            TsVertex::Pointer pV1 = VertexMap[std::pair<int, int>(VEdgeList[i].first, VEdgeList[i].second.first)];
            TsVertex::Pointer pV2 = VertexMap[std::pair<int, int>(VEdgeList[i].first, VEdgeList[i].second.second)];
            tmesh.AddVEdge(pV1, pV2);
        }

        tmesh.EndConstruct();
    }

    /// Export the topology mesh/knot coordinates mesh to matlab
    static void ExportMatlab(const TsMesh2D& tmesh,
        const std::string& fn, const std::string& mesh_type)
    {
        typedef TsMesh2D::edge_container_t edge_container_t;
        typedef TsMesh2D::cell_t cell_t;
        typedef TsMesh2D::anchor_t anchor_t;

        std::ofstream outfile(fn.c_str());

        outfile << "axis equal" << std::endl;
        outfile << "close all" << std::endl;
        outfile << "hold on" << std::endl << std::endl;

        // plot edges
        if(mesh_type == std::string("topology"))
        {
            for(edge_container_t::const_iterator it = tmesh.Edges().begin(); it != tmesh.Edges().end(); ++it)
            {
                if((*it)->EdgeType() == TsEdge::VIRTUAL_HORIZONTAL_EDGE || (*it)->EdgeType() == TsEdge::VIRTUAL_VERTICAL_EDGE)
                {
                    outfile << "line([" << (*it)->pV1()->Index1() << " " << (*it)->pV2()->Index1() << "],";
                    outfile << "[" << (*it)->pV1()->Index2() << " " << (*it)->pV2()->Index2() << "],'LineStyle',':');" << std::endl;
                }
                else
                {
                    outfile << "line([" << (*it)->pV1()->Index1() << " " << (*it)->pV2()->Index1() << "],";
                    outfile << "[" << (*it)->pV1()->Index2() << " " << (*it)->pV2()->Index2() << "]);" << std::endl;
                }
            }
        }
        else if(mesh_type == std::string("knots"))
        {
            for(edge_container_t::const_iterator it = tmesh.Edges().begin(); it != tmesh.Edges().end(); ++it)
            {
                if((*it)->EdgeType() == TsEdge::VIRTUAL_HORIZONTAL_EDGE || (*it)->EdgeType() == TsEdge::VIRTUAL_VERTICAL_EDGE)
                {
                    outfile << "line([" << (*it)->pV1()->pXi()->Value() << " " << (*it)->pV2()->pXi()->Value() << "],";
                    outfile << "[" << (*it)->pV1()->pEta()->Value() << " " << (*it)->pV2()->pEta()->Value() << "],'LineStyle',':');" << std::endl;
                }
                else
                {
                    outfile << "line([" << (*it)->pV1()->pXi()->Value() << " " << (*it)->pV2()->pXi()->Value() << "],";
                    outfile << "[" << (*it)->pV1()->pEta()->Value() << " " << (*it)->pV2()->pEta()->Value() << "]);" << std::endl;
                }
            }
        }
        outfile << std::endl;

        // find all anchors in the current topology mesh
        std::vector<anchor_t> Anchors;
        tmesh.FindAnchors(Anchors);

        // export the knot vectors for each anchors
        std::vector<double> Knots1;
        std::vector<double> Knots2;
        int cnt = 0;
        for(std::size_t i = 0; i < Anchors.size(); ++i)
        {
            tmesh.FindKnots<1, double>(Anchors[i].first, Anchors[i].second, Knots1, Knots2);
            outfile << "local_knots(" << ++cnt << ",:,:) = [";
            for(std::size_t i = 0; i < Knots1.size(); ++i)
                outfile << " " << Knots1[i];
            outfile << std::endl;
            for(std::size_t i = 0; i < Knots2.size(); ++i)
                outfile << " " << Knots2[i];
            outfile << "];" << std::endl;
        }

        outfile.close();
        std::cout << "Exported to " << fn << " completed!" << std::endl;

        std::cout << "Find cells in the T-splines topology mesh...";
        std::set<cell_t> cells;
        tmesh.FindCells(cells);
        std::cout << "OK!" << std::endl;

        std::cout << "Find cells in the extended T-splines topology mesh...";
        cells.clear();
        tmesh.FindCells(cells, true);
        std::cout << "OK!" << std::endl;
    }

    /// Export the T-mesh cells to mdpa format
    static void ExportMDPA(const TsMesh2D& tmesh, const std::string& fn, const std::vector<int>& division)
    {
        typedef TsMesh2D::cell_t cell_t;
        typedef TsMesh2D::anchor_container_t anchor_container_t;
        typedef TsMesh2D::cell_container_t cell_container_t;

        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for T-splines\n";
        outfile << "//(c) 2015 Hoang Giang Bui, Ruhr-University Bochum\n";

        std::time_t t = time(0);
        struct tm* now = std::localtime(&t);
        outfile << "//This file is created on " << now->tm_mday << "/" << now->tm_mon << "/" << (now->tm_year + 1900)
                << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl << std::endl;

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        for(anchor_container_t::const_iterator it = tmesh.Anchors().begin(); it != tmesh.Anchors().end(); ++it)
            outfile << (*it)->Id() << "\t" << (*it)->X() << "\t" << (*it)->Y() << "\t" << (*it)->Z() << std::endl;
        outfile << "End Nodes\n\n";

        // write elements
        outfile << "Begin Elements KinematicLinearGeo2dBezier\n";
        for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
        {
            const std::vector<std::size_t>& AnchorIds = (*it)->GetSupportedAnchors();
            outfile << (*it)->Id() << " 1";
            for(std::vector<std::size_t>::const_iterator it2 = AnchorIds.begin(); it2 != AnchorIds.end(); ++it2)
                outfile << " " << (*it2);
            outfile << std::endl;
        }
        outfile << "End Elements\n\n";

        // write weights
        outfile << "Begin ElementalData NURBS_WEIGHT\n";
        std::vector<double> AnchorWeights;
        for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
        {
            AnchorWeights = (*it)->GetAnchorWeights();
            outfile << (*it)->Id() << " [" << AnchorWeights.size() << "] (";
            for(std::size_t i = 0; i < AnchorWeights.size() - 1; ++i)
                outfile << AnchorWeights[i] << ",";
            outfile << AnchorWeights.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // read extraction operator from element in the form of triplet CSR
        // this requires that the Id of the cell must be unique
        std::vector<int> rowPtr;
        std::vector<int> colInd;
        std::vector<double> values;
        std::map<int, std::vector<int> > rowPtrMap;
        std::map<int, std::vector<int> > colIndMap;
        std::map<int, std::vector<double> > valuesMap;
        for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
        {
            rowPtr.clear();
            colInd.clear();
            values.clear();
            (*it)->GetExtractionOperator(rowPtr, colInd, values);
            rowPtrMap[(*it)->Id()] = rowPtr;
            colIndMap[(*it)->Id()] = colInd;
            valuesMap[(*it)->Id()] = values;
        }

        // write extraction operator
        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_ROWPTR\n";
        for(std::map<int, std::vector<int> >::iterator it = rowPtrMap.begin(); it != rowPtrMap.end(); ++it)
        {
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_COLIND\n";
        for(std::map<int, std::vector<int> >::iterator it = colIndMap.begin(); it != colIndMap.end(); ++it)
        {
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_VALUES\n";
        for(std::map<int, std::vector<double> >::iterator it = valuesMap.begin(); it != valuesMap.end(); ++it)
        {
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // write the degree
        outfile << "Begin ElementalData NURBS_DEGREE_1\n";
        for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
            outfile << (*it)->Id() << " " << tmesh.Order(0) << std::endl;
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData NURBS_DEGREE_2\n";
        for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
            outfile << (*it)->Id() << " " << tmesh.Order(1) << std::endl;
        outfile << "End ElementalData\n\n";

        if (division.size() != 0)
        {
            // write the division
            outfile << "Begin ElementalData NUM_DIVISION_1\n";
            for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
                outfile << (*it)->Id() << " " << division[0] << std::endl;
            outfile << "End ElementalData\n\n";

            outfile << "Begin ElementalData NUM_DIVISION_2\n";
            for(cell_container_t::const_iterator it = tmesh.Cells().begin(); it != tmesh.Cells().end(); ++it)
                outfile << (*it)->Id() << " " << division[1] << std::endl;
            outfile << "End ElementalData\n\n";
        }

        outfile.close();
        std::cout << "ExportMDPA completed" << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Testing
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TSplineUtils";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TSplineUtils& operator=(TSplineUtils const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    TSplineUtils(TSplineUtils const& rOther)
    {
    }

    ///@}

}; // Class TSplineUtils

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, TSplineUtils& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
        const TSplineUtils& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_TSPLINE_UTILS_H_INCLUDED  defined
