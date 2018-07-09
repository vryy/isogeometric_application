//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 May 2015 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_HN_MESH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_HN_MESH_H_INCLUDED

// System includes
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

// External includes
#include <omp.h>
#include "boost/progress.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/nurbs/knot.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/cell_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/nurbs/domain_manager_3d.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/hbsplines/deprecated_basis_function_manager.h"
#include "custom_utilities/hbsplines/deprecated_hb_basis_function.h"

namespace Kratos
{

/**
Hierarchical B-Splines mesh and refinement
*/
template<int TDim>
class DeprecatedHBMesh : public Patch<TDim>
{
private:
    template<class data_type>
    class unordered_pair
    {
        public:
            unordered_pair(data_type i, data_type j)
            {
                if(i < j)
                {
                    mi = i;
                    mj = j;
                }
                else
                {
                    mi = j;
                    mj = i;
                }
            }

            ~unordered_pair() {}

            bool operator<(const unordered_pair& rOther) const
            {
                if(mi == rOther.mi)
                    return mj < rOther.mj;
                else
                    return mi < rOther.mi;
            }

            data_type first() const {return mi;}
            data_type second() const {return mj;}

        private:
            data_type mi, mj;
    };

    template<class data_type>
    class unordered_tuple
    {
        public:
            unordered_tuple(data_type i, data_type j, data_type k, data_type l)
            {
                // TODO
            }

            ~unordered_tuple() {}

            bool operator<(const unordered_tuple& rOther) const
            {
                // TODO
            }

        private:
            data_type mi, mj, mk, ml;
    };

public:
    /// Const definition
    static const int NO_READ          = 0;
    static const int READ_PATCH       = 1;
    static const int READ_ORDER       = 2;
    static const int READ_NUMBER      = 3;
    static const int READ_KNOTS       = 4;
    static const int READ_COORDINATES = 5;
    static const int READ_WEIGHTS     = 6;

//    static const int ECHO_KNOT_VECTORS  = 0b0000000000000001;
//    static const int ECHO_CELLS         = 0b0000000000000010;
//    static const int ECHO_BASIS_FUNCS   = 0b0000000000000100;
//    static const int ECHO_REFIMENT      = 0b0000000000001000;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DeprecatedHBMesh);

    /// Type definition
    typedef Patch<TDim> BaseType;

    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    typedef HBCell<DeprecatedHBBasisFunction> CellType;
    typedef typename CellType::Pointer cell_t;
    typedef CellManager<CellType> cell_container_t;

    typedef DeprecatedHBBasisFunction::Pointer bf_t;
    typedef DeprecatedBasisFunctionManager<DeprecatedHBBasisFunction> bf_container_t;

    typedef DomainManager::Pointer domain_t;
    typedef std::map<std::size_t, domain_t> domain_container_t;

    /**************************************************************************
                           CONSTRUCTOR
    **************************************************************************/

    /// Default constructor
    DeprecatedHBMesh(const std::size_t& Id, const std::string& Name);

    /// Destructor
    ~DeprecatedHBMesh() {}

    /**************************************************************************
                           FUNDAMENTAL SUBROUTINES
    **************************************************************************/

    /// Get the name of this hierarchical mesh
    const std::string& Name() const {return mName;}

    /// Set the echo level for this mesh
    void SetEchoLevel(const int& Level) {mEchoLevel = Level;}

    /// Get the echo level
    const int& GetEchoLevel() const {return mEchoLevel;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const {return mBasisFuncs.size();}

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        if (i == 0) return mOrder1;
        else if (i == 1) return mOrder2;
        else if (i == 2) return mOrder3;
        else return 0;
    }

    /// Set the maximum level allowed in the hierarchical mesh
    void SetMaxLevels(const std::size_t& MaxLevels) {mMaxLevels = MaxLevels;}

    /**************************************************************************
                            INPUT SUBROUTINES
    **************************************************************************/

    /// Read and construct the first level from file
    void ReadMesh(const std::string& fn);

    /**************************************************************************
                            REFINEMENT SUBROUTINES
    **************************************************************************/

    /// Refine a basis function specified by the Id
    void Refine(const std::size_t& Id);

    /// Refine a list of bfs provided by python list. The list will be automatically sorted by the level.
    void RefineNodes(boost::python::list& pyList);

    /// Refine on a cuboid domain
    void RefineWindow(const double& Xi_min, const double& Xi_max,
            const double& Eta_min, const double& Eta_max,
            const double& Zeta_min, const double& Zeta_max);

    /// Perform additional refinement to ensure linear independence
    /// In this algorithm, every bf in each level will be checked. If the support domain of a bf contained in the domain_manager of that level, this bf will be refined. According to the paper of Vuong et al, this will produce a linear independent bases.
    void LinearDependencyRefine(const std::size_t& refine_cycle);

    /// Finalize the hierarchical mesh construction and compute the Bezier decomposition at each cell
    void BuildMesh();

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "Hierarchical B-Splines Patch" << TDim << "D";
        return ss.str();
    }

    /// Print information of this hierarchical mesh
    void PrintInfo(std::ostream& rOStream) const;

    /// Print data of this hierarchical mesh
    void PrintData(std::ostream& rOStream) const;

    /// Print knot vectors to standard output
    void PrintKnotVectors() const;

    /// Print cells to standard output; level = 0 will print cells of all the level
    void PrintCells(int level) const;

    /// Print basis functions to standard output; level = 0 will print bfs of all the level
    void PrintBasisFuncs(int level) const;

    /// Print local knot vector of a bf
    void PrintLocalKnotVectors(unsigned int Id) const;

    /// Print the refinement history
    void PrintRefinementHistory() const;

    /// Build the nested spaces K as in Definition 1, paper Vuong et al
    void BuildNestedSpace(std::size_t Level, std::map<std::size_t, std::set<std::size_t> >& rK);

    /// Check if the nested space condition is correct on all level
    void CheckNestedSpace();

    /**************************************************************************
                            AUXILIARY SUBROUTINES
    **************************************************************************/

    /// Validate the patch
    virtual bool Validate() const;

    /// Build the hierarchical boundary mesh
    void BuildBoundaryMesh(DeprecatedHBMesh& rMesh, std::string boundary_mesh_type) const;

    /**************************************************************************
                            MATLAB INTERFACE
    **************************************************************************/

    /// Export the cell topology to Matlab for visualization. Cell topology is the cell organization in the parameter space.
    void ExportCellTopology(std::string fn, bool cell_numbering) const;

    /// Export the cell geoology to Matlab for visualization. Cell geology is the cell organization in the physical space.
    void ExportCellGeology(std::string fn);

    /// Export the support domain of all levels to matlab for visualization
    void ExportSupportDomain(std::string fn, double z_distance);

    /// Export the mesh information (basis functions, Bezier elements) to Matlab for plotting and debugging
    void ExportMatlab(std::string fn) const;

    /**************************************************************************
                            KRATOS INTERFACE
    **************************************************************************/

    /// Export this hierarchical mesh to MDPA
    void ExportMDPA(std::string fn) const;
    void ExportMDPA2(std::string fn) const;

    /// Export this hierarchical mesh to post MDPA file for post-processing
    void ExportPostMDPA(std::string fn, int Division1, int Division2, int Division3);

    /// Export the cell geology as the post MDPA file
    void ExportCellGeologyAsPostMDPA(std::string fn);

    /**************************************************************************
                            DEBUG INTERFACE
    **************************************************************************/

private:

    std::string mName;
    unsigned int mEchoLevel;

    unsigned int mOrder1;
    unsigned int mOrder2;
    unsigned int mOrder3;

    knot_container_t mKnots1; // container for knot vector in u-direction
    knot_container_t mKnots2; // container for knot vector in v-direction
    knot_container_t mKnots3; // container for knot vector in v-direction

    std::size_t mLastLevel;
    std::size_t mMaxLevels;

    typename cell_container_t::Pointer mpCellManager;

    bf_container_t mBasisFuncs;

    domain_container_t mSupportDomains; // this domain manager manages the support of all bfs in each level

    std::vector<std::size_t> mRefinementHistory;

    /// Get the domain manager for support domain for level k
    domain_t GetSupportDomain(std::size_t Level)
    {
        domain_container_t::iterator it = mSupportDomains.find(Level);
        if(it != mSupportDomains.end())
            return it->second;
        else
        {
            domain_t p_domain;
            if(TDim == 2)
                p_domain = domain_t(new DomainManager2D(Level));
            else if(TDim == 3)
                p_domain = domain_t(new DomainManager3D(Level));
            mSupportDomains[Level] = p_domain;
            return p_domain;
        }
    }

    /// Generate the cell geological map for Matlab rendering and post processing
    void GenerateCellGeology(std::vector<unsigned int>& point_list, // list of node Id (of the active basis function) in the hierarchical mesh
                             std::map<unsigned int, double>& X_list, // global x-coordinate of the node in point_list
                             std::map<unsigned int, double>& Y_list, // global y-coordinate of the node in point_list
                             std::map<unsigned int, double>& Z_list, // global z-coordinate of the node in point_list
                             std::map<unsigned int, double>& xi_list, // local xi-coordinate of the node in point_list
                             std::map<unsigned int, double>& eta_list, // local eta-coordinate of the node in point_list
                             std::map<unsigned int, double>& zeta_list, // local eta-coordinate of the node in point_list
                             std::map<unsigned int, unsigned int>& cell_list, // cell index of the node in point_list
                             std::map<unsigned int, std::vector<std::vector<unsigned int> > >& Connectivities // element connectivities
                             );

};

/**
 * Template specific instantiation for null-D hierarchical B-Splines patch to terminate the compilation
 */
template<>
class DeprecatedHBMesh<0> : public Patch<0>
{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(DeprecatedHBMesh);

    /// Default constructor
    DeprecatedHBMesh(const std::size_t& Id, const std::string& Name) : Patch<0>(Id), mName(Name) {}

    /// Destructor
    ~DeprecatedHBMesh() {}

    /// Get the name of this hierarchical mesh
    const std::string& Name() const {return mName;}

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "Hierarchical B-Splines Patch0D";
    }

private:
    std::string mName;
};

/// output stream function
template<int TDim>
inline std::ostream& operator<<(std::ostream& rOStream, const DeprecatedHBMesh<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

// #undef DEBUG_LEVEL1
// #undef DEBUG_LEVEL2
// #undef ENABLE_PROFILING

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HN_MESH_2D_H_INCLUDED

