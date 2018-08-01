//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED

// System includes
#include <vector>
#include <tuple>

// External includes
#include <boost/any.hpp>
#include <boost/array.hpp>
#include <boost/enable_shared_from_this.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/array_1d.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/weighted_fespace.h"
#include "custom_utilities/control_grid_utility.h"
#include "isogeometric_application/isogeometric_application.h"

#define DEBUG_DESTROY

#define CONVERT_INDEX_IGA_TO_KRATOS(n) (n+1)
#define CONVERT_INDEX_KRATOS_TO_IGA(n) (n-1)

namespace Kratos
{


// Forward Declaration
template<int TDim> class MultiPatch;
template<int TDim> class PatchInterface;


/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
 */
template<int TDim>
class Patch : public boost::enable_shared_from_this<Patch<TDim> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;
    typedef Transformation<double> TransformationType;

    typedef GridFunction<TDim, double> DoubleGridFunctionType;
    typedef std::vector<typename DoubleGridFunctionType::Pointer> DoubleGridFunctionContainerType;

    typedef GridFunction<TDim, array_1d<double, 3> > Array1DGridFunctionType;
    typedef std::vector<typename Array1DGridFunctionType::Pointer> Array1DGridFunctionContainerType;

    typedef GridFunction<TDim, Vector> VectorGridFunctionType;
    typedef std::vector<typename VectorGridFunctionType::Pointer> VectorGridFunctionContainerType;

    typedef std::vector<typename Patch<TDim>::Pointer> NeighborPatchContainerType;

    typedef std::size_t vertex_t;
    typedef std::tuple<std::size_t, std::size_t, std::size_t, int> edge_t;
    //                  vertex1     vertex2     knot index   is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, int> face_t;
    //                  vertex1     vertex2     vertex3         vertex4     is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t> volume_t;

    typedef FESpace<TDim> FESpaceType;

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id), mpFESpace(NULL), mPrefix("Patch")
    {
    }

    /// Constructor with id and FESpace
    Patch(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace) : mId(Id), mpFESpace(pFESpace), mPrefix("Patch")
    {
        if (mpFESpace == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Invalid FESpace is provided", "")
    }

    /// Destructor
    virtual ~Patch()
    {
        #ifdef DEBUG_DESTROY
        std::cout << Type() << ", Id = " << Id()
                  << ", " << mpFESpace->Type()
                  << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Helper function to create new patch pointer
    static typename Patch<TDim>::Pointer Create(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace));
    }

    /// Get the working space dimension of the patch
    std::size_t WorkingSpaceDimension() const {return TDim;}

    /// Set the prefix for the patch
    void SetPrefix(const std::string& prefix) {mPrefix = prefix;}

    /// Get the prefix of the patch
    const std::string& Prefix() const {return mPrefix;}

    /// Set the Id of this patch
    void SetId(const std::size_t& Id) {mId = Id;}

    /// Get the name of the patch. The name is prefix + id
    std::string Name() const
    {
        std::stringstream ss;
        ss << mPrefix << "_" << mId;
        return ss.str();
    }

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Set the corresponding FESpace for the patch
    void SetFESpace(typename FESpace<TDim>::Pointer pFESpace) {mpFESpace = pFESpace;}

    /// Get the FESpace pointer
    typename FESpace<TDim>::Pointer pFESpace() {return mpFESpace;}

    /// Get the FESpace pointer
    typename FESpace<TDim>::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        assert(mpFESpace == NULL);
        return mpFESpace->TotalNumber();
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        assert(mpFESpace == NULL);
        if (i >= TDim) return 0;
        else return mpFESpace->Order(i);
    }

    /// Return true if this patch is a primary patch
    virtual bool IsPrimary() const
    {
        return true;
    }

    /// Enumerate the patch
    virtual void Enumerate()
    {
        std::size_t last = 0;
        last = this->pFESpace()->Enumerate(last);
    }

    /// Get the string representing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "Patch" << TDim << "D";
        return ss.str();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create the control point grid
    typename GridFunction<TDim, ControlPointType>::Pointer CreateControlPointGridFunction(typename ControlGrid<ControlPointType>::Pointer pControlPointGrid)
    {
        CheckSize(*pControlPointGrid, __FUNCTION__);
        pControlPointGrid->SetName("CONTROL_POINT");
        typename GridFunction<TDim, ControlPointType>::Pointer pNewGridFunc = GridFunction<TDim, ControlPointType>::Create(mpFESpace, pControlPointGrid);
        mpGridFunctions.push_back(pNewGridFunc);

        // create additional grid for control point coordinates, in order to compute the derivatives
        typedef typename ControlPointType::CoordinatesType CoordinatesType;
        ControlGrid<CoordinatesType>::Pointer pControlPointCoordinatesGrid = ControlGridUtility::CreateControlPointValueGrid<ControlPointType>(pControlPointGrid);
        pControlPointCoordinatesGrid->SetName("CONTROL_POINT_COORDINATES");
        typename FESpace<TDim>::Pointer pNewFESpace = WeightedFESpace<TDim>::Create(mpFESpace, this->GetControlWeights());
        typename GridFunction<TDim, CoordinatesType>::Pointer pNewCoordinatesGridFunc = GridFunction<TDim, CoordinatesType>::Create(pNewFESpace, pControlPointCoordinatesGrid);
        mpGridFunctions.push_back(pNewCoordinatesGridFunc);

        return pNewGridFunc;
    }

    /// Get the control point grid function
    GridFunction<TDim, ControlPointType>& ControlPointGridFunction() {return *(this->pGetGridFunction(CONTROL_POINT));}

    /// Get the control point grid function
    const GridFunction<TDim, ControlPointType>& ControlPointGridFunction() const {return *(this->pGetGridFunction(CONTROL_POINT));}

    /// Get the control point grid function pointer
    typename GridFunction<TDim, ControlPointType>::Pointer pControlPointGridFunction() {return this->pGetGridFunction(CONTROL_POINT);}

    /// Get the control point grid
    typename GridFunction<TDim, ControlPointType>::ConstPointer pControlPointGridFunction() const {return this->pGetGridFunction(CONTROL_POINT);}

    /// Get the control point weights vector
    std::vector<double> GetControlWeights() const
    {
        typename ControlGrid<ControlPointType>::ConstPointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        std::vector<double> Weights(pControlPointGrid->size());
        for (std::size_t i = 0; i < pControlPointGrid->size(); ++i)
            Weights[i] = (*pControlPointGrid)[i].W();
        return Weights;
    }

    /// Apply the homogeneous transformation to the patch by applying the homogeneous transformation to the control points grid. For DISPLACEMENT, access the grid function for DISPLACEMENT directly and transform it.
    void ApplyTransformation(const TransformationType& trans)
    {
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        ControlGridUtility::ApplyTransformation(*pControlPointGrid, trans);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function. This function will create the new FESpace based on the original FESpace of the control grid and the weights, and then assign to the new grid function.
    /// One must not use this function for the ControlPoint data type.
    template<typename TDataType>
    typename GridFunction<TDim, TDataType>::Pointer CreateGridFunction(typename ControlGrid<TDataType>::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename FESpace<TDim>::Pointer pNewFESpace = WeightedFESpace<TDim>::Create(mpFESpace, this->GetControlWeights());
        typename GridFunction<TDim, TDataType>::Pointer pNewGridFunc = GridFunction<TDim, TDataType>::Create(pNewFESpace, pControlGrid);
        mpGridFunctions.push_back(pNewGridFunc);
        return pNewGridFunc;
    }

    /// Create and add the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::Pointer CreateGridFunction(const TVariableType& rVariable,
            typename ControlGrid<typename TVariableType::Type>::Pointer pControlGrid)
    {
        pControlGrid->SetName(rVariable.Name());
        return this->CreateGridFunction<typename TVariableType::Type>(pControlGrid);
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::Pointer pGetGridFunction(const TVariableType& rVariable)
    {
        typedef typename GridFunction<TDim, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (std::size_t i = 0; i < mpGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(mpGridFunctions[i]);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                    return pGridFunc;
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        std::stringstream ss;
        ss << "The grid function with control grid " << rVariable.Name() << " does not exist in the database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, typename TVariableType::Type>::ConstPointer pGetGridFunction(const TVariableType& rVariable) const
    {
        typedef typename GridFunction<TDim, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (std::size_t i = 0; i < mpGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(mpGridFunctions[i]);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                    return pGridFunc;
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        std::stringstream ss;
        ss << "The grid function with control grid " << rVariable.Name() << " does not exist in the database";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }

    /// Filter out and get the underlying double grid functions
    DoubleGridFunctionContainerType DoubleGridFunctions() {return this->ExtractGridFunctions<DoubleGridFunctionContainerType>(mpGridFunctions);}
    DoubleGridFunctionContainerType DoubleGridFunctions() const {return this->ExtractGridFunctions<DoubleGridFunctionContainerType>(mpGridFunctions);}

    /// Filter out and get the underlying array_1d grid functions
    Array1DGridFunctionContainerType Array1DGridFunctions() {return this->ExtractGridFunctions<Array1DGridFunctionContainerType>(mpGridFunctions);}
    Array1DGridFunctionContainerType Array1DGridFunctions() const {return this->ExtractGridFunctions<Array1DGridFunctionContainerType>(mpGridFunctions);}

    /// Filter out and get the underlying Vector grid functions
    VectorGridFunctionContainerType VectorGridFunctions() {return this->ExtractGridFunctions<VectorGridFunctionContainerType>(mpGridFunctions);}
    VectorGridFunctionContainerType VectorGridFunctions() const {return this->ExtractGridFunctions<VectorGridFunctionContainerType>(mpGridFunctions);}

    /// Check if the grid function with name existed in the patch
    template<class TVariableType>
    bool HasGridFunction(const TVariableType& rVariable) const
    {
        std::vector<TVariableType*> var_list = this->ExtractVariables<TVariableType>();
        for (std::size_t i = 0; i < var_list.size(); ++i)
            if (*(var_list[i]) == rVariable)
                return true;
        return false;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Extract the Kratos variables from Grid functions. It is important that the Grid function has the same name and type as Kratos variable.
    template<class TVariableType>
    std::vector<TVariableType*> ExtractVariables() const
    {
        return this->ExtractVariables<TVariableType>(mpGridFunctions);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Validate the patch
    virtual bool Validate() const
    {
        if (Id() == 0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The patch must have an Id", "")
        }

        if (pControlPointGridFunction() != NULL)
            if (pControlPointGridFunction()->pControlGrid()->Size() != this->TotalNumber())
                KRATOS_THROW_ERROR(std::logic_error, "The control point grid is incompatible", "")

        DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The double variable grid is incompatible", (*it)->pControlGrid()->Name())
                return false;
            }
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The array_1d variable grid is incompatible", (*it)->pControlGrid()->Name())
                return false;
            }
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The vector variable grid is incompatible", (*it)->pControlGrid()->Name())
                return false;
            }
        }

        return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Construct the boundary patch based on side
    virtual typename Patch<TDim-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    {
        typename Patch<TDim-1>::Pointer pBPatch = typename Patch<TDim-1>::Pointer(new Patch<TDim-1>(-1));

        // construct the boundary FESpace
        typename FESpace<TDim-1>::Pointer pBFESpace = this->pFESpace()->ConstructBoundaryFESpace(side);
        pBPatch->SetFESpace(pBFESpace);

        // extract the local id associated with the boundary indices
        std::vector<std::size_t> local_ids = this->pFESpace()->LocalId(pBFESpace->FunctionIndices());

        // transfer the control values
//        typename ControlGrid<ControlPointType>::Pointer pBoundaryControlPointGrid = ControlGridUtility::ExtractSubGrid<ControlPointType>(this->pControlPointGridFunction()->pControlGrid(), local_ids);
        typename ControlGrid<ControlPointType>::Pointer pBoundaryControlPointGrid = ControlGridUtility::ExtractSubGrid<TDim, ControlPointType>(this->pControlPointGridFunction()->pControlGrid(), *(this->pFESpace()), *pBFESpace);
        pBPatch->CreateControlPointGridFunction(pBoundaryControlPointGrid);

        // transfer other values
        DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
//            typename ControlGrid<double>::Pointer pBoundaryDoubleControlGrid = ControlGridUtility::ExtractSubGrid<double>((*it)->pControlGrid(), local_ids);
            typename ControlGrid<double>::Pointer pBoundaryDoubleControlGrid = ControlGridUtility::ExtractSubGrid<TDim, double>((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<double>(pBoundaryDoubleControlGrid);
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
//            typename ControlGrid<array_1d<double, 3> >::Pointer pBoundaryArray1DControlGrid = ControlGridUtility::ExtractSubGrid<array_1d<double, 3> >((*it)->pControlGrid(), local_ids);
            typename ControlGrid<array_1d<double, 3> >::Pointer pBoundaryArray1DControlGrid = ControlGridUtility::ExtractSubGrid<TDim, array_1d<double, 3> >((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<array_1d<double, 3> >(pBoundaryArray1DControlGrid);
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
//            typename ControlGrid<Vector>::Pointer pBoundaryVectorControlGrid = ControlGridUtility::ExtractSubGrid<Vector>((*it)->pControlGrid(), local_ids);
            typename ControlGrid<Vector>::Pointer pBoundaryVectorControlGrid = ControlGridUtility::ExtractSubGrid<TDim, Vector>((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<Vector>(pBoundaryVectorControlGrid);
        }

        return pBPatch;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Search for the neighbor
    typename Patch<TDim>::Pointer pNeighbor(const BoundarySide& side)
    {
        for (std::size_t i = 0; i < this->NumberOfInterfaces(); ++i)
        {
            if (this->pInterface(i)->Side1() == side)
                return this->pInterface(i)->pPatch2();
        }
        return NULL;
    }

    typename Patch<TDim>::ConstPointer pNeighbor(const BoundarySide& side) const
    {
        for (std::size_t i = 0; i < this->NumberOfInterfaces(); ++i)
        {
            if (this->pInterface(i)->Side1() == side)
                return this->pInterface(i)->pPatch2();
        }
        return NULL;
    }

    /// Seach for the boundary side of the neighor patch, if it exists
    int FindBoundarySide(typename Patch<TDim>::ConstPointer pPatch) const
    {
        for (std::size_t i = 0; i < this->NumberOfInterfaces(); ++i)
        {
            if (this->pNeighbor(i)->pPatch2() == pPatch)
                return this->pNeighbor(i)->Side1();
        }
        return -1;
    }

    /// Add an interface to the patch
    void AddInterface(typename PatchInterface<TDim>::Pointer pInterface)
    {
        mpInterfaces.push_back(pInterface);
    }

    /// Get the number of interfaces
    std::size_t NumberOfInterfaces() const
    {
        return mpInterfaces.size();
    }

    /// Get the interface
    typename PatchInterface<TDim>::Pointer pInterface(const std::size_t& i)
    {
        return mpInterfaces[i];
    }

    typename PatchInterface<TDim>::ConstPointer pInterface(const std::size_t& i) const
    {
        return mpInterfaces[i];
    }

    /// Get/Set the parent multipatch
    void pSetParentMultiPatch(typename MultiPatch<TDim>::Pointer pPatch) {mpParentMultiPatch = pPatch;}
    MultiPatch<TDim>& ParentMultiPatch() {return *pParentMultiPatch();}
    const MultiPatch<TDim>& ParentMultiPatch() const {return *pParentMultiPatch();}
    typename MultiPatch<TDim>::Pointer pParentMultiPatch() {return mpParentMultiPatch.lock();}
    const typename MultiPatch<TDim>::Pointer pParentMultiPatch() const {return mpParentMultiPatch.lock();}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Generate topology data to visualize with GLVis
    void GenerateTopolgyData(std::size_t& starting_vertex_id,
            std::vector<vertex_t>& vertices,
            std::vector<edge_t>& edges,
            std::vector<face_t>& faces,
            std::vector<volume_t>& volumes,
            std::size_t& starting_knotv_id,
            std::vector<std::size_t>& knotv ) const
    {
        if (TDim == 1)
        {
            vertices.resize(2);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;

            knotv.resize(0);
            knotv[0] = starting_knotv_id++;

            edges.resize(1);
            edges[0] = std::make_tuple(vertices[0], vertices[1], knotv[0], 0);

            faces.resize(0);
            volumes.resize(0);
        }
        else if (TDim == 2)
        {
            /// Reference for edge mapping: Fig.2, Burstedde et al, SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
            /// Mapping for edges: table 2
            vertices.resize(4);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;
            vertices[2] = starting_vertex_id++;
            vertices[3] = starting_vertex_id++;

            knotv.resize(2);
            knotv[0] = starting_knotv_id++;
            knotv[1] = starting_knotv_id++;

            edges.resize(4);
            edges[0] = std::make_tuple(vertices[0], vertices[2], knotv[1], 1);
            edges[1] = std::make_tuple(vertices[1], vertices[3], knotv[1], 1);
            edges[2] = std::make_tuple(vertices[0], vertices[1], knotv[0], 1);
            edges[3] = std::make_tuple(vertices[2], vertices[3], knotv[0], 1);

            faces.resize(1);
            faces[0] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], 0);

            volumes.resize(0);
        }
        else if (TDim == 3)
        {
            /// Reference for edge mapping: Fig.2, Burstedde et al, SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
            /// Mapping for faces: table 2
            vertices.resize(8);
            vertices[0] = starting_vertex_id++;
            vertices[1] = starting_vertex_id++;
            vertices[2] = starting_vertex_id++;
            vertices[3] = starting_vertex_id++;
            vertices[4] = starting_vertex_id++;
            vertices[5] = starting_vertex_id++;
            vertices[6] = starting_vertex_id++;
            vertices[7] = starting_vertex_id++;

            knotv.resize(3);
            knotv[0] = starting_knotv_id++;
            knotv[1] = starting_knotv_id++;
            knotv[2] = starting_knotv_id++;

            edges.resize(12);
            edges[0] = std::make_tuple(vertices[0], vertices[1], knotv[0], 1);
            edges[1] = std::make_tuple(vertices[2], vertices[3], knotv[0], 1);
            edges[2] = std::make_tuple(vertices[4], vertices[5], knotv[0], 1);
            edges[3] = std::make_tuple(vertices[6], vertices[7], knotv[0], 1);
            edges[4] = std::make_tuple(vertices[0], vertices[2], knotv[1], 1);
            edges[5] = std::make_tuple(vertices[1], vertices[3], knotv[1], 1);
            edges[6] = std::make_tuple(vertices[4], vertices[6], knotv[1], 1);
            edges[7] = std::make_tuple(vertices[5], vertices[7], knotv[1], 1);
            edges[8] = std::make_tuple(vertices[0], vertices[4], knotv[2], 1);
            edges[9] = std::make_tuple(vertices[1], vertices[5], knotv[2], 1);
            edges[10] = std::make_tuple(vertices[2], vertices[6], knotv[2], 1);
            edges[11] = std::make_tuple(vertices[3], vertices[7], knotv[2], 1);

            faces.resize(6);
            faces[0] = std::make_tuple(vertices[0], vertices[2], vertices[4], vertices[6], 1);
            faces[1] = std::make_tuple(vertices[1], vertices[3], vertices[5], vertices[7], 1);
            faces[2] = std::make_tuple(vertices[0], vertices[1], vertices[4], vertices[5], 1);
            faces[3] = std::make_tuple(vertices[2], vertices[3], vertices[6], vertices[7], 1);
            faces[4] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], 1);
            faces[5] = std::make_tuple(vertices[4], vertices[5], vertices[6], vertices[7], 1);

            volumes.resize(1);
            volumes[0] = std::make_tuple(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5], vertices[6], vertices[7]);
        }
        else
        {
            std::stringstream ss;
            ss << __FUNCTION__ << " is not implemented for " << TDim;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two patches in terms of its parametric information. The grid function data, including control points, shall not be checked.
    virtual bool IsCompatible(const Patch<TDim>& rOtherPatch) const
    {
        return *(this->pFESpace()) == *(rOtherPatch.pFESpace());
    }

    /// Compare between two patches in terms of parametric information and control points.
    bool IsEquivalent(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsCompatible(rOtherPatch))
            return false;

        // TODO compare the control points

        return true;
    }

    /// Compare between two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const Patch<TDim>& rOtherPatch) const
    {
        if (!this->IsEquivalent(rOtherPatch))
            return false;

        // TODO compare the grid function values

        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<TDim>& rOther)
    {
        return (Id() == rOther.Id()) && this->IsSame(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Id = " << Id() << ", Addr = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        if (pFESpace() != NULL)
            rOStream << *pFESpace() << std::endl;

        if (pControlPointGridFunction() != NULL)
            rOStream << *(pControlPointGridFunction()->pControlGrid()) << std::endl;

        DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();

        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            rOStream << *((*it)->pControlGrid()) << std::endl;
        }

        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            rOStream << *((*it)->pControlGrid()) << std::endl;
        }

        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            rOStream << *((*it)->pControlGrid()) << std::endl;
        }

        rOStream << "Interfaces (" << this->NumberOfInterfaces() << "):" << std::endl;
        for (std::size_t i = 0; i < this->NumberOfInterfaces(); ++i)
        {
            rOStream << "  ";
            this->pInterface(i)->PrintInfo(rOStream);
            rOStream << std::endl;
        }
    }

private:

    std::size_t mId;
    std::string mPrefix;

    // FESpace contains the shape function information and various information with regards to the functional space.
    // Because the control point grid is in homogeneous coordinates, the FESpace shall be an unweighted spaces
    typename FESpace<TDim>::Pointer mpFESpace;

    // container to contain all the grid functions
    std::vector<boost::any> mpGridFunctions; // using boost::any to store pointer to grid function

    /**
     * interface data
     */
    std::vector<typename PatchInterface<TDim>::Pointer> mpInterfaces;

    /**
     * pointer to parent multipatch
     */
    typename MultiPatch<TDim>::WeakPointer mpParentMultiPatch;

    /// Empty Constructor for serializer
    Patch() : mId(0), mpFESpace(NULL) {}

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    /// Auxiliary
    template<class TGridFunctionType>
    void CheckSize(const TGridFunctionType& rGrid, const std::string& source) const
    {
        if (rGrid.Size() != this->TotalNumber())
        {
            std::stringstream ss;
            ss << "The size of grid function (" << rGrid.size() << ") is not compatible with the current number of control values (" << this->TotalNumber()
               << ") of patch " << Id() << ". Error at " << source;
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /// Helper to extract the grid functions out from boost::any
    template<class TContainerType>
    TContainerType ExtractGridFunctions(const std::vector<boost::any>& pGridFunctions) const
    {
        TContainerType GridFuncs;

        typedef typename TContainerType::value_type GridFunctionPointerType;

        for (std::size_t i = 0; i < pGridFunctions.size(); ++i)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(pGridFunctions[i]);
                GridFuncs.push_back(pGridFunc);
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }

        return GridFuncs;
    }

    /// Helper to extract the variables out from Grid functions. It is important that the variable is already registered to the Kratos kernel.
    template<class TVariableType>
    std::vector<TVariableType*> ExtractVariables(const std::vector<boost::any>& pGridFunctions) const
    {
        typedef GridFunction<TDim, typename TVariableType::Type> GridFunctionType;
        typedef std::vector<typename GridFunctionType::Pointer> GridFunctionContainerType;
        GridFunctionContainerType GridFuncs = this->ExtractGridFunctions<GridFunctionContainerType>(pGridFunctions);

        std::vector<TVariableType*> var_list;
        for (std::size_t i = 0; i < GridFuncs.size(); ++i)
        {
            const std::string& var_name = GridFuncs[i]->pControlGrid()->Name();

            if (KratosComponents<VariableData>::Has(var_name))
            {
                var_list.push_back(dynamic_cast<TVariableType*>(&KratosComponents<VariableData>::Get(var_name)));
            }
        }

        return var_list;
    }
};


/**
 * Template specific instantiation for null-D patch to terminate the compilation.
 * In fact, null-D patch is a vertex
 */
template<>
class Patch<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    // Type definitions
    typedef ControlPoint<double> ControlPointType;

    /// Default constructor
    Patch() : mId(0) {}

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<0>::Pointer pFESpace) {}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "Patch0D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Set the control point grid
    typename GridFunction<0, ControlPointType>::Pointer CreateControlPointGridFunction(typename ControlGrid<ControlPointType>::Pointer pControlPointGrid)
    {
        return NULL;
    }

    /// Create and add the grid function
    template<typename TDataType>
    typename GridFunction<0, TDataType>::Pointer CreateGridFunction(typename ControlGrid<TDataType>::Pointer pControlGrid)
    {
        return NULL;
    }

    /// Compare between two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const Patch<0>& rOtherPatch) const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<0>& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const Patch<0>& rPatch1, const BoundarySide& side1,
            const Patch<0>& rPatch2, const BoundarySide& side2)
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/**
 * Template specific instantiation for -1-D patch to terminate the compilation.
 */
template<>
class Patch<-1>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Default constructor
    Patch() : mId(0) {}

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<-1>::Pointer pFESpace) {}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "Patch<-1>D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<-1>& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const Patch<-1>& rPatch1, const BoundarySide& side1,
            const Patch<-1>& rPatch2, const BoundarySide& side2)
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<-1>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/**
 * Template specific instantiation for -2-D patch to terminate the compilation.
 */
template<>
class Patch<-2>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);

    /// Default constructor
    Patch() : mId(0) {}

    /// Constructor with id
    Patch(const std::size_t& Id) : mId(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpace<-2>::Pointer pFESpace) {}

    /// Get the Id of this patch
    const std::size_t& Id() const {return mId;}

    /// Get the number of basis functions defined over the patch
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "Patch<-2>D";
    }

    /// Validate the patch
    virtual bool Validate() const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const Patch<-2>& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const Patch<-2>& rPatch1, const BoundarySide& side1,
            const Patch<-2>& rPatch2, const BoundarySide& side2)
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-2>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<-2>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mId;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const Patch<TDim>& rThis)
{
    rOStream << "-------------Begin PatchInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End PatchInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED defined

