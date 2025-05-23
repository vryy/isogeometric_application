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
#include <list>
#include <tuple>
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include <memory> // for std::enable_shared_from_this
#endif

// External includes
#include <boost/any.hpp>
#include <boost/array.hpp>
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include <boost/enable_shared_from_this.hpp>
#endif

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/kratos_flags.h"
#include "containers/array_1d.h"
#include "utilities/indexed_object.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/weighted_fespace.h"
#include "custom_utilities/control_grid_utility.h"
#include "isogeometric_application_variables.h"

#define CONVERT_INDEX_IGA_TO_KRATOS(n) (n+1)
#define CONVERT_INDEX_KRATOS_TO_IGA(n) (n-1)

namespace Kratos
{

// Forward Declaration
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType> class MultiPatch;
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType> class PatchInterface;

/**
This class represents an isogeometric patch in parametric coordinates. An isogeometric patch can be a NURBS patch, a hierarchical BSplines patch, or a T-Splines patch.
 */
template<int TDim, typename TLocalCoordinateType = double, typename TCoordinateType = double, typename TDataType = double>
class Patch : public IndexedObject, public Flags
#ifdef SD_APP_FORWARD_COMPATIBILITY
    , public std::enable_shared_from_this<Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#else
    , public boost::enable_shared_from_this<Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> >
#endif
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const Patch> ConstPointer;
#endif

    /// Type definition
    typedef Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;
    typedef Patch<TDim-1, TLocalCoordinateType, TCoordinateType, TDataType> BoundaryPatchType;
    typedef MultiPatch<TDim, TLocalCoordinateType, TCoordinateType, TDataType> MultiPatchType;
    typedef PatchInterface<TDim, TLocalCoordinateType, TCoordinateType, TDataType> PatchInterfaceType;

    typedef TLocalCoordinateType LocalCoordinateType;
    typedef TCoordinateType CoordinateType;
    typedef TDataType DataType;
    typedef typename DataTypeToValueType<TCoordinateType>::value_type CoordinateValueType;
    typedef CoordinateValueType WeightType;
    typedef typename MatrixVectorTypeSelector<DataType>::VectorType VectorType;

    typedef ControlPoint<CoordinateType, WeightType> ControlPointType;
    typedef Transformation<CoordinateType> TransformationType;

    typedef std::map<std::string, boost::any> GridFunctionContainerType;

    typedef GridFunction<TDim, TLocalCoordinateType, TDataType> DoubleGridFunctionType;
    typedef std::vector<typename DoubleGridFunctionType::Pointer> DoubleGridFunctionContainerType;

    typedef GridFunction<TDim, TLocalCoordinateType, array_1d<TDataType, 3> > Array1DGridFunctionType;
    typedef std::vector<typename Array1DGridFunctionType::Pointer> Array1DGridFunctionContainerType;

    typedef GridFunction<TDim, TLocalCoordinateType, VectorType> VectorGridFunctionType;
    typedef std::vector<typename VectorGridFunctionType::Pointer> VectorGridFunctionContainerType;

    typedef GridFunction<TDim, TLocalCoordinateType, ControlPointType> ControlPointGridFunctionType;

    typedef std::vector<typename PatchType::Pointer> NeighborPatchContainerType;

//    typedef std::vector<typename PatchInterfaceType::Pointer> InterfaceContainerType;
    typedef std::list<typename PatchInterfaceType::Pointer> InterfaceContainerType;
    typedef typename InterfaceContainerType::iterator interface_iterator;
    typedef typename InterfaceContainerType::const_iterator interface_const_iterator;

    typedef std::size_t vertex_t;
    typedef std::tuple<std::size_t, std::size_t, std::size_t, int> edge_t;
    //                  vertex1     vertex2     knot index   is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, int> face_t;
    //                  vertex1     vertex2     vertex3         vertex4     is_boundary
    typedef std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t> volume_t;

    typedef FESpace<TDim, LocalCoordinateType> FESpaceType;
    typedef typename FESpaceType::BoundaryFESpaceType BoundaryFESpaceType;
    typedef WeightedFESpace<TDim, LocalCoordinateType, WeightType> WeightedFESpaceType;

    /// Constants
    static constexpr int Dim = TDim;
    static constexpr CoordinateValueType DISTANCE_TOLERANCE = 1e-13;

    /// Constructor with id
    Patch(std::size_t Id)
        : IndexedObject(Id), mPrefix("Patch"), mLayerIndex(Id), mpFESpace(NULL)
        , mLocalSearchMaxIters(30), mLocalSearchTolerance(1e-8)
    {
        this->Set(ACTIVE, true);
    }

    /// Constructor with id and FESpace
    Patch(std::size_t Id, typename FESpaceType::Pointer pFESpace)
        : IndexedObject(Id), mPrefix("Patch"), mLayerIndex(Id), mpFESpace(pFESpace)
        , mLocalSearchMaxIters(30), mLocalSearchTolerance(1e-8)
    {
        this->Set(ACTIVE, true);
        if (mpFESpace == nullptr)
            KRATOS_ERROR << "Invalid FESpace is provided";
        }

    /// Destructor
    ~Patch() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << Type() << ", Id = " << Id()
                  << ", " << mpFESpace->Type()
                  << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper function to create new patch pointer
    static typename PatchType::Pointer Create(std::size_t Id, typename FESpaceType::Pointer pFESpace)
    {
        return typename PatchType::Pointer(new PatchType(Id, pFESpace));
    }

    /// Get the working space dimension of the patch
    constexpr std::size_t WorkingSpaceDimension() {return TDim;}

    /// Set the prefix for the patch
    void SetPrefix(const std::string& prefix) {mPrefix = prefix;}

    /// Set the layer index
    void SetLayerIndex(int Index) {mLayerIndex = Index;}

    /// Get the layer index
    int LayerIndex() {return mLayerIndex;}

    /// Get the prefix of the patch
    const std::string& Prefix() const {return mPrefix;}

    /// Get the name of the patch. The name is prefix + id
    std::string Name() const
    {
        std::stringstream ss;
        ss << mPrefix << "_" << Id();
        return ss.str();
    }

    /// Set the local search max iteration
    void SetLocalSearchMaxIters(unsigned int max_iters)
    {
        mLocalSearchMaxIters = max_iters;
    }

    /// Set the local search tolerance
    void SetLocalSearchTolerance(CoordinateType tolerance)
    {
        mLocalSearchTolerance = tolerance;
    }

    /// Set the corresponding FESpace for the patch
    void SetFESpace(typename FESpaceType::Pointer pFESpace) {mpFESpace = pFESpace;}

    /// Get the FESpace pointer
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Get the FESpace pointer
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the number of basis functions defined over the patch
    virtual std::size_t TotalNumber() const
    {
        assert(mpFESpace != NULL);
        return mpFESpace->TotalNumber();
    }

    /// Get the order of the patch in specific direction
    virtual std::size_t Order(std::size_t i) const
    {
        assert(mpFESpace != NULL);
        if (i >= TDim) { return 0; }
        else { return mpFESpace->Order(i); }
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
        std::stringstream ss;
        ss << "Patch<" << TDim
           << ", " << DataTypeToString<LocalCoordinateType>::Get()
           << ", " << DataTypeToString<CoordinateType>::Get()
           << ", " << DataTypeToString<DataType>::Get()
           << ">";
        return ss.str();
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
    typename ControlPointGridFunctionType::Pointer CreateControlPointGridFunction(typename ControlGrid<ControlPointType>::Pointer pControlPointGrid)
    {
        CheckSize(*pControlPointGrid, __FUNCTION__);
        pControlPointGrid->SetName("CONTROL_POINT");
        typename ControlPointGridFunctionType::Pointer pNewGridFunc = ControlPointGridFunctionType::Create(mpFESpace, pControlPointGrid);
        mpGridFunctions["CONTROL_POINT"] = pNewGridFunc;

        // create additional grid for control point coordinates, in order to compute the derivatives
        typedef typename ControlPointType::CoordinatesType CoordinatesType;
        typename ControlGrid<CoordinatesType>::Pointer pControlPointCoordinatesGrid = ControlGridUtility::CreateControlPointValueGrid<ControlPointType>(pControlPointGrid);
        pControlPointCoordinatesGrid->SetName("CONTROL_POINT_COORDINATES");
        typename FESpaceType::Pointer pNewFESpace = WeightedFESpaceType::Create(mpFESpace, this->GetControlWeights());
        typename GridFunction<TDim, TLocalCoordinateType, CoordinatesType>::Pointer pNewCoordinatesGridFunc = GridFunction<TDim, TLocalCoordinateType, CoordinatesType>::Create(pNewFESpace, pControlPointCoordinatesGrid);
        mpGridFunctions["CONTROL_POINT_COORDINATES"] = pNewCoordinatesGridFunc;

        return pNewGridFunc;
    }

    /// Get the control point grid function
    ControlPointGridFunctionType& ControlPointGridFunction()
    {
        return *(this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT)));
    }

    /// Get the control point grid function
    const ControlPointGridFunctionType& ControlPointGridFunction() const
    {
        return *(this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT)));
    }

    /// Get the control point grid function pointer
    typename ControlPointGridFunctionType::Pointer pControlPointGridFunction()
    {
        return this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT));
    }

    /// Get the control point grid
    typename ControlPointGridFunctionType::ConstPointer pControlPointGridFunction() const
    {
        return this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT));
    }

    /// Get the control point weights vector
    std::vector<WeightType> GetControlWeights() const
    {
        typename ControlGrid<ControlPointType>::ConstPointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        std::vector<WeightType> Weights(pControlPointGrid->size());
        for (std::size_t i = 0; i < pControlPointGrid->size(); ++i)
        {
            Weights[i] = (*pControlPointGrid)[i].W();
        }
        return Weights;
    }

    /// Apply the homogeneous transformation to the patch by applying the homogeneous transformation to the control points grid. For DISPLACEMENT, access the grid function for DISPLACEMENT directly and transform it.
    void ApplyTransformation(const TransformationType& trans)
    {
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        ControlGridUtility::ApplyTransformation(*pControlPointGrid, trans);

        // transform the control point coordinates grid
        typedef typename ControlPointType::CoordinatesType CoordinatesType;
        typename ControlGrid<CoordinatesType>::Pointer pControlPointCoordinatesGrid = ControlGridUtility::CreateControlPointValueGrid<ControlPointType>(pControlPointGrid);
        pControlPointCoordinatesGrid->SetName("CONTROL_POINT_COORDINATES");
        typename FESpaceType::Pointer pNewFESpace = WeightedFESpaceType::Create(mpFESpace, this->GetControlWeights());
        typename GridFunction<TDim, TLocalCoordinateType, CoordinatesType>::Pointer pNewCoordinatesGridFunc = GridFunction<TDim, TLocalCoordinateType, CoordinatesType>::Create(pNewFESpace, pControlPointCoordinatesGrid);
        mpGridFunctions["CONTROL_POINT_COORDINATES"] = pNewCoordinatesGridFunc;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create and add the grid function. This function will create the new FESpace based on the original FESpace of the control grid and the weights, and then assign to the new grid function.
    /// One must not use this function for the ControlPoint data type.
    template<typename TOtherDataType>
    typename GridFunction<TDim, TLocalCoordinateType, TOtherDataType>::Pointer CreateGridFunction(typename ControlGrid<TOtherDataType>::Pointer pControlGrid)
    {
        CheckSize(*pControlGrid, __FUNCTION__);
        typename FESpaceType::Pointer pNewFESpace = WeightedFESpaceType::Create(mpFESpace, this->GetControlWeights());
        typename GridFunction<TDim, TLocalCoordinateType, TOtherDataType>::Pointer pNewGridFunc = GridFunction<TDim, TLocalCoordinateType, TOtherDataType>::Create(pNewFESpace, pControlGrid);
        mpGridFunctions[pControlGrid->Name()] = pNewGridFunc;
        return pNewGridFunc;
    }

    /// Create and add the grid function
    template<class TVariableType>
    typename GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type>::Pointer CreateGridFunction(const TVariableType& rVariable,
            typename ControlGrid<typename TVariableType::Type>::Pointer pControlGrid)
    {
        pControlGrid->SetName(rVariable.Name());
        return this->CreateGridFunction<typename TVariableType::Type>(pControlGrid);
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type>::Pointer pGetGridFunction(const TVariableType& rVariable)
    {
        typedef typename GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (GridFunctionContainerType::iterator it = mpGridFunctions.begin(); it != mpGridFunctions.end(); ++it)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(it->second);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                {
                    return pGridFunc;
                }
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        KRATOS_ERROR << "The grid function with control grid " << rVariable.Name() << " does not exist in the database of patch " << Id();
    }

    /// Get the grid function
    template<class TVariableType>
    typename GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type>::ConstPointer pGetGridFunction(const TVariableType& rVariable) const
    {
        typedef typename GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type>::Pointer GridFunctionPointerType;
        for (GridFunctionContainerType::const_iterator it = mpGridFunctions.begin(); it != mpGridFunctions.end(); ++it)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(it->second);
                if (pGridFunc->pControlGrid()->Name() == rVariable.Name())
                {
                    return pGridFunc;
                }
            }
            catch (boost::bad_any_cast& e)
            {
                continue;
            }
        }
        // shall not come here
        KRATOS_ERROR << "The grid function with control grid " << rVariable.Name() << " does not exist in the database of patch " << Id();
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
            {
                return true;
            }
        return false;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute a rough estimation of the local coordinates of a point by sampling technique
    void Predict(const array_1d<CoordinateType, 3>& point, array_1d<LocalCoordinateType, 3>& xi, const std::vector<int>& nsampling,
                 const array_1d<LocalCoordinateType, 3>& xi_min, const array_1d<LocalCoordinateType, 3>& xi_max) const
    {
        typename GridFunction<TDim, TLocalCoordinateType, array_1d<CoordinateType, 3> >::ConstPointer pGridFunc = this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT_COORDINATES));
        pGridFunc->Predict(point, xi, nsampling, xi_min, xi_max);
    }

    /// Compute a rough estimation of the local coordinates of a point by sampling technique
    void Predict(const array_1d<CoordinateType, 3>& point, array_1d<LocalCoordinateType, 3>& xi, const std::vector<int>& nsampling) const
    {
        typename GridFunction<TDim, TLocalCoordinateType, array_1d<CoordinateType, 3> >::ConstPointer pGridFunc = this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT_COORDINATES));
        pGridFunc->Predict(point, xi, nsampling);
    }

    /// Compute the local coordinates of a point
    int LocalCoordinates(const array_1d<CoordinateType, 3>& point, array_1d<LocalCoordinateType, 3>& xi) const
    {
        typename GridFunction<TDim, TLocalCoordinateType, array_1d<CoordinateType, 3> >::ConstPointer pGridFunc = this->pGetGridFunction(VARSEL(TCoordinateType, CONTROL_POINT_COORDINATES));
        int error_code = pGridFunc->LocalCoordinates(point, xi, mLocalSearchMaxIters, mLocalSearchTolerance);
        bool is_inside = this->pFESpace()->IsInside(std::vector<LocalCoordinateType> {xi[0], xi[1], xi[2]});
        if (!is_inside)
        {
            return 2;
        }
        else
        {
            return error_code;
        }
    }

    /// Check if the point is inside the patch
    /// This subroutine requires xi0, which is a prediction of the projected local point. This has to be determined, i.e. using a sampling technique.
    bool IsInside(const array_1d<CoordinateType, 3>& point, const array_1d<LocalCoordinateType, 3>& xi0) const
    {
        array_1d<LocalCoordinateType, 3> xi;
        noalias(xi) = xi0;
        int stat = this->LocalCoordinates(point, xi);

        if (stat == 0)
            return this->pFESpace()->IsInside(std::vector<LocalCoordinateType> {xi[0], xi[1], xi[2]});
        else
        {
            return false;
        }
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
            KRATOS_ERROR << "The patch must have an Id";
        }

        if (pControlPointGridFunction() != NULL)
            if (pControlPointGridFunction()->pControlGrid()->Size() != this->TotalNumber())
                KRATOS_ERROR << "The control point grid is incompatible";

                DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_ERROR << "The double variable grid " << (*it)->pControlGrid()->Name() << " is incompatible";
                return false;
            }
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_ERROR << "The array_1d variable grid " << (*it)->pControlGrid()->Name() << " is incompatible";
                return false;
            }
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Size() != this->TotalNumber())
            {
                KRATOS_ERROR << "The vector variable grid " << (*it)->pControlGrid()->Name() <<" is incompatible";
                return false;
            }
        }

        return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Construct the boundary patch based on side
    virtual typename BoundaryPatchType::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    {
        // construct the boundary FESpace
        typename BoundaryFESpaceType::Pointer pBFESpace = this->pFESpace()->ConstructBoundaryFESpace(side);

        return this->ConstructBoundaryPatch(pBFESpace);
    }

    /// Construct the boundary patch based on side and local relative configuration
    virtual typename BoundaryPatchType::Pointer ConstructBoundaryPatch(const BoundarySide& side,
            const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const
    {
        // construct the boundary FESpace
        typename BoundaryFESpaceType::Pointer pBFESpace = this->pFESpace()->ConstructBoundaryFESpace(side, local_parameter_map, directions);

        return this->ConstructBoundaryPatch(pBFESpace);
    }

    /// Construct the boundary patch based on boundary FESpace
    typename BoundaryPatchType::Pointer ConstructBoundaryPatch(const typename BoundaryFESpaceType::Pointer pBFESpace) const
    {
        // construct the boundary patch based on boundary FESpace
        typename BoundaryPatchType::Pointer pBPatch = typename BoundaryPatchType::Pointer(new BoundaryPatchType (-1));
        pBPatch->SetFESpace(pBFESpace);

        // transfer the control values
        // std::cout << "Control point grid function " << this->pControlPointGridFunction()->pControlGrid()->Name() << " will be constructed" << std::endl;
        typename ControlGrid<ControlPointType>::Pointer pBoundaryControlPointGrid = ControlGridUtility::ExtractSubGrid<TDim, ControlPointType>(this->pControlPointGridFunction()->pControlGrid(), *(this->pFESpace()), *pBFESpace);
        pBPatch->CreateControlPointGridFunction(pBoundaryControlPointGrid);

        // transfer other values
        DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            // std::cout << "Double grid function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<TDataType>::Pointer pBoundaryDoubleControlGrid = ControlGridUtility::ExtractSubGrid<TDim, TDataType>((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<TDataType>(pBoundaryDoubleControlGrid);
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") { continue; }
            // std::cout << "Array1D function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<array_1d<TDataType, 3> >::Pointer pBoundaryArray1DControlGrid = ControlGridUtility::ExtractSubGrid<TDim, array_1d<TDataType, 3> >((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<array_1d<TDataType, 3> >(pBoundaryArray1DControlGrid);
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            // std::cout << "Vector grid function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<VectorType>::Pointer pBoundaryVectorControlGrid = ControlGridUtility::ExtractSubGrid<TDim, VectorType>((*it)->pControlGrid(), *(this->pFESpace()), *pBFESpace);
            pBPatch->template CreateGridFunction<VectorType>(pBoundaryVectorControlGrid);
        }

        return pBPatch;
    }

    /// Construct the sliced patch on a specific direction
    virtual typename BoundaryPatchType::Pointer ConstructSlicedPatch(int idir, LocalCoordinateType xi) const
    {
        typename BoundaryPatchType::Pointer pSPatch = typename BoundaryPatchType::Pointer(new BoundaryPatchType (-1));

        // construct the sliced FESpace
        typename BoundaryFESpaceType::Pointer pSFESpace = this->pFESpace()->ConstructSlicedFESpace(idir, xi);
        pSPatch->SetFESpace(pSFESpace);

        // transfer the control values
        std::cout << "Control point grid function " << this->pControlPointGridFunction()->pControlGrid()->Name() << " will be constructed" << std::endl;
        typename ControlGrid<ControlPointType>::Pointer pSlicedControlPointGrid = ControlGridUtility::ComputeSlicedGrid<TDim, ControlPointType>(this->pControlPointGridFunction()->pControlGrid(), *(this->pFESpace()), idir, xi);
        pSPatch->CreateControlPointGridFunction(pSlicedControlPointGrid);

        // transfer other values
        DoubleGridFunctionContainerType DoubleGridFunctions_ = this->DoubleGridFunctions();
        for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            std::cout << "Double grid function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<TDataType>::Pointer pSlicedDoubleControlGrid = ControlGridUtility::ComputeSlicedGrid<TDim, TDataType>((*it)->pControlGrid(), *(this->pFESpace()), idir, xi);
            pSPatch->template CreateGridFunction<TDataType>(pSlicedDoubleControlGrid);
        }

        Array1DGridFunctionContainerType Array1DGridFunctions_ = this->Array1DGridFunctions();
        for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") { continue; }
            std::cout << "Array1D function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<array_1d<TDataType, 3> >::Pointer pSlicedArray1DControlGrid = ControlGridUtility::ComputeSlicedGrid<TDim, array_1d<TDataType, 3> >((*it)->pControlGrid(), *(this->pFESpace()), idir, xi);
            pSPatch->template CreateGridFunction<array_1d<TDataType, 3> >(pSlicedArray1DControlGrid);
        }

        VectorGridFunctionContainerType VectorGridFunctions_ = this->VectorGridFunctions();
        for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            std::cout << "Vector grid function " << (*it)->pControlGrid()->Name() << " will be constructed" << std::endl;
            typename ControlGrid<VectorType>::Pointer pSlicedVectorControlGrid = ControlGridUtility::ComputeSlicedGrid<TDim, VectorType>((*it)->pControlGrid(), *(this->pFESpace()), idir, xi);
            pSPatch->template CreateGridFunction<VectorType>(pSlicedVectorControlGrid);
        }

        return pSPatch;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Search for the neighbor
    typename PatchType::Pointer pNeighbor(const BoundarySide& side)
    {
        for (interface_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if ((*it)->Side1() == side)
            {
                return (*it)->pPatch2();
            }
        }
        return NULL;
    }

    /// Search for the neighbor
    typename PatchType::ConstPointer pNeighbor(const BoundarySide& side) const
    {
        for (interface_const_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if ((*it)->Side1() == side)
            {
                return (*it)->pPatch2();
            }
        }
        return NULL;
    }

    /// Search for the boundary side (in the current patch) of the neighbor patch, if it exists
    int FindBoundarySide(typename PatchType::ConstPointer pPatch) const
    {
        for (interface_const_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if ((*it)->pPatch2() == pPatch)
            {
                return (*it)->Side1();
            }
        }
        return -1;
    }

    /// Add an interface to the patch
    void AddInterface(typename PatchInterfaceType::Pointer pInterface)
    {
        if (&(*(pInterface->pPatch1())) != this)
        {
            std::cout << "WARNING: the patch 1 of the interface is not the same as this patch, skipped." << std::endl;
            return;
        }

        for (interface_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if ((*it) == pInterface)
            {
                std::cout << "WARNING: the interface exist in this patch, skipped." << std::endl;
                return;
            }
        }

        mpInterfaces.push_back(pInterface);
    }

    /// Remove an interface from the patch
    void RemoveInterface(typename PatchInterfaceType::Pointer pInterface)
    {
        if (&(*(pInterface->pPatch1())) != this)
        {
            std::cout << "WARNING: the patch 1 of the interface is not the same as this patch, skipped." << std::endl;
            return;
        }

        for (interface_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if ((*it) == pInterface)
            {
                std::cout << "REMARK: Interface " << (*it) << " is removed from patch " << Id() << std::endl;
                mpInterfaces.erase(it);
                return;
            }
        }
    }

    /// Remove all interfaces associated with patch
    void ClearInterface()
    {
        mpInterfaces.clear();
    }

    /// Get the number of interfaces
    std::size_t NumberOfInterfaces() const
    {
        return mpInterfaces.size();
    }

    /// Interface iterators
    interface_iterator InterfaceBegin() {return mpInterfaces.begin();}
    interface_const_iterator InterfaceBegin() const {return mpInterfaces.begin();}
    interface_iterator InterfaceEnd() {return mpInterfaces.end();}
    interface_const_iterator InterfaceEnd() const {return mpInterfaces.end();}

    /// Get the interface
    typename PatchInterfaceType::Pointer pInterface(std::size_t i)
    {
        std::size_t cnt = 0;
        for (interface_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if (cnt == i)
            {
                return (*it);
            }
            else
            {
                ++cnt;
            }
        }
        return NULL;
    }

    typename PatchInterfaceType::ConstPointer pInterface(std::size_t i) const
    {
        std::size_t cnt = 0;
        for (interface_const_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            if (cnt == i)
            {
                return (*it);
            }
            else
            {
                ++cnt;
            }
        }
        return NULL;
    }

    /// Get/Set the parent multipatch
    void pSetParentMultiPatch(typename MultiPatchType::Pointer pPatch) {mpParentMultiPatch = pPatch;}
    MultiPatchType& ParentMultiPatch() {return *pParentMultiPatch();}
    const MultiPatchType& ParentMultiPatch() const {return *pParentMultiPatch();}
    typename MultiPatchType::Pointer pParentMultiPatch() {return mpParentMultiPatch.lock();}
    const typename MultiPatchType::Pointer pParentMultiPatch() const {return mpParentMultiPatch.lock();}

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Generate topology data to visualize with GLVis
    void GenerateTopologyData(std::size_t& starting_vertex_id,
                              std::vector<vertex_t>& vertices,
                              std::vector<edge_t>& edges,
                              std::vector<face_t>& faces,
                              std::vector<volume_t>& volumes,
                              std::size_t& starting_knotv_id,
                              std::vector<std::size_t>& knotv ) const
    {
        if constexpr (TDim == 1)
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
        else if constexpr (TDim == 2)
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
        else if constexpr (TDim == 3)
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
            KRATOS_ERROR << "Not implemented for " << TDim;
        }
    }

    /// Get the bounding box of the patch
    /// The arrangement is [x_min, x_max, y_min, y_max, z_min, z_max]
    /// This function takes advantage of the convex hull properties of a NURBS patch
    void GetBoundingBox(std::vector<CoordinateType>& bounding_box) const
    {
        if (bounding_box.size() != 6)
        {
            bounding_box.resize(6);
        }

        CoordinateType& x_min = bounding_box[0];
        CoordinateType& x_max = bounding_box[1];
        CoordinateType& y_min = bounding_box[2];
        CoordinateType& y_max = bounding_box[3];
        CoordinateType& z_min = bounding_box[4];
        CoordinateType& z_max = bounding_box[5];

        x_min = 1.0e99; y_min = 1.0e99; z_min = 1.0e99;
        x_max = -1.0e99; y_max = -1.0e99; z_max = -1.0e99;

        typename ControlGrid<ControlPointType>::ConstPointer pControlPointGrid = pControlPointGridFunction()->pControlGrid();
        CoordinateType x, y, z;
        for (std::size_t i = 0; i < pControlPointGrid->size(); ++i)
        {
            x = (*pControlPointGrid)[i].X();
            y = (*pControlPointGrid)[i].Y();
            z = (*pControlPointGrid)[i].Z();

            if (x < x_min) { x_min = x; }
            if (x > x_max) { x_max = x; }
            if (y < y_min) { y_min = y; }
            if (y > y_max) { y_max = y; }
            if (z < z_min) { z_min = z; }
            if (z > z_max) { z_max = z; }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare two patches in terms of its parametric information. The grid function data, including control points, are not checked.
    virtual bool IsCompatible(const PatchType& rOtherPatch) const
    {
        return *(this->pFESpace()) == *(rOtherPatch.pFESpace());
    }

    /// Compare two patches in terms of parametric information and control points.
    bool IsEquivalent(const PatchType& rOtherPatch, const int echo_level = 0, const CoordinateValueType dist_tol = DISTANCE_TOLERANCE) const
    {
        if (!this->IsCompatible(rOtherPatch))
        {
            if (echo_level > 0)
                std::cout << "Patch equivalence check: The two patches are not compatible" << std::endl;
            return false;
        }

        // compare the control points

        typename ControlGrid<ControlPointType>::ConstPointer pThisControlPointGrid = pControlPointGridFunction()->pControlGrid();
        typename ControlGrid<ControlPointType>::ConstPointer pOtherControlPointGrid = rOtherPatch.pControlPointGridFunction()->pControlGrid();

        if (pThisControlPointGrid->size() != pOtherControlPointGrid->size())
        {
            if (echo_level > 0)
                std::cout << "Patch equivalence check: The two control grid have different size" << std::endl;
            return false;
        }

        for (std::size_t i = 0; i < pThisControlPointGrid->size(); ++i)
        {
            const ControlPointType& p1 = pThisControlPointGrid->GetData(i);
            const ControlPointType& p2 = pOtherControlPointGrid->GetData(i);

            const CoordinateType dist = p1.Distance(p2);
            if (std::abs(dist) > dist_tol)
            {
                if (echo_level > 0)
                   std::cout << "Patch equivalence check: The control point " << p1 << " and " << p2
                             << " are not the same. The difference = " << dist << " > " << dist_tol
                             << std::endl;
                return false;
            }
        }

        return true;
    }

    /// Compare two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const PatchType& rOtherPatch) const
    {
        if (!this->IsEquivalent(rOtherPatch))
        {
            return false;
        }

        // TODO compare the grid function values

        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const PatchType& rOther)
    {
        return (Id() == rOther.Id()) && this->IsSame(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Type() << ", Id = " << Id() << ", Addr = " << this;
    }

    void PrintData(std::ostream& rOStream) const override
    {
        if (pFESpace() != NULL)
        {
            rOStream << *pFESpace() << std::endl;
        }

        if (pControlPointGridFunction() != NULL)
        {
            rOStream << *(pControlPointGridFunction()->pControlGrid()) << std::endl;
        }

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
        for (interface_const_iterator it = InterfaceBegin(); it != InterfaceEnd(); ++it)
        {
            rOStream << "  ";
            (*it)->PrintInfo(rOStream);
            rOStream << std::endl;
        }
    }

private:

    std::string mPrefix;
    int mLayerIndex;

    // FESpace contains the shape function information and various information with regards to the functional space.
    // Because the control point grid is in homogeneous coordinates, the FESpace shall be an unweighted spaces
    typename FESpaceType::Pointer mpFESpace;

    // container to contain all the grid functions
    GridFunctionContainerType mpGridFunctions; // using boost::any to store pointer to grid function

    unsigned int mLocalSearchMaxIters;
    LocalCoordinateType mLocalSearchTolerance;

    /**
     * interface data
     */
    InterfaceContainerType mpInterfaces;

    /**
     * pointer to parent multipatch
     */
    typename MultiPatchType::WeakPointer mpParentMultiPatch;

    /// Empty Constructor for serializer
    Patch() : IndexedObject(0), mpFESpace(NULL) {}

    /// Serializer
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
    }

    void load(Serializer& rSerializer) override
    {
    }

    /// Auxiliary
    template<class TGridFunctionType>
    void CheckSize(const TGridFunctionType& rGrid, const std::string& source) const
    {
        if (rGrid.Size() != this->TotalNumber())
        {
            KRATOS_ERROR << "The size of grid function (" << rGrid.size() << ") is not compatible with the current number of control values (" << this->TotalNumber()
                         << ") of patch " << Id() << ". Error at " << source;
        }
    }

    /// Helper to extract the grid functions out from boost::any
    template<class TContainerType>
    TContainerType ExtractGridFunctions(const GridFunctionContainerType& pGridFunctions) const
    {
        TContainerType GridFuncs;

        typedef typename TContainerType::value_type GridFunctionPointerType;

        for (GridFunctionContainerType::const_iterator it = pGridFunctions.begin(); it != pGridFunctions.end(); ++it)
        {
            try
            {
                GridFunctionPointerType pGridFunc = boost::any_cast<GridFunctionPointerType>(it->second);
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
    std::vector<TVariableType*> ExtractVariables(const GridFunctionContainerType& pGridFunctions) const
    {
        typedef GridFunction<TDim, TLocalCoordinateType, typename TVariableType::Type> GridFunctionType;
        typedef std::vector<typename GridFunctionType::Pointer> GridFunctionVectorContainerType;
        GridFunctionVectorContainerType GridFuncs = this->ExtractGridFunctions<GridFunctionVectorContainerType>(pGridFunctions);

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
}; // class Patch

/**
 * Template specific instantiation for null-D patch to terminate the compilation.
 * In fact, null-D patch is a vertex
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class Patch<0, TLocalCoordinateType, TCoordinateType, TDataType> : public IndexedObject, public Flags
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const Patch> ConstPointer;
#endif

    /// Type definitions
    typedef TCoordinateType CoordinateType;
    typedef typename DataTypeToValueType<TCoordinateType>::value_type CoordinateValueType;
    typedef CoordinateValueType WeightType;

    typedef ControlPoint<CoordinateType, WeightType> ControlPointType;
    typedef FESpace<0, TLocalCoordinateType> FESpaceType;
    typedef GridFunction<0, TLocalCoordinateType, ControlPointType> ControlPointGridFunctionType;

    typedef Patch<0, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;

    /// Constants
    static constexpr CoordinateValueType DISTANCE_TOLERANCE = 1e-13;

    /// Default constructor
    Patch() : IndexedObject(0) {}

    /// Constructor with id
    Patch(std::size_t Id) : IndexedObject(Id) {}

    /// Destructor
    virtual ~Patch() {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpaceType::Pointer pFESpace) {mpFESpace = pFESpace;}

    /// Get the number of basis functions defined over the patch
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual std::size_t Order(std::size_t i) const
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
    typename ControlPointGridFunctionType::Pointer CreateControlPointGridFunction(typename ControlGrid<ControlPointType>::Pointer pControlPointGrid)
    {
        pControlPointGrid->SetName("CONTROL_POINT");
        typename ControlPointGridFunctionType::Pointer pNewGridFunc = ControlPointGridFunctionType::Create(mpFESpace, pControlPointGrid);
        mpControlPointGridFunc = pNewGridFunc;

        return pNewGridFunc;
    }

    /// Get the control point grid function
    ControlPointGridFunctionType& ControlPointGridFunction() {return *mpControlPointGridFunc;}

    /// Get the control point grid function
    const ControlPointGridFunctionType& ControlPointGridFunction() const {return *mpControlPointGridFunc;}

    /// Create and add the grid function
    template<typename TOtherDataType>
    typename GridFunction<0, TLocalCoordinateType, TOtherDataType>::Pointer CreateGridFunction(typename ControlGrid<TOtherDataType>::Pointer pControlGrid)
    {
        return NULL;
    }

    /// Get the FESpace pointer
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Get the FESpace pointer
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Compare between two patches in terms of parametric information and grid function data, including the control points.
    bool IsSame(const PatchType& rOtherPatch) const
    {
        return true;
    }

    /// Overload comparison operator
    virtual bool operator==(const PatchType& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const PatchType& rPatch1, const BoundarySide& side1,
                                           const PatchType& rPatch2, const BoundarySide& side2)
    {
        return true;
    }

    // /// Construct the boundary patch based on side
    // virtual typename Patch<-1>::Pointer ConstructBoundaryPatch(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Compare two patches in terms of parametric information and control points.
    bool IsEquivalent(const PatchType& rOtherPatch, const int echo_level = 0, const CoordinateValueType dist_tol = DISTANCE_TOLERANCE) const
    {
        // compare the control points

        typename ControlGrid<ControlPointType>::ConstPointer pThisControlPointGrid = ControlPointGridFunction().pControlGrid();
        typename ControlGrid<ControlPointType>::ConstPointer pOtherControlPointGrid = rOtherPatch.ControlPointGridFunction().pControlGrid();

        if (pThisControlPointGrid->size() != pOtherControlPointGrid->size())
        {
            return false;
        }

        for (std::size_t i = 0; i < pThisControlPointGrid->size(); ++i)
        {
            const ControlPointType& p1 = pThisControlPointGrid->GetData(i);
            const ControlPointType& p2 = pOtherControlPointGrid->GetData(i);

            const double dist = p1.Distance(p2);
            if (dist > dist_tol)
            {
                return false;
            }
        }

        return true;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Patch<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    typename ControlPointGridFunctionType::Pointer mpControlPointGridFunc;
    typename FESpaceType::Pointer mpFESpace;
};

/**
 * Template specific instantiation for -1-D patch to terminate the compilation.
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class Patch<-1, TLocalCoordinateType, TCoordinateType, TDataType> : public IndexedObject, public Flags
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const Patch> ConstPointer;
#endif

    /// Type definitions
    typedef TCoordinateType CoordinateType;
    typedef typename DataTypeToValueType<TCoordinateType>::value_type CoordinateValueType;
    typedef CoordinateValueType WeightType;

    typedef ControlPoint<CoordinateType, WeightType> ControlPointType;
    typedef FESpace<-1, TLocalCoordinateType> FESpaceType;

    typedef Patch<-1, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;

    /// Default constructor
    Patch() : IndexedObject(0) {}

    /// Constructor with id
    Patch(std::size_t Id) : IndexedObject(Id) {}

    /// Destructor
    ~Patch() override {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpaceType::Pointer pFESpace) {}

    /// Get the number of basis functions defined over the patch
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual std::size_t Order(std::size_t i) const
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
    virtual bool operator==(const PatchType& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const PatchType& rPatch1, const BoundarySide& side1,
                                           const PatchType& rPatch2, const BoundarySide& side2)
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
};

/**
 * Template specific instantiation for -2-D patch to terminate the compilation.
 */
template<typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
class Patch<-2, TLocalCoordinateType, TCoordinateType, TDataType> : public IndexedObject, public Flags
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Patch);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const Patch> ConstPointer;
#endif

    /// Type definitions
    typedef TCoordinateType CoordinateType;
    typedef typename DataTypeToValueType<TCoordinateType>::value_type CoordinateValueType;
    typedef CoordinateValueType WeightType;

    typedef ControlPoint<CoordinateType, WeightType> ControlPointType;
    typedef FESpace<-2, TLocalCoordinateType> FESpaceType;

    typedef Patch<-2, TLocalCoordinateType, TCoordinateType, TDataType> PatchType;

    /// Default constructor
    Patch() : IndexedObject(0) {}

    /// Constructor with id
    Patch(std::size_t Id) : IndexedObject(Id) {}

    /// Destructor
    ~Patch() override {}

    /// Set the FESpace for the patch
    void SetFESpace(typename FESpaceType::Pointer pFESpace) {}

    /// Get the number of basis functions defined over the patch
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the patch in specific direction
    virtual std::size_t Order(std::size_t i) const
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
    virtual bool operator==(const PatchType& rOther)
    {
        return Id() == rOther.Id();
    }

    /// Check the compatibility between boundaries of two patches
    static bool CheckBoundaryCompatibility(const PatchType& rPatch1, const BoundarySide& side1,
                                           const PatchType& rPatch2, const BoundarySide& side2)
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
};

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TCoordinateType, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const Patch<TDim, TLocalCoordinateType, TCoordinateType, TDataType>& rThis)
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PATCH_H_INCLUDED defined
