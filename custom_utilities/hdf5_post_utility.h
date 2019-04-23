//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Jan 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_HDF5_POST_UTILITY_H_INCLUDED )
#define  KRATOS_HDF5_POST_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>
#include "boost/progress.hpp"
#include "H5Cpp.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"


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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/*** Detail class definition.
 */
class HDF5PostUtility
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

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename Node<3>::DofsContainerType DofsContainerType;

    typedef typename Node<3>::Pointer NodeType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SerialSparseSpaceType;

    typedef UblasSpace<double, Matrix, Vector> SerialDenseSpaceType;

    typedef LinearSolver<SerialSparseSpaceType, SerialDenseSpaceType> LinearSolverType;

    typedef std::size_t IndexType;

    /// Pointer definition of HDF5PostUtility
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(HDF5PostUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HDF5PostUtility(const std::string h5_filename)
    {
        mpFile = Kratos::shared_ptr<H5::H5File>(new H5::H5File(h5_filename, H5F_ACC_TRUNC));
    }

    /// Constructor with access mode
    HDF5PostUtility(const std::string h5_filename, const std::string AccessMode)
    {
        unsigned int mode = H5F_ACC_TRUNC;

        if(AccessMode == std::string("Truncation"))
        {
            mode = H5F_ACC_TRUNC;
        }
        else if(AccessMode == std::string("Read-Only"))
        {
            mode = H5F_ACC_RDONLY;
        }
        else if(AccessMode == std::string("Read-Write"))
        {
            mode = H5F_ACC_RDWR;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "This access mode is not supported:", AccessMode)

        mpFile = Kratos::shared_ptr<H5::H5File>(new H5::H5File(h5_filename, mode));
    }

    /// Destructor.
    virtual ~HDF5PostUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Write all nodes of the model_part to the HDF5 datafile
     */
    void WriteNodes(ModelPart& r_model_part)
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Node_t {
	        int    id;
	        double x;
	        double y;
	        double z;
        } Node_t;

        Node_t* nodes = new Node_t[pNodes.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            nodes[cnt].id = (*it)->Id();
            nodes[cnt].x = (*it)->X0();
            nodes[cnt].y = (*it)->Y0();
            nodes[cnt].z = (*it)->Z0();
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype( sizeof(Node_t) );
        mtype.insertMember( "id", HOFFSET(Node_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember( "x",  HOFFSET(Node_t, x),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember( "y",  HOFFSET(Node_t, y),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember( "z",  HOFFSET(Node_t, z),  H5::PredType::NATIVE_DOUBLE);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet("Nodes", mtype, space));

        /*
         * Write data to the dataset;
         */
        dataset->write(nodes, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete nodes;
    }

    template<class TDataType>
    void WriteNodalResults(const Variable<TDataType>& rThisVariable, ModelPart& r_model_part)
    {
        WriteNodalResults_(rThisVariable, r_model_part);
    }

    template<class TDataType>
    void ReadNodalResults(const Variable<TDataType>& rThisVariable, ModelPart& r_model_part)
    {
        ReadNodalResults_(rThisVariable, r_model_part);
    }

    template<class TDataType>
    void WriteElementalData(const Variable<TDataType>& rThisVariable, ModelPart& r_model_part)
    {
        WriteElementalData_(rThisVariable, r_model_part);
    }

    template<class TDataType>
    void ReadElementalData(const Variable<TDataType>& rThisVariable, ModelPart& r_model_part, bool allow_unequal = true)
    {
        ReadElementalData_(rThisVariable, r_model_part, allow_unequal);
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
        buffer << "HDF5PostUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HDF5PostUtility";
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
    Kratos::shared_ptr<H5::H5File> mpFile;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /*****************************************************
           FUNCTIONS TO WRITE DATA TO NODES
    *****************************************************/
    void WriteNodalResults_(
        const Variable<double>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double v;
        } Data_t;

        Data_t* data = new Data_t[pNodes.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            data[cnt].id = (*it)->Id();
            data[cnt].v = (*it)->operator[](rThisVariable);
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t));
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v",  HOFFSET(Data_t, v),  H5::PredType::NATIVE_DOUBLE);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));

        /*
         * Write data to the dataset;
         */
        dataset->write(data, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete data;
    }

    void WriteNodalResults_(
        const Variable<array_1d<double, 3> >& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double d1;
	        double d2;
	        double d3;
        } Data_t;

        Data_t* data = new Data_t[pNodes.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            data[cnt].id = (*it)->Id();
            array_1d<double, 3>& v = (*it)->GetSolutionStepValue(rThisVariable);
            data[cnt].d1 = v[0];
            data[cnt].d2 = v[1];
            data[cnt].d3 = v[2];
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, d1),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, d2),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, d3),  H5::PredType::NATIVE_DOUBLE);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));

        /*
         * Write data to the dataset;
         */
        dataset->write(data, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete data;
    }

    void WriteNodalResults_(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();
        int len = (*(pNodes.ptr_begin()))->GetSolutionStepValue(rThisVariable).size();

        // I have to do this since it is not impossible AFAIK to create a struct with double pointer with HDF5/C++. If yes, is the performance OK?
        if(len == 3)
        {
            WriteNodalResults_Vector_3(rThisVariable, pModelPart);
        }
        else if(len == 6)
        {
            WriteNodalResults_Vector_6(rThisVariable, pModelPart);
        }
        else
            std::cout << "Vector length of size " << len
                      << " is not supported, skipping variable "
                      << rThisVariable.Name() << std::endl;
    }

    void WriteNodalResults_Vector_3(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double d1;
	        double d2;
	        double d3;
        } Data_t;

        Data_t* data = new Data_t[pNodes.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            data[cnt].id = (*it)->Id();
            Vector& v = (*it)->GetSolutionStepValue(rThisVariable);
            data[cnt].d1 = v[0];
            data[cnt].d2 = v[1];
            data[cnt].d3 = v[2];
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, d1),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, d2),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, d3),  H5::PredType::NATIVE_DOUBLE);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));

        /*
         * Write the attribute
         */
        int attr_data[] = { 3 };
        hsize_t dims[] = { 1 };
        H5::DataSpace attr_dataspace = H5::DataSpace (1, dims);
        H5::Attribute attribute = dataset->createAttribute("LEN", H5::PredType::STD_I32BE, attr_dataspace);
        attribute.write(H5::PredType::NATIVE_INT, attr_data);

        /*
         * Write data to the dataset;
         */
        dataset->write(data, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete data;
    }

    void WriteNodalResults_Vector_6(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double d1;
	        double d2;
	        double d3;
	        double d4;
	        double d5;
	        double d6;
        } Data_t;

        Data_t* data = new Data_t[pNodes.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            data[cnt].id = (*it)->Id();
            Vector& v = (*it)->GetSolutionStepValue(rThisVariable);
            data[cnt].d1 = v[0];
            data[cnt].d2 = v[1];
            data[cnt].d3 = v[2];
            data[cnt].d4 = v[3];
            data[cnt].d5 = v[4];
            data[cnt].d6 = v[5];
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, d1),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, d2),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, d3),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v4",  HOFFSET(Data_t, d4),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v5",  HOFFSET(Data_t, d5),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v6",  HOFFSET(Data_t, d6),  H5::PredType::NATIVE_DOUBLE);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));

        /*
         * Write the attribute
         */
        int attr_data[] = { 6 };
        hsize_t dims[] = { 1 };
        H5::DataSpace attr_dataspace = H5::DataSpace (1, dims);
        H5::Attribute attribute = dataset->createAttribute("LEN", H5::PredType::STD_I32BE, attr_dataspace);
        attribute.write(H5::PredType::NATIVE_INT, attr_data);

        /*
         * Write data to the dataset;
         */
        dataset->write(data, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete data;
    }

    /*****************************************************
                FUNCTIONS TO WRITE ELEMENTAL DATA
    *****************************************************/
    void WriteElementalData_(
        const Variable<bool>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        ElementsArrayType& pElements = r_model_part.Elements();

        typedef struct Data_t {
	        int id;
	        int v;
        } Data_t;

        Data_t* data = new Data_t[pElements.size()];

        /*
         * Initialize the data
         */
        int cnt = 0;
        for(typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            data[cnt].id = (*it)->Id();
            if((*it)->GetValue(rThisVariable) == true)
                data[cnt].v = 1;
            else
                data[cnt].v = 0;
            ++cnt;
        }

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t));
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v",  HOFFSET(Data_t, v),  H5::PredType::NATIVE_INT);

        /*
         * Create the data space.
         */
        hsize_t dim[] = {pElements.size()};   /* Dataspace dimensions */
        H5::DataSpace space(1, dim);

        /*
         * Create the dataset.
         */
        H5::DataSet* dataset;
        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));

        /*
         * Write data to the dataset;
         */
        dataset->write(data, mtype);

        /*
         * Release resources
         */
        delete dataset;
        delete data;
    }


    /*****************************************************
             FUNCTIONS TO READ DATA FROM NODES
    *****************************************************/
    int ReadNodalResults_(
        const Variable<double>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double v;
        } Data_t;

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v",  HOFFSET(Data_t, v),  H5::PredType::NATIVE_DOUBLE);

        try
        {
            /*
             * Open the dataset
             */
            H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

            /*
             * Get dataspace of the dataset.
             */
            H5::DataSpace dataspace = dataset.getSpace();

            /*
             * Get the number of dimensions in the dataspace.
             */
            int rank = dataspace.getSimpleExtentNdims();
//            KRATOS_WATCH(rank)

            /*
             * Get the dimension size of each dimension in the dataspace and
             * do the bound check.
             */
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

            if(dims_out[0] != pNodes.size())
                KRATOS_THROW_ERROR(std::logic_error, "The dimension of provided data is not compatible:", dims_out[0])

            /*
             * Output buffer
             */
            Data_t* data = new Data_t[pNodes.size()];

            /*
             * Read data from hyperslab in the file into the hyperslab in
             * memory and display the data.
             */
            dataset.read(data, mtype, dataspace, dataspace);

            double tmp;
            for(unsigned int i = 0; i < dims_out[0]; ++i)
            {
                tmp = data[i].v;
                pNodes[data[i].id].GetSolutionStepValue(rThisVariable) = tmp;
            }

            /*
             * Release memory
             */
            delete data;
        }
        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            error.printError();
            return -1;
        }
        // catch failure caused by the DataSpace operations
       catch(H5::DataSpaceIException error)
       {
          error.printError();
          return -2;
       }
       // catch failure caused by the DataSpace operations
       catch(H5::DataTypeIException error)
       {
          error.printError();
          return -3;
       }
    }

    int ReadNodalResults_(
        const Variable<array_1d<double, 3> >& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double v[3];
        } Data_t;

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, v[0]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, v[1]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, v[2]),  H5::PredType::NATIVE_DOUBLE);

        try
        {
            /*
             * Open the dataset
             */
            H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

            /*
             * Get dataspace of the dataset.
             */
            H5::DataSpace dataspace = dataset.getSpace();

            /*
             * Get the number of dimensions in the dataspace.
             */
            int rank = dataspace.getSimpleExtentNdims();
//            KRATOS_WATCH(rank)

            /*
             * Get the dimension size of each dimension in the dataspace and
             * do the bound check.
             */
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

            if(dims_out[0] != pNodes.size())
                KRATOS_THROW_ERROR(std::logic_error, "The dimension of provided data is not compatible:", dims_out[0])

            /*
             * Output buffer
             */
            Data_t* data = new Data_t[pNodes.size()];

            /*
             * Read data from hyperslab in the file into the hyperslab in
             * memory and display the data.
             */
            dataset.read(data, mtype, dataspace, dataspace);

            //display data
//            for(int i = 0; i < pNodes.size(); ++i)
//                std::cout << "node " << data[i].id
//                          << " d1: " << data[i].v[0]
//                          << " d2: " << data[i].v[1]
//                          << " d3: " << data[i].v[2]
//                          << std::endl;
            array_1d<double, 3> tmp;
            for(unsigned int i = 0; i < dims_out[0]; ++i)
            {
                tmp[0] = data[i].v[0];
                tmp[1] = data[i].v[1];
                tmp[2] = data[i].v[2];
                pNodes[data[i].id].GetSolutionStepValue(rThisVariable) = tmp;
            }

            /*
             * Release memory
             */
            delete data;
        }
        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            error.printError();
            return -1;
        }
        // catch failure caused by the DataSpace operations
       catch(H5::DataSpaceIException error)
       {
          error.printError();
          return -2;
       }
       // catch failure caused by the DataSpace operations
       catch(H5::DataTypeIException error)
       {
          error.printError();
          return -3;
       }
    }

    int ReadNodalResults_(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        /*
         * Open the dataset
         */
        H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

        /*
         * Read the attribute LEN
         */
        int len;
        H5::Attribute attr = dataset.openAttribute("LEN");
        attr.read(H5::PredType::NATIVE_INT, &len);

        if(len == 3)
        {
            ReadNodalResults_Vector_3(rThisVariable, pModelPart);
        }
        else if(len == 6)
        {
            ReadNodalResults_Vector_6(rThisVariable, pModelPart);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Vector length of this size is not supported:", len)
    }

    int ReadNodalResults_Vector_3(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double v[3];
        } Data_t;

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, v[0]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, v[1]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, v[2]),  H5::PredType::NATIVE_DOUBLE);

        try
        {
            /*
             * Open the dataset
             */
            H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

            /*
             * Get dataspace of the dataset.
             */
            H5::DataSpace dataspace = dataset.getSpace();

            /*
             * Get the number of dimensions in the dataspace.
             */
            int rank = dataspace.getSimpleExtentNdims();
//            KRATOS_WATCH(rank)

            /*
             * Get the dimension size of each dimension in the dataspace and
             * do the bound check.
             */
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

            if(dims_out[0] != pNodes.size())
                KRATOS_THROW_ERROR(std::logic_error, "The dimension of provided data is not compatible:", dims_out[0])

            /*
             * Output buffer
             */
            Data_t* data = new Data_t[pNodes.size()];

            /*
             * Read data from hyperslab in the file into the hyperslab in
             * memory and display the data.
             */
            dataset.read(data, mtype, dataspace, dataspace);

            Vector tmp(3);
            for(unsigned int i = 0; i < dims_out[0]; ++i)
            {
                tmp(0) = data[i].v[0];
                tmp(1) = data[i].v[1];
                tmp(2) = data[i].v[2];
                pNodes[data[i].id].GetSolutionStepValue(rThisVariable) = tmp;
            }

            /*
             * Release memory
             */
            delete data;
        }
        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            error.printError();
            return -1;
        }
        // catch failure caused by the DataSpace operations
       catch(H5::DataSpaceIException error)
       {
          error.printError();
          return -2;
       }
       // catch failure caused by the DataSpace operations
       catch(H5::DataTypeIException error)
       {
          error.printError();
          return -3;
       }
    }

    int ReadNodalResults_Vector_6(
        const Variable<Vector>& rThisVariable,
        ModelPart& r_model_part
    )
    {
        NodesArrayType& pNodes = r_model_part.Nodes();

        typedef struct Data_t {
	        int    id;
	        double v[6];
        } Data_t;

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v1",  HOFFSET(Data_t, v[0]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v2",  HOFFSET(Data_t, v[1]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v3",  HOFFSET(Data_t, v[2]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v4",  HOFFSET(Data_t, v[3]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v5",  HOFFSET(Data_t, v[4]),  H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("v6",  HOFFSET(Data_t, v[5]),  H5::PredType::NATIVE_DOUBLE);

        try
        {
            /*
             * Open the dataset
             */
            H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

            /*
             * Get dataspace of the dataset.
             */
            H5::DataSpace dataspace = dataset.getSpace();

            /*
             * Get the number of dimensions in the dataspace.
             */
            int rank = dataspace.getSimpleExtentNdims();
//            KRATOS_WATCH(rank)

            /*
             * Get the dimension size of each dimension in the dataspace and
             * do the bound check.
             */
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

            if(dims_out[0] != pNodes.size())
                KRATOS_THROW_ERROR(std::logic_error, "The dimension of provided data is not compatible:", dims_out[0])

            /*
             * Output buffer
             */
            Data_t* data = new Data_t[pNodes.size()];

            /*
             * Read data from hyperslab in the file into the hyperslab in
             * memory and display the data.
             */
            dataset.read(data, mtype, dataspace, dataspace);

            Vector tmp(6);
            for(unsigned int i = 0; i < dims_out[0]; ++i)
            {
                tmp(0) = data[i].v[0];
                tmp(1) = data[i].v[1];
                tmp(2) = data[i].v[2];
                tmp(3) = data[i].v[3];
                tmp(4) = data[i].v[4];
                tmp(5) = data[i].v[5];
                pNodes[data[i].id].GetSolutionStepValue(rThisVariable) = tmp;
            }

            /*
             * Release memory
             */
            delete data;
        }
        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            error.printError();
            return -1;
        }
        // catch failure caused by the DataSpace operations
       catch(H5::DataSpaceIException error)
       {
          error.printError();
          return -2;
       }
       // catch failure caused by the DataSpace operations
       catch(H5::DataTypeIException error)
       {
          error.printError();
          return -3;
       }
    }

    /*****************************************************
             FUNCTIONS TO READ ELEMENTAL DATA
    *****************************************************/
    int ReadElementalData_(
        const Variable<bool>& rThisVariable,
        ModelPart& r_model_part,
        bool allow_unequal = true
    )
    {
        ElementsArrayType& pElements = r_model_part.Elements();

        typedef struct Data_t {
	        int id;
	        int v;
        } Data_t;

        /*
         * Create the memory data type.
         */
        H5::CompType mtype(sizeof(Data_t) );
        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
        mtype.insertMember("v",  HOFFSET(Data_t, v),  H5::PredType::NATIVE_INT);

        try
        {
            /*
             * Open the dataset
             */
            H5::DataSet dataset = mpFile->openDataSet(rThisVariable.Name());

            /*
             * Get dataspace of the dataset.
             */
            H5::DataSpace dataspace = dataset.getSpace();

            /*
             * Get the number of dimensions in the dataspace.
             */
            int rank = dataspace.getSimpleExtentNdims();
//            KRATOS_WATCH(rank)

            /*
             * Get the dimension size of each dimension in the dataspace and
             * do the bound check.
             */
            hsize_t dims_out[rank];
            int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

            if(!allow_unequal)
                if(dims_out[0] != pElements.size())
                {
                    std::stringstream ss;
                    ss << "The dimension of provided data is not compatible.";
                    ss << " The number of elements is " << pElements.size();
                    ss << ", the dimension of read-out data is " << dims_out[0];
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }

            /*
             * Output buffer
             */
            Data_t* data = new Data_t[dims_out[0]];

            /*
             * Read data from hyperslab in the file into the hyperslab in
             * memory and display the data.
             */
            dataset.read(data, mtype, dataspace, dataspace);

            double tmp;
            for(unsigned int i = 0; i < dims_out[0]; ++i)
            {
                tmp = data[i].v;
                if(tmp == 0)
                    pElements[data[i].id].SetValue(rThisVariable, false);
                else
                    pElements[data[i].id].SetValue(rThisVariable, true);
            }

            /*
             * Release memory
             */
            delete data;
        }
        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            error.printError();
            return -1;
        }
        // catch failure caused by the DataSpace operations
       catch(H5::DataSpaceIException error)
       {
          error.printError();
          return -2;
       }
       // catch failure caused by the DataSpace operations
       catch(H5::DataTypeIException error)
       {
          error.printError();
          return -3;
       }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    HDF5PostUtility& operator=(HDF5PostUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    HDF5PostUtility(HDF5PostUtility const& rOther)
    {
    }

    ///@}

}; // Class HDF5PostUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, HDF5PostUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HDF5PostUtility& rThis)
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


//      void WriteNodalResults_(
//        const Variable<Vector>& rThisVariable,
//        ModelPart& r_model_part
//    )
//    {
//        NodesArrayType& pNodes = r_model_part.Nodes();
//
//        int len = (*(pNodes.ptr_begin()))->GetSolutionStepValue(rThisVariable).size();
//
//        typedef struct Data_t {
//	        int    id;
//	        double *v;
//        } Data_t;
//
//        typedef struct Data_buffer_t {
//	        int    id;
//	        hvl_t  v;
//        } Data_buffer_t;
//
//        Data_t* data = new Data_t[pNodes.size()];
//        Data_buffer_t* data_buffer = new Data_buffer_t[pNodes.size()];
//
//        /*
//         * Initialize the data
//         */
//        int cnt = 0;
//        for(typename NodesArrayType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
//        {
//            data[cnt].id = (*it)->Id();
//            data[cnt].v = new double[len];
//            Vector& v = (*it)->GetSolutionStepValue(rThisVariable);
//            for(int i = 0; i < len; ++i)
//                data[cnt].v[i] = v[i];
//
//            data_buffer[cnt].id = data[cnt].id;
//            data_buffer[cnt].v.len = len;
//            data_buffer[cnt].v.p = data[cnt].v;
//
//            ++cnt;
//        }
//
//        /*
//         * Create the memory data type.
//         */
//        H5::VarLenType vlen_tid(&H5::PredType::NATIVE_DOUBLE);
//        H5::CompType mtype(sizeof(Data_t) );
//        mtype.insertMember("id", HOFFSET(Data_t, id), H5::PredType::NATIVE_INT);
//        mtype.insertMember("v", HOFFSET(Data_t, v), vlen_tid);
//
//        /*
//         * Create the data space.
//         */
//        hsize_t dim[] = {pNodes.size()};   /* Dataspace dimensions */
//        H5::DataSpace space(1, dim);
//
//        /*
//         * Create the dataset.
//         */
//        H5::DataSet* dataset;
//        dataset = new H5::DataSet(mpFile->createDataSet(rThisVariable.Name(), mtype, space));
//
//        /*
//         * Write data to the dataset;
//         */
//        dataset->write(data_buffer, mtype);
//
//        /*
//         * Release resources
//         */
//        delete dataset;
//        delete data;
//        delete data_buffer;
//    }
//
