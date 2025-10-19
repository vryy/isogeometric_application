// see isogeometric_application/LICENCE.txt

//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Jan 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_PARTITIONING_PROCESS_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_PARTITIONING_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

// External includes
#include <parmetis.h>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/graph_coloring_process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/mpi_communicator.h"

extern "C"
{
    int METIS_PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize,
                           idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval,
                           idx_t *epart, idx_t *npart);
};

namespace Kratos
{

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

/** Detail class definition.
 */

class IsogeometricPartitioningProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsogeometricPartitioningProcess
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricPartitioningProcess);

    typedef std::size_t size_type;
    typedef std::size_t index_type;
    typedef boost::numeric::ublas::matrix<int> graph_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricPartitioningProcess(ModelPart& rModelPart, IO& rIO, size_type NumberOfPartitions)
        : mrModelPart(rModelPart), mrIO(rIO), mNumberOfPartitions(NumberOfPartitions)
    {
        KRATOS_TRY
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::stringstream log_filename;
        log_filename << "kratos_metis_" << rank << ".log";
        mLogFile.open(log_filename.str().c_str());
        KRATOS_CATCH("")
    }

    /// Copy constructor.
    IsogeometricPartitioningProcess(IsogeometricPartitioningProcess const& rOther)
        : mrModelPart(rOther.mrModelPart), mrIO(rOther.mrIO), mNumberOfPartitions(rOther.mNumberOfPartitions)
    {
        KRATOS_TRY
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::stringstream log_filename;
        log_filename << "kratos_metis_" << rank << ".log";
        mLogFile.open(log_filename.str().c_str());
        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~IsogeometricPartitioningProcess()
    {
        mLogFile.close();
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        KRATOS_TRY;

        int rank, number_of_processes;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

        if (number_of_processes < 2) // There is no need to partition it and just reading the input
        {
            mLogFile << "No partitioning" << std::endl;
            mrIO.ReadModelPart(mrModelPart);
            return;
        }

        // Set MPICommunicator as modelpart's communicator
        VariablesList* mVariables_List = &mrModelPart.GetNodalSolutionStepVariablesList();
        mrModelPart.SetCommunicator(Communicator::Pointer(new MPICommunicator(mVariables_List)));

        // if mNumberOfPartitions is not defined we set it to the number_of_processes
        if (mNumberOfPartitions == 0)
        {
            mNumberOfPartitions = static_cast<size_type>(number_of_processes);
        }

        if (rank == 0)
        {
            KRATOS_WATCH(mNumberOfPartitions);
        }

        // Read connectivities
        IO::ConnectivitiesContainerType elements_connectivities;
        int number_of_elements = mrIO.ReadElementsConnectivities(elements_connectivities);

        // Read nodes
        ModelPart::NodesContainerType temp_nodes;
        mrIO.ReadNodes(temp_nodes);
        int number_of_nodes = temp_nodes.size(); // considering sequential numbering!!

        mLogFile << rank << " : Reading nodes completed, number_of_nodes = " << number_of_nodes << ", number_of_elements = " << number_of_elements << std::endl;

        idx_t* epart = new idx_t[number_of_elements];
        idx_t* npart = new idx_t[number_of_nodes];

        graph_type domains_graph = boost::numeric::ublas::zero_matrix<int>(mNumberOfPartitions, mNumberOfPartitions);
        graph_type domains_colored_graph;
        int* coloring_send_buffer = NULL;

        // Adding interface meshes
        mrModelPart.GetMeshes().push_back(ModelPart::MeshType());

        int colors_number;
        if (rank == 0)
        {
            CallingMetis(number_of_nodes, number_of_elements, elements_connectivities, npart, epart);
            CalculateDomainsGraph(domains_graph, number_of_elements, elements_connectivities, npart, epart);
            GraphColoringProcess(mNumberOfPartitions, domains_graph, domains_colored_graph, colors_number).Execute();
            KRATOS_WATCH(colors_number);
            KRATOS_WATCH(domains_colored_graph);

            // Filling the sending buffer
            int buffer_index = 0;
            coloring_send_buffer = new int[mNumberOfPartitions * colors_number];
            for (unsigned int i = 0; i < mNumberOfPartitions; ++i)
                for (int j = 0; j < colors_number; j++)
                {
                    coloring_send_buffer[buffer_index++] = domains_colored_graph(i, j);
                }

            mLogFile << rank << " : colors_number = " << colors_number << std::endl;
            mLogFile << rank << " : coloring_send_buffer = [";
            for (size_type j = 0; j < mNumberOfPartitions * colors_number; ++j)
            {
                mLogFile << coloring_send_buffer[j] << " ,";
            }
            mLogFile << "]" << std::endl;
        }

        // Broadcasting partioning information
        MPI_Bcast(&colors_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(npart, number_of_nodes, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(epart, number_of_elements, MPI_INT, 0, MPI_COMM_WORLD);

        Communicator::NeighbourIndicesContainerType& neighbours_indices = mrModelPart.GetCommunicator().NeighbourIndices();
        if (neighbours_indices.size() != static_cast<unsigned int>(colors_number))
        {
            neighbours_indices.resize(colors_number, false);
        }
        for (int i = 0; i < colors_number; ++i)
        {
            neighbours_indices[i] = 0;
        }

        mLogFile << rank << " : neighbours_indices = [";
        for (unsigned int j = 0; j < neighbours_indices.size(); ++j)
        {
            mLogFile << neighbours_indices[j] << " ,";
        }
        mLogFile << "]" << std::endl;

        mLogFile << rank << " : colors_number = " << colors_number << std::endl;
        mrModelPart.GetCommunicator().SetNumberOfColors(colors_number);

        MPI_Scatter(coloring_send_buffer, colors_number, MPI_INT, &(neighbours_indices[0]), colors_number, MPI_INT, 0, MPI_COMM_WORLD);

        mLogFile << rank << " : [";
        for (int j = 0; j < colors_number; j++)
        {
            mLogFile << mrModelPart.GetCommunicator().NeighbourIndices()[j] << " ,";
        }
        mLogFile << "]" << std::endl;

        // Adding local, ghost and interface meshes to ModelPart if is necessary
        int number_of_meshes = ModelPart::Kratos_Ownership_Size + colors_number; // (all + local + ghost) + (colors_number for interfaces)
        if (mrModelPart.GetMeshes().size() < static_cast<unsigned int>(number_of_meshes))
            for (int i = mrModelPart.GetMeshes().size(); i < number_of_meshes; ++i)
            {
                mrModelPart.GetMeshes().push_back(ModelPart::MeshType());
            }

        for (ModelPart::NodeIterator i_node = temp_nodes.begin(); i_node != temp_nodes.end(); ++i_node)
        {
            i_node->SetSolutionStepVariablesList(&(mrModelPart.GetNodalSolutionStepVariablesList()));
        }

        // Adding nodes to model_part
        AddingNodes(temp_nodes, number_of_elements, elements_connectivities, npart, epart);

        // Adding properties to modelpart
        mLogFile << rank << " : Start reading Properties " << std::endl;
        mrIO.ReadProperties(mrModelPart.rProperties());
        mLogFile << rank << " : End adding Properties " << std::endl;

        // Adding elements to each partition mesh
        AddingElements(temp_nodes, npart, epart);

        // Adding conditions to each partition mesh
        ModelPart::ConditionsContainerType temp_conditions;
        AddingConditions(temp_nodes, npart, epart, temp_conditions);

        mrIO.ReadInitialValues(temp_nodes, mrModelPart.Elements(), temp_conditions);

        mLogFile << rank << " : start cleaning memory " << std::endl;
        delete [] epart;
        delete [] npart;
        if (rank == 0)
        {
            mLogFile << rank << " : deleting coloring_send_buffer " << std::endl;
            delete[] coloring_send_buffer;
        }
        mLogFile << rank << " : cleaning memory Finished" << std::endl;

        KRATOS_CATCH("")
    }

    void CalculateDomainsGraph(graph_type& rDomainsGraph, size_type NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idx_t* NPart, idx_t* EPart)
    {
        for (size_type i_element = 0; i_element < NumberOfElements; ++i_element)
            for (std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin();
                    i_node != ElementsConnectivities[i_element].end(); ++i_node)
            {
                size_type node_rank = NPart[*i_node - 1];
                size_type element_rank = EPart[i_element];
                if (node_rank != element_rank)
                {
                    rDomainsGraph(node_rank, element_rank) = 1;
                    rDomainsGraph(element_rank, node_rank) = 1;
                }
            }
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
        return "IsogeometricPartitioningProcess";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricPartitioningProcess";
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

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

    ModelPart& mrModelPart;

    IO& mrIO;

    size_type mNumberOfPartitions;

    std::ofstream mLogFile;

    size_type mDimension;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void CallingMetis(size_type NumberOfNodes, size_type NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idx_t* NPart, idx_t* EPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        mLogFile << rank << ": Calling Metis" << std::endl;

        idx_t ne = NumberOfElements;
        idx_t nn = NumberOfNodes;

        // fill up eptr and calculate the size of total connectivities
        idx_t* eptr = new idx_t[ne + 1];
        size_type connectivity_size = 0;
        index_type cnt = 0;
        eptr[0] = 0;
        for (IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin();
                i_connectivities != ElementsConnectivities.end(); ++i_connectivities)
        {
            connectivity_size += i_connectivities->size();
            eptr[++cnt] = connectivity_size;
        }

        // fill up eind
        idx_t* eind = new idx_t[connectivity_size];
        cnt = 0;
        for (IO::ConnectivitiesContainerType::iterator i_connectivities = ElementsConnectivities.begin();
                i_connectivities != ElementsConnectivities.end(); ++i_connectivities)
        {
            for (unsigned int j = 0; j < i_connectivities->size(); ++j)
            {
                eind[cnt++] = (*i_connectivities)[j] - 1;    //transform to zero-based indexing
            }
        }

//        mLogFile << rank << ": eptr:";
//        for(index_type i = 0; i < ne + 1; ++i)
//            mLogFile << " " << eptr[i];
//        mLogFile << std::endl;
//
//        mLogFile << rank << ": eind:";
//        for(index_type i = 0; i < ne; ++i)
//        {
//            for(index_type j = eptr[i]; j < eptr[i + 1]; ++j)
//                mLogFile << " " << eind[j];
//            mLogFile << std::endl;
//        }
//        mLogFile << std::endl;

        mLogFile << rank << ": Calling METIS_PartMeshDual" << std::endl;
        int status = 0;
        idx_t edgecut;
        idx_t ncommon = 4;
        idx_t nparts = mNumberOfPartitions;
        status = METIS_PartMeshDual(&ne, &nn, eptr, eind,
                                    NULL, NULL, &ncommon, &nparts,
                                    NULL, NULL, &edgecut, EPart, NPart);

        // release memory
        delete [] eptr;
        delete [] eind;
    }

    void AddingNodes(ModelPart::NodesContainerType& AllNodes, size_type NumberOfElements, IO::ConnectivitiesContainerType& ElementsConnectivities, idx_t* NPart, idx_t* EPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        Communicator& r_communicator = mrModelPart.GetCommunicator();

        // first adding the partition's nodes
        mLogFile << rank << " : Adding nodes to modelpart" << std::endl;
        for (ModelPart::NodeIterator i_node = AllNodes.begin(); i_node != AllNodes.end(); ++i_node)
            if (NPart[i_node->Id() - 1] == rank)
            {
                mrModelPart.AssignNode(*(i_node.base()));
                r_communicator.LocalMesh().Nodes().push_back(*(i_node.base()));
                i_node->GetSolutionStepValue(PARTITION_INDEX) = rank;
            }

        std::vector<int> interface_indices(mNumberOfPartitions, -1);
        vector<int>& neighbours_indices = r_communicator.NeighbourIndices();
        for (size_type i = 0; i < neighbours_indices.size(); i++)
            if (size_type(neighbours_indices[i]) < interface_indices.size())
            {
                interface_indices[neighbours_indices[i]] = i;
            }

        // now adding interface nodes which belongs to other partitions
        mLogFile << rank << " : Adding interface nodes to modelpart" << std::endl;
        for (size_type i_element = 0; i_element < NumberOfElements; ++i_element)
            if (EPart[i_element] == rank)
            {
                for (std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin();
                        i_node != ElementsConnectivities[i_element].end(); ++i_node)
                {
                    int node_partition = NPart[*i_node - 1];
                    if (node_partition != rank)
                    {

                        ModelPart::NodeType::Pointer p_node = AllNodes(*i_node);
                        // Giving model part's variables list to the node
                        p_node->SetSolutionStepVariablesList(&(mrModelPart.GetNodalSolutionStepVariablesList()));

                        //set buffer size
                        p_node->SetBufferSize(mrModelPart.GetBufferSize());
                        mrModelPart.Nodes().push_back(p_node);
                        r_communicator.GhostMesh().Nodes().push_back(p_node);
                        if (size_type(interface_indices[node_partition]) < neighbours_indices.size())
                        {
                            r_communicator.GhostMesh(interface_indices[node_partition]).Nodes().push_back(p_node);
                            r_communicator.InterfaceMesh(interface_indices[node_partition]).Nodes().push_back(p_node);
                        }
                        else
                        {
                            std::cout << rank << " : ERROR! Node #" << *i_node << " has not registered interface for partition #" << node_partition << std::endl;
                        }
                        r_communicator.InterfaceMesh().Nodes().push_back(p_node);
                        p_node->GetSolutionStepValue(PARTITION_INDEX) = NPart[*i_node - 1];
                    }
                }
            }
            else // adding the owned interface nodes
            {
                for (std::vector<std::size_t>::iterator i_node = ElementsConnectivities[i_element].begin();
                        i_node != ElementsConnectivities[i_element].end(); ++i_node)
                    if (NPart[*i_node - 1] == rank)
                    {
                        size_type mesh_index = interface_indices[EPart[i_element]];
                        if (mesh_index > neighbours_indices.size()) // Means the neighbour domain is not registered!!
                        {
                            //std::cout << rank << " : cannot find interface for element #" << i_element << " with rank " << EPart[i_element] << std::endl;
                            KRATOS_ERROR << "Cannot find the neighbour domain " << EPart[i_element];
                        }

                        ModelPart::NodeType::Pointer p_node = AllNodes(*i_node);
                        r_communicator.LocalMesh().Nodes().push_back(p_node);
                        r_communicator.LocalMesh(mesh_index).Nodes().push_back(p_node);
                        r_communicator.InterfaceMesh(mesh_index).Nodes().push_back(p_node);
                        r_communicator.InterfaceMesh().Nodes().push_back(p_node);
                    }
            }

        // After making push_back to the nodes list now we need to make unique and sort for all meshes in communicator
        mrModelPart.Nodes().Unique();
        r_communicator.LocalMesh().Nodes().Unique();
        r_communicator.GhostMesh().Nodes().Unique();
        r_communicator.InterfaceMesh().Nodes().Unique();
        for (size_type i = 0; i < r_communicator.LocalMeshes().size(); i++)
        {
            r_communicator.LocalMesh(i).Nodes().Unique();
        }
        for (size_type i = 0; i < r_communicator.GhostMeshes().size(); i++)
        {
            r_communicator.GhostMesh(i).Nodes().Unique();
        }
        for (size_type i = 0; i < r_communicator.InterfaceMeshes().size(); i++)
        {
            r_communicator.InterfaceMesh(i).Nodes().Unique();
        }

        mLogFile << rank << " : Nodes added to modelpart" << std::endl;

    }

    void AddingElements(ModelPart::NodesContainerType& AllNodes, idx_t* NPart, idx_t* EPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        mLogFile << rank << " : Reading elements" << std::endl;
        ModelPart::ElementsContainerType temp_elements;
        mrIO.ReadElements(AllNodes, mrModelPart.rProperties(), temp_elements);

        mLogFile << rank << " : Adding elements to modelpart" << std::endl;
        idx_t* epart_position = EPart;
        for (Kratos::ModelPart::ElementIterator i_element = temp_elements.begin();
                i_element != temp_elements.end(); i_element++)
        {
            if (*epart_position == rank)
            {
                mrModelPart.AddElement(*(i_element.base()));
                mrModelPart.GetCommunicator().LocalMesh().AddElement(*(i_element.base()));
                //                 mrModelPart.AddElement(*(i_element.base()),  ModelPart::Kratos_Local);
            }
            epart_position++;
        }
        mLogFile << rank << " : Elements added" << std::endl;
    }

    virtual void AddingConditions(ModelPart::NodesContainerType& AllNodes, idx_t* NPart, idx_t* EPart, ModelPart::ConditionsContainerType& AllConditions)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        mLogFile << rank << " : Reading conditions" << std::endl;
        mrIO.ReadConditions(AllNodes, mrModelPart.rProperties(), AllConditions);

        mLogFile << rank << " : Adding conditions to modelpart" << std::endl;
        for (Kratos::ModelPart::ConditionIterator i_condition = AllConditions.begin();
                i_condition != AllConditions.end(); i_condition++)
        {
            bool is_local = 1;
            // See if all of the condition nodes are in this partition as a local or even as a ghost
            // TODO: THIS IS DANGEROUSE AND MAY FAILE DUE TO THE MESH!!! MUST BE CHANGED!!
            for (ModelPart::ConditionType::GeometryType::iterator i_node = i_condition->GetGeometry().begin();
                    i_node != i_condition->GetGeometry().end(); ++i_node)
                if (mrModelPart.Nodes().find(i_node->Id()) == mrModelPart.Nodes().end())
                {
                    is_local = 0;
                }
            if (is_local)
            {
                mrModelPart.AddCondition(*(i_condition.base()));
                mrModelPart.GetCommunicator().LocalMesh().AddCondition(*(i_condition.base()));
            }
        }
        mLogFile << rank << " : Conditions added" << std::endl;
    }

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
    IsogeometricPartitioningProcess & operator=(IsogeometricPartitioningProcess const& rOther);

    /// Copy constructor.
    //IsogeometricPartitioningProcess(IsogeometricPartitioningProcess const& rOther);

    ///@}

}; // Class IsogeometricPartitioningProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  IsogeometricPartitioningProcess& rThis)
{
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const IsogeometricPartitioningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_PARTITIONING_PROCESS_INCLUDED defined 
