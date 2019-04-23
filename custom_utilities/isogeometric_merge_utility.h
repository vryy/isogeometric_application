//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 29 Jul 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_MERGE_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_MERGE_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#include "containers/model.h"
#include "containers/variables_list.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#include "isogeometric_application.h"


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

template<class PointType>
class LoosePointType
{
public:
    LoosePointType(PointType& rP) : mrP(rP) {}
    ~LoosePointType() {}
    bool operator<(const LoosePointType& rOther) const
    {
        if(mrP.X0() < rOther.mrP.X0() - 1.0e-6)
            return true;
        if(mrP.X0() > rOther.mrP.X0() + 1.0e-6)
            return false;
        if(mrP.Y0() < rOther.mrP.Y0() - 1.0e-6)
            return true;
        if(mrP.Y0() > rOther.mrP.Y0() + 1.0e-6)
            return false;
        if(mrP.Z0() < rOther.mrP.Z0() - 1.0e-6)
            return true;
        if(mrP.Z0() > rOther.mrP.Z0() + 1.0e-6)
            return false;
        return false;
    }
    PointType& GetPoint() const {return mrP;}
//    PointType GetPoint() {return mrP;}
private:
    PointType& mrP;
};


/// Short class definition.
/**
 * A utility to combine multiple model_part into one model_part and skipping the coincident nodes. It is useful to combine multiple FEM model_part arised from Bezier mesh in multipatch isogeometric analysis.
 */
class IsogeometricMergeUtility
{
public:
    ///@name Type Definitions
    ///@{

    typedef boost::numeric::ublas::vector<double> ValuesContainerType;

    typedef boost::numeric::ublas::matrix<double> ValuesArrayContainerType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef Element::GeometryType GeometryType;

    typedef GeometryType::PointType PointType;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename Node<3>::DofsContainerType DofsContainerType;

    typedef ModelPart::PropertiesIterator PropertiesIterator;

    typedef Node<3> NodeType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SerialSparseSpaceType;

    typedef UblasSpace<double, Matrix, Vector> SerialDenseSpaceType;

    typedef LinearSolver<SerialSparseSpaceType, SerialDenseSpaceType> LinearSolverType;

    typedef std::size_t IndexType;


    /// Pointer definition of IsogeometricMergeUtility
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(IsogeometricMergeUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsogeometricMergeUtility()
    {
    }

    /// Destructor.
    virtual ~IsogeometricMergeUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Add(Model::Pointer pModel, const std::string& ModelPartName)
    {
        mpModelContainer.push_back(pModel);
        mModelPartNames.push_back(ModelPartName);
    }

    void Export(ModelPart& r_model_part)
    {
        // firstly iterate through all the nodes in all model part to extract a non-repeated list of nodes
        typedef std::pair<int, int> KeyType;
        typedef std::set<LoosePointType<PointType> > SetType;
        typedef std::map<KeyType, SetType::iterator> MapType;
        SetType PointSet;
        MapType PointMap;

//        std::set<PointType> PointSetDummy;
//      // note: this does not work since node does not implement operator<
//        for(int i = 0; i < mpModelContainer.size(); ++i)
//        {
//            NodesContainerType ThisNodes = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).Nodes();
//            for(NodesContainerType::iterator it = ThisNodes.begin(); it != ThisNodes.end(); ++it)
//            {
//                PointSetDummy.insert(*it);
//            }
//        }
//        KRATOS_WATCH(PointSetDummy.size())

        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            NodesContainerType ThisNodes = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).Nodes();
            for(NodesContainerType::iterator it = ThisNodes.begin(); it != ThisNodes.end(); ++it)
            {
                LoosePointType<PointType> P(*it);
                std::pair<SetType::iterator, bool> iP = PointSet.insert(P);
                KeyType K(i, it->Id());
                PointMap[K] = iP.first;
            }
        }
        KRATOS_WATCH(PointSet.size())

        // create a unified variables_list
        VariablesList ThisVariablesList;
        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            VariablesList& tmp = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).GetNodalSolutionStepVariablesList();
            for(VariablesList::ptr_const_iterator it = tmp.ptr_begin(); it != tmp.ptr_end(); ++it)
                ThisVariablesList.Add(*(*it));
        }

        // add the unified variables_list to the new model_part
        for(VariablesList::ptr_const_iterator it = ThisVariablesList.ptr_begin(); it != ThisVariablesList.ptr_end(); ++it)
            r_model_part.GetNodalSolutionStepVariablesList().Add(*(*it));
        KRATOS_WATCH(r_model_part.GetNodalSolutionStepVariablesList())

        // create a maximum buffer
        int buffer_size = 0;
        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            int tmp_buff_size = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).GetBufferSize();
            if(tmp_buff_size > buffer_size)
                buffer_size = tmp_buff_size;
        }
        r_model_part.SetBufferSize(buffer_size);
        KRATOS_WATCH(r_model_part.GetBufferSize())

        // add nodes to the new model part
        int lastNode = 0;
        std::map<LoosePointType<PointType>, int> NodeIds; // map from reduced node set to new node set id
        for(SetType::iterator it = PointSet.begin(); it != PointSet.end(); ++it)
        {
            NodeType::Pointer pNewNode( new Node<3>( 0, it->GetPoint() ) );;
            pNewNode->SetId(++lastNode);
            NodeIds[*it] = lastNode;
            pNewNode->SetSolutionStepVariablesList(&(r_model_part.GetNodalSolutionStepVariablesList())); // to make sure it synchronized with model_part variables list
            pNewNode->SetBufferSize(r_model_part.GetBufferSize());
//            r_model_part.AddNode(NewNode); // for some reason it created segmentation fault error
            r_model_part.Nodes().push_back(pNewNode);
//            // TODO: transfer nodal data

//            NodeType temp_node;
//            temp_node.SetSolutionStepVariablesList(&r_model_part.GetNodalSolutionStepVariablesList());
//            temp_node.SetBufferSize(r_model_part.GetBufferSize());
//            temp_node.SetId(++lastNode);
//            NodeIds[*it] = lastNode;
//            temp_node.X() = it->GetPoint().X();
//            temp_node.Y() = it->GetPoint().Y();
//            temp_node.Z() = it->GetPoint().Z();
//            temp_node.X0() = temp_node.X();
//            temp_node.Y0() = temp_node.Y();
//            temp_node.Z0() = temp_node.Z();
//            r_model_part.Nodes().push_back(temp_node);
        }

        // add elements to the new model part
        int lastElement = 0;
        std::string NodeKey = std::string("Node");
        typedef typename KratosComponents<Element>::ComponentsContainerType ElementComponentsContainerType;
        ElementsContainerType NewElements;
        ElementComponentsContainerType ElementComponents = KratosComponents<Element>::GetComponents();
        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            ElementsContainerType pElements = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).Elements();
            for(ElementsContainerType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
            {
                std::string element_name;
                for(ElementComponentsContainerType::iterator cit = ElementComponents.begin(); cit != ElementComponents.end(); ++cit)
                    if(typeid(*(cit->second)).name() == typeid(*(*it)).name())
                        if(typeid(cit->second->GetGeometry()).name() == typeid((*it)->GetGeometry()).name())
                            element_name = cit->first;
//                KRATOS_WATCH(element_name)
                Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
                Properties::Pointer p_temp_properties = (*it)->pGetProperties();
                Element::NodesArrayType temp_element_nodes;
                for(int j = 0; j < (*it)->GetGeometry().size(); ++j)
                {
                    KeyType K(i, (*it)->GetGeometry()[j].Id());
                    SetType::iterator iP = PointMap[K];
                    temp_element_nodes.push_back(*(FindKey(r_model_part.Nodes(), NodeIds[*iP], NodeKey).base()));
                }
                Element::Pointer elem = r_clone_element.Create(++lastElement, temp_element_nodes, p_temp_properties);
                elem->Data() = (*it)->Data(); // transfer elemental data
                NewElements.push_back(elem);
            }
        }
        // add new elements
        for( ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it )
            r_model_part.Elements().push_back( *it );
        NewElements.clear();

        // add conditions to the new model part
        int lastCondition = 0;
        typedef typename KratosComponents<Condition>::ComponentsContainerType ConditionComponentsContainerType;
        ConditionsContainerType NewConditions;
        ConditionComponentsContainerType ConditionComponents = KratosComponents<Condition>::GetComponents();
        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            ConditionsContainerType pConditions = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).Conditions();
            for(ConditionsContainerType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
            {
                std::string condition_name;
                for(ConditionComponentsContainerType::iterator cit = ConditionComponents.begin(); cit != ConditionComponents.end(); ++cit)
                    if(typeid(*(cit->second)).name() == typeid(*(*it)).name())
                        if(typeid(cit->second->GetGeometry()).name() == typeid((*it)->GetGeometry()).name())
                            condition_name = cit->first;
//                KRATOS_WATCH(condition_name)
                Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
                Properties::Pointer p_temp_properties = (*it)->pGetProperties();
                Condition::NodesArrayType temp_condition_nodes;
                for(int j = 0; j < (*it)->GetGeometry().size(); ++j)
                {
                    KeyType K(i, (*it)->GetGeometry()[j].Id());
                    SetType::iterator iP = PointMap[K];
                    temp_condition_nodes.push_back(*(FindKey(r_model_part.Nodes(), NodeIds[*iP], NodeKey).base()));
                }
                Condition::Pointer cond = r_clone_condition.Create(++lastCondition, temp_condition_nodes, p_temp_properties);
                cond->Data() = (*it)->Data(); // transfer conditional data
                NewConditions.push_back(cond);
            }
        }
        // add new conditions
        for( ConditionsContainerType::ptr_iterator it = NewConditions.ptr_begin(); it != NewConditions.ptr_end(); ++it )
            r_model_part.Conditions().push_back( *it );
        NewConditions.clear();

        // add properties to the model_part
        int lastProperties = 0;
        for(int i = 0; i < mpModelContainer.size(); ++i)
        {
            for(PropertiesIterator it = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).PropertiesBegin(); it != mpModelContainer[i]->GetModelPart(mModelPartNames[i]).PropertiesEnd(); ++it)
            {
                Properties::Pointer tmp = mpModelContainer[i]->GetModelPart(mModelPartNames[i]).pGetProperties(it->Id());
                tmp->SetId(++lastProperties);
                r_model_part.AddProperties(tmp);
            }
        }
    }


    ///@}
    ///@name Access
    ///@{

    void DumpNodalVariablesList(ModelPart& r_model_part) const
    {
        NodesContainerType pNodes = r_model_part.Nodes();
        for(NodesContainerType::ptr_iterator it = pNodes.ptr_begin(); it != pNodes.ptr_end(); ++it)
        {
            KRATOS_WATCH((*it)->SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_WATCH((*it)->SolutionStepsDataHas(DISPLACEMENT_X))
            KRATOS_WATCH((*it)->GetSolutionStepValue(DISPLACEMENT_X))
        }
    }


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
        buffer << "IsogeometricMergeUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IsogeometricMergeUtility";
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
    std::vector<Model::Pointer> mpModelContainer;
    std::vector<std::string> mModelPartNames;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    //**********AUXILIARY FUNCTION**************************************************************
    //******************************************************************************************
    template<class TContainerType, class TKeyType>
    typename TContainerType::iterator FindKey(TContainerType& ThisContainer, TKeyType& ThisKey, std::string& ComponentName)
    {
        typename TContainerType::iterator i_result;
        if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            std::stringstream buffer;
            buffer << ComponentName << " #" << ThisKey << " is not found.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }


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
    IsogeometricMergeUtility& operator=(IsogeometricMergeUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    IsogeometricMergeUtility(IsogeometricMergeUtility const& rOther)
    {
    }

    ///@}

}; // Class IsogeometricMergeUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, IsogeometricMergeUtility& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const IsogeometricMergeUtility& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#undef DEBUG_LEVEL1
#undef DEBUG_LEVEL2
#undef DEBUG_MULTISOLVE
#undef DEBUG_GENERATE_MESH
#undef ENABLE_PROFILING

#endif // ISOGEOMETRIC_MERGE_UTILITY

