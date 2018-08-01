/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2015 $
//   Revision:            $Revision: 1.0 $
//
//


// Project includes
#include "bezier_model_part_io.h"
#include "custom_geometries/isogeometric_geometry.h"
#include "custom_geometries/geo_1d_bezier.h"
#include "custom_geometries/geo_2d_bezier.h"
#include "custom_geometries/geo_2d_bezier_3.h"
#include "custom_geometries/geo_3d_bezier.h"
#include "custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application.h"

namespace Kratos
{

    BezierModelPartIO::BezierModelPartIO(std::string const& Filename, const Flags Options)
    : ModelPartIO(Filename, Options)
    , mpBezierInfoContainer(new BezierInfoContainerType())
    {}

    BezierModelPartIO::~BezierModelPartIO()
    {}

    void BezierModelPartIO::ReadModelPart(ModelPart & rThisModelPart)
    {
        KRATOS_TRY

        Timer::Start("Reading Input");

        ModelPartIO::ResetInput();
        std::string word;
        while(true)
        {
            ModelPartIO::ReadWord(word);
            #if defined(KRATOS_SD_REF_NUMBER_2)
            if(mFile.eof())
            #elif defined(KRATOS_SD_REF_NUMBER_3)
            if(mpStream->eof())
            #endif
                break;
            ModelPartIO::ReadBlockName(word);
            if(word == "ModelPartData")
                ModelPartIO::ReadModelPartDataBlock(rThisModelPart);
            else if(word == "Table")
                ModelPartIO::ReadTableBlock(rThisModelPart.Tables());
            else if(word == "Properties")
                ModelPartIO::ReadPropertiesBlock(rThisModelPart.rProperties());
            else if(word == "Nodes")
                ModelPartIO::ReadNodesBlock(rThisModelPart);
            else if(word == "Elements")
                ModelPartIO::ReadElementsBlock(rThisModelPart);
            else if(word == "Conditions")
                ModelPartIO::ReadConditionsBlock(rThisModelPart);
            else if(word == "BezierBlock")
                this->ReadBezierBlock(rThisModelPart);
            else if(word == "NodalData")
                ModelPartIO::ReadNodalDataBlock(rThisModelPart);
            else if(word == "ElementalData")
                ModelPartIO::ReadElementalDataBlock(rThisModelPart.Elements());
            else if(word == "ConditionalData")
                ModelPartIO::ReadConditionalDataBlock(rThisModelPart.Conditions());
            else if(word == "CommunicatorData")
            {
                ModelPartIO::ReadCommunicatorDataBlock(rThisModelPart.GetCommunicator(), rThisModelPart.Nodes());
                //Adding the elements and conditions to the communicator
                rThisModelPart.GetCommunicator().LocalMesh().Elements() = rThisModelPart.Elements();
                rThisModelPart.GetCommunicator().LocalMesh().Conditions() = rThisModelPart.Conditions();
            }
            else if(word == "Mesh")
                ModelPartIO::ReadMeshBlock(rThisModelPart);

        }
        std::cout << "  [Total Lines Read : " << mNumberOfLines<<"]";
        std::cout << std::endl;

	    Timer::Stop("Reading Input");

        KRATOS_CATCH("")
    }

    void BezierModelPartIO::ReadBezierBlock(ModelPart & rThisModelPart)
    {
        std::string word;

        #if defined(KRATOS_SD_REF_NUMBER_2)
        while(!mFile.eof())
        #elif defined(KRATOS_SD_REF_NUMBER_3)
        while(!mpStream->eof())
        #endif
        {
            ModelPartIO::ReadWord(word);

            if(ModelPartIO::CheckEndBlock("BezierBlock", word))
                break;

            #if defined(KRATOS_SD_REF_NUMBER_2)
            if(mFile.eof())
            #elif defined(KRATOS_SD_REF_NUMBER_3)
            if(mpStream->eof())
            #endif
                break;

            ModelPartIO::ReadBlockName(word);
            if(word == "IsogeometricBezierData")
                this->ReadIsogeometricBezierDataBlock(*mpBezierInfoContainer);
            else if(word == "ElementsWithGeometry")
                this->ReadElementsWithGeometryBlock(rThisModelPart.Nodes(), rThisModelPart.rProperties(), *mpBezierInfoContainer, rThisModelPart.Elements());
            else if(word == "ConditionsWithGeometry")
                this->ReadConditionsWithGeometryBlock(rThisModelPart.Nodes(), rThisModelPart.rProperties(), *mpBezierInfoContainer, rThisModelPart.Conditions());
        }
    }

    void BezierModelPartIO::ReadIsogeometricBezierDataBlock(BezierInfoContainerType& rThisBezierInfo)
    {
        std::string word;
        SizeType number_of_geometries_read = 0;

        std::cout << "  [Reading Bezier Geometries : ";

        #if defined(KRATOS_SD_REF_NUMBER_2)
        while(!mFile.eof())
        #elif defined(KRATOS_SD_REF_NUMBER_3)
        while(!mpStream->eof())
        #endif
        {
            ModelPartIO::ReadWord(word);
            if(ModelPartIO::CheckEndBlock("IsogeometricBezierData", word))
                break;

            BezierInfo::Pointer p_temp_info = BezierInfo::Pointer(new BezierInfo());

            std::size_t id;
            ModelPartIO::ExtractValue(word, id);
            p_temp_info->SetId(id);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->n);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->local_space_dim);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->global_space_dim);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->p1);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->p2);

            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, p_temp_info->p3);

            ModelPartIO::ReadVectorialValue(p_temp_info->weights);

            ModelPartIO::ReadWord(word);
            if(word == std::string("Full"))
                p_temp_info->mat_type = 0;
            else if(word == std::string("MCSR"))
                p_temp_info->mat_type = 1;
            else if(word == std::string("CSR"))
                p_temp_info->mat_type = 2;
            else
                KRATOS_THROW_ERROR(std::logic_error, "Invalid matrix type", word)

            if(p_temp_info->mat_type == 0)
            {
                p_temp_info->C = ModelPartIO::ReadVectorialValue(p_temp_info->C);
            }
            else if(p_temp_info->mat_type == 1)
            {
                Matrix Temp;
                Temp = ModelPartIO::ReadVectorialValue(Temp);

                // check if the input is 2 rows
                if(Temp.size1() != 2)
                    KRATOS_THROW_ERROR(std::logic_error, "Invalid MCSR matrix for extraction operator found at geometry", p_temp_info->Id())

                // choose the best storage scheme based ratio between number of nonzeros and the full size of the matrix
                unsigned int size_ex_n = (unsigned int)(Temp(0, 0) - 1);
                unsigned int size_ex_nz = Temp.size2() - 1;
                if( ( (double)(size_ex_nz) ) / (size_ex_n * size_ex_n) < 0.2 )
                    p_temp_info->C = IsogeometricMathUtils::MCSR2CSR(Temp);
                else
                    p_temp_info->C = IsogeometricMathUtils::MCSR2MAT(Temp);
            }
            else if(p_temp_info->mat_type == 2)
            {
                Vector rowPtr;
                rowPtr = ModelPartIO::ReadVectorialValue(rowPtr);
                Vector colInd;
                colInd = ModelPartIO::ReadVectorialValue(colInd);
                Vector values;
                values = ModelPartIO::ReadVectorialValue(values);
                p_temp_info->C = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
            }

//            rThisBezierInfo.push_back(p_temp_info);
            rThisBezierInfo.insert(rThisBezierInfo.begin(), p_temp_info);
            ++number_of_geometries_read;
        }

        std::cout << number_of_geometries_read << " geometries read]" << std::endl;
    }

    void BezierModelPartIO::ReadElementsWithGeometryBlock(NodesContainerType& rThisNodes,
                                                          PropertiesContainerType& rThisProperties,
                                                          BezierInfoContainerType& rGeometryInfo,
                                                          ElementsContainerType& rThisElements)
    {
        KRATOS_TRY

        SizeType id;
        SizeType properties_id;
        SizeType geometry_id;
        SizeType node_id;
        SizeType number_of_read_elements = 0;

        std::string word;
        std::string element_name;

        ModelPartIO::ReadWord(element_name);
        std::cout << "  [Reading Elements : ";

        if(!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
        Element::NodesArrayType temp_element_nodes;

        #if defined(KRATOS_SD_REF_NUMBER_2)
        while(!mFile.eof())
        #elif defined(KRATOS_SD_REF_NUMBER_3)
        while(!mpStream->eof())
        #endif
        {
            ModelPartIO::ReadWord(word); // Reading the element id or End
            if(ModelPartIO::CheckEndBlock("ElementsWithGeometry", word))
                break;

            // Reading the element id
            ModelPartIO::ExtractValue(word, id);

            // Reading the properties id
            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, properties_id);
            Properties::Pointer p_temp_properties = *(ModelPartIO::FindKey(rThisProperties, properties_id, "Properties").base());

            // Reading the geometry id
            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, geometry_id);
            BezierInfo::Pointer p_temp_info = *(ModelPartIO::FindKey(rGeometryInfo, geometry_id, "BezierInfo").base());

            // Reading the connectivities
            temp_element_nodes.clear();
            SizeType number_of_nodes = p_temp_info->n;
            for(SizeType i = 0 ; i < number_of_nodes ; i++)
            {
                // Reading the node id;
                ModelPartIO::ReadWord(word);
                ModelPartIO::ExtractValue(word, node_id);
                temp_element_nodes.push_back(*(ModelPartIO::FindKey(rThisNodes, ModelPartIO::ReorderedNodeId(node_id), "Node").base()));
            }

            // create new geometry with bezier info
            typedef BezierModelPartIO::NodeType NodeType;
            typedef IsogeometricGeometry<NodeType> IsogeometricGeometryType;
            typename IsogeometricGeometryType::Pointer p_temp_geometry;

//            if(p_temp_info->local_space_dim == 2 && p_temp_info->global_space_dim == 2)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo2dBezier<NodeType>(temp_element_nodes));
//            else if(p_temp_info->local_space_dim == 3 && p_temp_info->global_space_dim == 3)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo3dBezier<NodeType>(temp_element_nodes));
//            else if(p_temp_info->local_space_dim == 2 && p_temp_info->global_space_dim == 3)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo2dBezier3<NodeType>(temp_element_nodes));
            p_temp_geometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(r_clone_element.GetGeometry().Create(temp_element_nodes));
            if (p_temp_geometry == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to IsogeometricGeometry is failed.", "")

            Vector dummy;
            int max_integration_method = (*p_temp_properties)[NUM_IGA_INTEGRATION_METHOD];
//            KRATOS_WATCH(max_integration_method)
            p_temp_geometry->AssignGeometryData(dummy,
                                                dummy,
                                                dummy,
                                                p_temp_info->weights,
                                                p_temp_info->C,
                                                p_temp_info->p1,
                                                p_temp_info->p2,
                                                p_temp_info->p3,
                                                max_integration_method);

            rThisElements.push_back(r_clone_element.Create(ReorderedElementId(id), p_temp_geometry, p_temp_properties));
            ++number_of_read_elements;
        }
        std::cout << number_of_read_elements << " elements read] [Type: " << element_name << "]" << std::endl;
        rThisElements.Unique();

        KRATOS_CATCH("")
    }

    void BezierModelPartIO::ReadConditionsWithGeometryBlock(NodesContainerType& rThisNodes,
                                                            PropertiesContainerType& rThisProperties,
                                                            BezierInfoContainerType& rGeometryInfo,
                                                            ConditionsContainerType& rThisConditions)
    {
        KRATOS_TRY

        SizeType id;
        SizeType properties_id;
        SizeType geometry_id;
        SizeType node_id;
        SizeType number_of_read_conditions = 0;

        std::string word;
        std::string condition_name;

        ModelPartIO::ReadWord(condition_name);
        std::cout << "  [Reading Conditions : ";

        if(!KratosComponents<Condition>::Has(condition_name))
        {
            std::stringstream buffer;
            buffer << "Condition " << condition_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the condition name and see if the application which containing it, is registered corectly.";
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
        Condition::NodesArrayType temp_condition_nodes;

        #if defined(KRATOS_SD_REF_NUMBER_2)
        while(!mFile.eof())
        #elif defined(KRATOS_SD_REF_NUMBER_3)
        while(!mpStream->eof())
        #endif
        {
            ModelPartIO::ReadWord(word); // Reading the condition id or End
            if(ModelPartIO::CheckEndBlock("ConditionsWithGeometry", word))
                break;

            // Reading the element id
            ModelPartIO::ExtractValue(word, id);

            // Reading the properties id
            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, properties_id);
            Properties::Pointer p_temp_properties = *(ModelPartIO::FindKey(rThisProperties, properties_id, "Properties").base());

            // Reading the geometry id
            ModelPartIO::ReadWord(word);
            ModelPartIO::ExtractValue(word, geometry_id);
            BezierInfo::Pointer p_temp_info = *(ModelPartIO::FindKey(rGeometryInfo, geometry_id, "BezierInfo").base());

            // Reading the connectivities
            temp_condition_nodes.clear();
            SizeType number_of_nodes = p_temp_info->n;
            for(SizeType i = 0 ; i < number_of_nodes ; i++)
            {
                // Reading the node id;
                ModelPartIO::ReadWord(word);
                ModelPartIO::ExtractValue(word, node_id);
                temp_condition_nodes.push_back(*(ModelPartIO::FindKey(rThisNodes, ModelPartIO::ReorderedNodeId(node_id), "Node").base()));
            }

            // create new geometry with bezier info
            typedef BezierModelPartIO::NodeType NodeType;
            typedef IsogeometricGeometry<NodeType> IsogeometricGeometryType;
            typename IsogeometricGeometryType::Pointer p_temp_geometry;

//            if(p_temp_info->local_space_dim == 2 && p_temp_info->global_space_dim == 2)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo2dBezier<NodeType>(temp_condition_nodes));
//            else if(p_temp_info->local_space_dim == 3 && p_temp_info->global_space_dim == 3)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo3dBezier<NodeType>(temp_condition_nodes));
//            else if(p_temp_info->local_space_dim == 2 && p_temp_info->global_space_dim == 3)
//                p_temp_geometry = IsogeometricGeometryType::Pointer(new Geo2dBezier3<NodeType>(temp_condition_nodes));
            p_temp_geometry = boost::dynamic_pointer_cast<IsogeometricGeometryType>(r_clone_condition.GetGeometry().Create(temp_condition_nodes));
            if (p_temp_geometry == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to IsogeometricGeometry is failed.", "")

            Vector dummy;
            int max_integration_method = (*p_temp_properties)[NUM_IGA_INTEGRATION_METHOD];
//            KRATOS_WATCH(max_integration_method)
            p_temp_geometry->AssignGeometryData(dummy,
                                                dummy,
                                                dummy,
                                                p_temp_info->weights,
                                                p_temp_info->C,
                                                p_temp_info->p1,
                                                p_temp_info->p2,
                                                p_temp_info->p3,
                                                max_integration_method);

            rThisConditions.push_back(r_clone_condition.Create(ReorderedConditionId(id), p_temp_geometry, p_temp_properties));
            ++number_of_read_conditions;

        }
        std::cout << number_of_read_conditions << " conditions read] [Type: " << condition_name << "]" << std::endl;
        rThisConditions.Unique();

        KRATOS_CATCH("")
    }

}

