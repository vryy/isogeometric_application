//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Dec 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes 

// External includes 

// Project includes 
#include "custom_conditions/dummy_isogeometric_condition.h"
#include "includes/deprecated_variables.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyIsogeometricCondition::DummyIsogeometricCondition()
{
}

DummyIsogeometricCondition::DummyIsogeometricCondition( IndexType NewId, 
                              GeometryType::Pointer pGeometry)
: Condition( NewId, pGeometry )
{
}

DummyIsogeometricCondition::DummyIsogeometricCondition( IndexType NewId, 
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Condition( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyIsogeometricCondition::~DummyIsogeometricCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer DummyIsogeometricCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyIsogeometricCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Condition::Pointer DummyIsogeometricCondition::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyIsogeometricCondition(NewId, pGeom, pProperties));
}

void DummyIsogeometricCondition::Initialize()
{
    KRATOS_TRY

//    ////////////////////Initialize geometry_data/////////////////////////////
//    #ifdef ENABLE_PROFILING
//    double start_compute = OpenMPUtils::GetCurrentTime();
//    #endif

//    // try to read the extraction operator from the elemental data
//    Matrix ExtractionOperator;
//    bool manual_initilization = false;
//    if( this->Has( EXTRACTION_OPERATOR ) )
//    {
//        ExtractionOperator = this->GetValue( EXTRACTION_OPERATOR );
//        manual_initilization = true;
//    }
//    else if( this->Has( EXTRACTION_OPERATOR_MCSR ) )
//    {
//        Matrix Temp = this->GetValue( EXTRACTION_OPERATOR_MCSR );

//        // make a simple check
//        if(Temp.size1() != 2)
//            KRATOS_THROW_ERROR(std::logic_error, "Invalid MCSR matrix for extraction operator found at element", this->Id())

//        // choose the best storage scheme based ratio between number of nonzeros and the full size of the matrix
//        unsigned int size_ex_n = (unsigned int)(Temp(0, 0) - 1);
//        unsigned int size_ex_nz = Temp.size2() - 1;
//        if( ( (double)(size_ex_nz) ) / (size_ex_n * size_ex_n) < 0.2 )
//            ExtractionOperator = IsogeometricMathUtils::MCSR2CSR(Temp);
//        else
//            ExtractionOperator = IsogeometricMathUtils::MCSR2MAT(Temp);

//        manual_initilization = true;
//    }
//    else if( this->Has( EXTRACTION_OPERATOR_CSR_ROWPTR )
//         && this->Has( EXTRACTION_OPERATOR_CSR_COLIND )
//         && this->Has( EXTRACTION_OPERATOR_CSR_VALUES ) )
//    {
//        Vector rowPtr = this->GetValue( EXTRACTION_OPERATOR_CSR_ROWPTR ); // must be 0-base
//        Vector colInd = this->GetValue( EXTRACTION_OPERATOR_CSR_COLIND ); // must be 0-base
//        Vector values = this->GetValue( EXTRACTION_OPERATOR_CSR_VALUES );
//        ExtractionOperator = IsogeometricMathUtils::Triplet2CSR(rowPtr, colInd, values);
//        manual_initilization = true;
//    }
////        else
////            KRATOS_THROW_ERROR(std::logic_error, "The extraction operator was not given for element", Id())
////        KRATOS_WATCH(ExtractionOperator)

//    // initialize the geometry
//    if(manual_initilization)
//    {
//        int num_integration_method = 2; // by default compute two integration rules
//        if( GetProperties().Has(NUM_IGA_INTEGRATION_METHOD) )
//            num_integration_method = GetProperties()[NUM_IGA_INTEGRATION_METHOD];
//        mpIsogeometricGeometry->AssignGeometryData(
//            this->GetValue(NURBS_KNOTS_1),
//            this->GetValue(NURBS_KNOTS_2),
//            this->GetValue(NURBS_KNOTS_3),
//            this->GetValue(NURBS_WEIGHT),
//            ExtractionOperator,
//            this->GetValue(NURBS_DEGREE_1),
//            this->GetValue(NURBS_DEGREE_2),
//            this->GetValue(NURBS_DEGREE_3),
//            num_integration_method
//        );
//    }

    KRATOS_CATCH("")
}

//************************************************************************************ 
//************************************************************************************
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void DummyIsogeometricCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector, 
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, 
                  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void DummyIsogeometricCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                          VectorType& rRightHandSideVector, 
                                          ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}
    //************************************************************************************
//************************************************************************************    /
/**
 * This function calculates all system contributions due to the contact problem
 * with regard to the current master and slave partners.
 * All Conditions are assumed to be defined in 2D/3D space and having 2/3 DOFs per node 
 */
void DummyIsogeometricCondition::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                  VectorType& rRightHandSideVector, 
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    rLeftHandSideMatrix.resize(0, 0, false);
    rRightHandSideVector.resize(0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.    
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void DummyIsogeometricCondition::EquationIdVector( EquationIdVectorType& rResult, 
                                      ProcessInfo& CurrentProcessInfo)
{
    rResult.resize(0);
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list for the current partners.
 * All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per Node.
 * All DOF are given Master first, Slave second
 */
//************************************************************************************
//************************************************************************************
void DummyIsogeometricCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    ConditionalDofList.resize(0);
}

} // Namespace Kratos

