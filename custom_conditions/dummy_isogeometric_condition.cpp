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
#include "includes/legacy_structural_app_vars.h"
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

GeometryData::IntegrationMethod DummyIsogeometricCondition::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

void DummyIsogeometricCondition::Initialize()
{
    KRATOS_TRY

    // integration rule
    if(this->Has( INTEGRATION_ORDER ))
    {
        if(this->GetValue(INTEGRATION_ORDER) == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(this->GetValue(INTEGRATION_ORDER) == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricCondition does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
    }
    else if(GetProperties().Has( INTEGRATION_ORDER ))
    {
        if(GetProperties()[INTEGRATION_ORDER] == 1)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 2)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 3)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 4)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if(GetProperties()[INTEGRATION_ORDER] == 5)
        {
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricCondition does not support for integration points", GetProperties()[INTEGRATION_ORDER])
    }
    else
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

    KRATOS_CATCH("")
}

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
void DummyIsogeometricCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    ConditionalDofList.resize(0);
}

} // Namespace Kratos

