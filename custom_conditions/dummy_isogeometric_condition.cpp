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
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "includes/legacy_structural_app_vars.h"
#endif
#include "custom_conditions/dummy_isogeometric_condition.h"

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

void DummyIsogeometricCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

#ifdef SD_APP_FORWARD_COMPATIBILITY
    // borrow the variable from other application
    const Variable<int>& INTEGRATION_ORDER_var = static_cast<const Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));
#else
    const Variable<int>& INTEGRATION_ORDER_var = INTEGRATION_ORDER;
#endif

    // integration rule
    if (this->Has( INTEGRATION_ORDER_var ))
    {
        if (this->GetValue(INTEGRATION_ORDER_var) == 1)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if (this->GetValue(INTEGRATION_ORDER_var) == 2)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if (this->GetValue(INTEGRATION_ORDER_var) == 3)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if (this->GetValue(INTEGRATION_ORDER_var) == 4)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if (this->GetValue(INTEGRATION_ORDER_var) == 5)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricCondition does not support for integration rule", this->GetValue(INTEGRATION_ORDER_var))
        }
    else if (GetProperties().Has( INTEGRATION_ORDER_var ))
    {
        if (GetProperties()[INTEGRATION_ORDER_var] == 1)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if (GetProperties()[INTEGRATION_ORDER_var] == 2)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if (GetProperties()[INTEGRATION_ORDER_var] == 3)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if (GetProperties()[INTEGRATION_ORDER_var] == 4)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if (GetProperties()[INTEGRATION_ORDER_var] == 5)
        {
            mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricCondition does not support for integration points", GetProperties()[INTEGRATION_ORDER_var])
        }
    else
    {
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();    // default method
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void DummyIsogeometricCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
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
        const ProcessInfo& rCurrentProcessInfo)
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
        const ProcessInfo& rCurrentProcessInfo,
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
        const ProcessInfo& CurrentProcessInfo) const
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
void DummyIsogeometricCondition::GetDofList( DofsVectorType& ConditionalDofList,
        const ProcessInfo& CurrentProcessInfo) const
{
    ConditionalDofList.resize(0);
}

} // Namespace Kratos
