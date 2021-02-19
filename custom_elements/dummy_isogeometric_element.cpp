//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jan 2021$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/dummy_isogeometric_element.h"
#include "includes/legacy_structural_app_vars.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "isogeometric_application/custom_utilities/isogeometric_math_utils.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyIsogeometricElement::DummyIsogeometricElement()
{
}

DummyIsogeometricElement::DummyIsogeometricElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

DummyIsogeometricElement::DummyIsogeometricElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyIsogeometricElement::~DummyIsogeometricElement()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer DummyIsogeometricElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyIsogeometricElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer DummyIsogeometricElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyIsogeometricElement(NewId, pGeom, pProperties));
}

GeometryData::IntegrationMethod DummyIsogeometricElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

void DummyIsogeometricElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricElement does not support for integration rule", this->GetValue(INTEGRATION_ORDER))
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
            KRATOS_THROW_ERROR(std::logic_error, "DummyIsogeometricElement does not support for integration points", GetProperties()[INTEGRATION_ORDER])
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
void DummyIsogeometricElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
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
void DummyIsogeometricElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
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
void DummyIsogeometricElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
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
void DummyIsogeometricElement::EquationIdVector( EquationIdVectorType& rResult,
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
void DummyIsogeometricElement::GetDofList( DofsVectorType& ConditionalDofList,
                                      const ProcessInfo& CurrentProcessInfo) const
{
    ConditionalDofList.resize(0);
}

} // Namespace Kratos

