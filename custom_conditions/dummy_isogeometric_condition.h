//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Dec 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED )
#define  KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_geometries/isogeometric_geometry.h"

namespace Kratos
{
/**
 * This condition does nothing. But it gives the isogeometric geometry for other condition that contains it.
 */
class DummyIsogeometricCondition : public Condition
{
    public:
        // Counted pointer of DummyIsogeometricCondition
        KRATOS_CLASS_POINTER_DEFINITION(DummyIsogeometricCondition);

        typedef GeometryData::IntegrationMethod IntegrationMethod;

        typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

        /**
         * Default constructor.
         */
        DummyIsogeometricCondition();
        DummyIsogeometricCondition( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyIsogeometricCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyIsogeometricCondition();

        /**
         * Operations.
         */

        virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

        virtual Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const;

        IntegrationMethod GetIntegrationMethod() const;

        /**
         * Calculates the local system contributions for this contact element
         */
        virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     ProcessInfo& rCurrentProcessInfo);

        virtual void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo);

        virtual void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);
        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();

        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        //virtual void PrintData(std::ostream& rOStream) const;

    private:

        IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

        IntegrationMethod mThisIntegrationMethod;

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyIsogeometricCondition

}  // namespace Kratos.


#endif // KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED defined

