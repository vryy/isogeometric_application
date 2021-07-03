//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 18, 2013 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "custom_utilities/control_point.h"
#include "custom_elements/dummy_isogeometric_element.h"
#include "custom_conditions/dummy_isogeometric_condition.h"

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
    class KratosIsogeometricApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{


        /// Pointer definition of KratosDiscontinuitiesApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosIsogeometricApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosIsogeometricApplication();

        /// Destructor.
        virtual ~KratosIsogeometricApplication()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register();



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
            return "KratosIsogeometricApplication";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            KRATOS_WATCH("in KratosIsogeometricApplication");
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
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

        // const DummyIsogeometricElement mDummyElementBezier;
        // const DummyIsogeometricElement mDummyElementBezier2D;
        // const DummyIsogeometricElement mDummyElementBezier3D;

        // const DummyIsogeometricCondition mDummyConditionBezier;
        // const DummyIsogeometricCondition mDummyConditionBezier2D;
        // const DummyIsogeometricCondition mDummyConditionBezier2D3;
        // const DummyIsogeometricCondition mDummyConditionBezier3D;


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
        KratosIsogeometricApplication& operator=(KratosIsogeometricApplication const& rOther);

        /// Copy constructor.
        KratosIsogeometricApplication(KratosIsogeometricApplication const& rOther);


        ///@}

    }; // Class KratosIsogeometricApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


}  // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_H_INCLUDED  defined


