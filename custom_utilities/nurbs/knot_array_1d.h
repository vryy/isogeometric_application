//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED

// System includes
#include <deque>
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/nurbs/knot.h"

namespace Kratos
{

/**
This container manages the knot array in 1D. It allows for easy insertion, extraction of knot from the array

Short description:
+   mpKnots is always sorted ascending.
+   the index of knot starts from 0.
+   this container stores the array of pointers to the knot, not the knot value itself.
 */
template<typename TDataType>
class KnotArray1D
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(KnotArray1D);

    /// Type definitions
    typedef Knot<TDataType> KnotType;
    typedef TDataType value_type; // this is to be in consistent with std::vector
    typedef typename KnotType::Pointer knot_t;
    typedef typename KnotType::ConstPointer const_knot_t;
    typedef std::deque<knot_t> knot_container_t;
    typedef typename knot_container_t::iterator iterator;
    typedef typename knot_container_t::const_iterator const_iterator;

    /// Default constructor
    KnotArray1D() {}

    /// Copy constructor
    KnotArray1D(const KnotArray1D& rOther)
    {
        for (auto it = rOther.mpKnots.begin(); it != rOther.mpKnots.end(); ++it)
        {
            knot_t p_knot = (*it)->Clone();
            this->mpKnots.push_back(p_knot);
        }
    }

    /// Destructor
    virtual ~KnotArray1D() {}

    /// Clear the internal knot container
    void clear()
    {
        mpKnots.clear();
    }

    /// Insert the knot to the array and return its pointer.
    /// This function creates the new knot regardless it is repetitive or not.
    knot_t pCreateKnot(const TDataType& k)
    {
        // insert to the correct location
        iterator it;
        for (it = mpKnots.begin(); it != mpKnots.end(); ++it)
            if (k < (*it)->Value())
            {
                break;
            }
        knot_t p_knot = knot_t(new KnotType(k));
        mpKnots.insert(it, p_knot);

        // update the index of the knot
        std::size_t index = 0;
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            (*it)->UpdateIndex(index);
            ++index;
        }

        return p_knot;
    }

    /// Normalize the knot vector so that the largest knot is 1.0
    void Normalize()
    {
        TDataType kmax = -1.0;
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            if ((*it)->Value() > kmax)
            {
                kmax = (*it)->Value();
            }
        }
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            (*it)->Value() /= kmax;
        }
    }

    /// Insert the knot to the array and return its pointer.
    /// In the case that the knot are repetitive within the tolerance, return the internal one.
    knot_t pCreateUniqueKnot(const TDataType& k, const TDataType& tol)
    {
        // insert to the correct location
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            if (fabs(k - (*it)->Value()) < tol)
            {
                return *it;
            }
        }
        return pCreateKnot(k);
    }

    /// Reverse this knot
    void Reverse()
    {
        TDataType maxv = mpKnots.back()->Value();
        std::reverse(mpKnots.begin(), mpKnots.end());
        // update the index and value of the knot
        std::size_t index = 0;
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            (*it)->UpdateIndex(index);
            (*it)->Value() = maxv - (*it)->Value();
            ++index;
        }
    }

    /// Create a clone of this knot vector
    KnotArray1D<TDataType> Clone(const BoundaryDirection dir = BoundaryDirection::_FORWARD_) const
    {
        KnotArray1D<TDataType> kvec;

        if (dir == BoundaryDirection::_FORWARD_)
        {
            for (std::size_t i = 0; i < mpKnots.size(); ++i)
            {
                kvec.pCreateKnot(mpKnots[i]->Value());
            }
        }
        else if (dir == BoundaryDirection::_REVERSED_)
        {
            for (std::size_t i = 0; i < mpKnots.size(); ++i)
            {
                kvec.pCreateKnot(mpKnots.back()->Value() - mpKnots[mpKnots.size() - 1 - i]->Value());
            }
        }

        return kvec;
    }

    /// Get the knot at index i
    const knot_t pKnotAt(std::size_t i) const
    {
        if (i >= 0 && i < mpKnots.size())
        {
            return mpKnots[i];
        }
        else
            KRATOS_ERROR << "Index " << i << " is out of range";
    }

    /// Get the knot at index i
    knot_t pKnotAt(std::size_t i)
    {
        if (i >= 0 && i < mpKnots.size())
        {
            return mpKnots[i];
        }
        else
            KRATOS_ERROR << "Index " << i << " is out of range";
    }

    /// Get the size of the knot vector
    std::size_t size() const {return mpKnots.size();}

    /// Get the number of knot spans
    std::size_t nspans() const
    {
        std::size_t nspan = 0;
        TDataType left = (*begin())->Value();
        for (const_iterator it = begin(); it != end(); ++it)
        {
            TDataType right = (*it)->Value();
            if (right != left)
            {
                ++nspan;
                left = right;
            }
        }
        return nspan;
    }

    /// Get the two knots bounded the span (the closest one)
    std::tuple<knot_t, knot_t> span(std::size_t i_span) const
    {
        std::size_t i = 0;
        knot_t left = *begin();
        for (const_iterator it = begin(); it != end(); ++it)
        {
            knot_t right = (*it);
            if (right->Value() != left->Value())
            {
                ++i;
                if (i == i_span)
                {
                    return std::make_tuple(left, right);
                }
                left = right;
            }
            else
            {
                left = right;    // move the knot
            }
        }
        // shall not come here
        KRATOS_ERROR << "the span index exceeds the number of span of the knot vector";
    }

    /// Return the values of the knot vector
    void GetValues(std::vector<TDataType>& r_values) const
    {
        if (r_values.size() != mpKnots.size())
        {
            r_values.resize(mpKnots.size());
        }

        for (std::size_t i = 0; i < mpKnots.size(); ++i)
        {
            r_values[i] = mpKnots[i]->Value();
        }
    }

    /// Get the knot repeated number
    /// For example, if the knot vector is [0 0 0 0.5 0.5 1 1 1] it should return [3 2 3]
    std::vector<unsigned int> GetKnotRepeatedNumber() const
    {
        std::vector<unsigned int> numbers{};
        if (size() < 1) return numbers;

        TDataType last_value;
        unsigned int cnt = 0;
        for (std::size_t i = 0; i < mpKnots.size(); ++i)
        {
            TDataType current_value = mpKnots[i]->Value();
            if (cnt == 0)
            {
                last_value = current_value;
                ++cnt;
                continue;
            }

            if (current_value > last_value)
            {
                numbers.push_back(cnt);
                last_value = current_value;
                cnt = 1;
            }
            else
            {
                ++cnt;
            }

            if (i == mpKnots.size() - 1)
                numbers.push_back(cnt);
        }

        return numbers;
    }

    /// Return the values of the knot vector
    std::vector<TDataType> GetValues() const
    {
        std::vector<TDataType> values;
        this->GetValues(values);
        return values;
    }

    /// Iterator
    iterator begin() {return mpKnots.begin();}

    /// Iterator
    iterator end() {return mpKnots.end();}

    /// Iterator
    const_iterator begin() const {return mpKnots.begin();}

    /// Iterator
    const_iterator end() const {return mpKnots.end();}

    /// Check if this knot vector is symmetric within a specified tolerance
    bool IsSymmetric(const TDataType& tol) const
    {
        std::vector<TDataType> knot_vec = this->GetValues();
        return IsSymmetric(knot_vec, tol);
    }

    /// Check if a knot vector is symmetric. A knot vector is symmetric if both x and 1-x are in the knot vector
    static bool IsSymmetric(const std::vector<TDataType>& knot_vec, const TDataType& tol, bool sorted = true)
    {
        std::vector<TDataType>* sorted_vec;
        if (sorted == false)
        {
            sorted_vec = new std::vector<TDataType>(knot_vec.size());
            std::copy(knot_vec.begin(), knot_vec.end(), sorted_vec->begin());
            std::sort(sorted_vec->begin(), sorted_vec->end());
        }
        else
        {
            sorted_vec = const_cast<std::vector<TDataType>*>(&knot_vec);
        }

        std::size_t size = sorted_vec->size();
        if (size % 2 == 0)
        {
            std::size_t half_size = size / 2;
            for (std::size_t i = 0; i < half_size; ++i)
            {
                if (fabs((*sorted_vec)[size - 1 - i] - (*sorted_vec)[i]) > tol)
                {
                    return false;
                }
            }
        }
        else
        {
            std::size_t half_size = (size - 1) / 2;
            if (fabs((*sorted_vec)[half_size] - 0.5) > tol)
            {
                return false;
            }
            for (std::size_t i = 0; i < half_size; ++i)
            {
                if (fabs((*sorted_vec)[size - 1 - i] - (*sorted_vec)[i]) > tol)
                {
                    return false;
                }
            }
        }

        if (sorted == false)
        {
            delete sorted_vec;
        }

        return true;
    }

    /// Check if a knot is inside the knot vector
    bool IsInside(const TDataType& knot) const
    {
        const TDataType& min = mpKnots.front()->Value();
        const TDataType& max = mpKnots.back()->Value();
        return (knot >= min) && (knot <= max);
    }

    /// Compute the knots in the respective direction
    static std::vector<TDataType> CloneKnots(const std::vector<TDataType>& knots, const BoundaryDirection direction)
    {
        if (direction == BoundaryDirection::_FORWARD_)
        {
            return knots;
        }
        else if (direction == BoundaryDirection::_REVERSED_)
        {
            return ReverseKnots(knots);
        }
        return std::vector<TDataType> {};
    }

    /// Compute the knots in the respective direction
    static std::vector<TDataType> CloneKnotsWithPivot(const TDataType& pivot, const std::vector<TDataType>& knots, const BoundaryDirection direction)
    {
        if (direction == BoundaryDirection::_FORWARD_)
        {
            return knots;
        }
        else if (direction == BoundaryDirection::_REVERSED_)
        {
            return ReverseKnotsWithPivot(pivot, knots);
        }
        return std::vector<TDataType> {};
    }

    /// Compute the reversed knots, i.e. 1-k
    static std::vector<TDataType> ReverseKnotsWithPivot(const TDataType& pivot, const std::vector<TDataType>& knots)
    {
        std::vector<TDataType> reversed_knots(knots.size());
        for (std::size_t i = 0; i < reversed_knots.size(); ++i)
        {
            reversed_knots[i] = pivot - knots[i];
        }
        return reversed_knots;
    }

    /// Compute the reversed knots, i.e. 1-k
    static std::vector<TDataType> ReverseKnots(const std::vector<TDataType>& knots)
    {
        return ReverseKnotsWithPivot(knots.back(), knots);
    }

    /// Check if a local knot vector is on the left side. It is on the left side if the outer left knots are repeated (p+1) times.
    /// It is assumed that the input local knot vector must have size (p+2).
    static bool IsOnLeft(const std::vector<knot_t>& knots, std::size_t p)
    {
        for (std::size_t i = 0; i < p; ++i)
            if (knots[i + 1]->Value() != knots[i]->Value())
            {
                return false;
            }
        return true;
    }

    /// Check if a local knot vector is on the right side. It is on the right side if the outer right knots are repeated (p+1) times.
    /// It is assumed that the input local knot vector must have size (p+2).
    static bool IsOnRight(const std::vector<knot_t>& knots, std::size_t p)
    {
        std::size_t last = knots.size() - 1;
        for (std::size_t i = 0; i < p; ++i)
            if (knots[last - i]->Value() != knots[last - 1 - i]->Value())
            {
                return false;
            }
        return true;
    }

    /// Compare the two knot vectors
    bool operator==(const KnotArray1D<TDataType>& rOther) const
    {
        if (this->size() != rOther.size())
        {
            return false;
        }

        for (std::size_t i = 0; i < 0; ++i)
        {
            if (this->pKnotAt(i)->Value() != rOther.pKnotAt(i)->Value())
            {
                return false;
            }
            else
            {
                KRATOS_ERROR << "The knot vector is different at loc " << i
                             << ": " << this->pKnotAt(i)->Value() << " != " << rOther.pKnotAt(i)->Value();
            }
        }

        return true;
    }

    /// Compare the two knot vectors
    bool operator==(const std::vector<knot_t>& rOther) const
    {
        if (this->size() != rOther.size())
        {
            return false;
        }

        for (std::size_t i = 0; i < 0; ++i)
        {
            if (this->pKnotAt(i)->Value() != rOther[i]->Value())
            {
                return false;
            }
            else
            {
                KRATOS_ERROR << "The knot vector is different at loc " << i
                             << ": " << this->pKnotAt(i)->Value() << " != " << rOther.pKnotAt(i)->Value();
            }
        }

        return true;
    }

    // overload assignment operator
    KnotArray1D& operator=(const KnotArray1D& rOther)
    {
        for (auto it = rOther.mpKnots.begin(); it != rOther.mpKnots.end(); ++it)
        {
            knot_t p_knot = (*it)->Clone();
            this->mpKnots.push_back(p_knot);
        }
        return *this;
    }

    // overload operator ()
    knot_t operator() (std::size_t i)
    {
        return pKnotAt(i);
    }

    // overload operator ()
    const_knot_t operator() (std::size_t i) const
    {
        return pKnotAt(i);
    }

    // overload operator []
    TDataType& operator[] (std::size_t i)
    {
        return pKnotAt(i)->Value();
    }

    // overload operator []
    const TDataType& operator[] (std::size_t i) const
    {
        return pKnotAt(i)->Value();
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        for (const_iterator it = begin(); it != end(); ++it)
        {
            rOStream << " (" << (*it)->Index() << "," << (*it)->Value() << ")";
        }
    }

private:

    knot_container_t mpKnots;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save( "size", mpKnots.size() );

        for (const_iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            rSerializer.save( "K", *(*it) );
        }
    }

    void load(Serializer& rSerializer)
    {
        std::size_t size;
        rSerializer.load( "size", size );

        KnotType knot;
        for (std::size_t i = 0; i < size; ++i)
        {
            rSerializer.load( "K", knot );
            knot_t p_knot = knot_t(new KnotType(knot));
            mpKnots.insert(mpKnots.end(), p_knot);
        }
    }
    ///@}
};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const KnotArray1D<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED
