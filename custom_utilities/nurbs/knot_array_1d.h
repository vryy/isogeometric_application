//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED

// System includes
#include <vector>
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
    typedef std::vector<knot_t> knot_container_t;
    typedef typename knot_container_t::iterator iterator;
    typedef typename knot_container_t::const_iterator const_iterator;

    struct ValueComparator
    {
    private:
        const TDataType epsilon_;

    public:

        explicit ValueComparator(TDataType tolerance) : epsilon_(tolerance)
        {}

        bool operator() (const TDataType a, const TDataType b) const {
            if (std::abs(a - b) <= epsilon_) {
                return false;
            }

            return a < b;
        }
    };

    /// Default constructor
    KnotArray1D() : mTol(1e-10) {}

    /// Copy constructor
    KnotArray1D(const KnotArray1D& rOther)
    {
        for (auto it = rOther.mpKnots.begin(); it != rOther.mpKnots.end(); ++it)
        {
            knot_t p_knot = (*it)->Clone();
            this->mpKnots.push_back(p_knot);
        }
        this->mTol = rOther.GetResolution();
    }

    /// Destructor
    virtual ~KnotArray1D() {}

    /// Clear the internal knot container
    void clear()
    {
        mpKnots.clear();
    }

    /// Set the resolution of the knot vector
    void SetResolution(const TDataType& rValue) {mTol = rValue;}

    /// Get the resolution of the knot vector
    TDataType GetResolution() const {return mTol;}

    /// Insert the knot to the array and return its pointer.
    /// This function creates the new knot regardless it is repetitive or not.
    /// After all the knot is inserted, it's important to call UpdateIndex() to index the internal array
    knot_t pCreateKnot(const TDataType& k)
    {
        // insert to the correct location
        iterator it;
        for (it = mpKnots.begin(); it != mpKnots.end(); ++it)
            if (k < (*it)->Value())
                break;
        knot_t p_knot = knot_t(new KnotType(k));
        mpKnots.insert(it, p_knot);

        return p_knot;
    }

    /// Update the index of the knot
    void UpdateIndex()
    {
        std::size_t index = 0;
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it, ++index)
            (*it)->UpdateIndex(index);
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
    knot_t pCreateUniqueKnot(const TDataType& k)
    {
        // insert to the correct location
        for (iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            if (std::abs(k - (*it)->Value()) < mTol)
            {
                return *it;
            }
        }
        return pCreateKnot(k);
    }

    /// Get the knot of specific value
    knot_t pGetKnot(const TDataType& k) const
    {
        // insert to the correct location
        for (const_iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            if (std::abs(k - (*it)->Value()) < mTol)
            {
                return *it;
            }
        }

        KRATOS_ERROR << "Knot value " << k << " does not exist";

        return nullptr;
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
        if (dir == BoundaryDirection::_FORWARD_)
        {
            KnotArray1D<TDataType> kvec(*this);
            return kvec;
        }
        else if (dir == BoundaryDirection::_REVERSED_)
        {
            KnotArray1D<TDataType> kvec(*this);
            kvec.Reverse();
            return kvec;
        }
    }

    /// Create a clone of this knot vector and insert the knot in the middle of each span.
    /// It is noted that knot span bounded by repetitive knots are not accounted
    KnotArray1D<TDataType> CloneAndRefineInTheMiddle() const
    {
        KnotArray1D<TDataType> kvec(*this);
        for (const_iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
        {
            const_iterator it2 = it + 1;
            if (it2 != mpKnots.end())
            {
                if (std::abs((*it2)->Value() - (*it)->Value()) > mTol)
                {
                    TDataType ins_knot = 0.5 * ((*it)->Value() + (*it2)->Value());
                    kvec.pCreateUniqueKnot(ins_knot);
                }
            }
        }
        kvec.UpdateIndex();
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
    template<typename TVectorType>
    void GetValues(TVectorType& r_values) const
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

    /// Return the inner values of the knot vector, assuming open one
    template<typename TVectorType>
    void GetInnerValues(TVectorType& r_values, const unsigned int order) const
    {
        int size = mpKnots.size() - 2*(order+1);
        if (size < 0)
            KRATOS_ERROR << "Invalid order " << order;

        if (r_values.size() != size)
            r_values.resize(size);

        for (std::size_t i = 0; i < size; ++i)
        {
            r_values[i] = mpKnots[i+order+1]->Value();
        }
    }

    /// Return the non-repeated values of the knot vector
    template<typename TVectorType>
    void GetNonRepeatedValues(TVectorType& r_values) const
    {
        ValueComparator comp(mTol);
        std::set<TDataType, ValueComparator> nr_values(comp);

        for (std::size_t i = 0; i < mpKnots.size(); ++i)
            nr_values.insert(mpKnots[i]->Value());

        r_values.clear();
        r_values.reserve(nr_values.size());
        for (auto v : nr_values)
            r_values.push_back(v);
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
    bool IsSymmetric() const
    {
        std::vector<TDataType> knot_vec = this->GetValues();
        return IsSymmetric(knot_vec, mTol);
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
                if (std::abs((*sorted_vec)[size - 1 - i] - (*sorted_vec)[i]) > tol)
                {
                    return false;
                }
            }
        }
        else
        {
            std::size_t half_size = (size - 1) / 2;
            if (std::abs((*sorted_vec)[half_size] - 0.5) > tol)
            {
                return false;
            }
            for (std::size_t i = 0; i < half_size; ++i)
            {
                if (std::abs((*sorted_vec)[size - 1 - i] - (*sorted_vec)[i]) > tol)
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

    /// Find the index number of a shape function providing the local knot vector
    std::size_t FindNumber(const std::vector<knot_t>& pLocalKnots) const
    {
        if (pLocalKnots.size() > mpKnots.size())
        {
            KRATOS_ERROR << "The local knot is larger than the global knot vector";
            return -1;
        }

        for (std::size_t i = 0; i < mpKnots.size() - pLocalKnots.size() + 1; ++i)
        {
            bool found = true;

            for (std::size_t j = 0; j < pLocalKnots.size(); ++j)
            {
                if (std::abs(pLocalKnots[j]->Value() - mpKnots[i+j]->Value()) > mTol)
                {
                    found = false;
                    break;
                }
            }

            if (found)
                return i;
        }

        KRATOS_ERROR << "Can't find the local knot vector within the global knot vector";

        return -1;
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
            knot_t p_knot = (*it)->Clone(); // be careful, this will create new knot in memory
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

    /// Print the value only
    void PrintValue(std::ostream& rOStream) const
    {
        rOStream << "[";
        for (const_iterator it = mpKnots.begin(); it != mpKnots.end(); ++it)
            rOStream << (*it)->Value() << ", ";
        rOStream << "]";
    }

private:

    knot_container_t mpKnots;
    TDataType mTol; // defines the resolution of the knot vector, in which two knots are considered
                    // non-repetitive if their difference is larger than mTol

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

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_ARRAY_1D_H_INCLUDED
