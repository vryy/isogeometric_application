//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Aug 2015 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_BASIS_FUNCTION_MANAGER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_BASIS_FUNCTION_MANAGER_H_INCLUDED

// System includes
#include <string>
#include <vector>
// #include <unordered_map>
#include <map>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/knot.h"

namespace Kratos
{

/**
  TODO
 */
template<class TBasisFuncType>
class DeprecatedBasisFunctionManager
{
public:
    /// Pointer definition
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(DeprecatedBasisFunctionManager);

    /// Type definitions
    typedef Knot<double> KnotType;
    typedef KnotType::Pointer knot_t;

    typedef typename TBasisFuncType::Pointer bf_t;
    struct bf_compare
    {
        bool operator() (const bf_t& lhs, const bf_t& rhs) const {return lhs->Id() < rhs->Id();}
    };
    typedef std::set<bf_t, bf_compare> bf_container_t;
    typedef typename bf_container_t::iterator iterator;

//    typedef std::unordered_map<std::size_t, bf_t> map_t;
    typedef std::map<std::size_t, bf_t> map_t;

    /// Default constructor
    DeprecatedBasisFunctionManager() : mLastId(0)
    {}

    /// Destructor
    virtual ~DeprecatedBasisFunctionManager()
    {}

    /// Check if the bf exists in the list; otherwise create new bf and return
    typename TBasisFuncType::Pointer CreateBf(unsigned int Level,
                                              const std::vector<knot_t>& rpKnots1,
                                              const std::vector<knot_t>& rpKnots2,
                                              const std::vector<knot_t>& rpKnots3)
    {
        // search in the current list of bfs, the one that has the same local knot vector with provided ones
        for(iterator it = mpBasisFuncs.begin(); it != mpBasisFuncs.end(); ++it)
            if((*it)->Contain(rpKnots1, rpKnots2, rpKnots3))
                return *it;

        // create the new bf and add the knot
        typename TBasisFuncType::Pointer p_bf = typename TBasisFuncType::Pointer(new TBasisFuncType(++mLastId, Level));
        p_bf->SetLocalKnotVectors(1, rpKnots1);
        p_bf->SetLocalKnotVectors(2, rpKnots2);
        p_bf->SetLocalKnotVectors(3, rpKnots3);
        mpBasisFuncs.insert(p_bf);
        function_map_is_created = false;

        return p_bf;
    }

    /// Iterators
    iterator begin() {return mpBasisFuncs.begin();}
    iterator begin() const {return mpBasisFuncs.begin();}
    iterator end() {return mpBasisFuncs.end();}
    iterator end() const {return mpBasisFuncs.end();}

    /// Get the number of basis functions in this container
    std::size_t size() const
    {
        return mpBasisFuncs.size();
    }

    /// Remove a basis function by its Id from the set
    void erase(bf_t p_bf)
    {
        for(iterator it = mpBasisFuncs.begin(); it != mpBasisFuncs.end(); ++it)
            if((*it) == p_bf)
            {
                mpBasisFuncs.erase(it);
                break;
            }
    }

    /// Get a basis function based on its Id
    bf_t get(const std::size_t& Id)
    {
        // create the index map if it's not created yet
        if(!function_map_is_created)
            CreateFunctionsMap();

        // return the bf if its Id exist in the list
        typename map_t::iterator it = mFunctionsMap.find(Id);
        if(it != mFunctionsMap.end())
            return it->second;
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Access index is not found:", Id)
    }

    /// Overload operator[]
    bf_t operator[](const std::size_t& Id)
    {
        return get(Id);
    }

    /// Reset all the Id of all the basis functions. Remarks: use it with care, you have to be responsible to the old indexing data of the basis functions before calling this function
    /// Temporary disable it for safety
//    std::size_t ReIndexing()
//    {
//        mLastId = 0;
//        for(iterator it = mpBasisFuncs.begin(); it != mpBasisFuncs.end(); ++it)
//            (*it)->SetId(++mLastId);
//        return mLastId;
//    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
    std::size_t mLastId;
    bf_container_t mpBasisFuncs;

    mutable map_t mFunctionsMap; // map from basis function id to the basis function. It's mainly used to search for the bf quickly. But it needs to be re-initialized whenever new bf is added to the set
    bool function_map_is_created;

    void CreateFunctionsMap()
    {
        mFunctionsMap.clear();
        for(iterator it = mpBasisFuncs.begin(); it != mpBasisFuncs.end(); ++it)
            mFunctionsMap[(*it)->Id()] = *it;
        function_map_is_created = true;
    }

};

/// output stream function
template<class TBasisFuncType>
inline std::ostream& operator <<(std::ostream& rOStream, const DeprecatedBasisFunctionManager<TBasisFuncType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_DEPRECATED_BASIS_FUNCTION_MANAGER_H_INCLUDED

