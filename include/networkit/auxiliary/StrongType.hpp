/*
 * StrongType.hpp
 *
 * Date: 28.06.2019
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef STRONGTYPE_HPP_
#define STRONGTYPE_HPP_

#include <algorithm>
#include <type_traits>

namespace Aux {

/**
 * Helper function to strongly type values passed. It's intended use is to define a named
 * specialization with the "using" keyword. Since two specializations of the same value
 * type are treated identically, the second "dummy" parameter can be used to tag the types
 * and make them distinguishable.
 *
 * \code
 * using Radius = StrongType<double, struct RadiusTag>;
 * using AverageDegree = StrongType<double, struct AverageDegreeTag>;
 *
 * // Two overloads that behave differently based on the type of parameter provided.
 * void findParameters(Radius radius);
 * void findParameters(AverageDegree avgDeg);
 * \endcode
 */
template <typename ValueType, typename UniqueTag = void>
class StrongType {
public:
    using value_type = ValueType;

    StrongType() = delete;
    StrongType(const StrongType&) = default;
    StrongType(StrongType&&) = default;
    StrongType& operator=(const StrongType&) = default;
    StrongType& operator=(StrongType&&) = default;

    ///! Copy data into wrapper
    explicit StrongType(value_type const& value) noexcept(std::is_nothrow_copy_constructible<value_type>::value)
        : value(value)
    {}

    ///! Moves data into wrapper
    explicit StrongType(value_type&& value) noexcept(std::is_nothrow_move_constructible<value_type>::value)
        : value(std::move(value))
    {}

    ///! Returns a reference to the value contained
    value_type& get() noexcept {
        return value;
    }

    ///! Returns a const reference to the value contained
    value_type const& get() const noexcept {
        return value;
    }

private:
    value_type value;

};

}

#endif // STRONGTYPE_HPP_
