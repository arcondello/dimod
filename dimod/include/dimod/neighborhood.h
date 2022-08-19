// Copyright 2021 D-Wave Systems Inc.
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#pragma once

#include <algorithm>
#include <utility>
#include <vector>

namespace dimod {

/**
 * Used internally by QuadraticModelBase to sparsely encode the neighborhood of
 * a variable.
 *
 * Internally, Neighborhoods keep two vectors, one of neighbors and the other
 * of biases. Neighborhoods are designed to make access more like a standard
 * library map.
 */
template <class Bias, class Index>
class Neighborhood {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    /// Exactly `pair<index_type, bias_type>`.
    using value_type = typename std::pair<index_type, bias_type>;

    /// A random access iterator to `pair<index_type, bias_type>`.
    using iterator = typename std::vector<value_type>::iterator;

    /// A random access iterator to `const pair<index_type, bias_type>.`
    using const_iterator = typename std::vector<value_type>::const_iterator;

    /**
     * Return a reference to the bias associated with `v`.
     *
     * This function automatically checks whether `v` is a variable in the
     * neighborhood and throws a `std::out_of_range` exception if it is not.
     */
    bias_type at(index_type v) const {
        auto it = this->lower_bound(v);
        if (it != this->cend() && it->first == v) {
            // it exists
            return it->second;
        } else {
            // it doesn't exist
            throw std::out_of_range("given variable has no interaction");
        }
    }

    value_type& back() { return this->neighborhood_.back(); }

    const value_type& back() const { return this->neighborhood_.back(); }

    /// Returns an iterator to the beginning.
    iterator begin() { return this->neighborhood_.begin(); }

    /// Returns an iterator to the end.
    iterator end() { return this->neighborhood_.end(); }

    /// Returns a const_iterator to the beginning.
    const_iterator cbegin() const { return this->neighborhood_.cbegin(); }

    /// Returns a const_iterator to the end.
    const_iterator cend() const { return this->neighborhood_.cend(); }

    /**
     * Insert a neighbor, bias pair at the end of the neighborhood.
     *
     * Note that this does not keep the neighborhood self-consistent and should
     * only be used when you know that the neighbor is greater than the current
     * last element.
     */
    void emplace_back(index_type v, bias_type bias) { this->neighborhood_.emplace_back(v, bias); }

    /// Returns whether the neighborhood is empty
    bool empty() const { return !this->size(); }

    /**
     * Erase an element from the neighborhood.
     *
     * Returns the number of element removed, either 0 or 1.
     */
    size_type erase(index_type v) {
        auto it = this->lower_bound(v);
        if (it != this->end() && it->first == v) {
            // is there to erase
            this->neighborhood_.erase(it);
            return 1;
        } else {
            return 0;
        }
    }

    /// Erase elements from the neighborhood.
    void erase(iterator first, iterator last) { this->neighborhood_.erase(first, last); }

    /// Return an iterator to the first element that does not come before `v`.
    iterator lower_bound(index_type v) {
        return std::lower_bound(this->begin(), this->end(), v, this->cmp);
    }

    /// Return an iterator to the first element that does not come before `v`.
    const_iterator lower_bound(index_type v) const {
        return std::lower_bound(this->cbegin(), this->cend(), v, this->cmp);
    }

    /**
     * Total bytes consumed by the biases and indices.
     *
     * If `capacity` is true, use the capacity of the underlying vectors rather
     * than the size.
     */
    size_type nbytes(bool capacity = false) const noexcept {
        // so there is no guaruntee that the compiler will not implement
        // pair as pointers or whatever, but this seems like a reasonable
        // assumption.
        if (capacity) {
            return this->neighborhood_.capacity() * sizeof(std::pair<index_type, bias_type>);
        } else {
            return this->neighborhood_.size() * sizeof(std::pair<index_type, bias_type>);
        }
    }

    /**
     * Return the bias at neighbor `v` or the default value.
     *
     * Return the bias of `v` if `v` is in the neighborhood, otherwise return
     * the `value` provided without inserting `v`.
     */
    bias_type get(index_type v, bias_type value = 0) const {
        auto it = this->lower_bound(v);

        if (it != this->cend() && it->first == v) {
            // it exists
            return it->second;
        } else {
            // it doesn't exist
            return value;
        }
    }

    /// Request that the neighborhood capacity be at least enough to contain `n`
    /// elements.
    void reserve(index_type n) { this->neighborhood_.reserve(n); }

    /// Return the size of the neighborhood.
    size_type size() const { return this->neighborhood_.size(); }

    /// Sort the neighborhood and sum the biases of duplicate variables.
    void sort_and_sum() {
        if (!std::is_sorted(this->begin(), this->end())) {
            std::sort(this->begin(), this->end());
        }

        // now remove any duplicates, summing the biases of duplicates
        size_type i = 0;
        size_type j = 1;

        // walk quickly through the neighborhood until we find a duplicate
        while (j < this->neighborhood_.size() &&
               this->neighborhood_[i].first != this->neighborhood_[j].first) {
            ++i;
            ++j;
        }

        // if we found one, move into de-duplication
        while (j < this->neighborhood_.size()) {
            if (this->neighborhood_[i].first == this->neighborhood_[j].first) {
                this->neighborhood_[i].second += this->neighborhood_[j].second;
                ++j;
            } else {
                ++i;
                this->neighborhood_[i] = this->neighborhood_[j];
                ++j;
            }
        }

        // finally resize to contain only the unique values
        this->neighborhood_.resize(i + 1);
    }

    /**
     * Access the bias of `v`.
     *
     * If `v` is in the neighborhood, the function returns a reference to
     * its bias. If `v` is not in the neighborhood, it is inserted and a
     * reference is returned to its bias.
     */
    bias_type& operator[](index_type v) {
        auto it = this->lower_bound(v);
        if (it == this->end() || it->first != v) {
            // it doesn't exist so insert
            it = this->neighborhood_.emplace(it, v, 0);
        }
        return it->second;
    }

 protected:
    std::vector<value_type> neighborhood_;

    static inline bool cmp(value_type ub, index_type v) { return ub.first < v; }
};

}  // namespace dimod
