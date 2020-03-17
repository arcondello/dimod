// Copyright 2019 D-Wave Systems Inc.
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

#ifndef DIMOD_SHAPEABLEBQM_H_
#define DIMOD_SHAPEABLEBQM_H_

// #include <algorithm>
#include <map>
#include <utility>
#include <vector>

namespace dimod {

template<class Neighborhood, class V, class B>
class cppShapeableBQM {  // temp name
 public:
    using bias_type = B;
    using variable_type = V;
    using size_type = std::size_t;

    using outvars_iterator = typename Neighborhood::iterator;
    using const_outvars_iterator = typename Neighborhood::const_iterator;

    // in the future we'd probably like to make this protected
    std::vector<std::pair<Neighborhood, B>> adj;

    cppShapeableBQM(): num_interactions_(0) {}

    size_type num_variables() const {
        return adj.size();
    }

    size_type num_interactions() const {
        return num_interactions_;
    }

    bias_type get_linear(variable_type v) const {
        return adj[v].second;
    }

 protected:
    size_type num_interactions_;
};


template<class V, class B>
class cppAdjMapBQM : public cppShapeableBQM<std::map<V, B>, V, B> {  // temp name
 public:
    cppAdjMapBQM() : cppShapeableBQM<std::map<V, B>, V, B>() {}

    template<class BQM>
    explicit cppAdjMapBQM(const BQM &bqm) : cppShapeableBQM<std::map<V, B>, V, B>() {
        adj.resize(bqm.num_variables());
    }
};

template<class V, class B>
class cppAdjVectorBQM : public cppShapeableBQM<std::vector<std::pair<V, B>>, V, B> {  // temp name
};


}  // namespace dimod

#endif  // DIMOD_SHAPEABLEBQM_H_
