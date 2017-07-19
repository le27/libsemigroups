//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2017 James D. Mitchell
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

// This file contains ...

#ifndef LIBSEMIGROUPS_SRC_SIMS_H_
#define LIBSEMIGROUPS_SRC_SIMS_H_

#include <iostream>
#include <unordered_map>
#include <vector>

#include "timer.h"

namespace libsemigroups {
  template <typename ElementType, typename PointType> class Orb {
   public:
    Orb(std::vector<ElementType*> gens, PointType seed)
        : _gens(gens), _orb({seed}), _map({std::make_pair(seed, 0)}) {}

    void enumerate() {
      for (size_t i = 0; i < _orb.size(); i++) {
        for (ElementType const* x : _gens) {
          PointType pt = (*x)[_orb[i]];
          if (_map.find(pt) == _map.end()) {
            _map.insert(std::make_pair(pt, _orb.size()));
            _orb.push_back(pt);
          }
        }
      }
    }

    size_t position(PointType pt) {
      auto it = _map.find(pt);
      if (it != _map.end()) {
        return (*it).second;
      } else {
        return -1;
      }
    }

    void reserve(size_t n) {
      _orb.reserve(n);
      _map.reserve(n);
    }

    size_t size() {
      return _orb.size();
    }

   private:
    std::vector<ElementType*>     _gens;
    std::vector<PointType>        _orb;
    std::unordered_map<PointType, size_t> _map;
    // TODO const?
  };
}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_SRC_SIMS_H_
