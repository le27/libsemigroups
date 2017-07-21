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
#include <unordered_set>
#include <algorithm>

#include "timer.h"
#include "semigroups.h"



namespace libsemigroups {

  template <typename ElementType, typename PointType> class Orb {
   public:
    Orb(std::vector<ElementType*> gens, PointType seed)
        : mover_map({}),
          _gens(gens), _orb({seed}),
          _map({std::make_pair(seed, 0)}) {}

    std::unordered_map<PointType, ElementType*>  mover_map;

    void enumerate() {
      for (size_t i = 0; i < _orb.size(); i++) {
        for (ElementType const* x : _gens) {
          PointType pt = (*x)[_orb[i]];
          if (_map.find(pt) == _map.end()) {
            _map.insert(std::make_pair(pt, _orb.size()));
            _orb.push_back(pt);
            ElementType* y = new ElementType(*x);
            if (mover_map.find(pt) == mover_map.end()) {
                mover_map.insert(std::make_pair(pt, y));
            }
            else {
                mover_map[pt]->redefine(mover_map[pt], y);
            }
            delete y;
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
      mover_map.reserve(n);
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

  template <typename ElementType, typename PointType> class DeltaU {
   public:
    DeltaU(std::unordered_set<ElementType*> strong_gen_set, PointType seed)
        : delta(), u(), new_gens({}), _strong_gen_set(strong_gen_set), _seed(seed) {}
    std::vector<PointType> delta;
    std::vector<ElementType*> new_gens;
    std::unordered_map<PointType, ElementType*> u;

    void enumerate() {
        for (ElementType* elt : _strong_gen_set) {
           if ((*elt)[_seed] == _seed)
               new_gens.push_back(elt);
        }
        Orb<ElementType, PointType> orb = Orb<ElementType, PointType>(new_gens, _seed);
        orb.enumerate();
        delta = orb._orb;
        u = orb.mover_map;
    }


    private:
     std::unordered_set<ElementType*>     _strong_gen_set;
     PointType                     _seed;
     // TODO const?
   };

  template <typename ElementType, typename PointType> class BSGS {
   public:
    BSGS(std::vector<ElementType*> gens)
        : base({}), strong_gen_set({}), _gens(gens) {}

    std::vector<PointType>        base;
    std::unordered_set<ElementType*> strong_gen_set;

    void enumerate_partial() {
        for (ElementType* elt : _gens){
            if (!(elt->is_identity()))
                strong_gen_set.insert(elt);
        }
        std::unordered_set<ElementType*> tempstrong_gen_set = strong_gen_set;
        bool baseelt_equal_base = true;

        for (ElementType* elt : tempstrong_gen_set) {
            for (PointType b: base) {
                baseelt_equal_base = true;
                if (std::find(base.begin(), base.end(), (*elt)[b]) == base.end()) {
                    baseelt_equal_base = false;
                }
            }
            if (baseelt_equal_base) {
                for (PointType point: *elt){
                    if (!(point == (*elt)[point])) {
                        base.push_back(point);
                        break;
                    }
                }
            }
            if (!(*elt == *(elt->inverse())))
                strong_gen_set.insert(elt->inverse());
        }
    }

    std::pair<ElementType*, size_t> strip(ElementType* g,
     std::vector<std::vector<PointType>>* delta,
     std::vector<std::unordered_map<PointType, ElementType*>>* u){
         ElementType* g_copy = new ElementType(*g);
         for (size_t l = 1; l <= base.size(); l++){
             if (std::find(delta[l].begin(), delta[l].end(),
                *g_copy[base[l]]) != delta[l].end())
                g_copy.redefine(g_copy, ((*u)[*g_copy[base[l]]])->inverse());
             else
                return std::make_pair(g_copy, l);
         }
         return make_pair(g_copy, base.size() + 1);
     }


     void enumerate(){
         this->enumerate_partial();
         for (size_t i = base.size(); i > 0; i = i - 1){
             DeltaU<ElementType, PointType> delta_u
              = DeltaU<ElementType, PointType>(strong_gen_set, base[i]);
             delta_u.enumerate();
             for (PointType pt : delta_u.delta) {
                 for (ElementType* elt : delta_u.new_gens){
                     ElementType* g = new ElementType(*((delta_u.u)[pt]));
                     g->redefine(g, elt);
                     g->redefine(g, (delta_u.u[(*elt)[pt]])->inverse());
                     if (!(g->is_identity())) {
                         strong_gen_set.insert(g);
                         strong_gen_set.insert(g->inverse());
                         bool g_fixes_base = true;
                         for (PointType b : base) {
                             if (!((*g)[b] == b)) {
                                 g_fixes_base = false;
                                 break;
                             }
                         }
                         if (g_fixes_base) {
                             for (PointType new_base_pt : *g) {
                                 if (!((*g)[new_base_pt] == new_base_pt)) {
                                     base.push_back(new_base_pt);
                                     break;
                                 }
                             }
                         }
                     }
                 }
             }
         }
     }


   private:
    std::vector<ElementType*>     _gens;
    // TODO const?
  };

  // template <typename ElementType, typename PointType> class StabChain {
  //  public:
  //   StabChain(std::vector<ElementType*> gens, std::vector<PointType> base)
  //       : _gens(gens), _base(base), _gen_chain({}) {}
  //
  //   void enumerate() {
  //     for (size_t i = 0; i < _orb.size(); i++) {
  //       for (ElementType const* x : _gens) {
  //         PointType pt = (*x)[_orb[i]];
  //         if (_map.find(pt) == _map.end()) {
  //           _map.insert(std::make_pair(pt, _orb.size()));
  //           _orb.push_back(pt);
  //         }
  //       }
  //     }
  //   }
  //
  //   size_t position(PointType pt) {
  //     auto it = _map.find(pt);
  //     if (it != _map.end()) {
  //       return (*it).second;
  //     } else {
  //       return -1;
  //     }
  //   }
  //
  //   void reserve(size_t n) {
  //     _orb.reserve(n);
  //     _map.reserve(n);
  //   }
  //
  //   size_t size() {
  //     return _orb.size();
  //   }
  //
  //  private:
  //   std::vector<ElementType*>     _gens;
  //   std::vector<PointType>        _base;
  //   std::vector<std::vector<ElementType*>>   _gen_chain;
  //   // TODO const?
  // };

}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_SRC_SIMS_H_
