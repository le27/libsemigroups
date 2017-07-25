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
//
#ifndef LIBSEMIGROUPS_SRC_SIMS_H_
#define LIBSEMIGROUPS_SRC_SIMS_H_

#include <algorithm>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "semigroups.h"
#include "timer.h"

namespace libsemigroups {
  //
  //   template <typename ElementType, typename PointType> class Orb {
  //    public:
  //     Orb(std::vector<ElementType*> gens, PointType seed)
  //         : mover_map({}),
  //           _gens(gens),
  //           _orb({seed}),
  //           _map({std::make_pair(seed, 0)}),
  //           _gen({nullptr}),
  //           _parent({-1}) {
  //       assert(!gens.empty());
  //       if (_tmp.degree() < _gens[0].degree()) {
  //         _tmp.copy(_gens[0]);
  //       }
  //     }
  //
  //     ~Orb();
  //
  //     void enumerate() {
  //       for (size_t i = 0; i < _orb.size(); i++) {
  //         for (ElementType const* x : _gens) {
  //           PointType pt = (*x)[_orb[i]];
  //           if (_map.find(pt) == _map.end()) {
  //             _map.insert(std::make_pair(pt, _orb.size()));
  //             _orb.push_back(pt);
  //             _gen.push_back(x);
  //             _parent.push_back(i);
  //           }
  //         }
  //       }
  //     }
  //
  //     // TODO change the name to find
  //     size_t position(PointType pt) {
  //       auto it = _map.find(pt);
  //       if (it != _map.end()) {
  //         return (*it).second;
  //       } else {
  //         return -1;
  //       }
  //     }
  //
  //     void reserve(size_t n) {
  //       _map.reserve(n);
  //       _orb.reserve(n);
  //       _gen.reserve(n);
  //     }
  //
  //     inline PointType operator[](size_t pos) const {
  //       assert(pos < _orb.size());
  //       return _orb[pos];
  //     }
  //
  //     inline PointType at(size_t pos) const {
  //       return _orb.at(pos);
  //     }
  //
  //     size_t size() {
  //       return _orb.size();
  //     }
  //
  //     Element const* mapper(size_t pos) {
  //       if (pos >= _mappers.size()) {
  //         size_t i = pos;
  //         assert(i < _orb.size());
  //         Element const* out = _gen[i].really_copy();
  //         Element const* tmp = out.really_copy();
  //         while (out != nullptr) {
  //           i = _parent[i];
  //           out->redefine(_gen[i], tmp);
  //           tmp.copy(out);
  //         }
  //         _mappers[i] = out;
  //       }
  //       return _mappers[pos];
  //     }
  //
  //    private:
  //     std::vector<ElementType*> _gens;
  //     std::unordered_map<PointType, size_t> _map;
  //     std::vector<PointType>    _orb;
  //     std::vector<ElementType*> _gen;
  //     std::vector<size_t> _parent;
  //     std::vector<ElementType*> _mappers;
  //     static ElementType* _tmp = new ElementType(new std::vector());
  //   };
  //
  //

  //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  // Schreier-Sims set up

  class PermColl {
   public:
    std::vector<Permutation<size_t>*>* gens;
    size_t                             nr_gens;
    size_t                             deg;
    size_t                             alloc_size;

    void really_delete() {
      for (Permutation<size_t>* x : *gens) {
        x->really_delete();
        delete x;
      }
      delete &nr_gens;
      delete &deg;
      delete &alloc_size;
    }

    PermColl* really_copy() {
      PermColl* out;
      for (Permutation<size_t>* x : *gens)
        out->gens->push_back(
            static_cast<Permutation<size_t>*>(x->really_copy()));
      out->nr_gens    = this->nr_gens;
      out->deg        = this->deg;
      out->alloc_size = this->alloc_size;
      return out;
    }

    void add_perm_coll(Permutation<size_t>* gen) {
      assert(this->nr_gens <= this->alloc_size);

      if (this->nr_gens == this->alloc_size) {
        this->gens->resize(this->nr_gens + 1 * sizeof(Permutation<size_t>));
        (this->alloc_size)++;
        // nr_ss_allocs++;
        // nr_ss_frees++;
      }
      (*(this->gens))[(this->nr_gens)] = gen;
    }
  };

  PermColl* new_perm_coll(size_t upper_bound, size_t deg) {
    // nr_new_perm_coll++;
    PermColl* coll = new PermColl;
    // nr_ss_allocs++;
    coll->gens = new std::vector<Permutation<size_t>*>;
    // nr_ss_allocs++;
    coll->nr_gens    = 0;
    coll->deg        = deg;
    coll->alloc_size = upper_bound;
    return coll;
  }

#define MAXVERTS 512

  static PermColl*            strong_gens[MAXVERTS];
  static Permutation<size_t>* transversal[MAXVERTS * MAXVERTS];
  static Permutation<size_t>* transversal_inv[MAXVERTS * MAXVERTS];
  static bool                 first_ever_call = true;
  static bool                 borbits[MAXVERTS * MAXVERTS];
  static size_t               orbits[MAXVERTS * MAXVERTS];
  static size_t               size_orbits[MAXVERTS];
  static size_t               base[MAXVERTS];
  static size_t               size_base;
  static size_t               deg;

  static inline void add_strong_gens(size_t const               pos,
                                     Permutation<size_t>* const value) {
    if (strong_gens[pos] == nullptr) {
      strong_gens[pos] = new_perm_coll(1, deg);
    }
    strong_gens[pos]->add_perm_coll(value);
  }

  static inline Permutation<size_t>* get_strong_gens(size_t const i,
                                                     size_t const j) {
    return (*(strong_gens[i]->gens))[j];
  }

  static inline Permutation<size_t>* get_transversal(size_t const i,
                                                     size_t const j) {
    return transversal[i * MAXVERTS + j];
  }

  static inline Permutation<size_t>* get_transversal_inv(size_t const i,
                                                         size_t const j) {
    return transversal_inv[i * MAXVERTS + j];
  }

  // TODO: Add back (or not) nr_ss_frees (used to debug membory leaks)
  static inline void set_transversal(size_t const               i,
                                     size_t const               j,
                                     Permutation<size_t>* const value) {
    // free the perm in this position if there is one already
    if (transversal[i * MAXVERTS + j] != nullptr) {
      delete transversal[i * MAXVERTS + j];
      // nr_ss_frees++;
      delete transversal_inv[i * MAXVERTS + j];
      // nr_ss_frees++;
    }
    transversal[i * MAXVERTS + j]     = value;
    transversal_inv[i * MAXVERTS + j] = value->inverse();
  }

  static bool perm_fixes_all_base_points(Permutation<size_t>* const x) {
    for (size_t i = 0; i < size_base; i++) {
      if ((*x)[base[i]] != base[i]) {
        return false;
      }
    }
    return true;
  }

  static inline void add_base_point(size_t const pt) {
    base[size_base]               = pt;
    size_orbits[size_base]        = 1;
    orbits[size_base * deg]       = pt;
    borbits[size_base * deg + pt] = true;
    set_transversal(size_base,
                    pt,
                    static_cast<Permutation<size_t>*>(
                        (*((strong_gens)[0]->gens))[0]->identity()));
    size_base++;
  }

  static inline void first_ever_init() {
    first_ever_call = false;
    std::memset((void*) size_orbits, 0, MAXVERTS * sizeof(size_t));
  }

  static void init_stab_chain() {
    if (first_ever_call) {
      first_ever_init();
    }

    std::memset((void*) borbits, false, deg * deg * sizeof(bool));
    size_base = 0;
  }

  static void free_stab_chain() {
    int i, j, k;

    std::memset((void*) size_orbits, 0, size_base * sizeof(size_t));

    // free the transversal
    // free the transversal_inv
    for (i = 0; i < (int) deg; i++) {
      for (j = 0; j < (int) deg; j++) {
        k = i * MAXVERTS + j;
        if (transversal[k] != nullptr) {
          transversal[k]->really_delete();
          transversal[k] = nullptr;
          // nr_ss_frees++;
          transversal_inv[k]->really_delete();
          transversal_inv[k] = nullptr;
          // nr_ss_frees++;
        }
      }
    }

    // free the strong_gens
    for (i = 0; i < (int) size_base; i++) {
      if (strong_gens[i] != nullptr) {
        strong_gens[i]->really_delete();
        strong_gens[i] = nullptr;
      }
    }
  }

  static void orbit_stab_chain(size_t const depth, size_t const init_pt) {
    size_t               i, j, pt, img;
    Permutation<size_t>* temp;
    Permutation<size_t>* x;

    assert(depth <= size_base);  // Should this be strict?

    for (i = 0; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->nr_gens; j++) {
        x   = get_strong_gens(depth, j);
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          temp->redefine(get_transversal(depth, pt), x);
          set_transversal(depth, img, temp);
        }
      }
    }
  }

  static void add_gen_orbit_stab_chain(size_t const               depth,
                                       Permutation<size_t>* const gen) {
    size_t               i, j, pt, img;
    Permutation<size_t>* x;
    Permutation<size_t>* temp;

    assert(depth <= size_base);

    // apply the new generator to existing points in orbits[depth]
    size_t nr = size_orbits[depth];
    for (i = 0; i < nr; i++) {
      pt  = orbits[depth * deg + i];
      img = (*gen)[pt];
      if (!borbits[depth * deg + img]) {
        orbits[depth * deg + size_orbits[depth]] = img;
        size_orbits[depth]++;
        borbits[depth * deg + img] = true;
        temp->redefine(get_transversal(depth, pt), x);
        set_transversal(depth, img, temp);
      }
    }

    for (i = nr; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->nr_gens; j++) {
        x   = get_strong_gens(depth, j);
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          temp->redefine(get_transversal(depth, pt), x);
          set_transversal(depth, img, temp);
        }
      }
    }
  }

  static void sift_stab_chain(Permutation<size_t>* g, size_t* depth) {
    size_t               beta;
    Permutation<size_t>* temp;
    assert(*depth == 0);

    for (; *depth < size_base; (*depth)++) {
      beta = (*g)[base[*depth]];
      if (!borbits[*depth * deg + beta]) {
        return;
      }
      temp = static_cast<Permutation<size_t>*>(g->really_copy());
      g->redefine(temp, get_transversal_inv(*depth, beta));
    }
  }

  static void schreier_sims_stab_chain(size_t const depth) {
    Permutation<size_t>* x;
    Permutation<size_t>* h;
    Permutation<size_t>* prod;
    bool                 escape, y;
    int                  i;
    size_t               j, jj, k, l, m, beta, betax;

    for (i = 0; i <= (int) depth; i++) {
      for (j = 0; j < strong_gens[i]->nr_gens; j++) {
        x = get_strong_gens(i, j);
        if (perm_fixes_all_base_points(x)) {
          for (k = 0; k < deg; k++) {
            if (k != (*x)[k]) {
              add_base_point(k);
              break;
            }
          }
        }
      }
    }

    for (i = depth + 1; i < (int) size_base + 1; i++) {
      beta = base[i - 1];
      // set up the strong generators
      for (j = 0; j < strong_gens[i - 1]->nr_gens; j++) {
        x = get_strong_gens(i - 1, j);
        if (beta == (*x)[beta]) {
          add_strong_gens(i,
                          static_cast<Permutation<size_t>*>(x->really_copy()));
        }
      }

      // find the orbit of <beta> under strong_gens[i - 1]
      orbit_stab_chain(i - 1, beta);
    }

    i = size_base - 1;

    while (i >= (int) depth) {
      escape = false;
      for (j = 0; j < size_orbits[i] && !escape; j++) {
        beta = orbits[i * deg + j];
        for (m = 0; m < strong_gens[i]->nr_gens && !escape; m++) {
          x = get_strong_gens(i, m);
          prod->redefine(get_transversal(i, beta), x);
          betax = (*x)[beta];
          if (!(prod == get_transversal(i, betax))) {
            y = true;
            h->redefine(prod, get_transversal_inv(i, betax));
            jj = 0;
            sift_stab_chain(h, &jj);
            if (jj < size_base) {
              y = false;
            } else if (!(h->is_identity())) {  // better method? IsOne(h)?
              y = false;
              for (k = 0; k < deg; k++) {
                if (k != (*h)[k]) {
                  add_base_point(k);
                  break;
                }
              }
            }

            if (!y) {
              for (l = i + 1; l <= jj; l++) {
                add_strong_gens(
                    l, static_cast<Permutation<size_t>*>(h->really_copy()));
                add_gen_orbit_stab_chain(l, h);
                // add generator to <h> to orbit of base[l]
              }
              i      = jj;
              escape = true;
            }
            h->really_delete();
            //          nr_ss_frees++;
          }
          prod->really_delete();
          //        nr_ss_frees++;
        }
      }
      if (!escape) {
        i--;
      }
    }
  }
  extern bool point_stabilizer(PermColl* gens, size_t const pt, PermColl* out) {
    init_stab_chain();
    strong_gens[0] = gens->really_copy();
    add_base_point(pt);
    schreier_sims_stab_chain(0);

    // The stabiliser we want is the PermColl pointed to by <strong_gens[1]>
    // UNLESS <strong_gens[1]> doesn't exists - this means that
    // <strong_gens[0]> is the stabilizer itself (????)
    if (out != nullptr) {
      out->really_delete();
    }
    if (strong_gens[1] == nullptr) {
      // this means that the stabilizer of pt under <gens> is trivial
      out = new_perm_coll(1, deg);
      out->add_perm_coll(static_cast<Permutation<size_t>*>(
          (*((strong_gens)[0]->gens))[0]->identity()));
      free_stab_chain();
      return true;
    }
    out = strong_gens[1]->really_copy();
    free_stab_chain();
    return false;
  }

}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_SRC_SIMS_H_
