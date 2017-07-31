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

// Schreier-Sims set up

#define MAXVERTS 512

  static std::vector<std::vector<Permutation<size_t>*>*> strong_gens(MAXVERTS);
  static std::vector<Permutation<size_t>*> transversal(MAXVERTS* MAXVERTS);
  static std::vector<Permutation<size_t>*> transversal_inv(MAXVERTS* MAXVERTS);
  static std::vector<bool>                 borbits(MAXVERTS* MAXVERTS);
  static std::vector<size_t>               orbits(MAXVERTS* MAXVERTS);
  static std::vector<size_t>               size_orbits(MAXVERTS);
  static std::vector<size_t>               base(MAXVERTS);
  static size_t                            size_base;
  static size_t                            deg;

  // deletes the contents of a contatainer
  template <typename T> static inline void really_delete_cont(T* cont) {
    for (Permutation<size_t>* x : *cont) {
      x->really_delete();
      delete x;
    }
  }

  // multiplies permutations, creating a new object.
  static Permutation<size_t>* prod_perms(Permutation<size_t>* x,
                                         Permutation<size_t>* y) {
    if (x == nullptr)
      return static_cast<Permutation<size_t>*>(y->really_copy());
    if (y == nullptr)
      return static_cast<Permutation<size_t>*>(x->really_copy());
    std::vector<size_t> z(deg);
    for (size_t i = 0; i < deg; i++) {
      z[i] = (*y)[(*x)[i]];
    }
    return new Permutation<size_t>(z);
  }

  // copies a vector of permutations
  static inline std::vector<Permutation<size_t>*>*
  really_copy_vec(std::vector<Permutation<size_t>*>* cont) {
    std::vector<Permutation<size_t>*>* out
        = new std::vector<Permutation<size_t>*>({});
    for (Permutation<size_t>* x : *cont) {
      out->push_back(static_cast<Permutation<size_t>*>(x->really_copy()));
    }
    return out;
  }

  // adds a permutaion to the posth position of strong_gens
  static inline void add_strong_gens(size_t const               pos,
                                     Permutation<size_t>* const value) {
    if (strong_gens[pos] == nullptr) {
      strong_gens[pos] = new std::vector<Permutation<size_t>*>({});
    }
    strong_gens[pos]->push_back(value);
  }

  // TODO: Add back (or not) nr_ss_frees (used to debug membory leaks)

  // stores a permutation in the transversal vector
  static inline void set_transversal(size_t const               i,
                                     size_t const               j,
                                     Permutation<size_t>* const value) {
    // // free the perm in this position if there is one already
    // if (transversal[i * MAXVERTS + j] != nullptr) {
    //   transversal[i * MAXVERTS + j]->really_delete();
    //   delete transversal[i * MAXVERTS + j];
    //   // nr_ss_frees++;
    //   transversal_inv[i * MAXVERTS + j]->really_delete();
    //   delete transversal_inv[i * MAXVERTS + j];
    //   // nr_ss_frees++;
    // }
    transversal[i * MAXVERTS + j]     = value;
    transversal_inv[i * MAXVERTS + j] = value->inverse();
  }

  // checks if x fixes all points in the base
  static bool perm_fixes_all_base_points(Permutation<size_t>* const x) {
    for (size_t i = 0; i < size_base; i++) {
      if ((*x)[base[i]] != base[i]) {
        return false;
      }
    }
    return true;
  }

  // adds a point to the base (and initialises orbit stuff)
  static inline void add_base_point(size_t const pt) {
    base[size_base]               = pt;
    size_orbits[size_base]        = 1;
    orbits[size_base * deg]       = pt;
    borbits[size_base * deg + pt] = true;
    set_transversal(size_base,
                    pt,
                    static_cast<Permutation<size_t>*>(
                        ((*(strong_gens[0]))[0]->identity())));
    size_base++;
  }

  // empties some variables
  static void free_stab_chain() {
    int i, j, k;

    // std::memset((void*) size_orbits, 0, size_base * sizeof(size_t));
    // free the transversal
    // free the transversal_inv
    for (i = 0; i < (int) deg; i++) {
      for (j = 0; j < (int) deg; j++) {
        k = i * MAXVERTS + j;
        if (transversal[k] != nullptr) {
          transversal[k]->really_delete();
          delete transversal[k];
          transversal[k] = nullptr;
          // nr_ss_frees++;
          transversal_inv[k]->really_delete();
          delete transversal[k];
          transversal_inv[k] = nullptr;
          // nr_ss_frees++;
        }
      }
    }

    // free the strong_gens
    for (i = 0; i < (int) size_base; i++) {
      if (strong_gens[i] != nullptr) {
        really_delete_cont(strong_gens[i]);
        delete strong_gens[i];
        strong_gens[i] = nullptr;
      }
    }
  }

  // calculates the orbit of init_pt in the `depth`th stabiliser
  static void orbit_stab_chain(size_t const depth) {
    size_t               i, j, pt, img;
    Permutation<size_t>* x;
    assert(depth <= size_base);  // Should this be strict?

    for (i = 0; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->size(); j++) {
        x   = (*(strong_gens[depth]))[j];
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          set_transversal(
              depth, img, prod_perms(transversal[depth * MAXVERTS + pt], x));
        }
      }
    }
  }

  // apply the new generator to existing points in orbits[depth]
  static void add_gen_orbit_stab_chain(size_t const               depth,
                                       Permutation<size_t>* const gen) {
    size_t               i, j, pt, img;
    Permutation<size_t>* x;

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
        set_transversal(
            depth, img, prod_perms(transversal[depth * MAXVERTS + pt], gen));
      }
    }

    for (i = nr; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->size(); j++) {
        x   = (*(strong_gens[depth]))[j];
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          set_transversal(
              depth, img, prod_perms(transversal[depth * MAXVERTS + pt], x));
        }
      }
    }
  }

  // changes g so that it fixes the current base
  static void sift_stab_chain(Permutation<size_t>* g, size_t* depth) {
    size_t beta;
    assert(*depth == 0);

    for (; *depth < size_base; (*depth)++) {
      beta = (*g)[base[*depth]];
      if (!borbits[*depth * deg + beta]) {
        return;
      }
      *g = *(prod_perms(g, transversal_inv[(*depth) * MAXVERTS + beta]));
    }
    assert(perm_fixes_all_base_points(g));
  }

  // builds a strong base and generating set for the group generated by
  // strong_gens[0]  and it's stabilisers up to depth
  static void schreier_sims_stab_chain(size_t const depth) {
    Permutation<size_t>* x;
    Permutation<size_t>* h;
    Permutation<size_t>* prod;
    bool                 escape, h_fixes_all_base_pts;
    int                  i;
    size_t               j, no_base_pts_h_fixes, k, l, m, beta, betax;

    // for (Permutation<size_t>* perm : *(strong_gens[1])) {
    //   assert((*perm)[base[0]] == base[0]);
    // }

    for (i = 0; i <= (int) depth; i++) {
      if (!(strong_gens[i] == nullptr)) {
        for (Permutation<size_t>* x : *(strong_gens[i])) {
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
    }

    for (i = depth + 1; i <= (int) size_base; i++) {
      beta = base[i - 1];
      // set up the strong generators
      for (Permutation<size_t>* y : *(strong_gens[i - 1])) {
        if (beta == (*y)[beta]) {
          assert((*y)[base[0]] == base[0]);
          add_strong_gens(i,
                          static_cast<Permutation<size_t>*>(y->really_copy()));
        }
      }

      // find the orbit of base[i -1] under strong_gens[i - 1]
      orbit_stab_chain(i - 1);
    }

    i = size_base - 1;

    while (i >= (int) depth) {
      escape = false;
      for (j = 0; j < size_orbits[i] && !escape; j++) {
        beta = orbits[i * deg + j];
        for (m = 0; m < strong_gens[i]->size() && !escape; m++) {
          x     = (*(strong_gens[i]))[m];
          prod  = prod_perms(transversal[i * MAXVERTS + beta], x);
          betax = (*x)[beta];
          if (!(prod == transversal[i * MAXVERTS + betax])) {
            h_fixes_all_base_pts = true;
            h = prod_perms(prod, transversal_inv[i * MAXVERTS + betax]);
            no_base_pts_h_fixes = 0;
            sift_stab_chain(h, &no_base_pts_h_fixes);
            if ((*h)[base[0]] == base[0])
              assert(no_base_pts_h_fixes >= 1);
            if (!((*h)[base[0]] == base[0]))
              assert(no_base_pts_h_fixes < 1);
            if (no_base_pts_h_fixes < size_base) {
              h_fixes_all_base_pts = false;
            } else if (!(h->is_identity())) {  // better method? IsOne(h)?
              h_fixes_all_base_pts = false;
              for (k = 0; k < deg; k++) {
                if (k != (*h)[k]) {
                  add_base_point(k);
                  break;
                }
              }
            }

            if (!h_fixes_all_base_pts) {
              for (l = i + 1; l <= no_base_pts_h_fixes; l++) {
                add_strong_gens(
                    l, static_cast<Permutation<size_t>*>(h->really_copy()));
                if (l >= 1) {
                  assert(no_base_pts_h_fixes >= 1);
                  assert((*h)[base[0]] == base[0]);
                }
                add_gen_orbit_stab_chain(l, h);
                // add generator to <h> to orbit of base[l]
              }
              i      = no_base_pts_h_fixes;
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
    for (Permutation<size_t>* perm : *(strong_gens[1])) {
      assert((*perm)[base[0]] == base[0]);
    }
  }

  // changes out into a generating set for the stabiliser of pt in the group
  // generated by gens
  extern bool point_stabilizer(std::vector<Permutation<size_t>*>* gens,
                               size_t const                       pt,
                               std::vector<Permutation<size_t>*>* out) {
    size_base      = 0;
    strong_gens[0] = really_copy_vec(gens);
    add_base_point(pt);
    schreier_sims_stab_chain(0);

    for (Permutation<size_t>* perm : *(strong_gens[1])) {
      assert((*perm)[pt] == pt);
    }

    // The stabiliser we want is the PermColl pointed to by
    // UNLESS <strong_gens[1]> doesn't exists - this means that
    // <strong_gens[0]> is the stabilizer itself (????)
    if (out != nullptr) {
      really_delete_cont(out);
    }
    if (strong_gens[1] == nullptr) {
      // this means that the stabilizer of pt under <gens> is trivial
      *out = std::vector<Permutation<size_t>*>({});
      out->push_back(static_cast<Permutation<size_t>*>(
          (*(strong_gens[0]))[0]->identity()));
      free_stab_chain();
      return true;
    }
    *out = *really_copy_vec(strong_gens[1]);
    free_stab_chain();
    return false;
  }

  // finds the size of the permutation group generated by gens
  extern size_t group_size(std::vector<Permutation<size_t>*>* gens) {
    std::unordered_set<size_t> orbset;
    deg = ((*gens)[0])->degree();
    std::vector<size_t>                orb;
    size_t                             out = 1;
    std::vector<Permutation<size_t>*>* gens1;
    std::vector<Permutation<size_t>*>* gens2 = really_copy_vec(gens);

    // No need to check for i = deg - 1, since a permutation that fixes
    // everything  except one point must also fix that point.
    for (size_t i = 0; (i <= (deg - 1)); i++) {
      gens1 = really_copy_vec(gens2);
      point_stabilizer(gens1, i, gens2);

      for (Permutation<size_t>* perm : *gens2) {
        assert((*perm)[i] == i);
      }
      // orbit
      orb    = {i};
      orbset = {i};

      // Cannot iterate through orbit the normal way as we are adding to it.
      for (size_t j = 0; j < orb.size(); j++) {
        for (Permutation<size_t>* perm : *gens1) {
          size_t pt = (*perm)[orb[j]];
          if (orbset.find(pt) == orbset.end()) {
            orbset.insert(pt);
            orb.push_back(pt);
          }
        }
      }
      out = out * orbset.size();
    }
    really_delete_cont(gens1);
    really_delete_cont(gens2);
    delete gens1;
    delete gens2;
    for (size_t i; i < size_base; i++) {
      base[i] = 0;
    }
    size_base = 0;
    for (size_t i; i < orbits.size(); i++) {
      orbits[i] = 0;
    }
    for (size_t i; i < size_orbits.size(); i++) {
      size_orbits[i] = 0;
    }
    for (size_t i; i < borbits.size(); i++) {
      borbits[i] = false;
    }
    free_stab_chain();
    return out;
  }
}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_SRC_SIMS_H_
