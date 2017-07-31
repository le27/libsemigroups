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
  static std::vector<size_t>               base;
  static size_t                            deg;

  // deletes the contents of a contatainer
  template <typename T> static inline void really_delete_cont(T* cont) {
    for (Permutation<size_t>* x : *cont) {
      x->really_delete();
      delete x;
    }
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
    transversal[i * deg + j]     = value;
    transversal_inv[i * deg + j] = value->inverse();
  }

  // checks if x fixes all points in the base
  static bool perm_fixes_all_base_points(Permutation<size_t>* const x) {
    for (size_t i = 0; i < base.size(); i++) {
      if ((*x)[base[i]] != base[i]) {
        return false;
      }
    }
    return true;
  }

  // adds a point to the base (and initialises orbit stuff)
  static inline void add_base_point(size_t const pt) {
    size_orbits[base.size()]            = 1;
    orbits[base.size() * deg]           = pt;
    borbits[base.size() * deg + pt]     = true;
    transversal[base.size() * deg + pt] = static_cast<Permutation<size_t>*>(
        ((*(strong_gens[0]))[0]->identity()));
    transversal_inv[base.size() * deg + pt] = static_cast<Permutation<size_t>*>(
        ((*(strong_gens[0]))[0]->identity()));
    base.push_back(pt);
  }

  // empties some variables
  static void free_stab_chain() {
    size_t i, j, k;
    for (i = 0; i < deg; i++) {
      for (j = 0; j < deg; j++) {
        k = i * deg + j;
        if (transversal[k] != nullptr) {
          transversal[k]->really_delete();
          delete transversal[k];
          transversal[k] = nullptr;
          transversal_inv[k]->really_delete();
          delete transversal[k];
          transversal_inv[k] = nullptr;
        }
      }
    }

    // free the strong_gens
    for (i = 0; i < MAXVERTS; i++) {
      if (strong_gens[i] != nullptr) {
        really_delete_cont(strong_gens[i]);
        delete strong_gens[i];
        strong_gens[i] = nullptr;
      }
    }

    for (i = 0; i < orbits.size(); i++) {
      orbits[i] = 0;
    }
    for (i = 0; i < size_orbits.size(); i++) {
      size_orbits[i] = 0;
    }
    for (i = 0; i < borbits.size(); i++) {
      borbits[i] = false;
    }
  }

  // calculates the orbit of init_pt in the `depth`th stabiliser
  static void orbit_stab_chain(size_t const depth) {
    size_t               pt, img;
    Permutation<size_t>* x;
    Permutation<size_t>* prod;
    assert(depth <= base.size());  // Should this be strict?

    for (size_t i = 0; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (size_t j = 0; j < strong_gens[depth]->size(); j++) {
        x   = (*(strong_gens[depth]))[j];
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          prod = new Permutation<size_t>(new std::vector<size_t>(deg));
          prod->redefine(transversal[depth * deg + pt], x);
          set_transversal(depth, img, prod);
        }
      }
    }
  }

  // apply the new generator to existing points in orbits[depth]
  static void add_gen_orbit_stab_chain(size_t const               depth,
                                       Permutation<size_t>* const gen) {
    size_t               i, j, pt, img;
    Permutation<size_t>* x;
    Permutation<size_t>* prod;
    assert(depth <= base.size());

    // apply the new generator to existing points in orbits[depth]
    for (i = 0; i < size_orbits[depth]; i++) {
      pt  = orbits[depth * deg + i];
      img = (*gen)[pt];
      if (!borbits[depth * deg + img]) {
        orbits[depth * deg + size_orbits[depth]] = img;
        size_orbits[depth]++;
        borbits[depth * deg + img] = true;
        prod = new Permutation<size_t>(new std::vector<size_t>(deg));
        prod->redefine(transversal[depth * deg + pt], gen);
        set_transversal(depth, img, prod);
      }
    }

    for (i = size_orbits[depth]; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->size(); j++) {
        x   = (*(strong_gens[depth]))[j];
        img = (*x)[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          prod = new Permutation<size_t>(new std::vector<size_t>(deg));
          prod->redefine(transversal[depth * deg + pt], x);
          set_transversal(depth, img, prod);
        }
      }
    }
  }

  // changes g so that it fixes the current base
  static void sift_stab_chain(Permutation<size_t>* g, size_t* depth) {
    size_t beta;
    assert(*depth == 0);

    for (; *depth < base.size(); (*depth)++) {
      beta = (*g)[base[*depth]];
      if (!borbits[*depth * deg + beta])
        return;
      g->redefine(g->really_copy(), transversal_inv[(*depth) * deg + beta]);
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

    for (i = depth + 1; i <= (int) base.size(); i++) {
      beta = base[i - 1];
      // set up the strong generators
      for (Permutation<size_t>* y : *(strong_gens[i - 1])) {
        if (beta == (*y)[beta]) {
          add_strong_gens(i,
                          static_cast<Permutation<size_t>*>(y->really_copy()));
        }
      }

      // find the orbit of base[i -1] under strong_gens[i - 1]
      orbit_stab_chain(i - 1);
    }

    i = base.size() - 1;

    while (i >= (int) depth) {
      escape = false;
      for (j = 0; j < size_orbits[i] && !escape; j++) {
        beta = orbits[i * deg + j];
        for (m = 0; m < strong_gens[i]->size() && !escape; m++) {
          x    = (*(strong_gens[i]))[m];
          prod = new Permutation<size_t>(std::vector<size_t>(deg));
          prod->redefine(transversal[i * deg + beta], x);
          betax = (*x)[beta];
          if (!(*prod == *(transversal[i * deg + betax]))) {
            h_fixes_all_base_pts = true;
            h = new Permutation<size_t>(std::vector<size_t>(deg));
            h->redefine(prod, transversal_inv[i * deg + betax]);
            no_base_pts_h_fixes = 0;
            sift_stab_chain(h, &no_base_pts_h_fixes);
            if (no_base_pts_h_fixes < base.size()) {
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
                add_gen_orbit_stab_chain(l, h);
                // add generator to <h> to orbit of base[l]
              }
              i      = no_base_pts_h_fixes;
              escape = true;
            }
            h->really_delete();
            delete h;
          }
          prod->really_delete();
          delete prod;
        }
      }
      if (!escape) {
        i--;
      }
    }
  }

  // changes out into a generating set for the stabiliser of pt in the group
  // generated by gens
  extern bool point_stabilizer(std::vector<Permutation<size_t>*>* const gens,
                               size_t const                             pt,
                               std::vector<Permutation<size_t>*>*       out) {
    base           = {};
    strong_gens[0] = really_copy_vec(gens);
    add_base_point(pt);
    schreier_sims_stab_chain(0);

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
  extern size_t group_size(std::vector<Permutation<size_t>*>* const gens) {
    std::unordered_set<size_t> orbset;
    deg = ((*gens)[0])->degree();
    std::vector<size_t>                orb;
    size_t                             out = 1;
    std::vector<Permutation<size_t>*>* gens1;
    std::vector<Permutation<size_t>*>* gens2 = really_copy_vec(gens);
    base                                     = {};
    size_t i;

    // No need to check for i = deg - 1, since a permutation that fixes
    // everything  except one point must also fix that point.
    for (i = 0; (i <= (deg - 1)); i++) {
      gens1 = really_copy_vec(gens2);
      point_stabilizer(gens1, i, gens2);
      orb    = {i};
      orbset = {i};

      // Cannot iterate through orbit the normal way as we are adding to it.
      for (size_t j = 0; j < orb.size(); j++) {
        for (Permutation<size_t>* perm : *gens1) {
          if (orbset.find((*perm)[orb[j]]) == orbset.end()) {
            orbset.insert((*perm)[orb[j]]);
            orb.push_back((*perm)[orb[j]]);
          }
        }
      }
      out = out * orbset.size();
      really_delete_cont(gens1);
      delete gens1;
    }
    really_delete_cont(gens2);
    delete gens2;
    for (i = 0; i < base.size(); i++) {
      base[i] = 0;
    }
    free_stab_chain();
    return out;
  }
}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_SRC_SIMS_H_
