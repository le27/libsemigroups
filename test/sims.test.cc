//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2016 James D. Mitchell
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

#include "../src/semigroups.h"
#include "../src/sims.h"
#include "catch.hpp"

#define SIMS_REPORT false

using namespace libsemigroups;

template <typename T> static inline void really_delete_cont(T cont) {
  for (Element* x : cont) {
    x->really_delete();
    delete x;
  }
}

// typedef  Orb<Permutation<u_int16_t>, u_int16_t> OrbPermInt;
//
// TEST_CASE("Sims 01: ", "[quick][sims][01][orb]") {
//   std::vector<Permutation<u_int16_t>*> gens
//       = {new Permutation<u_int16_t>({1, 0, 2}),
//          new Permutation<u_int16_t>({1, 2, 0})};
//   // Semigroup S = Semigroup(gens);
//   // S.set_report(SIMS_REPORT);
//   // REQUIRE(S.size() == 6);
//
//   Orb<Permutation<u_int16_t>, u_int16_t> o
//       = Orb<Permutation<u_int16_t>, u_int16_t>(gens, 0);
//
//   REQUIRE(o.size() == 1);
//   o.enumerate();
//   REQUIRE(o.size() == 3);
//   really_delete_cont(gens);
// }
//
// TEST_CASE("Sims 02: ", "[quick][sims][02][orb]") {
//   std::vector<u_int16_t> p;
//   for (size_t i = 0; i < 32768; i++) {
//     p.push_back(i + 1);
//   }
//   p.push_back(0);
//
//   std::vector<Permutation<u_int16_t>*> gens
//       = {new Permutation<u_int16_t>(p)};
//
//   Orb<Permutation<u_int16_t>, u_int16_t> o
//       = Orb<Permutation<u_int16_t>, u_int16_t>(gens, 2);
//
//   REQUIRE(o.size() == 1);
//   o.enumerate();
//   REQUIRE(o.size() == 32769);
//   really_delete_cont(gens);
// }
//
// TEST_CASE("Sims 03: ", "[quick][sims][03][orb]") {
//   std::vector<u_int32_t>* p = new std::vector<u_int32_t>();
//   p->reserve(32768001);
//   for (size_t i = 0; i < 32768000; i++) {
//     p->push_back(i + 1);
//   }
//   p->push_back(0);
//
//   std::vector<Permutation<u_int32_t>*> gens
//       = {new Permutation<u_int32_t>(p)};
//
//   Orb<Permutation<u_int32_t>, u_int32_t> o
//       = Orb<Permutation<u_int32_t>, u_int32_t>(gens, 29);
//
//   REQUIRE(o.size() == 1);
//   o.enumerate();
//   REQUIRE(o.size() == 32768001);
//   really_delete_cont(gens);
// }
//
// TEST_CASE("Sims 04: ", "[quick][sims][04][orb]") {
//   std::vector<u_int32_t>* p = new std::vector<u_int32_t>();
//   p->reserve(327681);
//   for (size_t i = 0; i < 327680; i++) {
//     p->push_back(i + 1);
//   }
//   p->push_back(0);
//
//   std::vector<Permutation<u_int32_t>*> gens
//       = {new Permutation<u_int32_t>(p)};
//
//   Orb<Permutation<u_int32_t>, u_int32_t> o
//       = Orb<Permutation<u_int32_t>, u_int32_t>(gens, 29);
//
//   REQUIRE(o.size() == 1);
//   o.enumerate();
//   REQUIRE(o.position(29) == 0);
//   REQUIRE(o.position(30) == 1);
//   REQUIRE(o.position(327681) == static_cast<size_t>(-1));
//   really_delete_cont(gens);
// }

// TEST_CASE("Sims 05: ", "[quick][sims][05][BSGS]") {
//   std::vector<Permutation<u_int16_t>*> gens
//       = {new Permutation<u_int16_t>({1, 0, 2}),
//          new Permutation<u_int16_t>({1, 2, 0})};
//   // Semigroup S = Semigroup(gens);
//   // S.set_report(SIMS_REPORT);
//   // REQUIRE(S.size() == 6);
//
//
//   BSGS<Permutation<u_int16_t>, u_int16_t> b
//       = BSGS<Permutation<u_int16_t>, u_int16_t>(gens);
//
//   b.enumerate_partial();
//   REQUIRE(b.base[0] == 1);
//   REQUIRE(b.strong_gen_set.size() == 3);
//
//
//   really_delete_cont(gens);
// }
//
// TEST_CASE("Sims 06: ", "[quick][sims][06][BSGS]") {
//   std::vector<Permutation<u_int16_t>*> gens
//       = {new Permutation<u_int16_t>({1, 0, 2}),
//          new Permutation<u_int16_t>({1, 2, 0})};
//   // Semigroup S = Semigroup(gens);
//   // S.set_report(SIMS_REPORT);
//   // REQUIRE(S.size() == 6);
//
//
//   BSGS<Permutation<u_int16_t>, u_int16_t> b
//       = BSGS<Permutation<u_int16_t>, u_int16_t>(gens);
//
//   b.enumerate();
//   REQUIRE(b.base[0] == 1);
//   REQUIRE(b.base[1] == 2);
//   REQUIRE(b.base.size() == 2);
//   REQUIRE(b.strong_gen_set.size() == 3);
//
//
//   really_delete_cont(gens);
// }

TEST_CASE("Sims 01: ", "[quick][sims][01][orb]") {
  // Semigroup S = Semigroup(gens);
  // S.set_report(SIMS_REPORT);
  // REQUIRE(S.size() == 6);

  PermColl* gens = new_perm_coll(3, 3);

  REQUIRE(group_size(gens) == 6);
  gens->really_delete();
}
