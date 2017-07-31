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

TEST_CASE("Sims 01: ", "[quick][sims][01][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({0, 2, 1}));
  gens->push_back(new Permutation<size_t>({1, 2, 0}));

  REQUIRE(group_size(gens) == 6);

  really_delete_cont(gens);
  delete gens;
}

TEST_CASE("Sims 02: ", "[quick][sims][02][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({0, 2, 1, 3}));
  gens->push_back(new Permutation<size_t>({1, 2, 3, 0}));

  REQUIRE(group_size(gens) == 24);

  really_delete_cont(gens);
  delete gens;
}

TEST_CASE("Sims 03: ", "[quick][sims][03][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({0, 2, 1, 3, 4}));
  gens->push_back(new Permutation<size_t>({1, 2, 3, 4, 0}));

  REQUIRE(group_size(gens) == 120);

  really_delete_cont(gens);
}
TEST_CASE("Sims 04: ", "[quick][sims][04][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({0, 4, 3, 2, 1}));
  gens->push_back(new Permutation<size_t>({1, 2, 3, 4, 0}));

  REQUIRE(group_size(gens) == 10);

  really_delete_cont(gens);
}

TEST_CASE("Sims 05: ", "[quick][sims][05][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({1, 0, 2, 3, 4, 5, 6, 7}));
  gens->push_back(new Permutation<size_t>({1, 2, 3, 4, 5, 6, 7, 0}));

  REQUIRE(group_size(gens) == 40320);

  really_delete_cont(gens);
}

TEST_CASE("Sims 06: ", "[quick][sims][06][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({1, 0, 2, 3, 4, 5, 6, 7}));
  gens->push_back(new Permutation<size_t>({0, 1, 3, 2, 4, 5, 6, 7}));
  gens->push_back(new Permutation<size_t>({0, 1, 2, 3, 4, 5, 7, 6}));

  REQUIRE(group_size(gens) == 8);

  really_delete_cont(gens);
}

TEST_CASE("Sims 07: ", "[quick][sims][07][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({1, 0, 2, 3, 4, 5, 6, 7}));
  gens->push_back(new Permutation<size_t>({0, 1, 3, 2, 4, 5, 6, 7}));
  gens->push_back(new Permutation<size_t>({0, 1, 2, 3, 5, 4, 6, 7}));
  gens->push_back(new Permutation<size_t>({0, 1, 2, 3, 4, 5, 7, 6}));

  REQUIRE(group_size(gens) == 16);

  really_delete_cont(gens);
}

TEST_CASE("Sims 08: ", "[quick][sims][08][group_size]") {

  std::vector<Permutation<size_t>*>* gens
      = new std::vector<Permutation<size_t>*>({});
  gens->push_back(new Permutation<size_t>({1, 0, 2, 3, 4, 5, 6, 7, 8}));
  gens->push_back(new Permutation<size_t>({1, 2, 3, 4, 5, 6, 7, 8, 0}));

  REQUIRE(group_size(gens) == 362880);

  really_delete_cont(gens);
}
