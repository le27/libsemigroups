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

// This file contains the declaration of the element class and its subclasses.

#ifndef LIBSEMIGROUPS_SRC_ELEMENTS_H_
#define LIBSEMIGROUPS_SRC_ELEMENTS_H_

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "blocks.h"
#include "recvec.h"
#include "semiring.h"

#define PP_UNDEFINED PartialTransformation<T, PartialPerm<T>>::UNDEFINED

namespace libsemigroups {

  //! Abstract base class for semigroup elements
  //!
  //! The Semigroup class consists of Element objects. Every derived class of
  //! Element implements the deleted methods of Element, which are used by the
  //! methods of the Semigroup class.
  class Element {
   public:
    //! This enum contains some different types of Element.
    //!
    //! This exists so that the type of a subclass of Element can be determined
    //! from a pointer to the base class. Currently, it is only necessary
    //! within **libsemigroups** to distiguish RWSE objects from other Element
    //! objects and so there are only two values in here.
    enum elm_t {
      //! Type for Element objects arising from a rewriting system RWS.
      RWSE = 0,
      //! Type for Element objects not arising from a rewriting system RWS.
      NOT_RWSE = 1
    };

    //! A constructor.
    //!
    //! The parameter \p hv is the hash value (for caching) of the element being
    //! created (defaults to libsemigroups::Element::UNDEFINED). This value
    //! should only be set if it is known *a priori* that \p hv is the correct
    //! value (such as when copying an Element object).
    //!
    //! The parameter \p type should be the type elm_t of the element being
    //! created (defaults to libsemigroups::Element::NOT_RWSE).
    explicit Element(size_t hv   = Element::UNDEFINED,
                     elm_t  type = Element::elm_t::NOT_RWSE)
        : _hash_value(hv), _type(type) {}

    //! A default destructor.
    //!
    //! This does not properly delete the underlying data of the object, this
    //! should be done using Element::really_delete.
    virtual ~Element() {}

    //! Returns the type libsemigroups::Element::elm_t of an Element object.
    //!
    //! There are currently only two types libsemigroups::Element::RWSE
    //! and libsemigroups::Element::NOT_RWSE distinguished by
    //! the return value of this method.
    elm_t get_type() const {
      return _type;
    }

    //! Returns \c true if \c this equals \p that.
    //!
    //! This method checks the mathematical equality of two Element objects in
    //! the same subclass of Element.
    virtual bool operator==(const Element& that) const = 0;

    //! Returns \c true if \c this is less than \p that.
    //!
    //! This method defines a total order on the set of objects in a given
    //! subclass of Element with a given Element::degree. The definition of
    //! this total order depends on the method for the operator < in the
    //! subclass.
    virtual bool operator<(const Element& that) const = 0;

    //!  Returns the approximate time complexity of multiplying two
    //! Elements in a given subclass.
    //!
    //! This method returns an integer which represents the approximate time
    //! complexity of multiplying two objects in the same subclass of Element
    //! which have the same Element::degree. For example, the approximate time
    //! complexity of multiplying two \f$3\times 3\f$ matrices over a common
    //! semiring is \f$O(3 ^ 3)\f$, and 27 is returned by
    //! MatrixOverSemiring::complexity.
    //!
    //! The returned value is used in, for example, Semigroup::fast_product and
    //! Semigroup::nridempotents to decide if it is better to multiply
    //! elements or follow a path in the Cayley graph.
    virtual size_t complexity() const = 0;

    //! Returns the degree of an Element.
    //!
    //! This method returns an integer which represents the size of the
    //! element, and is used to determine whether or not two elements are
    //! compatible for multiplication. For example, two Transformation objects
    //! of different degrees cannot be multiplied, or a Bipartition of degree
    //! 10 cannot be an element of a monoid of bipartitions of degree 3.
    //!
    //! See the relevant subclass for the particular meaning of the return value
    //! of this method for each subclass.
    virtual size_t degree() const = 0;

    //! Return the hash value of an Element.
    //!
    //! This method returns a hash value for an object in a subclass of
    //! Element. This value is only computed the first time this method is
    //! called.
    inline size_t hash_value() const {
      if (_hash_value == UNDEFINED) {
        this->cache_hash_value();
      }
      return this->_hash_value;
    }

    //! Returns the identity element.
    //!
    //! This method returns the multiplicative identity element for an object
    //! in a subclass of Element. The returned identity belongs to the same
    //! subclass and has the same degree as \c this.
    virtual Element* identity() const = 0;

    //! Returns a new element completely independent of \c this.
    //!
    //! If \p increase_deg_by is non-zero, then the degree (the size of the
    //! container) of the defining data of \c this will be increased by this
    //! amount.
    //!
    //! This method really copies an Element. To minimise the amount of copying
    //! when Elements are inserted in a std::unordered_map and other
    //! containers, an Element behaves somewhat like a pointer, in that the
    //! actual data in an Element is only copied when this method is called.
    //! Otherwise, if an Element is copied, then its defining data is only
    //! stored once.
    //!
    //! \sa Element::really_delete.
    virtual Element* really_copy(size_t increase_deg_by = 0) const = 0;

    //! Copy another Element into \c this.
    //!
    //! This method copies \p x into \c this by changing \c this in-place.
    virtual void copy(Element const* x) = 0;

    //! Deletes the defining data of an Element.
    //!
    //! This method really deletes an Element. To minimise the amount of
    //! copying when Elements are inserted in an std::unordered_map or other
    //! containers, an Element behaves somewhat like a pointer, in that the
    //! actual data in an Element is only deleted when this method is called.
    //! Otherwise, if an Element is deleted, then its defining data is not
    //! deleted, since it may be contained in several copies of the Element.
    //!
    //! \sa Element::really_copy.
    virtual void really_delete() = 0;

    //! Multiplies \p x and \p y and stores the result in \c this.
    //!
    //! Redefine \c this to be the product of \p x and \p y. This is in-place
    //! multiplication to avoid allocation of memory for products which do not
    //! need to be stored for future use.
    //!
    //! The implementation of this method in the Element base class simply
    //! calls the 3 parameter version with third parameter 0. Any subclass of
    //! Element can implement either a two or three parameter version of this
    //! method and the base class method implements the other method.
    virtual void redefine(Element const* x, Element const* y) {
      redefine(x, y, 0);
    }

    //! Multiplies \p x and \p y and stores the result in \c this.
    //!
    //! Redefine \c this to be the product of \p x and \p y. This is in-place
    //! multiplication to avoid allocation of memory for products which do not
    //! need to be stored for future use.
    //!
    //! The implementation of this method in the Element base class simply
    //! calls the 2 parameter version and ignores the third parameter \p
    //! thread_id. Any subclass of Element can implement either a two or three
    //! parameter version of this method and the base class method implements
    //! the other method.
    //!
    //! The parameter \p thread_id is required in some derived classes of
    //! Element because some temporary storage is required to find the product
    //! of \p x and \p y.
    //!
    //! Note that if different threads call this method, on a derived class of
    //! Element where static temporary storage is used in the redefine method
    //! with the same value of \p thread_id, then bad things may happen.
    virtual void
    redefine(Element const* x, Element const* y, size_t const& thread_id) {
      (void) thread_id;
      redefine(x, y);
    }

    //! Provides a call operator for comparing Elements via pointers.
    //!
    //! This struct provides a call operator for comparing const Element
    //! pointers (by comparing the Elements they point to). This is used by
    //! various methods of the Semigroup class.
    struct Equal {
      //! Returns \c true if \p x and \p y point to equal Element's.
      size_t operator()(Element const* x, Element const* y) const {
        return *x == *y;
      }
    };

    //!  Provides a call operator returning a hash value for an Element
    //! via a pointer.
    //!
    //! This struct provides a call operator for obtaining a hash value for the
    //! Element from a const Element pointer. This is used by
    //! various methods of the Semigroup class.
    struct Hash {
      //!  Returns the value of Element::hash_value applied to the
      //! Element pointed to by \p x.
      size_t operator()(Element const* x) const {
        return x->hash_value();
      }
    };

   protected:
    //! Calculate and cache a hash value.
    //!
    //! This method is used to compute and cache the hash value of \c this.
    virtual void cache_hash_value() const = 0;

    //! Reset the cached value used by Element::hash_value.
    //!
    //! This method is used to reset the cached hash value to
    //! libsemigroups::Element::UNDEFINED. This is required after running
    //! Element::redefine, Element::copy, or any other method that changes the
    //! defining data of \c this.
    void reset_hash_value() const {
      _hash_value = UNDEFINED;
    }

    //! UNDEFINED value.
    //!
    //! This variable is used to indicate that a value is undefined, such as,
    //! the cached hash value.
    static size_t const UNDEFINED;

    //! This data member holds a cached version of the hash value of an Element.
    //! It is stored here if it is ever computed. It is invalidated by
    //! libsemigroups::Element::redefine and sometimes by
    //! libsemigroups::Element::really_copy.
    mutable size_t _hash_value;

   private:
    elm_t _type;
  };

  //!  Abstract base class for elements using a vector to store their
  //! defining data.
  //!
  //! The template parameter \p S is the type entries in the vector containing
  //! the defining data.  The value of the template parameter \p S can be used
  //! to reduce (or increase) the amount of memory required by instances of
  //! this class.
  //!
  //! The template parameter \p T is the subclass of ElementWithVectorData used
  //! by certain methods to construct new instances of subclasses of
  //! ElementWithVectorData.
  //!
  //! For example, Transformation&lt;u_int128_t&gt; is a subclass of
  //! ElementWithVectorData&lt;u_int128_t,
  //! Transformation&lt;u_int128_t&gt;&gt;.

  template <typename S, class T> class ElementWithVectorData : public Element {
   public:
    //! A constructor.
    //!
    //! The parameter \p size is used to reserve the size of the underlying
    //! std::vector.
    //!
    //! Returns an object with an uninitialised vector which has size
    //! \p size (defaults to 0).
    explicit ElementWithVectorData(size_t size = 0)
        : Element(), _vector(new std::vector<S>()) {
      _vector->resize(size);
    }

    //! A constructor.
    //!
    //! The parameter \p vector should be a pointer to defining data of the
    //! element The parameter \p hv must be the hash value of the element being
    //! created (this defaults to Element::UNDEFINED). This should only be set
    //! if it is guaranteed that \p hv is the correct value.
    //!
    //! Returns an object whose defining data is stored in \p vector, which
    //! is not copied, and should be deleted using Element::really_delete.
    explicit ElementWithVectorData(std::vector<S>* vector,
                                   size_t          hv = Element::UNDEFINED)
        : Element(hv), _vector(vector) {}

    //! A constructor.
    //!
    //! The parameter \p vector should be a const reference to defining data of
    //! the element.
    //!
    //! Returns an object whose defining data is stored in \p vector, which
    //! is copied.
    explicit ElementWithVectorData(std::vector<S> const& vector)
        : ElementWithVectorData(new std::vector<S>(vector)) {}

    //! Returns the \p pos entry in the vector containing the defining data.
    //!
    //! This method returns the \p pos entry in the vector used to construct \c
    //! this. No checks are performed that \p pos in within the bounds of this
    //! vector.
    inline S operator[](size_t pos) const {
      return (*_vector)[pos];
    }

    //! Returns the \p pos entry in the vector containing the defining data.
    //!
    //! This method returns the \p pos entry in the vector used to construct \c
    //! this.
    inline S at(size_t pos) const {
      return _vector->at(pos);
    }

    //! Returns \c true if \c this equals \p that.
    //!
    //! This method checks that the underlying vectors of \c this and \p that
    //! are equal.
    bool operator==(Element const& that) const override {
      return *(static_cast<T const&>(that)._vector) == *(this->_vector);
    }

    //! Returns \c true if \c this is less than \p that.
    //!
    //! This method defines a total order on the set of objects in
    //! ElementWithVectorData of a given Element::degree, which is the
    //! short-lex order.
    bool operator<(const Element& that) const override {
      T const& ewvd = static_cast<T const&>(that);
      if (this->_vector->size() != ewvd._vector->size()) {
        return this->_vector->size() < ewvd._vector->size();
      }
      for (size_t i = 0; i < this->_vector->size(); i++) {
        if ((*this)[i] != ewvd[i]) {
          return (*this)[i] < ewvd[i];
        }
      }
      return false;
    }

    //! Returns a pointer to a copy of \c this.
    //!
    //! The size of the vector containing the defining data of \c this will be
    //! increased by \p increase_deg_by.  If \p increase_deg_by is not 0, then
    //! this method must be overridden by any subclass of ElementWithVectorData
    //! since there is no way of knowing how a subclass is defined by the data
    //! in the vector.
    Element* really_copy(size_t increase_deg_by) const override {
      assert(increase_deg_by == 0);
      (void) increase_deg_by;
      std::vector<S>* vector(new std::vector<S>(*_vector));
      return new T(vector, this->_hash_value);
    }

    //! Copy another Element into \c this.
    //!
    //! This method copies \p x into \c this by changing \c this in-place. This
    //! method asserts that the degrees of \c this and \p x are equal and then
    //! replaces the underlying vector of \c this with the underlying vector of
    //! \p x. Any method overriding this one must call
    //! Element::reset_hash_value.
    void copy(Element const* x) override {
      assert(x->degree() == this->degree());
      auto   xx  = static_cast<ElementWithVectorData const*>(x);
      size_t deg = _vector->size();
      for (size_t i = 0; i < deg; i++) {
        (*_vector)[i] = (*xx)[i];
      }
      this->reset_hash_value();
    }

    //! Deletes the defining data of an ElementWithVectorData.
    void really_delete() override {
      delete _vector;
    }

    //! Returns an iterator.
    //!
    //! This method returns an iterator pointing at the first entry in the
    //! vector that is the underlying defining data of \c this.
    inline typename std::vector<S>::iterator begin() const {
      return _vector->begin();
    }

    //! Returns an iterator.
    //!
    //! This method returns an iterator referring to the past-the-end element
    //! of the vector that is the underlying defining data of \c this.
    inline typename std::vector<S>::iterator end() const {
      return _vector->end();
    }

    //! Returns a const iterator.
    //!
    //! This method returns a const_iterator pointing at the first entry in the
    //! vector that is the underlying defining data of \c this.
    inline typename std::vector<S>::iterator cbegin() const {
      return _vector->cbegin();
    }

    //! Returns a const iterator.
    //!
    //! This method returns a const iterator referring to the past-the-end
    //! element of the vector that is the underlying defining data of \c this.
    inline typename std::vector<S>::iterator cend() const {
      return _vector->cend();
    }

   protected:
    //! The vector containing the defining data of \c this.
    //!
    //! The actual data defining of \c this is stored in _vector.
    std::vector<S>* _vector;
  };

  //! Abstract class for partial transformations.
  //!
  //! This is a template class for partial transformations, which is a subclass
  //! of ElementWithVectorData. For example, Transformation<u_int128_t> is a
  //! subclass of PartialTransformation<u_int128_t,
  //! Transformation<u_int128_t>>.
  //!
  //! The template parameter \p S is the type of image values, i.e. u_int16_t,
  //! and so on.  The value of the template parameter \p S can be used to
  //! reduce (or increase) the amount of memory required by instances of this
  //! class.
  //!
  //! The template parameter \p T is the subclass of PartialTransformation used
  //! by the PartialTransformation::identity method to construct an identity.
  //!
  //! This class is abstract since it does not implement all methods required by
  //! the Element class, it exists to provide common methods for its
  //! subclasses.
  //!
  //! A *partial transformation* \f$f\f$ is just a function defined on a subset
  //! of \f$\{0, 1, \ldots, n - 1\}\f$ for some integer \f$n\f$ called the
  //! *degree*  of *f*.
  //! A partial transformation is stored as a vector of the images of
  //! \f$\{0, 1, \ldots, n -1\}\f$, i.e. \f$\{(0)f, (1)f, \ldots, (n - 1)f\}\f$
  //! where the value PartialTransformation::UNDEFINED is used to indicate that
  //! \f$(i)f\f$ is, you guessed it, undefined (i.e. not among the points where
  //! \f$f\f$ is defined).

  template <typename S, typename T>
  class PartialTransformation : public ElementWithVectorData<S, T> {
   public:
    //! A constructor.
    //!
    //! Constructs a partial transformation with no data at all.
    PartialTransformation() : ElementWithVectorData<S, T>() {}

    //! A constructor.
    //!
    //! Constructs a partial transformation with list of images equal to
    //! \p vector, which is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p hv must be the hash value of the element
    //! being created (this defaults to Element::UNDEFINED). This should only
    //! be set if it is guaranteed that \p hv is the correct value. See
    //! Element::Element for more details.
    explicit PartialTransformation(std::vector<S>* vector,
                                   size_t          hv = Element::UNDEFINED)
        : ElementWithVectorData<S, T>(vector, hv) {}

    //! A constructor.
    //!
    //! Constructs a partial transformation with list of images equal to
    //! \p vector, which is copied into the constructed object.
    explicit PartialTransformation(std::vector<S> const& vector)
        : ElementWithVectorData<S, T>(vector) {}

    //! Returns the approximate time complexity of multiplying two
    //! partial transformations.
    //!
    //! The approximate time complexity of multiplying partial transformations
    //! is just their degree.
    size_t complexity() const override {
      return this->_vector->size();
    }

    //! Returns the degree of a partial transformation.
    //!
    //! The *degree* of a partial transformation is the number of points used
    //! in its definition, which is equal to the length of
    //! ElementWithVectorData::_vector.
    size_t degree() const override {
      return this->_vector->size();
    }

    //! Returns the rank of a partial transformation.
    //!
    //! The *rank* of a partial transformation is the number of its distinct
    //! image values, not including PartialTransformation::UNDEFINED. This
    //! method recomputes the return value every time it is called.
    size_t crank() const {
      _lookup.clear();
      _lookup.resize(degree(), false);
      size_t r = 0;
      for (auto const& x : *(this->_vector)) {
        if (x != UNDEFINED && !_lookup[x]) {
          _lookup[x] = true;
          r++;
        }
      }
      return r;
    }

    //! Find the hash value of a partial transformation.
    //!
    //! \sa Element::hash_value for more details.
    void cache_hash_value() const override {
      size_t seed = 0;
      size_t deg  = this->degree();
      for (auto const& val : *(this->_vector)) {
        seed *= deg;
        seed += val;
      }
      this->_hash_value = seed;
    }

    //! Returns the identity transformation with degrees of \c this.
    //!
    //! This method returns a new partial transformation with degree equal to
    //! the degree of \c this that fixes every value from *0* to the degree of
    //! \c this.
    Element* identity() const override {
      std::vector<S>* vector(new std::vector<S>());
      vector->reserve(this->degree());
      for (size_t i = 0; i < this->degree(); i++) {
        vector->push_back(i);
      }
      return new T(vector);
    }

    //! Undefined image value.
    //!
    //! This value is used to indicate that a partial transformation is not
    //! defined on a value.
    static S const UNDEFINED;

   private:
    // Used for determining rank, this is not thread safe FIXME
    static std::vector<bool> _lookup;
  };

  template <typename S, typename T>
  std::vector<bool> PartialTransformation<S, T>::_lookup = std::vector<bool>();

  template <typename S, typename T>
  S const PartialTransformation<S, T>::UNDEFINED
      = std::numeric_limits<S>::max();

  //! Template class for transformations.
  //!
  //! The value of the template parameter \p T can be used to reduce the amount
  //! of memory required by instances of this class; see PartialTransformation
  //! and ElementWithVectorData for more details.
  //!
  //! A *transformation* \f$f\f$ is just a function defined on the whole of
  //! \f$\{0, 1, \ldots, n - 1\}\f$ for some integer \f$n\f$ called the
  //! *degree* of \f$f\f$.
  //! A transformation is stored as a vector of the images of
  //! \f$\{0, 1, \ldots, n - 1\}\f$,
  //! i.e. \f$\{(0)f, (1)f, \ldots, (n - 1)f\}\f$.
  template <typename T>
  class Transformation : public PartialTransformation<T, Transformation<T>> {
   public:
    //! A constructor.
    //!
    //! Constructs a transformation with list of images equal to \p vector,
    //! \p vector is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p hv must be the hash value of the element
    //! being created (this defaults to Element::UNDEFINED). This should only
    //! be set if it is guaranteed that \p hv is the correct value. See
    //! Element::Element for more details.
    explicit Transformation(std::vector<T>* vector,
                            size_t          hv = Element::UNDEFINED)
        : PartialTransformation<T, Transformation<T>>(vector, hv) {}

    //! A constructor.
    //!
    //! Constructs a transformation with list of images equal to
    //! \p vector, which is copied into the constructed object.
    explicit Transformation(std::vector<T> const& vector)
        : PartialTransformation<T, Transformation<T>>(vector) {}

    //! Returns a pointer to a copy of \c this.
    //!
    //! See Element::really_copy for more details about this method.
    //!
    //! The copy returned by this method fixes all the values between the
    //! Transformation::degree of \c this and \p increase_deg_by.
    Element* really_copy(size_t increase_deg_by = 0) const override {
      std::vector<T>* vector_copy = new std::vector<T>(*this->_vector);
      if (increase_deg_by == 0) {
        return new Transformation<T>(vector_copy, this->_hash_value);
      } else {
        size_t n = vector_copy->size();
        vector_copy->reserve(n + increase_deg_by);
        for (size_t i = n; i < n + increase_deg_by; i++) {
          vector_copy->push_back(i);
        }
        return new Transformation<T>(vector_copy);
      }
    }

    //! Multiply \p x and \p y and stores the result in \c this.
    //!
    //! See Element::redefine for more details about this method.
    //!
    //! This method asserts that the degrees of \p x, \p y, and \c this, are
    //! all equal, and that neither \p x nor \p y equals \c this.
    void redefine(Element const* x, Element const* y) override {
      assert(x->degree() == y->degree());
      assert(x->degree() == this->degree());
      assert(x != this && y != this);
      Transformation<T> const* xx(static_cast<Transformation<T> const*>(x));
      Transformation<T> const* yy(static_cast<Transformation<T> const*>(y));
      size_t const             n = this->_vector->size();
      for (T i = 0; i < n; i++) {
        (*this->_vector)[i] = (*yy)[(*xx)[i]];
      }
      this->reset_hash_value();
    }
  };

template <typename T>
  class Permutation : public Transformation<T> {
   public:
    explicit Permutation(std::vector<T>* vector,
                            size_t          hv = Element::UNDEFINED)
        : Transformation<T>(vector, hv) {}

    explicit Permutation(std::vector<T> const& vector)
        : Transformation<T>(vector) {}

   Permutation* inverse() {
   size_t const             n = this->_vector->size();
   Permutation* id = static_cast<Permutation<T>*>(this->identity());
      for (T i = 0; i < n; i++) {
        (*id->_vector)[(*this->_vector)[i]] = i;
      }
    return id;
}
  };

  //! Template class for partial permutations.
  //!
  //! The value of the template parameter \p T can be used to reduce the amount
  //! of memory required by instances of this class; see PartialTransformation
  //! and ElementWithVectorData for more details.
  //!
  //! A *partial permutation* \f$f\f$ is just an injective partial
  //! transformation, which is stored as a vector of the images of
  //! \f$\{0, 1, \ldots, n - 1\}\f$, i.e.
  //! i.e. \f$\{(0)f, (1)f, \ldots, (n - 1)f\}\f$ where the value
  //! PartialTransformation::UNDEFINED is
  //! used to indicate that \f$(i)f\f$ is undefined (i.e. not among
  //! the points where \f$f\f$ is defined).
  template <typename T>
  class PartialPerm : public PartialTransformation<T, PartialPerm<T>> {
   public:
    //! A constructor.
    //!
    //! Constructs a partial perm with list of images equal to vector,
    //! vector is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p hv must be the hash value of the element
    //! being created (this defaults to Element::UNDEFINED). This should only
    //! be set if it is guaranteed that \p hv is the correct value. See
    //! Element::Element for more details.
    explicit PartialPerm(std::vector<T>* vector, size_t hv = Element::UNDEFINED)
        : PartialTransformation<T, PartialPerm<T>>(vector, hv) {}

    //! A constructor.
    //!
    //! Constructs a partial perm with list of images equal to vector,
    //! which is copied into the constructed object.
    explicit PartialPerm(std::vector<T> const& vector)
        : PartialTransformation<T, PartialPerm<T>>(vector) {}

    //! A constructor.
    //!
    //! Constructs a partial perm of degree \p deg such that
    //! \code (dom[i])f = ran[i] \endcode
    //! for all \c i and which is undefined on every other value in the range 0
    //! to (strictly less than \p deg). This method asserts that \p dom and \p
    //! ran have equal size and that \p deg is greater than or equal to the
    //! maximum value in \p dom or \p ran.
    explicit PartialPerm(std::vector<T> const& dom,
                         std::vector<T> const& ran,
                         size_t                deg)
        : PartialTransformation<T, PartialPerm<T>>() {
      assert(dom.size() == ran.size());
      assert(dom.empty() || deg >= *std::max_element(dom.begin(), dom.end()));

      this->_vector->resize(deg + 1, PP_UNDEFINED);
      for (size_t i = 0; i < dom.size(); i++) {
        (*this->_vector)[dom[i]] = ran[i];
      }
    }

    //! Returns \c true if \c this is less than \p that.
    //!
    //! This defines a total order on partial permutations that is equivalent to
    //! that used by GAP. It is not short-lex on the list of images.
    //!
    //! Returns \c true if something complicated is \c true and \c false if
    //! it is not.
    bool operator<(const Element& that) const override {
      auto pp_this = static_cast<const PartialPerm<T>*>(this);
      auto pp_that = static_cast<const PartialPerm<T>&>(that);

      size_t deg_this = pp_this->degree();
      for (auto it = pp_this->_vector->end() - 1;
           it >= pp_this->_vector->begin();
           it--) {
        if (*it == PP_UNDEFINED) {
          deg_this--;
        } else {
          break;
        }
      }
      size_t deg_that = pp_that.degree();
      for (auto it = pp_that._vector->end() - 1;
           it >= pp_that._vector->begin() && deg_that >= deg_this;
           it--) {
        if (*it == PP_UNDEFINED) {
          deg_that--;
        } else {
          break;
        }
      }

      if (deg_this != deg_that) {
        return deg_this < deg_that;
      }

      for (size_t i = 0; i < deg_this; i++) {
        if ((*pp_this)[i] != pp_that[i]) {
          return (*pp_this)[i] == PP_UNDEFINED
                 || (pp_that[i] != PP_UNDEFINED && (*pp_this)[i] < pp_that[i]);
        }
      }
      return false;
    }

    //! Returns a pointer to a copy of \c this.
    //!
    //! See Element::really_copy for more details about this method.
    //!
    //! The copy returned by this method is undefined on all the values between
    //! the PartialPerm::degree of \c this and \p increase_deg_by.
    Element* really_copy(size_t increase_deg_by = 0) const override {
      std::vector<T>* vector_copy = new std::vector<T>(*this->_vector);
      if (increase_deg_by == 0) {
        return new PartialPerm<T>(vector_copy, this->_hash_value);
      } else {
        size_t n = vector_copy->size();
        vector_copy->reserve(n + increase_deg_by);
        for (size_t i = n; i < n + increase_deg_by; i++) {
          vector_copy->push_back(PP_UNDEFINED);
        }
        return new PartialPerm<T>(vector_copy);
      }
    }

    //! Multiply \p x and \p y and stores the result in \c this.
    //!
    //! See Element::redefine for more details about this method.
    //!
    //! This method asserts that the degrees of \p x, \p y, and \c this, are
    //! all equal, and that neither \p x nor \p y equals \c this.
    void redefine(Element const* x, Element const* y) override {
      assert(x->degree() == y->degree());
      assert(x->degree() == this->degree());
      assert(x != this && y != this);
      PartialPerm<T> const* xx(static_cast<PartialPerm<T> const*>(x));
      PartialPerm<T> const* yy(static_cast<PartialPerm<T> const*>(y));
      size_t const          n = this->degree();
      for (T i = 0; i < n; i++) {
        (*this->_vector)[i]
            = ((*xx)[i] == PP_UNDEFINED ? PP_UNDEFINED : (*yy)[(*xx)[i]]);
      }
      this->reset_hash_value();
    }

    //! Returns the rank of a partial permutation.
    //!
    //! The *rank* of a partial permutation is the number of its distinct image
    //! values, not including PartialTransformation::UNDEFINED. This method
    //! involves slightly less work than PartialTransformation::crank since a
    //! partial permutation is injective, and so every image value occurs
    //! precisely once. This method recomputes the return value every time it
    //! is called.
    size_t crank() const {
      size_t nr_defined = 0;
      for (auto const& x : *this->_vector) {
        if (x != PP_UNDEFINED) {
          nr_defined++;
        }
      }
      return nr_defined;
    }
  };

  //! Class for square boolean matrices.
  //!
  //! A *boolean matrix* is a square matrix over the Boolean semiring, under
  //! the usual multiplication of matrices.
  class BooleanMat : public ElementWithVectorData<bool, BooleanMat> {
    // FIXME why is this not a subclass of MatrixOverSemiring?
   public:
    //! A constructor.
    //!
    //! Constructs a Boolean matrix defined by \p matrix, \p matrix is not
    //! copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p matrix should be a vector of boolean values of length
    //! \f$n ^  2\f$ for some integer \f$n\f$, so that the value in position
    //! \f$ni + j\f$ is the entry in row \f$i\f$ and column \f$j\f$ of
    //! the constructed matrix.
    //!
    //! The parameter \p hv must be the hash value of the element
    //! being created (this defaults to Element::UNDEFINED). This should only
    //! be set if it is guaranteed that \p hv is the correct value. See
    //! Element::Element for more details.
    explicit BooleanMat(std::vector<bool>* matrix,
                        size_t             hv = Element::UNDEFINED)
        : ElementWithVectorData<bool, BooleanMat>(matrix, hv) {}

    //! A constructor.
    //!
    //! Constructs a boolean matrix defined by matrix, which is copied into the
    //! constructed boolean matrix; see BooleanMat::BooleanMat for more
    //! details.
    explicit BooleanMat(std::vector<std::vector<bool>> const& matrix);

    //!  Returns the approximate time complexity of multiplying two
    //! boolean matrices.
    //!
    //! See Element::complexity for more details.
    //!
    //! The approximate time complexity of multiplying boolean matrices is
    //! \f$n ^ 3\f$ where \f$n\f$ is the dimension of the matrix.
    size_t complexity() const override;

    //! Returns the dimension of the boolean matrix.
    //!
    //! The *dimension* of a boolean matrix is just the number of rows (or,
    //! equivalently columns).
    //!
    //! See Element::degree for more details.
    size_t degree() const override;

    //! Find the hash value of a boolean matrix.
    //!
    //! See Element::hash_value or Element::cache_hash_value for more details.
    void cache_hash_value() const override;

    //! Returns the identity boolean matrix with dimension of \c this.
    //!
    //! This method returns a new boolean matrix with dimension equal to that
    //! of \c this, where the main diagonal consists of the value \c true and
    //! every other entry is \c false.
    Element* identity() const override;

    //! Multiplies \p x and \p y and stores the result in \c this.
    //!
    //! This method asserts that the dimensions of \p x, \p y, and \c this, are
    //! all equal, and that neither \p x nor \p y equals \c this.
    void redefine(Element const* x, Element const* y) override;
  };

  //! Class for bipartitions.
  //!
  //! A *bipartition* is a partition of the set
  //! \f$\{0, ..., 2n - 1\}\f$ for some integer \f$n\f$; see [TODO](TODO) for
  //! more details.
  // FIXME Reference the Semigroups in GAP manual.

  //! The Bipartition class is more complex (i.e. has more methods) than
  //! strictly required by the algorithms for a Semigroup object because the
  //! extra methods are used in the GAP package [Semigroups package for
  //! GAP](https://gap-packages.github.io/Semigroups/).

  class Bipartition : public ElementWithVectorData<u_int32_t, Bipartition> {
    // TODO(JDM) add more explanation to the doc here
   public:
    //! A constructor.
    //!
    //! Constructs a uninitialised bipartition of degree \p degree.
    explicit Bipartition(size_t degree)
        : ElementWithVectorData<u_int32_t, Bipartition>(2 * degree),
          _nr_blocks(Bipartition::UNDEFINED),
          _nr_left_blocks(Bipartition::UNDEFINED),
          _trans_blocks_lookup(),
          _rank(Bipartition::UNDEFINED) {}

    //! A constructor.
    //!
    //! The parameter \p blocks must have length *2n* for some positive integer
    //! *n*, consist of non-negative integers, and have the property that if
    //! *i*, *i > 0*, occurs in \p blocks, then *i - 1* occurs earlier in
    //! blocks.  None of this is checked.
    //!
    //! The parameter \p blocks is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p hv can be the hash value of the element being created,
    //! if it is known (defaults to  Element::UNDEFINED).
    explicit Bipartition(std::vector<u_int32_t>* blocks,
                         size_t                  hv = Element::UNDEFINED)
        : ElementWithVectorData<u_int32_t, Bipartition>(blocks, hv),
          _nr_blocks(Bipartition::UNDEFINED),
          _nr_left_blocks(Bipartition::UNDEFINED),
          _trans_blocks_lookup(),
          _rank(Bipartition::UNDEFINED) {}

    //! A constructor.
    //!
    //! The parameter \p blocks must have length *2n* for some positive integer
    //! *n*, consist of non-negative integers, and have the property that if
    //! *i*, *i > 0*, occurs in \p blocks, then *i - 1* occurs earlier in
    //! blocks.  None of this is checked.
    //!
    //! The parameter \p blocks is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    explicit Bipartition(std::vector<u_int32_t> const& blocks)
        : Bipartition(new std::vector<u_int32_t>(blocks)) {}

    // TODO(JDM) another constructor that accepts an actual partition

    //! Returns the approximate time complexity of multiplication.
    //!
    //! In the case of a Bipartition of degree *n* the value *2n ^ 2* is
    //! returned.
    size_t complexity() const override;

    //! Returns the degree of the bipartition.
    //!
    //! A bipartition is of degree *n* if it is a partition of
    //! \f$\{0, \ldots, 2n -  1\}\f$.
    size_t degree() const override;

    void cache_hash_value() const override;

    //! Returns an identity bipartition.
    //!
    //! The *identity bipartition* of degree \f$n\f$ has blocks \f$\{i, -i\}\f$
    //! for all \f$i\in \{0, \ldots, n - 1\}\f$. This method returns a new
    //! identity bipartition of degree equal to the degree of \c this.
    Element* identity() const override;

    //! Multiply \p x and \p y and stores the result in \c this.
    //!
    //! This method redefines \c this to be the product (as defined at the top
    //! of this page) of the parameters  \p x and \p y. This method asserts
    //! that the degrees of \p x, \p y, and \c this, are all equal, and that
    //! neither \p x nor  \p y equals \c this.
    //!
    //! The parameter \p thread_id is required since some temporary storage is
    //! required to find the product of \p x and \p y.
    //! Note that if different threads call this method with the same value of
    //! \p thread_id then bad things will happen.
    void redefine(Element const* x,
                  Element const* y,
                  size_t const&  thread_id) override;

    //! Returns the number of transverse blocks.
    //!
    //! The *rank* of a bipartition is the number of blocks containing both
    //! positive and negative values.  This value is cached after it is first
    //! computed.
    size_t rank();

    //! Returns the index of the block containing \p pos.
    //!
    //! The parameter \p pos must be a value between *0* and *2n - 1* where *n*
    //! is the degree of this. This method asserts that pos is in the correct
    //! range of values.

    // FIXME remove this it is redundant: the method for [] or at of
    // ElementWithVectorData makes this unnecessary
    u_int32_t block(size_t pos) const;

    //! Returns the number of blocks in a bipartition.
    //!
    //! This method differs for Bipartition::nr_blocks in that the number of
    //! blocks is not cached if it has not been previously computed.
    u_int32_t const_nr_blocks() const;

    //! Returns the number of blocks in a bipartition.
    //!
    //! This value is cached the first time it is computed.
    u_int32_t nr_blocks();

    //! Returns the number of blocks containing a positive integer.
    //!
    //! The *left blocks* of a bipartition is the partition of
    //! \f$\{0, \ldots, n - 1\}\f$ induced by the bipartition. This method
    //! returns the number of blocks in this partition.
    u_int32_t nr_left_blocks();

    //! Returns the number of blocks containing a negative integer.
    //!
    //! The *right blocks* of a bipartition is the partition of
    //! \f$\{n, \ldots, 2n - 1\}\f$ induced by the bipartition. This method
    //! returns the number of blocks in this partition.
    u_int32_t nr_right_blocks();

    //! Returns \c true if the block with index \p index is transverse.
    //!
    //! A block of a biparition is *transverse* if it contains integers less
    //! than and greater than \f$n\f$, which is the degree of the bipartition.
    //! This method asserts that the parameter \p index is less than the number
    //! of blocks in the bipartition.
    bool is_transverse_block(size_t index);

    //! Return the left blocks of a bipartition
    //!
    //! The *left blocks* of a bipartition is the partition of
    //! \f$\{0, \ldots, n - 1\}\f$ induced by the bipartition. This method
    //! returns a Blocks object representing this partition.
    Blocks* left_blocks();

    //! Return the left blocks of a bipartition
    //!
    //! The *right blocks* of a bipartition is the partition of
    //! \f$\{n, \ldots, 2n - 1\}\f$ induced by the bipartition. This method
    //! returns a Blocks object representing this partition.
    Blocks* right_blocks();

    //! Set the cached number of blocks
    //!
    //! This method sets the cached value of the number of blocks of \c this
    //! to \p nr_blocks. It asserts that either there is no existing cached
    //! value or \p nr_blocks equals the existing cached value.
    inline void set_nr_blocks(size_t nr_blocks) {
      assert(_nr_blocks == Bipartition::UNDEFINED || _nr_blocks == nr_blocks);
      _nr_blocks = nr_blocks;
    }

    //! Set the cached number of left blocks
    //!
    //! This method sets the cached value of the number of left blocks of
    //! \c this to \p nr_left_blocks. It asserts that either there is no
    //! existing cached value or \p nr_left_blocks equals the existing cached
    //! value.
    inline void set_nr_left_blocks(size_t nr_left_blocks) {
      assert(_nr_left_blocks == Bipartition::UNDEFINED
             || _nr_left_blocks == nr_left_blocks);
      _nr_left_blocks = nr_left_blocks;
    }

    //! Set the cached rank
    //!
    //! This method sets the cached value of the rank of \c this to \p rank.
    //! It asserts that either there is no existing cached value or
    //! \p rank equals the existing cached value.
    inline void set_rank(size_t rank) {
      assert(_rank == Bipartition::UNDEFINED || _rank == rank);
      _rank = rank;
    }

   private:
    u_int32_t fuseit(std::vector<u_int32_t>& fuse, u_int32_t pos);
    void init_trans_blocks_lookup();

    static std::vector<std::vector<u_int32_t>> _fuse;
    static std::vector<std::vector<u_int32_t>> _lookup;

    size_t            _nr_blocks;
    size_t            _nr_left_blocks;
    std::vector<bool> _trans_blocks_lookup;
    size_t            _rank;

    static u_int32_t const UNDEFINED;
  };

  //! Class for square matrices over a Semiring.
  //!
  //! This class is abstract since it does not implement all methods required by
  //! the Element class, it exists to provide common methods for its
  //! subclasses.
  class MatrixOverSemiring
      : public ElementWithVectorData<int64_t, MatrixOverSemiring> {
   public:
    //! A constructor.
    //!
    //! Constructs a matrix defined by \p matrix, \p matrix is not copied,
    //! and should be deleted using ElementWithVectorData::really_delete.
    //!
    //! The parameter \p matrix should be a vector of integer values of length
    //! \f$n ^ 2\f$ for some integer \f$n\f$, so that the value in position
    //! \f$in + j\f$ is the entry in the \f$i\f$th row and \f$j\f$th column of
    //! the constructed matrix.
    //!
    //! The parameter \p hv must be the hash value of the element being
    //! created. It is the responsibility of the caller to ensure that \p hv is
    //! the correct value. See Element::Element for more details.
    //!
    //! The parameter \p semiring should be a pointer to a Semiring, which
    //! is the semiring over which the matrix is defined (this defaults to \c
    //! nullptr).
    explicit MatrixOverSemiring(std::vector<int64_t>* matrix,
                                size_t                hv,
                                Semiring*             semiring = nullptr)
        : ElementWithVectorData<int64_t, MatrixOverSemiring>(matrix, hv),
          _semiring(semiring) {}

    //! A constructor.
    //!
    //! Constructs a matrix defined by \p matrix, \p matrix is not copied,
    //! and should be deleted using ElementWithVectorData::really_delete.
    //!
    //! The parameter \p matrix should be a vector of integer values of length
    //! \f$n ^ 2\f$ for some integer \f$n\f$, so that the value in position
    //! \f$in + j\f$ is the entry in the \f$i\f$th row and \f$j\f$th column of
    //! the constructed matrix.
    //!
    //! The parameter \p semiring should be a pointer to a Semiring, which
    //! is the semiring over which the matrix is defined (this defaults to \c
    //! nullptr).
    explicit MatrixOverSemiring(std::vector<int64_t>* matrix,
                                Semiring*             semiring = nullptr)
        : ElementWithVectorData<int64_t, MatrixOverSemiring>(matrix),
          _semiring(semiring) {}

    //! A constructor.
    //!
    //! Constructs a matrix defined by \p matrix, which is copied into the
    //! constructed object.
    //!
    //! The parameter \p matrix should be a vector of integer values of length
    //! \f$n ^ 2\f$ for some integer \f$n\f$, so that the value in position
    //! \f$in + j\f$ is the entry in the \f$i\f$th row and \f$j\f$th column of
    //! the constructed matrix.
    //!
    //! The parameter \p semiring should be a pointer to a Semiring, which
    //! is the semiring over which the matrix is defined.
    explicit MatrixOverSemiring(std::vector<std::vector<int64_t>> const& matrix,
                                Semiring* semiring);

    //! Returns a pointer to the Semiring over which the matrix is defined.
    Semiring* semiring() const;

    //! Returns the approximate time complexity of multiplying two matrices.
    //!
    //! The approximate time complexity of multiplying matrices is \f$n ^ 3\f$
    //! where \f$n\f$ is the dimension of the matrix.
    size_t complexity() const override;

    //! Returns the dimension of the matrix.
    //!
    //! The *dimension* of a matrix is just the number of rows (or,
    //! equivalently columns).
    //!
    size_t degree() const override;

    //! Find the hash value of a matrix over a semiring.
    //!
    //! \sa Element::hash_value.
    void cache_hash_value() const override;

    //! Returns the identity matrix with dimension of \c this.
    //!
    //! This method returns a new matrix with dimension equal to that of \c
    //! this, where the main diagonal consists of the value Semiring::one and
    //! every other entry Semiring::zero.
    Element* identity() const override;

    //! Returns a pointer to a copy of \c this.
    //!
    //! The parameter \p increase_deg_by must be 0, since it does not make
    //! sense to increase the degree of a matrix. The semiring of the \c this
    //! and the returned copy are equal, so if the method
    //! Element::really_delete is called on either, then the other will become
    //! invalid.
    Element* really_copy(size_t increase_deg_by = 0) const override;

    //! Multiply \p x and \p y and stores the result in \c this.
    //!
    //! This method asserts that the degrees of \p x, \p y, and \c this, are
    //! all equal, and that neither \p x nor \p y equals \c this. It does not
    //! currently verify that \c x, \c y, and \c this are defined over the same
    //! semiring.

    // FIXME verify that x, y, and this are defined over the same semiring.
    void redefine(Element const* x, Element const* y) override;

   private:
    //! a function applied after redefinition
    virtual void after() {}

    Semiring* _semiring;
  };

  //! Class for projective max-plus matrices.
  //!
  //! These matrices belong to the quotient of the monoid of all max-plus
  //! matrices by the congruence where two matrices are related if they differ
  //! by a scalar multiple. \sa MaxPlusSemiring and MatrixOverSemiring.
  //!
  //! Matrices in this class are modified when constructed to be in a normal
  //! form which is obtained by subtracting the maximum finite entry in
  //! the matrix from the every finite entry.
  class ProjectiveMaxPlusMatrix : public MatrixOverSemiring {
   public:
    //! A constructor.
    //!
    //! See MatrixOverSemiring::MatrixOverSemiring for details about this
    //! method.
    //!
    //! The parameter \p matrix is converted into its normal form when
    //! when the object is constructed.
    explicit ProjectiveMaxPlusMatrix(std::vector<int64_t>* matrix,
                                     size_t    hv       = Element::UNDEFINED,
                                     Semiring* semiring = nullptr)
        : MatrixOverSemiring(matrix, hv, semiring) {
      after();  // this is to put the matrix in normal form
    }

    //! A constructor.
    //!
    //! See MatrixOverSemiring::MatrixOverSemiring for details about this
    //! method.
    //!
    //! The parameter \p matrix is converted into its normal form when
    //! when the object is constructed.
    explicit ProjectiveMaxPlusMatrix(std::vector<int64_t>* matrix,
                                     Semiring*             semiring = nullptr)
        : MatrixOverSemiring(matrix, semiring) {
      after();  // this is to put the matrix in normal form
    }

    //! A constructor.
    //!
    //! See MatrixOverSemiring::MatrixOverSemiring for details about this
    //! method.
    //!
    //! The copy of the parameter \p matrix is converted into its normal form
    //! when the object is constructed.
    explicit ProjectiveMaxPlusMatrix(
        std::vector<std::vector<int64_t>> const& matrix,
        Semiring*                                semiring = nullptr)
        : MatrixOverSemiring(matrix, semiring) {
      after();  // this is to put the matrix in normal form
    }

    //! Find the hash value of a projective max-plus matrix.
    //!
    //! \sa Element::hash_value.
    void cache_hash_value() const override;

   private:
    //! a function applied after redefinition
    void after() override;
  };

  //! Class for partitioned binary relations (PBR).
  //!
  //! Partitioned binary relations (PBRs) are a generalisation of bipartitions,
  //! which were introduced by
  //! [Martin and Mazorchuk](https://arxiv.org/abs/1102.0862).
  class PBR : public ElementWithVectorData<std::vector<u_int32_t>, PBR> {
   public:
    //! Constructs a PBR defined by the vector pointed to by \p vector, which
    //! is not copied, and should be deleted using
    //! ElementWithVectorData::really_delete.
    //!
    //! The parameter \p vector should be a pointer to a vector of vectors of
    //! non-negative integer values of length \f$2n\f$ for some integer
    //! \f$n\f$, the vector in position \f$i\f$ is the list of points adjacent
    //! to \f$i\f$ in the PBR.
    explicit PBR(std::vector<std::vector<u_int32_t>>* vector,
                 size_t                               hv = Element::UNDEFINED)
        : ElementWithVectorData<std::vector<u_int32_t>, PBR>(vector, hv) {}

    //! Constructs a PBR defined by \p vector, where \p vector is copied into
    //! the PBR being constructed.
    //!
    //! The parameter \p vector should be a vector of vectors of non-negative
    //! integer values of length \f$2n\f$ for some integer \f$n\f$, the vector
    //! in position \f$i\f$ is the list of points adjacent to \f$i\f$ in the
    //! PBR.
    explicit PBR(std::vector<std::vector<u_int32_t>> const& vector)
        : PBR(new std::vector<std::vector<u_int32_t>>(vector)) {}

    //! Returns the approximate time complexity of multiplying PBRs.
    //!
    //! The approximate time complexity of multiplying PBRs is \f$2n ^ 3\f$
    //! where \f$n\f$ is the degree.
    size_t complexity() const override;

    //! Returns the degree of a PBR.
    //!
    //! The *degree* of a PBR is half the number of points in the PBR, which is
    //! also half the length of the underlying vector
    //! ElementWithVectorData::_vector.
    size_t degree() const override;

    //! Find the hash value of a PBR.
    //!
    //! \sa Element::hash_value.
    void cache_hash_value() const override;

    //! Returns the identity PBR with degree equal to that of \c this.
    //!
    //! This method returns a new PBR with degree equal to the degree of \c
    //! this where every value is adjacent to its negative. Equivalently,
    //! \f$i\f$ is adjacent \f$i + n\f$ and vice versa for every \f$i\f$ less
    //! than the degree \f$n\f$.
    Element* identity() const override;

    //! Multiply \p x and \p y and stores the result in \c this.
    //!
    //! This method redefines \c this to be the product
    //! of the parameters  \p x and \p y. This method asserts
    //! that the degrees of \p x, \p y, and \c this, are all equal, and that
    //! neither \p x nor  \p y equals \c this.
    //!
    //! The parameter \p thread_id is required since some temporary storage is
    //! required to find the product of \p x and \p y.  Note that if different
    //! threads call this method with the same value of \p thread_id then bad
    //! things will happen.
    void redefine(Element const* x,
                  Element const* y,
                  size_t const&  thread_id) override;

   private:
    void unite_rows(RecVec<bool>& out,
                    RecVec<bool>& tmp,
                    size_t const& vertex1,
                    size_t const& vertex2);

    void x_dfs(std::vector<bool>& x_seen,
               std::vector<bool>& y_seen,
               RecVec<bool>&      tmp,
               u_int32_t const&   n,
               u_int32_t const&   i,
               PBR const* const   x,
               PBR const* const   y,
               size_t const&      adj);

    void y_dfs(std::vector<bool>& x_seen,
               std::vector<bool>& y_seen,
               RecVec<bool>&      tmp,
               u_int32_t const&   n,
               u_int32_t const&   i,
               PBR const* const   x,
               PBR const* const   y,
               size_t const&      adj);

    static std::vector<std::vector<bool>> _x_seen;
    static std::vector<std::vector<bool>> _y_seen;
    static std::vector<RecVec<bool>>      _out;
    static std::vector<RecVec<bool>>      _tmp;
  };

  // Schreier-Sims set up

  u_int16_t nr_new_perm_coll;

  struct perm_coll {
    Permutation<u_int16_t>* gens;
    u_int16_t nr_gens;
    u_int16_t deg;
    u_int16_t alloc_size;
  };

  typedef struct perm_coll PermColl;

  PermColl* new_perm_coll(u_int16_t upper_bound) {
    nr_new_perm_coll++;
    PermColl* coll = malloc(sizeof(PermColl));
    nr_ss_allocs++;
    coll->gens = malloc(upper_bound * sizeof(Perm));
    nr_ss_allocs++;
    coll->nr_gens    = 0;
    coll->deg        = deg;
    coll->alloc_size = upper_bound;
    return coll;
  }

  static inline void add_strong_gens(u_int16_t const pos, Permutation<u_int16_t> const value) {
    if (strong_gens[pos] == NULL) {
      strong_gens[pos] = new_perm_coll(1);
    }
    add_perm_coll(strong_gens[pos], value);
  }

  static inline Permutation<u_int16_t> get_strong_gens(u_int16_t const i, u_int16_t const j) {
    return strong_gens[i]->gens[j];
  }

  static inline Permutation<u_int16_t> get_transversal(u_int16_t const i, u_int16_t const j) {
    return transversal[i * MAXVERTS + j];
  }

  static inline Permutation<u_int16_t> get_transversal_inv(u_int16_t const i, u_int16_t const j) {
    return transversal_inv[i * MAXVERTS + j];
  }

  static inline void
  set_transversal(u_int16_t const i, u_int16_t const j, Permutation<u_int16_t> const value) {
    // free the perm in this position if there is one already
    if (transversal[i * MAXVERTS + j] != NULL) {
      free(transversal[i * MAXVERTS + j]);
      nr_ss_frees++;
      free(transversal_inv[i * MAXVERTS + j]);
      nr_ss_frees++;
    }
    transversal[i * MAXVERTS + j]     = value;
    transversal_inv[i * MAXVERTS + j] = invert_perm(value);
  }

  static bool perm_fixes_all_base_points(Perm const x) {
    u_int16_t i;

    for (i = 0; i < size_base; i++) {
      if (x[base[i]] != base[i]) {
        return false;
      }
    }
    return true;
  }

  static inline void add_base_point(u_int16_t const pt) {
    base[size_base]               = pt;
    size_orbits[size_base]        = 1;
    orbits[size_base * deg]       = pt;
    borbits[size_base * deg + pt] = true;
    set_transversal(size_base, pt, id_perm());
    size_base++;
  }

  static inline void first_ever_init() {
    first_ever_call = false;
    memset((void*) size_orbits, 0, MAXVERTS * sizeof(u_int16_t));
  }

  static void init_stab_chain() {
    if (first_ever_call) {
      first_ever_init();
    }

    memset((void*) borbits, false, deg * deg * sizeof(bool));
    size_base = 0;
  }

  /*static void init_endos_base_points() {
    u_int16_t  i;

    for (i = 0; i < deg - 1; i++) {
      add_base_point(i);
    }
  }*/

  static void free_stab_chain() {
    int i, j, k;

    memset((void*) size_orbits, 0, size_base * sizeof(u_int16_t));

    // free the transversal
    // free the transversal_inv
    for (i = 0; i < (int) deg; i++) {
      for (j = 0; j < (int) deg; j++) {
        k = i * MAXVERTS + j;
        if (transversal[k] != NULL) {
          free(transversal[k]);
          transversal[k] = NULL;
          nr_ss_frees++;
          free(transversal_inv[k]);
          transversal_inv[k] = NULL;
          nr_ss_frees++;
        }
      }
    }

    // free the strong_gens
    for (i = 0; i < (int) size_base; i++) {
      if (strong_gens[i] != NULL) {
        free_perm_coll(strong_gens[i]);
        strong_gens[i] = NULL;
      }
    }
  }

  static void orbit_stab_chain(u_int16_t const depth, u_int16_t const init_pt) {
    u_int16_t i, j, pt, img;
    Perm  x;

    assert(depth <= size_base);  // Should this be strict?

    for (i = 0; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->nr_gens; j++) {
        x   = get_strong_gens(depth, j);
        img = x[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          set_transversal(depth, img, prod_perms(get_transversal(depth, pt), x));
        }
      }
    }
  }

  static void add_gen_orbit_stab_chain(u_int16_t const depth, Permutation<u_int16_t> const gen) {
    u_int16_t i, j, pt, img;
    Permutation<u_int16_t>  x;

    assert(depth <= size_base);

    // apply the new generator to existing points in orbits[depth]
    u_int16_t nr = size_orbits[depth];
    for (i = 0; i < nr; i++) {
      pt  = orbits[depth * deg + i];
      img = gen[pt];
      if (!borbits[depth * deg + img]) {
        orbits[depth * deg + size_orbits[depth]] = img;
        size_orbits[depth]++;
        borbits[depth * deg + img] = true;
        set_transversal(depth, img, prod_perms(get_transversal(depth, pt), gen));
      }
    }

    for (i = nr; i < size_orbits[depth]; i++) {
      pt = orbits[depth * deg + i];
      for (j = 0; j < strong_gens[depth]->nr_gens; j++) {
        x   = get_strong_gens(depth, j);
        img = x[pt];
        if (!borbits[depth * deg + img]) {
          orbits[depth * deg + size_orbits[depth]] = img;
          size_orbits[depth]++;
          borbits[depth * deg + img] = true;
          set_transversal(depth, img, prod_perms(get_transversal(depth, pt), x));
        }
      }
    }
  }

  static void sift_stab_chain(Permutation<u_int16_t>* g, u_int16_t* depth) {
    u_int16_t beta;

    assert(*depth == 0);

    for (; *depth < size_base; (*depth)++) {
      beta = (*g)[base[*depth]];
      if (!borbits[*depth * deg + beta]) {
        return;
      }
      prod_perms_in_place(*g, get_transversal_inv(*depth, beta));
    }
  }

  static void schreier_sims_stab_chain(u_int16_t const depth) {
    Permutation<u_int16_t>  x, h, prod;
    bool  escape, y;
    int   i;
    u_int16_t j, jj, k, l, m, beta, betax;

    for (i = 0; i <= (int) depth; i++) {
      for (j = 0; j < strong_gens[i]->nr_gens; j++) {
        x = get_strong_gens(i, j);
        if (perm_fixes_all_base_points(x)) {
          for (k = 0; k < deg; k++) {
            if (k != x[k]) {
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
        if (beta == x[beta]) {
          add_strong_gens(i, copy_perm(x));
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
          x     = get_strong_gens(i, m);
          prod  = prod_perms(get_transversal(i, beta), x);
          betax = x[beta];
          if (!eq_perms(prod, get_transversal(i, betax))) {
            y  = true;
            h  = prod_perms(prod, get_transversal_inv(i, betax));
            jj = 0;
            sift_stab_chain(&h, &jj);
            if (jj < size_base) {
              y = false;
            } else if (!is_one(h)) {  // better method? IsOne(h)?
              y = false;
              for (k = 0; k < deg; k++) {
                if (k != h[k]) {
                  add_base_point(k);
                  break;
                }
              }
            }

            if (!y) {
              for (l = i + 1; l <= jj; l++) {
                add_strong_gens(l, copy_perm(h));
                add_gen_orbit_stab_chain(l, h);
                // add generator to <h> to orbit of base[l]
              }
              i      = jj;
              escape = true;
            }
            free(h);
            nr_ss_frees++;
          }
          free(prod);
          nr_ss_frees++;
        }
      }
      if (!escape) {
        i--;
      }
    }
  }

  extern bool point_stabilizer(PermColl* gens, u_int16_t const pt, PermColl** out) {
    init_stab_chain();

    strong_gens[0] = copy_perm_coll(gens);
    add_base_point(pt);
    schreier_sims_stab_chain(0);

    // The stabiliser we want is the PermColl pointed to by <strong_gens[1]>
    // UNLESS <strong_gens[1]> doesn't exists - this means that <strong_gens[0]>
    // is the stabilizer itself (????)
    if (*out != NULL) {
      free_perm_coll(*out);
    }
    if (strong_gens[1] == NULL) {
      // this means that the stabilizer of pt under <gens> is trivial
      *out = new_perm_coll(1);
      add_perm_coll(*out, id_perm());
      free_stab_chain();
      return true;
    }
    *out = copy_perm_coll(strong_gens[1]);
    free_stab_chain();
    return false;
  }

}  // namespace libsemigroups
#endif  // LIBSEMIGROUPS_SRC_ELEMENTS_H_
