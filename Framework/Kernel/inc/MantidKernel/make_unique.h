#ifndef Mantid_make_unique_h
#define Mantid_make_unique_h

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

// implementation of memory::make_unique for platforms that don't currently have
// it.
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm

namespace Mantid {

namespace Kernel {

#if __cplusplus > 201103L // C++14

using std::make_unique;

#else  // C++11

template <class T> struct _Unique_if {
  typedef std::unique_ptr<T> _Single_object;
};

template <class T> struct _Unique_if<T[]> {
  typedef std::unique_ptr<T[]> _Unknown_bound;
};

template <class T, size_t N> struct _Unique_if<T[N]> {
  struct __invalid_type {};
};

template <class T>
typename _Unique_if<T>::_Unknown_bound make_unique(size_t n) {
  typedef typename std::remove_extent<T>::type U;
  return std::unique_ptr<T>(new U[n]());
}

template <class T, class... Args>
inline typename _Unique_if<T>::_Single_object make_unique(Args &&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <class T, class... Args>
inline typename _Unique_if<T>::_Known_bound make_unique(Args &&...) = delete;
#endif // __cplusplus == 201402L
} // namespace Kernel
} // namespace Mantid
#endif // Mantid_make_unique_h
