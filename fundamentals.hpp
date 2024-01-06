#pragma once
#include<vector>
#include<iostream>
#include<iomanip>
#include<initializer_list>
#include<cmath>
#include<numeric>
#include<functional>
#include<type_traits>

// Here we define out fundamental structures:
// 1. TensorNd, Tensor3d, Tensor2d: in order to be able to include this header
//    file to multiple locations we NEED to split the implementation of each
//    method into a different folder ('fundamentals.cpp')
// 2. Vector template class: since this is a template (rule) we can keep the
//    implementation in the same folder!!
//
namespace ueg {
namespace fundamentals {
    typedef std::vector<double> dvec_t;

    struct TensorNd {
      std::vector<size_t> dimensions;
      std::vector<double> values;
      size_t n_axes;

      // Constructors
      TensorNd() {};
      TensorNd(const std::initializer_list<size_t> &dimensions_);
      TensorNd(const TensorNd &tensor);

      // Routine to set all elements to zero
      void clear() { for (auto &element : this->values) element = 0.0; };

      // Routine to get the total size of the tensor
      const size_t size() const;
      
      // Routine to rescale all values of the tensor
      void rescale(double prefactor);

      // Inplace setting from operation ( no copies involved )
      void multiply_by(const TensorNd &other);
      void subtract_by(const TensorNd &other);

      void set_from_product(const double scalar, const TensorNd &tensor_nd);
      void set_from_product(double scalar, const dvec_t &vec);
      void set_from_product(const dvec_t &vec, double scalar);
      void set_from_product(const dvec_t &vec1, const dvec_t &vec2);
      void set_from_division(const dvec_t &vec1, const dvec_t &vec2);
      void set_from_addition(const dvec_t &vec1, const dvec_t &vec2);
      void set_from_subtraction(const dvec_t &vec1, const dvec_t &vec2);
   
      // Generic way to access an element of the tensor
      const size_t globalIndex(const std::vector<size_t> &indices) const;

      // Operator overloads ( getters/setters )
      double  operator[] (const std::vector<size_t> &indices) const { return values[globalIndex(indices)]; };
      double &operator[] (const std::vector<size_t> &indices) { return values[globalIndex(indices)]; };

      // Routine to find the maximum and minimum element of the tensor
      std::tuple<double, double> get_min_max_value() const;

      // Write to and read from binary stream in native binary encoding
      void write_tensor_elements(std::ostream &stream);
      void read_tensor_elements(std::istream &stream);
    };

    // --------------------------------------------------------------------------------------------

    struct Tensor2d : public TensorNd {
        // Member variables as specifications of 2D tensor
      size_t dim1, dim2;
      Tensor2d() : TensorNd() {};
      Tensor2d(const std::initializer_list<size_t> &dimensions_);
      Tensor2d(size_t i, size_t j);
      Tensor2d(const Tensor2d &tensor);

      double  operator() (const size_t &i, const size_t &j) const { return this->values[dim1*j + i]; };
      double &operator() (const size_t &i, const size_t &j) { return this->values[dim1*j + i]; };
    };
    
    // --------------------------------------------------------------------------------------------
    
    struct Tensor3d : public TensorNd {
      // Member variables as specifications of 3D tensor
      size_t dim1, dim2, dim3;
      
      Tensor3d() : TensorNd() {};
      Tensor3d(const std::initializer_list<size_t> &dimensions_);
      Tensor3d(size_t i, size_t j, size_t k);
      Tensor3d(const Tensor3d &tensor);

      double  operator() (const size_t &i, const size_t &j, const size_t &k) const { return this->values[dim1*(dim2*k + j) + i]; };
      double &operator() (const size_t &i, const size_t &j, const size_t &k) { return this->values[dim1*(dim2*k + j) + i]; };

      void print_elements() const;
    };

    struct Tensor5d : public TensorNd {
      // Member variables as specifications of 3D tensor
      size_t dim1, dim2, dim3, dim4, dim5;
      
      Tensor5d() : TensorNd() {};
      Tensor5d(const std::initializer_list<size_t> &dimensions_);
      Tensor5d(size_t i, size_t j, size_t k, size_t a, size_t b);
      Tensor5d(const Tensor5d &tensor);

      double  operator() (const size_t &i, const size_t &j, const size_t &k, const size_t &a, const size_t &b) const { return this->values[dim1*(dim2*(dim3*(dim4*a+b) +k ) + j) + i]; };
      double &operator() (const size_t &i, const size_t &j, const size_t &k, const size_t &a, const size_t &b) { return this->values[dim1*(dim2*(dim3*(dim4*a+b) +k) + j) + i]; };

      void print_elements() const;
    };

    // --------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------

    // Definition of 3 components Vector class, (x,y,z)
    template <typename DataType>
    struct Vector {
      // Member variables
      DataType x,y,z;
      
      // Constructors
      Vector() {};
      Vector(DataType x_, DataType y_, DataType z_) : x(x_), y(y_), z(z_) {};
      
      // Initialize from another vector and allow for type-casting
      template<typename T>
      Vector(const Vector<T> &v) {
        if (std::is_same<DataType,T>::value) { x = v.x ; y = v.y ; z = v.z ;
        } else {
          x = static_cast<int>(v.x >= 0 ? v.x+0.5 : v.x-0.5); // Ternary operators
          y = static_cast<int>(v.y >= 0 ? v.y+0.5 : v.y-0.5);
          z = static_cast<int>(v.z >= 0 ? v.z+0.5 : v.z-0.5);
        };
      }

      // Plain dot product between two vectors ( of potentially different type )
      template<typename T>
      auto dot(const Vector<T> &rhs) const { return x*rhs.x + y*rhs.y + z*rhs.z; };

      // Plain cross product between two vectors ( of potentially different type )
      template<typename T>
      auto cross(const Vector<T> &rhs) const { return Vector<decltype(x*rhs.x)>(y*rhs.z-z*rhs.y, z*rhs.x-x*rhs.z, x*rhs.y-y*rhs.x); };

      // Squared length of the vector
      const DataType squaredLength() const { return x*x + y*y + z*z; }

      // Norm of the vector
      const double norm() const { return std::sqrt(this->squaredLength()); }
      
      // Distance between two vectors | v1 - v2 | ** 2
      const double distance(const Vector<DataType> &vec) const { return (*this - vec).norm(); };
      
      // Rounding function (isn't obsolete to have it as a template??)
      template <typename Type>
      Vector<Type> round() { return Vector<Type>(std::round(x), std::round(y), std::round(z)); }
      
      // Overloads of unary operators
      Vector<DataType> operator - () const { return Vector<DataType>(-x,-y,-z); };
      Vector<DataType> operator + () const { return Vector<DataType>(+x,+y,+z); };
    };
    
    // Comparison operators between vectors
    template<typename Lhs, typename Rhs>
    inline bool operator == (const Vector<Lhs> &lhs, const Vector<Rhs> &rhs) {
      return ( (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z) );
    };

    template<typename Lhs, typename Rhs>
    inline bool operator != (const Vector<Lhs> &lhs, const Vector<Rhs> &rhs) {
      return ( (lhs.x != rhs.x) || (lhs.y != rhs.y) || (lhs.z != rhs.z) );
    };

    template<typename DataType>
    inline bool operator < (const Vector<DataType> &lhs, const Vector<DataType> &rhs) {
      double length_difference = lhs.squaredLength() - rhs.squaredLength();
      if ( std::abs(length_difference) > 1e-15 ) return length_difference < 0;

      // At this point we are sure that we have to do with 'degenerate' states
      // and we sort them 'lexicographically'
      if ( lhs.x != rhs.x ) {
        return lhs.x < rhs.x;
      }
      else {
        if ( lhs.y != rhs.y )
          return lhs.y < rhs.y;
        else
          return lhs.z < rhs.z;
      }
    };

    // Addition / Subtraction between two vectors
    template<typename Lhs, typename Rhs>
    inline auto operator + (const Vector<Lhs> &lhs, const Vector<Rhs> &rhs) {
      return Vector<decltype(lhs.x+rhs.x)>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
    };

    template<typename Lhs, typename Rhs>
    inline auto operator - (const Vector<Lhs> &lhs, const Vector<Rhs> &rhs) {
      return Vector<decltype(lhs.x-rhs.x)>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
    };

    // Multiplication / Division with scalars
    template<typename DataType, typename ScalarType>
    inline auto operator * (const Vector<DataType> &vec, const ScalarType &scalar) {
      static_assert(std::is_scalar<ScalarType>::value, "Multiplication is defined only for scalars");
      return Vector<decltype(vec.x * scalar)>(vec.x*scalar, vec.y*scalar, vec.z*scalar);
    };

    template<typename ScalarType, typename DataType>
    inline auto operator * (const ScalarType &scalar, const Vector<DataType> &vec) {
      static_assert(std::is_scalar<ScalarType>::value, "Multiplication is defined only for scalars");
      return Vector<decltype(vec.x * scalar)>(vec.x*scalar, vec.y*scalar, vec.z*scalar);
    };
    
    template<typename DataType, typename ScalarType>
    inline auto operator / (const Vector<DataType> &vec, const ScalarType &scalar) {
      if ( scalar == 0 ) throw std::invalid_argument("(!) Division by 0 encountered");
      if (std::is_same<DataType,ScalarType>::value && std::is_same<DataType,int>::value)
        throw std::invalid_argument("Division of integer vector with integer scalar is ill-defined!");
      else
        return Vector<decltype(vec.x / scalar)>(vec.x/scalar, vec.y/scalar, vec.z/scalar);
    };

    // Overload for the std::cout operator (printing out)
    template <typename DataType>
    std::ostream &operator << (std::ostream &out, const Vector<DataType> &vec) {
      out << std::showpos  << "(" << vec.x << "," << vec.y << "," << vec.z << ")" << std::noshowpos;
      return out;
    };

  }; // end namespace fundamentals
}; // end namepsace ueg

