#include "fundamentals.hpp"

namespace ueg {
  namespace fundamentals {
  
    // Implementation of TensorNd constructors / methods
    TensorNd::TensorNd(const std::initializer_list<size_t> &dimensions_) :
      dimensions(dimensions_), n_axes(dimensions.size())
    { 
      size_t size = std::accumulate(dimensions.begin(), dimensions.end(), 1UL, std::multiplies<size_t>());
      values.resize(size);
    };

    TensorNd::TensorNd(const TensorNd &tensor) :
        dimensions(tensor.dimensions), values(tensor.values), n_axes(tensor.n_axes) {};

    const size_t TensorNd::globalIndex(const std::vector<size_t> &indices) const {
      size_t index = 0;
      for (int dim(n_axes-1); dim >= 0; --dim) {
        index *= dimensions[dim];
        index += indices[dim];
      };
      return index;
    };

    const size_t TensorNd::size() const {
      return std::accumulate(dimensions.begin(), dimensions.end(), 1UL, std::multiplies<size_t>());
    };

    void TensorNd::rescale(double prefactor) {
      for (auto &element : values) element *= prefactor;
    };

    void TensorNd::multiply_by(const TensorNd &other) {
      if (this->size() != other.size())
        throw std::invalid_argument("Size of tensor does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] *= other.values[index];
    }

    void TensorNd::subtract_by(const TensorNd &other) {
      if (this->size() != other.size())
        throw std::invalid_argument("Size of tensor does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] -= other.values[index];
    }

    void TensorNd::set_from_product(const double scalar, const TensorNd &tensor_nd) {
      set_from_product(scalar, tensor_nd.values);
    };

    void TensorNd::set_from_product(const dvec_t &vec, double scalar) {
      if (this->size() != vec.size())
        throw std::invalid_argument("Size of tensor does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] = scalar * vec[index];
    };

    void TensorNd::set_from_product(double scalar, const dvec_t &vec) {
      this->set_from_product(vec,scalar);
    };

    void TensorNd::set_from_product(const dvec_t &vec1, const dvec_t &vec2) {
      if ( (this->size() != vec1.size()) || (this->size() != vec2.size()) )
        throw std::invalid_argument("Size of tensors does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] = vec1[index] * vec2[index];
    };

    void TensorNd::set_from_division(const dvec_t &vec1, const dvec_t &vec2) {
      if ( (this->size() != vec1.size()) || (this->size() != vec2.size()) )
        throw std::invalid_argument("Size of tensors does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] = vec1[index] / vec2[index];
    };

    void TensorNd::set_from_addition(const dvec_t &vec1, const dvec_t &vec2) {
      if ( (this->size() != vec1.size()) || (this->size() != vec2.size()) )
        throw std::invalid_argument("Size of tensors does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] = vec1[index] + vec2[index];
    };

    void TensorNd::set_from_subtraction(const dvec_t &vec1, const dvec_t &vec2) {
      if ( (this->size() != vec1.size()) || (this->size() != vec2.size()) )
        throw std::invalid_argument("Size of tensors does not match");
      for (size_t index(0); index < this->size(); ++index)
        this->values[index] = vec1[index] - vec2[index];
    };
    
    std::tuple<double, double> TensorNd::get_min_max_value() const{
      double minimum, maximum;
      auto vals = this->values;
      vals.erase(std::remove(vals.begin(), vals.end(), 0),vals.end());
      vals.shrink_to_fit();
      if (vals.size() != 0) {
        maximum = *std::max_element(vals.begin(), vals.end());
        minimum = *std::min_element(vals.begin(), vals.end());
      } else {
        maximum = 0.0;
        minimum = 0.0;
      };
      
      return std::make_tuple(minimum, maximum);
    };

    void TensorNd::write_tensor_elements(std::ostream &stream) {
      stream.write(reinterpret_cast<const char *>(values.data()), sizeof(double)*values.size());
    }

    void TensorNd::read_tensor_elements(std::istream &stream) {
      stream.read(reinterpret_cast<char *>(values.data()), sizeof(double)*values.size());
    }

    // Implementation of Tensor3d and Tensor2d methods / constructors
    Tensor2d::Tensor2d(const std::initializer_list<size_t> &dimensions_) : TensorNd(dimensions_) {
      dim1 = this->dimensions[0]; dim2 = this->dimensions[1];
    };

    Tensor2d::Tensor2d(size_t i, size_t j) : TensorNd({i,j}) {
      dim1 = i; dim2 = j;
    };

    Tensor2d::Tensor2d(const Tensor2d &tensor) : TensorNd({tensor.dim1,tensor.dim2}) {
      dim1 = tensor.dim1; dim2 = tensor.dim2;
      this->values = tensor.values;
    };

    Tensor3d::Tensor3d(const std::initializer_list<size_t> &dimensions_) : TensorNd(dimensions_) {
      dim1 = this->dimensions[0]; dim2 = this->dimensions[1]; dim3 = this->dimensions[2];
    };
    
    Tensor3d::Tensor3d(size_t i, size_t j, size_t k) : TensorNd({i,j,k}) {
      dim1 = i; dim2 = j; dim3 = k;
    };
      
    Tensor3d::Tensor3d(const Tensor3d &tensor) : TensorNd({tensor.dim1,tensor.dim2,tensor.dim3}) {
      dim1 = tensor.dim1; dim2 = tensor.dim2; dim3 = tensor.dim3;
      this->values = tensor.values;
    };

    void Tensor3d::print_elements() const {
      for (size_t i(0); i < dim1; ++i)
        for (size_t j(0); j < dim2; ++j)
          for (size_t k(0); k < dim3; ++k)
              std::cout << std::setw(10) << i << std::setw(10) << j << std::setw(10) << k
                        << std::setw(20) << this->values[globalIndex({i,j,k})] << std::endl;
    };

    Tensor5d::Tensor5d(const std::initializer_list<size_t> &dimensions_) : TensorNd(dimensions_) {
      dim1 = this->dimensions[0]; dim2 = this->dimensions[1]; dim3 = this->dimensions[2]; dim4 = this->dimensions[3]; dim5 = this->dimensions[4];
    };

    Tensor5d::Tensor5d(size_t i, size_t j, size_t k, size_t a, size_t b) : TensorNd({i,j,k,a,b}) {
      dim1 = i; dim2 = j; dim3 = k; dim4 = a; dim5 = b;
    };

    Tensor5d::Tensor5d(const Tensor5d &tensor) : TensorNd({tensor.dim1,tensor.dim2,tensor.dim3,tensor.dim4,tensor.dim5}) {
      dim1 = tensor.dim1; dim2 = tensor.dim2; dim3 = tensor.dim3; dim4 = tensor.dim4; dim5 = tensor.dim5;
      this->values = tensor.values;
    };
  
  }; // end namespace fundamentals
}; // end namespace ueg
