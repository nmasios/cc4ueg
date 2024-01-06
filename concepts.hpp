#pragma once
#include "fundamentals.hpp"

// The 'concept' structures defined above inherit all functionality and
// structure from a LocalTensor and (optimally) adds new functionality.
// In the case that no additional functions are written you should consider
// the corresponding structure as a pure LocalTensor!! It's in general a
// good practice to have multiple classes that define a specific mathematical
// or physical concept so that if you need to add something in the future
// you do not mess up the rest.

namespace ueg {
  namespace concepts {
    typedef ueg::fundamentals::Tensor3d tensor3d_t;
     
    struct Residuum : public tensor3d_t {
      Residuum() : tensor3d_t() {};
      Residuum(const std::initializer_list<size_t> &dimensions): tensor3d_t(dimensions) {};
      void info() const { std::cout << "Hello from Residuum" << std::endl; };

      // A good example of 'extended' functionality is the function
      void symmetrize() { std::cout << "Write the symmetrization of the Residuum" << std::endl; }; 

      void info() { 
        double minimum,maximum;
        std::tie(minimum,maximum) = get_min_max_value();

        std::cout << "- Information of the Residuum" << std::endl;
        std::cout << "\t- Total elements: " << this->size() << std::endl;
        std::cout << "\t- Minimum value (non-zero): " << minimum << std::endl; 
        std::cout << "\t- Maximum value (non-zero): " << maximum << std::endl;
        std::cout << "\t- Dimensions              : " << dimensions[0] << "," 
                  << dimensions[1] << "," << dimensions[2] << std::endl;
      };

      void print() const { 
        std::cout << "- Residuum elements. Using the conversion (i,j,a,value)" << std::endl;
        this->print_elements(); 
      };
    }; // end struct Residuum 


    struct Propagator : public tensor3d_t {
      Propagator() : tensor3d_t() {};
      Propagator(const std::initializer_list<size_t> &dimensions): tensor3d_t(dimensions) {};
      void info() { 
        double maximum, minimum; 
        auto vals = this->values; 
        vals.erase(std::remove(vals.begin(), vals.end(), 0),vals.end());
        vals.shrink_to_fit();
        
        if ( vals.size() != 0 ) { 
          maximum = *std::max_element(vals.begin(), vals.end());
          minimum = *std::min_element(vals.begin(), vals.end());
        } else { 
          maximum = 0.0;
          minimum = 0.0; 
        };

        std::cout << "- Information of the Propagator" << std::endl;
        std::cout << "\t- Total elements: " << this->size() << std::endl;
        std::cout << "\t- Minimum value (non-zero): " << minimum << std::endl; 
        std::cout << "\t- Maximum value (non-zero): " << maximum << std::endl;
        std::cout << "\t- Dimensions              : " << dimensions[0] << "," 
                  << dimensions[1] << "," << dimensions[2] << std::endl;
      };
    
      void print() const { 
        std::cout << "- Propagator elements. Using the conversion (i,j,a,value)" << std::endl;
        this->print_elements(); 
      };
    }; // end struct Propagator


    struct CoulombInteraction : public tensor3d_t {
      CoulombInteraction() : tensor3d_t() {};
      CoulombInteraction(const std::initializer_list<size_t> &dimensions): tensor3d_t(dimensions) {};
      void info() { 
        double maximum, minimum; 
        auto vals = this->values; 
        vals.erase(std::remove(vals.begin(), vals.end(), 0),vals.end());
        vals.shrink_to_fit();

        if ( vals.size() != 0 ) { 
          maximum = *std::max_element(vals.begin(), vals.end());
          minimum = *std::min_element(vals.begin(), vals.end());
        } else { 
          maximum = 0.0;
          minimum = 0.0; 
        };

        std::cout << "- Information of the Coulomb Interaction" << std::endl;
        std::cout << "\t- Total elements: " << this->size() << std::endl;
        std::cout << "\t- Minimum value (non-zero): " << minimum << std::endl; 
        std::cout << "\t- Maximum value (non-zero): " << maximum << std::endl;
        std::cout << "\t- Dimensions              : " << dimensions[0] << "," 
                  << dimensions[1] << "," << dimensions[2] << std::endl;
      };
      void print() const { 
        std::cout << "- Coulomb Interaction elements. Using the conversion (i,j,k,value)" << std::endl;
        this->print_elements(); 
      };
    
    }; // end struct CoulombInteraction
    
    struct Amplitude : public tensor3d_t {
      Amplitude() : tensor3d_t() {}; 
      Amplitude(const std::initializer_list<size_t> &dimensions): tensor3d_t(dimensions) {};
      void info() { 
        double maximum, minimum; 
        auto vals = this->values; 
        vals.erase(std::remove(vals.begin(), vals.end(), 0),vals.end());
        vals.shrink_to_fit();

        if ( vals.size() != 0 ) { 
          maximum = *std::max_element(vals.begin(), vals.end());
          minimum = *std::min_element(vals.begin(), vals.end());
        } else { 
          maximum = 0.0;
          minimum = 0.0; 
        };

        std::cout << "- Information of the Amplitudes" << std::endl;
        std::cout << "\t- Total elements: " << this->size() << std::endl;
        std::cout << "\t- Minimum value (non-zero): " << minimum << std::endl; 
        std::cout << "\t- Maximum value (non-zero): " << maximum << std::endl;
        std::cout << "\t- Dimensions              : " << dimensions[0] << "," 
                  << dimensions[1] << "," << dimensions[2] << std::endl;
      };
      
      void print() const { 
        std::cout << "- Amplitude elements. Using the conversion (i,j,k,value)" << std::endl;
        this->print_elements(); 
      };

    
    }; // end struct CoulombInteraction
    
  
  }; // end namespace concepts
}; // end namespace ueg 

