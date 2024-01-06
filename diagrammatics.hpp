#include "constructors.hpp"
#include "fundamentals.hpp"

namespace ueg {
  namespace diagrammatics {
    typedef ueg::constructors::System system_t;
    typedef ueg::fundamentals::Vector<double> dVec_t;
    typedef ueg::fundamentals::Vector<int> iVec_t;
    typedef ueg::fundamentals::Tensor2d tensor2d_t;
    
    struct HartreeFock : public system_t {
      size_t n_occupied, n_virtuals;
      std::vector<double> en_occ, en_vir;

      // Constructor
      HartreeFock(const system_t &system, std::string system_description) : system_t(system) {
        n_occupied = occupied.size(); n_virtuals = virtuals.size();
        en_occ.reserve(n_occupied); en_vir.reserve(n_virtuals);

        if (system_description == "free") {
          for (size_t index(0); index < n_occupied; ++index) en_occ.push_back(getKineticOccupied(index));
          for (size_t index(0); index < n_virtuals; ++index) en_vir.push_back(getKineticVirtual(index));
        }
        else if (system_description == "interacting") {
          for (size_t index(0); index < n_occupied; ++index) en_occ.push_back(getEnergyOccupied(index));
          for (size_t index(0); index < n_virtuals; ++index) en_vir.push_back(getEnergyVirtual(index));
        }
      };
        
      // Calculates the Hartree-Fock eigenenergy of occupied and virtual orbitals
      // using Alavi's & Spencer's or madelung's schemes.
      // @params index: the corresponding index in 'occupied' orbitals
      double getEnergyOccupied(const size_t &index) const;
      double getEnergyVirtual(const size_t &index) const;
 
      // Get free-electron eigenenergies of occupied and virtual orbitals
      double getKineticOccupied(const size_t &index) const;
      double getKineticVirtual(const size_t &index) const;

      // Calculate the Propagator (total, kinetic or exchange) given the indices.
      double getPropagator(const size_t &i, const size_t &j, const size_t &a, const size_t &b) const;
      
      // Calculate th Hartree-Fock energy of your system.
      void getHFenergy() const;

      // Trial function to calculate Triple's Propagator of a certain diagram.
      double getTriplesPropagator(const size_t &i, const size_t &j, const size_t &k, const size_t &a, const size_t &c, const size_t &d) const;
        
      // Returns the coulomb interaction in the virtual or occupied orbitals' subspace
      tensor2d_t getCoulombInteraction(std::string subspace) const;

      // Returns the coulomb interaction in virtual orbitals' subspace in a matrix
      tensor2d_t getCoulombInteractionVirtuals() const;

      // Returns the coulomb interaction in occupied orbitals' subspace in a matrix
      tensor2d_t getCoulombInteractionOccupied() const;

      // Returns the coulomb interaction between the occupied and virtuals in a matrix
      tensor2d_t getCoulombInteractionMixed() const;

      // Calculates the coulomb interaction in virtual orbitals' subspace on the fly
      double getCoulombInteractionVirtuals(const size_t &k, const size_t &l) const;

      // Calculates the coulomb interaction in occupied orbitals' subspace on the fly
      double getCoulombInteractionOccupied(const size_t &k, const size_t &l) const; 

      // For more expressive and easy to read-code (helpers inline == no overhead)
      inline double getVijji(const dVec_t &vec) const;
      inline double getViiii() const;

    }; // end of class HartreeFock
  
  }; // end of namespace diagrammatics
}; // end of namespace ueg
