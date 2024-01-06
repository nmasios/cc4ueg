#pragma once
#include "fundamentals.hpp"
#include "utility.hpp"
#include <exception>
#include <stdexcept>
#include <map>
#include <limits>
#include <algorithm>
#include <string>
#include <iostream>

namespace ueg {
  namespace constructors {

    // Type-definitions for code readability
    typedef ueg::fundamentals::Vector<double> dVec_t;
    typedef ueg::fundamentals::Vector<int> iVec_t;
    typedef std::map<iVec_t, std::vector<std::pair<size_t,size_t>>> mapping_t;
    typedef std::map<iVec_t, std::vector<std::tuple<size_t,size_t,size_t>>> mapping_3d_t;
    typedef std::tuple<dVec_t, dVec_t, dVec_t, dVec_t, dVec_t, dVec_t, double> lattice_properties_t;

    // The routine is called to retrieve the definitions of provided Lattices
    lattice_properties_t lattice_definition(std::string lattice_type);

    // The class defines the abstraction of a Lattice, together with
    // all the properties that comes with.
    struct Lattice
    {

      dVec_t a1,a2,a3;              // Primitive vectors
      dVec_t b1,b2,b3;              // Reciprocal vectors
      std::pair<int,int> grid;      // Grid limits along an axis
      std::string type;             // Bravais lattice type
      double volume;                // Volume of the unit-cell
      double constant = 1;          // Lattice constant
      double boxVolume = 1;         // Volume of the box L^3
      double radius = std::numeric_limits<double>::infinity();
      
      // Constructors + Copy
      Lattice() {};
      Lattice(const Lattice &other);
      Lattice(std::string celltype, int gridAxisSize);
      
      // Routine to construct all possible states (shifted or not) within a spherical shell.
      // The shift can be generated from random numbers and is a double array of three elements (x,y,z-directions).
      std::vector<dVec_t> getAllStates(double normalization, dVec_t shift = {0.0, 0.0, 0.0}) const;
  
      // Routine to reset the lattice constant and recalculate the boxVolume
      void reset_constant(size_t n_electrons, double wigner_seitz);
      void reset_constant(double value);
    
    }; // end class Lattice
    
    // The class defines the fundamental properties of the problem at hand
    // and will be passed in all structures such as Correlators or defined
    // diagrammatic methodologies such as Hartree-Fock.
    struct System
    {
      Lattice lattice;                // The lattice of the system
      std::vector<dVec_t> occupied;   // The occupied orbitals
      std::vector<dVec_t> virtuals;   // The virtual orbitals
      double wigner_seitz_radius;     // The Wigner-Seitz radius
      double madelung;                // The Madelung constant
      dVec_t shift;                   // The k-mesh shift

      // Constructors (full or partially interactive)
      System();
      System(std::string celltype, int gridAxisSize);
      System(std::string celltype, size_t numVirtuals, size_t numElectrons, double WignerSeitzRadius, dVec_t shift);
      System(std::string celltype, int gridAxisSize, size_t numElectrons, double WignerSeitzRadius);
      System(const Lattice &lat, size_t numElectrons, double WignerSeitzRadius);
      System(const System &sys);
      
      // Routine to renormalize and return either the "virtuals" or the "occupied" states
      // into their integer form.
      std::vector<iVec_t> transform_to_integer_notation(std::string orbital_type) const;
      
      // Routine to create the particle-hole and hole-hole 'mappings'
      mapping_t get_ph_mapping() const;
      mapping_t get_pp_mapping() const;
      mapping_t get_hh_mapping() const;
      mapping_3d_t get_hhh_mapping() const;
      mapping_3d_t get_ppp_mapping() const;
      
      // Routine that obtains the minimum and maximum momentum transfer of the given system. 
      std::tuple<double,double> get_min_max_momentum_transfer() const;
      
      // Routines to evaluate (or converge) the madelung constant
      inline double get_madelung_constant(int limit) const;
      double calculate_madelung_constant(double threshold=1e-5) const;

      // Print out basic information
      void information() const;

    }; // end class System
    
    // Simple routine to return the available electron filling given the grid and the lattice-type
    std::vector<size_t> lattice_available_electron_filling(std::string celltype, int gridAxisSize);

  }; // end namespace constructors
}; // end namespace ueg
