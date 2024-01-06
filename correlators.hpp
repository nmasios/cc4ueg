#pragma once
#include "fundamentals.hpp"
//#include "./constructors.hpp"
#include "concepts.hpp"
#include "diagrammatics.hpp"
#include <fstream>
#include <string>
#include <numeric>
#include <deque>
#include <chrono>
#include <armadillo>
//#include "./utility.hpp"


namespace ueg {
  namespace correlators {
    // Type-definitions for simplicity
    typedef ueg::constructors::System system_t;
    typedef ueg::fundamentals::Vector<int> iVec_t;
    typedef ueg::fundamentals::Vector<double> dVec_t;
    typedef ueg::fundamentals::Tensor2d tensor2d_t;
    typedef ueg::fundamentals::Tensor5d tensor5d_t;
    typedef std::map<iVec_t, std::vector<std::pair<size_t, size_t>>> mapping_t;
    typedef std::map<iVec_t, std::vector<std::tuple<size_t, size_t, size_t>>> mapping_3d_t;
    typedef std::vector<std::tuple<double,double,double,double>> vector_dtuples_t;
    typedef std::deque<ueg::concepts::Amplitude> deque_t;
    typedef std::tuple<ueg::concepts::Amplitude, ueg::concepts::Amplitude> tensor_tuple_t;

    // Definition of the Correlator
    template<class StartingPoint>
    struct Correlator {
      StartingPoint refpoint;                                            // The starting point of the diagrammatics
      ueg::concepts::Residuum residuum;                                  // The residuum structure
      ueg::concepts::Amplitude amplitude;                                // The amplitude (hole-hole-particle) structure
      ueg::concepts::Propagator propagator;                              // The propagator structure
      ueg::concepts::Propagator propagator_kin;                          // The kinetic propagator structure
      ueg::concepts::CoulombInteraction interaction;                     // The interaction structure
      tensor2d_t virtualSpaceCoulomb, occupiedSpaceCoulomb;              // Coulomb interaction in a subspace
      tensor2d_t mixedSpaceCoulomb;                                      // Coulomb interaction between occupied and virtuals
      tensor5d_t triple_residuum, triple_conjugate;                      // The triples' residuum and conjugate residuum structure
      mapping_t mappingPHvertex, mappingPP, mappingHH;                   // PH, PP and HH mappings
      mapping_3d_t mappingHHH, mappingPPP;                               // PPP and HHH mappings for the triples
      deque_t stored_amplitudes, stored_errors;                          // Stores a successive sequence of amplitudes and errors
      bool verbosity = true;                                             // If True then 'message' function is enabled...
      double kFermi, q_min, q_max, lattice_constant;                     // Minimum-maximum momentum transfer and kFermi
  
      // Constructor
      Correlator(const system_t &system, std::string system_description);
      
      // Routine to initialize the interaction and the propagator
      void initialize();
      
      // Routine to request the evaluation of the energy (direct, exchange or total)
      double getEnergy(std::string type) const;

      // Routine to request the evaluation of the energy (direct, exchange or total)
      double getDecomposedEnergy(ueg::concepts::Amplitude received_amplitude, std::string type) const;

      // Routine to evaluate the Direct energy of the system
      double getDirectEnergy() const;

      // Routine to evaluate the Direct energy of the system
      double getDecomposedDirectEnergy(ueg::concepts::Amplitude received_amplitude) const;
      
      // Routine to evaluate the Exchange energy of the system
      double getExchangeEnergy() const;

      // Routine to evaluate the Exchange energy of the system
      double getDecomposedExchangeEnergy(ueg::concepts::Amplitude received_amplitude) const;
      
      // Routine to symmetrize the 'residuum' structure
      void symmetrize();

      // Routine to request a specific 'diagram' to be added into the residuum
      void add_to_residuum(std::string diagram);

      // Routine to request a sequence of 'diagrams' to be added into the residuum
      void add_to_residuum(std::vector<std::string> diagrams);
  
      // Adds the Direct Ring Coupled Cluster into the residuum
      void _direct_ring_cc();
      void _rings();
      void _ring_quadratic();

      // Adds ladder diagrams into the residuum
      void _pp_ladder();
      void _hh_ladder();
      void _ph_vertex_ladder();
      void _ph_vertices_ladder();

      // Adds the Xijkl_Tklab quadratic term into the residuum
      void _Xijkl_Tklab_quadratic();
     
      // Adds the Xciak_Tkjcb quadratic term into the residuum
      void _Xciak_Tkjcb_quadratic();

      // Adds the Xcibk_Tkjac quadratic term into the residuum
      void _Xcibk_Tkjac_quadratic();

      // Adds the Fock term into the residuum
      void _F_term();
      
      // Adds the lambda terms into the residuum
      void _Lca_Tijcb();
      void _Lik_Tkjab();

      // Adds Xicak_Tkjcb terms into the residuum.
      // Term: "all", "linear", "quadratics", "non antisymmetrized integeral", "antisymmetrized integral"
      void _Xicak_Tkjcb();
      void _Xicak_Tkjcb_linear();
      void _Xicak_Tkjcb_non_antisym();
      void _Xicak_Tkjcb_antisym();
      void _Xicak_Tkjcb_antisym_2();
      void _Xicak_Tkjcb_quadratics();
      
      // Adds Xicak_Tkjbc terms into the residuum.
      // Terms: "all", "linear", "quadratics", "non antisymmetrized integral", "antisymmetrized integral"
      void _Xicak_Tkjbc();
      void _Xicak_Tkjbc_linear();
      void _Xicak_Tkjbc_non_antisym();
      void _Xicak_Tkjbc_antisym();
      void _Xicak_Tkjbc_quadratics();

      // Calculates the Triples energy with PPP and HHH momentum conservation maps.
      // The user can choose between direct, single/double exchange, or all terms.
      // It contains two different implementations for the sake of condintentiality.
      void _T(std::string option);

      // Calculates the Complete Triples energy with PPP and HHH
      // momentum conservation maps and the denominator from Piecuch's paper.
      void _cT();
 
      // Calculates the (T) energy and structure factor based on a triangular 
      // algorithm for memory efficiency (On the fly) and returns the SF 
      // results in a vector of tuples of doubles as (momentum transfer q,
      // triples structure factor, FT Coulomb interaction).
      vector_dtuples_t _OTF_CCDT(bool structure_factor);

      // Calculates both (T) and (cT) energies and structure factors and
      // returns the results in a vector of tuples of doubles as (momentum transfer q,
      // triples structure factor, cr-triples structure factor, FT Coulomb interaction).
      vector_dtuples_t _OTF_CCDcT(bool structure_factor);

      // Routine that composes the body of a one shot calculation (MP2/MP3/MP4)
      ueg::concepts::Amplitude one_shot_methods(std::string method);
      
      // Routine that composes the body of a full self-consistency cycle
      ueg::concepts::Residuum self_consistent_cycle_body(std::string method);
      
      // Calculates the amplitude, from the residuum, and the error as the difference of the current and preceding amplitudes
      tensor_tuple_t calculate_errors_and_amplitudes();

      // Returns the extrapolated amplitude using the DIIS technique
      ueg::concepts::Amplitude diis_extrapolation(tensor_tuple_t amplitude_and_error);

      // Returns the decomposed amplitude tensor for the optional and selected channel decomosition.
      ueg::concepts::Amplitude decomposition(std::string diagram);
      
      // Calculates the structure factor and returns the results in a vector of
      // tuples of doubles as (momentum transfer q, direct structure factor
      // exchange structure factor, FT Coulomb interaction).
      vector_dtuples_t structure_factor(ueg::concepts::Amplitude amplitude_received);
      
      // Routine to evaluate the Coulomb Integral
      inline double coulombIntegral(const iVec_t &vec) const {
        return (4*M_PI/( ((2*M_PI/this->refpoint.lattice.constant) * vec).squaredLength() )/refpoint.lattice.boxVolume);
      };

    }; // end class Correlator

  }; // end namespace correlators
}; // end namespace ueg

#include "./correlators.ipp"

