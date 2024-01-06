#pragma once
#include<iterator>
#include<armadillo>
#include<valarray>
#include <chrono>

namespace ueg {
  namespace correlators {
    // ------------------------------------------------------------------------------------------------------
    //                               IMPLEMENTATION OF THE CORRELATOR STRUCTURE
    // ------------------------------------------------------------------------------------------------------
    template<class StartingPoint> Correlator<StartingPoint>::Correlator(const system_t &system, std::string system_description) : refpoint(system, system_description){
      // Construction of the basic concepts
      interaction = ueg::concepts::CoulombInteraction({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals});
      propagator = ueg::concepts::Propagator({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals});
      residuum = ueg::concepts::Residuum({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals});
      amplitude = ueg::concepts::Amplitude({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals});
      //triple_residuum = ueg::fundamentals::Tensor5d({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});
      //triple_conjugate = ueg::fundamentals::Tensor5d({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});
      interaction.clear(); propagator.clear(); residuum.clear(); amplitude.clear(); 
      //triple_residuum.clear(); triple_conjugate.clear();

      // Get the PH and HH mappings (from the reference point)
      mappingPHvertex = refpoint.get_ph_mapping();
      mappingPP = refpoint.get_pp_mapping();
      mappingHH = refpoint.get_hh_mapping();

      // Get the HHH and PPP mappings (from the reference point)
      //mappingHHH = refpoint.get_hhh_mapping();
      //mappingPPP = refpoint.get_ppp_mapping();

      // Get the minimum-maximum momentum transfer and kFermi (not sure if it's working)
      //std::tie(q_min,q_max) = refpoint.get_min_max_momentum_transfer();
      //kFermi = refpoint.occupied.back().norm();

      // Construct the coulomb interaction in the subspaces (occupied + virtuals)
      occupiedSpaceCoulomb = refpoint.getCoulombInteraction("occupied");
      mixedSpaceCoulomb    = refpoint.getCoulombInteraction("mixed");
    }
      
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::initialize() {


      for (const auto [kvector, connections] : this->mappingPHvertex) 
        for (auto const &[i,a] : connections)
          for (auto const &[j,b] : this->mappingPHvertex[-kvector]) {
            this->interaction(i,j,a) = this->mixedSpaceCoulomb(i,a);
            this->propagator(i,j,a) = this->refpoint.getPropagator(i,j,a,b); }

      ueg::utility::message("\t- Interaction and Propagator are setup", this->verbosity);
    }
      
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> double Correlator<StartingPoint>::getEnergy(std::string type) const {
      if (type == "total") return this->getDirectEnergy() + this->getExchangeEnergy();
      else if (type == "direct") return this->getDirectEnergy();
      else return this->getExchangeEnergy();
    }

    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> double Correlator<StartingPoint>::getDecomposedEnergy(ueg::concepts::Amplitude received_amplitude, std::string type) const {
      if (type == "total") return this->getDecomposedDirectEnergy(received_amplitude) + this->getDecomposedExchangeEnergy(received_amplitude);
      else if (type == "direct") return this->getDecomposedDirectEnergy(received_amplitude);
      else return this->getDecomposedExchangeEnergy(received_amplitude);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> double Correlator<StartingPoint>::getDirectEnergy() const {
      double energy = 0.0;
      for (size_t index(0); index < amplitude.size(); ++index)
        energy += amplitude.values[index] * interaction.values[index];
      energy *= std::pow(-1,2) * std::pow(-2,2) / 2.0 / (2 * refpoint.n_occupied);

      if (verbosity) {
        std::cout.precision(12);
        std::cout << "- Direct energy in [Ha/N]  : " << std::fixed << energy << std::endl;
      }
      return energy;
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> double Correlator<StartingPoint>::getDecomposedDirectEnergy(ueg::concepts::Amplitude received_amplitude) const {
      double energy = 0.0;
      for (size_t index(0); index < received_amplitude.size(); ++index)
        energy += received_amplitude.values[index] * interaction.values[index];
      energy *= std::pow(-1,2) * std::pow(-2,2) / 2.0 / (2 * refpoint.n_occupied);

      if (verbosity) {
        std::cout.precision(12);
        std::cout << "- Direct energy in [Ha/N]  : " << std::fixed << energy << std::endl;
      }
      return energy;
    }
    
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> double Correlator<StartingPoint>::getExchangeEnergy() const {
      double energy = 0.0;
      for (size_t a(0); a < refpoint.n_virtuals; ++a)
        for (size_t j(0); j < refpoint.n_occupied; ++j)
          for (size_t i(0); i < refpoint.n_occupied; ++i)
            energy += this->amplitude(i,j,a) * this->interaction(j,i,a);
       
      energy *= std::pow(-1,2) * std::pow(-2,1) / 2.0 / (2*refpoint.n_occupied);
      
      if (verbosity) {
        std::cout.precision(12);
        std::cout << "- Exchange energy in [Ha/N]: " << std::fixed << energy << std::endl;
      }
      return energy;
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> double Correlator<StartingPoint>::getDecomposedExchangeEnergy(ueg::concepts::Amplitude received_amplitude) const {
      double energy = 0.0;

      for (size_t a(0); a < refpoint.n_virtuals; ++a)
        for (size_t j(0); j < refpoint.n_occupied; ++j)
          for (size_t i(0); i < refpoint.n_occupied; ++i)
            energy += received_amplitude(i,j,a) * this->interaction(j,i,a);
       
      energy *= std::pow(-1,2) * std::pow(-2,1) / 2.0 / (2*refpoint.n_occupied);
      
      if (verbosity) {
        std::cout.precision(12);
        std::cout << "- Exchange energy in [Ha/N]: " << std::fixed << energy << std::endl;
      }
      return energy;
    }
      
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::symmetrize() {
      auto tempres = this->residuum;
      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[j,b] : mappingPHvertex[-kvector])
          for (auto const &[i,a] : connections)
            residuum(i,j,a) = tempres(i,j,a) + tempres(j,i,b);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::add_to_residuum(std::string diagram) {
      if      (diagram == "Rings")                         this->_rings();
      else if (diagram == "Ring-Quadratic")                this->_ring_quadratic();
      else if (diagram == "Direct-Ring-CC")                this->_direct_ring_cc();
      else if (diagram == "PP-ladder")                     this->_pp_ladder();
      else if (diagram == "HH-ladder")                     this->_hh_ladder();
      else if (diagram == "PH-vertex-ladder")              this->_ph_vertex_ladder();
      else if (diagram == "PH-vertices-ladder")            this->_ph_vertices_ladder();
      else if (diagram == "Xijkl_Tklab-quadratic")         this->_Xijkl_Tklab_quadratic();
      else if (diagram == "Xciak_Tkjcb-quadratic")         this->_Xciak_Tkjcb_quadratic();
      else if (diagram == "Xcibk_Tkjac-quadratic")         this->_Xcibk_Tkjac_quadratic();
      else if (diagram == "Lca_Tijcb")                     this->_Lca_Tijcb();
      else if (diagram == "Lik_Tkjab")                     this->_Lik_Tkjab();
      else if (diagram == "Xicak_Tkjcb")                   this->_Xicak_Tkjcb();
      else if (diagram == "Xicak_Tkjcb-linear")            this->_Xicak_Tkjcb_linear();
      else if (diagram == "Xicak_Tkjcb-non-antisym")       this->_Xicak_Tkjcb_non_antisym();
      else if (diagram == "Xicak_Tkjcb-antisym")           this->_Xicak_Tkjcb_antisym();
      else if (diagram == "Xicak_Tkjcb-antisym-2")         this->_Xicak_Tkjcb_antisym_2();
      else if (diagram == "Xicak_Tkjcb-quadratics")        this->_Xicak_Tkjcb_quadratics();
      else if (diagram == "Xicak_Tkjbc")                   this->_Xicak_Tkjbc();
      else if (diagram == "Xicak_Tkjbc-linear")            this->_Xicak_Tkjbc_linear();
      else if (diagram == "Xicak_Tkjbc-non-antisym")       this->_Xicak_Tkjbc_non_antisym();
      else if (diagram == "Xicak_Tkjbc-antisym")           this->_Xicak_Tkjbc_antisym();
      else if (diagram == "Xicak_Tkjbc-quadratics")        this->_Xicak_Tkjbc_quadratics();
      else throw std::invalid_argument("Unrecognized diagram provided: " + diagram);
    }

    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::add_to_residuum(std::vector<std::string> diagrams) {
      for (auto diagram : diagrams) this->add_to_residuum(diagram);
    }

    // ------------------------------------------------------------------------------------------------------
    //        IMPLEMENTATION OF THE CLOSED-SHELL ENERGY EXPRESSION AND AMPLITUDE EQUATIONS FOR THE UEG
    //                      HIRATA et al. CHEMICAL PHYSICS LETTERS 345 (2001) 475-480
    // ------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------
    //                                       DIRECT RING COUPLED-CLUSTER
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::_rings() {

      double contraction = 0.0, Vicak = 0.0;
      std::vector<double> summation(refpoint.n_occupied, 0.0);

      for (auto const [kvector, connections] : mappingPHvertex) {
        Vicak = 2.0 * this->coulombIntegral(kvector);

        for (size_t j(0); j < refpoint.n_occupied; ++j) {
          contraction = 0.0;
          for (auto const &[k,c] : connections)
            contraction +=  this->amplitude(k,j,c);
          summation[j] =  contraction * Vicak;
        }

        for (auto const &[i,a] : connections)
          for (size_t j(0); j < refpoint.n_occupied; ++j)
            this->residuum(i,j,a) += summation[j];
      }
      ueg::utility::message("\t- Ring term is added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_ring_quadratic() {
      auto tempres = residuum; tempres.clear();
      double contraction = 0.0, Vilad = 0.0;
      std::vector<double> summation(refpoint.n_occupied, 0.0);

      // Calculate the contraction
      for (auto const [kvector, connections] : mappingPHvertex) {
        Vilad = std::sqrt(this->coulombIntegral(kvector));
        
        for (size_t j(0); j < refpoint.n_occupied; ++j) {
          contraction = 0.0;
          for (auto const &[l,d] : connections)
            contraction += this->amplitude(l,j,d);
          summation[j] = Vilad * contraction;
        }

        for (auto const &[i,a] : connections)
          for (size_t j(0); j < refpoint.n_occupied; ++j)
            tempres(i,j,a) += summation[j];
      }

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[j,b] : mappingPHvertex[-kvector])
          for (auto const &[i,a] : connections)
            this->residuum(i,j,a) += 2.0 * tempres(i,j,a) * tempres(j,i,b);
      
      ueg::utility::message("\t- Quadratic-ring is added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_direct_ring_cc() {
      _rings();
      _ring_quadratic();
      ueg::utility::message("\t- Direct Ring Coupled Cluster terms added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    //                                    IMPLEMENTATION OF LADDER DIAGRAMS
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::_pp_ladder() {
      
      double Vabcd_inner_shell, Vabcd_outer_shell, Vabcd_inner_outer, Vabcd_outer_inner;
      double inner_shell_amplitude, outer_shell_amplitude;
      size_t Nv = refpoint.n_virtuals;

      for (size_t a(0); a < refpoint.n_virtuals/2; ++a) 
        for (size_t c(0); c < refpoint.n_virtuals/2; ++c) {
           Vabcd_inner_shell = this-> refpoint.getCoulombInteractionVirtuals(a,c);
           Vabcd_inner_outer = this-> refpoint.getCoulombInteractionVirtuals(a,Nv/2+c);
           Vabcd_outer_inner = this-> refpoint.getCoulombInteractionVirtuals(Nv/2+a,c);
           Vabcd_outer_shell = this-> refpoint.getCoulombInteractionVirtuals(Nv/2+a,Nv/2+c);
           for (size_t i(0); i < refpoint.n_occupied; ++i) 
             for (size_t j(0); j < refpoint.n_occupied; ++j) {
               inner_shell_amplitude = amplitude(i,j,c);
               outer_shell_amplitude = amplitude(i,j,Nv/2+c);
               this->residuum(i,j,a) += inner_shell_amplitude * Vabcd_inner_shell;
               this->residuum(i,j,a) += outer_shell_amplitude * Vabcd_inner_outer;
               this->residuum(i,j,Nv/2+a) += inner_shell_amplitude * Vabcd_outer_inner;
               this->residuum(i,j,Nv/2+a) += outer_shell_amplitude * Vabcd_outer_shell;
            }
        }

      ueg::utility::message("\t- Particle-Particle ladder diagram added into the Residuum", this->verbosity);
    }
    
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_hh_ladder() {

      for (auto const [kvector, connections] : mappingHH) 
        for (auto const &[i,k] : connections) 
          for (auto const &[j,l] : mappingHH[-kvector]) 
            for (size_t a(0); a < refpoint.n_virtuals; ++a) 
              this->residuum(i,j,a) += amplitude(k,l,a) * this->occupiedSpaceCoulomb(i,k); 

      ueg::utility::message("\t- Hole-Hole ladder diagram added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_ph_vertex_ladder() {
      
      for (auto const [kvector, connections ] : mappingPHvertex) 
        for (auto const &[i,a] : connections) 
          for (auto const &[j,b] : mappingPHvertex[-kvector]) 
            for (auto const &[k,c] : mappingPHvertex[-kvector]) 
              this->residuum(i,j,a) += -amplitude(i,k,a) * this->occupiedSpaceCoulomb(j,k); 

      ueg::utility::message("\t- Particle-Hole ladder diagram added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::_ph_vertices_ladder() {
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      
      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[j,b] : mappingPHvertex[-kvector])
            for (auto const &[k,c] : mappingPHvertex[occ[i]-vir[b]])
              this->residuum(i,j,a) += -amplitude(k,j,a) * this->occupiedSpaceCoulomb(k,i);

      ueg::utility::message("\t- Particle-hole ladder diagram added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    //                             IMPLEMENTATION OF TERM Xijkl_Tklab (Quadratic)
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xijkl_Tklab_quadratic() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_occupied, refpoint.n_occupied});

      for (size_t i(0); i < refpoint.n_occupied; ++i)
        for (size_t j(0); j < refpoint.n_occupied; ++j)
          for (size_t k(0); k < refpoint.n_occupied; ++k)
            for (size_t c(0); c < refpoint.n_virtuals; ++c)
              chi(i,j,k) += amplitude(i,j,c) * this->mixedSpaceCoulomb(k,c);

      for (auto const [kvector, connections] : mappingHH)
        for (auto const &[i,k] : connections)
          for (auto const &[j,l] : mappingHH[-kvector])
            for (size_t a(0); a < refpoint.n_virtuals; ++a)
              this->residuum(i,j,a) += amplitude(k,l,a) * chi(i,j,k);

      ueg::utility::message("\t- Xijkl_Tklab Quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    //                             IMPLEMENTATION OF TERM Xciak_Tkjcb (Quadratic)
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xciak_Tkjcb_quadratic() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,d] : connections)
          for (auto const &[l,c] : mappingPHvertex[-kvector])
            for ( auto const &[i,a] : mappingPHvertex[vir[d]-occ[l]])
              chi(i,c,k) += amplitude(i,l,d) * this->mixedSpaceCoulomb(k,d);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,c] : connections)
          for (auto const &[j,b] : mappingPHvertex[-kvector])
            for (auto const &[i,a] : connections)
              this->residuum(i,j,a) += 0.5 * amplitude(k,j,c) * chi(i,c,k);

      ueg::utility::message("\t- Xciak_Tkjcb Quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    //                             IMPLEMENTATION OF TERM Xcibk_Tkjac (Quadratic)
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xcibk_Tkjac_quadratic() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,d] : connections)
          for (auto const &[l,c] : mappingPHvertex[-kvector])
            for (auto const &[i,b] : mappingPHvertex[vir[d]-occ[l]])
              chi(i,c,k) += amplitude(i,l,d) * this->mixedSpaceCoulomb(k,d);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[j,c] : connections)
          for (auto const &[k,a] : mappingPHvertex[-kvector])
            for (auto const &[i,b] : mappingPHvertex[occ[k]-vir[c]])
              this->residuum(i,j,a) += 0.5 * amplitude(j,k,c) * chi(i,c,k);

      ueg::utility::message("\t- Xcibk_Tkjac Quadratic term added into the Residuum", this->verbosity);
    }


    // ------------------------------------------------------------------------------------------------------
    //                                       IMPLEMENTATION OF λ TERMS
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::_Lca_Tijcb() {
      std::vector<double> kappa_ca(this->refpoint.n_virtuals, 0.0);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,d] : connections)
          for (auto const &[k,c] : mappingPHvertex[-kvector])
            kappa_ca[c] += -amplitude(k,l,c) * (2.0*this->mixedSpaceCoulomb(l,d) - this->mixedSpaceCoulomb(l,c));

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[j,b] : mappingPHvertex[-kvector])
            this->residuum(i,j,a) += amplitude(i,j,a) * kappa_ca[a];
      
      ueg::utility::message("\t- Λca_Tijcb term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Lik_Tkjab() {
      std::vector<double> kappa_ik(this->refpoint.n_occupied, 0.0);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,d] : connections)
          for (auto const &[i,c] : mappingPHvertex[-kvector])
            kappa_ik[i] += amplitude(i,l,c) * (2.0*this->mixedSpaceCoulomb(l,d) - this->mixedSpaceCoulomb(l,c));
      
      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[j,b] : mappingPHvertex[-kvector])
            this->residuum(i,j,a) += -amplitude(i,j,a) * kappa_ik[i];
      
      ueg::utility::message("\t- Λik_Tkjab term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------
    //                                   IMPLEMENTATION OF TERM Xicak_Tkjcb
    // ------------------------------------------------------------------------------------------------------
    
    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb_linear() {

      double contraction = 0.0, Vicak = 0.0;
      std::vector<double> summation(refpoint.n_occupied, 0.0);

      for (auto const [kvector, connections] : mappingPHvertex) {
        Vicak = 2.0 * this->coulombIntegral(kvector);

        for (size_t j(0); j < refpoint.n_occupied; ++j) {
          contraction = 0.0;
          for (auto const &[k,c] : connections)
            contraction +=  this->amplitude(k,j,c);
          summation[j] =  contraction * Vicak;
        }

        for (auto const &[i,a] : connections)
          for (size_t j(0); j < refpoint.n_occupied; ++j)
            this->residuum(i,j,a) += summation[j];
      }
      ueg::utility::message("\t- Xicak_Tkjcb linear term is added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb_non_antisym() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,c] : connections)
          for (auto const &[l,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : connections)
              chi(c,j,k) += amplitude(l,j,b) * this->mixedSpaceCoulomb(k,c);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[k,c] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[-kvector])
              this->residuum(i,j,a) += -amplitude(i,k,a) * chi(c,j,k);

      ueg::utility::message("\t- Xicak_Tkjcb non-antisymmetrized quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb_antisym() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,c] : connections)
          for (auto const &[k,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[vir[d]-occ[l]])
              chi(c,j,k) += amplitude(l,j,d) * (2.0 * this->mixedSpaceCoulomb(k,c) - this->mixedSpaceCoulomb(l,c));

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[k,c] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[-kvector])
              this->residuum(i,j,a) += amplitude(i,k,a) * chi(c,j,k);
      
      ueg::utility::message("\t- Xicak_Tkjcb antisymmetrized quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb_antisym_2() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,c] : connections)
          for (auto const &[k,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[vir[d]-occ[l]])
              chi(c,j,k) += -amplitude(l,j,d) * this->mixedSpaceCoulomb(l,c);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[k,c] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[-kvector])
              this->residuum(i,j,a) += amplitude(i,k,a) * chi(c,j,k);
      
      ueg::utility::message("\t- Xicak_Tkjcb antisymmetrized quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb_quadratics() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,c] : connections)
          for (auto const &[k,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[vir[d]-occ[l]])
              chi(c,j,k) += amplitude(l,j,d) * (2.0 * this->mixedSpaceCoulomb(k,c) - this->mixedSpaceCoulomb(l,c))
                             - amplitude(l,j,b) * this->mixedSpaceCoulomb(k,c);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[k,c] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[-kvector])
              this->residuum(i,j,a) += amplitude(i,k,a) * chi(c,j,k);
      
      ueg::utility::message("\t- Xicak_Tkjcb quadratic terms added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjcb() {
      _Xicak_Tkjcb_linear(); 
      _Xicak_Tkjcb_quadratics();
      ueg::utility::message("\t- Xicak_Tkjcb terms added into the Residuum", this->verbosity);
    }

  
    // ------------------------------------------------------------------------------------------------------
    //                                   IMPLEMENTATION OF TERM Xicak_Tkjbc
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjbc_linear() {
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      
      for (auto const [kvector, connections] : mappingPHvertex) 
        for (auto const &[i,a] : connections) 
          for (auto const &[k,c] : connections) 
            for (auto const &[j,b] : mappingPHvertex[vir[c]-occ[k]]) 
              this->residuum(i,j,a) += -amplitude(k,j,b) * this->mixedSpaceCoulomb(i,a);

      ueg::utility::message("\t- Xicak_Tkjbc linear term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjbc_non_antisym() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,c] : connections)
          for (auto const &[l,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : connections)
              chi(c,j,k) += amplitude(l,j,b) * this->mixedSpaceCoulomb(k,c);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,c] : connections)
          for (auto const &[k,a] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[occ[k]-vir[c]])
              this->residuum(i,j,a) += 0.5 * amplitude(i,k,c) * chi(c,j,k);

      ueg::utility::message("\t- Xicak_Tkjbc non-antisymmetrized quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjbc_antisym() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[l,c] : connections)
          for (auto const &[k,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[vir[d]-occ[l]])
              chi(c,j,k) += amplitude(l,j,d) * (0.5 * this->mixedSpaceCoulomb(l,c) - this->mixedSpaceCoulomb(l,d));

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,c] : connections)
          for (auto const &[k,a] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[occ[k]-vir[c]])
              this->residuum(i,j,a) += amplitude(i,k,c) * chi(c,j,k);

      ueg::utility::message("\t- Xicak_Tkjbc antisymmetrized quadratic term added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjbc_quadratics() {
      auto chi = ueg::fundamentals::Tensor3d({refpoint.n_virtuals, refpoint.n_occupied, refpoint.n_occupied});
      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");
      chi.clear();

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[k,c] : connections)
          for (auto const &[l,d] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : connections)
              chi(c,j,k) += amplitude(l,j,d) * (0.5 * this->mixedSpaceCoulomb(l,c) - this->mixedSpaceCoulomb(l,d))
                             + 0.5 * amplitude(l,j,b) * this->mixedSpaceCoulomb(k,c);

      for (auto const [kvector, connections] : mappingPHvertex)
        for (auto const &[i,c] : connections)
          for (auto const &[k,a] : mappingPHvertex[-kvector])
            for (auto const &[j,b] : mappingPHvertex[occ[k]-vir[c]])
              this->residuum(i,j,a) += amplitude(i,k,c) * chi(c,j,k);

      ueg::utility::message("\t- Xicak_Tkjbc quadratic terms added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> void Correlator<StartingPoint>::_Xicak_Tkjbc() {
      _Xicak_Tkjbc_linear();
      _Xicak_Tkjbc_quadratics();
      ueg::utility::message("\t- Xicak_Tkjbc terms added into the Residuum", this->verbosity);
    }

    // ------------------------------------------------------------------------------------------------------


    // ------------------------------------------------------------------------------------------------------
    //                        IMPLEMENTATION OF (T) and (cT) (normal and on the fly)                                    
    //                             https://doi.org/10.1103/PhysRevLett.131.186401
    // ------------------------------------------------------------------------------------------------------
 
    template<class StartingPoint> void Correlator<StartingPoint>::_T(std::string option) {

      auto occ = this->refpoint.transform_to_integer_notation("occupied");
      auto vir = this->refpoint.transform_to_integer_notation("virtuals");

      // Initialize a hole-particle-particle amplitude to facilitate the contractions.
      auto amplitude_hpp = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});
      amplitude_hpp.clear();

      // Fill up the HPP amplitude with the corresponding values.
      for (const auto [kvector, connections] : this->mappingPHvertex)
        for (auto const &[i,a] : connections)
          for (auto const &[j,b] : this->mappingPHvertex[-kvector])
            amplitude_hpp(i,a,b) = this->amplitude(i,j,a);
      
      // Calculate the (T) residuum by inserting the parrticle-particle and hole-hole parts.
      for (auto const [kvector,connections] : mappingHHH)
      for (auto const &[i,j,k] : connections)
      for (size_t c(0); c < refpoint.n_virtuals; ++c)
      for (auto const &[a,b]: mappingPP[kvector-vir[c]])
        triple_residuum(i,j,k,a,b) +=  ( amplitude(i,j,a) - amplitude_hpp(i,a,b) ) * this->mixedSpaceCoulomb(k,c);
      
      // Compute the energy. The user can distinguish between, direct/single-exchange/double-exchange/all terms.
      double energy(0.0);
      for (auto const [kvector,connections] : mappingHHH)
      for (auto const &[i,j,k] : connections)
      for (size_t c(0); c < refpoint.n_virtuals; ++c)
      for (auto const &[a,b]: mappingPP[kvector-vir[c]]){

            if(option == "direct")
              triple_conjugate(i,j,k,a,b) +=  8.0 * triple_residuum(i,j,k,a,b);
            else if (option == "bc exchange")
              triple_conjugate(i,j,k,a,b) += -4.0 * triple_residuum(i,j,k,a,c);
            else if (option == "ba exchange")
              triple_conjugate(i,j,k,a,b) += -4.0 * triple_residuum(i,j,k,b,a);
            else if (option == "ca exchange")
              triple_conjugate(i,j,k,a,b) += -4.0 * triple_residuum(i,j,k,c,b);
            else if (option == "bca exchange")
              triple_conjugate(i,j,k,a,b) +=  2.0 * triple_residuum(i,j,k,b,c);
            else if (option == "cab exchange")
              triple_conjugate(i,j,k,a,b) +=  2.0 * triple_residuum(i,j,k,c,a);
            else if (option == "all")
              triple_conjugate(i,j,k,a,b) +=   8.0 * triple_residuum(i,j,k,a,b) - 4.0 * triple_residuum(i,j,k,b,a) 
                                             - 4.0 * triple_residuum(i,j,k,a,c) - 4.0 * triple_residuum(i,j,k,c,b)
                                             + 2.0 * triple_residuum(i,j,k,c,a) + 2.0 * triple_residuum(i,j,k,b,c);

            energy += ( triple_residuum(i,j,k,a,b) + triple_residuum(j,i,k,b,a) +
                        triple_residuum(i,k,j,a,c) + triple_residuum(k,j,i,c,b) +
                        triple_residuum(k,i,j,c,a) + triple_residuum(j,k,i,b,c) ) 
                      * triple_conjugate(i,j,k,a,b) * this->refpoint.getTriplesPropagator(i,j,k,a,b,c);
            
      }
      
      // Print energy info (in eV and Hartree/electron).
      if (option == "all") 
        printf("[INFO]: Perturbative Triples Correction Energy in eV (total (T)) calculated: %+10.10f \n", energy*27.211386);
      else 
        printf("[INFO]: Perturbative Triples Correction Energy in eV (%s term) calculated: %+10.10f \n", option.c_str(), energy*27.211386);
      std::cout << " " << std::endl; 
      printf("[INFO]: Perturbative Triples Correction Energy in Hartree/electron calculated: %+10.10f \n", energy/(2*refpoint.n_occupied));
      std::cout << " " << std::endl;
      triple_residuum.clear(); triple_conjugate.clear();
      
      // Different implementation of perturbative triples.
      //for (auto const [kvector,connections] : mappingHHH)
      //for (auto const &[i,j,k] : connections)
      //for (auto const &[a,b,c]: mappingPPP[kvector])
      //  triple_residuum(i,j,k,a,b) += ( ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) ) * this->mixedSpaceCoulomb(k,c) 
      //                                + ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) ) * this->mixedSpaceCoulomb(j,b)
      //                                + ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) ) * this->mixedSpaceCoulomb(i,a) );
      //energy = 0.0;
      //for (auto const [kvector,connections] : mappingHHH)
      //for (auto const &[i,j,k] : connections)
      //for (auto const &[a,b,c]: mappingPPP[kvector])
      //  energy += 1./3 * ( 4.0 * triple_residuum(i,j,k,a,b) - 6.0 * triple_residuum(i,j,k,a,c) + 2.0 * triple_residuum(i,j,k,b,c) )
      //                 * triple_residuum(i,j,k,a,b) * this->refpoint.getTriplesPropagator(i,j,k,a,b,c);
 
      //printf("[INFO]: Perturbative Triples Correction Energy in eV calculated: %+10.14f \n", energy*27.211386);
      //std::cout << " " << std::endl;
      //printf("[INFO]: Perturbative Triples Correction Energy in Hartree/electron calculated: %+10.14f \n", energy/(2*refpoint.n_occupied));
      //std::cout << " " << std::endl;
      //triple_residuum.clear(); triple_conjugate.clear();

    }
      
   
    // ------------------------------------------------------------------------------------------------------
 
    template<class StartingPoint> void Correlator<StartingPoint>::_cT() {

    auto occ = this->refpoint.transform_to_integer_notation("occupied");
    auto vir = this->refpoint.transform_to_integer_notation("virtuals");

    // Create a hole-particle-particle amplitude to facilitate the contractions.
    auto amplitude_hpp = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});

    // Clear the amplitude.
    amplitude_hpp.clear();

    // Define the denominator from Piecuch's paper.   
    double contraction(0.0), denominator; 

    //  Initialize the HPP amplitude with the corresponding values and calculate the denominator.
    for (const auto [kvector, connections] : this->mappingPHvertex)
      for (auto const &[i,a] : connections)
        for (auto const &[j,b] : this->mappingPHvertex[-kvector]){
          amplitude_hpp(i,a,b) = this->amplitude(i,j,a);
          contraction += ( 2.0 * this->amplitude(i,j,a) - this->amplitude(j,i,a) ) * this->amplitude(i,j,a); }
    
    denominator = 1.0 + contraction/(2.0*refpoint.n_occupied);

    // Create an intermediate hole-particle-particle tensor (Intermediate Particle Tensor).
    auto _IP_tensor = ueg::fundamentals::Tensor3d(refpoint.n_occupied,refpoint.n_virtuals,refpoint.n_virtuals);

    // Create an intermediate hole-hole-particle tensor (Intermediate Hole Tensor).
    auto _IH_tensor = ueg::fundamentals::Tensor3d(refpoint.n_occupied,refpoint.n_occupied,refpoint.n_virtuals);
  
    // Set the elements of 2 tensors to zero.
    _IP_tensor.clear(); _IH_tensor.clear();

    // Calculate the hole-particle-particle intermediate quantity.
    for (auto const [kvector, connections] : mappingPHvertex)
      for (auto const &[k,c] : connections){
        for (size_t b(0); b < refpoint.n_virtuals; ++b){
          _IP_tensor(k,b,c) = this->mixedSpaceCoulomb(k,c);
          for (auto const &[m,d]: mappingPHvertex[-kvector])
            _IP_tensor(k,b,c) += this->mixedSpaceCoulomb(m,d) * ( 2.0 * amplitude(k,m,c) - amplitude(k,m,d) )
                                 - amplitude(m,k,d) * this->refpoint.getCoulombInteractionVirtuals(b,d);
          for (auto const &[m,d]: mappingPHvertex[vir[b]-occ[k]])
            _IP_tensor(k,b,c) -= amplitude(k,m,d) * this->refpoint.getCoulombInteractionVirtuals(c,d);}
        for(size_t m(0); m < refpoint.n_occupied; ++m)
          for (auto const &[n,b]: mappingPHvertex[vir[c]-occ[m]])
            _IP_tensor(k,b,c) += amplitude(m,n,c) * this->refpoint.getCoulombInteractionOccupied(k,m);}

    // Calculate the hole-hole-particle intermediate quantity.
    for (auto const [kvector, connections] : mappingPHvertex)
      for (auto const &[k,c] : connections){
        for (size_t j(0); j < refpoint.n_occupied; ++j){
          _IH_tensor(j,k,c) = this->mixedSpaceCoulomb(k,c);
          for (auto const &[l,e] : mappingPHvertex[-kvector])
            _IH_tensor(j,k,c) += this->mixedSpaceCoulomb(l,e) * ( 2.0 * amplitude(k,l,c) - amplitude(k,l,e) )
                                 - amplitude(l,k,e) * this->refpoint.getCoulombInteractionOccupied(l,j);
          for (auto const &[l,e] : mappingPHvertex[vir[c]-occ[j]])
            _IH_tensor(j,k,c) -= amplitude(l,j,c) * this->refpoint.getCoulombInteractionOccupied(l,k);}
        for(size_t e(0); e < refpoint.n_virtuals; ++e)
          for (auto const &[j,f] : mappingPHvertex[vir[e]-occ[k]])
            _IH_tensor(j,k,c) += amplitude(k,j,e) * this->refpoint.getCoulombInteractionVirtuals(c,e);}


      // Straight-out implementation of the complete renormalized triples.
      double energy = 0.0;
      for (auto const [kvector,connections] : mappingHHH)
      for (auto const &[i,j,k] : connections)
      for (auto const &[a,b,c]: mappingPPP[kvector]){
        triple_conjugate(i,j,k,a,b) += (  amplitude(i,j,a) * _IP_tensor(k,b,c) + amplitude(j,i,b) * _IP_tensor(k,a,c) + amplitude(i,k,a) * _IP_tensor(j,c,b)
                                        + amplitude(k,i,c) * _IP_tensor(j,a,b) + amplitude(j,k,b) * _IP_tensor(i,c,a) + amplitude(k,j,c) * _IP_tensor(i,b,a)
                                        - amplitude_hpp(i,a,b) * _IH_tensor(j,k,c) - amplitude_hpp(j,b,a) * _IH_tensor(i,k,c) - amplitude_hpp(i,a,c) * _IH_tensor(k,j,b)
                                        - amplitude_hpp(k,c,a) * _IH_tensor(i,j,b) - amplitude_hpp(k,c,b) * _IH_tensor(j,i,a) - amplitude_hpp(j,b,c) * _IH_tensor(k,i,a) );

        triple_residuum(i,j,k,a,b) += ( ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) ) * this->mixedSpaceCoulomb(k,c) 
                                      + ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) ) * this->mixedSpaceCoulomb(j,b)
                                      + ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) ) * this->mixedSpaceCoulomb(i,a) )
                                      * this->refpoint.getTriplesPropagator(i,j,k,a,b,c); }

      for (auto const [kvector,connections] : mappingHHH)
      for (auto const &[i,j,k] : connections)
      for (auto const &[a,b,c]: mappingPPP[kvector])
        energy += ( 4.0 * triple_residuum(i,j,k,a,b) - 6.0 * triple_residuum(i,j,k,a,c) + 2.0 * triple_residuum(i,j,k,b,c) )
                  * triple_conjugate(i,j,k,a,b)  / 3.0;

      printf("[INFO]: Complete (T) Correction Energy in eV calculated: %+10.10f \n", energy*27.211386);
      std::cout << " " << std::endl;
      printf("[INFO]: Complete (T) Correction Energy in Hartree/electron calculated: %+10.10f \n", energy/(2*refpoint.n_occupied));
      std::cout << " " << std::endl;
      printf("[INFO]: Complete (T) Denominator is (per electron): %+10.10f \n", denominator);
      std::cout << " " << std::endl;

      triple_residuum.clear(); triple_conjugate.clear();

    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> vector_dtuples_t Correlator<StartingPoint>::_OTF_CCDT(bool structure_factor) {

    auto occ = this->refpoint.transform_to_integer_notation("occupied");
    auto vir = this->refpoint.transform_to_integer_notation("virtuals");

    // Create a hole-particle-particle amplitude to facilitate the contractions.
    auto amplitude_hpp = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});

    // Create a matrix of size (No x Nv) where its elements (pairs of an occupied and a virtual orbital)
    // with the same momentum transfer will correspond to the same integer value.
    auto _momentum_index = ueg::fundamentals::Tensor2d({refpoint.n_virtuals, refpoint.n_occupied}); 
    
    // Clear the amplitude and matrix.
    amplitude_hpp.clear(); _momentum_index.clear();

    // Create a vector to store the Perturbative Triples' Structure Factor results.
    std::vector<double> _structure_factor(mappingPHvertex.size());

    // Create a vector of tuples (each tuple will contain 3 elements) to return the structure factors results.
    vector_dtuples_t structure_factor_results(this->mappingPHvertex.size());
    
    // Create a map that associates a unique momentum transfer vector to a position in the structure factor vector.
    std::map <iVec_t, size_t> _position;
    
    // Initialize the matrix and the map.
    size_t index(0);
    for (auto const [kvector, connections]: mappingPHvertex){
      _position[kvector] = index;
      for (auto const &[i,a] : connections)
        _momentum_index(a,i) = index;
      index++; }

    // Initialize the HPP amplitude with the corresponding values.
    for (auto const [kvector, connections] : this->mappingPHvertex)
      for (auto const &[i,a] : connections)
        for (auto const &[j,b] : this->mappingPHvertex[-kvector])
          amplitude_hpp(i,a,b) = this->amplitude(i,j,a);

    // Create a 2D sliced residuum in the virtual space for the on the fly calculation.
    auto triple_residuum_ijk  = ueg::fundamentals::Tensor2d(refpoint.n_virtuals,refpoint.n_virtuals);

    // Create three 'bare' (will not contain the Coulomb interaction) triples residuums of size Nq x Nq, 
    // where q is a unique momentum transfer vector.
    auto bare_residuum_ab  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());
    auto bare_residuum_ca  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());
    auto bare_residuum_bc  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());
    
    // Initialize Perturbative Triples energy to 0 and all the needed terms (allocated once and for all). B stands for bottom and T for top parts.
    double energy(0.0), B1(0.0), B2(0.0), B3(0.0), B4(0.0), B5(0.0), B6(0.0), factor(0.0);
    double T1_1(0.0), T1_2(0.0), T1_3(0.0), T2_1(0.0), T2_2(0.0), T2_3(0.0), T3_1(0.0), T3_2(0.0), T3_3(0.0);
    double T4_1(0.0), T4_2(0.0), T4_3(0.0), T5_1(0.0), T5_2(0.0), T5_3(0.0), T6_1(0.0), T6_2(0.0), T6_3(0.0);
    double _renorm_factor(2*M_PI/this->refpoint.lattice.constant), constant(4*M_PI/this->refpoint.lattice.boxVolume);
    double R1(0.0), R2(0.0), R3(0.0), Q1(0.0), Q2(0.0), _norm_kvector, _coulomb_interaction;
    
    // First indices comparison and flag 'identical'.
    size_t identical(0);
  
    // Define a 3D-kvector with coordinates integer numbers.
    iVec_t _kvector;

    // Slice the ijk-space in an "upper triangular" form.
    for (size_t i(0); i < refpoint.n_occupied; ++i)
    for (size_t j(i); j < refpoint.n_occupied; ++j)
    for (size_t k(j); k < refpoint.n_occupied; ++k){

        // Evaluate once the equality.
        identical = size_t(i==j) + size_t(j==k);

        // Check if we need to terminate this i,j,k combination.
        if (identical>1) continue;

        // Perform the rescaling of the residuums. 
        triple_residuum_ijk.rescale(0.0);
        bare_residuum_ab.rescale(0.0);
        bare_residuum_ca.rescale(0.0);
        bare_residuum_bc.rescale(0.0);

        // Particle and Hole Contractions for the Perturbative Triples ( R = [T(particle) - T(hole)] * V  ).
        // Search and obtain in the particle-particle map the a,b indices based on
        // the momentum conservation (ki + kj + kk - kc) = (ka + kb).
        // Initialize the three 'bare' residuums by removing the corresponding Coulomb kernel.
        for (size_t c(0); c < refpoint.n_virtuals; ++c){
          _kvector = occ[i] + occ[j] + occ[k] - vir[c];
          for (auto const &[a,b]: mappingPP[_kvector]){
            triple_residuum_ijk(a,b) += ( ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) ) * this->mixedSpaceCoulomb(k,c) 
                                        + ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) ) * this->mixedSpaceCoulomb(j,b) 
                                        + ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) ) * this->mixedSpaceCoulomb(i,a) );
            bare_residuum_ab(_momentum_index(a,i),_momentum_index(b,j)) += ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) );
            bare_residuum_ca(_momentum_index(c,k),_momentum_index(a,i)) += ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) );
            bare_residuum_bc(_momentum_index(b,j),_momentum_index(c,k)) += ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) );   }}
       
        // Energy contraction.
        for (size_t c(0); c < refpoint.n_virtuals; ++c){
          _kvector = occ[i] + occ[j] + occ[k] - vir[c];
          for (auto const &[a,b]: mappingPP[_kvector]){

            // Slice the abc-space in a "lower triangular" way.
            if (a > b || b > c) continue;

            // Calculate the factor without using the if statement.
            factor = ( 1.0 - 0.5 * int(a==b || b==c) ) * this->refpoint.getTriplesPropagator(i,j,k,a,b,c);
     
            // Evaluate the terms (Bottom/Top).
            B1 = triple_residuum_ijk(a,b);
            B2 = triple_residuum_ijk(b,c);
            B3 = triple_residuum_ijk(c,a);
             
            if (structure_factor) {
            R1 = (2.0 * B1 - B2 - B3);
            R2 = (2.0 * B2 - B1 - B3);
            R3 = (2.0 * B3 - B1 - B2);
            T1_1 = bare_residuum_ab(_momentum_index(a,i),_momentum_index(b,j));
            T1_2 = bare_residuum_ca(_momentum_index(c,k),_momentum_index(a,i));
            T1_3 = bare_residuum_bc(_momentum_index(b,j),_momentum_index(c,k));
            T2_1 = bare_residuum_ab(_momentum_index(b,i),_momentum_index(c,j));
            T2_2 = bare_residuum_ca(_momentum_index(a,k),_momentum_index(b,i));
            T2_3 = bare_residuum_bc(_momentum_index(c,j),_momentum_index(a,k));
            T3_1 = bare_residuum_ab(_momentum_index(c,i),_momentum_index(a,j));
            T3_2 = bare_residuum_ca(_momentum_index(b,k),_momentum_index(c,i));
            T3_3 = bare_residuum_bc(_momentum_index(a,j),_momentum_index(b,k)); }

            switch (identical) { 
              case 0:  // all indices are distinct.

                // Only now we access the elements that we need.
                B4 = triple_residuum_ijk(a,c);
                B5 = triple_residuum_ijk(b,a);
                B6 = triple_residuum_ijk(c,b);

                // Evaluate the (T) energy.
                energy += 2.0 * factor * ( 4.0 * ( B1*B1 + B2*B2 + B3*B3 + B4*B4 + B5*B5 + B6*B6 ) + 2.0 * (B1*B2 + B1*B3 + B2*B3 + B4*B5 + B4*B6 + B5*B6)
                                          -4.0 * (B4+B5+B6) * (B1+B2+B3) );

                if (structure_factor) {
                Q1 = ((B1 + B2 + B3) - 2.0 * (B4 + B5 +B6));
                Q2 = ((B4 + B5 + B6) - 2.0 * (B1 + B2 +B3));
                T4_1 = bare_residuum_ab(_momentum_index(a,i),_momentum_index(c,j));
                T4_2 = bare_residuum_ca(_momentum_index(b,k),_momentum_index(a,i));
                T4_3 = bare_residuum_bc(_momentum_index(c,j),_momentum_index(b,k));
                T5_1 = bare_residuum_ab(_momentum_index(b,i),_momentum_index(a,j));
                T5_2 = bare_residuum_ca(_momentum_index(c,k),_momentum_index(b,i));
                T5_3 = bare_residuum_bc(_momentum_index(a,j),_momentum_index(c,k));
                T6_1 = bare_residuum_ab(_momentum_index(c,i),_momentum_index(b,j));
                T6_2 = bare_residuum_ca(_momentum_index(a,k),_momentum_index(c,i));
                T6_3 = bare_residuum_bc(_momentum_index(b,j),_momentum_index(a,k)); 

                // Evaluate the (T) structure factor.

                _structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( 3.0 * ( B1 * T1_3 + B4 * T4_3 ) + Q1 * T1_3 + Q2 * T4_3);
                _structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( 3.0 * ( B1 * T1_2 + B6 * T6_2 ) + Q1 * T1_2 + Q2 * T6_2);
                _structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( 3.0 * ( B1 * T1_1 + B5 * T5_1 ) + Q1 * T1_1 + Q2 * T5_1);
                _structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( 3.0 * ( B2 * T2_3 + B5 * T5_3 ) + Q1 * T2_3 + Q2 * T5_3);
                _structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( 3.0 * ( B2 * T2_2 + B4 * T4_2 ) + Q1 * T2_2 + Q2 * T4_2);
                _structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( 3.0 * ( B2 * T2_1 + B6 * T6_1 ) + Q1 * T2_1 + Q2 * T6_1); 
                _structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( 3.0 * ( B3 * T3_3 + B6 * T6_3 ) + Q1 * T3_3 + Q2 * T6_3);
                _structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( 3.0 * ( B3 * T3_2 + B5 * T5_2 ) + Q1 * T3_2 + Q2 * T5_2);
                _structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( 3.0 * ( B3 * T3_1 + B4 * T4_1 ) + Q1 * T3_1 + Q2 * T4_1); }

                break;

              case 1: // two indices are identical.

                // Evaluate the energy and the structure factor.
                energy += 2.0 * factor * 2.0 * ( B1*B1 + B2*B2 + B3*B3 - B1*B2 - B1*B3 - B3*B2 );
  
                if (structure_factor) {
                _structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( R1 * T1_3);
                _structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( R1 * T1_2);
                _structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( R1 * T1_1);
                _structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( R2 * T2_3);
                _structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( R2 * T2_2); 
                _structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( R2 * T2_1);
                _structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( R3 * T3_3);
                _structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( R3 * T3_2);
                _structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( R3 * T3_1); }

                break;
             
            } 
          } //ab-loop
        } //c-loop
      } //k-loop

      // Print energy info (in eV and Hartree/electron).
      printf("[INFO]: Perturbative triples correction energy in Hartree/electron calculated: %+10.8f \n", energy/(2*refpoint.n_occupied));
      std::cout << " " << std::endl;
      
      // Return the structure factor results.
      if(structure_factor) {
        size_t _index = 0;
        for (auto const &[kvector,connections]: mappingPHvertex){
          _norm_kvector = _renorm_factor * (kvector).norm();
          _coulomb_interaction = constant/std::pow(_norm_kvector, 2.0);
          structure_factor_results[_index]=std::make_tuple(_norm_kvector, _structure_factor[_position[kvector]], _coulomb_interaction, 0.00000000);
          ++_index; }}
      
      return structure_factor_results;
      
    }

    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> vector_dtuples_t Correlator<StartingPoint>::_OTF_CCDcT(bool structure_factor) {

    auto occ = this->refpoint.transform_to_integer_notation("occupied");
    auto vir = this->refpoint.transform_to_integer_notation("virtuals");

    // Create a hole-particle-particle amplitude to facilitate the contractions.
    auto amplitude_hpp = ueg::fundamentals::Tensor3d({refpoint.n_occupied, refpoint.n_virtuals, refpoint.n_virtuals});

    // Create a matrix of size (No x Nv) where its elements (pairs of an occupied and a virtual orbital)
    // with the same momentum transfer will correspond to the same integer value.
    auto _momentum_index = ueg::fundamentals::Tensor2d({refpoint.n_virtuals, refpoint.n_occupied}); 
    
    // Clear the amplitude and matrix.
    amplitude_hpp.clear(); _momentum_index.clear();

    // Create a vector to store the Perturbative and Complete Renormalized Triples' Structure Factor results.
    std::vector<double> _pt_structure_factor(mappingPHvertex.size());
    std::vector<double> _cr_structure_factor(mappingPHvertex.size()); 
    
    // Create a vector of tuples (each tuple will contain 4 elements) to return the structure factors results.
    vector_dtuples_t structure_factor_results(this->mappingPHvertex.size());
    
    // Create a map that associates a unique momentum transfer vector to a position in the structure factor vector.
    std::map <iVec_t, size_t> _position;

    // Initialize the matrix and the map.
    size_t index(0);
    for (auto const [kvector, connections]: mappingPHvertex){
      _position[kvector] = index;
      for (auto const &[i,a] : connections)
        _momentum_index(a,i) = index;
      index++; }

    // Define the denominator from Piecuch's paper.   
    double contraction(0.0), denominator; 

    //  Initialize the HPP amplitude with the corresponding values and calculate the denominator.
    for (const auto [kvector, connections] : this->mappingPHvertex)
      for (auto const &[i,a] : connections)
        for (auto const &[j,b] : this->mappingPHvertex[-kvector]){
          amplitude_hpp(i,a,b) = this->amplitude(i,j,a);
          contraction += ( 2.0 * this->amplitude(i,j,a) - this->amplitude(j,i,a) ) * this->amplitude(i,j,a); }
    
    denominator = 1.0 + contraction/(2.0*refpoint.n_occupied);

    // Create an intermediate hole-particle-particle tensor (Intermediate Particle Tensor).
    auto _IP_tensor = ueg::fundamentals::Tensor3d(refpoint.n_occupied,refpoint.n_virtuals,refpoint.n_virtuals);

    // Create an intermediate hole-hole-particle tensor (Intermediate Hole Tensor).
    auto _IH_tensor = ueg::fundamentals::Tensor3d(refpoint.n_occupied,refpoint.n_occupied,refpoint.n_virtuals);

    // Set the elements of 2 tensors to zero.
    _IP_tensor.clear(); _IH_tensor.clear();

    // Calculate the hole-particle-particle intermediate quantity.
    for (auto const [kvector, connections] : mappingPHvertex)
      for (auto const &[k,c] : connections){
        for (size_t b(0); b < refpoint.n_virtuals; ++b){
          _IP_tensor(k,b,c) = this->mixedSpaceCoulomb(k,c);
          for (auto const &[m,d]: mappingPHvertex[-kvector])
            _IP_tensor(k,b,c) += this->mixedSpaceCoulomb(m,d) * ( 2.0 * amplitude(k,m,c) - amplitude(k,m,d) )
                                 - amplitude(m,k,d) * this->refpoint.getCoulombInteractionVirtuals(b,d);
          for (auto const &[m,d]: mappingPHvertex[vir[b]-occ[k]])
            _IP_tensor(k,b,c) -= amplitude(k,m,d) * this->refpoint.getCoulombInteractionVirtuals(c,d);}
        for(size_t m(0); m < refpoint.n_occupied; ++m)
          for (auto const &[n,b]: mappingPHvertex[vir[c]-occ[m]])
            _IP_tensor(k,b,c) += amplitude(m,n,c) * this->refpoint.getCoulombInteractionOccupied(k,m);}

    // Calculate the hole-hole-particle intermediate quantity.
    for (auto const [kvector, connections] : mappingPHvertex)
      for (auto const &[k,c] : connections){
        for (size_t j(0); j < refpoint.n_occupied; ++j){
          _IH_tensor(j,k,c) = this->mixedSpaceCoulomb(k,c);
          for (auto const &[l,e] : mappingPHvertex[-kvector])
            _IH_tensor(j,k,c) += this->mixedSpaceCoulomb(l,e) * ( 2.0 * amplitude(k,l,c) - amplitude(k,l,e) )
                                 - amplitude(l,k,e) * this->refpoint.getCoulombInteractionOccupied(l,j);
          for (auto const &[l,e] : mappingPHvertex[vir[c]-occ[j]])
            _IH_tensor(j,k,c) -= amplitude(l,j,c) * this->refpoint.getCoulombInteractionOccupied(l,k);}
        for(size_t e(0); e < refpoint.n_virtuals; ++e)
          for (auto const &[j,f] : mappingPHvertex[vir[e]-occ[k]])
            _IH_tensor(j,k,c) += amplitude(k,j,e) * this->refpoint.getCoulombInteractionVirtuals(c,e);}

    // Initialize Perturbative Triples energy to 0 and all the terms (allocated once and for all).
    double tr_energy(0.0), cr_energy(0.0), factor(1.0), T1(0.0), T2(0.0), T3(0.0), T4(0.0), T5(0.0), T6(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0);
    double B1_1(0.0), B1_2(0.0), B1_3(0.0), B2_1(0.0), B2_2(0.0), B2_3(0.0), B3_1(0.0), B3_2(0.0), B3_3(0.0);
    double B4_1(0.0), B4_2(0.0), B4_3(0.0), B5_1(0.0), B5_2(0.0), B5_3(0.0), B6_1(0.0), B6_2(0.0), B6_3(0.0);
    double _renorm_factor(2*M_PI/this->refpoint.lattice.constant), constant(4*M_PI/this->refpoint.lattice.boxVolume);
    double R1(0.0), R2(0.0), R3(0.0), Q1(0.0), Q2(0.0), C1(0.0), C2(0.0), C3(0.0), X1(0.0), X2(0.0),  _norm_kvector, _coulomb_interaction;
    
    // Create a 2D sliced residuum in the virtual space for the on the fly calculation.
    auto triple_residuum_ijk = ueg::fundamentals::Tensor2d(refpoint.n_virtuals,refpoint.n_virtuals);

    // Create three 'bare' (will not contain the Coulomb interaction) triples residuums of size Nq x Nq, 
    // where q is a unique momentum transfer vector.
    auto bare_residuum_ab  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());
    auto bare_residuum_ca  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());
    auto bare_residuum_bc  = ueg::fundamentals::Tensor2d(mappingPHvertex.size(),mappingPHvertex.size());

    // Create a 2D sliced quantity (M in Piecuch paper) in the virtual space for the on the fly calculation.
    auto _M_ijk = ueg::fundamentals::Tensor2d(refpoint.n_virtuals,refpoint.n_virtuals);

    // First indices comparison and flag 'identical'.
    size_t identical(0);

    // Define a 3D-kvector with coordinates integer numbers.
    iVec_t kvector;  

    // Slice the ijk-space in an "upper triangular" form. 
    for (size_t i(0); i < refpoint.n_occupied; ++i)
    for (size_t j(i); j < refpoint.n_occupied; ++j)
    for (size_t k(j); k < refpoint.n_occupied; ++k){

        // Evaluate once the equality.
        identical = size_t(i==j) + size_t(j==k);

        // Check if we need to terminate this i,j,k combination.
        if (identical>1) continue;

        // Perform the rescaling of the residuums and M quantity.
        triple_residuum_ijk.rescale(0.0);
        bare_residuum_ab.rescale(0.0);
        bare_residuum_ca.rescale(0.0);
        bare_residuum_bc.rescale(0.0);
        _M_ijk.rescale(0.0);
      
        // Particle and Hole Contractions for the Renormalized Perturbative Triples
        // ( R = [T(particle) - T(hole)] * V and M = [T(hole)*I(particle) - T(particle)*I(hole)] ).
        // Search and obtain in the particle-particle map the a,b indices based
        // on the momentum conservation (ki + kj + kk - kc) = (ka + kb).
        // Initialize also the three 'bare' residuums by removing the corresponding Coulomb kernel.
        for (size_t c(0); c < refpoint.n_virtuals; ++c){
          kvector = occ[i] + occ[j] + occ[k] - vir[c];
          for (auto const &[a,b]: mappingPP[kvector]){
            triple_residuum_ijk(a,b) += ( ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) ) * this->mixedSpaceCoulomb(k,c) 
                                        + ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) ) * this->mixedSpaceCoulomb(j,b)
                                        + ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) ) * this->mixedSpaceCoulomb(i,a) );
            bare_residuum_ab(_momentum_index(a,i),_momentum_index(b,j)) += ( amplitude(i,j,a) + amplitude(j,i,b) - amplitude_hpp(i,a,b) - amplitude_hpp(j,b,a) );
            bare_residuum_ca(_momentum_index(c,k),_momentum_index(a,i)) += ( amplitude(k,i,c) + amplitude(i,k,a) - amplitude_hpp(k,c,a) - amplitude_hpp(i,a,c) );
            bare_residuum_bc(_momentum_index(b,j),_momentum_index(c,k)) += ( amplitude(j,k,b) + amplitude(k,j,c) - amplitude_hpp(k,c,b) - amplitude_hpp(j,b,c) );
            _M_ijk(a,b) += (  amplitude(i,j,a) * _IP_tensor(k,b,c) + amplitude(j,i,b) * _IP_tensor(k,a,c) + amplitude(i,k,a) * _IP_tensor(j,c,b)
                            + amplitude(k,i,c) * _IP_tensor(j,a,b) + amplitude(j,k,b) * _IP_tensor(i,c,a) + amplitude(k,j,c) * _IP_tensor(i,b,a)
                            - amplitude_hpp(i,a,b) * _IH_tensor(j,k,c) - amplitude_hpp(j,b,a) * _IH_tensor(i,k,c) - amplitude_hpp(i,a,c) * _IH_tensor(k,j,b)
                            - amplitude_hpp(k,c,a) * _IH_tensor(i,j,b) - amplitude_hpp(k,c,b) * _IH_tensor(j,i,a) - amplitude_hpp(j,b,c) * _IH_tensor(k,i,a) ); }}

        // Energy contraction.
        for (size_t c(0); c < refpoint.n_virtuals; ++c){
          kvector = occ[i] + occ[j] + occ[k] - vir[c];
          for (auto const &[a,b]: mappingPP[kvector]){

            // Slice the abc-space in an "lower triangular" way.
            if (a > b || b > c) continue;

            // Calculate the factor without using the if statement;
            factor = ( 1.0 - 0.5 * int(a==b || b==c) ) * this->refpoint.getTriplesPropagator(i,j,k,a,b,c);
     
            // Evaluate the terms (amplitude T and M quantity terms) for the contraction.
            T1 = triple_residuum_ijk(a,b);
            T2 = triple_residuum_ijk(b,c);
            T3 = triple_residuum_ijk(c,a);
            M1 = _M_ijk(a,b);
            M2 = _M_ijk(b,c);
            M3 = _M_ijk(c,a);

            //if(structure_factor){
            B1_1 = bare_residuum_ab(_momentum_index(a,i),_momentum_index(b,j));
            B1_2 = bare_residuum_ca(_momentum_index(c,k),_momentum_index(a,i));
            B1_3 = bare_residuum_bc(_momentum_index(b,j),_momentum_index(c,k));
            B2_1 = bare_residuum_ab(_momentum_index(b,i),_momentum_index(c,j));
            B2_2 = bare_residuum_ca(_momentum_index(a,k),_momentum_index(b,i));
            B2_3 = bare_residuum_bc(_momentum_index(c,j),_momentum_index(a,k));
            B3_1 = bare_residuum_ab(_momentum_index(c,i),_momentum_index(a,j));
            B3_2 = bare_residuum_ca(_momentum_index(b,k),_momentum_index(c,i));
            B3_3 = bare_residuum_bc(_momentum_index(a,j),_momentum_index(b,k));
            R1 = (2.0 * M1 - M2 - M3);
            R2 = (2.0 * M2 - M1 - M3);
            R3 = (2.0 * M3 - M1 - M2);
            C1 = (2.0 * T1 - T2 - T3);
            C2 = (2.0 * T2 - T1 - T3);
            C3 = (2.0 * T3 - T1 - T2); 

            switch (identical) { 

              case 0:  // all indices are distinct.

                // Only now we access the elements that we need (of the amplitude and M quantity).
                T4 = triple_residuum_ijk(a,c);
                T5 = triple_residuum_ijk(b,a);
                T6 = triple_residuum_ijk(c,b);
                M4 = _M_ijk(a,c);
                M5 = _M_ijk(b,a);
                M6 = _M_ijk(c,b); 

                // Evaluate the triples and complete renormalized energies.
                tr_energy += 2.0 * factor * ( 4.0 * ( T1*T1 + T2*T2 + T3*T3 + T4*T4 + T5*T5 + T6*T6 ) + 2.0 * (T1*T2 + T1*T3 + T2*T3 + T4*T5 + T4*T6 + T5*T6)
                                             -4.0 * (T4+T5+T6) * (T1+T2+T3) );
                cr_energy += 2.0 * factor * ( 3.0 * ( T1*M1 + T2*M2 + T3*M3 + T4*M4 + T5*M5 + T6*M6 ) + ((M1+M2+M3) - 2.0 * (M4+M5+M6)) * (T1+T2+T3)
                                            + ((M4+M5+M6) - 2.0 * (M1+M2+M3)) * (T4+T5+T6) );

                if (structure_factor) { 

                B4_1 = bare_residuum_ab(_momentum_index(a,i),_momentum_index(c,j));
                B4_2 = bare_residuum_ca(_momentum_index(b,k),_momentum_index(a,i));
                B4_3 = bare_residuum_bc(_momentum_index(c,j),_momentum_index(b,k));
                B5_1 = bare_residuum_ab(_momentum_index(b,i),_momentum_index(a,j));
                B5_2 = bare_residuum_ca(_momentum_index(c,k),_momentum_index(b,i));
                B5_3 = bare_residuum_bc(_momentum_index(a,j),_momentum_index(c,k));
                B6_1 = bare_residuum_ab(_momentum_index(c,i),_momentum_index(b,j));
                B6_2 = bare_residuum_ca(_momentum_index(a,k),_momentum_index(c,i));
                B6_3 = bare_residuum_bc(_momentum_index(b,j),_momentum_index(a,k));
                Q1 = ((M1 + M2 + M3) - 2.0 * (M4 + M5 + M6));
                Q2 = ((M4 + M5 + M6) - 2.0 * (M1 + M2 + M3));
                X1 = ((T1 + T2 + T3) - 2.0 * (T4 + T5 + T6));
                X2 = ((T4 + T5 + T6) - 2.0 * (T1 + T2 + T3));
            
                // Evaluate the triples and complete renormalized structure factors.
                _pt_structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( 3.0 * ( T1 * B1_3 + T4 * B4_3 ) + X1 * B1_3 + X2 * B4_3);
                _pt_structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( 3.0 * ( T1 * B1_2 + T6 * B6_2 ) + X1 * B1_2 + X2 * B6_2);
                _pt_structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( 3.0 * ( T1 * B1_1 + T5 * B5_1 ) + X1 * B1_1 + X2 * B5_1);
                _pt_structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( 3.0 * ( T2 * B2_3 + T5 * B5_3 ) + X1 * B2_3 + X2 * B5_3);
                _pt_structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( 3.0 * ( T2 * B2_2 + T4 * B4_2 ) + X1 * B2_2 + X2 * B4_2);
                _pt_structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( 3.0 * ( T2 * B2_1 + T6 * B6_1 ) + X1 * B2_1 + X2 * B6_1); 
                _pt_structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( 3.0 * ( T3 * B3_3 + T6 * B6_3 ) + X1 * B3_3 + X2 * B6_3);
                _pt_structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( 3.0 * ( T3 * B3_2 + T5 * B5_2 ) + X1 * B3_2 + X2 * B5_2);
                _pt_structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( 3.0 * ( T3 * B3_1 + T4 * B4_1 ) + X1 * B3_1 + X2 * B4_1);

                _cr_structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( 3.0 * ( M1 * B1_3 + M4 * B4_3 ) + Q1 * B1_3 + Q2 * B4_3);
                _cr_structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( 3.0 * ( M1 * B1_2 + M6 * B6_2 ) + Q1 * B1_2 + Q2 * B6_2);
                _cr_structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( 3.0 * ( M1 * B1_1 + M5 * B5_1 ) + Q1 * B1_1 + Q2 * B5_1);
                _cr_structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( 3.0 * ( M2 * B2_3 + M5 * B5_3 ) + Q1 * B2_3 + Q2 * B5_3);
                _cr_structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( 3.0 * ( M2 * B2_2 + M4 * B4_2 ) + Q1 * B2_2 + Q2 * B4_2);
                _cr_structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( 3.0 * ( M2 * B2_1 + M6 * B6_1 ) + Q1 * B2_1 + Q2 * B6_1); 
                _cr_structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( 3.0 * ( M3 * B3_3 + M6 * B6_3 ) + Q1 * B3_3 + Q2 * B6_3);
                _cr_structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( 3.0 * ( M3 * B3_2 + M5 * B5_2 ) + Q1 * B3_2 + Q2 * B5_2);
                _cr_structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( 3.0 * ( M3 * B3_1 + M4 * B4_1 ) + Q1 * B3_1 + Q2 * B4_1); } 

                break;

              case 1: // two indices are identical.

                // Evaluate the triples and complete renormalized energies.
                tr_energy += 2.0 * factor * 2.0 * ( T1*T1 + T2*T2 + T3*T3 - T1*T2 - T1*T3 - T3*T2 );
                cr_energy += 2.0 * factor * ( 3.0 * ( T1*M1 + T2*M2 + T3*M3 ) - (T1+T2+T3) * (M1+M2+M3) );

                // Evaluate the triples and complete renormalized structure factors.
                if (structure_factor){
                _pt_structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( C1 * B1_3 );
                _pt_structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( C1 * B1_2 );
                _pt_structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( C1 * B1_1 );
                _pt_structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( C2 * B2_3 );
                _pt_structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( C2 * B2_2 ); 
                _pt_structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( C2 * B2_1 );
                _pt_structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( C3 * B3_3 );
                _pt_structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( C3 * B3_2 );
                _pt_structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( C3 * B3_1 );

                _cr_structure_factor[_momentum_index(a,i)] += 2.0 * factor * ( R1 * B1_3 );
                _cr_structure_factor[_momentum_index(b,j)] += 2.0 * factor * ( R1 * B1_2 );
                _cr_structure_factor[_momentum_index(c,k)] += 2.0 * factor * ( R1 * B1_1 );
                _cr_structure_factor[_momentum_index(b,i)] += 2.0 * factor * ( R2 * B2_3 );
                _cr_structure_factor[_momentum_index(c,j)] += 2.0 * factor * ( R2 * B2_2 ); 
                _cr_structure_factor[_momentum_index(a,k)] += 2.0 * factor * ( R2 * B2_1 );
                _cr_structure_factor[_momentum_index(c,i)] += 2.0 * factor * ( R3 * B3_3 );
                _cr_structure_factor[_momentum_index(a,j)] += 2.0 * factor * ( R3 * B3_2 );
                _cr_structure_factor[_momentum_index(b,k)] += 2.0 * factor * ( R3 * B3_1 ); }

                break;
            }
          } //ab-loop
        } //c-loop
      } //k-loop

      // Print energy info (in eV and Hartree/electron).
      printf("[INFO]: Perturbative triples correction energy in Hartree/electron calculated: %+10.8f \n", tr_energy/(2*refpoint.n_occupied));
      std::cout << " " << std::endl;
      printf("[INFO]: Complete (T) correction energy in Hartree/electron calculated: %+10.8f \n", cr_energy/(2*refpoint.n_occupied));
      std::cout << " " << std::endl;
      printf("[INFO]: Complete (T) denominator (per electron): %+10.8f \n", denominator);
      std::cout << " " << std::endl;

      // Return the structure factor results.
      if(structure_factor) {
        size_t _index = 0;
        for (auto const &[kvector,connections]: mappingPHvertex){
          _norm_kvector = _renorm_factor * (kvector).norm();
          _coulomb_interaction = constant/std::pow(_norm_kvector, 2.0);
          structure_factor_results[_index]=std::make_tuple(_norm_kvector, _pt_structure_factor[_position[kvector]], _cr_structure_factor[_position[kvector]],_coulomb_interaction);
          ++_index; }}
      
      return structure_factor_results;
        
    }

    // ------------------------------------------------------------------------------------------------------
    //                                   IMPLEMENTATION OF ONE-SHOT METHODS                                 
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> ueg::concepts::Amplitude Correlator<StartingPoint>::one_shot_methods(std::string method) {
      
      // Distinguish between the one-shot methods.
      if (method == "MP3"){
        this->add_to_residuum(ueg::utility::diagrams::smp3::identifiers); this->symmetrize();
        this->add_to_residuum(ueg::utility::diagrams::nsmp3::identifiers); }
      else if (method == "MP4"){
        this->add_to_residuum(ueg::utility::diagrams::smp4::identifiers); this->symmetrize();
        this->add_to_residuum(ueg::utility::diagrams::nsmp4::identifiers); }
      else { throw std::invalid_argument("Provided method wasn't found!"+method); exit(1); }

      // Calculate the amplitude based on the re-calculated residuum of the provided method.
      this->amplitude.set_from_product(this->residuum.values, this->propagator.values); 

      // Return the amplitude.
      return this->amplitude; 
    }


    // ------------------------------------------------------------------------------------------------------
    //                                 IMPLEMENTATION OF SELF CONSISTENT CYCLE
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> ueg::concepts::Residuum Correlator<StartingPoint>::self_consistent_cycle_body(std::string method) {
 
      // Calculate Residuum according to given method.
      if (method == "CCD") {
        this->residuum.set_from_product(0.5, this->interaction.values);
        this->add_to_residuum(ueg::utility::diagrams::symmetrizable::identifiers); this->symmetrize(); 
        this->add_to_residuum(ueg::utility::diagrams::unsymmetrizable::identifiers); }
      else if (method == "drCCD") {
        this->residuum.set_from_product(0.5, this->interaction.values);
        this->add_to_residuum(ueg::utility::diagrams::drccd::identifiers); this->symmetrize(); } 
      else if (method ==  "ldrCCD") {
        this->residuum.set_from_product(0.5,interaction.values);
        this->add_to_residuum(ueg::utility::diagrams::ldrccd::identifiers); this->symmetrize(); }
      else { throw std::invalid_argument("Provided method wasn't found! " +method); exit(1); }
 
      return this->residuum;

    }

    // ------------------------------------------------------------------------------------------------------
    //                                         IMPLEMENTATION OF DIIS                                       
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> tensor_tuple_t Correlator<StartingPoint>::calculate_errors_and_amplitudes() {
  
      // Copy the current amplitude as the one that received from previous iteration.
      auto preceding_amplitude = this->amplitude;

      // The error that occur for the first iteration is the amplitude itself.
      ueg::concepts::Amplitude computed_error = this->amplitude;
      
      // Initialize the empty list of stored amplitudes and errors with the foregoing 'values'.
      if(stored_amplitudes.empty()) {
        stored_amplitudes.push_back(this->amplitude);
        stored_errors.push_back(this->amplitude);
      }
      // If we have the initialized 'values' calculate the new amplitude and error.
      else {
        this->amplitude.set_from_product(this->residuum.values, this->propagator.values);
        computed_error.set_from_subtraction(this->amplitude.values, preceding_amplitude.values);
      }

      return std::make_tuple(this->amplitude, computed_error);
    }

    template<class StartingPoint> ueg::concepts::Amplitude Correlator<StartingPoint>::diis_extrapolation(tensor_tuple_t amplitude_and_error) {
      
      // Set the number of amplitudes to be stored.
      constexpr size_t amplitudes_to_store = 4;

      ueg::concepts::Amplitude computed_error;
      
      // Get the received amplitude and error and store them to the corresponding deques.
      std::tie(this->amplitude,computed_error) = amplitude_and_error;
      stored_amplitudes.push_back(this->amplitude);
      stored_errors.push_back(computed_error);

      // If the size is not equal to the given number of stored amplitudes, diis extrapolation cannot start. Proceed the self consistency with current amplitude.
      if (!stored_amplitudes.empty() && stored_amplitudes.size() < amplitudes_to_store) {
        return this->amplitude;
      }
      // If the size is greater, delete the first elements at the stored amplitudes and errors.
      else if (stored_amplitudes.size() > amplitudes_to_store) {
        stored_amplitudes.pop_front();  stored_errors.pop_front();
      }      
      
      // Declare the matrices of the DIIS extrapolation.
      arma::vec coefficients(stored_amplitudes.size()+1);
      arma::mat cc_diis_matrix(stored_amplitudes.size()+1, stored_amplitudes.size()+1);

      // A*X = B, with B = [0 0 ... 0 -1], X to be determined and A the CC DIIS matrix. Set B matrix.
      arma::vec right_hand_matrix(stored_amplitudes.size()+1, arma::fill::zeros);
      right_hand_matrix.back() = -1;

      // Set the moduli of the coefficients and lagrange multiplier.
      arma::vec lagrange_multiplier_coeff(stored_amplitudes.size()+1);
      for (size_t i(0); i <= stored_amplitudes.size(); ++i) {
        lagrange_multiplier_coeff[i] = -1;
      }
      lagrange_multiplier_coeff.back() = 0;

      // Set the subblock of the CC DIIS matrix as the inner product of the error 'vectors'. Symmetric for real values.
      double init_value = (0.0);
      for (size_t row(0); row < stored_errors.size(); ++row) {
        for (size_t column(0); column < stored_errors.size(); ++column) {
          if(column <= row) {
            auto element = std::inner_product(std::begin(stored_errors.at(row).values), std::end(stored_errors.at(row).values), std::begin(stored_errors.at(column).values), init_value);
            cc_diis_matrix(row,column) = element;
            cc_diis_matrix(column,row) = element;
          }
        }
      }

      // Set the last column and row of the CC DIIS matrix with the moduli of the coefficients and lagrange multiplier.
      cc_diis_matrix.col(stored_errors.size()) = lagrange_multiplier_coeff;
      cc_diis_matrix.row(stored_errors.size()) = lagrange_multiplier_coeff.t();

      // Solve the system of linear equations.
      coefficients = arma::solve(cc_diis_matrix, right_hand_matrix);
      
      // Clear the amplitude for the DIIS extrapolation.
      this->amplitude.clear();
 
      // Declare scaled amplitudes.
      ueg::concepts::Amplitude scaled_amplitude = this->amplitude;
   
      // Return the DIIS amplitude as a linear combination of the computed scaled amplitudes.
      for (size_t i(0); i < stored_errors.size(); ++i) {
        scaled_amplitude.set_from_product(coefficients[i], stored_amplitudes[i]);
        this->amplitude.set_from_addition(this->amplitude.values, scaled_amplitude.values);
      }

      return this->amplitude;
    }


    // ------------------------------------------------------------------------------------------------------
    //                                IMPLEMENTATION OF CHANNELS DECOMPOSITION                              
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> ueg::concepts::Amplitude Correlator<StartingPoint>::decomposition(std::string diagram) {
      // Use of lambda function so that we don't need to declare new types
      auto find_diagram = [] (std::vector<std::string> diagram_class, std::string diagram_choice) {
        return ( std::find(diagram_class.begin(), diagram_class.end(), diagram_choice) != diagram_class.end() );
      };

      if (diagram == "MP2") {
              return this->amplitude; }
      else {

        // Create copy of amplitudes
        auto temporal_amplitude = this->amplitude;
      
        { // Creating a 'scope' ensures that all unnecessary objects will be 'deleted' when it ends
 
          auto temporal_residuum = this->residuum; residuum.clear();

          // Check the class of the requested diagram and proceed accordingly
          if (diagram == "CCD") {
            this->residuum.set_from_product(0.5,interaction.values);
            this->add_to_residuum(ueg::utility::diagrams::symmetrizable::identifiers); this->symmetrize();
            this->add_to_residuum(ueg::utility::diagrams::unsymmetrizable::identifiers); }
          else if (diagram ==  "drCCD") {
            this->residuum.set_from_product(0.5,interaction.values);
            this->add_to_residuum(ueg::utility::diagrams::drccd::identifiers); this->symmetrize();}
          else if (diagram ==  "ldrCCD") {
            this->residuum.set_from_product(0.5,interaction.values);
            this->add_to_residuum(ueg::utility::diagrams::ldrccd::identifiers); this->symmetrize(); }
          else if ( find_diagram(ueg::utility::diagrams::symmetrizable::identifiers, diagram) ) {
            this->add_to_residuum(diagram); this->symmetrize(); }
          else if ( find_diagram(ueg::utility::diagrams::unsymmetrizable::identifiers, diagram) ) {
            this->add_to_residuum(diagram); }
          else { throw std::invalid_argument("Provided diagram wasn't found!"+diagram); exit(1); }

          // Update the amplitude based on the re-calculated residuum of the provided diagram
          this->amplitude.set_from_product(this->residuum.values, this->propagator.values);
        
          // Reset the value of the residuum to the initial one
          this->residuum = temporal_residuum;
          
        }

        // Swap amplitudes
        auto decomposed_amplitude = this->amplitude;
        this->amplitude = temporal_amplitude;
 
        return decomposed_amplitude;
      }

    }

    // ------------------------------------------------------------------------------------------------------
    //                                   IMPLEMENTATION OF STRUCTURE FACTOR                                 
    // ------------------------------------------------------------------------------------------------------

    template<class StartingPoint> vector_dtuples_t Correlator<StartingPoint>::structure_factor(ueg::concepts::Amplitude amplitude_received) {

      size_t _index = 0;
      double _norm_kvector, _direct_structure_factor, _exchange_structure_factor, _coulomb_interaction;
      double _renorm_factor(2*M_PI/this->refpoint.lattice.constant);
      vector_dtuples_t structure_factor_data(this->mappingPHvertex.size());

      for (auto const [kvector, connections] : mappingPHvertex) {
        _norm_kvector = _renorm_factor * kvector.norm();
        // Set-reset temporal values
        _direct_structure_factor = 0.0;  _exchange_structure_factor = 0.0;
        for (auto const &[i,a] : connections) {
          _coulomb_interaction = mixedSpaceCoulomb(i,a);
          for (auto const &[j,b] : mappingPHvertex[-kvector]) {
            _direct_structure_factor += 2.0 * amplitude_received(i,j,a);
            _exchange_structure_factor -= amplitude_received(j,i,a); }}
   
        // Store results
        structure_factor_data[_index]=std::make_tuple(_norm_kvector, _direct_structure_factor, _exchange_structure_factor,_coulomb_interaction);
        ++_index; }

      return structure_factor_data;
      
    }
     
  } // end namespace correlators
}  // end namespace ueg
