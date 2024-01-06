#include "diagrammatics.hpp"

namespace ueg {
namespace diagrammatics {

// -------------------------------------------------------------------------------
// Implentation of HartreeFock methods (usually we put them after the declaration)
// -------------------------------------------------------------------------------
//

//double HartreeFock::getPropagator(const size_t &i, const size_t &j, const size_t &a, const size_t &b) const {
double HartreeFock::getPropagator(const size_t &i, const size_t &j, const size_t &a, const size_t &b) const {
  return 1.0/(en_occ[i] + en_occ[j] - en_vir[a] - en_vir[b]);
};

double HartreeFock::getTriplesPropagator(const size_t &i, const size_t &j, const size_t &k, const size_t &a, const size_t &c, const size_t &d) const {
  return 1.0/(en_occ[i] + en_occ[j] + en_occ[k] - en_vir[a] - en_vir[c] - en_vir[d]);
};

void HartreeFock::getHFenergy() const {

  double kin_energy(0.0), exchange_energy(0.0), total_energy, hartree_energy, spins(2);

  for (size_t i(0); i < n_occupied; ++i){
    kin_energy += 0.5 * occupied[i].squaredLength();
    for(size_t j(0); j < n_occupied; ++j)
      i!=j ? exchange_energy+=getVijji(occupied[i]-occupied[j]) : exchange_energy+=getViiii(); }
   
   // Fock: 1 closed diagram, 1 interaction line, 2 holes, 1 loop, left-right symmetric.  
   total_energy = (spins * kin_energy + 0.5 * std::pow(-1,1+1+2) * std::pow(-spins,1) * exchange_energy)/(2.0 * n_occupied);
   hartree_energy = (spins * kin_energy)/(2.0 * n_occupied);
   std::cout << std::string(153,'.') << std::endl;
   std::cout << "[INFO]: Hartree Energy calculated in Hartree/electron: " << hartree_energy << std::endl;
   std::cout << std::string(153,'.') << std::endl;
   std::cout << "[INFO]: Hartree Fock Energy calculated in Hartree/electron: " << total_energy << std::endl;
   std::cout << std::string(153,'.') << std::endl;
   std::cout << "[INFO]: HOMO-LUMO-BandGap in Hartree: " << "[" << en_occ.back() << "," << en_vir.front() << "," << en_vir.front()-en_occ.back() << "]" << std::endl;

};

double HartreeFock::getKineticOccupied(const size_t &index) const {
  auto orbital = occupied.at(index);
  return 0.5*orbital.squaredLength();
};

double HartreeFock::getKineticVirtual(const size_t &index) const {
  auto orbital = virtuals.at(index);
  return 0.5*orbital.squaredLength();
};

double HartreeFock::getVijji(const dVec_t &vec) const {
  double Vijji;
  Vijji = 4.0*M_PI/vec.squaredLength()/this->lattice.boxVolume;
  return Vijji;
};

double HartreeFock::getViiii() const {
  double Viiii;
  Viiii = -this->madelung;
  return Viiii;
};

double HartreeFock::getEnergyOccupied(const size_t &index) const {
  double exch_energy = 0.0;
  auto orbital = occupied.at(index);

  for (size_t indexPrime(0); indexPrime < n_occupied; ++indexPrime)
    index!=indexPrime ? exch_energy+=getVijji(orbital-occupied[indexPrime]) : exch_energy+=getViiii(); 

  return 0.5 * orbital.squaredLength() - exch_energy;
};

double HartreeFock::getEnergyVirtual(const size_t &index) const {
  double exch_energy = 0.0;
  auto orbital = virtuals.at(index);

  for (size_t indexPrime(0); indexPrime < n_occupied; ++indexPrime)
    exch_energy += getVijji(orbital-occupied[indexPrime]);

  return 0.5 * orbital.squaredLength() - exch_energy;
};

tensor2d_t HartreeFock::getCoulombInteractionVirtuals() const {
  auto result = tensor2d_t({n_virtuals, n_virtuals});
  for (size_t k(0); k < n_virtuals; ++k)
    for (size_t l(0); l < n_virtuals; ++l){
      if ( k!=l )
        result(k,l) = (4*M_PI/(virtuals[k]-virtuals[l]).squaredLength())/this->lattice.boxVolume;
      else
        result(k,l) = -this->madelung; }
  
  return result;
};

tensor2d_t HartreeFock::getCoulombInteractionOccupied() const {
  auto result = tensor2d_t({n_occupied, n_occupied});
  for (size_t q(0); q < n_occupied; ++q)
    for (size_t r(0); r < n_occupied; ++r){
      if ( q!=r )
        result(q,r) = 4*M_PI/((occupied[q]-occupied[r]).squaredLength())/this->lattice.boxVolume;
      else
        result(q,r) = -this->madelung; }

  return result;
};

tensor2d_t HartreeFock::getCoulombInteractionMixed() const {
  auto result = tensor2d_t({n_occupied, n_virtuals});
  for (size_t k(0); k < n_occupied; ++k)
    for (size_t l(0); l < n_virtuals; ++l)
      result(k,l) = (4*M_PI/(occupied[k]-virtuals[l]).squaredLength())/this->lattice.boxVolume;
  return result;
};

double HartreeFock::getCoulombInteractionVirtuals(const size_t &k, const size_t &l) const {
  double result;
  if( k!=l )
    result = (4*M_PI/(virtuals[k]-virtuals[l]).squaredLength())/this->lattice.boxVolume;
  else
    result = -this->madelung;
  return result;
};

double HartreeFock::getCoulombInteractionOccupied(const size_t &k, const size_t &l) const {
  double result;
  if( k!=l )
    result = (4*M_PI/(occupied[k]-occupied[l]).squaredLength())/this->lattice.boxVolume;
  else
    result = -this->madelung;
  return result;
};
      
tensor2d_t HartreeFock::getCoulombInteraction(std::string subspace) const {
  if (subspace == "virtuals")
    return getCoulombInteractionVirtuals();
  else if (subspace == "occupied")
    return getCoulombInteractionOccupied();
  else if (subspace == "mixed")
    return getCoulombInteractionMixed();
};

};
};
