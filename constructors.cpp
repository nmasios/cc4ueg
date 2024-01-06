#include "constructors.hpp"

namespace ueg {
namespace constructors {

std::vector<size_t> get_available_electron_filling(std::vector<dVec_t> states) {
  std::vector<double> squaredLength, temporal;
  std::vector<size_t> available_filling;
  
  std::transform(states.cbegin(), states.cend(), std::back_inserter(squaredLength),
          [](const dVec_t &state) { return state.squaredLength(); });

  // Make a copy of the vector that we have just created
  std::copy(squaredLength.begin(), squaredLength.end(), std::back_inserter(temporal));

  // Find the unique value of 'square-lengths'
  std::vector<double>::iterator it = std::unique(squaredLength.begin(), squaredLength.end(),
          [](const double &x, const double &y) { return std::abs(x-y) < 1e-12; } );
  squaredLength.resize(std::distance(squaredLength.begin(),it));

  // At this point 'squaredLength' contains only the unique values, so we just need to
  // find in which index these values first appear at 'temporal'
  std::transform(squaredLength.begin()+1,squaredLength.end(),
          std::back_inserter(available_filling), [&temporal](const double &v) {
          return 2*std::distance(temporal.begin(), std::find(temporal.begin(),temporal.end(),v));
          });

  return available_filling;
};

lattice_properties_t lattice_definition(std::string lattice_type) {

  dVec_t a1,a2,a3,b1,b2,b3;

  if (lattice_type == "sc") {
    a1 = dVec_t(+1.0, +0.0, +0.0);
    a2 = dVec_t(+0.0, +1.0, +0.0);
    a3 = dVec_t(+0.0, +0.0, +1.0);
  } else if (lattice_type == "bcc") {
    a1 = dVec_t(-0.5, +0.5, +0.5);
    a2 = dVec_t(+0.5, -0.5, +0.5);
    a3 = dVec_t(+0.5, +0.5, -0.5);
  } else if (lattice_type == "fcc") {
    a1 = dVec_t(+0.0, +0.5, +0.5);
    a2 = dVec_t(+0.5, +0.0, +0.5);
    a3 = dVec_t(+0.5, +0.5, +0.0);
  } else {
    throw std::invalid_argument("Unknown lattice-type provided: "+lattice_type);
  }

  // Calculate cell-volume
  double volume = std::abs(a1.dot(a2.cross(a3)));

  // calculate reciprocal vectors
  b1 = 2*M_PI * a2.cross(a3) / volume;
  b2 = 2*M_PI * a3.cross(a1) / volume;
  b3 = 2*M_PI * a1.cross(a2) / volume;

  return std::make_tuple(a1,a2,a3,b1,b2,b3,volume);
};

// ------------------------------------------------------------------------------------------------
//                               Lattice structure definitions
// ------------------------------------------------------------------------------------------------

Lattice::Lattice(const Lattice &other) :
  a1(other.a1), a2(other.a2), a3(other.a3), b1(other.b1), b2(other.b2), b3(other.b3),
  grid(other.grid), type(other.type), volume(other.volume), constant(other.constant),
  boxVolume(other.boxVolume), radius(other.radius) {};

// ------------------------------------------------------------------------------------------------

Lattice::Lattice(std::string celltype, int gridAxisSize) : type(celltype) {
  std::tie(a1,a2,a3,b1,b2,b3,volume) = lattice_definition(type);
  grid = std::make_pair(-gridAxisSize,+gridAxisSize);

  // Find the maximum radius of a sphere that fits inside the lattice
  double temporal_radius;
  for (int m2 = grid.first; m2 <= grid.second; ++m2) {
    for (int m3 = grid.first; m3 <= grid.second; ++m3) {
      temporal_radius = (grid.first*b1 + m2*b2 + m3*b3).squaredLength();
      if (temporal_radius < radius) radius = temporal_radius;
    }
  }
};

// ------------------------------------------------------------------------------------------------

std::vector<dVec_t> Lattice::getAllStates(double normalization, dVec_t shift) const {

  std::vector<dVec_t> states;
  dVec_t shiftedVec = shift.x * b1 + shift.y * b2 + shift.z * b3; 

  for (int m1 = grid.first; m1 <= grid.second; ++m1)
  for (int m2 = grid.first; m2 <= grid.second; ++m2)
  for (int m3 = grid.first; m3 <= grid.second; ++m3){
    dVec_t gridpoint = m1*b1 + m2*b2 + m3*b3;

    // Check if vector is within the sphere
    if (gridpoint.squaredLength() < radius + 1e-12)
      states.push_back((gridpoint-shiftedVec)/normalization);
    
  }

  std::stable_sort(states.begin(), states.end());
  return states;
};
  
// ------------------------------------------------------------------------------------------------

void Lattice::reset_constant(double value) {
  constant = value; boxVolume = std::pow(constant,3)*volume;
};

void Lattice::reset_constant(size_t n_electrons, double wigner_seitz) {
  double value = std::pow(4*M_PI*n_electrons/(3.0*this->volume), 1.0/3.0)*wigner_seitz;
  this->reset_constant(value);
};

// ------------------------------------------------------------------------------------------------
//                               System structure definitions
// ------------------------------------------------------------------------------------------------

System::System() {
  auto celltype = ueg::utility::ask_user_input<std::string>("Specify the lattice type");
  int gridAxisSize = ueg::utility::ask_user_input<int>("Specify the number of grid-points along an (positive) axis");

  // Construct the lattice
  lattice = Lattice(celltype,gridAxisSize);
  auto states = lattice.getAllStates(lattice.constant);
  
  // Obtain the available number of electons that can fill the system
  auto closedShellElectrons = get_available_electron_filling(states);

  // Parsing the total-number of electrons and the Wigner-Seitz radius from the command line
  auto numElectrons = ueg::utility::ask_user_input<size_t>("Specify the number of electrons", closedShellElectrons);
  if (std::find(closedShellElectrons.begin(), closedShellElectrons.end(), numElectrons) == closedShellElectrons.end())
    throw std::invalid_argument("The number of electrons provided is not valid: " + std::to_string(numElectrons));

  auto wigner_seitz_radius = ueg::utility::ask_user_input<double>("Provide the Wigner-Seitz radius");
  if (wigner_seitz_radius <= 0.) { throw std::invalid_argument(" Wigner seitz radius can't be negative or zero! ");}

  // Calculate the lattice-constant and rescale the states with the radius
  lattice.reset_constant(numElectrons, wigner_seitz_radius);
  
  // Rescale all states based on (new) lattice-constant
  for (auto &state : states) state = state / lattice.constant;
  
  // Obtain the occupied and unoccupied orbitals
  occupied.insert(occupied.end(), states.begin(), states.begin()+numElectrons/2);
  virtuals.insert(virtuals.end(), states.begin()+numElectrons/2, states.end());
};

// ------------------------------------------------------------------------------------------------

System::System(std::string celltype, size_t numVirtuals, size_t numElectrons, double WignerSeitzRadius, dVec_t shift_) {
  
  // Check if the provided Wigner Seitz radius is correct.
  if ( WignerSeitzRadius <= 0.) { throw std::invalid_argument(" Wigner seitz radius can't be negative or zero! ");}
  
  // Set the attribute equal to the constructor parameter.  
  shift = shift_;
  auto numOccupied = numElectrons / 2;

  //  Take as maximum length of grid the radius of a sphere which ensues from the input of occupied and virtual orbitals.
  size_t maxGridLength = std::pow(5*(numOccupied+numVirtuals),0.3333);

  // Construct the lattice and get the states (k-points).
  lattice = Lattice(celltype, maxGridLength);
  auto states = lattice.getAllStates(lattice.constant,shift);
  
  // Check if everything went well with the shifted mesh.
  if (states.size() < numOccupied+numVirtuals) 
     throw std::invalid_argument("Something is wrong with the shifted mesh! Try a larger outer or a smaller inner box for the states."); 

  // Check if the number of electrons is valid (because of the random shift you might have a virtual state inside the Fermi "sphere").
  if (std::abs(states[numOccupied-1].squaredLength() - states[numOccupied].squaredLength()) < 1e-12) {
    throw std::invalid_argument("No differentiation between the occupied and virtuals. Try a different shift."); }

  // Calculate the lattice-constant and rescale the states with the radius
  lattice.reset_constant(numElectrons, WignerSeitzRadius);
  
  // Rescale all states based on (new) lattice-constant
  for (auto &state : states) state = state / lattice.constant;
  
  // Obtain the occupied and unoccupied orbitals
  occupied.insert(occupied.end(), states.begin(), states.begin()+numOccupied);
  virtuals.insert(virtuals.end(), states.begin()+numOccupied, states.begin()+numOccupied+numVirtuals);

  wigner_seitz_radius = WignerSeitzRadius;
};

// ------------------------------------------------------------------------------------------------

System::System(std::string celltype, int gridAxisSize) : lattice(celltype, gridAxisSize) {
  auto states = lattice.getAllStates(lattice.constant);
  
  // Obtain the available number of electons that can fill the system
  auto closedShellElectrons = get_available_electron_filling(states);

  // Parsing the total-number of electrons from the command line
  auto numElectrons = ueg::utility::ask_user_input<size_t>("Specify the number of electrons", closedShellElectrons);
  if (std::find(closedShellElectrons.begin(), closedShellElectrons.end(), numElectrons) == closedShellElectrons.end())
    throw std::invalid_argument("The number of electrons provided is not valid: " + std::to_string(numElectrons));
  auto wigner_seitz_radius = ueg::utility::ask_user_input<double>("Provide the Wigner-Seitz radius");
  if (wigner_seitz_radius <= 0.) { throw std::invalid_argument(" Wigner seitz radius can't be negative or zero! ");}

  // Calculate the lattice-constant and rescale the states with the radius
  lattice.reset_constant(numElectrons, wigner_seitz_radius);
  
  // Rescale all states based on (new) lattice-constant
  for (auto &state : states) state = state / lattice.constant;
  
  // Obtain the occupied and unoccupied orbitals
  occupied.insert(occupied.end(), states.begin(), states.begin()+numElectrons/2);
  virtuals.insert(virtuals.end(), states.begin()+numElectrons/2, states.end());
};

// ------------------------------------------------------------------------------------------------

System::System(std::string celltype, int gridAxisSize, size_t numElectrons, double WignerSeitzRadius) {

  if ( WignerSeitzRadius <= 0.) { throw std::invalid_argument(" Wigner seitz radius can't be negative or zero! ");}
  lattice = Lattice(celltype, gridAxisSize);
  auto states = lattice.getAllStates(lattice.constant);

  // Check if the provided number of electrons is valid
  auto closedShellElectrons = get_available_electron_filling(states);
  if (std::find(closedShellElectrons.begin(), closedShellElectrons.end(), numElectrons) == closedShellElectrons.end())
    throw std::invalid_argument("The number of electrons provided is not valid...");

  // Calculate the lattice-constant and rescale the states with the radius
  lattice.reset_constant(numElectrons, WignerSeitzRadius);
  
  // Rescale all states based on (new) lattice-constant
  for (auto &state : states) state = state / lattice.constant;

  
  // Obtain the occupied and unoccupied orbitals
  occupied.insert(occupied.end(), states.begin(), states.begin()+numElectrons/2);
  virtuals.insert(virtuals.end(), states.begin()+numElectrons/2, states.end());

  wigner_seitz_radius = WignerSeitzRadius;
};

// ------------------------------------------------------------------------------------------------

System::System(const Lattice &lat, size_t numElectrons, double WignerSeitzRadius) : lattice(lat) {

  if ( WignerSeitzRadius <= 0.) { throw std::invalid_argument(" Wigner seitz radius can't be negative or zero! ");}
  auto states = lattice.getAllStates(lattice.constant);

  // Check if the provided number of electrons is valid
  auto closedShellElectrons = get_available_electron_filling(states);
  if ( std::find(closedShellElectrons.begin(), closedShellElectrons.end(), numElectrons) == closedShellElectrons.end() )
    throw std::invalid_argument("The number of electrons provided is not valid: " + std::to_string(numElectrons));

  // Calculate the lattice-constant and rescale the states with the radius
  lattice.reset_constant(numElectrons, WignerSeitzRadius);
  
  // Rescale all states based on (new) lattice-constant
  for (auto &state : states) state = state / lattice.constant;
  
  // Obtain the occupied and unoccupied orbitals
  occupied.insert(occupied.end(), states.begin(), states.begin()+numElectrons/2);
  virtuals.insert(virtuals.end(), states.begin()+numElectrons/2, states.end());
  wigner_seitz_radius = WignerSeitzRadius;
};

// ------------------------------------------------------------------------------------------------

System::System(const System &sys) :
  lattice(sys.lattice), occupied(sys.occupied), virtuals(sys.virtuals),
  wigner_seitz_radius(sys.wigner_seitz_radius), madelung(sys.madelung), shift(sys.shift) {};

// ------------------------------------------------------------------------------------------------

std::vector<iVec_t> System::transform_to_integer_notation(std::string orbital_type) const {
  std::vector<iVec_t> result;
  auto coefficient = lattice.constant / (2*M_PI);

  if (orbital_type == "occupied") {
    result.reserve(occupied.size());
    for (auto const &state : occupied) 
      result.push_back( (coefficient * (state + (shift/coefficient))).round<int>() ); 
  } else if (orbital_type == "virtuals") {
    result.reserve(virtuals.size());
    for (auto const &state : virtuals) 
      result.push_back( (coefficient * (state + (shift/coefficient))).round<int>() ); 
  };

  return result;
};

// ------------------------------------------------------------------------------------------------

mapping_t System::get_ph_mapping() const {
  mapping_t mapping;

  auto occ = this->transform_to_integer_notation("occupied");
  auto vir = this->transform_to_integer_notation("virtuals");

  for (size_t i(0); i < occ.size(); ++i)
    for (size_t a(0); a < vir.size(); ++a)
      mapping[occ[i]-vir[a]].push_back(std::make_pair(i,a));

  return mapping;
};

// ------------------------------------------------------------------------------------------------

mapping_t System::get_pp_mapping() const {
  mapping_t mapping;

  auto vir = this->transform_to_integer_notation("virtuals");

  for (size_t a(0); a < vir.size(); ++a)
    for (size_t b(0); b < vir.size(); ++b)
      mapping[vir[a]+vir[b]].push_back(std::make_pair(a,b));

  return mapping;
};

// ------------------------------------------------------------------------------------------------

mapping_t System::get_hh_mapping() const {
  mapping_t mapping;

  auto occ = this->transform_to_integer_notation("occupied");

  for (size_t i(0); i < occ.size(); ++i)
    for (size_t k(0); k < occ.size(); ++k)
      mapping[occ[i]-occ[k]].push_back(std::make_pair(i,k));
    
  return mapping;
};

// ------------------------------------------------------------------------------------------------

mapping_3d_t System::get_hhh_mapping() const {
  mapping_3d_t mapping;

  auto occ = this->transform_to_integer_notation("occupied");

  for (size_t i(0); i < occ.size(); ++i)
    for (size_t j(0); j < occ.size(); ++j)
      for (size_t k(0); k < occ.size(); ++k)
        mapping[occ[i]+occ[j]+occ[k]].push_back(std::make_tuple(i,j,k));
    
  return mapping;
};

// ------------------------------------------------------------------------------------------------

mapping_3d_t System::get_ppp_mapping() const {
  mapping_3d_t mapping;

  auto vir = this->transform_to_integer_notation("virtuals");

  for (size_t a(0); a < vir.size(); ++a)
    for (size_t b(0); b < vir.size(); ++b)
      for (size_t c(0); c < vir.size(); ++c)
        mapping[vir[a]+vir[b]+vir[c]].push_back(std::make_tuple(a,b,c));

  return mapping;
};

// ------------------------------------------------------------------------------------------------

std::tuple<double,double> System::get_min_max_momentum_transfer() const {
  
  std::tuple<double,double> min_max_momentum_transfer;
  std::vector<double> momentum_transfer;

  for (auto const [kvector, connections]: this-> get_ph_mapping())
    momentum_transfer.push_back( (2*M_PI*kvector/this->lattice.constant).norm() );

  double q_min = *std::min_element(momentum_transfer.begin(), momentum_transfer.end());
  double q_max = *std::max_element(momentum_transfer.begin(), momentum_transfer.end());
  min_max_momentum_transfer = std::make_tuple(q_min,q_max);
  
  return min_max_momentum_transfer;
};

// ------------------------------------------------------------------------------------------------

inline double System::get_madelung_constant(int limit) const {
  double kappa = std::pow(lattice.boxVolume,-1.0/3.0);
  double term2 = -M_PI / (kappa*kappa*lattice.boxVolume);
  double term4 = -2 * kappa/std::sqrt(M_PI);
  double boxLength = 1.0/kappa;
  double recipsum = 0.0;
  double realsum = 0.0;

  for (int l1=-limit; l1 < limit; ++l1) {
    for (int l2=-limit; l2 < limit; ++l2) {
      for (int l3=-limit; l3 < limit; ++l3) {
        double n2 = l1*l1 + l2*l2 + l3*l3;
        double modr = boxLength * std::sqrt(n2);
        double k2 = kappa*kappa*n2;
        if (n2 > 0)
          recipsum += 1.0/(M_PI*k2)*std::exp(-M_PI*M_PI*k2/kappa/kappa)/lattice.boxVolume;
        if (modr > 0)
          realsum += erfc(kappa*modr)/modr;
      }
    }
  }
  return realsum + term2 + term4 + recipsum;
};

// ------------------------------------------------------------------------------------------------

double System::calculate_madelung_constant(double threshold) const {
  std::cout << "- Evaluating the Madelung constant (iterative convergence)" << std::endl;
  bool run_convergence = true;
  int limit = 1;
  auto reference_point = this->get_madelung_constant(limit);
  
  while (run_convergence) {
    limit += 1;
    auto temporal = this->get_madelung_constant(limit);
    if (std::abs(temporal-reference_point) < threshold)
      run_convergence = false;
    else
      std::cout << "\t- Converging the Madelung constant = " << reference_point << " (Limit = " << limit << ")" << std::endl;
    reference_point = temporal;
  };
  std::cout << "\t- Madelung constant converged to the value: " << reference_point << std::endl;

  return reference_point;
};

void System::information() const {
  std::cout << "\n( System information )" << std::endl;
  std::cout << "- Lattice constant in Bohr units: " << lattice.constant << std::endl;
  std::cout << "- Lattice constant in Angstrom  : " << 0.529*lattice.constant << std::endl;
  std::cout << "- Box-volume in Bohr units      : " << lattice.boxVolume << std::endl;
  std::cout << "- Box-volume in cubic Angstrom  : " << std::pow(0.529,3) * lattice.boxVolume << std::endl;
  std::cout << "- Cell volume                   : " << lattice.volume << std::endl;
  std::cout << "- Madelung constant             : " << madelung << std::endl;
  std::cout << "- Occupied orbitals             : " << occupied.size() << std::endl;
  std::cout << "- Virtuals orbitals             : " << virtuals.size() << std::endl;
  std::cout << "                                  " << std::endl;
};

std::vector<size_t> lattice_available_electron_filling(std::string celltype, int gridAxisSize) {
  auto lattice = Lattice(celltype, gridAxisSize);
  auto states = lattice.getAllStates(1); // 1 is the normalization
  return get_available_electron_filling(states);
};

}; // end namespace constructors
}; // end namespace ueg

