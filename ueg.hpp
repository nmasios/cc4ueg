#pragma once
#include <iostream>
#include <chrono>
#include <iomanip>
#include <iterator>
#include <stdio.h>
#include <tuple>
//#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

// Note that I don't need to write relative paths (already set from cmake)
#include "correlators.hpp"

// Minimize your code by simplifying notations
typedef ueg::diagrammatics::HartreeFock hartree_fock_t;
typedef ueg::correlators::Correlator<hartree_fock_t> correlator_t;

// On each iteration the correlator is updated and we need to keep
// track of the direct and exchange energy and the elapsed time
typedef std::tuple<double,double,double> iteration_results_t;
typedef ueg::fundamentals::Vector<double> dVec_t;

// Print-functions for nice-looks forward declaration (implementation in the end)
void print_preamble();
void print_system_info(ueg::constructors::System &system);
void print_iteration_results(int inum, double d, double e, double t, double diff, double dt, std::string unit);
void print_energy_results(correlator_t &correlator, int iteration_number, double direct, double exchange, double total, std::string method);
void check_memory_requirements(ueg::constructors::System &system, double memory);
iteration_results_t self_consistent_cycle(correlator_t &correlator, std::string method_option);
void energy_decomposition(correlator_t &correlator, std::vector<std::string> &diagrams);
void calculate_mp2_energy(correlator_t &correlator,std::string method_option);
void run_self_consistency(correlator_t &correlator, double &precision, int &max_iterations, std::string option);
void run_triples(correlator_t &correlator, bool structure_factor, std::string system_parameters);
void run_complete_renormalized(correlator_t &correlator, bool structure_factor, std::string system_parameters);
void structure_factor_calculation(correlator_t &correlator, std::string system_parameters, std::vector<std::string> diagrams);
dVec_t kmesh_shift(const std::vector<double> &shift_info);
std::tuple<size_t,size_t,size_t> grid_decomposition(const std::vector<size_t> &grid_info);
std::string parameters_concatenation(ueg::constructors::System &system, size_t n_electrons, std::string method);
std::vector<int> pseudoStates(const int &outerLimit);
void checkStatesConstruction(size_t n_elec, size_t innerLimit, std::vector<int> pseudoStates);
size_t numVirtuals(size_t n_elec, size_t innerLimit, std::vector<int> pseudoStates);
std::tuple<bool,bool> beyond_ccd_info(const std::vector<bool> &info);

// Use of lambda function to track the given method.
auto find_method = [] (std::vector<std::string> method_class, std::string method_option) {
  return ( std::find( method_class.begin(), method_class.end(), method_option) != method_class.end() );
};

// Definition of different methods that are used in the calculations.
std::vector<std::string> self_consistent_methods = {"CCD", "drCCD", "ldrCCD"};
std::vector<std::string> one_shot_methods = {"MP3", "MP4"};

// First I define everything that happens during one only iteration
template<typename TimeUnit>
iteration_results_t self_consistent_cycle(correlator_t &correlator, std::string method_option) {
  
  // Set the timer.
  auto t1 = std::chrono::high_resolution_clock::now();

  if ( find_method(self_consistent_methods, method_option) ) {

    // Calculate the contribution of the corresponding diagrams to the residuum.
    correlator.self_consistent_cycle_body(method_option);
   
    // Calculate the amplitude and compute the "error" with the preceding one.
    auto pair_tensor = correlator.calculate_errors_and_amplitudes();
  
    // When the desired number of amplitudes and errors will be stored the DIIS extrapolation will take place.
    correlator.diis_extrapolation(pair_tensor);

  } else if ( find_method (one_shot_methods, method_option) ) {
   
    // Calculate the one shot method energy.
    correlator.one_shot_methods(method_option);

  } else { throw std::invalid_argument ("Class of method wasn't found! " + method_option); exit(1); }

  // Get total elapsed time
  auto elapsed_time = std::chrono::duration_cast<TimeUnit>(std::chrono::high_resolution_clock::now()-t1).count();

  return std::make_tuple(correlator.getEnergy("direct"),correlator.getEnergy("exchange"),elapsed_time);
};

// Routine to obtain the direct, exchange and total energies of the desired decomposed channels.
template<typename TimeUnit>
void energy_decomposition(correlator_t &correlator, std::vector<std::string> &diagrams) {
 
  // Declare the decomposed energies.
  double d_direct, d_exchange, d_total_energy;
  
  // Loop through the list of the diagrams.
  for (auto const &diagram : diagrams) {
    
    // Brute force implementation of channels' decomposition.
    if(diagram == "channels") continue;

    // Obtain the amplitude of the corresponding diagram/diagrams from the decomposition.
    auto d_amplitude = correlator.decomposition(diagram);

    // Return the Direct, Exchange and total energy of the corresponding diagram/diagrams.
    auto results = std::make_tuple( correlator.getDecomposedEnergy(d_amplitude, "direct"), 
                                    correlator.getDecomposedEnergy(d_amplitude, "exchange"),
                                    correlator.getDecomposedEnergy(d_amplitude, "total") );
    // Unpack the results.
    std::tie(d_direct,d_exchange,d_total_energy)= results;
   
    // Print the results.
    printf("[INFO]: Decomposed energy contribution in Hartree/electron is = ( %s | %zu | %+10.08f | %+10.08f | %+10.14f )\n\n",
            diagram.c_str(), correlator.refpoint.n_virtuals, d_direct, d_exchange, d_total_energy);}
};

// Calculate MP2 energy.
template<typename TimeUnit>
void calculate_mp2_energy(correlator_t &correlator, std::string method_option) {

  // Set the timer and time-unit that we use for the elapsed time.
  auto tunit = ueg::utility::time_unit_as_string<TimeUnit>::unit;
  auto t1 = std::chrono::high_resolution_clock::now();

  if ( method_option == "MP2" || find_method(self_consistent_methods, method_option) || find_method(one_shot_methods, method_option) ) {

  // Declaration of variables to hold.
  double Direct, Exchange, ElapsedTime, difference, total_energy;

  // Provide the Hartree-Fock energy of the system.
  //correlator.refpoint.getHFenergy();

  // First we need to initialize the interaction and the propagator.
  correlator.initialize();

  // Since MP2 must be enabled we will do it directly to the amplitude.
  correlator.amplitude.set_from_product(correlator.interaction.values, correlator.propagator.values);

  // Return MP2 direct, exchange and total energy.
  auto results = std::make_tuple(correlator.getEnergy("direct"),correlator.getEnergy("exchange"),correlator.getEnergy("total"));

  // Unpack the results.
  std::tie(Direct,Exchange,total_energy)= results;

  // Print out the results of MP2 calculation.
  std::cout << std::string(153,'.') << std::endl;
  printf("[INFO]: MP2 energy calculated. Energy contributions in Hartree/electron are ( D | E | T | O ) = ( %+10.08f | %+10.08f | %+10.14f | %zu )\n",
                                                                                        Direct, Exchange, total_energy, correlator.refpoint.n_virtuals);
  // Get total elapsed time and print the results.
  std::cout << "" << std::endl;
  auto elapsed_time = std::chrono::duration_cast<TimeUnit>(std::chrono::high_resolution_clock::now()-t1).count();
  printf("[INFO]: Execution time of calculating MP2 energy: %.4f%s \n", double(elapsed_time), tunit.c_str());
  std::cout << std::string(153,'.') << "\n" << std::endl;}

  else {throw std::invalid_argument ("Class of method wasn't found! " + method_option); exit(1);}

};

// Define the full self-consistency
template<typename TimeUnit>
void run_self_consistency(correlator_t &correlator, double &precision, int &max_iterations, std::string option) {

  // Declaration of variables to hold
  double Direct, Exchange, ElapsedTime, difference, total_energy;
  int iteration_number = 1;

  // The time-unit that we use for the elapsed time
  auto tunit = ueg::utility::time_unit_as_string<TimeUnit>::unit;

  // Terminate if the option is the MP2 calculation 
  if (option == "MP2") {return;}

  if ( find_method (one_shot_methods, option) ) {

    // Unpack the results
    std::tie(Direct,Exchange,ElapsedTime) = self_consistent_cycle<TimeUnit>(correlator, option);
 
    // Get the total energy.
    total_energy = Direct + Exchange;

    // Print out the results of the one shot calculation
    print_energy_results(correlator, iteration_number, Direct, Exchange, total_energy, option);
    
    return; }

  // Store the amplitude and the error 
  auto pair_tensor = correlator.calculate_errors_and_amplitudes();

  // Calculate the total-energy
  print_preamble();

  // The whole convergence
  while ( (iteration_number == 1) || (std::abs(difference) > precision) ) {
    // Unpack the results
    std::tie(Direct,Exchange,ElapsedTime) = self_consistent_cycle<TimeUnit>(correlator, option);
    
    // Update the convergence criteria
    difference = ( Direct + Exchange ) - total_energy;
    total_energy += difference;

    // Print-out information on the loop
    print_iteration_results(iteration_number,Direct,Exchange,total_energy,difference,ElapsedTime,tunit);
    
    // Check if we have exceeded maximum_iterations
    if (iteration_number == max_iterations) break;
    else iteration_number++;
  }

  // Closing message (can be greped easily)
  if ( (iteration_number == max_iterations) && (std::abs(difference) > precision) ) {
    std::cout << std::string(153,'.') << std::endl;
    std::cout << "[INFO]: Unable to reach convergence..." << std::endl;  exit(1); 
  } else print_energy_results(correlator, iteration_number, Direct, Exchange, total_energy,option);

};

// Calculates the (T) energy and (if required) the (T) structure factor. Writes the results 
// in a txt file as (momentum transfer q, structure factor, FT Coulomb interaction).
template<typename TimeUnit>
void run_triples(correlator_t &correlator, bool structure_factor, std::string system_parameters) {

  // Set the timer and time-unit that we use for the elapsed time.
  auto tunit = ueg::utility::time_unit_as_string<TimeUnit>::unit;
  auto t1 = std::chrono::high_resolution_clock::now();
  
  // Start the evaluation.
  if(structure_factor) {
    auto sf_results = correlator._OTF_CCDT(structure_factor);
    ueg::utility::write_to_file("SF" + system_parameters, sf_results, ".txt"); }
  else 
    correlator._OTF_CCDT(structure_factor);

  // Get total elapsed time and print the results.
  auto elapsed_time = std::chrono::duration_cast<TimeUnit>(std::chrono::high_resolution_clock::now()-t1).count();
  printf("[INFO]: Execution time of calculating perturbative triples correction energy: %.4f%s \n", double(elapsed_time), tunit.c_str());
  std::cout << "" << std::endl;

};

// Calculates the (T) and CR-(T) energies and (if required) the corresponding structure factors. Writes the
// results in a txt file as (momentum transfer q, (T) structure factor, CR-(T) structure factor, FT COulomb interaction).
template<typename TimeUnit>
void run_complete_renormalized(correlator_t &correlator, bool structure_factor, std::string system_parameters) {

  // Set the timer and time-unit that we use for the elapsed time.
  auto tunit = ueg::utility::time_unit_as_string<TimeUnit>::unit;
  auto t1 = std::chrono::high_resolution_clock::now();

  // Start the evaluation.
  if(structure_factor) {
    auto sf_results = correlator._OTF_CCDcT(structure_factor);
    ueg::utility::write_to_file("SF" + system_parameters, sf_results, ".txt"); }
  else 
    correlator._OTF_CCDcT(structure_factor);

  // Get total elapsed time and print the results.
  auto elapsed_time = std::chrono::duration_cast<TimeUnit>(std::chrono::high_resolution_clock::now()-t1).count();
  printf("[INFO]: Execution time of calculating (T) and CR-(T) energy: %.4f%s \n", double(elapsed_time), tunit.c_str());
  std::cout << "" << std::endl;

};

// How the preamble of self-consistent results is printed. 
void print_preamble() {
  std::cout << "\n " << std::setw(99) << "Information on the self-consistency procedure\n"
      << std::string(153,'-') << "\n" << std::string(23,' ')
      << "Direct-energy          Exchange-energy           Total-Energy             Difference                Elapsed-time\n"
      << std::string(153,'-') << std::endl;
};

// Return some system information. It can be modified according to the user's taste.
void print_system_info(ueg::constructors::System &system) {
  std::cout << "\n- General information of the system" << std::endl;
  std::cout << "\t- Lattice type                  : " << system.lattice.type << std::endl;
  std::cout << "\t- Lattice constant in Bohr units: " << system.lattice.constant << std::endl;
  std::cout << "\t- Lattice constant in Angstrom  : " << 0.529*system.lattice.constant << std::endl;
  std::cout << "\t- Box-volume in Bohr units      : " << system.lattice.boxVolume << std::endl;
  std::cout << "\t- Box-volume in cubic Angstrom  : " << std::pow(0.529,3) * system.lattice.boxVolume << std::endl;
  std::cout << "\t- Cell volume                   : " << system.lattice.volume << std::endl;
  std::cout << "\t- WS radius in Borh units       : " << system.wigner_seitz_radius << std::endl;
  std::cout << "\t- Madelung constant             : " << system.madelung << std::endl;
  std::cout << "\t- Occupied orbitals             : " << system.occupied.size() << std::endl;
  std::cout << "\t- Virtual  orbitals             : " << system.virtuals.size() << std::endl;
  std::cout << "\t- K-mesh shift                  : " << system.shift << std::endl;
  std::cout << "                                  " << std::endl;
  std::cout << std::string(153,'.') << std::endl;
};

// Formation of each self-consistent cycle results.
void print_iteration_results(int inum, double d, double e, double t, double diff, double dt, std::string unit) {
  printf("Iteration-%03d     %4s%+14.12f         %+14.12f         %+14.12f         %+14.12f         % 14.4f", inum, "", d, e, t, diff, dt);
  std::cout << unit << std::endl;
};
 
// Output of corresponding calculation (can be grepped easily).
void print_energy_results(correlator_t &correlator, int iteration_number, double direct, double exchange, double total, std::string method) {
    std::cout << std::string(153,'.') << std::endl;
    printf("[INFO]: %s Energy calculated. Energy contributions in Hartree/electron are ( D | E | T | O ) = ( %+10.08f | %+10.08f | %+10.14f | %zu )\n",
            method.c_str(), direct, exchange, total, correlator.refpoint.n_virtuals);
    std::cout << std::string(153,'.') << "\n" << std::endl;
};

// Functions that checks (approximately) the memory requirements according to the number of
// tensors which are used during a calculation and their size (No*No*Nv).
void check_memory_requirements(ueg::constructors::System &system, double memory) {

  // Calculate the size of a tensor used in the calculations (usually No*No*Nv).
   auto tensor_size = system.virtuals.size() * std::pow(system.occupied.size(),2);

  // Calculate approximately the requisite memory of the calculation and compare to the memory capacity.
  if (!ueg::utility::is_memory_enough(tensor_size, 12, memory)) {
    throw std::invalid_argument("--> Memory requirements exceed available resources. Please try a different calculation."); exit(1);
  };
};

// Function that checks the parsed components of the shifted k-mesh (the norm of each
// component shoulb be smaller than or equal to 0.5) and returns the values in a dVec_t type.
dVec_t kmesh_shift(const std::vector<double> &shift_info){
  
  dVec_t shift;
  double shift_length, threshold(0.5);
  double x_norm = std::abs(shift_info[0]), y_norm = std::abs(shift_info[1]), z_norm = std::abs(shift_info[2]);
  
  if ( shift_info.size() == 3 && x_norm <= 0.5 && y_norm <= threshold && z_norm <= threshold)
    shift = { shift_info[0], shift_info[1], shift_info[2] };
  else throw std::invalid_argument("The number of elements in the k-mesh shift vector should be 3 and each component in the interval [-0.5,0.5]. Please try again.");

  return shift;
}

// Function that 'decomposes' a grid. User can give an initial starting point, 
// a final one and the desirable step for continuous calculations.
std::tuple<size_t,size_t,size_t> grid_decomposition(const std::vector<size_t> &grid_info) {

  // Define the parameters.
  size_t grid_start, grid_end, grid_step;
  
  // Distinguish between the following cases.
  if      (grid_info.size() == 1) { grid_start = grid_info[0]; grid_end = grid_info[0]; grid_step = 1; }
  else if (grid_info.size() == 2) { grid_start = grid_info[0]; grid_end = grid_info[1]; grid_step = 1; }
  else if (grid_info.size() == 3) { grid_start = grid_info[0]; grid_end = grid_info[1]; grid_step = grid_info[2]; }
  else { throw std::invalid_argument("Wrong grid initialization."); }
  
  // Return a tuple with the analyzed grid parameters.
  return std::make_tuple(grid_start,grid_end,grid_step);
};

// Build pseudo-states which are included in a very large box. The limits are given by the user
// and the states are projected to an integer number.
std::vector<int> pseudoStates(const int &outerLimit) {

  // Define a vector which will contain all the possible states inside a large box. 
  std::vector<int> maximumStates;

  // Build the pseudo-states.
  for(int m1 = -outerLimit; m1 <= outerLimit; ++m1)
  for(int m2 = -outerLimit; m2 <= outerLimit; ++m2)
  for(int m3 = -outerLimit; m3 <= outerLimit; ++m3) 
    maximumStates.push_back(m1*m1 + m2*m2 + m3*m3);
  
  // Sort the pseudo-states.
  std::sort(maximumStates.begin(), maximumStates.end());

  return maximumStates; 
  
};

// Check if everything went well during the construction of pseudo-states. The user will parse a secondary
// point on which a smaller box of states will be constructed.
void checkStatesConstruction(size_t n_elec, size_t innerLimit, std::vector<int> pseudoStates) {

  size_t No(n_elec/2), Nv(innerLimit*n_elec/2);

  // Check if the given number of electrons is correct.
  if (pseudoStates[n_elec/2] == pseudoStates[n_elec/2 - 1])
    throw std::invalid_argument("The number of electrons (closed-shell filling) is not correct! Please try again.");

  // Check if the number of the built pseudo-states is smaller than the number 
  // of the states built by the user in the secondary smaller box.
  if (pseudoStates.size() < No+Nv)
    throw std::invalid_argument("The number of the approximate inner states is outside the larger boundary. Please try again.");
};

// Returns the number of virtual orbitals/states according to the inner box boundaries.
size_t numVirtuals(size_t n_elec, size_t innerLimit, std::vector<int> pseudoStates) {
  
  size_t Nv_init(innerLimit*n_elec/2), Nv_final;
  
  // Search for all the states that form the corresponding shell.
  for (size_t num(Nv_init); num < pseudoStates.size(); ++num) {
    if (pseudoStates[num] != pseudoStates[num+1]) {
      Nv_final = num+1-n_elec/2; break; }
  }

  return Nv_final; 
};

// Function that concatenates in a string the chosen system's parameters. 
// It can be extended according to the user's taste. 
// Usage: File name of strructure factor results.
std::string parameters_concatenation(ueg::constructors::System &system, size_t n_electrons, std::string method) {

  // Define the variables for the desired system's information.
  std::string full_info, ldash("_");

  // Get the total number of orbitals.
  auto orbitals = system.occupied.size() + system.virtuals.size();
  
  // Concatenate data.
  full_info.append( ldash + std::to_string(n_electrons) + ldash + std::to_string(orbitals) + ldash + method );
 
  return full_info; 
};

// Function for calculating the structure factor of a method (MP2, MP3, MP4, DRCCD, CCD).
// Capability of calculating the structure factor of different channels through decomposition.
// Returns a txt file as "SF_+system's information_+method/decomposed channel".
void structure_factor_calculation(correlator_t &correlator, std::string system_parameters, std::vector<std::string> diagrams) {

  // Define the variables for the full system's information.
  std::string full_info, ldash("_");
   
  // Start the structure factor's calculation.  
  for(auto const &diagram: diagrams) {

    // Brute force implementation of channels' decomposed structure factors.
    if (diagram == "channels") continue;

    std::cout << "[INFO]: Starting the calculation of the Structure Factor for "+diagram << "\n" << std::endl;
    auto decomposed_amplitude = correlator.decomposition(diagram);
    auto structure_factor_results = correlator.structure_factor(decomposed_amplitude);
    full_info = system_parameters + ldash + diagram;
    ueg::utility::write_to_file("SF" + full_info, structure_factor_results, ".txt"); }

};

// Returns in a tuple boolean values that are parsed for the calculation of the energy 
// and structure factor of either perturbative triples or complete renormalized triples.
std::tuple<bool,bool> beyond_ccd_info(const std::vector<bool> &info){

  // Define the parameters.
  bool energy, structure_factor;
  
  // Distinguish between the following cases.
  if      (info.size() == 2) { energy = info[0]; structure_factor = info[1]; }
  else { throw std::invalid_argument("Only 2 boolean arguments for the calculation. Please try again."); }
  
  // Return a tuple with the analyzed requirements.
  return std::make_tuple(energy,structure_factor);

};
