// --------------------------------------------------------------------------------------
// SCOPE: This is a simple program that can be used to make a full-convergent calcution
//        It allows for the following:
//        - Setting once and for all the lattice-type               (from command-line)
//        - Setting once and for all the number of electrons        (from command-line)
//        - Setting once and for all the Wigner-Seitz radius        (from command-line)
//        - Setting multiple grid-axis dimensions  to scan          (from command-line)
//        - Setting the requested convergence precision             (from command-line)
//        - Setting the maximum-iterations allowed                  (from command-line)
//        - Submission of the calculation of the structure factor   (from command-line)
//        - Choose the desirable diagram for the calculation        (from command-line)
// --------------------------------------------------------------------------------------

#include "ueg.hpp"
#include "parser.hpp"
#include <typeinfo>

// CHANGE THE TIME-UNITS HERE IF YOU WANT
#define time_unit_t std::chrono::milliseconds

// DEFINE THE MAXIMUM AVAILABLE MEMORY
#define max_memory 200

int main(int argc, char* argv[]) {

  // Please place this function always in the beginning of your programs
  //ueg::version::print_information();
  //print_application_version();

  // Here we parse the command-line arguments (see SCOPE)
  auto [cell, grid_L, grid_S, n_elec, radius, method, decomposition, precision, max_iterations, structure_factor, diagrams, shift_, triples_, crtriples_] = parse_arguments(argc,argv);

  // Decompose the parsing values of the grid
  size_t grid_start, grid_end, grid_step;
  std::tie(grid_start,grid_end,grid_step) = grid_decomposition(grid_S);
 
  // Now we have everything to start the main-part of the program
  for (size_t grid_axis = grid_start; grid_axis <= grid_end; grid_axis += grid_step) {
    std::cout << "\n\n" << std::string(100,'=') << " New calculation \n\n" << std::endl;
    
    // Get the states that are contained in the large box. The boundaries are defined by the user (grid_L).
    // The final states for the calculation will be set up in the system constructor.
    auto getStates = pseudoStates(grid_L);
    
    // Check if the occupied and total states for the calculation are reserved.
    checkStatesConstruction(n_elec,grid_axis,getStates);
 
    // Obtain the correct number of the virtuals (full shell).
    auto Nv = numVirtuals(n_elec,grid_axis,getStates);
    
    // Obtain the shifted k-mesh.
    auto shift = kmesh_shift(shift_);

    // Create the corresponding system + and set the madelung constant
    auto system = ueg::constructors::System(cell,Nv,n_elec,radius,shift);
    //auto system = ueg::constructors::System(cell,grid_axis,n_elec,radius);
    system.madelung = system.calculate_madelung_constant(1e-8);
    
    // Print some necessary info of the system and the memory requirements for the calculation.
    print_system_info(system); check_memory_requirements(system, max_memory);
    
    // Construct the Correlator object.
    auto correlator = ueg::correlators::Correlator<hartree_fock_t>(system, "interacting");
    correlator.verbosity = false;
     
    // Get some parameters of the system (number of electrons, total number of orbitals and chosen method) 
    // in a concatenated string form.
    auto system_parameters = parameters_concatenation(system, n_elec, method);
 
    // Calculate the MP2 Energy of the given system.
    calculate_mp2_energy<time_unit_t>(correlator, method);
 
    // Now we have everything to run a self-consistent calculation.
    run_self_consistency<time_unit_t>(correlator, precision, max_iterations, method);
  
    // Define boolean variables and obtain the values that are parsed for the
    // perturbative triples energy and structure factor calculation. 
    //correlator._triples("all");

    bool triples, triples_sf; 
    std::tie(triples,triples_sf) = beyond_ccd_info(triples_);

    // Get the parameters of the system for naming the structure factor output file.
    auto triples_parameters = parameters_concatenation(system, n_elec, "triples");
    
    // Perturbative triples energy and structure factor (upon request).
    if (triples)
      run_triples<time_unit_t>(correlator, triples_sf, triples_parameters);
    
    // Define boolean variables and obtain the values that are parsed for the complete
    // renormalized and perturbative triples energies and structure factors calculation. 
    bool crtriples, crtriples_sf;
    std::tie(crtriples,crtriples_sf) = beyond_ccd_info(crtriples_);
    
    // Get the parameters of the system for naming the structure factor output file.
    auto crtriples_parameters = parameters_concatenation(system, n_elec, "crtriples");

    // Perturbative triples and Complete renormalized triples energies and structure factors (upon request).
    if (crtriples)
      run_complete_renormalized<time_unit_t>(correlator, crtriples_sf, crtriples_parameters);

    // Store the CCD diagrams in a vector of strings for presumptive decomposition. 
    auto ccd_diagrams = ueg::utility::diagrams::all_diagrams::identifiers;

    // Observables (upon request) - Structure Factor and method/channel decomposition.
    if (structure_factor) {
      if (std::find(diagrams.begin(), diagrams.end(), "channels") != diagrams.end()) 
        structure_factor_calculation(correlator,system_parameters,ccd_diagrams);
      structure_factor_calculation(correlator,system_parameters,diagrams);}
    
    if (decomposition){
      if (std::find(diagrams.begin(), diagrams.end(), "channels") != diagrams.end()) 
        energy_decomposition<time_unit_t>(correlator,ccd_diagrams); 
      energy_decomposition<time_unit_t>(correlator, diagrams); }
  };
 
  return 0;
}
