#pragma once 
#include <iostream>
#include <string>
#include <chrono>
#include <utility> 
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <fstream>
#include <iomanip>
#include <armadillo>

// IMPORTANT NOTE: 
// - As you will see all functions that are not 'templates' are defined as inline. This is 
//   to assure the compiler that he might receive multiple definitions, thus allowing us 
//   to include this header file into multiple other files

namespace ueg {
  namespace utility {

    // Type-definition
    typedef std::vector<std::tuple<double,double,double,double>> vector_tuples_t;

    // Template routine to request a value from a user 
    template<typename DataType>
    DataType ask_user_input(std::string message, const std::vector<DataType>& available_options = std::vector<DataType>()) { 
      DataType user_input;
      if ( !available_options.empty() ) { 
        bool provide_option = true; 
        std::cout << "[USER-INPUT]: " << message << std::endl;
        std::cout << " --> Available options are: \n\t";
        std::for_each(available_options.begin(),available_options.begin()+20,[](const DataType& el) { std::cout << el << ", "; });
        if ( available_options.size() > 20 ) std::cout << " ... "; 
        while ( provide_option ) { 
          std::cout << "\n--> Select from the above options: "; std::cin >> user_input; 
          if ( std::find(available_options.begin(), available_options.end(), user_input) != available_options.end() )
            provide_option = false;
        }; 
      } else {   
        std::cout << "[USER-INPUT]: " << message << ": "; std::cin >> user_input;
      };
      return user_input; 
    };

    // Routine to transform the energy units from [Ha/N] to [eV]
    inline double hartree_to_ev(int n_occupied, double energy) { return 2*27.211386*n_occupied*energy; };  
    
    
    // Routine to print messages 
    inline void message(std::string msg, bool verbosity=true) { 
      if ( verbosity ) std::cout << msg << std::endl; 
    }; 


    // Routine to check if the memory limit is compromised
    inline bool is_memory_enough(size_t n_elements, size_t n_objects, double threshold) {
      double memory =  n_objects*n_elements*sizeof(double)/std::pow(1024.0,3);
      std::cout << "[INFO]: Estimate for the memory required: " << memory <<  " GB" << std::endl;
      return memory < threshold;
    };
    
    // Routine to write a vector of 4-component tuples in a file
    inline void write_to_file(std::string filename, vector_tuples_t results, std::string suffix=".txt") {
      auto has_extension = [] (std::string n, std::string s) {
         return ( n.size() >= s.size() && (n.compare(n.size()-s.size(), s.size(), s))==0 );
      };
      
      if ( !has_extension(filename,suffix) ) { filename = filename + suffix; };
    
      // Print name of file
      std::cout << "[INFO]: Writing results to the file: " << filename << "\n" << std::endl;
   
      // Open - Write results - Close 
      std::ofstream output; output.open(filename);
      for (auto const &[t1, t2, t3, t4] : results)
        output << std::fixed << std::setprecision(8) << t1 << "\t" << t2 << "\t"
               << t3 << "\t" << t4 << std::endl;
      output.close();
    };

    inline void store_to_file(std::string filename, arma::vec results, std::string suffix=".txt") {
      auto has_extension = [] (std::string n, std::string s) {
         return ( n.size() >= s.size() && (n.compare(n.size()-s.size(), s.size(), s))==0 );
      };
      
      if ( !has_extension(filename,suffix) ) { filename = filename + suffix; };
    
      // Print name of file
      std::cout << "-Writing results to the file: " << filename << std::endl;
   
      // Open - Write results - Close 
      std::ofstream output; output.open(filename);
      for (auto const &data : results)
        output << std::fixed << std::setprecision(8) << data << "\t" << std::endl;
      output.close();
    };

    // Function to time a callback ( variadic template quite difficult ) 
    template<typename TimeUnit, typename FunctionType, typename... Args>
    double get_execution_time(FunctionType callback_function, Args&&... args) {
      std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      callback_function(std::forward<Args>(args)...);
      return std::chrono::duration_cast<TimeUnit>(std::chrono::high_resolution_clock::now()-t1).count();
    };

    
    // Time convensions to string format 
    template<typename TimeUnit> struct time_unit_as_string {}; 
    template<> struct time_unit_as_string<std::chrono::seconds> { inline static const std::string unit = " s";       }; 
    template<> struct time_unit_as_string<std::chrono::milliseconds> { inline static const std::string unit = " ms";       }; 
    template<> struct time_unit_as_string<std::chrono::microseconds> { inline static const std::string unit = " \u03BCs"; }; 
    template<> struct time_unit_as_string<std::chrono::nanoseconds>  { inline static const std::string unit = " ns";       };

  }; // end namespace utility
}; // end namespace ueg 


// #define is just a preprocessor directive. Here what it does is to 
// substitute all names on the left with the expressions on the right 
// Quite useful and old ( sure there is a better way to do it )
#define R       "Rings"
#define QR      "Ring-Quadratic"
#define DRCC    "Direct-Ring-CC" 
#define PPL     "PP-ladder"
#define HHL     "HH-ladder"
#define VPHL    "PH-vertex-ladder"
#define VVPHL   "PH-vertices-ladder"
#define HHQ     "Xijkl_Tklab-quadratic"
#define VVPHLQ  "Xcibk_Tkjac-quadratic"
#define VPHLQ   "Xciak_Tkjcb-quadratic"
#define XTBC    "Xicak_Tkjbc"
#define XTCB    "Xicak_Tkjcb"
#define LIK     "Lik_Tkjab"
#define LCA     "Lca_Tijcb"
#define XTCBL   "Xicak_Tkjcb-linear"
#define XTCBNA  "Xicak_Tkjcb-non-antisym"
#define XTCBA2  "Xicak_Tkjcb-antisym-2"
#define XTCBA   "Xicak_Tkjcb-antisym"
#define XTCBQ   "Xicak_Tkjcb-quadratics"
#define XTBCL   "Xicak_Tkjbc-linear"
#define XTBCNA  "Xicak_Tkjbc-non-antisym"
#define XTBCA   "Xicak_Tkjbc-antisym"
#define XTBCQ   "Xicak_Tkjbc-quadratics"

// This is really just to minimize the linewidth of the declarations 
#define ISC_SVEC_T inline static const std::vector<std::string> 

namespace ueg { 
  namespace utility { 
    namespace diagrams { 
      // Inside here place all the convensions you want to use. 
      // Example: auto selected_diagrams = ueg::utility::diagrams::symmetrizable::identifiers; 
      struct ccd { ISC_SVEC_T identifiers { PPL, HHL, VPHL, VVPHL, HHQ, VVPHLQ, VPHLQ, XTBC, XTCB, LIK, LCA }; };
      struct drccd { ISC_SVEC_T identifiers {DRCC}; };
      struct ldrccd { ISC_SVEC_T identifiers {R}; };
      struct smp3 {ISC_SVEC_T identifiers {R, VPHL, VVPHL, XTBCL}; };
      struct nsmp3 { ISC_SVEC_T identifiers {PPL, HHL}; };
      struct smp4 { ISC_SVEC_T identifiers {QR, LIK, LCA, XTCBNA, XTCBA2, XTBCNA, XTBCA, VPHLQ, VVPHLQ }; };
      struct nsmp4 { ISC_SVEC_T identifiers {HHQ}; };
      struct symmetrizable { ISC_SVEC_T identifiers { VVPHL, VPHL, VVPHLQ, VPHLQ, XTCBL, XTBCL, XTBCNA, XTBCA, XTCBNA, XTCBA2, QR, LIK, LCA }; }; 
      struct unsymmetrizable { ISC_SVEC_T identifiers { PPL, HHL, HHQ }; }; 
      struct all_diagrams { ISC_SVEC_T identifiers { VPHL, VVPHL, VVPHLQ, VPHLQ, XTBCL, XTBCNA, XTBCA, XTCBL, XTCBNA, XTCBA2, QR, LIK, LCA, PPL, HHL, HHQ }; }; 
    
    }; // end namespace diagrams 
  }; // end namespace utility
}; // end namespace ueg 

// For safety always undefine the 'definitions'
#undef R
#undef QR
#undef DRCC
#undef PPL
#undef HHL
#undef VPHL
#undef VVPHL
#undef HHQ
#undef VVPHLQ
#undef VPHLQ
#undef XTBC
#undef XTCB
#undef LIK
#undef LCA
#undef XTCBL
#undef XTCBNA
#undef XTCBA
#undef XTCBQ
#undef XTBCL
#undef XTBCNA
#undef XTBCA
#undef XTBCQ
#undef ISC_SVEC_T
