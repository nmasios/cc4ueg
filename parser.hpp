#pragma once 
#include<iostream>
#include<vector>
#include<string>
#include<tuple>
#include "cxxopts.hpp"

// Create your own beautiful parser using the cxxopts!!
// ----------------------------------------------------
// The only things that you have to keep in mind are:
// - which parameters you would like to parse from the command line 
// - what is the type of those parameters 
// - in which sequence you want them to be returned ( and build a tuple ) 
//
// The parser below for example is made under the following convensions:
// - Parameters to parse: cell, grid, electrons, radim, threshold, maxiter
// - Types of them: string, vector<int>, size_t, double, double, int 
// - Sequence to return: same as above 
//
// For convenience I (always) define the following type 
#include "fundamentals.hpp"

typedef ueg::fundamentals::Vector<double> dVec_t;
typedef std::tuple< std::string, int, std::vector<size_t>, size_t, double, std::string, bool, double, int, bool, std::vector<std::string>,
                    std::vector<double>, std::vector<bool>, std::vector<bool> >  arguments_t; 


// Now I define the function that I am going to call in the main program 
arguments_t parse_arguments(int argc, char* argv[]) { 
  // First I create the parser object ( using the library )
  cxxopts::Options parser("Skeleton-example", "Application for the Structure Factor"); 

  // Now we need to define the named-arguments that the parser should expect 
  // to receive together with the types that they need to be returned. 
  parser.add_options()
    ("c,cell", "The type of the lattice-cell that we use", cxxopts::value<std::string>())
    ("G,Grid", "The grid dimensions to build the boundaries of the large box. Suggested value:30", cxxopts::value<int>())
    ("g,grid", "The smallest grid dimensions that we want for the calculation", cxxopts::value<std::vector<size_t>>())
    ("e,electrons", "The number of electrons in the system", cxxopts::value<size_t>()) 
    ("r,radius", "The Wigner-Seitz radius of the system", cxxopts::value<double>())
    ("m,method", "Chosen method (MP2,DRCC,CCD,DCC)", cxxopts::value<std::string>())
    ("d,decomposition", "Energy decomposition to the different channels of CCD", cxxopts::value<bool>())
    ("p,precision", "The convergence limit that we require", cxxopts::value<double>())
    ("i,iterations"  , "The maximum-allowed iterations allowed", cxxopts::value<int>())
    ("f,structureFactor" , "Calculate Structure Factor (true/false)", cxxopts::value<bool>())
    ("t,terms"  , "Choose term/diagrams for the calculation of the SF ", cxxopts::value<std::vector<std::string>>()->default_value("none"))
    ("s,shift"  , "Choose shift for a shifted k-mesh calculation  ", cxxopts::value<std::vector<double>>())
    ("T,Triples", "Calculate (T) energy and structure factor. Tuple of boolean values (energy,SF)", cxxopts::value<std::vector<bool>>())
    ("C,CRtriples", "Calculate CR-(T) energy and structure factor. Tuple of boolean values (energy,SF)", cxxopts::value<std::vector<bool>>())
    ("h,help"     , "Print-out help message")
    ; // Here the command for adding options ends 

  // Parse the arguments and store them into a variable
  auto arguments = parser.parse(argc, argv); 

  // In case that you don't remember the names of the arguments just type 'help'
  if ((arguments.count("help") == 1) || (arguments.count("h") == 1)) {
    std::cout << parser.help() << std::endl; exit(1); 
  };
  // Now we have everything and we just need to be careful with the 
  // sequence that we return so that it matches our 'argument_t' definitions
  return std::make_tuple( arguments["cell"].as<std::string>(), arguments["Grid"].as<int>(),
         arguments["grid"].as<std::vector<size_t>>(), arguments["electrons"].as<size_t>(),
         arguments["radius"].as<double>(), arguments["method"].as<std::string>(),
         arguments["decomposition"].as<bool>(), arguments["precision"].as<double>(),
         arguments["iterations"].as<int>(), arguments["structureFactor"].as<bool>(),
         arguments["terms"].as<std::vector<std::string>>(), arguments["shift"].as<std::vector<double>>(),
         arguments["Triples"].as<std::vector<bool>>(), arguments["CRtriples"].as<std::vector<bool>>() );
}; 

