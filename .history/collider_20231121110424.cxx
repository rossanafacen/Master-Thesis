// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License


#include "collider.h"

#include <cmath>
#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"
#include "nucleus.h"
#include "common.h"

#include <iostream>

double PHI_A;
double THETA_A;

double PHI_B;
double THETA_B;

namespace trento {

namespace {

// Helper functions for Collider ctor.

// Create one nucleus from the configuration.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  const auto& nucleon_dmin = var_map["nucleon-min-dist"].as<double>();
  return Nucleus::create(species, nucleon_dmin);
}

// Determine the maximum impact parameter.  If the configuration contains a
// non-negative value for bmax, use it; otherwise, fall back to the minimum-bias
// default.
double determine_bmax(const VarMap& var_map,
    const Nucleus& A, const Nucleus& B, const NucleonCommon& nc) {
  auto bmax = var_map["b-max"].as<double>();
  if (bmax < 0.)
    bmax = A.radius() + B.radius() + nc.max_impact();
  return bmax;
}

// Determine the asymmetry parameter (Collider::asymmetry_) for a pair of
// nuclei.  It's just rA/(rA+rB), falling back to 1/2 if both radii are zero
// (i.e. for proton-proton).
double determine_asym(const Nucleus& A, const Nucleus& B) {
  double rA = A.radius();
  double rB = B.radius();
  double sum = rA + rB;
  if (sum < 0.1)
    return 0.5;
  else
    return rA/sum;
}

}  // unnamed namespace

// Lots of members to initialize...
// Several helper functions are defined above.
Collider::Collider(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_common_(var_map),
      nevents_(var_map["number-events"].as<int>()),
      with_ncoll_(var_map["ncoll"].as<bool>()),
      bmin_(var_map["b-min"].as<double>()),
      bmax_(determine_bmax(var_map, *nucleusA_, *nucleusB_, nucleon_common_)),
      centmin_(var_map["cent-min"].as<double>()),
      centmax_(var_map["cent-max"].as<double>()),
      asymmetry_(determine_asym(*nucleusA_, *nucleusB_)),
      event_(var_map),
      output_(var_map) {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

// See header for explanation.
Collider::~Collider() = default;

void Collider::run_events() {
  // The main event loop.
  for (int n = 0; n < nevents_; ++n) {
    // Sampling the impact parameter also implicitly prepares the nuclei for
    // event computation, i.e. by sampling nucleon positions and participants.
    //auto collision_attr = sample_collision();
    //double b = std::get<0>(collision_attr);
    double b = sample_collision();
    //int ncoll = std::get<1>(collision_attr);

    // Pass the prepared nuclei to the Event.  It computes the entropy profile
    // (thickness grid) and other event observables.
    event_.compute(*nucleusA_, *nucleusB_, nucleon_common_);
 //Check if even is in a given centrality class
    if ( event_.multiplicity() < centmax_ && centmin_ <= event_.multiplicity()){
    // Write event data.
    output_(n, b, event_);
    }
    else{
        --n;
    }
  }
}

//std::tuple<double, int> Collider::sample_collision() {
double Collider::sample_collision() {
  // Sample impact parameters until at least one nucleon-nucleon pair
  // participates.  The bool 'collision' keeps track -- it is effectively a
  // logical OR over all possible participant pairs.
  // Returns the sampled impact parameter b, and binary collision number ncoll.
  double b;
  //int ncoll = 0;
  bool collision = false;

  do {
    // Sample b from P(b)db = 2*pi*b.
    b = bmin_ + (bmax_ - bmin_) * std::sqrt(random::canonical<double>());

    // Offset each nucleus depending on the asymmetry parameter (see header).
    nucleusA_->sample_nucleons(asymmetry_ * b);
    THETA_A = THETA;
    PHI_A = PHI;
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);
    THETA_B = THETA;
    PHI_B = PHI;


    // Check each nucleon-nucleon pair.
    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
        auto new_collision = nucleon_common_.participate(A, B);
        if (with_ncoll_) {
          if (new_collision && (!collision) ) event_.clear_TAB();
          // WK: to calculate binary collision denstiy, each collision 
          // contribute independently its Tpp. Therefore, if one pair collide, 
          // it calls the event object to accumulate Tpp to the Ncoll density
          // Ncoll density = Sum Tpp		
          if (new_collision) event_.accumulate_TAB(A, B, nucleon_profile_);
          //if (new_collision) event_.compute_ncoll();
        }
        
        
        collision = new_collision || collision;
      }
    }
  } while (!collision);

  //return std::make_tuple(b, ncoll);
  return b;
}

}  // namespace trento
