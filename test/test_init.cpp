// test_init.cpp
//
// Copyright (c) 2017 Ben Lindsay <benjlindsay@gmail.com>

#include "test_init.hpp"
#include "catch.hpp"

TEST_CASE("Bare bones input", "[input]") {
  const char* yaml_1 =
      "dim: 3\n"
      "lx:  [1,2,3]\n"
      "nx:  [10,10,10]\n"
      "rho_0: 10.0\n";
  YAML::Node input = YAML::Load(yaml_1);

  // Test sim_plan Input
  std::cout.setstate(std::ios_base::failbit);
  Sim_Plan* sim_plan = Sim_Plan_Factory::New_Sim_Plan(input);
  std::cout.clear();
  Single_Sim_Plan* single_sim_plan;
  Canonical_Sim* sim;

  REQUIRE(dynamic_cast<Single_Sim_Plan*>(sim_plan) != NULL);
  single_sim_plan = dynamic_cast<Single_Sim_Plan*>(sim_plan);
  REQUIRE(dynamic_cast<Canonical_Sim*>(single_sim_plan->sim) != NULL);
  sim = dynamic_cast<Canonical_Sim*>(single_sim_plan->sim);

  SECTION("dx calculation") {
    REQUIRE(sim->dx[0] == 0.1);
    REQUIRE(sim->dx[1] == 0.2);
    REQUIRE(sim->dx[2] == 0.3);
  }

  const char* yaml_2 =
      "dim: 3\n"
      "dx:  [.1,.2,.3]\n"
      "nx:  [10,10,10]\n"
      "rho_0: 10.0\n";
  input = YAML::Load(yaml_2);

  // Test sim_plan Input
  std::cout.setstate(std::ios_base::failbit);
  sim_plan = Sim_Plan_Factory::New_Sim_Plan(input);
  std::cout.clear();

  REQUIRE(dynamic_cast<Single_Sim_Plan*>(sim_plan) != NULL);
  single_sim_plan = dynamic_cast<Single_Sim_Plan*>(sim_plan);
  REQUIRE(dynamic_cast<Canonical_Sim*>(single_sim_plan->sim) != NULL);
  sim = dynamic_cast<Canonical_Sim*>(single_sim_plan->sim);

  SECTION("Lx calculation") {
    REQUIRE(sim->Lx[0] == 1);
    REQUIRE(sim->Lx[1] == 2);
    REQUIRE(sim->Lx[2] == 3);
  }

  const char* yaml_3 =
      "dim: 3\n"
      "dx:  [.1,.1,.1]\n"
      "lx:  [0.9,1.0,1.05]\n"
      "rho_0: 10.0\n";
  input = YAML::Load(yaml_3);

  // Test sim_plan Input
  std::cout.setstate(std::ios_base::failbit);
  sim_plan = Sim_Plan_Factory::New_Sim_Plan(input);
  std::cout.clear();

  REQUIRE(dynamic_cast<Single_Sim_Plan*>(sim_plan) != NULL);
  single_sim_plan = dynamic_cast<Single_Sim_Plan*>(sim_plan);
  REQUIRE(dynamic_cast<Canonical_Sim*>(single_sim_plan->sim) != NULL);
  sim = dynamic_cast<Canonical_Sim*>(single_sim_plan->sim);

  SECTION("Nx calculation") {
    REQUIRE(sim->Nx[0] == 9);
    REQUIRE(sim->Nx[1] == 11);
    REQUIRE(sim->Nx[2] == 11);
  }
}
