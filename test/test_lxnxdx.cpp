// test_lxnxdx.cpp
//
// Copyright (c) 2019 Ben Lindsay <benjlindsay@gmail.com>

#include "test_lxnxdx.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

class LxNxDxTest : public ::testing::Test {
 public:
  Canonical_Sim *sim;
  void SetUp(const std::string var_to_drop) {
    std::cout << var_to_drop << std::endl;
    const char *yaml_1 =
        "dim: 3\n"
        "lx:  [1.1,2.2,3.3]\n"
        "nx:  [11,11,11]\n"
        "dx:  [0.1,0.2,0.3]\n"
        "rho_0: 10.0\n";
    YAML::Node input = YAML::Load(yaml_1);
    input.remove(var_to_drop);
    Sim_Plan *sim_plan = Sim_Plan_Factory::New_Sim_Plan(input);
    Single_Sim_Plan *single_sim_plan =
        dynamic_cast<Single_Sim_Plan *>(sim_plan);
    sim = dynamic_cast<Canonical_Sim *>(single_sim_plan->sim);
  }
  void TearDown() { delete sim; }
};

TEST_F(LxNxDxTest, CalcuatesLxRight) {
  SetUp("lx");
  ASSERT_EQ(sim->Lx[0], 1.1);
  ASSERT_EQ(sim->Lx[1], 2.2);
  ASSERT_EQ(sim->Lx[2], 3.3);
}

TEST_F(LxNxDxTest, CalcuatesNxRight) {
  SetUp("nx");
  ASSERT_EQ(sim->Nx[0], 11);
  ASSERT_EQ(sim->Nx[1], 11);
  ASSERT_EQ(sim->Nx[2], 11);
}

TEST_F(LxNxDxTest, CalculatesDxRight) {
  SetUp("dx");
  ASSERT_EQ(sim->dx[0], 0.1);
  ASSERT_EQ(sim->dx[1], 0.2);
  ASSERT_EQ(sim->dx[2], 0.3);
}
