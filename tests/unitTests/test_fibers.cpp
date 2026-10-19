// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "Fibers.h"
#include "Parameters.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

namespace {

mshType make_single_element_mesh()
{
  mshType mesh;
  mesh.gnEl = 1;
  mesh.gnNo = 1;
  mesh.nG = 1;
  mesh.eNoN = 1;
  mesh.name = "unit-test-mesh";
  mesh.N.resize(1, 1);
  mesh.N(0, 0) = 1.0;
  mesh.gIEN.resize(1, 1);
  mesh.gIEN(0, 0) = 0;
  return mesh;
}

} // namespace

TEST(Fibers, ReadsOrthonormalConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.6, 0.8, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;

  EXPECT_FALSE(fibers.has_fibers());
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.nfibers, 1);

  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 1);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 0.6);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.8);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);

  fibers.clear();
  EXPECT_FALSE(fibers.has_fibers());
}

TEST(Fibers, ReadsMultipleOrthonormalConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  ASSERT_EQ(fibers.nfibers, 2);
  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 3);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);
}

TEST(Fibers, RejectsNonUnitConstantDirection)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{3.0, 4.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), std::runtime_error);
}

TEST(Fibers, RejectsNonOrthogonalConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.6, 0.8, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), std::runtime_error);
}

TEST(Fibers, LoadsCellDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.nfibers, 1);

  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 1);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
}

TEST(Fibers, LoadsPointDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.nfibers, 1);

  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 1);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
}

TEST(Fibers, RejectsCellTagWithPointData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), std::runtime_error);
}

TEST(Fibers, RejectsPointTagWithCellData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), std::runtime_error);
}

TEST(Fibers, RejectsMixedInputModes)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), std::runtime_error);
}
