/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTCAPSULEBASEDDIVISIONRULES_HPP_
#define TESTCAPSULEBASEDDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "FixedCentreBasedDivisionRule.hpp"
#include "SmartPointers.hpp"
#include "CapsuleBasedDivisionRule.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCapsuleBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestFixedCentreBasedDivisionRule()
    {
    	EXIT_IF_PARALLEL;

		/* We then create a couple of cells at the base of each germarium.
		 * (we put two cells in each crypt to set off delta-notch patterning) */
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u,  false,  0.0, 0.0));
		nodes.push_back(new Node<2>(1u,  false,  1.0, 0.0));
		nodes.push_back(new Node<2>(2u,  false,  2.0, 0.0));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		//Create cells
		std::vector<CellPtr> cells;
		auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
		CellsGenerator<NoCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, 3u, p_diff_type);

		// Create cell population
		NodeBasedCellPopulation<2> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);
        c_vector<double, 2> expected_parent_location;
        expected_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        c_vector<double, 2> expected_daughter_location;
        expected_daughter_location[0] = 1.2;
        expected_daughter_location[1] = 3.4;

        // Set the division rule for our population to be the random direction division rule
        typedef CapsuleBasedDivisionRule<2,2> CapsuleRule;
        MAKE_PTR_ARGS(CapsuleRule, p_division_rule_to_set, (expected_daughter_location));

        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[0], 1.2, 1e-6);
        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[1], 3.4, 1e-6);

        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        // Check that the division rule returns the correct pair of vectors
        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_location;
        parent_location = positions.first;
        TS_ASSERT_DELTA(parent_location[0], expected_parent_location[0], 1e-6);
        TS_ASSERT_DELTA(parent_location[1], expected_parent_location[1], 1e-6);

        c_vector<double, 2> daughter_location;
        daughter_location = positions.second;
        TS_ASSERT_DELTA(daughter_location[0], expected_daughter_location[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location[1], expected_daughter_location[1], 1e-6);
    }

    void TestFixedCentreBasedDivisionRule3d()
    {
    	EXIT_IF_PARALLEL;

		/* We then create a couple of cells at the base of each germarium.
		 * (we put two cells in each crypt to set off delta-notch patterning) */
		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<2>(0u,  false,  0.0, 0.0, 0.0));
		nodes.push_back(new Node<2>(1u,  false,  1.0, 0.0, 0.0));
		nodes.push_back(new Node<2>(2u,  false,  2.0, 0.0, 0.0));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

		//Create cells
		std::vector<CellPtr> cells;
		auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
		CellsGenerator<NoCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, 3u, p_diff_type);

		// Create cell population
		NodeBasedCellPopulation<3> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);
        c_vector<double, 3> expected_parent_location;
        expected_parent_location = cell_population.GetLocationOfCellCentre(p_cell0);

        c_vector<double, 3> expected_daughter_location;
        expected_daughter_location[0] = 1.2;
        expected_daughter_location[1] = 3.4;
        expected_daughter_location[2] = 0.0;

        // Set the division rule for our population to be the random direction division rule
        typedef CapsuleBasedDivisionRule<3,3> CapsuleRule;
        MAKE_PTR_ARGS(CapsuleRule, p_division_rule_to_set, (expected_daughter_location));

        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[0], 1.2, 1e-6);
        TS_ASSERT_DELTA(p_division_rule_to_set->rGetDaughterLocation()[1], 3.4, 1e-6);

        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        // Check that the division rule returns the correct pair of vectors
        std::pair<c_vector<double, 3>, c_vector<double, 2> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 3> parent_location;
        parent_location = positions.first;
        TS_ASSERT_DELTA(parent_location[0], expected_parent_location[0], 1e-6);
        TS_ASSERT_DELTA(parent_location[1], expected_parent_location[1], 1e-6);
        TS_ASSERT_DELTA(parent_location[2], expected_parent_location[2], 1e-6);

        c_vector<double, 3> daughter_location;
        daughter_location = positions.second;
        TS_ASSERT_DELTA(daughter_location[0], expected_daughter_location[0], 1e-6);
        TS_ASSERT_DELTA(daughter_location[1], expected_daughter_location[1], 1e-6);
        TS_ASSERT_DELTA(daughter_location[2], expected_daughter_location[2], 1e-6);
    }
};

#endif /*TESTCAPSULEBASEDDIVISIONRULES_HPP_*/
