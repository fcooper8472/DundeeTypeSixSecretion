/*

Copyright (c) 2005-2015, University of Oxford.
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


#ifndef TESTCAPSULESIMULATION_HPP_
#define TESTCAPSULESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cycle/UniformCellCycleModel.hpp>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Cell.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"


// Header files included in this project
#include "TypeSixSecretionEnumerations.hpp"
#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "CapsuleForce.hpp"
#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "SquareBoundaryCondition.hpp"

// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsuleSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleCapsuleSimulation() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        const unsigned num_nodes = 20u;
        auto p_rand_gen = RandomNumberGenerator::Instance();

        // Create some capsules
        std::vector<Node<2>*> nodes;

        nodes.push_back(new Node<2>(0u, Create_c_vector(1.0, 0.0)));
        for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
        {
            c_vector<double, 2> safe_location;

            bool safe = false;
            while(!safe)
            {
                safe = true;
                safe_location = Create_c_vector(10.0 * p_rand_gen->ranf(), 10.0 * p_rand_gen->ranf());

                for(auto&& p_node : nodes)
                {
                    if(norm_2(p_node->rGetLocation() - safe_location) < 2.0)
                    {
                        safe = false;
                    }
                }
            }

            nodes.push_back(new Node<2>(node_idx, safe_location));
        }

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        mesh.GetNode(0u)->AddNodeAttribute(0.0);
        mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 1.0;

        for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
        {
            mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
            mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = 2.0 * M_PI * p_rand_gen->ranf();
            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = p_rand_gen->NormalRandomDeviate(2.0, 0.5);
            mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
        }

        //Create cells
        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        // Create cell population
        NodeBasedCellPopulation<2> population(mesh, cells);

        population.AddCellWriter<CapsuleOrientationWriter>();
        population.AddCellWriter<CapsuleScalingWriter>();

        // Create simulation
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory("CapsuleSimulation");
        simulator.SetDt(1.0/1200.0);
        simulator.SetSamplingTimestepMultiple(1u);

        auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
        simulator.SetNumericalMethod(p_numerical_method);

        /*
         * We now create a force law and pass it to the simulation
         * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
         */
        auto p_calsule_force = boost::make_shared<CapsuleForce<2>>();
        simulator.AddForce(p_calsule_force);

//        auto p_boundary_condition = boost::make_shared<SquareBoundaryCondition>(&population);
//        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* We then set an end time and run the simulation */
        simulator.SetEndTime(1.0);
        simulator.Solve();
    }
};

#endif /*TESTCAPSULESIMULATION_HPP_*/