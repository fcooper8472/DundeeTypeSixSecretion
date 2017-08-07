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


#ifndef TESTRODCELLSIMULATION_HPP_
#define TESTRODCELLSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cycle/UniformCellCycleModel.hpp>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellAncestorWriter.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedCentreBasedDivisionRule.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

// Header files included in this project
#include "DifferentialAdhesionGermariumForce.hpp"
#include "DrosophilaOogenesisEnumerations.hpp"
#include "DrosophilaOogenesisSimulationModifier.hpp"
#include "GermariumBoundaryCondition.hpp"
#include "GermariumDivisionRule.hpp"

// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestOogenesis : public AbstractCellBasedTestSuite
{
public:

    void Test3dOogenesis() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        /* First we set up some numbers that will define the crypt and villus geometry */
        double germarium_radius = 1.0;

        /* We then create a couple of cells at the base of each germarium.
         * (we put two cells in each crypt to set off delta-notch patterning) */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.0, 0.0, 0.0));

        nodes.push_back(new Node<3>(1u,  false,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, -1.0, 0.0));

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        mesh.GetNode(0u)->SetRadius(0.7);
        mesh.GetNode(1u)->SetRadius(0.3);
        mesh.GetNode(2u)->SetRadius(0.3);

        /*
         * Next we have to create the cells that will be associated with these nodes.
         * So we make an empty vector in which to store the cells and then loop over
         * each node, adding cells as we go.
         */
        std::vector<CellPtr> cells;

        // Set up the proliferative types
        auto p_stem_type = boost::make_shared<StemCellProliferativeType>();

        // Set up the mutation state
        auto p_state = boost::make_shared<WildTypeCellMutationState>();

        // Set up cell labels
        auto p_aggregate_label = boost::make_shared<CellLabel>(TYPE_AGGREGATE);
        auto p_follicle_label = boost::make_shared<CellLabel>(TYPE_FOLLICLE);

        // The "aggregate stem cell"
        {
            /* We then create a cell with a mutation state (Wild Type in this case) and a cell cycle model */
            cells.push_back(boost::make_shared<Cell>(p_state, new FixedG1GenerationalCellCycleModel()));

            auto p_agg_ccm = static_cast<FixedG1GenerationalCellCycleModel*>(cells.back()->GetCellCycleModel());
            p_agg_ccm->SetDimension(3);
            p_agg_ccm->SetMaxTransitGenerations(1u);
            p_agg_ccm->SetStemCellG1Duration(30.0);
            p_agg_ccm->SetMDuration(0.0);
            p_agg_ccm->SetG2Duration(0.0);
            p_agg_ccm->SetSDuration(0.0);

            // Set the params we want for the cell we have just created
            cells.back()->SetBirthTime(0.0);
            cells.back()->SetCellProliferativeType(p_stem_type);
            cells.back()->AddCellProperty(p_aggregate_label);
        }

        // The two follicle stem cells
        for(auto&& i : {1, 2})
        {
            /* We then create a cell with a mutation state (Wild Type in this case) and a cell cycle model */
            cells.push_back(boost::make_shared<Cell>(p_state, new FixedG1GenerationalCellCycleModel()));

            auto p_fol_ccm = static_cast<FixedG1GenerationalCellCycleModel*>(cells.back()->GetCellCycleModel());
            p_fol_ccm->SetDimension(3);
            p_fol_ccm->SetMaxTransitGenerations(1u);
            p_fol_ccm->SetStemCellG1Duration(2.0);
            p_fol_ccm->SetMDuration(0.0);
            p_fol_ccm->SetG2Duration(0.0);
            p_fol_ccm->SetSDuration(0.0);

            // Set the params we want for the cell we have just created
            cells.back()->SetBirthTime(0.0);
            cells.back()->SetCellProliferativeType(p_stem_type);
            cells.back()->AddCellProperty(p_follicle_label);
        }

        /*
         * We now create a cell population, which keeps track of a mesh and cells and the association between them.
         * In this case we need a `NodeBasedCellPopulation` in three dimensions.
         */
        NodeBasedCellPopulation<3> germarium(mesh, cells);

        c_vector<double, 3> aggregate_sep = 0.5 * unit_vector<double>(3, 0);
        c_vector<double, 3> follicle_sep = 0.2 * unit_vector<double>(3, 0);

        auto p_division_rule = boost::make_shared<GermariumDivisionRule<3, 3>>(aggregate_sep, follicle_sep);
        germarium.SetCentreBasedDivisionRule(p_division_rule);

        /* We limit the absolute movement that cells can make to cause error messages if numerics become unstable */
        germarium.SetAbsoluteMovementThreshold(10);

        /* We then instruct the cell population to output some useful information for plotting in VTK format in e.g. paraview */
        germarium.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        germarium.AddCellWriter<CellRadiusWriter>();

        /*
         * We now set up our cell-based simulation class.
         * We have sub-classed the main simulator class `OffLatticeSimulation`
         * which is in the core of Chaste with a new simulation class `SimplifiedDeltaNotchOffLatticeSimulation`
         * which can be found in this project's `src` folder.
         *
         * The reason for this is that this class performs calculations of average Delta levels surrounding each cell
         * at the end of each timestep. You will see a minimum number of methods have been overridden, and the class
         * is fairly simple.
         */
        OffLatticeSimulation<3> simulator(germarium);
        simulator.SetOutputDirectory("DrosophilaOogenesis");
        simulator.SetDt(1.0/120.0);
        /* We limit the output to every 120 time steps (1 hour) to reduce output file sizes */
        simulator.SetSamplingTimestepMultiple(50);

        // Add simulation modifier
        auto p_mod = boost::make_shared<DrosophilaOogenesisSimulationModifier<3>>();
        simulator.AddSimulationModifier(p_mod);

        /*
         * We now create a force law and pass it to the simulation
         * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
         */
        auto p_linear_force = boost::make_shared<DifferentialAdhesionGermariumForce<3>>();
        p_linear_force->SetMeinekeSpringStiffness(30.0); // default is 15.0;
        p_linear_force->SetMeinekeDivisionRestingSpringLength(0.5);
        p_linear_force->SetCutOffLength(2.0);
        p_linear_force->SetFollicleAggregateSpringConstantMultiplier(0.00001);
        simulator.AddForce(p_linear_force);

        /*
         * The most complex part of this problem definition is that of the boundary condition that limits
         * cell locations to a 2D surface in 3D space. This has been defined in a separate class
         * `MultipleCryptGeometryBoundaryCondition` which can be found in this project's `src` folder.
         */
        auto p_boundary_condition = boost::make_shared<GermariumBoundaryCondition>(&germarium, 0.0);
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* We then set an end time and run the simulation */
        simulator.SetEndTime(100.0);
        simulator.Solve(); // to 250 hours
    }
};

#endif /*TESTRODCELLSIMULATION_HPP_*/