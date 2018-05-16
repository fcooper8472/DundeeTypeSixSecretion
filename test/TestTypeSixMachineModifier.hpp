
#ifndef TESTTYPESIXMACHINEMODIFIER_HPP_
#define TESTTYPESIXMACHINEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "TypeSixMachineProperty.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "TypeSixMachineModifier.hpp"
#include "TypeSixSecretionEnumerations.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestTypeSixMachineModifier : public AbstractCellBasedTestSuite
{
public:

    void TestTypeSixMachineModifierException1() throw(Exception)
    {
        // Create a simple 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetBirthTime(-10.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_type);
            cells.push_back(p_cell);
            
	        std::vector<double>& attributes = mesh.GetNode(i)->rGetNodeAttributes();
	        attributes.resize(NA_VEC_LENGTH);
	        attributes[NA_THETA] = 1.23;
	        attributes[NA_LENGTH] = 2.34;
	        attributes[NA_RADIUS] = 3.45;
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulation(cell_population);
        simulation.SetOutputDirectory("TestTypeSixMachineModifierException1");
        simulation.SetEndTime(simulation.GetDt()/2.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulation.AddForce(p_force);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
        simulation.AddSimulationModifier(p_modifier);

        // Test the correct exception is thrown if the output directory has not been passed to the modifier
        TS_ASSERT_THROWS_THIS(simulation.Solve(),
            "SetOutputDirectory() must be called on TypeSixMachineModifier");
    }

    void TestTypeSixMachineModifierException2() throw(Exception)
    {
        // Create a simple 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetBirthTime(-10.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_type);
            cells.push_back(p_cell);
            
	        std::vector<double>& attributes = mesh.GetNode(i)->rGetNodeAttributes();
	        attributes.resize(NA_VEC_LENGTH);
	        attributes[NA_THETA] = 1.23;
	        attributes[NA_LENGTH] = 2.34;
	        attributes[NA_RADIUS] = 3.45;
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulation(cell_population);
        simulation.SetOutputDirectory("TestNodeBasedSimulationWithTypeSixMachineModifier");
        simulation.SetEndTime(simulation.GetDt()/2.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulation.AddForce(p_force);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
        simulation.AddSimulationModifier(p_modifier);        
        p_modifier->SetOutputDirectory("TestTypeSixMachineModifierException2");
        
        // Test the correct exception is thrown if any cell does not have a TypeSixMachineProperty
        TS_ASSERT_THROWS_THIS(simulation.Solve(),
            "TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
    }

    void TestNodeBasedSimulationWithTypeSixMachineModifier() throw(Exception)
    {
        // Create a simple 2D NodeBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetBirthTime(-10.0);

            MAKE_PTR(TypeSixMachineProperty, p_property);
            std::vector<double> machine_angles;
            machine_angles.push_back(0.0);
	        p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(i%4, machine_angles));

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_type);
	        p_cell->AddCellProperty(p_property);
            cells.push_back(p_cell);
            
	        std::vector<double>& attributes = mesh.GetNode(i)->rGetNodeAttributes();
	        attributes.resize(NA_VEC_LENGTH);
	        attributes[NA_THETA] = 0.0;
	        attributes[NA_LENGTH] = 0.1;
	        attributes[NA_RADIUS] = 0.5;
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create a simulation
        OffLatticeSimulation<2> simulation(cell_population);
        simulation.SetOutputDirectory("TestNodeBasedSimulationWithTypeSixMachineModifier");
        simulation.SetEndTime(10*simulation.GetDt());

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulation.AddForce(p_force);

        // Create a volume-tracking modifier and pass it to the simulation
        MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
        p_modifier->SetOutputDirectory("TestNodeBasedSimulationWithTypeSixMachineModifier");
        simulation.AddSimulationModifier(p_modifier);

        // Run simulation
        simulation.Solve();
        //TS_ASSERT_THROWS_NOTHING(simulation.Solve());

        ///\todo Test something
    }
///\todo test archiving and parameter output method
};

#endif /*TESTTYPESIXMACHINEMODIFIER_HPP_*/
