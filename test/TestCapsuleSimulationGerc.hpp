
#ifndef TESTCAPSULESIMULATIONGERC_HPP_
#define TESTCAPSULESIMULATIONGERC_HPP_

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
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellIdWriter.hpp"

// Header files included in this project
#include "TypeSixSecretionEnumerations.hpp"
#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "CapsuleForce.hpp"
#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "SquareBoundaryCondition.hpp"
#include "CapsuleBasedDivisionRule.hpp"
#include "TypeSixMachineModifier.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "MachineStateCountWriter.hpp"



// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsuleSimulationGerc : public AbstractCellBasedTestSuite
{
public:



	void xTestSingleCapsuleSimulationWithDivisionAndMachinesKillerGerc() throw (Exception)
                		  {
		EXIT_IF_PARALLEL;

		//auto p_rand_gen = RandomNumberGenerator::Instance();

		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 100.0);
		c_vector<double, 4> domain_size;
		domain_size[0] = -1000.0;
		domain_size[1] = 1000.0;
		domain_size[2] = -1000.0;
		domain_size[3] = 1000.0;
		mesh.SetInitialBoxCollection(domain_size, 10.0);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
		{
			mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
			mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_THETA] = 0.0;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
		}

		// Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.6);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);



			double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			std::vector<double> machine_angles;
             machine_angles.push_back(rand_angle);
             MAKE_PTR(TypeSixMachineProperty, p_property);
             p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
             p_cell->AddCellProperty(p_property);

			//double birth_time = -RandomNumberGenerator::Instance()->ranf();
			p_cell->SetBirthTime(-0.9);
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			cells.push_back(p_cell);


		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<MachineStateCountWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKillerGerc");
		double dt = 1.0/1200.0;
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(10);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);


		//
		MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
		simulator.AddCellKiller(p_killer);
		/*
		 * We now create a capsuleforce law and pass it to the simulation
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		p_capsule_force->SetYoungModulus(200.0);

		simulator.AddForce(p_capsule_force);
		//

		MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		p_modifier->SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKillerGerc");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(8.20527000050075); // was 1.0075
		simulator.Solve();
		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
                		  }

	void xTestSingleCapsuleSimulationWithDivisionAndMachinesNoKillerGerc() throw (Exception)
                		  {
		EXIT_IF_PARALLEL;

		//auto p_rand_gen = RandomNumberGenerator::Instance();

		// Create some capsules
		std::vector<Node<2>*> nodes;
		double initial_pos=20.0;
		nodes.push_back(new Node<2>(0u, Create_c_vector(initial_pos, initial_pos)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 100.0);
		c_vector<double, 4> domain_size;
		domain_size[0] = -1000.0;
		domain_size[1] = 1000.0;
		domain_size[2] = -1000.0;
		domain_size[3] = 1000.0;
		mesh.SetInitialBoxCollection(domain_size, 2.0*initial_pos);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
		{
			mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
			mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_THETA] = 0.0;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
		}

		// Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.6);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);



			std::vector<double> machine_angles;
			             machine_angles.push_back(0.0);
			             MAKE_PTR(TypeSixMachineProperty, p_property);
			             p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
			             p_cell->AddCellProperty(p_property);

			//double birth_time = -RandomNumberGenerator::Instance()->ranf();
			p_cell->SetBirthTime(-0.9);
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			cells.push_back(p_cell);


		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<MachineStateCountWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesNoKillerGerc");
		double dt = 1.0/1200.0;
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(1);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);


		//
		//MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
		//simulator.AddCellKiller(p_killer);
		/*
		 * We now create a capsuleforce law and pass it to the simulation
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		p_capsule_force->SetYoungModulus(200.0);

		simulator.AddForce(p_capsule_force);
		//

		MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		p_modifier->SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesNoKillerGerc");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(8.20527000050075); // was 1.0075
		simulator.Solve();
		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
    }



};

#endif /*TESTCAPSULESIMULATIONGerc_HPP_*/
