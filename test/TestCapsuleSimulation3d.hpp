
#ifndef TESTCAPSULESIMULATION3D_HPP_
#define TESTCAPSULESIMULATION3D_HPP_

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
#include "DiffusionForce.hpp"


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

class TestCapsuleSimulation3d : public AbstractCellBasedTestSuite
{
public:



	void xTestSmallSymmetric3dCapsuleSimulation() throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<3>*> nodes;

		nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 0.0)));
		nodes.push_back(new Node<3>(1u, Create_c_vector(4.0, 6.0, 0.0)));
		nodes.push_back(new Node<3>(2u, Create_c_vector(7.0, 4.0, 0.0)));
		nodes.push_back(new Node<3>(3u, Create_c_vector(7.0, 6.0, 0.0)));
		nodes.push_back(new Node<3>(4u, Create_c_vector(7.0, 5.0, 2.0)));
		nodes.push_back(new Node<3>(5u, Create_c_vector(4.0, 5.0, -2.0)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(3u)->AddNodeAttribute(0.0);
		mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(4u)->AddNodeAttribute(0.0);
		mesh.GetNode(4u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_PHI] = 0.25 * M_PI;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(5u)->AddNodeAttribute(0.0);
		mesh.GetNode(5u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_PHI] = 0.75 * M_PI;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		//Create cells
		std::vector<CellPtr> cells;
		auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
		CellsGenerator<NoCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

		// Create cell population
		NodeBasedCellPopulation<3> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();

		// Create simulation
		OffLatticeSimulation<3> simulator(population);
		simulator.SetOutputDirectory("TestSmallSymmetric3dCapsuleSimulation");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(1u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
		simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
		simulator.AddForce(p_capsule_force);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestSmallSymmetric3dCapsuleSimulationWithMachineProperties() throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<3>*> nodes;

		nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 0.0)));
		nodes.push_back(new Node<3>(1u, Create_c_vector(4.0, 6.0, 0.0)));
		nodes.push_back(new Node<3>(2u, Create_c_vector(7.0, 4.0, 0.0)));
		nodes.push_back(new Node<3>(3u, Create_c_vector(7.0, 6.0, 0.0)));
		nodes.push_back(new Node<3>(4u, Create_c_vector(7.0, 5.0, 2.0)));
		nodes.push_back(new Node<3>(5u, Create_c_vector(4.0, 5.0, -2.0)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(3u)->AddNodeAttribute(0.0);
		mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(4u)->AddNodeAttribute(0.0);
		mesh.GetNode(4u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_PHI] = 0.25 * M_PI;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(5u)->AddNodeAttribute(0.0);
		mesh.GetNode(5u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_PHI] = 0.75 * M_PI;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(5u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		//Create cells
		 // Create cells
				std::vector<CellPtr> cells;
				MAKE_PTR(WildTypeCellMutationState, p_state);
				MAKE_PTR(DifferentiatedCellProliferativeType, p_type);
				for (unsigned i=0; i<mesh.GetNumNodes(); i++)
				{
					UniformCellCycleModel* p_model = new UniformCellCycleModel();
					p_model->SetMinCellCycleDuration(1.0);
					p_model->SetMaxCellCycleDuration(1.01);
					CellPtr p_cell(new Cell(p_state, p_model));
					p_cell->SetCellProliferativeType(p_type);

					double angle_1 = 0.0;
					double angle_2 = M_PI;

					std::vector<double> machine_angles1;
					machine_angles1.push_back(angle_1);
					machine_angles1.push_back(angle_2);



					MAKE_PTR(TypeSixMachineProperty, p_property);
					p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_angles1));


					p_cell->AddCellProperty(p_property);
					p_cell->SetBirthTime(-0.9);
					mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

					cells.push_back(p_cell);
				}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
		population.SetCentreBasedDivisionRule(p_division_rule);


		// Create simulation
		OffLatticeSimulation<3> simulator(population);
		simulator.SetOutputDirectory("TestSmallSymmetric3dCapsuleSimulationWithMachineProperties");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(1u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
		simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
		simulator.AddForce(p_capsule_force);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestSmallSymmetric3dCapsuleSimulationWithMachinePropertiesAndDivision() throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<3>*> nodes;

		nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 2.0)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;




		//Create cells
		 // Create cells
				std::vector<CellPtr> cells;
				MAKE_PTR(WildTypeCellMutationState, p_state);
				MAKE_PTR(TransitCellProliferativeType, p_type);
				for (unsigned i=0; i<mesh.GetNumNodes(); i++)
				{
					UniformCellCycleModel* p_model = new UniformCellCycleModel();
					p_model->SetMinCellCycleDuration(1.0);
					p_model->SetMaxCellCycleDuration(1.01);
					CellPtr p_cell(new Cell(p_state, p_model));
					p_cell->SetCellProliferativeType(p_type);

					double angle_1 = 0.0;
					double angle_2 = M_PI;

					std::vector<double> machine_angles1;
					machine_angles1.push_back(angle_1);
					machine_angles1.push_back(angle_2);



					MAKE_PTR(TypeSixMachineProperty, p_property);
					p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_angles1));


					p_cell->AddCellProperty(p_property);
					p_cell->SetBirthTime(-0.9);
					mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

					cells.push_back(p_cell);
				}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
		population.SetCentreBasedDivisionRule(p_division_rule);


		// Create simulation
		OffLatticeSimulation<3> simulator(population);
		simulator.SetOutputDirectory("TestSmallSymmetric3dCapsuleSimulationWithMachinePropertiesAndDivision");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(100u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
		simulator.SetNumericalMethod(p_numerical_method);


		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
		simulator.AddForce(p_capsule_force);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(0.7);
		simulator.Solve();

		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2u);
	}

	void TestSmallSymmetric3dCapsuleSimulationWithMachinePropertiesAndDivisionAndModifier() throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<3>*> nodes;

		nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 0.0)));
		nodes.push_back(new Node<3>(1u, Create_c_vector(4.0, 6.0,-1.0)));
		nodes.push_back(new Node<3>(2u, Create_c_vector(7.0, 4.0, 2.0)));
		nodes.push_back(new Node<3>(3u, Create_c_vector(7.0, 6.0, -1.0)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = -0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_PHI] = -0.25 * M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = -0.25 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_PHI] = -0.25 * M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(3u)->AddNodeAttribute(0.0);
		mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_PHI] = -0.25 * M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;




//Create cells
 // Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.5);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);

			double angle_1 = 0.5;
			double angle_2 = M_PI;

			std::vector<double> machine_angles1;
			machine_angles1.push_back(angle_1);
			machine_angles1.push_back(angle_2);



			MAKE_PTR(TypeSixMachineProperty, p_property);
			p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_angles1));


			double birth_time = -RandomNumberGenerator::Instance()->ranf();
			p_cell->AddCellProperty(p_property);
			p_cell->SetBirthTime(birth_time);
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			cells.push_back(p_cell);
		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
				population.SetCentreBasedDivisionRule(p_division_rule);


		// Create simulation
		OffLatticeSimulation<3> simulator(population);
		simulator.SetOutputDirectory("Test3dCapsuleWithMachineDivisionAndModifier");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(1u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
		simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
		simulator.AddForce(p_capsule_force);

		MAKE_PTR(TypeSixMachineModifier<3>, p_modifier);
		p_modifier->SetOutputDirectory("Test3dCapsuleWithMachineDivisionAndModifier");
		p_modifier->Setk_1(0.0);
		simulator.AddSimulationModifier(p_modifier);


		unsigned num_machines=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());

		TS_ASSERT_EQUALS(num_machines,2u);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(2.0); // was 1.0075
		simulator.SetSamplingTimestepMultiple(20);

		simulator.Solve();

		//		unsigned num_machines2=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());
		//
		//		TS_ASSERT_EQUALS(num_machines2,2u);
		//
		//		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());

		/* We then set an end time and run the simulation */


		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2u);
	}

    void TestMachinesWithModifiersAndDivision3d()
    {
        EXIT_IF_PARALLEL;
        // Create some capsules
        std::vector<Node<3>*> nodes;

        nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 2.0)));

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 150.5);

        mesh.GetNode(0u)->AddNodeAttribute(0.0);
        mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;

        mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


        //Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetMinCellCycleDuration(1.0);
            p_model->SetMaxCellCycleDuration(1.01);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_type);


            p_cell->SetBirthTime(-0.6);
            mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


            double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
            double azimuthal_coordinate = M_PI ;


            std::vector<double> machine_coordinates;
            machine_coordinates.push_back(vertical_coordinate);
            machine_coordinates.push_back(azimuthal_coordinate);

            MAKE_PTR(TypeSixMachineProperty, p_property);
            p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1, machine_coordinates));

            p_cell->AddCellProperty(p_property);

            cells.push_back(p_cell);
        }

        // Create cell population
        NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CapsuleOrientationWriter>();
        population.AddCellWriter<CapsuleScalingWriter>();

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
                 population.SetCentreBasedDivisionRule(p_division_rule);

        // Create simulation
        OffLatticeSimulation<3> simulator(population);
        simulator.SetOutputDirectory("TestMachineModifierWithDivision3d");
        simulator.SetDt(1.0/1200.0);
        simulator.SetSamplingTimestepMultiple(30u);

        auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
        simulator.SetNumericalMethod(p_numerical_method);

        /*
         * We now create a force law and pass it to the simulation
         * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
         */
        auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
        simulator.AddForce(p_capsule_force);

        MAKE_PTR(TypeSixMachineModifier<3>, p_modifier);
        p_modifier->SetOutputDirectory("TestMachineModifierWithDivision3d");
        p_modifier->Setk_1(0.0);
        p_modifier->Setk_5(0.0);
        p_modifier->Setk_2(0.0);


        simulator.AddSimulationModifier(p_modifier);

    	unsigned num_machines=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());

    	TS_ASSERT_EQUALS(num_machines,1u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),1u);


        /* We then set an end time and run the simulation */
        simulator.SetEndTime(0.35);
        simulator.Solve();
    	unsigned num_machines2=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());

    	TS_ASSERT_EQUALS(num_machines2,1u);

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),1u);

        simulator.SetEndTime(0.5);
        simulator.Solve();
        unsigned num_machines3=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());
        TS_ASSERT_EQUALS(num_machines3,1u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2u);
    }

    void TestMachinesWithModifiersAndDivisionVis3d()
       {
           EXIT_IF_PARALLEL;
           // Create some capsules
           std::vector<Node<3>*> nodes;

           nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 1.0)));

           /*
            * We then convert this list of nodes to a `NodesOnlyMesh`,
            * which doesn't do very much apart from keep track of the nodes.
            */
           NodesOnlyMesh<3> mesh;
           mesh.ConstructNodesWithoutMesh(nodes, 150.5);

           mesh.GetNode(0u)->AddNodeAttribute(0.0);
           mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
           mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
   		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;

           mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
           mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


           //Create cells
           std::vector<CellPtr> cells;
           MAKE_PTR(WildTypeCellMutationState, p_state);
           MAKE_PTR(TransitCellProliferativeType, p_type);
           for (unsigned i=0; i<mesh.GetNumNodes(); i++)
           {
               UniformCellCycleModel* p_model = new UniformCellCycleModel();
               p_model->SetMinCellCycleDuration(1.0);
               p_model->SetMaxCellCycleDuration(1.01);
               CellPtr p_cell(new Cell(p_state, p_model));
               p_cell->SetCellProliferativeType(p_type);


               p_cell->SetBirthTime(-0.6);
               mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


               double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
               double azimuthal_coordinate = M_PI ;


               std::vector<double> machine_coordinates;
               machine_coordinates.push_back(vertical_coordinate);
               machine_coordinates.push_back(azimuthal_coordinate);

               MAKE_PTR(TypeSixMachineProperty, p_property);
               p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1, machine_coordinates));

               p_cell->AddCellProperty(p_property);

               cells.push_back(p_cell);
           }

           // Create cell population
           NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

           population.AddCellWriter<CellIdWriter>();
           population.AddCellWriter<CapsuleOrientationWriter>();
           population.AddCellWriter<CapsuleScalingWriter>();

           boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
                    population.SetCentreBasedDivisionRule(p_division_rule);

           // Create simulation
           OffLatticeSimulation<3> simulator(population);
           simulator.SetOutputDirectory("TestMachineModifierWithDivision3dVis");
           simulator.SetDt(1.0/1200.0);
           simulator.SetSamplingTimestepMultiple(30u);

           auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
           simulator.SetNumericalMethod(p_numerical_method);

           /*
            * We now create a force law and pass it to the simulation
            * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
            */
           auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
           simulator.AddForce(p_capsule_force);

           MAKE_PTR(TypeSixMachineModifier<3>, p_modifier);
           p_modifier->SetOutputDirectory("TestMachineModifierWithDivision3dVis");
           p_modifier->Setk_1(0.0);
           p_modifier->Setk_5(0.0);
           p_modifier->Setk_2(0.0);


           simulator.AddSimulationModifier(p_modifier);

       	unsigned num_machines=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());

       	TS_ASSERT_EQUALS(num_machines,1u);
           TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),1u);



           simulator.SetEndTime(2.5);
           simulator.Solve();
           unsigned num_machines3=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());
           TS_ASSERT_EQUALS(num_machines3,1u);
           TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2u);
       }

    void TestMachinesWithModifiersAndDivisionVisManyMachines3d()
          {
              EXIT_IF_PARALLEL;
              // Create some capsules
              std::vector<Node<3>*> nodes;

              nodes.push_back(new Node<3>(0u, Create_c_vector(4.0, 4.0, 1.0)));

              /*
               * We then convert this list of nodes to a `NodesOnlyMesh`,
               * which doesn't do very much apart from keep track of the nodes.
               */
              NodesOnlyMesh<3> mesh;
              mesh.ConstructNodesWithoutMesh(nodes, 150.5);

              mesh.GetNode(0u)->AddNodeAttribute(0.0);
              mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
              mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25 * M_PI;
      		mesh.GetNode(0u)->rGetNodeAttributes()[NA_PHI] = 0.5 * M_PI;

              mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
              mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


              //Create cells
              std::vector<CellPtr> cells;
              MAKE_PTR(WildTypeCellMutationState, p_state);
              MAKE_PTR(TransitCellProliferativeType, p_type);
              for (unsigned i=0; i<mesh.GetNumNodes(); i++)
              {
                  UniformCellCycleModel* p_model = new UniformCellCycleModel();
                  p_model->SetMinCellCycleDuration(1.0);
                  p_model->SetMaxCellCycleDuration(1.5);
                  CellPtr p_cell(new Cell(p_state, p_model));
                  p_cell->SetCellProliferativeType(p_type);

                  double birth_time = -1.0*RandomNumberGenerator::Instance()->ranf();
                  p_cell->SetBirthTime(birth_time);
                  mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


                  double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
                  double azimuthal_coordinate = M_PI ;


                  std::vector<double> machine_coordinates;
                  machine_coordinates.push_back(vertical_coordinate);
                  machine_coordinates.push_back(azimuthal_coordinate);

                  MAKE_PTR(TypeSixMachineProperty, p_property);
                  p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1, machine_coordinates));

                  p_cell->AddCellProperty(p_property);

                  cells.push_back(p_cell);
              }

              // Create cell population
              NodeBasedCellPopulationWithCapsules<3> population(mesh, cells);

              population.AddCellWriter<CellIdWriter>();
              population.AddCellWriter<CapsuleOrientationWriter>();
              population.AddCellWriter<CapsuleScalingWriter>();

              boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule(new CapsuleBasedDivisionRule<3,3>());
                       population.SetCentreBasedDivisionRule(p_division_rule);

              // Create simulation
              OffLatticeSimulation<3> simulator(population);
              simulator.SetOutputDirectory("TestMachineModifierWithDivision3dManyVis");
              simulator.SetDt(1.0/1200.0);
              simulator.SetSamplingTimestepMultiple(30u);

              auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<3,3>>();
              simulator.SetNumericalMethod(p_numerical_method);

              /*
               * We now create a force law and pass it to the simulation
               * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
               */
              auto p_capsule_force = boost::make_shared<CapsuleForce<3>>();
              simulator.AddForce(p_capsule_force);

              MAKE_PTR(DiffusionForce<3>, p_force_diff);
                        p_force_diff->SetAbsoluteTemperature(0.150);
                                   simulator.AddForce(p_force_diff);


              MAKE_PTR(TypeSixMachineModifier<3>, p_modifier);
              p_modifier->SetOutputDirectory("TestMachineModifierWithDivision3dManyVis");
              //p_modifier->Setk_1(0.0);
              //p_modifier->Setk_5(0.0);
              //p_modifier->Setk_2(0.0);


              simulator.AddSimulationModifier(p_modifier);

          	unsigned num_machines=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());

          	TS_ASSERT_EQUALS(num_machines,1u);
              TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),1u);



              simulator.SetEndTime(10.5);
              simulator.Solve();
              unsigned num_machines3=p_modifier->GetTotalNumberOfMachines(simulator.rGetCellPopulation());
              TS_ASSERT_EQUALS(num_machines3,1u);
              TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2u);
          }





};

#endif /*TESTCAPSULESIMULATION3D_HPP_*/
