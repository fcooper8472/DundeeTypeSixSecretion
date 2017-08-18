
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

// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsuleSimulation : public AbstractCellBasedTestSuite
{
public:

    void xTestSmallSymmetricCapsuleSimulation() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        // Create some capsules
        std::vector<Node<2>*> nodes;

        nodes.push_back(new Node<2>(0u, Create_c_vector(4.0, 4.0)));
        nodes.push_back(new Node<2>(1u, Create_c_vector(4.0, 6.0)));
        nodes.push_back(new Node<2>(2u, Create_c_vector(7.0, 4.0)));
        nodes.push_back(new Node<2>(3u, Create_c_vector(7.0, 6.0)));

        /*
         * We then convert this list of nodes to a `NodesOnlyMesh`,
         * which doesn't do very much apart from keep track of the nodes.
         */
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 150.5);

        mesh.GetNode(0u)->AddNodeAttribute(0.0);
        mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.25 * M_PI;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

        mesh.GetNode(1u)->AddNodeAttribute(0.0);
        mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_ANGLE] = -0.25 * M_PI;
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

        mesh.GetNode(2u)->AddNodeAttribute(0.0);
        mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(2u)->rGetNodeAttributes()[NA_ANGLE] = -0.25 * M_PI;
        mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

        mesh.GetNode(3u)->AddNodeAttribute(0.0);
        mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(3u)->rGetNodeAttributes()[NA_ANGLE] = 0.25 * M_PI;
        mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

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
        simulator.SetOutputDirectory("TestSmallSymmetricCapsuleSimulation");
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

        /* We then set an end time and run the simulation */
        simulator.SetEndTime(100.0/1200.0);
        simulator.Solve();
    }




    void NoTestLongerCapsuleSimulation() throw (Exception)
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
        mesh.ConstructNodesWithoutMesh(nodes, 10.0);

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
        simulator.SetOutputDirectory("TestLongerCapsuleSimulation");
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
        simulator.SetEndTime(150.0/1200.0);
        simulator.Solve();
    }

    void NoTestSingleCapsuleSimulationWithDivision() throw (Exception)
       {
   	    EXIT_IF_PARALLEL;
   	    std::cout<< "Ahhhhh 1" << "\n";

              const unsigned num_nodes = 1u;
              //auto p_rand_gen = RandomNumberGenerator::Instance();

              // Create some capsules
              std::vector<Node<2>*> nodes;

              nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
              for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
              {
                  c_vector<double, 2> safe_location;

                  bool safe = false;
                  while(!safe)
                  {
                      safe = true;
                      safe_location = Create_c_vector(3.0, 3.0);

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
              mesh.ConstructNodesWithoutMesh(nodes, 100.0);
              c_vector<double, 4> domain_size;
              domain_size[0] = -1000.0;
              domain_size[1] = 1000.0;
              domain_size[2] = -1000.0;
              domain_size[3] = 1000.0;
              mesh.SetInitialBoxCollection(domain_size, 10.0);

              mesh.GetNode(0u)->AddNodeAttribute(0.0);
              mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
              mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
              mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
              mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

              for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
              {
                  mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
                  mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                  mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
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



                  //MAKE_PTR(TypeSixMachineProperty, p_property);
   				//p_property->rGetMachineData().emplace_back(std::pair<unsigned, double>(i%4, 0.0));
   				//	p_cell->AddCellProperty(p_property);

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

              boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
              population.SetCentreBasedDivisionRule(p_division_rule);

              // Create simulation
              OffLatticeSimulation<2> simulator(population);
              simulator.SetOutputDirectory("TestSingleCapsuleWithDivision");
              double dt = 1.0/1200.0;
              simulator.SetDt(dt);
              simulator.SetSamplingTimestepMultiple(10);

              auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
              simulator.SetNumericalMethod(p_numerical_method);


   //
   //           MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
   //           simulator.AddCellKiller(p_killer);
              /*
               * We now create a force law and pass it to the simulation
               * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
               */
              auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
              simulator.AddForce(p_capsule_force);
   //
              auto p_boundary_condition = boost::make_shared<SquareBoundaryCondition>(&population);
              simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

              /* We then set an end time and run the simulation */
              simulator.SetEndTime(5.527000050075); // was 1.0075
              simulator.Solve();
              PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
          }

    void TestSingleCapsuleSimulationWithDivisionAndMachines() throw (Exception)
           {
       	    EXIT_IF_PARALLEL;

                  const unsigned num_nodes = 1u;
                  //auto p_rand_gen = RandomNumberGenerator::Instance();

                  // Create some capsules
                  std::vector<Node<2>*> nodes;

                  nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
                  for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
                  {
                      c_vector<double, 2> safe_location;

                      bool safe = false;
                      while(!safe)
                      {
                          safe = true;
                          safe_location = Create_c_vector(3.0, 3.0);

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
                  mesh.ConstructNodesWithoutMesh(nodes, 100.0);
                  c_vector<double, 4> domain_size;
                  domain_size[0] = -1000.0;
                  domain_size[1] = 1000.0;
                  domain_size[2] = -1000.0;
                  domain_size[3] = 1000.0;
                  mesh.SetInitialBoxCollection(domain_size, 10.0);

                  mesh.GetNode(0u)->AddNodeAttribute(0.0);
                  mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                  mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
                  mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
                  mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

                  for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
                  {
                      mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
                      mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                      mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
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



                      double rand_angle = 2*3.147*RandomNumberGenerator::Instance()->ranf();
                      MAKE_PTR(TypeSixMachineProperty, p_property);
       				 p_property->rGetMachineData().emplace_back(std::pair<unsigned, double>(i%4, rand_angle));
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

                  boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
                  population.SetCentreBasedDivisionRule(p_division_rule);

                  // Create simulation
                  OffLatticeSimulation<2> simulator(population);
                  simulator.SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachines");
                  double dt = 1.0/1200.0;
                  simulator.SetDt(dt);
                  simulator.SetSamplingTimestepMultiple(10);

                  auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
                  simulator.SetNumericalMethod(p_numerical_method);


       //
       //           MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
       //           simulator.AddCellKiller(p_killer);
                  /*
                   * We now create a force law and pass it to the simulation
                   * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
                   */
                  auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
                  simulator.AddForce(p_capsule_force);
       //
                  auto p_boundary_condition = boost::make_shared<SquareBoundaryCondition>(&population);
                  simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

                  MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
                  p_modifier->SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachines");
                  simulator.AddSimulationModifier(p_modifier);

                  /* We then set an end time and run the simulation */
                  simulator.SetEndTime(4.0527000050075); // was 1.0075
                  simulator.Solve();
                  PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
              }


    void TestSingleCapsuleSimulationWithDivisionAndMachinesKiller() throw (Exception)
               {
           	    EXIT_IF_PARALLEL;

                      const unsigned num_nodes = 1u;
                      //auto p_rand_gen = RandomNumberGenerator::Instance();

                      // Create some capsules
                      std::vector<Node<2>*> nodes;

                      nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
                      for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
                      {
                          c_vector<double, 2> safe_location;

                          bool safe = false;
                          while(!safe)
                          {
                              safe = true;
                              safe_location = Create_c_vector(3.0, 3.0);

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
                      mesh.ConstructNodesWithoutMesh(nodes, 100.0);
                      c_vector<double, 4> domain_size;
                      domain_size[0] = -1000.0;
                      domain_size[1] = 1000.0;
                      domain_size[2] = -1000.0;
                      domain_size[3] = 1000.0;
                      mesh.SetInitialBoxCollection(domain_size, 10.0);

                      mesh.GetNode(0u)->AddNodeAttribute(0.0);
                      mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                      mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
                      mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
                      mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

                      for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
                      {
                          mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
                          mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
                          mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
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



                          double rand_angle = 2*3.147*RandomNumberGenerator::Instance()->ranf();
                          MAKE_PTR(TypeSixMachineProperty, p_property);
           				 p_property->rGetMachineData().emplace_back(std::pair<unsigned, double>(i%4, rand_angle));
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

                      boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
                      population.SetCentreBasedDivisionRule(p_division_rule);

                      // Create simulation
                      OffLatticeSimulation<2> simulator(population);
                      simulator.SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKiller");
                      double dt = 1.0/1200.0;
                      simulator.SetDt(dt);
                      simulator.SetSamplingTimestepMultiple(10);

                      auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
                      simulator.SetNumericalMethod(p_numerical_method);


           //
                      MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
                      simulator.AddCellKiller(p_killer);
                      /*
                       * We now create a force law and pass it to the simulation
                       * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
                       */
                      auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
                      simulator.AddForce(p_capsule_force);
           //
                      auto p_boundary_condition = boost::make_shared<SquareBoundaryCondition>(&population);
                      simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

                      MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
                      p_modifier->SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKiller");
                      simulator.AddSimulationModifier(p_modifier);

                      /* We then set an end time and run the simulation */
                      simulator.SetEndTime(4.20527000050075); // was 1.0075
                      simulator.Solve();
                      PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
                  }


    void NoTestLongerCapsuleSimulationWithDivisionAndMachine() throw (Exception)
    {


           const unsigned num_nodes = 1u;
           auto p_rand_gen = RandomNumberGenerator::Instance();

           // Create some capsules
           std::vector<Node<2>*> nodes;

           nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
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
           mesh.ConstructNodesWithoutMesh(nodes, 100.0);
           c_vector<double, 4> domain_size;
           domain_size[0] = -1000.0;
           domain_size[1] = 1000.0;
           domain_size[2] = -1000.0;
           domain_size[3] = 1000.0;
           mesh.SetInitialBoxCollection(domain_size, 10.0);

           mesh.GetNode(0u)->AddNodeAttribute(0.0);
           mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
           mesh.GetNode(0u)->rGetNodeAttributes()[NA_ANGLE] = 0.0;
           mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
           mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

           for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
           {
               mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
               mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
               mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_ANGLE] = 2.0 * M_PI * p_rand_gen->ranf();
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
        	   p_model->SetMaxCellCycleDuration(2.0);
               CellPtr p_cell(new Cell(p_state, p_model));
               p_cell->SetCellProliferativeType(p_type);

               MAKE_PTR(TypeSixMachineProperty, p_property);
				p_property->rGetMachineData().emplace_back(std::pair<unsigned, double>(i%4, 0.0));
					p_cell->AddCellProperty(p_property);

//               double birth_time = -RandomNumberGenerator::Instance()->ranf();
               p_cell->SetBirthTime(0.0);
               cells.push_back(p_cell);


           }

           // Create cell population
           NodeBasedCellPopulation<2> population(mesh, cells);

           population.AddCellWriter<CellIdWriter>();
           population.AddCellWriter<CapsuleOrientationWriter>();
           population.AddCellWriter<CapsuleScalingWriter>();

           boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
           population.SetCentreBasedDivisionRule(p_division_rule);

           // Create simulation
           OffLatticeSimulation<2> simulator(population);
           simulator.SetOutputDirectory("TestLongerCapsuleSimulationWithDivision");
           double dt = 1.0/1200.0;
           simulator.SetDt(dt);
           simulator.SetSamplingTimestepMultiple(10);

           auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
           simulator.SetNumericalMethod(p_numerical_method);

           MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
           p_modifier->SetOutputDirectory("TestLongerCapsuleSimulationWithDivision");
           simulator.AddSimulationModifier(p_modifier);
//
//           MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
//           simulator.AddCellKiller(p_killer);
           /*
            * We now create a force law and pass it to the simulation
            * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
            */
           auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
           simulator.AddForce(p_capsule_force);
//
           auto p_boundary_condition = boost::make_shared<SquareBoundaryCondition>(&population);
           simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

           /* We then set an end time and run the simulation */
           simulator.SetEndTime(3.0075); // was 1.0075
           std::cout<< "BEfore solve " << "\n";
           simulator.Solve();
           PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
       }
};

#endif /*TESTCAPSULESIMULATION_HPP_*/
