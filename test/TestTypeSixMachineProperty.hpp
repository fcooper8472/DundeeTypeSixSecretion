
#ifndef TESTTYPESIXMACHINEPROPERTY_HPP_
#define TESTTYPESIXMACHINEPROPERTY_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "TypeSixMachineProperty.hpp"

class TestTypeSixMachineProperty : public AbstractCellBasedTestSuite
{
public:

    void TestProperty() throw (Exception)
    {
        // Create cell property and test that the data structure is empty
        MAKE_PTR(TypeSixMachineProperty, p_property);        
        TS_ASSERT_EQUALS(p_property->rGetMachineData().empty(), true);

        // Add the property to a cell and test that the data structure is still empty 
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        UniformCellCycleModel* p_model = new UniformCellCycleModel();

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->AddCellProperty(p_property);
        
        TS_ASSERT_EQUALS(p_cell->rGetCellPropertyCollection().HasProperty(p_property), true);
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);

        CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<TypeSixMachineProperty>();
        boost::shared_ptr<TypeSixMachineProperty> p_property_from_cell = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        TS_ASSERT_EQUALS(p_property_from_cell->rGetMachineData().empty(), true);

        // Test that we can add some data to the data structure
        std::set<std::pair<unsigned, double> >& data = p_property_from_cell->rGetMachineData();
        std::pair<unsigned, double> machine(4, 0.5);
        
        data.insert(machine);
        TS_ASSERT_EQUALS(data.empty(), false);

        collection = p_cell->rGetCellPropertyCollection().GetProperties<TypeSixMachineProperty>();
        p_property_from_cell = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        TS_ASSERT_EQUALS(p_property_from_cell->rGetMachineData().empty(), false);

        // Test that we can recover the data from the data structure
        std::set<std::pair<unsigned, double> >& data_from_cell = p_property_from_cell->rGetMachineData();
        TS_ASSERT_EQUALS(data_from_cell.size(), 1u);

        std::pair<unsigned, double> data_pair = *(data_from_cell.begin());
        unsigned data_1 = data_pair.first;
        TS_ASSERT_EQUALS(data_1, 4u);
        double data_2 = data_pair.second;
        TS_ASSERT_DELTA(data_2, 0.5, 1e-6);
    }

    ///\todo test that property can be used in a simulation with different values for each cell
    void TestSimulationWithProperty() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a 'nodes only' mesh
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create some cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Associate the mesh and cells in a cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Pass the cell population to a simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestTypeSixMachineProperty");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        // Pass a force law to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        // Run the simulation
        simulator.Solve();
    }
};

#endif /* TESTTYPESIXMACHINEPROPERTY_HPP_ */
