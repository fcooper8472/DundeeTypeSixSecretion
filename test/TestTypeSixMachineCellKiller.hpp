
#ifndef TESTTYPESIXMACHINECELLKILLER_HPP_
#define TESTTYPESIXMACHINECELLKILLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestTypeSixMachineCellKiller : public AbstractCellBasedTestSuite
{
public:

    void TestCellKiller() throw(Exception)
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
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create cell killer
        TypeSixMachineCellKiller<2> cell_killer(&cell_population);

        TS_ASSERT_EQUALS(cell_killer.GetIdentifier(), "TypeSixMachineCellKiller-2");

        ///\todo Test something
    }

///\todo test archiving and parameter output method
};

#endif /*TESTTYPESIXMACHINECELLKILLER_HPP_*/
