
#include "TypeSixMachineModifier.hpp"
#include "TypeSixMachineProperty.hpp"

template<unsigned DIM>
TypeSixMachineModifier<DIM>::TypeSixMachineModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
TypeSixMachineModifier<DIM>::~TypeSixMachineModifier()
{
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    ///\todo Make sure the cell population is updated?
    //rCellPopulation.Update();

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::set<std::pair<unsigned, double> >& data = p_property->rGetMachineData();

        ///\todo change parameters to be member variables
        double k_1 = 1.0;
        double k_2 = 1.0;
        double k_3 = 1.0;
        double k_4 = 1.0;
        double k_5 = 1.0;
        double k_6 = 1.0;
        double k_7 = 1.0;
        double k_8 = 1.0;

	    double dt = SimulationTime::Instance()->GetTimeStep();
        assert((k_1 + k_2 + k_3 + k_4 + k_5 + k_6 + k_7 + k_8)*dt <= 1.0);

        for (std::set<std::pair<unsigned, double> >::iterator it = data.begin();
             it != data.end();
             )
        {
	        unsigned old_state = (*it).first;
	        unsigned new_state = UNSIGNED_UNSET;
		    double r = RandomNumberGenerator::Instance()->ranf();
		
		    switch (old_state)
		    {
		        case 1:
		            if (r < k_2*dt)
		            {
		                new_state = 0;
		            }
		            else if (r < (k_2 + k_3)*dt)
		            {
		                new_state = 2;
		            }
		            break;
		        case 2:
		            if (r < k_4*dt)
		            {
		                new_state = 1;
		            }
		            else if (r < (k_4 + k_5)*dt)
		            {
		                new_state = 3;
		            }
		            break;
		        case 3:
		            if (r < k_6*dt)
		            {
		                new_state = 2;
		            }
		            else if (r < (k_6 + k_7)*dt)
		            {
		                new_state = 4;
		            }
		            break;
		        case 4:
		            if (r < k_8*dt)
		            {
		                new_state = 8;
		            }
		            break;
		        default:
		            NEVER_REACHED;
		    }
		    data.first = new_state;
		    
		    if (new_state == 0)
		    {
		        data.erase(it);
		    }
		    else
		    {
		        ++it;
		    }
		}
		
		// Create a machine?
	    double r = RandomNumberGenerator::Instance()->ranf();
	
        if (r < k_1*dt)
        {
            double theta = 2*M_PI*RandomNumberGenerator::Instance()->ranf();
        	std::pair<unsigned, double> machine(1, theta);
            data.insert(machine);
        }
    }
}

template<unsigned DIM>
void TypeSixMachineModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class TypeSixMachineModifier<1>;
template class TypeSixMachineModifier<2>;
template class TypeSixMachineModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeSixMachineModifier)

