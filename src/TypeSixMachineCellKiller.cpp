
#include "TypeSixMachineCellKiller.hpp"

template<unsigned DIM>
TypeSixMachineCellKiller<DIM>::TypeSixMachineCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    // Iterate over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a TypeSixMachineProperty");
        }
        boost::shared_ptr<TypeSixMachineProperty> p_property = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
        std::vector<std::pair<unsigned, double> >& r_data = p_property->rGetMachineData();

		// Create a new vector to store all pairs less any we might throw away
		std::vector<std::pair<unsigned, double> > new_data;
		new_data.reserve(r_data.size() + 1);

        // Iterate over machines in this cell
        for (auto& r_pair : r_data)
        {
            // If this machine is ready to kill a cell...
            unsigned state = r_pair.first;
            if (state == 4)
            {
                // ...check if any neighbouring cells are close enough to kill...
                
                
                double theta = r_pair.second;
                
                if (///\todo)
                {                   
                    ///\todo Kill the neighbouring cell
                    
                    // Note: In this case we don't store this machine's data in new_data,
                    // since the machine is assumed to be destroyed upon killing the 
                    // neighbouring cell.
                }
                else
                {
		            new_data.emplace_back(std::pair<unsigned, double>(r_pair));                    
                }
            }
            else
            {
		        new_data.emplace_back(std::pair<unsigned, double>(r_pair));
		    }
        }

		r_data = new_data;
    }
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class TypeSixMachineCellKiller<1>;
template class TypeSixMachineCellKiller<2>;
template class TypeSixMachineCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeSixMachineCellKiller)
