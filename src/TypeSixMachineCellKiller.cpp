
#include "TypeSixMachineCellKiller.hpp"

template<unsigned DIM>
TypeSixMachineCellKiller<DIM>::TypeSixMachineCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
}

template<unsigned DIM>
void TypeSixMachineCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
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
