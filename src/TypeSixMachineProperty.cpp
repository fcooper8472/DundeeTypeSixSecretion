#include "TypeSixMachineProperty.hpp"

TypeSixMachineProperty::TypeSixMachineProperty()
    : AbstractCellProperty()
{
}

TypeSixMachineProperty::~TypeSixMachineProperty()
{
}

std::vector<std::pair<unsigned, std::vector<double>> >& TypeSixMachineProperty::rGetMachineData()
{
    return mMachineData;
}

void TypeSixMachineProperty::SetMachineData(std::vector<std::pair<unsigned, std::vector<double>> > machineData)
{
     mMachineData=machineData;
}

unsigned TypeSixMachineProperty::GetNumMachineFiresInThisTimeStep()
{
    return mNumMachineFiresInThisTimeStep;
}

void TypeSixMachineProperty::SetNumMachineFiresInThisTimeStep(unsigned numMachineFiresInThisTimeStep)
{
     mNumMachineFiresInThisTimeStep=numMachineFiresInThisTimeStep;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TypeSixMachineProperty)
