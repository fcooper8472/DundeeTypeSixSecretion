#include "TypeSixMachineProperty.hpp"

TypeSixMachineProperty::TypeSixMachineProperty()
    : AbstractCellProperty()
{
}

TypeSixMachineProperty::~TypeSixMachineProperty()
{
}

std::vector<std::pair<unsigned, double> >& TypeSixMachineProperty::rGetMachineData()
{
    return mMachineData;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TypeSixMachineProperty)
