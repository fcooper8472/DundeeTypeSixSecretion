
#ifndef TYPESIXMACHINEPROPERTY_HPP_
#define TYPESIXMACHINEPROPERTY_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"
#include "PetscTools.hpp"
#include <set>

class TypeSixMachineProperty : public AbstractCellProperty
{
private:

    std::vector<std::pair<unsigned, double> > mMachineData;

public:

    /**
     * Constructor.
     */
    TypeSixMachineProperty();

    /**
     * Destructor.
     */
    virtual ~TypeSixMachineProperty();

    /**
     * @return #mMachineData
     */
    std::vector<std::pair<unsigned, double> >& rGetMachineData();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TypeSixMachineProperty)

#endif /* TYPESIXMACHINEPROPERTY_HPP_ */
