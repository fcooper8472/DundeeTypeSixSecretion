
#ifndef CAPSULEFORCE_HPP_
#define CAPSULEFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law between two capsules (cylinder with hemispherical caps)
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class CapsuleForce : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestCapsuleForce;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

    /**
     * Calculate the shortest distance between two capsules.
     * @param rNodeA the node at the centre of mass of the first capsule
     * @param rNodeB the node at the centre of mass of the second capsule
     * @return the shortest distance
     */
    double CalculateDistanceBetweenCapsules(Node<SPACE_DIM>& rNodeA, Node<SPACE_DIM>& rNodeB);

    /**
     * Calculate the overlap between two capsules. This is the sum of the radii minus the shortest distance between
     * the capsules.
     * @param rNodeA the node at the centre of mass of the first capsule
     * @param rNodeB the node at the centre of mass of the second capsule
     * @param distance the shortest distance between the two capsules
     * @return the overlap
     */
    double CalculateOverlapBetweenCapsules(Node<SPACE_DIM>& rNodeA, Node<SPACE_DIM>& rNodeB, const double distance);

public:

    /**
     * Constructor.
     */
    CapsuleForce();

    /**
     * Destructor.
     */
    virtual ~CapsuleForce() = default;

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by AddForceContribution()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleForce)

#endif /*CAPSULEFORCE_HPP_*/
