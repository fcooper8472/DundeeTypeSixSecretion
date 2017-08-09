
#ifndef CAPSULEFORCE_HPP_
#define CAPSULEFORCE_HPP_

#include "AbstractForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law between two capsules (cylinder with hemispherical caps)
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class CapsuleForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
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
     * @param shortestDistance the shortest distance between the two capsules
     * @return the overlap
     */
    double CalculateOverlapBetweenCapsules(Node<SPACE_DIM>& rNodeA, Node<SPACE_DIM>& rNodeB, const double shortestDistance);

    /**
     * Calculate the direction of the force between two overlapping capsules, and the distance along the rod from each
     * capsule's centre of mass to the point at which the force will act.
     * @param rNodeA the node at the centre of mass of the first capsule
     * @param rNodeB the node at the centre of mass of the second capsule
     * @param shortestDistance the shortest distance between the two capsules
     * @param rVecAToB filled in as a unit vector from the contact point on capsule A to that on capsule B
     * @param rContactDistA filled in as the distance from the centre of mass of capsule A to contact point
     * @param rContactDistB filled in as the distance from the centre of mass of capsule B to contact point
     */
    void CalculateForceDirectionAndContactPoints(Node<SPACE_DIM>& rNodeA,
                                                 Node<SPACE_DIM>& rNodeB,
                                                 const double shortestDistance,
                                                 c_vector<double, SPACE_DIM>& rVecAToB,
                                                 double& rContactDistA,
                                                 double& rContactDistB);

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
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

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
