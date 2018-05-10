/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef FORWARDEULERNUMERICALMETHODFORCAPSULES_HPP_
#define FORWARDEULERNUMERICALMETHODFORCAPSULES_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractNumericalMethod.hpp"

/**
 * Implements forward Euler time stepping specific for Capsules,
 * where the position and angle are both updated each time step.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class ForwardEulerNumericalMethodForCapsules : public AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestNumericalMethodForCapsules;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM> >(*this);
    }

    /**
     * Calculate the 3D mass per unit volume of a 2D capsule (a cylinder capped by two hemispheres).
     * @param length the length of the cylinder
     * @param radius the radius of the cylinder and each hemisphere
     * @return the mass per unit volume of the capsule
     */
    double CalculateMassOfCapsule(const double length, const double radius);

    /**
     * Calculate the 3D moment of inertia of a 2D capsule (a cylinder capped by two hemispheres), about its centre.
     * @param length the length of the cylinder
     * @param radius the radius of the cylinder and each hemisphere
     * @return the moment of inertia of the capsule
     */
    double CalculateMomentOfInertiaOfCapsule(const double length, const double radius);

public:

    /**
     * Constructor.
     */
    ForwardEulerNumericalMethodForCapsules();

    /**
     * Destructor.
     */
    virtual ~ForwardEulerNumericalMethodForCapsules() = default;

    /**
     * Overridden UpdateAllNodePositions() method.
     *
     * @param dt Time step size
     */
    void UpdateAllNodePositions(double dt);

    /**
     * Overridden OutputNumericalMethodParameters() method.
     *
     * @param rParamsFile Reference to the parameter output filestream
     */
    virtual void OutputNumericalMethodParameters(out_stream& rParamsFile);

    bool mAxialCapsuleGrowth;

    void SetAxialCapsuleGrowth(bool);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethodForCapsules)

#endif /*FORWARDEULERNUMERICALMETHODFORCAPSULES_HPP_*/
