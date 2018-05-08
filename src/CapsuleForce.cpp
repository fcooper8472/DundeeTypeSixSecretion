
#include "CapsuleForce.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"

#ifdef CHASTE_VTK
#include <vtkLine.h>
#endif // CHASTE_VTK

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <cmath>

#include "Debug.hpp"



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM, SPACE_DIM>::SetYoungModulus(double youngModulus)

{
	mYoungModulus=youngModulus;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM, SPACE_DIM>::GetYoungModulus()

{
	return mYoungModulus;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CapsuleForce<ELEMENT_DIM, SPACE_DIM>::CapsuleForce()
        : AbstractForce<ELEMENT_DIM, SPACE_DIM>(),
          mYoungModulus(100.0)
{
    // Has to be either element and space dimensions are both 2 or both 3.
    assert((ELEMENT_DIM == 2u && SPACE_DIM == 2u) || (ELEMENT_DIM == 3u && SPACE_DIM == 3u));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateDistanceBetweenCapsules(
        Node<SPACE_DIM>& rNodeA,
        Node<SPACE_DIM>& rNodeB)
{
    auto location_a = rNodeA.rGetLocation();
    auto location_b = rNodeB.rGetLocation();

    if (SPACE_DIM==2u)
    {
        using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
        using geom_segment = boost::geometry::model::segment<geom_point>;

        const double angle_theta_a = rNodeA.rGetNodeAttributes()[NA_THETA];
        const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];

        const double angle_theta_b = rNodeB.rGetNodeAttributes()[NA_THETA];
        const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

        geom_point capsule_a_end_1(location_a[0] + 0.5 * length_a * cos(angle_theta_a),
                                   location_a[1] + 0.5 * length_a * sin(angle_theta_a));

        geom_point capsule_a_end_2(location_a[0] - 0.5 * length_a * cos(angle_theta_a),
                                   location_a[1] - 0.5 * length_a * sin(angle_theta_a));

        geom_point capsule_b_end_1(location_b[0] + 0.5 * length_b * cos(angle_theta_b),
                                   location_b[1] + 0.5 * length_b * sin(angle_theta_b));

        geom_point capsule_b_end_2(location_b[0] - 0.5 * length_b * cos(angle_theta_b),
                                   location_b[1] - 0.5 * length_b * sin(angle_theta_b));

        geom_segment capsule_axis_a(capsule_a_end_1, capsule_a_end_2);
        geom_segment capsule_axis_b(capsule_b_end_1, capsule_b_end_2);

        return boost::geometry::distance(capsule_axis_a, capsule_axis_b);
    }
    else
    {
        EXCEPT_IF_NOT(SPACE_DIM==3u);

        const double angle_theta_a = rNodeA.rGetNodeAttributes()[NA_THETA];
        const double angle_phi_a = rNodeA.rGetNodeAttributes()[NA_PHI];
        const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];

        //PRINT_3_VARIABLES(length_a,angle_theta_a,angle_phi_a);

        const double angle_theta_b = rNodeB.rGetNodeAttributes()[NA_THETA];
        const double angle_phi_b = rNodeB.rGetNodeAttributes()[NA_PHI];
        const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

        //PRINT_3_VARIABLES(length_b,angle_theta_b,angle_phi_b);

/*      // THIS IS A BOOST IMPLEMENTATION THAT THROWS A `NOT YET IMPLEMENTED' ERROR ON 1.58.
 *      using geom_point = boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
        using geom_segment = boost::geometry::model::segment<geom_point>;

        geom_point capsule_a_end_1(location_a[0] + 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a),
                                   location_a[1] + 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a),
                                   location_a[2] + 0.5 * length_a * cos(angle_phi_a));

        geom_point capsule_a_end_2(location_a[0] - 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a),
                                   location_a[1] - 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a),
                                   location_a[2] - 0.5 * length_a * cos(angle_phi_a));

        geom_point capsule_b_end_1(location_b[0] + 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b),
                                   location_b[1] + 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b),
                                   location_b[2] + 0.5 * length_b * cos(angle_phi_b));

        geom_point capsule_b_end_2(location_b[0] - 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b),
                                   location_b[1] - 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b),
                                   location_b[2] - 0.5 * length_b * cos(angle_phi_b));

        geom_segment capsule_axis_a(capsule_a_end_1, capsule_a_end_2);
        geom_segment capsule_axis_b(capsule_b_end_1, capsule_b_end_2);

        return boost::geometry::distance(capsule_axis_a, capsule_axis_b);

        */

        c_vector<double,3u> segment_a_point_1;
        segment_a_point_1[0] = location_a[0] + 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
        segment_a_point_1[1] = location_a[1] + 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
        segment_a_point_1[2] = location_a[2] + 0.5 * length_a * cos(angle_phi_a);
//        PRINT_VECTOR(segment_a_point_1);

        c_vector<double,3u> segment_a_point_2;
        segment_a_point_2[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
        segment_a_point_2[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
        segment_a_point_2[2] = location_a[2] - 0.5 * length_a * cos(angle_phi_a);
//        PRINT_VECTOR(segment_a_point_2);

        c_vector<double,3u> segment_b_point_1;
        segment_b_point_1[0] = location_b[0] + 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
        segment_b_point_1[1] = location_b[1] + 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
        segment_b_point_1[2] = location_b[2] + 0.5 * length_b * cos(angle_phi_b);
//        PRINT_VECTOR(segment_b_point_1);

        c_vector<double,3u> segment_b_point_2;
        segment_b_point_2[0] = location_b[0] - 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
        segment_b_point_2[1] = location_b[1] - 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
        segment_b_point_2[2] = location_b[2] - 0.5 * length_b * cos(angle_phi_b);
//        PRINT_VECTOR(segment_b_point_2);
        // Copyright 2001 softSurfer, 2012 Dan Sunday
        // This code may be freely used, distributed and modified for any purpose
        // providing that this copyright notice is included with it.
        // SoftSurfer makes no warranty for this code, and cannot be held
        // liable for any real or imagined damage resulting from its use.
        // Users of this code must verify correctness for their application.


        //            Vector   u = S1.P1 - S1.P0;
        //            Vector   v = S2.P1 - S2.P0;
        //            Vector   w = S1.P0 - S2.P0;
        c_vector<double, 3u> u = segment_a_point_2 - segment_a_point_1;
        c_vector<double, 3u> v = segment_b_point_2 - segment_b_point_1;
        c_vector<double, 3u> w = segment_a_point_1 - segment_b_point_1;

//        PRINT_VECTOR(u);
//        PRINT_VECTOR(v);
//        PRINT_VECTOR(w);

        double    a = inner_prod(u,u);         // always >= 0
        double    b = inner_prod(u,v);
        double    c = inner_prod(v,v);         // always >= 0
        double    d = inner_prod(u,w);
        double    e = inner_prod(v,w);
        double    D = a*c - b*b;        // always >= 0
        double    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
        double    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

        const double SMALL_NUM = 1e-9;

        //PRINT_VARIABLE(D);

        if (D < SMALL_NUM)
        {   // parallel lines
            sN = 0;
            sD = 1;
            tN = e;
            tD = c;
        }
        else
        {
            sN = (b * e - c * d);
            tN = (a * e - b * d);
            if (sN < 0)
            {
                sN = 0;
                tN = e;
                tD = c;
            }
            else if (sN > sD)
            {
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN < 0)
        {
            tN = 0;
            if (-d < 0)
            {
                sN = 0;
            }
            else if (-d > a)
            {
                sN = sD;
            }
            else
            {
                sN = -d;
                sD = a;
            }
        }
        else if (tN > tD)
        {
            tN = tD;
            if ((-d + b) < 0)
            {
                sN = 0;
            }
            else if ((-d + b) > a)
            {
                sN = sD;
            }
            else
            {
                sN = (-d + b);
                sD = a;
            }
        }

        if (std::fabs(sN) < SMALL_NUM)
        {
            sc = 0;
        }
        else
        {
            sc = sN / sD;
        }

        if (std::fabs(tN) < SMALL_NUM)
        {
            tc = 0;
        }
        else
        {
            tc = tN / tD;
        }

        // get the difference of the two closest points
        c_vector<double, 3u> dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)

        return norm_2(dP);   // return the closest distance
     }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateOverlapBetweenCapsules(
        Node<SPACE_DIM>& rNodeA,
        Node<SPACE_DIM>& rNodeB,
        const double shortestDistance)
{
    const double radius_a = rNodeA.rGetNodeAttributes()[NA_RADIUS];
    const double radius_b = rNodeB.rGetNodeAttributes()[NA_RADIUS];

    return radius_a + radius_b - shortestDistance;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM, SPACE_DIM>::CalculateForceAttributes(Node<SPACE_DIM>& rNodeA,
                                                                    Node<SPACE_DIM>& rNodeB,
                                                                    c_vector<double, SPACE_DIM>& rForceDirection,
                                                                    double& rOverlap,
                                                                    double& rContactDistA,
                                                                    double& rContactDistB)
{
    EXCEPT_IF_NOT(SPACE_DIM == 2);

#ifdef CHASTE_VTK

    const auto& location_a = rNodeA.rGetLocation();
    const auto& location_b = rNodeB.rGetLocation();

    const double angle_theta_a = rNodeA.rGetNodeAttributes()[NA_THETA];
    const double radius_a = rNodeA.rGetNodeAttributes()[NA_RADIUS];
    const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];

    const double angle_theta_b = rNodeB.rGetNodeAttributes()[NA_THETA];
    const double radius_b = rNodeB.rGetNodeAttributes()[NA_RADIUS];
    const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

    // Calculate the ends of capsule a, with a default value of z = 0
    double capsule_a_end_1[3] = {location_a[0] + 0.5 * length_a * cos(angle_theta_a),
                                 location_a[1] + 0.5 * length_a * sin(angle_theta_a), 0.0};
    double capsule_a_end_2[3] = {location_a[0] - 0.5 * length_a * cos(angle_theta_a),
                                 location_a[1] - 0.5 * length_a * sin(angle_theta_a), 0.0};

    // Calculate the ends of capsule b, with a default value of z = 0
    double capsule_b_end_1[3] = {location_b[0] + 0.5 * length_b * cos(angle_theta_b),
                                 location_b[1] + 0.5 * length_b * sin(angle_theta_b), 0.0};
    double capsule_b_end_2[3] = {location_b[0] - 0.5 * length_b * cos(angle_theta_b),
                                 location_b[1] - 0.5 * length_b * sin(angle_theta_b), 0.0};

    // Create two arrays to be filled by the vtk algorithm
    double closest_point_a[3];
    double closest_point_b[3];

    // Create two doubles to be filled by the vtk algorithm
    double parametric_a;
    double parametric_b;

    const double shortest_dist_squared = vtkLine::DistanceBetweenLineSegments(capsule_a_end_1, capsule_a_end_2,
                                                                              capsule_b_end_1, capsule_b_end_2,
                                                                              closest_point_a, closest_point_b,
                                                                              parametric_a, parametric_b);

    // Fill in the force direction
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        rForceDirection[dim] = closest_point_b[dim] - closest_point_a[dim];
    }
    rForceDirection /= norm_2(rForceDirection);

    // Fill in overlap and contact distances
    rOverlap = radius_a + radius_b - sqrt(shortest_dist_squared);
    rContactDistA = (parametric_a - 0.5) * length_a;
    rContactDistB = (parametric_b - 0.5) * length_b;

#else  // CHASTE_VTK not defined
    EXCEPTION("VTK is required for this calculation");
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM, SPACE_DIM>::CalculateForceDirectionAndContactPoints(Node<SPACE_DIM>& rNodeA,
                                                                                   Node<SPACE_DIM>& rNodeB,
                                                                                   const double shortestDistance,
                                                                                   c_vector<double, SPACE_DIM>& rVecAToB,
                                                                                   double& rContactDistA,
                                                                                   double& rContactDistB)
{


    EXCEPT_IF_NOT(SPACE_DIM == 2);

    using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using geom_segment = boost::geometry::model::segment<geom_point>;

    rContactDistA = DOUBLE_UNSET;
    rContactDistB = DOUBLE_UNSET;

    auto location_a = rNodeA.rGetLocation();
    auto location_b = rNodeB.rGetLocation();

    const double angle_theta_a = rNodeA.rGetNodeAttributes()[NA_THETA];
    const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];

    const double angle_theta_b = rNodeB.rGetNodeAttributes()[NA_THETA];
    const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

    geom_point capsule_a_end_1(location_a[0] + 0.5 * length_a * cos(angle_theta_a),
                               location_a[1] + 0.5 * length_a * sin(angle_theta_a));

    geom_point capsule_a_end_2(location_a[0] - 0.5 * length_a * cos(angle_theta_a),
                               location_a[1] - 0.5 * length_a * sin(angle_theta_a));

    geom_point capsule_b_end_1(location_b[0] + 0.5 * length_b * cos(angle_theta_b),
                               location_b[1] + 0.5 * length_b * sin(angle_theta_b));

    geom_point capsule_b_end_2(location_b[0] - 0.5 * length_b * cos(angle_theta_b),
                               location_b[1] - 0.5 * length_b * sin(angle_theta_b));

    geom_segment capsule_axis_a(capsule_a_end_1, capsule_a_end_2);
    geom_segment capsule_axis_b(capsule_b_end_1, capsule_b_end_2);

    c_vector<double, SPACE_DIM> end_point;
    c_vector<double, SPACE_DIM> centre_capsule;
    double angle_for_calculation;
    bool point_in_capsule_a;

    if (shortestDistance < 1e-12)
    {
    	EXCEPTION("Centre lines of capsules have collided - something is wrong!");
    }

    if (fabs(boost::geometry::distance(capsule_a_end_1, capsule_axis_b) - shortestDistance) < 1e-12)
    {
        end_point[0] = capsule_a_end_1.get<0>();
        end_point[1] = capsule_a_end_1.get<1>();
        centre_capsule = location_b;
        rContactDistA = 0.5 * length_a;
        angle_for_calculation = angle_theta_b;
        point_in_capsule_a = true;
    }
    else if(fabs(boost::geometry::distance(capsule_a_end_2, capsule_axis_b) - shortestDistance) < 1e-12)
    {
        end_point[0] = capsule_a_end_2.get<0>();
        end_point[1] = capsule_a_end_2.get<1>();
        centre_capsule = location_b;
        rContactDistA = -0.5 * length_a;
        angle_for_calculation = angle_theta_b;
        point_in_capsule_a = true;
    }
    else if(fabs(boost::geometry::distance(capsule_b_end_1, capsule_axis_a) - shortestDistance) < 1e-12)
    {
        end_point[0] = capsule_b_end_1.get<0>();
        end_point[1] = capsule_b_end_1.get<1>();
        centre_capsule = location_a;
        rContactDistB = 0.5 * length_b;
        angle_for_calculation = angle_theta_a;
        point_in_capsule_a = false;
    }
    else if(fabs(boost::geometry::distance(capsule_b_end_2, capsule_axis_a) - shortestDistance) < 1e-12)
    {
        end_point[0] = capsule_b_end_2.get<0>();
        end_point[1] = capsule_b_end_2.get<1>();
        centre_capsule = location_a;
        rContactDistB = -0.5 * length_b;
        angle_for_calculation = angle_theta_a;
        point_in_capsule_a = false;
    }
    else
    {
        NEVER_REACHED;
    }

    double cos_theta = cos(angle_for_calculation);
    double sin_theta = sin(angle_for_calculation);

    auto delta = end_point - centre_capsule;
    double other_length = (delta[0] * cos_theta + delta[1] * sin_theta) / (cos_theta * cos_theta + sin_theta * sin_theta);

    if (rContactDistA == DOUBLE_UNSET)
    {
        if (fabs(other_length) > 0.5 * length_a)
        {
            NEVER_REACHED;
            /*
             * We should never have to cap "other_length" in Capsule A: if we did, the implication would be that one of
             * the endpoints of Capsule A is involved in the closest distance, and "other_length" would then have to be
             * in Capsule B instead.
             */
            //other_length = copysign(0.5 * length_a, other_length);
        }

        rContactDistA = other_length;
    }
    else
    {
        if (fabs(other_length) > 0.5 * length_b)
        {
            other_length = copysign(0.5 * length_b, other_length);
        }

        rContactDistB = other_length;
    }

    c_vector<double, SPACE_DIM> cos_sin_theta;
    cos_sin_theta[0] = cos_theta;
    cos_sin_theta[1] = sin_theta;

    rVecAToB = delta - other_length * cos_sin_theta;

    if (point_in_capsule_a)
    {
        rVecAToB *= -1.0;
    }

    rVecAToB /= norm_2(rVecAToB);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceMagnitude(const double overlap,
                                                                    const double radiusA,
                                                                    const double radiusB)
{
    const double effective_radius = 2.0 * radiusA * radiusB / (radiusA + radiusB);
    const double force = 4.0 * mYoungModulus * pow(overlap, 1.5) * sqrt(effective_radius) / 3.0;

    // Horrific hack to stop explosions after division and before appropriate length is set!
    if (overlap > radiusA)
    {
    	return 0.0;
    }

    return force;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    auto p_cell_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);

    if (p_cell_population == nullptr)
    {
        EXCEPTION("Capsule force only works with AbstractCentreBasedCellPopulation");
    }

    // Set all applied angles back to zero
    for (auto iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
         iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
         ++iter)
    {
        iter->rGetNodeAttributes()[NA_APPLIED_THETA] = 0.0;
    }

    // Calculate force and applied angle contributions from each pair
    for (auto& node_pair : p_cell_population->rGetNodePairs())
    {
        Node<SPACE_DIM>& r_node_a = *(node_pair.first);
        Node<SPACE_DIM>& r_node_b = *(node_pair.second);

        double min_dist = CalculateDistanceBetweenCapsules(r_node_a, r_node_b);
        double overlap = CalculateOverlapBetweenCapsules(r_node_a, r_node_b, min_dist);

        if (overlap > 0.0)
        {
            c_vector<double, SPACE_DIM> force_direction_a_to_b;
            double contact_dist_a;
            double contact_dist_b;

            CalculateForceDirectionAndContactPoints(r_node_a,
                                                    r_node_b,
                                                    min_dist,
                                                    force_direction_a_to_b,
                                                    contact_dist_a,
                                                    contact_dist_b);

            const double radius_a = r_node_a.rGetNodeAttributes()[NA_RADIUS];
            const double radius_b = r_node_b.rGetNodeAttributes()[NA_RADIUS];
            double force_magnitude = CalculateForceMagnitude(overlap, radius_a, radius_b);

            c_vector<double, SPACE_DIM> force_a_b = force_direction_a_to_b * force_magnitude;
            c_vector<double, SPACE_DIM> force_b_a = -1.0 * force_a_b;

            const double angle_theta_a = r_node_a.rGetNodeAttributes()[NA_THETA];
            const double angle_theta_b = r_node_b.rGetNodeAttributes()[NA_THETA];

            c_vector<double, SPACE_DIM> torque_vec_a;
            torque_vec_a[0] = contact_dist_a * cos(angle_theta_a);
            torque_vec_a[1] = contact_dist_a * sin(angle_theta_a);

            c_vector<double, SPACE_DIM> torque_vec_b;
            torque_vec_b[0] = contact_dist_b * cos(angle_theta_b);
            torque_vec_b[1] = contact_dist_b * sin(angle_theta_b);

            // Calculate the 2D cross product of two vectors
            auto cross_product = [](c_vector<double, SPACE_DIM> a, c_vector<double, SPACE_DIM> b) -> double
            {
                return a[0] * b[1] - b[0] * a[1];
            };

            r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product(torque_vec_a, force_b_a);
            r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product(torque_vec_b, force_a_b);

            r_node_b.AddAppliedForceContribution(force_a_b);
            r_node_a.AddAppliedForceContribution(force_b_a);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CapsuleForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class CapsuleForce<1,1>;
template class CapsuleForce<1,2>;
template class CapsuleForce<2,2>;
template class CapsuleForce<1,3>;
template class CapsuleForce<2,3>;
template class CapsuleForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleForce)
