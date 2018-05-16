
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
double CapsuleForce<ELEMENT_DIM,SPACE_DIM>::CalculateOverlapBetweenCapsules(
        Node<SPACE_DIM>& rNodeA,
        Node<SPACE_DIM>& rNodeB,
		const double& rShortestDistance)
{
    const double radius_a = rNodeA.rGetNodeAttributes()[NA_RADIUS];
    const double radius_b = rNodeB.rGetNodeAttributes()[NA_RADIUS];

    return radius_a + radius_b - rShortestDistance;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CapsuleForce<ELEMENT_DIM, SPACE_DIM>::CalculateForceDirectionAndContactPoints(Node<SPACE_DIM>& rNodeA,
                                                                                     Node<SPACE_DIM>& rNodeB,
                                                                                     c_vector<double, SPACE_DIM>& rVecAToB,
                                                                                     double& rContactDistA,
                                                                                     double& rContactDistB)
{


	auto location_a = rNodeA.rGetLocation();
	auto location_b = rNodeB.rGetLocation();

	const double angle_theta_a = rNodeA.rGetNodeAttributes()[NA_THETA];
	const double angle_theta_b = rNodeB.rGetNodeAttributes()[NA_THETA];

	double angle_phi_a;
	double angle_phi_b;

	const double length_a = rNodeA.rGetNodeAttributes()[NA_LENGTH];
	const double length_b = rNodeB.rGetNodeAttributes()[NA_LENGTH];

	if (SPACE_DIM==3u)
	{
		angle_phi_a = rNodeA.rGetNodeAttributes()[NA_PHI];
		angle_phi_b = rNodeB.rGetNodeAttributes()[NA_PHI];
	}

	c_vector<double,SPACE_DIM> segment_a_point_1;
	c_vector<double,SPACE_DIM> segment_a_point_2;
	c_vector<double,SPACE_DIM> segment_b_point_1;
	c_vector<double,SPACE_DIM> segment_b_point_2;

	// The line segments are backwards. If the centre point of a is location_a, then would segment_a_point_1 actually be
	// calculated by subtracting half the length? Like
	//segment_a_point_1[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
	//segment_a_point_1[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
	//segment_a_point_1[2] = location_a[2] - 0.5 * length_a * cos(angle_phi_a);
	//And segment_a_point_2 be the one added? This might account for the differing signs in the test for parallels?

	if (SPACE_DIM==3u)
	{
		segment_a_point_1[0] = location_a[0] + 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
		segment_a_point_1[1] = location_a[1] + 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
		segment_a_point_1[2] = location_a[2] + 0.5 * length_a * cos(angle_phi_a);
	//        PRINT_VECTOR(segment_a_point_1);

		segment_a_point_2[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
		segment_a_point_2[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
		segment_a_point_2[2] = location_a[2] - 0.5 * length_a * cos(angle_phi_a);
	//        PRINT_VECTOR(segment_a_point_2);

		segment_b_point_1[0] = location_b[0] + 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
		segment_b_point_1[1] = location_b[1] + 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
		segment_b_point_1[2] = location_b[2] + 0.5 * length_b * cos(angle_phi_b);
	//        PRINT_VECTOR(segment_b_point_1);

		segment_b_point_2[0] = location_b[0] - 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
		segment_b_point_2[1] = location_b[1] - 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
		segment_b_point_2[2] = location_b[2] - 0.5 * length_b * cos(angle_phi_b);
	}
	else
	{
		segment_a_point_1[0] = location_a[0] + 0.5 * length_a * cos(angle_theta_a);
		segment_a_point_1[1] = location_a[1] + 0.5 * length_a * sin(angle_theta_a);
	//        PRINT_VECTOR(segment_a_point_1);

		segment_a_point_2[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a);
		segment_a_point_2[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a);
	//        PRINT_VECTOR(segment_a_point_2);

		segment_b_point_1[0] = location_b[0] + 0.5 * length_b * cos(angle_theta_b);
		segment_b_point_1[1] = location_b[1] + 0.5 * length_b * sin(angle_theta_b);
	//        PRINT_VECTOR(segment_b_point_1);

		segment_b_point_2[0] = location_b[0] - 0.5 * length_b * cos(angle_theta_b);
		segment_b_point_2[1] = location_b[1] - 0.5 * length_b * sin(angle_theta_b);
	}
//        PRINT_VECTOR(segment_b_point_2);
	// Copyright 2001 softSurfer, 2012 Dan Sunday
	// This code may be freely used, distributed and modified for any purpose
	// providing that this copyright notice is included with it.
	// SoftSurfer makes no warranty for this code, and cannot be held
	// liable for any real or imagined damage resulting from its use.
	// Users of this code must verify correctness for their application.

	/* Terminology in this section:
	 * Segment 1 given by
	 *  P(s) = P_0 + s(P_1-P_0) = P_0 + s*vector_u
	 *
	 * Segment 2 given by
	 *  Q(t) = Q_0 + t(Q_1-Q_0) = Q_0 + t*vector_v
	 *
	 * w(s,t) = P(s)-Q(t) a generic vector between the two lines. Find w(s,t) s.t. it has a minimum length.
	 *
	 * closest points are
	 * P(sc)
	 * Q(tc)
	 * that minimise w(sc,tc) = minimum distance.
	 *
	 * Explained further here: http://geomalgorithms.com/a07-_distance.html
	 */


	//            Vector   u = S1.P1 - S1.P0;
	//            Vector   v = S2.P1 - S2.P0;
	//            Vector   w = S1.P0 - S2.P0;
	c_vector<double, SPACE_DIM> u = segment_a_point_2 - segment_a_point_1; // distance along here parameterised by 0 <= s <= 1
	c_vector<double, SPACE_DIM> v = segment_b_point_2 - segment_b_point_1; // distance along here parameterised by 0 <= t <= 1
	c_vector<double, SPACE_DIM> w = segment_a_point_1 - segment_b_point_1;

//        PRINT_VECTOR(u);
//        PRINT_VECTOR(v);
//        PRINT_VECTOR(w);

	double    a = inner_prod(u,u);         // always >= 0
	double    b = inner_prod(u,v);
	double    c = inner_prod(v,v);         // always >= 0
	double    d = inner_prod(u,w);
	double    e = inner_prod(v,w);
	double    D = a*c - b*b;        // always >= 0
	double    sc, sN, sD = D;       // sc = sN[umerator] / sD[eterminant], default sD (for closest point not on end of segment) = D >= 0
	double    tc, tN, tD = D;       // tc = tN[umerator] / tD[eterminant], default tD (for closest point not on end of segment) = D >= 0

	const double SMALL_NUM = 1e-9;

	//PRINT_VARIABLE(D);

	if (D < SMALL_NUM)
	{   // parallel lines
		sN = 0.0; // force using point P0 on segment S1
		sD = 1.0; // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
		//std::cout << "TRACE: parallel lines\n" << std::flush;

		// Special case that we have added for parallel lines that
		if (std::fabs(tN) < SMALL_NUM)
		{
			//std::cout << "TRACE: small sN and tN\n" << std::flush;
			sN = 0.5;
			tN = 0.5*tD;
		}
	}
	else
	{   // get the closest points on the infinite lines drawn by the segments
		sN = (b * e - c * d);
		tN = (a * e - b * d);

		// Work out if these are off the ends of the segments
		if (sN < 0) // sc < 0 => the s=0 on segment A is a closest point
		{
			sN = 0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) // sc > 1  => the s=1 on Segment A is closest point
		{
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0) // tc < 0 => the t=0 on Segment B is closest point
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
	else if (tN > tD)// tc > 1  => the t=1 on segment B is closest point
	{
		tN = tD;
		if ((-d + b) < 0)
		{
			//std::cout << "TRACE: setting sN=0" << std::endl;
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

	// finally do the division by determinant to get closest point parameters sc and tc
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

	// vector between the two closest points
	rVecAToB = -(w + (sc * u) - (tc * v));  // =  Segment_A(sc) - Segment_B(tc)

	rContactDistA = -(sc-0.5)*length_a;
	rContactDistB = -(tc-0.5)*length_b;

	double shortest_distance = norm_2(rVecAToB);   // return the closest distance

	// A zero vector is OK if they are intersecting...
	// We produce a unit vector if not.
	if (shortest_distance>SMALL_NUM)
	{
		rVecAToB /= shortest_distance; // Turn rVecAToB into a unit vector.
	}


    return CalculateOverlapBetweenCapsules(rNodeA, rNodeB, shortest_distance); // return the overlap
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
        if (SPACE_DIM==3u)
        {
        	iter->rGetNodeAttributes()[NA_APPLIED_PHI] = 0.0;
        }
    }

    // Calculate force and applied angle contributions from each pair
    for (auto& node_pair : p_cell_population->rGetNodePairs())
    {
        Node<SPACE_DIM>& r_node_a = *(node_pair.first);
        Node<SPACE_DIM>& r_node_b = *(node_pair.second);

        c_vector<double, SPACE_DIM> force_direction_a_to_b;
		double contact_dist_a;
		double contact_dist_b;



		double overlap = CalculateForceDirectionAndContactPoints(r_node_a,
																 r_node_b,
																 force_direction_a_to_b,
																 contact_dist_a,
																 contact_dist_b);

		if (overlap > 0.0)
		{

            const double radius_a = r_node_a.rGetNodeAttributes()[NA_RADIUS];
            const double radius_b = r_node_b.rGetNodeAttributes()[NA_RADIUS];
            double force_magnitude = CalculateForceMagnitude(overlap, radius_a, radius_b);


            c_vector<double, SPACE_DIM> force_a_b = force_direction_a_to_b * force_magnitude;
            c_vector<double, SPACE_DIM> force_b_a = -1.0 * force_a_b;

            const double angle_theta_a = r_node_a.rGetNodeAttributes()[NA_THETA];
            const double angle_theta_b = r_node_b.rGetNodeAttributes()[NA_THETA];

            c_vector<double, SPACE_DIM> torque_vec_a;
            c_vector<double, SPACE_DIM> torque_vec_b;
            c_vector<double, SPACE_DIM> cross_torque_vec;

            if (SPACE_DIM==2u)
            {
				torque_vec_a[0] = contact_dist_a * cos(angle_theta_a);
				torque_vec_a[1] = contact_dist_a * sin(angle_theta_a);
				torque_vec_b[0] = contact_dist_b * cos(angle_theta_b);
				torque_vec_b[1] = contact_dist_b * sin(angle_theta_b);
            }
            else
            {
            	const double angle_phi_a = r_node_a.rGetNodeAttributes()[NA_PHI];
            	const double angle_phi_b = r_node_b.rGetNodeAttributes()[NA_PHI];

        		torque_vec_a[0] = contact_dist_a * cos(angle_theta_a) * sin(angle_phi_a);
				torque_vec_a[1] = contact_dist_a * sin(angle_theta_a) * sin(angle_phi_a);
				torque_vec_a[2] = contact_dist_a * cos(angle_phi_a);

				torque_vec_b[0] = contact_dist_b * cos(angle_theta_b) * sin(angle_phi_b);
				torque_vec_b[1] = contact_dist_b * sin(angle_theta_b) * sin(angle_phi_b);
				torque_vec_b[2] = contact_dist_b * cos(angle_phi_b);
            }

			// Calculate the 2D cross product of two vectors
        	//cross_torque_vec = torque_vec_a[0]*torque_vec_b[1]-torque_vec_b[0]*torque_vec_a[1];
			auto cross_product = [](c_vector<double, SPACE_DIM> u, c_vector<double, SPACE_DIM> v) -> double
			{
				return u[0] * v[1] - u[1] * v[0];
			};

			// Calculate the 3D cross product of two vectors
			auto cross_product_3d = [](c_vector<double, SPACE_DIM> u, c_vector<double, SPACE_DIM> v) -> c_vector<double, SPACE_DIM>
			{
            	//Standard cross product where u={u1,u2,u3} and v={v1,v2,v3}
            	//uxv = {u2v3-u3v2,u3v1-u1v3,u1v2-u2v1}
//            	cross_torque_vec[0] = torque_vec_a[1]*torque_vec_b[2] - torque_vec_a[2]*torque_vec_b[1];
//            	cross_torque_vec[1] = torque_vec_a[2]*torque_vec_b[0] - torque_vec_a[0]*torque_vec_b[2];
//            	cross_torque_vec[2] = torque_vec_a[0]*torque_vec_b[1] - torque_vec_a[1]*torque_vec_b[0];
					c_vector<double, SPACE_DIM> c;
					c[0] = u[1]*v[2]-u[2]*v[1];
					c[1] = u[2]*v[0]-u[0]*v[2];
					c[2] = u[0]*v[1]-u[1]*v[0];
					return c;
			};


			if (SPACE_DIM==2u)
			{
				r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product(torque_vec_a, force_b_a);
				r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product(torque_vec_b, force_a_b);
			}
			else
			{
				r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product_3d(torque_vec_a, force_b_a)[2];
				r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product_3d(torque_vec_b, force_a_b)[2];
				r_node_a.rGetNodeAttributes()[NA_APPLIED_PHI] -= cross_product_3d(torque_vec_a, force_b_a)[0];
				r_node_b.rGetNodeAttributes()[NA_APPLIED_PHI] -= cross_product_3d(torque_vec_b, force_a_b)[0];
			}

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
