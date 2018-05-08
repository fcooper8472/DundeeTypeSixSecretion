
#ifndef _TESTCAPSULEFORCE_HPP_
#define _TESTCAPSULEFORCE_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "CapsuleForce.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "Node.hpp"
#include "NodeAttributes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsuleForce : public AbstractCellBasedTestSuite
{
public:

    void TestDistanceBetweenTwoCapsules2d() throw(Exception)
    {
        CapsuleForce<2, 2> force;

        // Two horizontal rods a distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{1.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{1.0, 2.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.25;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.25, 1e-6);
        }

        // One horizontal, one vertical, distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{1.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.75;

            Node<2> node_b(0u, std::vector<double>{1.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.05;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.2, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{0.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), 1.0, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{4.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 4.0;
            attributes_a[NA_RADIUS] = 0.12;

            Node<2> node_b(0u, std::vector<double>{4.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = atan(0.75);
            attributes_b[NA_LENGTH] = 4.0;
            attributes_b[NA_RADIUS] = 0.13;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.55, 1e-6);
        }
    }

    void TestDistanceBetweenTwoCapsules3d() throw(Exception)
    {
        CapsuleForce<3, 3> force;

        // Two horizontal rods a distance 2 from each other
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{1.0, 0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_PHI] = 1.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<3> node_b(0u, std::vector<double>{1.0, 2.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_PHI] = 1.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.25;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 2.0, 1e-12);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.25, 1e-6);
        }

        // One horizontal, one vertical, distance 2 from each other
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{1.0, 0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_PHI] = M_PI/2.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.75;

            Node<3> node_b(0u, std::vector<double>{1.0, 3.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_PHI] = M_PI/2.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.05;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 2.0, 1e-9);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.2, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_PHI] = M_PI/2.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_PHI] = M_PI/2.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 0.0, 1e-9);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), 1.0, 1e-6);
        }

        // 3,4,5 triangle angles for capsule B above capsule A which is horizontal.
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{4.0, 0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_PHI] = M_PI/2.0;
            attributes_a[NA_LENGTH] = 4.0;
            attributes_a[NA_RADIUS] = 0.12;

            Node<3> node_b(0u, std::vector<double>{4.0, 3.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = atan(0.75);
            attributes_b[NA_PHI] = M_PI/2.0;
            attributes_b[NA_LENGTH] = 4.0;
            attributes_b[NA_RADIUS] = 0.13;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 1.8, 1e-9);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.55, 1e-6);
        }

        // Intersecting cross shape, distance 0 from each other
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 1.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_PHI] = M_PI/2.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_PHI] = M_PI/2.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 1.0, 1e-9);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), 0.0, 1e-6);
        }

        // 3,4,5 triangle angles for capsule B above capsule A which is horizontal.
        {
            MARK;
            // index, {x, y, z}
            Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 1.0;
            attributes_a[NA_PHI] = M_PI/2.0;
            attributes_a[NA_LENGTH] = 4.0;
            attributes_a[NA_RADIUS] = 0.12;

            Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 1.0;
            attributes_b[NA_PHI] = M_PI/2.0-atan(0.75);
            attributes_b[NA_LENGTH] = 4.0;
            attributes_b[NA_RADIUS] = 0.13;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);
            TS_ASSERT_DELTA(d, 1.8, 1e-9);

            TS_ASSERT_DELTA(force.CalculateOverlapBetweenCapsules(node_a, node_b, d), -1.55, 1e-6);
        }
    }

    void TestCalculateForceDirectionAndContactPoints() throw(Exception)
    {
        CapsuleForce<2, 2> force;

        // Two horizontal rods a distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{10.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;

            Node<2> node_b(0u, std::vector<double>{0.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 2.0;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            double contact_dist_a;
            double contact_dist_b;
            c_vector<double, 2> vec_a_to_b;

            force.CalculateForceDirectionAndContactPoints(node_a, node_b, d, vec_a_to_b, contact_dist_a, contact_dist_b);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], -8.0 / sqrt(73.0), 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], 3.0 / sqrt(73.0), 1e-6);

            // Swap a->b for coverage
            force.CalculateForceDirectionAndContactPoints(node_b, node_a, d, vec_a_to_b, contact_dist_b, contact_dist_a);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 8.0 / sqrt(73.0), 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], -3.0 / sqrt(73.0), 1e-6);
        }

        // Two horizontal rods a distance 2 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{10.0, 0.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;

            Node<2> node_b(0u, std::vector<double>{0.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 18.0;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            double contact_dist_a;
            double contact_dist_b;
            c_vector<double, 2> vec_a_to_b;

            force.CalculateForceDirectionAndContactPoints(node_a, node_b, d, vec_a_to_b, contact_dist_a, contact_dist_b);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);

            // Swap a->b
            force.CalculateForceDirectionAndContactPoints(node_b, node_a, d, vec_a_to_b, contact_dist_b, contact_dist_a);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);
        }

        // 3-4-5 triangle-based orientation
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{4.0, 3.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = atan(0.75);
            attributes_a[NA_LENGTH] = 4.0;

            Node<2> node_b(0u, std::vector<double>{4.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 4.0;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            double contact_dist_a;
            double contact_dist_b;
            c_vector<double, 2> vec_a_to_b;

            force.CalculateForceDirectionAndContactPoints(node_a, node_b, d, vec_a_to_b, contact_dist_a, contact_dist_b);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);

            // Swap a->b for coverage
            force.CalculateForceDirectionAndContactPoints(node_b, node_a, d, vec_a_to_b, contact_dist_b, contact_dist_a);

            TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);
        }

        // 3-4-5 triangle-based orientation, with the horizontal capsule rotated 180 from the test above
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{4.0, 3.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = atan(0.75) + M_PI;
            attributes_a[NA_LENGTH] = 4.0;

            Node<2> node_b(0u, std::vector<double>{4.0, 0.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 4.0;

            const double d = force.CalculateDistanceBetweenCapsules(node_a, node_b);

            double contact_dist_a;
            double contact_dist_b;
            c_vector<double, 2> vec_a_to_b;

            force.CalculateForceDirectionAndContactPoints(node_a, node_b, d, vec_a_to_b, contact_dist_a, contact_dist_b);

            TS_ASSERT_DELTA(contact_dist_a, 0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);

            // Swap a->b for coverage
            force.CalculateForceDirectionAndContactPoints(node_b, node_a, d, vec_a_to_b, contact_dist_b, contact_dist_a);

            TS_ASSERT_DELTA(contact_dist_a, 0.5 * attributes_a[NA_LENGTH], 1e-6);
            TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);
        }
    }

    void TestCalculateForceMagnitude() throw(Exception)
    {
        CapsuleForce<2, 2> force;

        // Hand calculate a trivial case
        {
            double overlap = 1.0;
            double radiusA = 4.0;
            double radiusB = 4.0;

            TS_ASSERT_DELTA(force.CalculateForceMagnitude(overlap, radiusA, radiusB), 800.0 / 3.0, 1e-6);
        }

        // Hand calculate a nontrivial case
        {
            double overlap = 1.23;
            double radiusA = 2.34;
            double radiusB = 3.45;

            TS_ASSERT_DELTA(force.CalculateForceMagnitude(overlap, radiusA, radiusB), 303.731332875, 1e-6);
        }
    }

    void TestAddForceContribution() throw(Exception)
    {
        // Create two nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u,  false,  0.0, 0.0));
        nodes.push_back(new Node<2>(1u,  false,  2.0, 0.0));

        // Create mesh with massive interaction distance so all nodes interact with each other
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1e6);

        mesh.GetNode(0u)->AddNodeAttribute(0.0);
        mesh.GetNode(0u)->ClearAppliedForce();
        mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.25;
        mesh.GetNode(0u)->rGetNodeAttributes()[NA_APPLIED_THETA] = 0.0;

        mesh.GetNode(1u)->AddNodeAttribute(0.0);
        mesh.GetNode(1u)->ClearAppliedForce();
        mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.25;
        mesh.GetNode(1u)->rGetNodeAttributes()[NA_APPLIED_THETA] = 0.0;

        //Create cells
        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        // Create cell population
        NodeBasedCellPopulation<2> population(mesh, cells);

        MARK;
        CapsuleForce<2, 2> force;
        MARK;
        force.AddForceContribution(population);
        MARK;

        // Nodes 0 and 1 are too far apart to interact with each other, so no applied force or angle
        {
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetNodeAttributes()[NA_APPLIED_THETA], 0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetNodeAttributes()[NA_APPLIED_THETA], 0.0, 1e-6);
        }

        MARK;
    }
};

#endif /*_TESTCAPSULEFORCE_HPP_*/
