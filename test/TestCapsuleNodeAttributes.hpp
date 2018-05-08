
#ifndef _TESTCAPSULENODEATTRIBUTES_HPP_
#define _TESTCAPSULENODEATTRIBUTES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "OutputFileHandler.hpp"
#include "Node.hpp"
#include "NodeAttributes.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCapsuleNodeAttributes : public CxxTest::TestSuite
{
public:

    void TestAttributes2d() throw(Exception)
    {
        // index, {x, y}
        Node<2> node(0u, std::vector<double>{0.0, 0.0});

        // Create node attributes
        node.AddNodeAttribute(0.0);

        std::vector<double>& attributes = node.rGetNodeAttributes();

        attributes.resize(NA_VEC_LENGTH);
        attributes[NA_THETA] = 1.23;
        attributes[NA_LENGTH] = 2.34;
        attributes[NA_RADIUS] = 3.45;

        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_THETA], 1.23, 1e-6);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_LENGTH], 2.34, 1e-6);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_RADIUS], 3.45, 1e-6);
    }

    void TestAttributes3d() throw(Exception)
    {
        // index, {x, y, z}
        Node<3> node(0u, std::vector<double>{0.0, 0.0, 0.0});

        // Create node attributes
        node.AddNodeAttribute(0.0);

        std::vector<double>& attributes = node.rGetNodeAttributes();

        attributes.resize(NA_VEC_LENGTH);
        attributes[NA_THETA] = 1.23;
        attributes[NA_PHI] = 0.34;
        attributes[NA_LENGTH] = 2.34;
        attributes[NA_RADIUS] = 3.45;

        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_THETA], 1.23, 1e-6);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_PHI], 0.34, 1e-6);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_LENGTH], 2.34, 1e-6);
        TS_ASSERT_DELTA(node.rGetNodeAttributes()[NA_RADIUS], 3.45, 1e-6);
    }
};

#endif /*_TESTCAPSULENODEATTRIBUTES_HPP_*/
