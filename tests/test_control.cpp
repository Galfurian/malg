#include "malg/control/control.hpp"
#include "malg/io.hpp"

int test_c2d()
{
    malg::Matrix<double> A = {
        { -5.48328, 0.96393, -7.51204, 5.44579, -7.60848 },
        { 1.34257, -8.9994, -9.63156, 5.80513, 4.864 },
        { 8.5414, -0.685415, -4.95847, 7.58094, -0.0461616 },
        { 8.53525, 0.564686, 6.03278, -0.992654, 3.41276 },
        { 4.14144, 2.8026, -8.80536, -4.99359, -6.67707 }
    };
    malg::Matrix<double> B = {
        { 4.14299 },
        { 4.56572 },
        { 5.34455 },
        { 4.78882 },
        { -9.34907 }
    };
    malg::Matrix<double> C = {
        { -6.75743, 1.61057, -8.59706, -3.852, -5.76971 },
        { 5.33681, -2.14648, -1.37983, -0.30633, 4.17382 },
        { -1.38097, 5.1317, -2.39843, 5.62033, -2.95562 },
        { -8.32778, 2.35308, -2.10552, -3.82267, 5.30287 },
        { -7.79774, 1.97377, -5.66731, 0.177432, 0.132991 }
    };
    malg::Matrix<double> D = {
        { 0 },
        { 0 },
        { 0 },
        { 0 },
        { 0 }
    };
    // Reference.
    malg::Matrix<double> ref_A = {
        { 0.994508, 0.000950468, -0.00742783, 0.00542156, -0.00755056 },
        { 0.00132652, 0.991053, -0.00957357, 0.00573143, 0.00483101 },
        { 0.00852862, -0.000674497, 0.995048, 0.00757982, -6.70223e-05 },
        { 0.00854067, 0.000568664, 0.00596532, 0.999047, 0.00336857 },
        { 0.00405941, 0.00278427, -0.008798, -0.00498858, 0.993328 }
    };
    malg::Matrix<double> ref_B = {
        { 0.00416226 },
        { 0.00451353 },
        { 0.00536587 },
        { 0.00480567 },
        { -0.0093385 }
    };
    // Create the state-space model.
    malg::control::StateSpace<double> sys{ A, B, C, D };
    // Perform the discretization.
    auto dsys = malg::control::c2d(sys, 1e-03);
    // Set the tollerance.
    malg::feq::tolerance() = 1e-04;
    // Check the operations.
    if (!malg::all(dsys.A == ref_A)) {
        return 1;
    }
    if (!malg::all(dsys.B == ref_B)) {
        return 1;
    }
    return 0;
}

int test_acker()
{
    malg::Matrix<double> A = {
        { -5.48328, 0.96393, -7.51204, 5.44579, -7.60848 },
        { 1.34257, -8.9994, -9.63156, 5.80513, 4.864 },
        { 8.5414, -0.685415, -4.95847, 7.58094, -0.0461616 },
        { 8.53525, 0.564686, 6.03278, -0.992654, 3.41276 },
        { 4.14144, 2.8026, -8.80536, -4.99359, -6.67707 }
    };
    malg::Matrix<double> B = {
        { 4.14299 },
        { 4.56572 },
        { 5.34455 },
        { 4.78882 },
        { -9.34907 }
    };
    malg::Matrix<double> C = {
        { -6.75743, 1.61057, -8.59706, -3.852, -5.76971 },
        { 5.33681, -2.14648, -1.37983, -0.30633, 4.17382 },
        { -1.38097, 5.1317, -2.39843, 5.62033, -2.95562 },
        { -8.32778, 2.35308, -2.10552, -3.82267, 5.30287 },
        { -7.79774, 1.97377, -5.66731, 0.177432, 0.132991 }
    };
    malg::Matrix<double> D = {
        { 0 },
        { 0 },
        { 0 },
        { 0 },
        { 0 }
    };
    // Reference.
    malg::Matrix<double> ref_K = { { -3.69417, -5.54803, 21.2759, -6.09239, -1.67523 } };
    // Create the state-space model.
    malg::control::StateSpace<double> sys{ A, B, C, D };
    // Select the poles.
    malg::Matrix<double> poles = { { -19.8672, -18.54, -13.3962, -19.5758, -15.2926 } };
    // Apply ackerman and compute the gain.
    auto K = malg::control::acker(sys.A, sys.B, poles);
    // Set the tollerance.
    malg::feq::tolerance() = 1e-04;
    // Check the operations.
    if (!malg::all(K == ref_K)) {
        {
            return 1;
        }
    }
    return 0;
}

int test_polyreduce()
{
    // Input.
    malg::Vector<int> a = { 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 };
    malg::Vector<int> b = { 0, 1, 1, 1, 1, 0, 1, 0, 0, 0 };
    malg::Vector<int> c = { 0, 1, 1, 0, 1, 1, 0, 1, 0, 0 };
    malg::Vector<int> d = { 0, 0, 0, 1, 0, 1, 0, 1, 1, 1 };
    malg::Vector<int> e = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };
    // Reference.
    malg::Vector<int> ref_a = { 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 };
    malg::Vector<int> ref_b = { 1, 1, 1, 1, 0, 1, 0, 0, 0 };
    malg::Vector<int> ref_c = { 1, 1, 0, 1, 1, 0, 1, 0, 0 };
    malg::Vector<int> ref_d = { 1, 0, 1, 0, 1, 1, 1 };
    malg::Vector<int> ref_e = { 1, 0 };
    // Check the operations.
    if (!malg::all(malg::control::polyreduce(a) == ref_a)) {
        return 1;
    }
    if (!malg::all(malg::control::polyreduce(b) == ref_b)) {
        return 1;
    }
    if (!malg::all(malg::control::polyreduce(c) == ref_c)) {
        return 1;
    }
    if (!malg::all(malg::control::polyreduce(d) == ref_d)) {
        return 1;
    }
    if (!malg::all(malg::control::polyreduce(e) == ref_e)) {
        return 1;
    }
    return 0;
}

int test_poly()
{
    // Input.
    malg::Vector<double> a = { -4.31365916, -5.43236374, 8.30901341, -4.1052015, -8.75671849, -8.20968452, -7.9560314, 4.93680757, 9.39615916, -4.92315493 };
    malg::Vector<double> b = { -6.59672341, -0.663223138, 2.65719019, -6.36291406, 9.07896614, 0.975672407, 6.18531446, -5.29307772, -5.14025403, -1.46744895 };
    malg::Vector<double> c = { 2.44049723, -9.49302124, -8.37147845, 1.51190319, -7.40282487, 8.66589357, -5.72460565, 2.24800701, -4.27068812, 6.31589424 };
    malg::Vector<double> d = { 1.18684296, -8.83575517, 9.88988867, 2.2718604, 7.85940986, 2.26749284, -5.22560682, 2.30190354, -4.51284885, 3.30124804 };
    malg::Vector<double> e = { 7.1682439, -8.06318395, -0.528348051, 2.64482488, 4.27461831, 7.59099286, 8.16414622, -5.98537147, 3.02603983, 1.75029028 };
    // Reference.
    malg::Vector<double> ref_a = { 1, 21.0548336, -17.4897007, -3268.17296, -19192.7843, 107667.157, 1353231.96, 2417940.15, -17396408.5, -82999934.9, -104405414 };
    malg::Vector<double> ref_b = { 1, 6.62649812, -112.911725, -902.680926, 2500.57233, 30199.3369, 19059.1166, -210355.545, -212728.833, 195535.603, 161817.207 };
    malg::Vector<double> ref_c = { 1, 14.0804231, -98.0180687, -1916.89923, 887.608435, 78625.4753, 66775.9288, -1145145.79, -234697.692, 6666982.31, -6529736.95 };
    malg::Vector<double> ref_d = { 1, -10.5044355, -101.444433, 1251.43984, 1027.49399, -37680.6611, 74338.4234, 262291.876, -1195069.98, 1627862.7, -752486.363 };
    malg::Vector<double> ref_e = { 1, -20.0422528, 43.7125343, 1560.5739, -11639.3617, -302.944468, 275213.098, -1066030.37, 1436090.15, -187362.574, -678295.446 };
    // Set the tollerance, specifically for div.
    malg::feq::tolerance() = 1e-06;
    // Check the operations.
    if (!malg::all(malg::control::poly(a) == ref_a)) {
        return 1;
    }
    if (!malg::all(malg::control::poly(b) == ref_b)) {
        return 1;
    }
    if (!malg::all(malg::control::poly(c) == ref_c)) {
        return 1;
    }
    if (!malg::all(malg::control::poly(d) == ref_d)) {
        return 1;
    }
    if (!malg::all(malg::control::poly(e) == ref_e)) {
        return 1;
    }
    return 0;
}

int main(int, char *[])
{
    if (test_c2d()) {
        return 1;
    }
    if (test_acker()) {
        return 1;
    }
    if (test_polyreduce()) {
        return 1;
    }
    if (test_poly()) {
        return 1;
    }
    return 0;
}
