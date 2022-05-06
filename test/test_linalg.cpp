#include "malg/matrix.hpp"
#include "malg/linalg.hpp"
#include "malg/math.hpp"
#include "malg/io.hpp"

int test_expm()
{
    // Input.
    malg::Matrix<double> a = {
        { -15.5458944, -97.1885676, -2.25545699, 70.1863391, -35.2873173 },
        { 66.4465788, 43.8484354, 12.583171, -79.0772638, 95.624119 },
        { 99.8136885, 59.3016462, 84.9035347, -23.2639528, -50.9771463 },
        { 91.4942282, -22.5984618, 17.8706132, 44.7630693, -9.99790168 },
        { 39.0343824, -19.6871119, 98.6161879, -60.2191987, 98.6344944 }
    };
    malg::Matrix<double> b = {
        { 63.3507647, -39.7977986, -55.8888088, -61.3599192, -64.8363302 },
        { -36.84952, -71.5332789, 46.5774608, -74.3500775, 70.8524849 },
        { -21.2086081, 65.8513081, -53.2621921, 54.9914801, -7.02221337 },
        { 18.1242303, -83.9546014, -35.9837116, 32.0856732, -77.3596723 },
        { -69.7916162, 90.902458, -2.02615437, 70.279005, 85.7506663 }
    };
    malg::Matrix<double> c = {
        { -54.1910383, 67.3419784, -4.97696063, 7.51835372, 70.7248749 },
        { -83.3971181, -90.3452801, -89.8086915, 73.2110995, -4.93238127 },
        { -39.2404358, 83.4662255, 86.2141902, 82.8031647, -82.6533288 },
        { -25.5733975, -50.9251651, 70.8319844, 74.3622488, -8.98154252 },
        { 29.7198651, -47.0061289, 63.3928138, 20.0981107, -56.0546803 }
    };
    malg::Matrix<double> d = {
        { 65.1809585, -42.8504736, 68.1575291, -25.2619994, -71.1127418 },
        { -79.749854, -63.0507856, 58.6389908, 27.6389739, -5.2120999 },
        { -0.49745321, 60.8086088, -53.7952418, -54.9404498, 62.5375388 },
        { 55.0985219, -14.0161005, 96.1509587, -17.6416804, 54.0363257 },
        { 64.5736433, 6.62711579, 12.3644458, -68.3727389, 2.99360322 }
    };
    malg::Matrix<double> e = {
        { -52.1272295, 88.780273, 78.3526637, -13.276907, 34.6016492 },
        { -64.1692363, 31.8109079, -7.50080231, -66.57565, -88.2098535 },
        { 18.360975, -68.0999633, -66.1906277, 57.7654585, 1.02118668 },
        { -84.1868939, -99.0993309, -92.1614526, 66.0547531, -84.7551219 },
        { -2.78259091, -75.5274384, 74.9108197, 20.548016, -24.8396294 }
    };
    // Reference.
    malg::Matrix<double> ref_a = {
        { -6.8641161e+38, -5.28956323e+38, -1.17870565e+39, 8.60190225e+38, -1.64630967e+37 },
        { 1.99255725e+38, 7.7264104e+37, 4.6416759e+37, 8.30193829e+37, -3.95544728e+38 },
        { -6.20626954e+38, -5.01862772e+38, -1.20254272e+39, 9.16562108e+38, -1.98642391e+38 },
        { -9.92990851e+38, -8.77721929e+38, -2.10767045e+39, 1.70844095e+39, -5.69710177e+38 },
        { -3.66179457e+38, -4.0415798e+38, -1.13267554e+39, 1.01542957e+39, -6.89858307e+38 }
    };
    malg::Matrix<double> ref_b = {
        { -2.74079999e+52, 5.73787112e+50, 3.39221923e+51, 3.83171483e+52, 1.26271261e+52 },
        { 5.01839137e+52, -2.03093079e+52, -1.25544954e+52, -4.43370007e+52, -4.64752456e+52 },
        { -2.01144126e+51, 3.35846821e+51, 1.34128045e+51, -1.63442568e+51, 4.94843928e+51 },
        { -7.82244008e+52, 3.57085589e+52, 2.09037888e+52, 6.36785887e+52, 7.73565517e+52 },
        { 4.92769625e+52, -1.13857557e+52, -9.50929294e+51, -5.50080382e+52, -3.52588356e+52 }
    };
    malg::Matrix<double> ref_c = {
        { -2.26319519e+54, 2.8723809e+53, 4.14939661e+54, 5.34457492e+54, -2.99957483e+54 },
        { 1.02862586e+54, -1.30550175e+53, -1.88590745e+54, -2.42911792e+54, 1.3633116e+54 },
        { -1.15267357e+55, 1.46293945e+54, 2.1133395e+55, 2.72205874e+55, -1.52772092e+55 },
        { -1.40615333e+55, 1.78464853e+54, 2.57807542e+55, 3.32065563e+55, -1.86367582e+55 },
        { -6.12971162e+54, 7.77965005e+53, 1.1238361e+55, 1.44754209e+55, -8.12414626e+54 }
    };
    malg::Matrix<double> ref_d = {
        { -536745180, -150592812, -858208588, -108725257, 823064724 },
        { 1.72332308e+09, -111399009, 1.5766801e+09, 264337160, -2.04796303e+09 },
        { 1.13733153e+09, -63466680.6, 1.06047066e+09, 175884570, -1.36163042e+09 },
        { -74309435.1, -111544364, -298521411, -27972724.2, 204606147 },
        { -759985590, 86666024.2, -620935752, -111225953, 865629943 }
    };
    malg::Matrix<double> ref_e = {
        { 1.12298542e+44, 6.60836243e+44, 3.62974711e+44, -6.33814012e+44, -2.44414379e+42 },
        { 6.27078362e+44, 3.69012903e+45, 2.02686147e+45, -3.5392361e+45, -1.36481708e+43 },
        { -2.62626957e+44, -1.54546452e+45, -8.48870718e+44, 1.48226899e+45, 5.71599627e+42 },
        { -6.45505798e+43, -3.79856781e+44, -2.08642318e+44, 3.64324072e+44, 1.40492384e+42 },
        { -5.29652557e+44, -3.11681346e+45, -1.71195886e+45, 2.98936395e+45, 1.15277277e+43 }
    };
    // Set the tollerance, specifically for div.
    feq::tolerance() = 1e-06;
    // Check the operations.
    if (!malg::all(malg::linalg::expm(a) == ref_a))
        return 1;
    if (!malg::all(malg::linalg::expm(b) == ref_b))
        return 1;
    if (!malg::all(malg::linalg::expm(c) == ref_c))
        return 1;
    if (!malg::all(malg::linalg::expm(d) == ref_d))
        return 1;
    if (!malg::all(malg::linalg::expm(e) == ref_e))
        return 1;
    return 0;
}

int main(int, char *[])
{
    test_expm();
    return 0;
}
