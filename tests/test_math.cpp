#include "malg/matrix.hpp"
#include "malg/io.hpp"
#include "malg/linalg.hpp"
#include "malg/math.hpp"

int test_int()
{
    // Input.
    const malg::Matrix<int> a = {
        { -77, 43, 81, -60, -33 },
        { -95, -40, 14, 41, 13 },
        { 0, 41, 93, -89, -17 },
        { -60, 1, -74, -53, 69 },
        { -9, 10, -89, 30, -72 }
    };
    const malg::Matrix<int> b = {
        { 67, 38, 34, 88, 28 },
        { 0, -12, -12, 42, -69 },
        { 20, -54, 15, -27, -20 },
        { 77, 16, 54, -18, 18 },
        { -8, 93, 88, -40, -9 }
    };
    // Reference.
    const malg::Matrix<int> sum = {
        { -10, 81, 115, 28, -5 },
        { -95, -52, 2, 83, -56 },
        { 20, -13, 108, -116, -37 },
        { 17, 17, -20, -71, 87 },
        { -17, 103, -1, -10, -81 }
    };
    const malg::Matrix<int> sub = {
        { -144, 5, 47, -148, -61 },
        { -95, -28, 26, -1, 82 },
        { -20, 95, 78, -62, 3 },
        { -137, -15, -128, -35, 51 },
        { -1, -83, -177, 70, -63 }
    };
    const malg::Matrix<int> mul = {
        { -7895, -11845, -8063, -4757, -7526 },
        { -3032, -2021, 818, -11676, 441 },
        { -4857, -8519, -5399, 1493, -6138 },
        { -10133, 7273, 48, -5046, -1844 },
        { 503, -1872, -6477, 4371, 2026 }
    };
    const malg::Matrix<double> div = {
        { -0.724101333, 0.959172965, -3.50660292, 3.52737442, -3.30035297 },
        { -3.87515054, 3.01491411, -6.44053984, 7.70150764, -2.09035033 },
        { 0.609014897, 0.368501698, -1.70527853, -0.311900352, -2.16473747 },
        { 1.87805156, 1.93870341, 2.98797968, -0.462939925, 1.27233358 },
        { 2.98660079, -4.03525548, 5.63094744, -3.59205073, 1.89608685 }
    };
    // Set the tollerance, specifically for div.
    malg::feq::tolerance() = 1e-06;
    // Check the operations.
    if (!malg::all((a + b) == sum)) {
        return 1;
    }
    if (!malg::all((a - b) == sub)) {
        return 1;
    }
    if (!malg::all((a * b) == mul)) {
        return 1;
    }
    if (!malg::all(malg::linalg::div(a, b) == div)) {
        return 1;
    }
    return 0;
}

int test_double()
{
    // Input.
    const malg::Matrix<double> a = {
        { -37.2013365, -68.6581972, 41.3047901, 39.1891392, -83.7775724 },
        { -29.4343112, 19.2121221, -92.9539817, 79.2170016, 87.7778305 },
        { -17.6058458, -4.69375854, 22.4131617, 67.8548823, -72.9728209 },
        { 82.8826396, -23.414967, -63.2745926, 99.0710991, -39.0158564 },
        { -32.7507545, -20.9186668, -0.269948579, 27.4472885, -97.0661874 }
    };
    const malg::Matrix<double> b = {
        { -37.4315122, -46.8778636, 98.621834, 64.875136, -82.3723058 },
        { -19.5194821, -29.6733099, 27.642446, -66.8425277, 18.3564081 },
        { -56.4430537, -85.474696, 44.795422, 44.2905858, -2.45158714 },
        { 27.0582983, -64.0869579, 11.8315765, 58.6659616, -95.9791052 },
        { -33.592816, 37.7442705, 6.9947062, 90.3429215, 9.45861974 }
    };
    // Reference.
    const malg::Matrix<double> sum = {
        { -74.6328488, -115.536061, 139.926624, 104.064275, -166.149878 },
        { -48.9537934, -10.4611878, -65.3115357, 12.374474, 106.134239 },
        { -74.0488995, -90.1684545, 67.2085837, 112.145468, -75.424408 },
        { 109.940938, -87.501925, -51.4430161, 157.737061, -134.994962 },
        { -66.3435705, 16.8256037, 6.72475762, 117.79021, -87.6075677 }
    };
    const malg::Matrix<double> sub = {
        { 0.230175709, -21.7803337, -57.3170439, -25.6859969, -1.40526653 },
        { -9.91482914, 48.885432, -120.596428, 146.059529, 69.4214224 },
        { 38.8372079, 80.7809375, -22.3822602, 23.5642965, -70.5212337 },
        { 55.8243412, 40.6719909, -75.1061691, 40.4051375, 56.9632487 },
        { 0.84206149, -58.6629372, -7.26465478, -62.8956329, -106.524807 }
    };
    const malg::Matrix<double> mul = {
        { 4276.02224, -5422.91529, -3838.80925, -1264.38308, -2850.97902 },
        { 5168.13944, 6991.27723, -4984.46639, 5266.71917, -3767.79738 },
        { 3772.96607, -8054.07327, -569.653942, -2447.5528, -5893.7473 },
        { 4917.38457, -5603.97532, 5591.64172, 6426.96416, -16979.7172 },
        { 5652.87203, -3243.6204, -4174.47982, -7897.43065, -1238.05296 }
    };
    const malg::Matrix<double> div = {
        { -1.75792991, 7.4526154, -1.77961696, 4.36930452, 5.24532028 },
        { -1.84669495, 1.44407952, 0.432614559, 1.16256939, 2.3043699 },
        { -1.14684653, 5.19297263, -1.67231909, 3.19086747, 4.16458535 },
        { 1.58260119, -9.06219024, 2.47852257, -3.32018408, -5.80381381 },
        { -4.05907923, 15.0937256, -4.35863652, 8.57273315, 10.9560582 }
    };
    // Set the tollerance, specifically for div.
    malg::feq::tolerance() = 1e-06;
    // Check the operations.
    if (!malg::all((a + b) == sum)) {
        return 1;
    }
    if (!malg::all((a - b) == sub)) {
        return 1;
    }
    if (!malg::all((a * b) == mul)) {
        return 1;
    }
    if (!malg::all(malg::linalg::div(a, b) == div)) {
        return 1;
    }
    return 0;
}

int test_complex()
{
    using namespace std::complex_literals;
    // Input.
    const malg::Matrix<std::complex<double>> a = {
        { std::complex<double>(6.21972047, 55.6003688), std::complex<double>(-62.0717653, -9.29039483), std::complex<double>(-79.2333101, +25.5526935), std::complex<double>(70.1014075, -15.4969014), std::complex<double>(32.7315868, +79.1097553) },
        { std::complex<double>(7.51023105, +3.64680008), std::complex<double>(63.7374363, +43.92121), std::complex<double>(36.9210945, -56.4619748), std::complex<double>(-20.2259046, +73.9420527), std::complex<double>(48.0314307, +22.4845855) },
        { std::complex<double>(10.9792993, +19.9594746), std::complex<double>(79.8537074, +55.1074495), std::complex<double>(83.2460263, -63.3195366), std::complex<double>(-58.0035035, +95.0898121), std::complex<double>(98.1701286, -18.5235781) },
        { std::complex<double>(-46.0580658, -92.6997876), std::complex<double>(-97.3461191, -8.6252516), std::complex<double>(87.8804822, -57.8865619), std::complex<double>(3.50403728, +47.7014354), std::complex<double>(-15.613705, +24.7798349) },
        { std::complex<double>(67.6734711, +56.103517), std::complex<double>(-26.3088772, +48.460602), std::complex<double>(70.6932588, -79.7588588), std::complex<double>(-85.9500325, +35.6230139), std::complex<double>(-11.8113297, +78.4160521) }
    };
    const malg::Matrix<std::complex<double>> b = {
        { std::complex<double>(37.6463194, +55.5178134), std::complex<double>(-18.4074263, -65.6354342), std::complex<double>(11.291236, +19.5497491), std::complex<double>(26.5698575, -96.6589965), std::complex<double>(-31.1807715, -81.1488153) },
        { std::complex<double>(27.0583953, +15.5557727), std::complex<double>(41.347285, -49.2771551), std::complex<double>(-86.1100143, -8.4328456), std::complex<double>(-91.8965878, +97.7787111), std::complex<double>(-46.1050965, +47.456248) },
        { std::complex<double>(25.6042146, +53.7801993), std::complex<double>(-40.6112622, +9.19292961), std::complex<double>(2.85852944, +2.69556217), std::complex<double>(97.6854439, -0.892334903), std::complex<double>(12.801631, -12.0367261) },
        { std::complex<double>(-45.3993398, +85.3966666), std::complex<double>(74.6956836, -29.4508223), std::complex<double>(51.8387618, -51.7090829), std::complex<double>(81.4358985, -69.3236742), std::complex<double>(58.2171779, -7.16924139) },
        { std::complex<double>(-69.3605792, +48.1518686), std::complex<double>(39.8427618, +51.2287156), std::complex<double>(-27.0086001, +30.983493), std::complex<double>(90.1640195, -84.6278121), std::complex<double>(-59.2749554, -68.5078407) }
    };
    // Reference.
    const malg::Matrix<std::complex<double>> sum = {
        { std::complex<double>(43.8660398, +111.118182), std::complex<double>(-80.4791916, -74.9258291), std::complex<double>(-67.942074, +45.1024427), std::complex<double>(96.671265, -112.155898), std::complex<double>(1.55081528, -2.03905996) },
        { std::complex<double>(34.5686264, +19.2025728), std::complex<double>(105.084721, -5.3559451), std::complex<double>(-49.1889198, -64.8948204), std::complex<double>(-112.122492, +171.720764), std::complex<double>(1.92633415, +69.9408335) },
        { std::complex<double>(36.5835139, +73.7396739), std::complex<double>(39.2424452, +64.3003792), std::complex<double>(86.1045557, -60.6239744), std::complex<double>(39.6819404, +94.1974772), std::complex<double>(110.97176, -30.5603042) },
        { std::complex<double>(-91.4574056, -7.30312099), std::complex<double>(-22.6504355, -38.0760739), std::complex<double>(139.719244, -109.595645), std::complex<double>(84.9399358, -21.6222388), std::complex<double>(42.6034729, +17.6105935) },
        { std::complex<double>(-1.68710806, +104.255386), std::complex<double>(13.5338845, +99.6893176), std::complex<double>(43.6846587, -48.7753658), std::complex<double>(4.21398695, -49.0047983), std::complex<double>(-71.0862851, +9.90821134) }
    };
    const malg::Matrix<std::complex<double>> sub = {
        { std::complex<double>(-31.4265989, +0.0825553932), std::complex<double>(-43.6643391, +56.3450394), std::complex<double>(-90.5245461, +6.00294437), std::complex<double>(43.5315499, +81.1620951), std::complex<double>(63.9123582, +160.258571) },
        { std::complex<double>(-19.5481643, -11.9089726), std::complex<double>(22.3901514, +93.1983651), std::complex<double>(123.031109, -48.0291292), std::complex<double>(71.6706833, -23.8366585), std::complex<double>(94.1365272, -24.9716626) },
        { std::complex<double>(-14.6249152, -33.8207247), std::complex<double>(120.46497, +45.9145199), std::complex<double>(80.3874968, -66.0150988), std::complex<double>(-155.688947, +95.982147), std::complex<double>(85.3684976, -6.48685207) },
        { std::complex<double>(-0.658726072, -178.096454), std::complex<double>(-172.041803, +20.8255707), std::complex<double>(36.0417204, -6.177479), std::complex<double>(-77.9318612, +117.02511), std::complex<double>(-73.830883, +31.9490763) },
        { std::complex<double>(137.03405, +7.95164839), std::complex<double>(-66.151639, -2.76811362), std::complex<double>(97.7018589, -110.742352), std::complex<double>(-176.114052, +120.250826), std::complex<double>(47.4636257, +146.923893) }
    };
    const malg::Matrix<std::complex<double>> mul = {
        { std::complex<double>(-15729.3783, +393.534741), std::complex<double>(5524.72444, +1083.43283), std::complex<double>(3452.0503, -3618.4312), std::complex<double>(18715.5456, -3531.4511), std::complex<double>(14363.4158, -11811.2754) },
        { std::complex<double>(-4706.7783, -1056.76389), std::complex<double>(5349.15112, +10222.821), std::complex<double>(-4065.73103, +1566.41267), std::complex<double>(3668.84828, +1404.63678), std::complex<double>(-7122.24205, -1064.22327) },
        { std::complex<double>(-5258.73983, +3691.5543), std::complex<double>(7654.7563, +13694.4559), std::complex<double>(-6436.41659, +6535.30201), std::complex<double>(6721.99659, -2260.14522), std::complex<double>(-14499.1048, -1752.82647) },
        { std::complex<double>(1933.13282, -8887.39193), std::complex<double>(-12948.3454, +15975.7617), std::complex<double>(12311.3721, +826.480396), std::complex<double>(12419.5692, -5272.63418), std::complex<double>(2408.46411, +2959.91843) },
        { std::complex<double>(1970.03643, -6433.90808), std::complex<double>(-8259.55967, +9426.04534), std::complex<double>(-1965.60586, +1775.14428), std::complex<double>(12776.1178, -3001.52339), std::complex<double>(2624.64344, -13744.6046) }
    };
    const malg::Matrix<std::complex<double>> div = {
        { std::complex<double>(-2.61040517, -0.51938165), std::complex<double>(-0.86673903, -1.64601849), std::complex<double>(0.103491678, +0.796748593), std::complex<double>(0.314309467, -2.28212603), std::complex<double>(1.47281764, -0.731910495) },
        { std::complex<double>(0.546504938, +0.924958461), std::complex<double>(0.502346129, +1.28819315), std::complex<double>(0.94023401, -0.537185686), std::complex<double>(-0.280741161, +1.18529653), std::complex<double>(-1.2076058, +0.0729997121) },
        { std::complex<double>(1.15521647, +2.4762995), std::complex<double>(0.0113777518, +2.59992251), std::complex<double>(1.62600223, -0.265827624), std::complex<double>(-1.43559827, +1.79784449), std::complex<double>(-2.2968276, -0.378406349) },
        { std::complex<double>(0.477587767, -2.40854546), std::complex<double>(0.179103912, -0.638206495), std::complex<double>(-0.891307585, +0.400488878), std::complex<double>(1.53100862, +0.408404799), std::complex<double>(0.0222926676, +1.29257084) },
        { std::complex<double>(-0.223354037, +0.0842494207), std::complex<double>(-0.0380571152, +1.29062814), std::complex<double>(1.74043017, +0.253183147), std::complex<double>(-0.360988464, +0.553631548), std::complex<double>(-1.07908079, -0.559970571) }
    };

    // Set the tollerance, specifically for div.
    malg::feq::tolerance() = 1e-06;
    // Check the operations.
    if (!malg::all((a + b) == sum)) {
        return 1;
    }
    if (!malg::all((a - b) == sub)) {
        return 1;
    }
    if (!malg::all((a * b) == mul)) {
        return 1;
    }
    if (!malg::all(malg::linalg::div(a, b) == div)) {
        return 1;
    }
    return 0;
}

int main(int, char *[])
{
    if (test_int()) {
        return 1;
    }
    if (test_double()) {
        return 1;
    }
    if (test_complex()) {
        return 1;
    }
    return 0;
}
