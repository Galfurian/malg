#include "malg/control/control.hpp"
#include "malg/io.hpp"

int test_view()
{
    malg::Matrix<int> a{
        { 10, 1, 2, 6, -3 },
        { -4, 4, -1, -4, 7 },
        { -4, -2, 0, 9, -9 },
        { -4, -8, -4, -1, 7 },
        { -4, -9, -6, -6, 8 },
        { -5, 0, -1, 8, 10 }
    };
    // References.
    malg::Matrix<int> ref1 = { { 1, 2, 6 }, { 4, -1, -4 }, { -2, 0, 9 } };
    malg::Matrix<int> ref2 = { { 10, 1 }, { -4, 4 }, { -4, -2 }, { -4, -8 }, { -4, -9 } };
    malg::Matrix<int> ref3 = { { 2, 6, -3 }, { -1, -4, 7 }, { 0, 9, -9 }, { -4, -1, 7 }, { -6, -6, 8 }, { -1, 8, 10 } };
    // Views.
    auto view1 = a({ 0, 3 }, { 1, 4 });
    auto view2 = a({ 0, 5 }, { 0, 2 });
    auto view3 = a(malg::Range::all(), { 2, 5 });
    // Check the operations.
    if (!malg::all(ref1 == view1))
        return 1;
    if (!malg::all(ref2 == view2))
        return 1;
    if (!malg::all(ref3 == view3))
        return 1;
    return 0;
}

int main(int, char *[])
{
    if (test_view())
        return 1;
    return 0;
}
