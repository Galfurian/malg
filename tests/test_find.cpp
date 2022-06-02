#include "malg/utility.hpp"
#include "malg/matrix.hpp"
#include "malg/io.hpp"

int main(int, char *[])
{
    malg::Matrix<int> a{
        { 10, 1, 2, 6, -3 },
        { -4, 4, -1, -4, 7 },
        { -4, -2, 0, 9, -9 },
        { -4, -8, -4, -1, 7 },
        { -4, -9, -6, -6, 8 },
        { -5, 0, -1, 8, 10 }
    };
    auto indices = malg::utility::find(a < -5);
#ifdef ROW_MAJOR
    malg::Matrix<std::size_t> control{ { 14},{16},{21},{22},{23 } };
#else
    malg::Matrix<std::size_t> control{ { 9},{10},{16},{22},{26 } };
#endif
    // Check the operations.
    if (!malg::all(indices == control))
        return 1;
    return 0;
}
