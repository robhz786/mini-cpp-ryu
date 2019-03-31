#include "ryu/ryu.hpp"

#include <cstdio>

int main()
{
    auto v = ryu::to_decimal(1.12345e+8);

    std::printf(" %lu\n %d\n", v.mantissa, v.exponent);
    
    return 0;
}
