# mini-cpp-ryu

This repository is just a sandbox and will soon be deleted. My intention is just study the Ryu algoritm and do some experiments with the code of https://github.com/ulfjack/ryu, which I may later contibute to, or create a fork from.

This is C++11 header-only library. The header `ryu.hpp` contains the following function:


    enum class float_state
    {
        normal, nan, infinity
    }; 
    
    struct decimal_double
    {
        uint64_t mantissa;
        int32_t exponent;
        bool negative;
        float_state state;
    };
    
    inline decimal_double to_decimal(double f);

