// Copyright 2018 Ulf Adams
//
// The contents of this file may be used under the terms of the Apache License,
// Version 2.0.
//
//    (See accompanying file LICENSE-Apache or copy at
//     http://www.apache.org/licenses/LICENSE-2.0)
//
// Alternatively, the contents of this file may be used under the terms of
// the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE-Boost or copy at
//     https://www.boost.org/LICENSE_1_0.txt)
//
// Unless required by applicable law or agreed to in writing, this software
// is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.
#ifndef RYU_D2S_HPP
#define RYU_D2S_HPP

#include <cassert>
#include <cstdint>

// ABSL avoids uint128_t on Win32 even if __SIZEOF_INT128__ is defined.
// Let's do the same for now.
#if defined(__SIZEOF_INT128__) && !defined(_MSC_VER) && !defined(RYU_ONLY_64_BIT_OPS)
#define HAS_UINT128
#elif defined(_MSC_VER) && !defined(RYU_ONLY_64_BIT_OPS) && defined(_M_X64)
#define HAS_64_BIT_INTRINSICS
#endif

#include "ryu/d2s_full_table.hpp"
#include "ryu/d2s_intrinsics.hpp"

namespace ryu {

// Returns e == 0 ? 1 : ceil(log_2(5^e)).
inline std::int32_t pow5bits(const std::int32_t e) {
  // This approximation works up to the point that the multiplication overflows at e = 3529.
  // If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
  // than 2^9297.
  assert(e >= 0);
  assert(e <= 3528);
  return (std::int32_t) (((((std::uint32_t) e) * 1217359) >> 19) + 1);
}

// Returns floor(log_10(2^e)).
inline std::uint32_t log10Pow2(const std::int32_t e) {
  // The first value this approximation fails for is 2^1651 which is just greater than 10^297.
  assert(e >= 0);
  assert(e <= 1650);
  return (((std::uint32_t) e) * 78913) >> 18;
}

// Returns floor(log_10(5^e)).
inline std::uint32_t log10Pow5(const std::int32_t e) {
  // The first value this approximation fails for is 5^2621 which is just greater than 10^1832.
  assert(e >= 0);
  assert(e <= 2620);
  return (((std::uint32_t) e) * 732923) >> 20;
}

inline std::uint64_t double_to_bits(const double d) {
  union { double v; std::uint64_t bits;};
  v = d;
  return bits;
}

inline std::uint32_t pow5Factor(std::uint64_t value) {
  std::uint32_t count = 0;
  for (;;) {
    assert(value != 0);
    const std::uint64_t q = div5(value);
    const std::uint32_t r = (std::uint32_t) (value - 5 * q);
    if (r != 0) {
      break;
    }
    value = q;
    ++count;
  }
  return count;
}

// Returns true if value is divisible by 5^p.
inline bool multipleOfPowerOf5(const std::uint64_t value, const std::uint32_t p) {
  // I tried a case distinction on p, but there was no performance difference.
  return pow5Factor(value) >= p;
}

// Returns true if value is divisible by 2^p.
inline bool multipleOfPowerOf2(const std::uint64_t value, const std::uint32_t p) {
  // return __builtin_ctzll(value) >= p;
  return (value & ((1ull << p) - 1)) == 0;
}


// We need a 64x128-bit multiplication and a subsequent 128-bit shift.
// Multiplication:
//   The 64-bit factor is variable and passed in, the 128-bit factor comes
//   from a lookup table. We know that the 64-bit factor only has 55
//   significant bits (i.e., the 9 topmost bits are zeros). The 128-bit
//   factor only has 124 significant bits (i.e., the 4 topmost bits are
//   zeros).
// Shift:
//   In principle, the multiplication result requires 55 + 124 = 179 bits to
//   represent. However, we then shift this value to the right by j, which is
//   at least j >= 115, so the result is guaranteed to fit into 179 - 115 = 64
//   bits. This means that we only need the topmost 64 significant bits of
//   the 64x128-bit multiplication.
//
// There are several ways to do this:
// 1. Best case: the compiler exposes a 128-bit type.
//    We perform two 64x64-bit multiplications, add the higher 64 bits of the
//    lower result to the higher result, and shift by j - 64 bits.
//
//    We explicitly cast from 64-bit to 128-bit, so the compiler can tell
//    that these are only 64-bit inputs, and can map these to the best
//    possible sequence of assembly instructions.
//    x64 machines happen to have matching assembly instructions for
//    64x64-bit multiplications and 128-bit shifts.
//
// 2. Second best case: the compiler exposes intrinsics for the x64 assembly
//    instructions mentioned in 1.
//
// 3. We only have 64x64 bit instructions that return the lower 64 bits of
//    the result, i.e., we have to use plain C.
//    Our inputs are less than the full width, so we have three options:
//    a. Ignore this fact and just implement the intrinsics manually.
//    b. Split both into 31-bit pieces, which guarantees no internal overflow,
//       but requires extra work upfront (unless we change the lookup table).
//    c. Split only the first factor into 31-bit pieces, which also guarantees
//       no internal overflow, but requires extra work since the intermediate
//       results are not perfectly aligned.
#if defined(HAS_UINT128)

// Best case: use 128-bit type.
inline std::uint64_t mulShift(const std::uint64_t m, const std::uint64_t* const mul, const std::int32_t j) {

  typedef __uint128_t uint128_t;

  const uint128_t b0 = ((uint128_t) m) * mul[0];
  const uint128_t b2 = ((uint128_t) m) * mul[1];
  return (std::uint64_t) (((b0 >> 64) + b2) >> (j - 64));
}

inline std::uint64_t mulShiftAll(const std::uint64_t m, const std::uint64_t* const mul, const std::int32_t j,
  std::uint64_t* const vp, std::uint64_t* const vm, const std::uint32_t mmShift) {
//  m <<= 2;
//  uint128_t b0 = ((uint128_t) m) * mul[0]; // 0
//  uint128_t b2 = ((uint128_t) m) * mul[1]; // 64
//
//  uint128_t hi = (b0 >> 64) + b2;
//  uint128_t lo = b0 & 0xffffffffffffffffull;
//  uint128_t factor = (((uint128_t) mul[1]) << 64) + mul[0];
//  uint128_t vpLo = lo + (factor << 1);
//  *vp = (std::uint64_t) ((hi + (vpLo >> 64)) >> (j - 64));
//  uint128_t vmLo = lo - (factor << mmShift);
//  *vm = (std::uint64_t) ((hi + (vmLo >> 64) - (((uint128_t) 1ull) << 64)) >> (j - 64));
//  return (std::uint64_t) (hi >> (j - 64));
  *vp = mulShift(4 * m + 2, mul, j);
  *vm = mulShift(4 * m - 1 - mmShift, mul, j);
  return mulShift(4 * m, mul, j);
}

#elif defined(HAS_64_BIT_INTRINSICS)

inline std::uint64_t mulShift(const std::uint64_t m, const std::uint64_t* const mul, const std::int32_t j) {
  // m is maximum 55 bits
  std::uint64_t high1;                                   // 128
  const std::uint64_t low1 = umul128(m, mul[1], &high1); // 64
  std::uint64_t high0;                                   // 64
  umul128(m, mul[0], &high0);                       // 0
  const std::uint64_t sum = high0 + low1;
  if (sum < high0) {
    ++high1; // overflow into high1
  }
  return shiftright128(sum, high1, j - 64);
}

inline std::uint64_t mulShiftAll(const std::uint64_t m, const std::uint64_t* const mul, const std::int32_t j,
  std::uint64_t* const vp, std::uint64_t* const vm, const std::uint32_t mmShift) {
  *vp = mulShift(4 * m + 2, mul, j);
  *vm = mulShift(4 * m - 1 - mmShift, mul, j);
  return mulShift(4 * m, mul, j);
}

#else // !defined(HAS_UINT128) && !defined(HAS_64_BIT_INTRINSICS)

inline std::uint64_t mulShiftAll(std::uint64_t m, const std::uint64_t* const mul, const std::int32_t j,
  std::uint64_t* const vp, std::uint64_t* const vm, const std::uint32_t mmShift) {
  m <<= 1;
  // m is maximum 55 bits
  std::uint64_t tmp;
  const std::uint64_t lo = umul128(m, mul[0], &tmp);
  std::uint64_t hi;
  const std::uint64_t mid = tmp + umul128(m, mul[1], &hi);
  hi += mid < tmp; // overflow into hi

  const std::uint64_t lo2 = lo + mul[0];
  const std::uint64_t mid2 = mid + mul[1] + (lo2 < lo);
  const std::uint64_t hi2 = hi + (mid2 < mid);
  *vp = shiftright128(mid2, hi2, (std::uint32_t) (j - 64 - 1));

  if (mmShift == 1) {
    const std::uint64_t lo3 = lo - mul[0];
    const std::uint64_t mid3 = mid - mul[1] - (lo3 > lo);
    const std::uint64_t hi3 = hi - (mid3 > mid);
    *vm = shiftright128(mid3, hi3, (std::uint32_t) (j - 64 - 1));
  } else {
    const std::uint64_t lo3 = lo + lo;
    const std::uint64_t mid3 = mid + mid + (lo3 < lo);
    const std::uint64_t hi3 = hi + hi + (mid3 < mid);
    const std::uint64_t lo4 = lo3 - mul[0];
    const std::uint64_t mid4 = mid3 - mul[1] - (lo4 > lo3);
    const std::uint64_t hi4 = hi3 - (mid4 > mid3);
    *vm = shiftright128(mid4, hi4, (std::uint32_t) (j - 64));
  }

  return shiftright128(mid, hi, (std::uint32_t) (j - 64 - 1));
}

#endif // HAS_64_BIT_INTRINSICS


struct d2d_step3_result
{
  std::uint64_t vr;
  std::uint64_t vp;
  std::uint64_t vm;
  const std::int32_t e10;
  bool vmIsTrailingZeros;
  bool vrIsTrailingZeros;
  const bool acceptBounds;
};


inline d2d_step3_result d2d_step_3
    ( std::int32_t e2
    , std::uint64_t m2
    , const std::uint64_t mv
    , const std::uint32_t mmShift
    , const bool acceptBounds )
{
  // Step 3: Convert to a decimal power base using 128-bit arithmetic.

  std::uint64_t vr, vp, vm;
  std::int32_t e10;
  bool vmIsTrailingZeros = false;
  bool vrIsTrailingZeros = false;

  if (e2 >= 0) {
    // I tried special-casing q == 0, but there was no effect on performance.
    // This expression is slightly faster than max(0, log10Pow2(e2) - 1).
    const std::uint32_t q = log10Pow2(e2) - (e2 > 3);
    e10 = (std::int32_t) q;
    constexpr std::int32_t double_pow5_inv_bitcount = 122;
    const std::int32_t k = double_pow5_inv_bitcount + pow5bits((std::int32_t) q) - 1;
    const std::int32_t i = -e2 + (std::int32_t) q + k;
    vr = mulShiftAll(m2, double_pow5_inv_split(q), i, &vp, &vm, mmShift);
    if (q <= 21) {
      // This should use q <= 22, but I think 21 is also safe. Smaller values
      // may still be safe, but it's more difficult to reason about them.
      // Only one of mp, mv, and mm can be a multiple of 5, if any.
      const std::uint32_t mvMod5 = (std::uint32_t) (mv - 5 * div5(mv));
      if (mvMod5 == 0) {
        vrIsTrailingZeros = multipleOfPowerOf5(mv, q);
      } else if (acceptBounds) {
        // Same as min(e2 + (~mm & 1), pow5Factor(mm)) >= q
        // <=> e2 + (~mm & 1) >= q && pow5Factor(mm) >= q
        // <=> true && pow5Factor(mm) >= q, since e2 >= q.
        vmIsTrailingZeros = multipleOfPowerOf5(mv - 1 - mmShift, q);
      } else {
        // Same as min(e2 + 1, pow5Factor(mp)) >= q.
        vp -= multipleOfPowerOf5(mv + 2, q);
      }
    }
  } else {
    // This expression is slightly faster than max(0, log10Pow5(-e2) - 1).
    const std::uint32_t q = log10Pow5(-e2) - (-e2 > 1);
    e10 = (std::int32_t) q + e2;
    const std::int32_t i = -e2 - (std::int32_t) q;
    constexpr std::int32_t double_pow5_bitcount = 121;
    const std::int32_t k = pow5bits(i) - double_pow5_bitcount;
    const std::int32_t j = (std::int32_t) q - k;
    vr = mulShiftAll(m2, double_pow5_split(i), j, &vp, &vm, mmShift);
    if (q <= 1) {
      // {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
      // mv = 4 * m2, so it always has at least two trailing 0 bits.
      vrIsTrailingZeros = true;
      if (acceptBounds) {
        // mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
        vmIsTrailingZeros = mmShift == 1;
      } else {
        // mp = mv + 2, so it always has at least one trailing 0 bit.
        --vp;
      }
    } else if (q < 63) { // TODO(ulfjack): Use a tighter bound here.
      // We need to compute min(ntz(mv), pow5Factor(mv) - e2) >= q - 1
      // <=> ntz(mv) >= q - 1 && pow5Factor(mv) - e2 >= q - 1
      // <=> ntz(mv) >= q - 1 (e2 is negative and -e2 >= q)
      // <=> (mv & ((1 << (q - 1)) - 1)) == 0
      // We also need to make sure that the left shift does not overflow.
      vrIsTrailingZeros = multipleOfPowerOf2(mv, q - 1);
    }
  }

  return d2d_step3_result{ vr, vp, vm, e10, vmIsTrailingZeros
                         , vrIsTrailingZeros, acceptBounds };
}

struct d2d_step4_result {
  std::uint64_t dec_mantissa;
  std::int32_t dec_exponent;
};

inline d2d_step4_result d2d_step4
    ( std::uint64_t vr
    , std::uint64_t vp
    , std::uint64_t vm
    , const std::int32_t e10
    , bool vmIsTrailingZeros
    , bool vrIsTrailingZeros
    , const bool acceptBounds )
{
  // Step 4: Find the shortest decimal representation in the interval of valid representations.
  std::int32_t removed = 0;
  uint8_t lastRemovedDigit = 0;
  std::uint64_t dec_mantissa;
  // On average, we remove ~2 digits.
  if (vmIsTrailingZeros || vrIsTrailingZeros) {
    // General case, which happens rarely (~0.7%).
    for (;;) {
      const std::uint64_t vpDiv10 = div10(vp);
      const std::uint64_t vmDiv10 = div10(vm);
      if (vpDiv10 <= vmDiv10) {
        break;
      }
      const std::uint32_t vmMod10 = (std::uint32_t) (vm - 10 * vmDiv10);
      const std::uint64_t vrDiv10 = div10(vr);
      const std::uint32_t vrMod10 = (std::uint32_t) (vr - 10 * vrDiv10);
      vmIsTrailingZeros &= vmMod10 == 0;
      vrIsTrailingZeros &= lastRemovedDigit == 0;
      lastRemovedDigit = (uint8_t) vrMod10;
      vr = vrDiv10;
      vp = vpDiv10;
      vm = vmDiv10;
      ++removed;
    }
    if (vmIsTrailingZeros) {
      for (;;) {
        const std::uint64_t vmDiv10 = div10(vm);
        const std::uint32_t vmMod10 = (std::uint32_t) (vm - 10 * vmDiv10);
        if (vmMod10 != 0) {
          break;
        }
        const std::uint64_t vpDiv10 = div10(vp);
        const std::uint64_t vrDiv10 = div10(vr);
        const std::uint32_t vrMod10 = (std::uint32_t) (vr - 10 * vrDiv10);
        vrIsTrailingZeros &= lastRemovedDigit == 0;
        lastRemovedDigit = (uint8_t) vrMod10;
        vr = vrDiv10;
        vp = vpDiv10;
        vm = vmDiv10;
        ++removed;
      }
    }
    if (vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0) {
      // Round even if the exact number is .....50..0.
      lastRemovedDigit = 4;
    }
    // We need to take vr + 1 if vr is outside bounds or we need to round up.
    dec_mantissa = vr + ((vr == vm && (!acceptBounds || !vmIsTrailingZeros)) || lastRemovedDigit >= 5);
  } else {
    // Specialized for the common case (~99.3%). Percentages below are relative to this.
    bool roundUp = false;
    const std::uint64_t vpDiv100 = div100(vp);
    const std::uint64_t vmDiv100 = div100(vm);
    if (vpDiv100 > vmDiv100) { // Optimization: remove two digits at a time (~86.2%).
      const std::uint64_t vrDiv100 = div100(vr);
      const std::uint32_t vrMod100 = (std::uint32_t) (vr - 100 * vrDiv100);
      roundUp = vrMod100 >= 50;
      vr = vrDiv100;
      vp = vpDiv100;
      vm = vmDiv100;
      removed += 2;
    }
    // Loop iterations below (approximately), without optimization above:
    // 0: 0.03%, 1: 13.8%, 2: 70.6%, 3: 14.0%, 4: 1.40%, 5: 0.14%, 6+: 0.02%
    // Loop iterations below (approximately), with optimization above:
    // 0: 70.6%, 1: 27.8%, 2: 1.40%, 3: 0.14%, 4+: 0.02%
    for (;;) {
      const std::uint64_t vpDiv10 = div10(vp);
      const std::uint64_t vmDiv10 = div10(vm);
      if (vpDiv10 <= vmDiv10) {
        break;
      }
      const std::uint64_t vrDiv10 = div10(vr);
      const std::uint32_t vrMod10 = (std::uint32_t) (vr - 10 * vrDiv10);
      roundUp = vrMod10 >= 5;
      vr = vrDiv10;
      vp = vpDiv10;
      vm = vmDiv10;
      ++removed;
    }
    // We need to take vr + 1 if vr is outside bounds or we need to round up.
    dec_mantissa = vr + (vr == vm || roundUp);
  }

  return {dec_mantissa, e10 + removed};
}

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


inline decimal_double to_decimal(double f)
{
  constexpr int mantissa_bits = 52;
  constexpr int exponent_bits = 11;
  constexpr int bias = 1023;
  
  const std::uint64_t bits = double_to_bits(f);
  const bool ieeeSign = bits >> 63;
  const std::uint64_t ieeeMantissa = bits & 0xFFFFFFFFFFFFFull;
  const std::uint32_t ieeeExponent = static_cast<std::uint32_t>((bits << 1) >> 53);

  if (ieeeExponent == 0x7FF) {
    if (ieeeMantissa == 0) {
      return {0, 0, ieeeSign, float_state::infinity};
    }
    else {
      return {0, 0, ieeeSign, float_state::nan};
    }
  }

  std::int32_t e2;
  std::uint64_t m2;
  if (ieeeExponent == 0) {
    // We subtract 2 so that the bounds computation has 2 additional bits.
    e2 = 1 - bias - mantissa_bits - 2;
    m2 = ieeeMantissa;
  } else {
    e2 = (std::int32_t) ieeeExponent - bias - mantissa_bits - 2;
    m2 = (1ull << mantissa_bits) | ieeeMantissa;
  }
  const bool even = (m2 & 1) == 0;
  const bool acceptBounds = even;


  // Step 2: Determine the interval of valid decimal representations.
  const std::uint64_t mv = 4 * m2;
  // Implicit bool -> int conversion. True is 1, false is 0.
  const std::uint32_t mmShift = ieeeMantissa != 0 || ieeeExponent <= 1;
  // We would compute mp and mm like this:
  // std::uint64_t mp = 4 * m2 + 2;
  // std::uint64_t mm = mv - 1 - mmShift;

  auto res3 = d2d_step_3(e2, m2, mv, mmShift, acceptBounds );
  auto res4 = d2d_step4( res3.vr, res3.vp, res3.vm, res3.e10
                       , res3.vmIsTrailingZeros, res3.vrIsTrailingZeros
                       , res3.acceptBounds );

  return {res4.dec_mantissa, res4.dec_exponent, ieeeSign, float_state::normal};
}

} // namespace ryu

#endif //  RYU_D2S_HPP
