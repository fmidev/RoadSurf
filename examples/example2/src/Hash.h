#pragma once

#include <smartmet/newbase/NFmiPoint.h>

namespace std
{
/*
 * Hash function for coordinates
 */

template <>
struct hash<NFmiPoint>
{
  using argument_type = NFmiPoint;
  using result_type = std::size_t;

  result_type operator()(argument_type const& p) const noexcept
  {
    const result_type h1 = std::hash<double>{}(p.X());
    const result_type h2 = std::hash<double>{}(p.Y());
    return h1 ^ (h2 << 1);
  }
};

}  // namespace std
