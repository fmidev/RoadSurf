#include "InputData.h"
#include <fmt/format.h>

InputData::InputData(int pSize)
    : tair(pSize, -9999.9),
      tdew(pSize, -9999.9),
      VZ(pSize, -9999.9),
      Rhz(pSize, -9999.9),
      prec(pSize, -9999.9),
      SW(pSize, -9999.9),
      LW(pSize, -9999.9),
      SW_dir(pSize, -9999.9),
      LW_net(pSize, -9999.9),
      TSurfObs(pSize, -9999.9),
      PrecPhase(pSize, -9999),
      local_horizons(360, 0.0),
      Depth(pSize, -9999.9),
      year(pSize, -9999),
      month(pSize, -9999),
      day(pSize, -9999),
      hour(pSize, -9999),
      minute(pSize, -9999),
      second(pSize, -9999)
{
}

InputPointers InputData::pointers() const
{
  return InputPointers{static_cast<int>(tair.size()),
                      tair.data(),
                      tdew.data(),
                      VZ.data(),
                      Rhz.data(),
                      prec.data(),
                      SW.data(),
                      LW.data(),
                      SW_dir.data(),
                      LW_net.data(),
                      TSurfObs.data(),
                      PrecPhase.data(),
                      local_horizons.data(),
                      Depth.data(),
                      year.data(),
                      month.data(),
                      day.data(),
                      hour.data(),
                      minute.data(),
                      second.data()};
}

std::ostream& operator<<(std::ostream& out, const InputData& data)
{
  out << "i       date     time   T  TDEW VZ   RH  RR    SW    LW  SW_dir LW_net TSurf PrecPhase "
         "Depth\n";
  for (std::size_t i = 0; i < data.tair.size(); i++)
    out << fmt::format(
        "{} {:02}.{:02}.{:04} {:02}:{:02}:{:02} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} "
        "{:.1f} {:.1f} {:.1f} {} {:.4f}"
        "\n",
        i,
        data.day[i],
        data.month[i],
        data.year[i],
        data.hour[i],
        data.minute[i],
        data.second[i],
        data.tair[i],
        data.tdew[i],
        data.VZ[i],
        data.Rhz[i],
        data.prec[i],
        data.SW[i],
        data.LW[i],
        data.SW_dir[i],
        data.LW_net[i],
        data.TSurfObs[i],
        data.PrecPhase[i],
        data.Depth[i]);

  return out;
}
