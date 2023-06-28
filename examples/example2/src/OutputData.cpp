#include "OutputData.h"
#include <fmt/format.h>

OutputData::OutputData(int pSize)
    : TsurfOut(pSize, -9999.0),
      SnowOut(pSize, -9999.0),
      WaterOut(pSize, -9999.0),
      IceOut(pSize, -9999.0),
      DepositOut(pSize, -9999.0),
      Ice2Out(pSize, -9999.0)
{
}

OutputPointers OutputData::pointers() const
{
  return OutputPointers{static_cast<int>(TsurfOut.size()),
                            TsurfOut.data(),
                            SnowOut.data(),
                            WaterOut.data(),
                            IceOut.data(),
                            DepositOut.data(),
                            Ice2Out.data()};
}

std::ostream& operator<<(std::ostream& out, const OutputData& data)
{
  out << "i  Tr   Snow Water Ice Depo Ice2\n";

  for (std::size_t i = 0; i < data.TsurfOut.size(); i++)
  {
    out << fmt::format(
        "{} {:.1f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n",
        i,
        data.TsurfOut[i],
        data.SnowOut[i],
        data.WaterOut[i],
        data.IceOut[i],
        data.DepositOut[i],
        data.Ice2Out[i]);
  }
  return out;
}
