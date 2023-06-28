#pragma once

#include <json/json.h>
#include <memory>

class NFmiQueryData;
class NFmiFastQueryInfo;
struct InputData;
struct OutputData;
struct InputSettings;

std::unique_ptr<NFmiQueryData> create_coordinates_querydata(const InputSettings& pSettings,
                                                            const Json::Value& pJson);

std::unique_ptr<NFmiQueryData> create_grid_querydata(const InputSettings& pSettings,
                                                     const Json::Value& pJson);

struct WriteStride
{
  std::size_t start;  // first index of simulation data written out
  std::size_t step;   // index step size for writes
};

WriteStride get_write_stride(NFmiFastQueryInfo& qi, const InputSettings& pSettings);

void write_timeseries(NFmiFastQueryInfo& pQ,
                      long pLocIndex,
                      const InputData& pInput,
                      const OutputData& pOutput,
                      const WriteStride& pStride);
