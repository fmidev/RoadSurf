#pragma once

#include "InputData.h"
#include "LocalParameters.h"
#include <json/json.h>

//Defines generic source class for input data
class GenericSource
{
 public:
  GenericSource() = default;

  GenericSource(const GenericSource& pOther) = delete;
  GenericSource& operator=(const GenericSource& pOther) = delete;
  GenericSource(GenericSource&& pOther) = delete;
  GenericSource& operator=(GenericSource&& pOther) = delete;

  virtual ~GenericSource() = default; //Destructor
  //Get input data for point
  virtual void GetWeather(InputData& pData,int pointID) const = 0;
  //Get vector with location data
  virtual std::vector<LocalParameters> GetLocations() const = 0;
  //Get vector with pointIDs
  virtual std::vector<int> GetPointIDs() const = 0;
  //Get index of latest observation
  virtual int GetLatestObsIndex(int statPos) const = 0;

  bool is_observation = false;

 private:
};  // class GenericSource
