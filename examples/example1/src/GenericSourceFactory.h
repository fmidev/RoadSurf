#pragma once
#include "GenericSource.h"
#include <string>
#include <time.h>
//>! Used to create input data entities
namespace GenericSourceFactory
{
GenericSource* create(const time_t& pWallClock,
                   const time_t& pStartTime,
                   const time_t& pEndTime,
                   const Json::Value& pSource,
                   const Json::Value& pConfig,
                   const int SimLen,
                   const int DTSecs);
}  // namespace GenericSourceFactory
