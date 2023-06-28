#pragma once
#include "DataSource.h"
#include <string>
//>! Used to create input data entities
namespace DataSourceFactory
{
DataSource* create(const NFmiMetTime& pWallClock,
                   const NFmiMetTime& pStartTime,
                   const NFmiMetTime& pEndTime,
                   const Json::Value& pSource,
                   const Json::Value& pConfig);
}  // namespace DataSourceFactory
