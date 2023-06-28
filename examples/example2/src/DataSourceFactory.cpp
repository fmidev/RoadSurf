#include "DataSourceFactory.h"
#include "AsciiSource.h"
#include "QueryDataSource.h"
#include "RoadSurfSource.h"
#include "SmartMetSource.h"
#include <json/json.h>
#include <stdexcept>

namespace DataSourceFactory
{
// ----------------------------------------------------------------------
/*!
 * \brief Create a new data source
 */
// ----------------------------------------------------------------------

DataSource* create(const NFmiMetTime& pWallClock,
                   const NFmiMetTime& pStartTime,
                   const NFmiMetTime& pEndTime,
                   const Json::Value& pSource,
                   const Json::Value& pConfig)
{
  if (!pSource.isObject())
    throw std::logic_error("Data sources must be defined with JSON objects");

  Json::Value nulljson;

  auto jname = pSource.get("name", nulljson);
  if (jname.isNull())
    throw std::logic_error("Encountered an element in 'sources' without a name");
  std::cout << "Reading data source named '" << jname.asString() << "'" << std::endl;

  auto jtype = pSource.get("type", nulljson);
  if (jtype.isNull())
    throw std::logic_error("Data source '" + jname.asString() + "' has no type");

  auto type = jtype.asString();

  if (type == "file")
    return new QueryDataSource(pWallClock, pStartTime, pEndTime, pSource);
  if (type == "directory")
    return new QueryDataSource(pWallClock, pStartTime, pEndTime, pSource);
  if (type == "smartmet")
    return new SmartMetSource(pWallClock, pStartTime, pEndTime, pSource, pConfig);
  if (type == "ascii")
    return new AsciiSource(pWallClock, pStartTime, pEndTime, pSource);
  if (type == "RoadSurf")
    return new RoadSurfSource(pWallClock, pStartTime, pEndTime, pSource, pConfig);
  throw std::runtime_error("Unknown data source type '" + type + "'");
}
}  // namespace DataSourceFactory
