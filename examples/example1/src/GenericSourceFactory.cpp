#include "GenericSourceFactory.h"
#include "JsonSource.h"
#include <json/json.h>
#include <stdexcept>

namespace GenericSourceFactory
{
// ----------------------------------------------------------------------
/*!
 * \brief Create a new data source
 */
// ----------------------------------------------------------------------

GenericSource* create(const time_t& pWallClock,
                   const time_t& pStartTime,
                   const time_t& pEndTime,
                   const Json::Value& pSource,
                   const Json::Value& pConfig,
                   const int SimLen,
                   const int DTSecs)
{
  //Check if soruce is given
  if (!pSource.isObject())
    throw std::logic_error("Data sources must be defined with JSON objects");

  Json::Value nulljson;

  //Get data source name
  auto jname = pSource.get("name", nulljson);
  if (jname.isNull())
    throw std::logic_error("Encountered an element in 'sources' without a name");
  std::cout << "Reading data source named '" << jname.asString() << "'" << std::endl;

  //Get data source type, e.g. "json"
  auto jtype = pSource.get("type", nulljson);
  if (jtype.isNull())
    throw std::logic_error("Data source '" + jname.asString() + "' has no type");

  auto type = jtype.asString();
  //if type is json, create new json source (Input data is read at this point)
  if (type == "json")
    return new JsonSource(pWallClock, pStartTime, pEndTime, pSource, pConfig,SimLen,DTSecs);
  throw std::runtime_error("Unknown data source type '" + type + "'");
}
}  // namespace GenericSourceFactory
