#include "PointMode.h"

PointMode pointmode(const Json::Value& pJson)
{
  Json::Value nulljson;
  auto json = pJson.get("points", nulljson);

  bool has_latlon = (json.isMember("longitude") || json.isMember("latitude"));
  bool has_coordinates = (json.isMember("keyword") || json.isMember("filename"));
  bool has_grid = json.isMember("projection");

  if (has_latlon && !has_coordinates && !has_grid)
    return PointMode::Coordinate;

  if (has_coordinates && !has_latlon && !has_grid)
    return PointMode::Coordinates;

  if (has_grid && !has_latlon && !has_coordinates)
    return PointMode::Grid;

  throw std::runtime_error(
      "Multiple points definitions, use only one of keyword/filename, projection or longitude & "
      "latitude");
}
