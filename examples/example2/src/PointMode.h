#pragma once
#include <json/json.h>

enum class PointMode
{
  Coordinate,   // single latitude + longitude for reseacrh runs
  Coordinates,  // multiple coordinates
  Grid          // projected NxM grid
};

PointMode pointmode(const Json::Value& pJson);
