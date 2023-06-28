#pragma once
#include <json/json.h>
#include <smartmet/newbase/NFmiStation.h>
#include <list>
std::list<NFmiStation> read_filename_locations(const Json::Value& pJson);
std::list<NFmiStation> read_filename_locations(const std::string& pFilename);
