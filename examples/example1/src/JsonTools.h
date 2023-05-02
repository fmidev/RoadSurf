#pragma once
#include <json/json.h>

Json::Value read_json(const std::string& filename);

bool override(double* pValue, const Json::Value& pJson, const char* pName);
bool override(int* pValue, const Json::Value& pJson, const char* pName);
