#include "JsonTools.h"
#include <fstream>
#include <memory>
#include <json/json.h>
#include <json/reader.h>
#include <json/value.h>
//Set pvalue to double given in json
bool override(double* pValue, const Json::Value& pJson, const char* pName)
{
  //Check if json is given
  if (pJson.isNull())
    return false;

  Json::Value nulljson;
  //Check if value is giben for pName in json
  auto tmp = pJson.get(pName, nulljson);
  if (tmp.isNull())
    return false;
  //Set pvalue to value given in json
  *pValue = tmp.asDouble();
  return true;
}
//Set pvalue to int given in json
bool override(int* pValue, const Json::Value& pJson, const char* pName)
{
  //Check if json is given
  if (pJson.isNull())
    return false;

  //Check that value is given for pName in json
  Json::Value nulljson;
  auto tmp = pJson.get(pName, nulljson);
  if (tmp.isNull())
    return false;

  //Set pvalue to inreger given in json
  *pValue = tmp.asInt();
  return true;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read a file into a string
 */
// ----------------------------------------------------------------------

std::string read_file(const std::string& pFilename)
{
  std::ifstream in(pFilename.c_str());
  if (!in)
    throw std::logic_error("Failed to open '" + pFilename + "' for reading");

  // A fast way to read an entire file into a string.
  std::string content;
  content.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());

  return content;
}

// ----------------------------------------------------------------------
/*!
 * \brief Read the JSON configuration file
 */
// ----------------------------------------------------------------------

Json::Value read_json(const std::string& pFilename)
{
  auto text = read_file(pFilename);

  Json::Value json;
 // Json::Reader reader;
 // if (reader.parse(text, json))
 //   return json;
  JSONCPP_STRING err;
  const auto jsonLength = static_cast<int>(text.length());
  Json::CharReaderBuilder builder;
  const std::unique_ptr<Json::CharReader> reader(builder.newCharReader());
  if (reader->parse(text.c_str(), text.c_str() + jsonLength, &json,
                       &err)) {
     return json;
  }

  throw std::runtime_error("failed to read "+pFilename);
}
