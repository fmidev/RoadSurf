#include "JsonTools.h"
#include <fstream>

bool override(double* pValue, const Json::Value& pJson, const char* pName)
{
  if (pJson.isNull())
    return false;

  Json::Value nulljson;
  auto tmp = pJson.get(pName, nulljson);
  if (tmp.isNull())
    return false;

  *pValue = tmp.asDouble();
  return true;
}

bool override(int* pValue, const Json::Value& pJson, const char* pName)
{
  if (pJson.isNull())
    return false;

  Json::Value nulljson;
  auto tmp = pJson.get(pName, nulljson);
  if (tmp.isNull())
    return false;

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

  // A fast way to read an entire file into a string. Not sure if there are faster ways.
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
  Json::Reader reader;
  if (reader.parse(text, json))
    return json;

  throw std::logic_error(reader.getFormattedErrorMessages());
}
