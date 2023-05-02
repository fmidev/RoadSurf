#include "SkyView.h"
#include "JsonTools.h"
#include "SkyVariables.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

//------------------------------------------------------------------------
//
// Read sky view factors and local horizon angles
//
//------------------------------------------------------------------------
SkyView::SkyView(const Json::Value& pJson )
{
  // check if parameters are given in json
  if (pJson.isNull())
    return;
  Json::Value nulljson;

  int id = 0;
  double lat = 0;
  double lon = 0;
  std::string name;

  // make default local horizon array
  std::vector<double> default_horizons(360);
  fill(default_horizons.begin(), default_horizons.end(), 0.0);

  // Read sky view values if filename is given in config file
  auto tmp = pJson.get("sky_view_file", nulljson);
  if (!tmp.isNull())
  {
    //Get file name
    std::string sky_view_file = tmp.asString();
    std::ifstream in(sky_view_file);
    //Check if file exists
    if (!in)
      throw std::runtime_error("Failed to open '" + sky_view_file + "'");
    double sky_view;
    //Read in data
    while (in.good())
    {
      in >> id >> name >> lat >> lon >> sky_view;
      SkyVariables skyVar;  // sky view factor and horizons
                            // are saved in this structure
      // check if sky view factor is on right range
      if (sky_view >= 0.0 && sky_view <= 1.0)
        skyVar.sky_view = sky_view;
      else
        skyVar.sky_view = 1.0;
      //Insert dafault horizon angles in case they are not given
      skyVar.local_horizons = default_horizons;
      mData.insert(std::make_pair(id, skyVar));
    }
  }
  // Read local horizon values if filename is given in config file
  auto tmp2 = pJson.get("local_horizon_file", nulljson);
  if (!tmp2.isNull())
  {
    //Get file name
    std::string local_horizon_file = tmp2.asString();
    std::ifstream in(local_horizon_file);
    //Check if file exists
    if (!in)
      throw std::runtime_error("Failed to open '" + local_horizon_file + "'");
    std::string line;
    std::string value;
    //Read in data
    while (std::getline(in, line))
    {
      std::istringstream iss(line);
      std::vector<double> local_horizons;
      int l = 1;
      //Read line
      while (iss >> value)
      {
        if (l == 1) //pointID
          id = std::stoi(value);
        else if (l == 2) //name
          name = value;
        else if (l == 3) //latitude
          lat = std::stod(value);
        else if (l == 4)  //longitude
          lon = std::stod(value);
        else //local horizon
          local_horizons.push_back(std::stod(value));
        l++;
      }
      //If sky view factors are already saved to mdata
      if (!mData.empty())
      {
	//Find the place in mData with same pointID
        std::map<int, SkyVariables>::iterator it;
        it = mData.find(id);
	//If found, set local horizons to madata
        if (it != mData.end())
        {
          it->second.local_horizons = local_horizons;
        }
        else
        {
	  //Else give sky view factor a default value and save
	  //local horizon angles
          SkyVariables skyVar;
          skyVar.local_horizons = local_horizons;
          skyVar.sky_view = 1.0;
          mData.insert(std::make_pair(id, skyVar));
        }
      }
      else
      {
	//If mdata is empty give sky view factr a default value
	//and save local horizon angels
        SkyVariables skyVar;
        skyVar.sky_view = 1.0;
        skyVar.local_horizons = local_horizons;
        mData.insert(std::make_pair(id, skyVar));
      }
    }
  }

}
//Get sky view factor and local horizon angles for point
void SkyView::GetSkyVariables(LocalParameters& lParameters,
                              InputData& pData,
                              int pointID)
{

  // Check if there is data for point
  auto myData = mData.find(pointID);
  if (myData == mData.end())
    return;
  //Save sky view factor and locan horizon angle to input data
  const auto& SkyVar = myData->second;
  lParameters.sky_view = SkyVar.sky_view;
  pData.local_horizons = SkyVar.local_horizons;
}
