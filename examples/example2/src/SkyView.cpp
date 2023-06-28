#include "SkyView.h"
#include "JsonTools.h"
#include "SkyVariables.h"
#include "Weather.h"
#include <fstream>
#include <iostream>
#include <string>

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
    std::string sky_view_file = tmp.asString();
    std::ifstream in(sky_view_file);
    if (!in)
      throw std::runtime_error("Failed to open '" + sky_view_file + "'");
    double sky_view;
    while (in.good())
    {
      in >> id >> name >> lat >> lon >> sky_view;
      // add station to station tree
      mStationTree.insert(Station(lon, lat, id));
      // check if sky view factor is on right range
      SkyVariables skyVar;  // sky view factor and horizons
                            // are saved in this structure
      if (sky_view >= 0.0 && sky_view <= 1.0)
        skyVar.sky_view = sky_view;
      else
        skyVar.sky_view = 1.0;
      //
      skyVar.local_horizons = default_horizons;
      mData.insert(std::make_pair(id, skyVar));
    }
  }
  // Read local horizon values if filename is given in config file
  auto tmp2 = pJson.get("local_horizon_file", nulljson);
  if (!tmp2.isNull())
  {
    std::string local_horizon_file = tmp2.asString();
    // std::string sky_view_file,std::string local_horizons_file
    std::ifstream in(local_horizon_file);
    if (!in)
      throw std::runtime_error("Failed to open '" + local_horizon_file + "'");
    std::string line;
    std::string value;
    while (std::getline(in, line))
    {
      std::istringstream iss(line);
      std::vector<double> local_horizons;
      int l = 1;
      while (iss >> value)
      {
        if (l == 1)
          id = std::stoi(value);
        else if (l == 2)
          name = value;
        else if (l == 3)
          lat = std::stod(value);
        else if (l == 4)
          lon = std::stod(value);
        else
          local_horizons.push_back(std::stod(value));
        l++;
      }
      if (!mData.empty())
      {
        std::map<int, SkyVariables>::iterator it;
        it = mData.find(id);
        if (it != mData.end())
        {
          it->second.local_horizons = local_horizons;
        }
        else
        {
          SkyVariables skyVar;
          skyVar.local_horizons = local_horizons;
          skyVar.sky_view = 1.0;
          mData.insert(std::make_pair(id, skyVar));
        }
      }
      else
      {
        mStationTree.insert(Station(lon, lat, id));
        SkyVariables skyVar;
        skyVar.sky_view = 1.0;
        skyVar.local_horizons = local_horizons;
        mData.insert(std::make_pair(id, skyVar));
      }
    }
  }

  // std::map<int,SkyVariables>:: iterator tt = mData.begin();
  // while (tt!=mData.end())
  // {
  //   std::cout<<tt->first<<" "<<tt->second.sky_view<<"
  //   "<<tt->second.local_horizons.size()<<std::endl; auto lh=tt->second.local_horizons; for
  //   (std::size_t i=0; i<lh.size(); ++i)
  //      std::cout << i<<" "<<lh[i]<<std::endl;
  //  tt++;
  // }
}
void SkyView::GetSkyVariables(LocalParameters& lParameters,
                              InputData& pData,
                              const NFmiPoint& pLatLon)
{
  if (is_missing(pLatLon.X()) || is_missing(pLatLon.Y()))
    return;
  double mRadius = 0.1;  // km
  // No station near enough?
  auto nearest_station =
      mStationTree.nearest(Station{pLatLon.X(), pLatLon.Y()}, Station::ChordLength(mRadius));
  if (!nearest_station)
    return;

  // No data for the station?!?!
  auto myData = mData.find(nearest_station->ID());
  if (myData == mData.end())
    return;

  const auto& SkyVar = myData->second;
  lParameters.sky_view = SkyVar.sky_view;
  pData.local_horizons = SkyVar.local_horizons;
}
