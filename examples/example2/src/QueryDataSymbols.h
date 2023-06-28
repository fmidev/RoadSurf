// ======================================================================
/*!
 * \file
 * \brief Interface of a symbol value extractor for expressions
 */
// ======================================================================

#pragma once

#include <smartmet/newbase/NFmiFastQueryInfo.h>
#include <stx-exparser/AnyScalar.h>
#include <stx-exparser/ExpressionParser.h>

class QueryDataSymbols : public stx::SymbolTable
{
 public:
  QueryDataSymbols() = default;

  using paramlist_type = std::vector<stx::AnyScalar>;

  stx::AnyScalar lookupVariable(const std::string& varname) const override;

  stx::AnyScalar processFunction(const std::string& funcname,
                                 const paramlist_type& paramlist) const override;

  // Additional methods:

  void setData(NFmiFastQueryInfo& qinfo) { mInfo = &qinfo; }
  void setLocation(const NFmiPoint& latlon) { mLatLon = latlon; }

 private:
  NFmiFastQueryInfo* mInfo = nullptr;  // not owner
  NFmiPoint mLatLon;

};  // class QueryDataSymbols

// ======================================================================
