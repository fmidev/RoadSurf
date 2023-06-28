
#include "QueryDataSymbols.h"
#include <smartmet/newbase/NFmiEnumConverter.h>
#include <smartmet/newbase/NFmiParameterName.h>
#include <stdexcept>

using namespace std;

// ----------------------------------------------------------------------
/*!
 * \brief Name converter
 */
// ----------------------------------------------------------------------

FmiParameterName toenum(const string& varname)
{
  static NFmiEnumConverter converter;
  return static_cast<FmiParameterName>(converter.ToEnum(varname));
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish variable value
 */
// ----------------------------------------------------------------------

stx::AnyScalar QueryDataSymbols::lookupVariable(const string& varname) const
{
  // Constants

  if (varname == "PI")
    return {3.14159265358979323846};

  // Querydata variables

  FmiParameterName p = toenum(varname);

  if (p == kFmiBadParameter)
    throw runtime_error("Unrecognized variable name: " + varname);

  if (!mInfo->Param(p))
    throw runtime_error("Parameter " + varname + " is not available in the querydata");

  double value = mInfo->InterpolatedValue(mLatLon);
  return {value};
}

// ----------------------------------------------------------------------
/*!
 * \brief Establish function value
 */
// ----------------------------------------------------------------------

stx::AnyScalar QueryDataSymbols::processFunction(const string& funcname,
                                                 const paramlist_type& paramlist) const
{
  if (funcname == "missing")
  {
    if (paramlist.size() != 1)
      throw runtime_error("missing function takes exactly one argument");
    return {paramlist[0] == 32700};
  }

  throw runtime_error("Unrecognized function: " + funcname);
}
