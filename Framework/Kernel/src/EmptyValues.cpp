//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------
#include "MantidKernel/EmptyValues.h"

namespace Mantid {

/**
 * Returns what we consider an "empty" integer within a property
 * @returns An flag value
 */
int32_t EMPTY_INT() { return emptyValue<int32_t>(); }

/**
 * Returns what we consider an "empty" long within a property
 * @returns An flag value
 */
int64_t EMPTY_LONG() { return emptyValue<int64_t>(); }

/**
 * Returns what we consider an "empty" double within a property
 * @returns An flag value
 */
double EMPTY_DBL() { return emptyValue<double>(); }

} // namespace Mantid
