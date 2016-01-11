// Includes
#include "MantidMDAlgorithms/Quantification/Models/QCoordinate.h"
#include <cmath>

namespace Mantid {
namespace MDAlgorithms {
DECLARE_FOREGROUNDMODEL(QCoordinate)

namespace {
const char *QCOORD_ATTR = "Coord";
}

QCoordinate::QCoordinate() : ForegroundModel(), m_coord(0) {}

/**
 * Initialize the model
 */
void QCoordinate::init() {
  declareAttribute(QCOORD_ATTR, API::IFunction::Attribute("H"));
}

/**
 * Called when an attribute is set from the Fit string
 * @param name :: The name of the attribute
 * @param attr :: The value of the attribute
 */
void QCoordinate::setAttribute(const std::string &name,
                               const API::IFunction::Attribute &attr) {
  if (name == QCOORD_ATTR) {
    static std::map<std::string, size_t> coords = {
        {"H", 0}, {"K", 1}, {"L", 2}, {"En", 3}, {"Unity", 4}};
    auto it = coords.find(attr.asString());
    if (it != coords.end()) {
      m_coord = it->second;
    } else
      throw std::invalid_argument(
          "Unknown coordinate name passed to QCoordinate model");
  } else
    ForegroundModel::setAttribute(name, attr);
}

/**
 * Calculates the scattering intensity
 * @param exptSetup :: Details of the current experiment
 * @param point :: The axis values for the current point
 * @return The weight contributing from this point
 */
double
QCoordinate::scatteringIntensity(const API::ExperimentInfo &exptSetup,
                                 const std::vector<double> &point) const {
  UNUSED_ARG(exptSetup);
  double weight(0.0);
  if (m_coord < 3) {
    double qh, qk, ql;
    toHKL(exptSetup, point[0], point[1], point[2], qh, qk, ql);
    switch (m_coord) {
    case 0:
      weight = qh;
      break;
    case 1:
      weight = qk;
      break;
    case 2:
      weight = ql;
      break;
    }
  } else if (m_coord == 3)
    weight = point[3];
  else {
    weight = 1.0;
  }
  return weight;
}
}
}
