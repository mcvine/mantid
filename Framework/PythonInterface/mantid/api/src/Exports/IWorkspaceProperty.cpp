#include "MantidAPI/IWorkspaceProperty.h"
#include <boost/python/class.hpp>

void export_IWorkspaceProperty() {
  using namespace boost::python;
  using Mantid::API::IWorkspaceProperty;

  class_<IWorkspaceProperty, boost::noncopyable>("IWorkspaceProperty", no_init)
      .def("isOptional", &IWorkspaceProperty::isOptional, arg("self"),
           "Is the input workspace property optional")
      .def("isLocking", &IWorkspaceProperty::isLocking, arg("self"),
           "Will the workspace be locked when starting an algorithm");
}
