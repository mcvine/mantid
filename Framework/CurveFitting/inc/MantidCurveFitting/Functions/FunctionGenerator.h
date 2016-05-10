#ifndef MANTID_CURVEFITTING_FUNCTIONGENERATOR_H_
#define MANTID_CURVEFITTING_FUNCTIONGENERATOR_H_

#include "MantidAPI/IFunction.h"

namespace Mantid {
namespace CurveFitting {
namespace Functions {
/**
FunctionGenerator is a partial implementation of IFunction that defines a
function consisting of two parts: the source and the target. The source
function generates the target function which in turn is used to calculate
the output. Concrete functions subclassing FunctionGenerator must implement
the following virtual methods:
  
  name()
  category()
  buildTargetFunction()
  updateTargetFunction()

Parameters and attributes of both source and target functions become
parameters (attributes) of FunctionGenerator without changing names.
Virtual method isSourceName(name) decides to which function a parameter
belongs to. By default if a name has the signature of a composite function
(f<number>.name) then it is attributed to the target function.

Copyright &copy; 2007-8 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport FunctionGenerator : public API::IFunction {
public:
  /// Constructor
  FunctionGenerator(API::IFunction_sptr source);

  /// @name Overrides implementing composition of two functions:
  /// m_source and m_target.
  //@{
  /// Set i-th parameter
  virtual void setParameter(size_t, const double &value,
                            bool explicitlySet = true) override;
  /// Set i-th parameter description
  virtual void setParameterDescription(size_t,
                                       const std::string &description) override;
  /// Get i-th parameter
  virtual double getParameter(size_t i) const override;
  /// Set parameter by name.
  virtual void setParameter(const std::string &name, const double &value,
                            bool explicitlySet = true) override;
  /// Set description of parameter by name.
  virtual void setParameterDescription(const std::string &name,
                                       const std::string &description) override;
  /// Get parameter by name.
  virtual double getParameter(const std::string &name) const override;
  /// Total number of parameters
  virtual size_t nParams() const override;
  /// Returns the index of parameter name
  virtual size_t parameterIndex(const std::string &name) const override;
  /// Returns the name of parameter i
  virtual std::string parameterName(size_t i) const override;
  /// Returns the description of parameter i
  virtual std::string parameterDescription(size_t i) const override;
  /// Checks if a parameter has been set explicitly
  virtual bool isExplicitlySet(size_t i) const override;
  /// Get the fitting error for a parameter
  virtual double getError(size_t i) const override;
  /// Set the fitting error for a parameter
  virtual void setError(size_t i, double err) override;

  /// Check if a declared parameter i is fixed
  virtual bool isFixed(size_t i) const override;
  /// Removes a declared parameter i from the list of active
  virtual void fix(size_t i) override;
  /// Restores a declared parameter i to the active status
  virtual void unfix(size_t i) override;

  /// Return parameter index from a parameter reference.
  virtual size_t
  getParameterIndex(const API::ParameterReference &ref) const override;
  /// Tie a parameter to other parameters (or a constant)
  virtual API::ParameterTie *tie(const std::string &parName,
                                 const std::string &expr,
                                 bool isDefault = false) override;
  /// Apply the ties
  virtual void applyTies() override;
  /// Remove all ties
  virtual void clearTies() override;
  // Unhide base class function: removeTie(string).
  using IFunction::removeTie;
  /// Removes i-th parameter's tie
  virtual bool removeTie(size_t i) override;
  /// Get the tie of i-th parameter
  virtual API::ParameterTie *getTie(size_t i) const override;

  /// Add a constraint to function
  virtual void addConstraint(API::IConstraint *ic) override;
  /// Get constraint of i-th parameter
  virtual API::IConstraint *getConstraint(size_t i) const override;
  /// Remove a constraint
  virtual void removeConstraint(const std::string &parName) override;

  /// Set up the function for a fit.
  virtual void setUpForFit() override;

  /// Build target function.
  virtual void buildTargetFunction() const = 0;

protected:
  /// Declare a new parameter
  virtual void declareParameter(const std::string &name, double initValue = 0,
                                const std::string &description = "") override;

  /// Add a new tie. Derived classes must provide storage for ties
  virtual void addTie(API::ParameterTie *tie) override;
  //@}

public:
  /** @name Attributes */
  //@{
  /// Returns the number of attributes associated with the function
  virtual size_t nAttributes() const override;
  /// Returns a list of attribute names
  virtual std::vector<std::string> getAttributeNames() const override;
  /// Return a value of attribute attName
  virtual Attribute getAttribute(const std::string &name) const override;
  /// Set a value to attribute attName
  virtual void setAttribute(const std::string &name,
                            const Attribute &) override;
  /// Check if attribute attName exists
  virtual bool hasAttribute(const std::string &name) const override;
  //@}

  /// Evaluate the function
  void function(const API::FunctionDomain &domain,
                API::FunctionValues &values) const override;

protected:
  /// overwrite IFunction base class method, which declare function parameters
  void init() override;

  /// Test if a name (parameter's or attribute's) belongs to m_source
  virtual bool isSourceName(const std::string &aName) const;
  /// Update target function.
  virtual void updateTargetFunction() const = 0;
  /// Update target function if necessary.
  void checkTargetFunction() const;
  /// Function that calculates peak centres and intensities.
  API::IFunction_sptr m_source;
  /// Function that actually callculates the spectrum.
  mutable API::IFunction_sptr m_target;
  /// Cached number of parameters in m_source.
  size_t m_nOwnParams;
  /// Flag indicating that updateTargetFunction() is required.
  mutable bool m_dirty;
};

} // namespace Functions
} // namespace CurveFitting
} // namespace Mantid

#endif /*MANTID_CURVEFITTING_FUNCTIONGENERATOR_H_*/