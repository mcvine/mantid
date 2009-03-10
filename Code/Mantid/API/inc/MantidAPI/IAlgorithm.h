#ifndef MANTID_KERNEL_IALGORITHM_H_
#define MANTID_KERNEL_IALGORITHM_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include <string>
#include <vector>
#include "MantidKernel/INamedInterface.h"
#include "MantidKernel/Property.h"
#include "MantidKernel/IPropertyManager.h"

#include <Poco/ActiveResult.h>

namespace Poco
{
    class AbstractObserver;
}

namespace Mantid
{
namespace API
{
// Declaration of the interface ID ( interface id, major version, minor version)
// RJT: Have not yet imported the code for this (in IInterface.h in Gaudi)
//static const InterfaceID IID_IAlgorithm("IAlgorithm", 3 , 0); 

/** @class IAlgorithm IAlgorithm.h Kernel/IAlgorithm.h

 IAlgorithm is the interface implemented by the Algorithm base class.
 Concrete algorithms, derived from the Algorithm base class are controlled 
 via this interface.

 @author Russell Taylor, Tessella Support Services plc
 @author Based on the Gaudi class of the same name (see http://proj-gaudi.web.cern.ch/proj-gaudi/)
 @date 11/09/2007
 
 Copyright &copy; 2007 STFC Rutherford Appleton Laboratories

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

 File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>.    
 Code Documentation is available at: <http://doxygen.mantidproject.org>
 */

/**
     As we have multiple interfaces to the same logical algorithm (Algorithm & AlgorithmProxy)
     we need a way of uniquely identifying managed algorithms. It can be AlgorithmID.
 */
typedef void* AlgorithmID;

class DLLExport IAlgorithm : virtual public Kernel::INamedInterface, virtual public Kernel::IPropertyManager
{
public:
  // Retrieve interface ID
  //    static const InterfaceID& interfaceID() { return IID_IAlgorithm; }

  /// Virtual destructor (always needed for abstract classes)
  virtual ~IAlgorithm() {}

  /// function to return a name of the algorithm, must be overridden in all algorithms
  virtual const std::string name() const = 0;

  /// function to return a version of the algorithm, must be overridden in all algorithms
  virtual const int version() const = 0;
  
  /// function to return a category of the algorithm.
  virtual const std::string category() const = 0;

  /// Algorithm ID. Unmanaged algorithms return 0 (or NULL?) values. Managed ones have non-zero.
  virtual AlgorithmID getAlgorithmID()const = 0;

  /** Initialization method invoked by the framework. This method is responsible
   *  for any bookkeeping of initialization required by the framework itself.
   *  It will in turn invoke the init() method of the derived algorithm,
   *  and of any sub-algorithms which it creates.
   */
  virtual void initialize() = 0;

  /// System execution. This method invokes the exec() method of a concrete algorithm.
  virtual bool execute() = 0;

  /// Asynchronous execution.
  virtual Poco::ActiveResult<bool> executeAsync() = 0;

  /// Check whether the algorithm is initialized properly
  virtual bool isInitialized() const = 0;
  /// Check whether the algorithm has already been executed
  virtual bool isExecuted() const = 0;

  /// Raises the cancel flag. interuption_point() method if called inside exec() checks this flag
  /// and if true terminates the algorithm.
  virtual void cancel()const = 0;

  /// True if the algorithm is running asynchronously.
  virtual bool isRunningAsync() = 0;

  /// True if the algorithm is running.
  virtual bool isRunning() = 0;

  /// To query whether algorithm is a child. Default to false
  virtual bool isChild() const = 0;

  /** To set whether algorithm is a child.
   *  @param isChild True - the algorithm is a child algorithm.  False - this is a full managed algorithm.
   */
  virtual void setChild(const bool isChild) = 0;

  /// Add an observer for a notification
  virtual void addObserver(const Poco::AbstractObserver& observer)const = 0;

  /// Remove an observer
  virtual void removeObserver(const Poco::AbstractObserver& observer)const = 0;

};

typedef boost::shared_ptr<IAlgorithm> IAlgorithm_sptr;
typedef boost::shared_ptr<const IAlgorithm> IAlgorithm_const_sptr;

} // namespace API
} // namespace Mantid

#endif /*MANTID_KERNEL_IALGORITHM_H_*/
