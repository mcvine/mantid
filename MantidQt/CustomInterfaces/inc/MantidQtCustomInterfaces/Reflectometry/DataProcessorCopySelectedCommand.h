#ifndef MANTID_CUSTOMINTERFACES_DATAPROCESSORCOPYSELECTEDCOMMAND_H
#define MANTID_CUSTOMINTERFACES_DATAPROCESSORCOPYSELECTEDCOMMAND_H

#include "MantidQtCustomInterfaces/Reflectometry/DataProcessorCommandBase.h"

namespace MantidQt {
namespace CustomInterfaces {
/** @class DataProcessorCopySelectedCommand

DataProcessorCopySelectedCommand defines the action "Copy Selected"

Copyright &copy; 2011-14 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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

File change history is stored at: <https://github.com/mantidproject/mantid>.
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DataProcessorCopySelectedCommand : public DataProcessorCommandBase {
public:
  DataProcessorCopySelectedCommand(DataProcessorPresenter *tablePresenter)
      : DataProcessorCommandBase(tablePresenter){};
  virtual ~DataProcessorCopySelectedCommand(){};

  void execute() override {
    m_presenter->notify(DataProcessorPresenter::CopySelectedFlag);
  };
  std::string name() override { return std::string("Copy Selected"); }
  std::string icon() override { return std::string("://copy.png"); }
};
}
}
#endif /*MANTID_CUSTOMINTERFACES_DATAPROCESSORCOPYSELECTEDCOMMAND_H*/