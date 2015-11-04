#ifndef MANTIDABOUT_H
#define MANTIDABOUT_H

//----------------------
// Includes
//----------------------
#include "ui_MantidAbout.h"
#include <MantidQtAPI/MantidDialog.h>

	/** 
    This class implements About MantidPlot dialog for mnatid help menu.

    @author Sofia Antony, ISIS, RAL
    @date 13/01/2010

    Copyright &copy; 2009 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge National Laboratory & European Spallation Source

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

class MantidAbout : public MantidQt::API::MantidDialog
	{

	Q_OBJECT
	
	public:
		/// constructor
		MantidAbout(QWidget* parent = 0);
		///destructor
		~MantidAbout() {}
	private:
		/// form generated by QT Designer
		Ui::MantidAbout m_uiForm;
	};


#endif /* MANTIDABOUT_H */