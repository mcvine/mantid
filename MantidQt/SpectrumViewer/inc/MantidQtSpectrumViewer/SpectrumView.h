#ifndef  SPECTRUM_VIEW_H
#define  SPECTRUM_VIEW_H

#include <QMainWindow>
#include <QMdiSubWindow>
#include <QtGui>
#include <QList>

#include "MantidAPI/MatrixWorkspace_fwd.h"
#include "MantidQtAPI/WorkspaceObserver.h"
#include "MantidQtSpectrumViewer/GraphDisplay.h"
#include "MantidQtSpectrumViewer/SpectrumDataSource.h"
#include "MantidQtSpectrumViewer/MatrixWSDataSource.h"
#include "MantidQtSpectrumViewer/DllOptionSV.h"

/**
    @class SpectrumView

    This is the QMainWindow for the SpectrumView data viewer.  Data is
    displayed in an SpectrumView, by constructing the SpectrumView object and
    specifying a particular data source.

    @author Dennis Mikkelson
    @date   2012-04-03

    Copyright © 2012 ORNL, STFC Rutherford Appleton Laboratories

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

    Code Documentation is available at
                 <http://doxygen.mantidproject.org>
 */

namespace Ui
{
class SpectrumViewer; // Forward declaration of UI file
}

namespace MantidQt
{
namespace SpectrumView
{

// Forward declarations
class EModeHandler;
class RangeHandler;
class SliderHandler;
class SpectrumDisplay;
class SVConnections;
class MatrixWSDataSource;

class EXPORT_OPT_MANTIDQT_SPECTRUMVIEWER SpectrumView : public QMainWindow, public MantidQt::API::WorkspaceObserver
{
  Q_OBJECT

public:
  /// Construct a SpectrumView to display data from the specified data source
  SpectrumView( QWidget * parent = 0 );

  ~SpectrumView() override;
  void renderWorkspace(Mantid::API::MatrixWorkspace_const_sptr wksp);
  QList<boost::shared_ptr<SpectrumDisplay>> getSpectrumDisplays() const { return m_spectrumDisplay; }
  bool isTrackingOn() const;

signals:
  void spectrumDisplayChanged(SpectrumDisplay*);

protected slots:
  void closeWindow();
  void changeSpectrumDisplay(int tab);
  void respondToTabCloseReqest(int tab);
  void changeTracking(bool on);

protected:
  void resizeEvent(QResizeEvent *event) override;
  void
  preDeleteHandle(const std::string &wsName,
                  const boost::shared_ptr<Mantid::API::Workspace> ws) override;
  void afterReplaceHandle(
      const std::string &wsName,
      const boost::shared_ptr<Mantid::API::Workspace> ws) override;

  void dragMoveEvent(QDragMoveEvent *de) override;
  void dragEnterEvent(QDragEnterEvent *de) override;
  void dropEvent(QDropEvent *de) override;

private:
  void updateHandlers();
  void loadSettings();
  void saveSettings() const;

  QList<MatrixWSDataSource_sptr> m_dataSource;
  QList<boost::shared_ptr<SpectrumDisplay>> m_spectrumDisplay;
  boost::shared_ptr<GraphDisplay> m_hGraph;
  boost::shared_ptr<GraphDisplay> m_vGraph;
  boost::shared_ptr<SVConnections> m_svConnections;


  Ui::SpectrumViewer *m_ui;
  SliderHandler      *m_sliderHandler;
  RangeHandler       *m_rangeHandler;
  EModeHandler       *m_emodeHandler;

signals:
  void needToClose();
  void needToUpdate();
};

} // namespace SpectrumView
} // namespace MantidQt

#endif   // SPECTRUM_VIEW_H
