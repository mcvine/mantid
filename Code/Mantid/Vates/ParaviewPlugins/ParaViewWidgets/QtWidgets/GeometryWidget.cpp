#include "GeometryWidget.h"
#include "DimensionWidget.h"
#include "MantidVatesAPI/GeometryPresenter.h"
#include <QLabel>
#include <QGridLayout>

GeometryWidget::GeometryWidget(Mantid::VATES::GeometryPresenter* pPresenter, bool readOnlyLimits) : m_widgetFactory(readOnlyLimits), m_pPresenter(pPresenter) 
{
  QGridLayout* headerLayout = new QGridLayout();
  QVBoxLayout* bodyLayout = new QVBoxLayout();
  
  headerLayout->addWidget(new QLabel("Geometry"), 0, 0, 1, 2, Qt::AlignCenter); 
  bodyLayout->addLayout(headerLayout);
  
  this->setLayout(bodyLayout);
  m_pPresenter->acceptView(this);
}

GeometryWidget::~GeometryWidget()
{
  delete m_pPresenter;
}

void GeometryWidget::addDimensionView(Mantid::VATES::DimensionView* dimView)
{
  DimensionWidget* dimWidget = dynamic_cast<DimensionWidget*>(dimView); //TODO. design should not need capability queries!
  if(dimWidget != NULL)
  {
    QLayout* layout = this->layout();
    layout->addWidget(dimWidget);
  }
}

std::string GeometryWidget::getGeometryXMLString() const
{
  return m_pPresenter->getGeometryXML();
}

const Mantid::VATES::DimensionViewFactory& GeometryWidget::getDimensionViewFactory()
{
  return m_widgetFactory;
}

void GeometryWidget::raiseModified()
{
  emit valueChanged();
}
