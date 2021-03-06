set ( SRC_FILES
    src/AlgorithmSelectorWidget.cpp
    src/CatalogHelper.cpp
    src/CatalogSearch.cpp
    src/CatalogSelector.cpp
    src/CheckboxHeader.cpp
    src/ColorBarWidget.cpp
    src/DataSelector.cpp
    src/DiagResults.cpp
    src/DisplayCurveFit.cpp
    src/DoubleDialogEditor.cpp
    src/DoubleSpinBox.cpp
    src/ErrorCurve.cpp
    src/FilenameDialogEditor.cpp
    src/FindDialog.cpp
    src/FindReplaceDialog.cpp
    src/FitPropertyBrowser.cpp
    src/FormulaDialogEditor.cpp
    src/FitOptionsBrowser.cpp
    src/FunctionBrowser.cpp
    src/HintingLineEdit.cpp
    src/InstrumentSelector.cpp
    src/IndirectInstrumentConfig.cpp
    src/InputController.cpp
    src/MantidHelpWindow.cpp
    src/MWDiag.cpp
    src/MWRunFiles.cpp
    src/MWView.cpp
    src/MessageDisplay.cpp
    src/MultifitSetupDialog.cpp
    src/MuonFitPropertyBrowser.cpp
    src/MuonSequentialFitDialog.cpp
    src/PeakPicker.cpp
    src/PeriodicTableWidget.cpp
    src/PreviewPlot.cpp
    src/ProcessingAlgoWidget.cpp
    src/PropertyHandler.cpp
    src/RangeSelector.cpp
    src/RenameParDialog.cpp
    src/SafeQwtPlot.cpp
    src/SaveWorkspaces.cpp
    src/ScriptEditor.cpp
    src/SelectFunctionDialog.cpp
    src/SelectWorkspacesDialog.cpp
    src/SequentialFitDialog.cpp
    src/SlicingAlgorithmDialog.cpp
    src/SlitCalculator.cpp
    src/StringDialogEditor.cpp
    src/StringEditorFactory.cpp
    src/UserFunctionDialog.cpp
    src/WorkspaceEditorFactory.cpp
    src/WorkspaceSelector.cpp
    src/pqHelpWindow.cxx
    src/pythonCalc.cpp
    src/LineEditWithClear.cpp
    src/InstrumentView/BinDialog.cpp
    src/InstrumentView/CollapsiblePanel.cpp
    src/InstrumentView/ColorMapWidget.cpp
    src/InstrumentView/CompAssemblyActor.cpp
    src/InstrumentView/ComponentActor.cpp
    src/InstrumentView/DetXMLFile.cpp
    src/InstrumentView/GLActor.cpp
    src/InstrumentView/GLActorCollection.cpp
    src/InstrumentView/GLActorVisitor.cpp
    src/InstrumentView/GLColor.cpp
    src/InstrumentView/GLObject.cpp
    src/InstrumentView/ICompAssemblyActor.cpp
    src/InstrumentView/InstrumentActor.cpp
    src/InstrumentView/InstrumentTreeModel.cpp
    src/InstrumentView/InstrumentTreeWidget.cpp
    src/InstrumentView/InstrumentWidget.cpp
    src/InstrumentView/InstrumentWidgetMaskTab.cpp
    src/InstrumentView/InstrumentWidgetPickTab.cpp
    src/InstrumentView/InstrumentWidgetRenderTab.cpp
    src/InstrumentView/InstrumentWidgetTab.cpp
    src/InstrumentView/InstrumentWidgetTreeTab.cpp
    src/InstrumentView/MantidGLWidget.cpp
    src/InstrumentView/MaskBinsData.cpp
    src/InstrumentView/ObjCompAssemblyActor.cpp
    src/InstrumentView/ObjComponentActor.cpp
    src/InstrumentView/OneCurvePlot.cpp
    src/InstrumentView/OpenGLError.cpp
    src/InstrumentView/PanelsSurface.cpp
    src/InstrumentView/PeakMarker2D.cpp
    src/InstrumentView/PeakOverlay.cpp
    src/InstrumentView/Projection3D.cpp
    src/InstrumentView/ProjectionSurface.cpp
    src/InstrumentView/RectangularDetectorActor.cpp
    src/InstrumentView/StructuredDetectorActor.cpp
    src/InstrumentView/RectF.cpp
    src/InstrumentView/RotationSurface.cpp
    src/InstrumentView/SampleActor.cpp
    src/InstrumentView/Shape2D.cpp
    src/InstrumentView/Shape2DCollection.cpp
    src/InstrumentView/SimpleWidget.cpp
    src/InstrumentView/UCorrectionDialog.cpp
    src/InstrumentView/UnwrappedCylinder.cpp
    src/InstrumentView/UnwrappedSphere.cpp
    src/InstrumentView/UnwrappedSurface.cpp
    src/InstrumentView/Viewport.cpp
    src/InstrumentView/XIntegrationControl.cpp
)

# Header files with Q_OBJECT that qmake will "moc"
set ( MOC_FILES 
    inc/MantidQtMantidWidgets/AlgorithmSelectorWidget.h
    inc/MantidQtMantidWidgets/CheckboxHeader.h
    inc/MantidQtMantidWidgets/ColorBarWidget.h
    inc/MantidQtMantidWidgets/DataSelector.h
    inc/MantidQtMantidWidgets/DiagResults.h
    inc/MantidQtMantidWidgets/DisplayCurveFit.h
    inc/MantidQtMantidWidgets/DoubleDialogEditor.h
    inc/MantidQtMantidWidgets/DoubleSpinBox.h
    inc/MantidQtMantidWidgets/FilenameDialogEditor.h
    inc/MantidQtMantidWidgets/FormulaDialogEditor.h
    inc/MantidQtMantidWidgets/FindReplaceDialog.h
    inc/MantidQtMantidWidgets/FindDialog.h
    inc/MantidQtMantidWidgets/FitPropertyBrowser.h
    inc/MantidQtMantidWidgets/FitOptionsBrowser.h
    inc/MantidQtMantidWidgets/FunctionBrowser.h
    inc/MantidQtMantidWidgets/HintingLineEdit.h
    inc/MantidQtMantidWidgets/CatalogSearch.h
    inc/MantidQtMantidWidgets/CatalogSelector.h
    inc/MantidQtMantidWidgets/InstrumentSelector.h
    inc/MantidQtMantidWidgets/IndirectInstrumentConfig.h
    inc/MantidQtMantidWidgets/InputController.h
    inc/MantidQtMantidWidgets/MantidHelpWindow.h
    inc/MantidQtMantidWidgets/MessageDisplay.h
    inc/MantidQtMantidWidgets/MultifitSetupDialog.h
    inc/MantidQtMantidWidgets/MuonFitPropertyBrowser.h
    inc/MantidQtMantidWidgets/MuonSequentialFitDialog.h
    inc/MantidQtMantidWidgets/MWDiag.h
    inc/MantidQtMantidWidgets/MWRunFiles.h
    inc/MantidQtMantidWidgets/MWView.h
    inc/MantidQtMantidWidgets/PeakPicker.h
    inc/MantidQtMantidWidgets/PeriodicTableWidget.h
    inc/MantidQtMantidWidgets/pqHelpWindow.h
    inc/MantidQtMantidWidgets/PreviewPlot.h
    inc/MantidQtMantidWidgets/PropertyHandler.h
    inc/MantidQtMantidWidgets/ProcessingAlgoWidget.h
    inc/MantidQtMantidWidgets/pythonCalc.h
    inc/MantidQtMantidWidgets/RangeSelector.h
    inc/MantidQtMantidWidgets/RenameParDialog.h
    inc/MantidQtMantidWidgets/SafeQwtPlot.h
    inc/MantidQtMantidWidgets/SaveWorkspaces.h
    inc/MantidQtMantidWidgets/ScriptEditor.h
    inc/MantidQtMantidWidgets/SelectFunctionDialog.h
    inc/MantidQtMantidWidgets/SelectWorkspacesDialog.h
    inc/MantidQtMantidWidgets/SequentialFitDialog.h
    inc/MantidQtMantidWidgets/SlicingAlgorithmDialog.h
    inc/MantidQtMantidWidgets/SlitCalculator.h
    inc/MantidQtMantidWidgets/StringDialogEditor.h
    inc/MantidQtMantidWidgets/StringEditorFactory.h
    inc/MantidQtMantidWidgets/UserFunctionDialog.h
    inc/MantidQtMantidWidgets/WorkspaceEditorFactory.h
    inc/MantidQtMantidWidgets/WorkspaceSelector.h
    inc/MantidQtMantidWidgets/LineEditWithClear.h
    inc/MantidQtMantidWidgets/InstrumentView/BinDialog.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentTreeModel.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentTreeWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetRenderTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetPickTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetMaskTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetTreeTab.h
    inc/MantidQtMantidWidgets/InstrumentView/ColorMapWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/MantidGLWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/OneCurvePlot.h
    inc/MantidQtMantidWidgets/InstrumentView/CollapsiblePanel.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentActor.h
    inc/MantidQtMantidWidgets/InstrumentView/PeakOverlay.h
    inc/MantidQtMantidWidgets/InstrumentView/ProjectionSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/Projection3D.h
    inc/MantidQtMantidWidgets/InstrumentView/Shape2DCollection.h
    inc/MantidQtMantidWidgets/InstrumentView/UnwrappedSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/UCorrectionDialog.h
	inc/MantidQtMantidWidgets/InstrumentView/XIntegrationControl.h
)

# Add the include files are NOT already in MOC_FILES
set ( INC_FILES
    ${MOC_FILES}
    inc/MantidQtMantidWidgets/AlgorithmHintStrategy.h
    inc/MantidQtMantidWidgets/CatalogHelper.h
    inc/MantidQtMantidWidgets/ErrorCurve.h
    inc/MantidQtMantidWidgets/WidgetDllOption.h
    inc/MantidQtMantidWidgets/HintStrategy.h
    inc/MantidQtMantidWidgets/InstrumentView/BinDialog.h
    inc/MantidQtMantidWidgets/InstrumentView/CollapsiblePanel.h
    inc/MantidQtMantidWidgets/InstrumentView/ColorMapWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/CompAssemblyActor.h
    inc/MantidQtMantidWidgets/InstrumentView/ComponentActor.h
    inc/MantidQtMantidWidgets/InstrumentView/DetXMLFile.h
    inc/MantidQtMantidWidgets/InstrumentView/GLActor.h
    inc/MantidQtMantidWidgets/InstrumentView/GLActorCollection.h
    inc/MantidQtMantidWidgets/InstrumentView/GLActorVisitor.h
    inc/MantidQtMantidWidgets/InstrumentView/GLColor.h
    inc/MantidQtMantidWidgets/InstrumentView/GLObject.h
    inc/MantidQtMantidWidgets/InstrumentView/ICompAssemblyActor.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentActor.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentTreeModel.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentTreeWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetMaskTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetPickTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetRenderTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetTreeTab.h
    inc/MantidQtMantidWidgets/InstrumentView/InstrumentWidgetTypes.h
    inc/MantidQtMantidWidgets/InstrumentView/MantidGLWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/MaskBinsData.h
    inc/MantidQtMantidWidgets/InstrumentView/ObjCompAssemblyActor.h
    inc/MantidQtMantidWidgets/InstrumentView/ObjComponentActor.h
    inc/MantidQtMantidWidgets/InstrumentView/OneCurvePlot.h
    inc/MantidQtMantidWidgets/InstrumentView/OpenGLError.h
    inc/MantidQtMantidWidgets/InstrumentView/PanelsSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/PeakMarker2D.h
    inc/MantidQtMantidWidgets/InstrumentView/PeakOverlay.h
    inc/MantidQtMantidWidgets/InstrumentView/Projection3D.h
    inc/MantidQtMantidWidgets/InstrumentView/ProjectionSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/RectangularDetectorActor.h
    inc/MantidQtMantidWidgets/InstrumentView/StructuredDetectorActor.h
    inc/MantidQtMantidWidgets/InstrumentView/RectF.h
    inc/MantidQtMantidWidgets/InstrumentView/RotationSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/SampleActor.h
    inc/MantidQtMantidWidgets/InstrumentView/Shape2D.h
    inc/MantidQtMantidWidgets/InstrumentView/Shape2DCollection.h
    inc/MantidQtMantidWidgets/InstrumentView/SimpleWidget.h
    inc/MantidQtMantidWidgets/InstrumentView/UCorrectionDialog.h
    inc/MantidQtMantidWidgets/InstrumentView/UnwrappedCylinder.h
    inc/MantidQtMantidWidgets/InstrumentView/UnwrappedSphere.h
    inc/MantidQtMantidWidgets/InstrumentView/UnwrappedSurface.h
    inc/MantidQtMantidWidgets/InstrumentView/Viewport.h
    inc/MantidQtMantidWidgets/InstrumentView/XIntegrationControl.h
)

# QtDesigner UI files to process
set ( UI_FILES 
    inc/MantidQtMantidWidgets/DataSelector.ui
    inc/MantidQtMantidWidgets/CatalogSearch.ui
    inc/MantidQtMantidWidgets/CatalogSelector.ui
    inc/MantidQtMantidWidgets/ColorBarWidget.ui
    inc/MantidQtMantidWidgets/DisplayCurveFit.ui
    inc/MantidQtMantidWidgets/IndirectInstrumentConfig.ui
    inc/MantidQtMantidWidgets/MWDiag.ui
    inc/MantidQtMantidWidgets/MWRunFiles.ui
    inc/MantidQtMantidWidgets/MWView.ui
    inc/MantidQtMantidWidgets/MultifitSetupDialog.ui
    inc/MantidQtMantidWidgets/MuonSequentialFitDialog.ui
    inc/MantidQtMantidWidgets/ProcessingAlgoWidget.ui
    inc/MantidQtMantidWidgets/RenameParDialog.ui
    inc/MantidQtMantidWidgets/SelectFunctionDialog.ui
    inc/MantidQtMantidWidgets/SequentialFitDialog.ui
    inc/MantidQtMantidWidgets/SlicingAlgorithmDialog.ui
    inc/MantidQtMantidWidgets/SlitCalculator.ui
    inc/MantidQtMantidWidgets/UserFunctionDialog.ui
    inc/MantidQtMantidWidgets/PeriodicTableWidget.ui
    inc/MantidQtMantidWidgets/PreviewPlot.ui
    inc/MantidQtMantidWidgets/pqHelpWindow.ui
    inc/MantidQtMantidWidgets/InstrumentView/UCorrectionDialog.ui
)

# Python unit tests
set ( TEST_PY_FILES
  test/MWRunFilesTest.py
)

set ( TEST_FILES
  AlgorithmHintStrategyTest.h
)

find_package (Qt4 REQUIRED QtHelp QtWebKit QtNetwork QUIET)
include(${QT_USE_FILE})

include_directories ( ../../QtPropertyBrowser/src )
include_directories ( ../../MantidQt/API/inc )
include_directories ( ../../MantidQt/MantidWidgets/inc )
include_directories ( ../../MantidQt/SliceViewer/inc )
include_directories ( ../../MantidQt/SpectrumViewer/inc )
include_directories ( ../../MantidQt/Factory/inc )
qt4_wrap_cpp ( MOCCED_FILES ${MOC_FILES} )

set ( ALL_SRC ${SRC_FILES} ${MOCCED_FILES} )

qt4_wrap_ui ( UI_HDRS ${UI_FILES} ) 
qt4_add_resources ( RES_FILES ${PROJECT_SOURCE_DIR}/images/MantidWidgets.qrc )

include_directories ( ${QSCINTILLA_INCLUDE_DIR} )
add_definitions ( -DQSCINTILLA_DLL )     # Will only have an effect on Windows (as is desired)


# For Windows:
add_definitions ( -DIN_MANTIDQT_MANTIDWIDGETS -DQT_QTPROPERTYBROWSER_IMPORT )
# Use a precompiled header where they are supported
enable_precompiled_headers( inc/MantidQtMantidWidgets/PrecompiledHeader.h ALL_SRC )
add_library ( MantidWidgets ${ALL_SRC} ${INC_FILES} ${UI_HDRS} ${RES_FILES} )

if (OSX_VERSION VERSION_GREATER 10.8)
  set_target_properties ( MantidWidgets PROPERTIES INSTALL_RPATH "@loader_path/../MacOS")
endif ()

target_link_libraries ( MantidWidgets LINK_PRIVATE ${TCMALLOC_LIBRARIES_LINKTIME} 
            MantidQtAPI
            QtPropertyBrowser
            ${CORE_MANTIDLIBS}
            ${QT_LIBRARIES} ${QWT_LIBRARIES}
            ${QSCINTILLA_LIBRARIES}
            ${POCO_LIBRARIES}
            ${Boost_LIBRARIES}
            ${ZLIB_LIBRARIES}
            ${OPENGL_LIBRARIES}
            )

###########################################################################
# Testing
###########################################################################

add_subdirectory ( test )

###########################################################################
# Installation settings
###########################################################################

install ( TARGETS MantidWidgets ${SYSTEM_PACKAGE_TARGET} DESTINATION ${LIB_DIR} )
