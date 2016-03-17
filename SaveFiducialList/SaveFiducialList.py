import os
import unittest
import numpy
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *

#
# SaveFiducialList
#

class SaveFiducialList(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "SaveFiducialList" # TODO make this more human readable by adding spaces
    self.parent.categories = ["TestModules"]
    self.parent.contributors = ["Davide Scorza (Vicomtech)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    Module to save a fiducial list into a txt file, saving only the translation coordinates x, y, z.
    The file will be generated in a predefined folder and the text file named with the name of the fiducial list.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# CompareDistancesWidget

class SaveFiducialListWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)   

    # Markup selector
    self.markupSelectorLabel = qt.QLabel()
    self.markupSelectorLabel.setText( "Markup list: " )
    self.markupSelector = slicer.qMRMLNodeComboBox()
    self.markupSelector.nodeTypes = ( "vtkMRMLMarkupsFiducialNode", "" )
    self.markupSelector.noneEnabled = False
    self.markupSelector.selectNodeUponCreation = True
    self.markupSelector.setMRMLScene( slicer.mrmlScene )
    self.markupSelector.setToolTip( "Pick the markup list to save" )
    parametersFormLayout.addRow(self.markupSelectorLabel, self.markupSelector)    

    self.savePathLabel = qt.QLabel()
    self.savePathLabel.setText("Choose the saving Path: [Actually is FIXED, to modify]")
    parametersFormLayout.addRow(self.savePathLabel)
    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Save")
    self.applyButton.toolTip = "Save the fiducials in a text file"
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)    
    
    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.markupSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)
    

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.markupSelector.currentNode() 

  def onApplyButton(self):
    #slicer.app.processEvents()
    self.logic = SaveFiducialListLogic()
    print("Run the algorithm")
    self.logic.run(self.markupSelector.currentNode())

#
# SaveFiducialListLogic
#

class SaveFiducialListLogic(ScriptedLoadableModuleLogic):

  def run(self,markupSelector):
    """
    Run the actual algorithm
    """
    self.info={}
    
    FiducialNumber = markupSelector.GetNumberOfFiducials()
    #FidMatrix = numpy.ndarray((FiducialNumber, 3))
    FidCoords = numpy.empty(3)
    
    #save these metrics in a txt file
    
    savePath = "C:/MyProjects/FiducialLists/" + markupSelector.GetName()
    print savePath
    f = open(savePath + ".txt", "w")
    
    for fidIndex in xrange(FiducialNumber):
        markupSelector.GetNthFiducialPosition(fidIndex,FidCoords)
        f.write(str(FidCoords[0]) + "\t" + str(FidCoords[1]) + "\t" + str(FidCoords[2]) + "\n")
        
    f.close()
    self.delayDisplay('List Saved\n')
    
    return True