import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

# vtkModel
class vtkModel(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "vtkModel" # TODO make this more human readable by adding spaces
    self.parent.categories = ["TestModules"]
    self.parent.dependencies = []
    self.parent.contributors = ["Davide Scorza"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# vtkModelWidget
#

class vtkModelWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Parameters Area
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    
    # Apply Button
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Generate Model"
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Add vertical spacer
    self.layout.addStretch(1)

  def cleanup(self):
    pass

  def onApplyButton(self):
    logic = vtkModelLogic()
    logic.run()

# vtkModelLogic
class vtkModelLogic(ScriptedLoadableModuleLogic):
  """
  """

  def run(self):
    """
    Run the actual algorithm
    """
    #generating Nodes for displaying a new model
    modelNode = slicer.vtkMRMLModelNode()
    dispNode  = slicer.vtkMRMLModelDisplayNode()
    transform = slicer.vtkMRMLLinearTransformNode()
    
    #Display node characteristics
    dispNode.SetVisibility(True)
    dispNode.SetSliceIntersectionVisibility(True)
    dispNode.SetOpacity(1)
    dispNode.SetColor(1, 1, 0)
    dispNode.SetScene(slicer.mrmlScene)

    #generate sphere data
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(10,10,10)
    sphere.SetRadius(40)
    sphere.Update()
    
    #adding necessary nodes to the Scene
    slicer.mrmlScene.AddNode(dispNode)
    slicer.mrmlScene.AddNode(transform)
    slicer.mrmlScene.AddNode(modelNode)
    
    #model node name and associations!
    modelNode.SetName("SphereModelNode")
    modelNode.SetScene(slicer.mrmlScene)
    modelNode.SetAndObserveTransformNodeID(transform.GetID())
    modelNode.SetAndObserveDisplayNodeID(dispNode.GetID())

    apd = vtk.vtkAppendPolyData()
    apd.AddInputData(sphere.GetOutput())
    apd.Update()
    
    #adding model node poly data! Here there are  sphere's data!!!    
    modelNode.SetAndObservePolyData(apd.GetOutput())
    #update the scene
    slicer.mrmlScene.Modified()
    

    return True