# coding=utf-8
import os
import unittest
import vtk, qt, ctk, slicer, numpy
from slicer.ScriptedLoadableModule import *
import logging

#
# surgeonEyeView
#

class surgeonEyeView(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "surgeonEyeView" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Segmentation"]
    self.parent.dependencies = []
    self.parent.contributors = ["Giuseppe De Luca (POLIMI)"]
    self.parent.helpText = """
    This module will be used to verify the correctness of the trajectories implemented.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# surgeonEyeViewWidget
#

class surgeonEyeViewWidget(ScriptedLoadableModuleWidget):
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

    # fiducials
    self.FiducialsSelector = slicer.qMRMLNodeComboBox()
    self.FiducialsSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.FiducialsSelector.selectNodeUponCreation = True
    self.FiducialsSelector.addEnabled = False
    self.FiducialsSelector.removeEnabled = False
    self.FiducialsSelector.noneEnabled = False
    self.FiducialsSelector.showHidden = False
    self.FiducialsSelector.showChildNodeTypes = False
    self.FiducialsSelector.setMRMLScene( slicer.mrmlScene )
    self.FiducialsSelector.setToolTip( "Pick the fiducials." )
    parametersFormLayout.addRow("Fiducials: ", self.FiducialsSelector)


    # Compute Distance Button
    self.distanceButton = qt.QPushButton("Compute Distance")
    self.distanceButton.toolTip = "Compute Distance between F1 and F2"
    self.distanceButton.enabled = False
    parametersFormLayout.addRow(self.distanceButton)
    self.distanceButton.connect('clicked(bool)', self.ondistanceButton)

    # Apply Button
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Results
    self.distanceValueLabel = qt.QLabel()
    parametersFormLayout.addRow("Distance between F1 and F2: ", self.distanceValueLabel)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.FiducialsSelector.currentNode()
    self.distanceButton.enabled = self.FiducialsSelector.currentNode()

  def onApplyButton(self):
    pass

    logic = surgeonEyeViewLogic()
    self.distanceValue, self.F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText('%.3f' % self.distanceValue + ' millimeters')
    logic.run(self.F)


  def ondistanceButton(self):

    logic = surgeonEyeViewLogic()
    self.distanceValue, F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText('%.3f' % self.distanceValue + ' millimeters')




#
# surgeonEyeViewLogic
#

class surgeonEyeViewLogic(ScriptedLoadableModuleLogic):


  def calcDistance(self, Fiducials):

    F = numpy.zeros((2, 3))
    FidCoords = numpy.empty(3)
    fiducialCount = Fiducials.GetNumberOfFiducials()
    if fiducialCount == 0:
      print 'No points found'
    elif fiducialCount == 1:
      print 'There is only one point'

    for x in range(fiducialCount):
      Fiducials.GetNthFiducialPosition(x, FidCoords)
      F[x] = FidCoords
    x = numpy.power(F[0, 0] - F[1, 0], 2)
    y = numpy.power(F[0, 1] - F[1, 1], 2)
    z = numpy.power(F[0, 2] - F[1, 2], 2)

    distanceValue = numpy.sqrt((x + y + z))
    print 'Distance Value between fiducials Computed: ',distanceValue

    if distanceValue < 0:
      slicer.util.delayDisplay("Error Computing Distance")
      return 1
    else:
      return distanceValue, F


  def run(self, F ):
    """
    Run the actual algorithm
    """

    logging.info('Processing started')

    P0 = F[0]
    P1 = F[1]
    # New Support Point
    P2 = [P0[0]-20, P0[1]+10,P0[2]+30]


    # Compute the 3 new axis of LCS (Local Coordinate System)
    Az = P1 - P0
    Az = Az / numpy.linalg.norm(Az)
    Ax = numpy.cross(Az, P2 - P0)
    Ax = Ax / numpy.linalg.norm(Ax)
    Ay = numpy.cross(Az, Ax)

    #Rototranslation Matrix
    MRT = numpy.zeros((4, 4))
    MRT[:3, 0] = Ax
    MRT[:3, 1] = Ay
    MRT[:3, 2] = Az
    MRT[:3, 3] = P0
    MRT[3, :3] = 0
    MRT[3, 3] = 1

    det = numpy.linalg.det(MRT)
    print "Determinant: ", det

    self.createNewLinearTransform(MRT)

    #TODO Check why is not working correctly
    if det == 1:
      logging.info('Processing completed')
    else:
      logging.info("Error computing the rotation")
      slicer.util.delayDisplay("Error computing the rotation")
      return False

  def createNewLinearTransform(self, numpyMatrix):

    linearTransformNode = slicer.vtkMRMLLinearTransformNode()
    transformMatrix = vtk.vtkMatrix4x4()
    numpyMatrix = numpyMatrix.ravel()
    numpyMatrix = numpyMatrix.squeeze()

    transformMatrix.DeepCopy([numpyMatrix[0], numpyMatrix[1], numpyMatrix[2], numpyMatrix[3],
                              numpyMatrix[4], numpyMatrix[5], numpyMatrix[6], numpyMatrix[7],
                              numpyMatrix[8], numpyMatrix[9], numpyMatrix[10], numpyMatrix[11],
                              numpyMatrix[12], numpyMatrix[13], numpyMatrix[14], numpyMatrix[15]])


    slicer.mrmlScene.AddNode(linearTransformNode)
    linearTransformNode.SetAndObserveMatrixTransformToParent(transformMatrix)
    linearTransformNode.SetName('TransformToLCS')
    self.moveSliceToNewReferenceFrame(transformMatrix)
    return transformMatrix


  def moveSliceToNewReferenceFrame(self, transformMatrix):

    redSlice = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
    redSlice.SetSliceToRAS(transformMatrix)
    redSlice.UpdateMatrices()
    lm = slicer.app.layoutManager()
    lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

  def RotateVolume(self, linearTransformNode):

    scene = slicer.mrmlScene
    inode = scene.GetNodeByID('vtkMRMLScalarVolumeNode1')
    event = vtk.vtkIntArray()
    event.InsertNextValue(slicer.vtkMRMLTransformableNode.TransformModifiedEvent)
    inode.SetAndObserveNodeReferenceID('transform', 'vtkMRMLLinearTransformNode4', event)
