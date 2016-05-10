import os
import unittest
import vtk, qt, ctk, slicer, numpy
from slicer.ScriptedLoadableModule import *
import logging

#
# MIPonPlane
#

class MIPonPlane(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "MIPonPlane" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Segmentation"]
    self.parent.dependencies = []
    self.parent.contributors = ["Giuseppe De Luca (POLIMI)"]
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# MIPonPlaneWidget
#

class MIPonPlaneWidget(ScriptedLoadableModuleWidget):
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

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene(slicer.mrmlScene)
    self.inputSelector.setToolTip("Pick the input to the algorithm.")
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # fiducials
    #
    self.FiducialsSelector = slicer.qMRMLNodeComboBox()
    self.FiducialsSelector.nodeTypes = (("vtkMRMLMarkupsFiducialNode"), "")
    self.FiducialsSelector.selectNodeUponCreation = True
    self.FiducialsSelector.addEnabled = False
    self.FiducialsSelector.removeEnabled = False
    self.FiducialsSelector.noneEnabled = False
    self.FiducialsSelector.showHidden = False
    self.FiducialsSelector.showChildNodeTypes = False
    self.FiducialsSelector.setMRMLScene(slicer.mrmlScene)
    self.FiducialsSelector.setToolTip("Pick the fiducials.")
    parametersFormLayout.addRow("Fiducials: ", self.FiducialsSelector)


    #
    # Set Slice value
    #
    self.numberOfSlicesWidget = ctk.ctkSliderWidget()
    self.numberOfSlicesWidget.singleStep = 1
    self.numberOfSlicesWidget.minimum = 1
    self.numberOfSlicesWidget.maximum = 1000
    self.numberOfSlicesWidget.value = 50
    self.numberOfSlicesWidget.setToolTip("Set the number of slices.")
    parametersFormLayout.addRow("Number of Slices", self.numberOfSlicesWidget)

    #
    # Set SpacingFraction
    #
    self.spacingFractionWidget = ctk.ctkSliderWidget()
    self.spacingFractionWidget.singleStep = 0.05
    self.spacingFractionWidget.minimum = 0.1
    self.spacingFractionWidget.maximum = 1
    self.spacingFractionWidget.value = 0.5
    self.spacingFractionWidget.setToolTip(
      "Set Spacing Fraction to reduce interpolation artifacts")
    parametersFormLayout.addRow("Spacing Fraction", self.spacingFractionWidget)

    # Compute Distance Button
    self.distanceButton = qt.QPushButton("Compute Distance")
    self.distanceButton.toolTip = "Compute Distance between F1 and F2"
    self.distanceButton.enabled = False
    parametersFormLayout.addRow(self.distanceButton)
    self.distanceButton.connect('clicked(bool)', self.ondistanceButton)


    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    #
    # Apply MIP
    #
    self.applyMIP = qt.QPushButton("Apply MIP")
    self.applyMIP.toolTip = "MIP algorithm."
    self.applyMIP.enabled = True
    parametersFormLayout.addRow(self.applyMIP)
    self.applyMIP.connect('clicked(bool)', self.onApplyMIP)

    #
    # Remove MIP
    #
    self.removeMIP = qt.QPushButton("Remove MIP")
    self.removeMIP.toolTip = "MIP algorithm."
    self.removeMIP.enabled = True
    parametersFormLayout.addRow(self.removeMIP)
    self.removeMIP.connect('clicked(bool)', self.onRemoveMIP)

    #
    # Export MIP Button
    #
    self.ExportMIPButton = qt.QPushButton("Export MIP Image")
    self.ExportMIPButton.toolTip = "Export MIP for MATLAB."
    self.ExportMIPButton.enabled = False
    parametersFormLayout.addRow(self.ExportMIPButton)
    self.ExportMIPButton.connect('clicked(bool)', self.onExportMIPButton)

    # Results
    self.distanceValueLabel = qt.QLabel()
    parametersFormLayout.addRow("Distance between F1 and F2: ", self.distanceValueLabel)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)


    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.FiducialsSelector.currentNode()
    self.distanceButton.enabled = self.FiducialsSelector.currentNode()
    self.ExportMIPButton.enabled = self.inputSelector.currentNode()

  def onApplyButton(self):

    logic = MIPonPlaneLogic()
    self.distanceValue, self.F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText('%.3f' % self.distanceValue + ' millimeters')
    self.MRT = logic.run(self.F)

  def ondistanceButton(self):
    logic = MIPonPlaneLogic()
    self.distanceValue, F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText('%.3f' % self.distanceValue + ' millimeters')

  def onApplyMIP(self):

    logic = MIPonPlaneLogic()
    Slab = int(self.numberOfSlicesWidget.value)
    SpacingFraction = self.spacingFractionWidget.value
    self.reslice = logic.MipOnPlane(Slab, SpacingFraction)

  def onRemoveMIP(self):

    logic = MIPonPlaneLogic()
    logic.RemoveMIP()

  def onExportMIPButton(self):
    logic = MIPonPlaneLogic()
    logic.ExportMIP(self.inputSelector.currentNode(), self.MRT)


#
# MIPonPlaneLogic
#

class MIPonPlaneLogic(ScriptedLoadableModuleLogic):
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
    print 'Distance Value between fiducials Computed: ', distanceValue

    if distanceValue < 0:
      slicer.util.delayDisplay("Error Computing Distance")
      return 1
    else:
      return distanceValue, F

  def run(self, F):
    """
    Run the actual algorithm
    """

    logging.info('Processing started')

    P0 = F[0]
    P1 = F[1]
    # New Support Point
    P2 = [P0[0] - 20, P0[1] + 10, P0[2] + 30]

    # Compute the 3 new axis of LCS (Local Coordinate System)
    Az = P1 - P0
    Az = Az / numpy.linalg.norm(Az)
    Ax = numpy.cross(Az, P2 - P0)
    Ax = Ax / numpy.linalg.norm(Ax)
    Ay = numpy.cross(Az, Ax)

    # Rototranslation Matrix
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

    return MRT

    # # TODO Check why is not working correctly
    # if det == 1:
    #   logging.info('Processing completed')
    # else:
    #   logging.info("Error computing the rotation")
    #   slicer.util.delayDisplay("Error computing the rotation")
    #   return False

  def createNewLinearTransform(self, numpyMatrix):

    linearTransformNode = slicer.vtkMRMLLinearTransformNode()
    transformMatrix = vtk.vtkMatrix4x4()
    numpyMatrix = numpyMatrix.ravel()
    numpyMatrix = numpyMatrix.squeeze()

    transformMatrix.DeepCopy([numpyMatrix[0], numpyMatrix[1], numpyMatrix[2], numpyMatrix[3],
                              numpyMatrix[4], numpyMatrix[5], numpyMatrix[6], numpyMatrix[7],
                              numpyMatrix[8], numpyMatrix[9], numpyMatrix[10], numpyMatrix[11],
                              numpyMatrix[12], numpyMatrix[13], numpyMatrix[14], numpyMatrix[15]])

    TransfNode = slicer.util.getNode('TransformToLCS')
    if not TransfNode:
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

   #UNUSED#

  def RotateVolume(self, linearTransformNode):

    scene = slicer.mrmlScene
    inode = scene.GetNodeByID('vtkMRMLScalarVolumeNode1')
    event = vtk.vtkIntArray()
    event.InsertNextValue(slicer.vtkMRMLTransformableNode.TransformModifiedEvent)
    inode.SetAndObserveNodeReferenceID('transform', 'vtkMRMLLinearTransformNode4', event)


  def MipOnPlane(self, Slab, SpacingFraction):


    sliceNode = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
    appLogic = slicer.app.applicationLogic()
    sliceLogic = appLogic.GetSliceLogic(sliceNode)
    sliceLayerLogic = sliceLogic.GetBackgroundLayer()
    reslice = sliceLayerLogic.GetReslice()
    reslice.SetSlabModeToMax()
    reslice.SetSlabNumberOfSlices(Slab)
    reslice.SetSlabSliceSpacingFraction(SpacingFraction)
    sliceNode.Modified()
    return reslice



  def RemoveMIP(self):


    sliceNode = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
    appLogic = slicer.app.applicationLogic()
    sliceLogic = appLogic.GetSliceLogic(sliceNode)
    sliceLayerLogic = sliceLogic.GetBackgroundLayer()
    reslice = sliceLayerLogic.GetReslice()
    reslice.SetSlabModeToMean() #default value
    reslice.SetSlabNumberOfSlices(1) #default value
    reslice.SetSlabSliceSpacingFraction(1) #default value
    sliceNode.Modified()


  def ExportMIP(self, inputVolume,  MRT):
    """
    This function is to create a new volume containing the Xray and display in the Slicer Scene
    """

    #extract input volume information

    inputVolumeData = inputVolume.GetImageData()

    #reslice Filter
    resliceFilter = vtk.vtkImageReslice()
    resliceFilter.SetInputData(inputVolumeData)

    # #The Output dimensionality need to be fixed on '2' in order to extract a single slice.
    resliceFilter.SetOutputDimensionality(2)
    Spacing = inputVolume.GetSpacing()
    resliceFilter.SetOutputSpacing(Spacing[0],Spacing[1],Spacing[2]) # FIXME recuperare spacing del volume in input
    resliceFilter.SetSlabModeToMax()
    resliceFilter.SetSlabNumberOfSlices(50)
    resliceFilter.SetSlabSliceSpacingFraction(0.5)
    # #set reslices axes origin in order to get the slice I want
    resliceFilter.SetResliceAxesOrigin(MRT[:3, 3])
    # #direction cosine for standard axial slice
    resliceFilter.SetResliceAxesDirectionCosines(MRT[:3, 0],MRT[:3, 1],MRT[:3, 2])

    # #Update the filter
    resliceFilter.Update()

    #Create a Volume and fill it with the resliceFilter Output
    Origin = inputVolume.GetOrigin()
    self.volumeNode = slicer.vtkMRMLScalarVolumeNode()
    self.volumeNode.SetSpacing(Spacing[0],Spacing[1],Spacing[2])
    self.volumeNode.SetImageDataConnection(resliceFilter.GetOutputPort())
    #self.volumeNode.SetIJKToRASMatrix(transformMatrix)
    self.volumeNode.SetOrigin(Origin)
    # self.imageData = self.volumeNode.GetImageData()
    # self.imageData.SetExtent(inputVolumeExtent[0], inputVolumeExtent[1], inputVolumeExtent[2], inputVolumeExtent[3], 0, 0)

    # Add volume to scene
    slicer.mrmlScene.AddNode(self.volumeNode)
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    self.volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    self.volumeNode.CreateDefaultStorageNode()
    self.volumeNode.SetName('MIP')
    self.name = self.volumeNode.GetName()

    print "MIP extracted"
