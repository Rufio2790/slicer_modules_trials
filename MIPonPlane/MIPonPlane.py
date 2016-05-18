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
    # self.VTKtransform = logic.numpyMatrixtoVTKtransform(self.MRT)
    # logic.slice(self.inputSelector.currentNode(), 'sagittal', 80, self.VTKtransform)

  def ondistanceButton(self):
    logic = MIPonPlaneLogic()
    self.distanceValue, F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText('%.3f' % self.distanceValue + ' millimeters')

  def onApplyMIP(self):

    logic = MIPonPlaneLogic()
    Slab = int(self.numberOfSlicesWidget.value)
    SpacingFraction = self.spacingFractionWidget.value
    self.reslice = logic.MipOnPlane(Slab, SpacingFraction)
    #Origin = self.inputSelector.currentNode.GetOutputOrigin()
    #print Origin


  def onRemoveMIP(self):

    logic = MIPonPlaneLogic()
    logic.RemoveMIP()

  def onExportMIPButton(self):
    logic = MIPonPlaneLogic()
    self.distanceValue, self.F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.MRT = logic.run(self.F)
    self.VTKtransform = logic.numpyMatrixtoVTKtransform(self.MRT)
    logic.slice(self.inputSelector.currentNode(), 'oblique', 80, self.VTKtransform)

    #logic.ExportMIP()
    #Origin = self.inputSelector.currentNode().GetOrigin()
    #logic.resliceVolume(self.inputSelector.currentNode(), self.VTKtransform)
    #logic.slice(self.inputSelector.currentNode(), 'oblique', 50)

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
    #lm = slicer.app.layoutManager()
    #lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpRedSliceView)

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
    # imageData = reslice.GetOutput()
    #self.resliceVolume(input)
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

  def ExctactSlice(self, inputVolume):

    inputVolumeData = inputVolume.GetImageData()




  def ExportMIP(self):
    """
    This function is to create a new volume containing the Xray and display in the Slicer Scene
    """

    #extract input volume information

    # inputVolumeData = inputVolume.GetImageData()
    #
    # #reslice Filter
    # resliceFilter = vtk.vtkImageReslice()
    # resliceFilter.SetInputData(inputVolumeData)
    #
    # # #The Output dimensionality need to be fixed on '2' in order to extract a single slice.
    # resliceFilter.SetOutputDimensionality(2)
    # Spacing = inputVolume.GetSpacing()
    # resliceFilter.SetOutputSpacing(Spacing[0],Spacing[1],Spacing[2])
    # resliceFilter.SetSlabModeToMax()
    # resliceFilter.SetSlabNumberOfSlices(50)
    # resliceFilter.SetSlabSliceSpacingFraction(0.5)
    # # # #set reslices axes origin in order to get the slice I want
    # resliceFilter.SetResliceAxesOrigin(MRT[:3, 3])
    # # #direction cosine for standard axial slice
    # resliceFilter.SetResliceAxesDirectionCosines(MRT[:3, 0],MRT[:3, 1],MRT[:3, 2])
    # #
    # # # #Update the filter
    # resliceFilter.Update()
    #
    # #Create a Volume and fill it with the resliceFilter Output
    # Origin = inputVolume.GetOrigin()
    # self.volumeNode = slicer.vtkMRMLScalarVolumeNode()
    # self.volumeNode.SetSpacing(Spacing[0],Spacing[1],Spacing[2])
    # inputVolume.SetImageDataConnection(resliceFilter.GetOutputPort())
    # #self.volumeNode.SetIJKToRASMatrix(transformMatrix)
    # self.volumeNode.SetOrigin(Origin)
    # # self.imageData = self.volumeNode.GetImageData()
    # # self.imageData.SetExtent(inputVolumeExtent[0], inputVolumeExtent[1], inputVolumeExtent[2], inputVolumeExtent[3], 0, 0)
    #
    #
    # # Add volume to scene
    # slicer.mrmlScene.AddNode(self.volumeNode)
    # displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    # slicer.mrmlScene.AddNode(displayNode)
    # colorNode = slicer.util.getNode('Grey')
    # displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    # self.volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    # self.volumeNode.CreateDefaultStorageNode()
    # self.volumeNode.SetName('MIP')
    # self.name = self.volumeNode.GetName()
    #
    # print "MIP extracted"


    layoutName = 'Red'
    imagePathPattern = '/tmp/Mip.png'
    widget = slicer.app.layoutManager().sliceWidget(layoutName)
    view = widget.sliceView()
    logic = widget.sliceLogic()
    sliceNode = slicer.mrmlScene.GetNodeByID('vtkMRMLSliceNodeRed')
    offset = sliceNode.GetSliceOffset()
    logic.SetSliceOffset(offset)
    image = qt.QPixmap.grabWidget(view).toImage()
    image.save(imagePathPattern)


  def volumeNodeFromVolArray(self, volArray, shape, vtkDataType, dimensions, volName, spacing):

      volArray = volArray.reshape(shape)

      # make a vtkImage data of the vol data
      image = vtk.vtkImageData()
      image.SetDimensions(dimensions)
      image.AllocateScalars(vtk.VTK_FLOAT, 1)
      imageArray = vtk.util.numpy_support.vtk_to_numpy(image.GetPointData().GetScalars()).reshape(shape)
      imageArray[:] = volArray

      # make a slicer volume node for the image
      volumeNode = slicer.vtkMRMLScalarVolumeNode()
      volumeNode.SetName(volName)
      volumeNode.SetSpacing(*spacing)
      volumeNode.SetAndObserveImageData(image)
      slicer.mrmlScene.AddNode(volumeNode)
      volumeNode.CreateDefaultDisplayNodes()

      # make this volume visible
      applicationLogic = slicer.app.applicationLogic()
      selectionNode = applicationLogic.GetSelectionNode()
      selectionNode.SetReferenceActiveVolumeID( volumeNode.GetID() )
      applicationLogic.PropagateVolumeSelection(1)

      return volumeNode

  def resliceVolume(self, VolumeNode, transformMatrix):

      image = VolumeNode.GetImageData()
      resliceFilter = vtk.vtkImageReslice()
      resliceFilter.SetInputData(image)

      # The Output dimensionality need to be fixed on '2' in order to extract a single slice.
      resliceFilter.SetOutputDimensionality(2)
      #resliceFilter.SetOutputSpacing(1.0, 1.0, 1.0)
      # set the output extent to 1024,1024,1
      #inputVolumeExtent = image.GetExtent()
      #resliceFilter.SetOutputExtent(inputVolumeExtent[0], inputVolumeExtent[1], inputVolumeExtent[2], inputVolumeExtent[3], 0, 0)
      # set reslices axes origin in order to get the slice I want
      Origin = image.GetOrigin()
      print Origin
      #resliceFilter.SetResliceAxesOrigin(Origin[0],Origin[1],Origin[2]+80)
      #resliceFilter.SetResliceAxesOrigin(MRT[:3, 3])
      # direction cosine for standard axial slice
      teta=90
      #resliceFilter.SetResliceAxesDirectionCosines(1.0, 0.0, 0.0, 0.0, numpy.cos(teta), -numpy.sin(teta), 0.0, numpy.sin(teta), numpy.cos(teta))
      #resliceFilter.SetResliceAxesDirectionCosines(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
      #resliceFilter.SetResliceAxesDirectionCosines(MRT[:3, 0],MRT[:3, 1],MRT[:3, 2])
      (xMin, xMax, yMin, yMax, zMin, zMax) = image.GetExtent()
      (xSpacing, ySpacing, zSpacing) = image.GetSpacing()
      (x0, y0, z0) = image.GetOrigin()

      center = [x0 + xSpacing * 0.5 * (xMin + xMax),
                y0 + ySpacing * 0.5 * (yMin + yMax),
                z0 + zSpacing * 0.5 * (zMin + zMax)]

      sagittal = vtk.vtkMatrix4x4()
      sagittal.DeepCopy((0, 0, -1, center[0],
                         1, 0, 0, center[1],
                         0, -1, 0, center[2],
                         0, 0, 0, 1))
      resliceFilter.SetResliceAxes(sagittal)

      resliceFilter.SetSlabModeToMax()
      resliceFilter.SetSlabNumberOfSlices(10)
      print transformMatrix
      #resliceFilter.SetResliceAxes(transformMatrix)
      #resliceFilter.SetInterpolationModeToLinear()
      resliceFilter.SetSlabModeToMax()
      resliceFilter.SetSlabNumberOfSlices(50)
      resliceFilter.SetSlabSliceSpacingFraction(0.5)

      # Update the filter
      resliceFilter.Update()

      OutImage = vtk.vtkImageData()
      OutImage = resliceFilter.GetOutput()
      volumeNode = slicer.vtkMRMLScalarVolumeNode()
      volumeNode.SetName('ResVolume')
      volumeNode.SetAndObserveImageData(OutImage)
      slicer.mrmlScene.AddNode(volumeNode)
      volumeNode.CreateDefaultDisplayNodes()

  def RotateresliceVolume(self, VolumeNode, MRT, Origin):

        image = VolumeNode.GetImageData()
        resliceFilter = vtk.vtkImageReslice()
        resliceFilter.SetInputData(image)

        # The Output dimensionality need to be fixed on '2' in order to extract a single slice.
        resliceFilter.SetOutputDimensionality(2)
        #resliceFilter.SetOutputSpacing(1.0, 1.0, 1.0)
        # set the output extent to 1024,1024,1
        # resliceFilter.SetOutputExtent(inputVolumeExtent[0], inputVolumeExtent[1], inputVolumeExtent[2], inputVolumeExtent[3], 0, 0)
        # set reslices axes origin in order to get the slice I want
        #resliceFilter.SetResliceAxesOrigin(Origin[0],Origin[1],Origin[2]+50)
        # direction cosine for standard axial slice
        (xMin, xMax, yMin, yMax, zMin, zMax) = image.GetExtent()
        (xSpacing, ySpacing, zSpacing) = image.GetSpacing()
        (x0, y0, z0) = image.GetOrigin()

        center = [x0 + xSpacing * 0.5 * (xMin + xMax),
                  y0 + ySpacing * 0.5 * (yMin + yMax),
                  z0 + zSpacing * 0.5 * (zMin + zMax)]

        sagittal = vtk.vtkMatrix4x4()
        sagittal.DeepCopy((0, 0, -1, center[0],
                           1, 0, 0, center[1],
                           0, -1, 0, center[2],
                           0, 0, 0, 1))
        resliceFilter.SetResliceAxesDirectionCosines(sagittal)
        resliceFilter.SetSlabModeToMax()
        resliceFilter.SetSlabNumberOfSlices(10)

        # Update the filter
        resliceFilter.Update()

        OutImage = vtk.vtkImageData()
        OutImage = resliceFilter.GetOutput()
        volumeNode = slicer.vtkMRMLScalarVolumeNode()
        volumeNode.SetName('ResVolume')
        volumeNode.SetAndObserveImageData(OutImage)
        slicer.mrmlScene.AddNode(volumeNode)
        volumeNode.CreateDefaultDisplayNodes()

  def numpyMatrixtoVTKtransform(self, numpyMatrix):
        transformMatrix = vtk.vtkMatrix4x4()
        numpyMatrix = numpyMatrix.ravel()
        numpyMatrix = numpyMatrix.squeeze()

        transformMatrix.DeepCopy([numpyMatrix[0], numpyMatrix[1], numpyMatrix[2], numpyMatrix[3],
                                  numpyMatrix[4], numpyMatrix[5], numpyMatrix[6], numpyMatrix[7],
                                  numpyMatrix[8], numpyMatrix[9], numpyMatrix[10], numpyMatrix[11],
                                  numpyMatrix[12], numpyMatrix[13], numpyMatrix[14], numpyMatrix[15]])


        return transformMatrix

  def slice(self, Volume, Orientation, slice, VTKtransform):
      # reader is the input VTK volume
      # Calculate the center of the volume
      ImageSet = Volume.GetImageData()
      (xMin, xMax, yMin, yMax, zMin, zMax) = ImageSet.GetExtent()
      Spacing = Volume.GetSpacing()
      Origin = Volume.GetOrigin()

      center = [Origin[0] + Spacing[0] * 0.5 * (xMin + xMax),
                Origin[1] + Spacing[1] * 0.5 * (yMin + yMax),
                Origin[2] + Spacing[2] * 0.5 * (zMin + zMax)]
      print center


      # Extract a slice in the desired orientation
      #reslice.SetBlendModeToMax()
      #reslice.SetResliceAxesOrigin(10, 20, 20)
      if Orientation == 'axial':

          axial = vtk.vtkMatrix4x4()
          axial.DeepCopy((1, 0, 0, center[0],
                          0, 1, 0, center[1],
                          0, 0, 1, slice,
                          0, 0, 0, 1))

          AXreslice = vtk.vtkImageReslice()
          AXreslice.SetInterpolationModeToLinear()
          AXreslice.SetInputData(ImageSet)
          AXreslice.SetOutputDimensionality(2)
          AXreslice.SetResliceAxes(axial)
          AXreslice.SetSlabModeToMax()
          AXreslice.SetSlabNumberOfSlices(50)
          AXreslice.SetSlabSliceSpacingFraction(0.5)
          AXreslice.Update()
          self.ViewReslice(AXreslice,'axial', Spacing, Origin)

      elif Orientation == 'sagittal':

          sagittal = vtk.vtkMatrix4x4()
          sagittal.DeepCopy((0, 0, -1, slice,
                             1, 0, 0, center[1],
                             0, -1, 0, center[2],
                             0, 0, 0, 1))

          SAreslice = vtk.vtkImageReslice()
          SAreslice.SetInterpolationModeToLinear()
          SAreslice.SetInputData(ImageSet)
          SAreslice.SetOutputDimensionality(2)
          SAreslice.SetResliceAxes(sagittal)
          SAreslice.SetSlabModeToMax()
          SAreslice.SetSlabNumberOfSlices(50)
          SAreslice.SetSlabSliceSpacingFraction(0.5)
          SAreslice.Update()
          self.ViewReslice(SAreslice,'sagittal', Spacing, Origin)

      elif Orientation == 'coronal':

          coronal = vtk.vtkMatrix4x4()
          coronal.DeepCopy((1, 0, 0, center[0],
                            0, 0, 1, slice,
                            0, -1, 0, center[2],
                            0, 0, 0, 1))

          COreslice = vtk.vtkImageReslice()
          COreslice.SetInterpolationModeToLinear()
          COreslice.SetInputData(ImageSet)
          COreslice.SetOutputDimensionality(2)
          COreslice.SetResliceAxes(coronal)
          COreslice.SetSlabModeToMax()
          COreslice.SetSlabNumberOfSlices(50)
          COreslice.SetSlabSliceSpacingFraction(0.5)
          COreslice.Update()
          self.ViewReslice(COreslice,'coronal', Spacing, Origin)

      elif Orientation == 'oblique':
          matrix = VTKtransform
          oblique = vtk.vtkMatrix4x4()
          oblique.DeepCopy((1, 0, 0, center[0], # matrix.GetElement(0,3)
                            0, 0.866025404, -0.5, center[1], #matrix.GetElement(1,3)
                            0, 0.5, 0.866025404, center[2],  #matrix.GetElement(2,3)
                            0, 0, 0, 1))
          OBreslice = vtk.vtkImageReslice()
          OBreslice.SetInterpolationModeToLinear()
          OBreslice.SetInputData(ImageSet)
          OBreslice.SetOutputDimensionality(2)
          # d = vtk.vtkMatrixToHomogeneousTransform() #Provo a convertire da una vtk4X4 ad un'abstract transform per poi
          # d.SetInput(oblique)                        #usare SetResliceTransform
          # oblique = d.MakeTransform()
          OBreslice.SetResliceAxes(oblique)
          OBreslice.SetSlabModeToMax()
          OBreslice.SetSlabNumberOfSlices(50)
          OBreslice.SetSlabSliceSpacingFraction(0.5)
          OBreslice.Update()
          self.ViewReslice(OBreslice, 'oblique', Spacing, Origin)

      #reslice.AutoCropOutputOn()
      # reslice.SetSlabModeToMax()
      # reslice.SetSlabNumberOfSlices(50)
      # reslice.SetSlabSliceSpacingFraction(0.5)
      # reslice.Update()
      # NewImage = reslice.GetOutput()
      # volumeNode = slicer.vtkMRMLScalarVolumeNode()
      # volumeNode.SetName('ResVolume2')
      # volumeNode.SetAndObserveImageData(NewImage)
      # slicer.mrmlScene.AddNode(volumeNode)
      # volumeNode.CreateDefaultDisplayNodes()



  def ViewReslice (self, reslice, Orientation, Spacing, Origin):

      if Orientation == 'axial':
        Image = reslice.GetOutput()
        volumeNode = slicer.vtkMRMLScalarVolumeNode()
        volumeNode.SetName('AxialSlice')
        volumeNode.SetAndObserveImageData(Image)
        volumeNode.SetSpacing(Spacing[0],Spacing[1],1.0)
        volumeNode.SetOrigin(Origin)
        slicer.mrmlScene.AddNode(volumeNode)
        volumeNode.CreateDefaultDisplayNodes()

      elif Orientation == 'coronal':
          Image = reslice.GetOutput()
          volumeNode = slicer.vtkMRMLScalarVolumeNode()
          volumeNode.SetName('CoronalSlice')
          volumeNode.SetAndObserveImageData(Image)
          volumeNode.SetSpacing(Spacing[0],Spacing[2], 1.0)
          volumeNode.SetOrigin(Origin)
          slicer.mrmlScene.AddNode(volumeNode)
          volumeNode.CreateDefaultDisplayNodes()

      elif Orientation == 'sagittal':
          Image = reslice.GetOutput()
          volumeNode = slicer.vtkMRMLScalarVolumeNode()
          volumeNode.SetName('SagittalSlice')
          volumeNode.SetAndObserveImageData(Image)
          volumeNode.SetSpacing(Spacing[1],Spacing[2],1.0)
          volumeNode.SetOrigin(Origin)
          slicer.mrmlScene.AddNode(volumeNode)
          volumeNode.CreateDefaultDisplayNodes()

      elif Orientation == 'oblique':
          Image = reslice.GetOutput()
          volumeNode = slicer.vtkMRMLScalarVolumeNode()
          volumeNode.SetName('ObliqueSlice')
          volumeNode.SetAndObserveImageData(Image)
          volumeNode.SetSpacing(Spacing[0], Spacing[1], 1.0)
          volumeNode.SetOrigin(Origin)
          slicer.mrmlScene.AddNode(volumeNode)
          volumeNode.CreateDefaultDisplayNodes()