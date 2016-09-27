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
    self.spacingFractionWidget.setToolTip("Set Spacing Fraction to reduce interpolation artifacts")
    parametersFormLayout.addRow("Spacing Fraction", self.spacingFractionWidget)

    # Compute Distance Button
    self.distanceButton = qt.QPushButton("Compute Distance")
    self.distanceButton.toolTip = "Compute Distance between EP and TP"
    self.distanceButton.enabled = False
    parametersFormLayout.addRow(self.distanceButton)
    self.distanceButton.connect('clicked(bool)', self.ondistanceButton)


    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply Roto-Translation")
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
    parametersFormLayout.addRow("Distance between EP and TP: ", self.distanceValueLabel)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    self.segmentationSection = SegmentationParametersWidget(self)

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
    Slab = int(self.numberOfSlicesWidget.value)
    SpacingFraction = self.spacingFractionWidget.value
    self.reslice = logic.MipOnPlane(Slab, SpacingFraction)
    logic.ExctactSlice(self.reslice)


#
# MIPonPlaneLogic
#


class SegmentationParametersWidget(object):

  def __init__(self, parent):

    self.parent = parent

    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)












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

    v1 = F[1] - F[0]  # electrode direction = z (from entry to target)
    v2 = [v1[0], - v1[1], 0]  # add initialization to maintain visibilty in red view

    Az = v1 / numpy.linalg.norm(v1)
    Ax = v2 - numpy.sum(v2 * Az) * Az
    Ax = Ax / numpy.linalg.norm(Ax)
    Ay = numpy.cross(Az, Ax)

    # Rototranslation Matrix
    MRT = numpy.zeros((4, 4))
    MRT[:3, 0] = Ax
    MRT[:3, 1] = Ay
    MRT[:3, 2] = Az
    MRT[:3, 3] = F[0]
    MRT[3, :3] = 0
    MRT[3, 3] = 1

    det = numpy.linalg.det(MRT)
    print "Determinant: ", det

    self.createNewLinearTransform(MRT)

    return MRT

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

  def ExctactSlice(self, reslice):

      Image = reslice.GetOutput()
      volumeNode = slicer.vtkMRMLScalarVolumeNode()
      volumeNode.SetName('MIPSlice')
      volumeNode.SetAndObserveImageData(Image)
      slicer.mrmlScene.AddNode(volumeNode)
      volumeNode.CreateDefaultDisplayNodes()

      imm = slicer.util.array('MIPSlice')
      img = numpy.array(imm.squeeze())
      #img = numpy.fliplr(img)

      row_range_img = numpy.ceil(img.shape[1])
      col_range_img = numpy.ceil(img.shape[0])

      center_x = int(col_range_img / 2)
      center_y = int(row_range_img / 2)
      offset = 20
      size_crop = 2 * offset
      crop_img = img[center_x - offset : center_x + offset, center_y - offset : center_y + offset]

      no_MRF = GMM_Segmentation(crop_img, nCluster = 3, likelihoodTreshold = 1, maxIter = 100, beta = 1.5, nNeighbors = None)
      seg_no_MRF = no_MRF.compute()
      seg_no_MRF_vol = numpy.zeros([1, size_crop, size_crop])
      seg_no_MRF_vol[0, :, :] = seg_no_MRF
      vassel_vol = numpy.zeros([1, size_crop, size_crop])
      vassel_vol[0, :, :] = crop_img
      dd = [size_crop, size_crop, 1]

      vasVolumeNode = self.volumeNodeFromVolArray(imm, img.shape, None, [img.shape[1], img.shape[0], 1], 'Vassel_tot', volumeNode.GetSpacing())

      dimension = [no_MRF.nCol, no_MRF.nRow, 1]

      not_segemntedVolumeNode = self.volumeNodeFromVolArray(vassel_vol, crop_img.shape, None, dd, 'Vassel_crop',
                                                            volumeNode.GetSpacing())
      segmentedVolumeNode = self.volumeNodeFromVolArray(seg_no_MRF_vol, crop_img.shape, None, dimension, 'segNoMRF',
                                                        not_segemntedVolumeNode.GetSpacing())


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


  def numpyMatrixtoVTKtransform(self, numpyMatrix):
        transformMatrix = vtk.vtkMatrix4x4()
        numpyMatrix = numpyMatrix.ravel()
        numpyMatrix = numpyMatrix.squeeze()

        transformMatrix.DeepCopy([numpyMatrix[0], numpyMatrix[1], numpyMatrix[2], numpyMatrix[3],
                                  numpyMatrix[4], numpyMatrix[5], numpyMatrix[6], numpyMatrix[7],
                                  numpyMatrix[8], numpyMatrix[9], numpyMatrix[10], numpyMatrix[11],
                                  numpyMatrix[12], numpyMatrix[13], numpyMatrix[14], numpyMatrix[15]])


        return transformMatrix

  def slice(self, Volume, Orientation, slice, VTKtransform=None):
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
          #matrix = VTKtransform
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



class GMM_Segmentation (object):
    """
    img = np array, nCluster = n classi

    """
    def __init__(self, img, nCluster, likelihoodTreshold=None, maxIter= None, beta = None, nNeighbors = None):
        self.img = img #np array
        self.nCluster = nCluster #int, numb classes
        self.beta = beta #MRF parameter
        self.nNeighbors = nNeighbors #4 or 8, MRF paraneter
        self.MRFflag = False #no MRF computation
        self.exit_conditions = False

        #likelihood and n max iteration initialization in case of None
        if likelihoodTreshold is None:
            self.likelihoodTreshold = 130
            print('Absent loglikelihood value, computation with default value'), self.likelihoodTreshold
        else:
            self.likelihoodTreshold = likelihoodTreshold
        if maxIter is None:
            self.maxIter = 300
            print('Absent number max iterations, computation with default value'), self.maxIter
        else:
            self.maxIter = maxIter
        if (self.beta is not None or self.nNeighbors is not None):
            self.MRFflag = True #avvio MRF in caso di chiamata con almeno uno dei 2 parametri caratteristici
            if self.beta is None:
                self.beta = 1.5
                print('Absent beta value, computation with default value'), self.beta
            else:
                self.beta = beta
            if self.nNeighbors is None:
                self.nNeighbors = 4
                print('Absent number of neighbors, computation with default value'), self.nNeighbors
            else:
                self.nNeighbors = nNeighbors


        #initialization to enter in the first Exp-Max step
        self.loglikelihood = 1000
        self.diff_loglikelihood = 20000
        self.iter = 0
        self.seg_img = img
        self.inMRF = False

        #IMAGES PARAMETERS
        self.minInt = numpy.amin(self.img)
        self.maxInt = numpy.amax(self.img)
        self.nRow = self.img.shape[0]
        self.nCol = self.img.shape[1]
        self.nPixel = self.nRow * self.nCol

        #INITIALIZATION (mean, variance)
        #mu = (1 : nCluster) * maxVal / (nCluster + 1);
        self.mu = (numpy.arange(self.nCluster) + numpy.ones(self.nCluster)) * (self.maxInt / (self.nCluster + 1))
        #vari = ones(1, nCluster) * maxVal;
        self.vari = numpy.ones([self.nCluster]) * self.maxInt

        # INITIALIZATION (probmatrix, sum of each pixel distr ~ sum_gauss, proportion of each cluster ~ pr)
        self.prob_matrix = numpy.zeros([self.nRow, self.nRow, self.nCluster])
        self.sum_gauss = numpy.zeros([self.nRow, self.nCol])
        self.pr = numpy.ones([self.nCluster]) / self.nCluster
        self.distr = numpy.zeros([self.nRow, self.nCol, self.nCluster])

    #CALLING HERE THE PROPER FUNCTION
    def compute(self):
        self.seg_img = self.expectation_maximization_steps()
        if self.MRFflag is True:
            self.inMRF = True
            self.exit_conditions = False
            self.loglikelihood = 20000
            self.seg_img = self.expectation_maximization_steps()
            return self.seg_img
        else:
            return self.seg_img

    #EXPECTATION MAXIMIZATION STEPS (BOTH SIMPLE ~inMRF=False~ AND MRF ~inMRF=True)
    def expectation_maximization_steps (self):
        while(self.exit_conditions == False):
            self.iter += 1
            if self.inMRF is True:
                self.matrix_prior = self.matrix_prior_computation()

            #EXPECTATION STEP
            loglikelihood_old = self.loglikelihood
            self.distr = self.distr_computation()
            self.sum_gauss = self.sum_gauss_computation()
            self.loglikelihood = self.loglikelihood_computation()
            self.prob_matrix = self.prob_matrix_computation()
            self.diff_loglikelihood = abs(loglikelihood_old - self.loglikelihood)
            #print 'Diff_loglikelihood', self.diff_loglikelihood
            #MAXIMIZATION STEP
            self.pr = self.pr_computation()
            self.mu = self.mu_comp()
            self.vari = self.vari_comp()
            #UPDATE EXIT CONDITIONS
            self.exit_conditions = self.exit_conditions_computation()
        else:
            if self.iter == self.maxIter:
                print('ERROR: Max iter')
                return -1
            return self.seg_img

    def exit_conditions_computation(self):   #VERIFY EXIT CONDITIONS
        if self.diff_loglikelihood <= self.likelihoodTreshold:
            print '\nConvergenza raggiunta in ', self.iter, ' iterazioni\n'
            self.seg_img = numpy.argmax(self.prob_matrix, 2)
            self.exit_conditions = True
        elif self.iter == self.maxIter:
            #display di errore per maxIter raggiunto
            self.exit_conditions = True
        else:
            self.exit_conditions = False
        return self.exit_conditions

    def distr_computation(self):
        #costn(l, c, z) = pr(z) / sqrt(2 * pi * vari(z));
        #espon(l, c, z) = -((img(l, c) - mu(z)). ^ 2 / (2 * vari(z)));
        #distr(l, c, z) = costn(l, c, z) * exp(espon(l, c, z));
        costn = numpy.zeros([self.nRow, self.nCol, self.nCluster])
        espon = numpy.zeros([self.nRow, self.nCol, self.nCluster])

        for l in range(self.nRow):
            for c in range(self.nCol):
                for z in range(self.nCluster):
                    if self.inMRF is False:
                        costn[l, c, z] = self.pr[z] / numpy.sqrt(2 * numpy.pi * self.vari[z])
                    else:
                        costn[l, c, z] = self.matrix_prior[l, c, z] / numpy.sqrt(2 * numpy.pi *self.vari[z])
                    espon[l, c, z] = -((self.img[l, c] - self.mu[z])**2 / (2 * self.vari[z]) + numpy.exp(-13))
                    self.distr[l, c, z] = costn[l, c, z] * numpy.exp(espon[l, c, z])
        return self.distr

    def sum_gauss_computation(self):
        # somma lungo la terza direzione (0 righe, 1 colonne, 2 terza dim)
        self.sum_gauss = numpy.sum(self.distr, axis = 2) + numpy.exp(-13)
        #sum_gauss: nR * nC
        return self.sum_gauss

    def prob_matrix_computation(self):
        prob_matrix = numpy.zeros([self.nRow, self.nCol, self.nCluster])
        #self.sum_gauss = self.sum_gauss_computation()
        for l in range(self.nRow):
            for c in range(self.nCol):
                for z in range(self.nCluster):
                    #prob_matrix(l, c, z) = distr(l, c, z). / sum_gauss(l, c);
                    prob_matrix[l, c, z] = self.distr[l, c, z] / self.sum_gauss[l, c]
        return prob_matrix

    def loglikelihood_computation(self):
        return (numpy.sum(numpy.sum(numpy.log(self.sum_gauss))))

    def pr_computation(self):
        pr_sum_lines = numpy.sum(self.prob_matrix_computation(), 0)
        pr_sum = numpy.sum(pr_sum_lines, 0)
        pr = pr_sum / self.nPixel
        pr = pr / numpy.sum(pr)
        return pr

    def mu_comp(self):
        s_row = numpy.sum(self.prob_matrix, 0)
        s_row_line = numpy.sum(s_row, 0)
        for cl in range(self.nCluster):
            self.mu[cl] = numpy.sum(numpy.sum(self.img[:, :] * self.prob_matrix[:, :, cl])) / s_row_line[cl]
        return self.mu

    def vari_comp(self):
        s_row = numpy.sum(self.prob_matrix, 0)
        s_row_line = numpy.sum(s_row, 0)
        for cl in range(self.nCluster):
            stnde = numpy.zeros(self.img.shape)
            stnde[:, :] = self.img[:, :] - self.mu[cl]
            self.vari[cl] = + numpy.exp(-13) + numpy.sum(numpy.sum(stnde[:, :]**2 * self.prob_matrix[:, :, cl])) / s_row_line[cl]
        return self.vari

    #da chiamare all'inizio o a exp-max completata
    def original_intensities(self):
        x_intensities = numpy.arange(self.minInt, self.maxInt + 1)
        y_intensities = numpy.zeros(x_intensities.shape)
        for intensity in range(self.minInt, self.maxInt + 1):
            if intensity != 0:
                array = (self.img == intensity)
                y_intensities[intensity - self.minInt] = sum(sum(array))
        #x_intensities contiene i valori di ogni intensita' nel range dell'immagine
        #y_intensities[i] contiene il numero di pixel per ogni intensita'
        y_intensities = y_intensities / self.nPixel
        #y_intensities contiene ora la frequenza di ogni intensita'
        return [x_intensities, y_intensities]

    #MRF only
    def matrix_prior_computation(self):
        self.matrix_prior = numpy.zeros([self.nRow, self.nCol, self.nCluster])
        for row in range(self.nRow):
            for col in range(self.nCol):
                neighborhood = self.find_neighborhood(row, col)
                cl_frequence = numpy.zeros(self.nCluster)
                for nb in range(neighborhood.size):
                    for cl in range(self.nCluster):
                        if neighborhood[nb] == cl:
                            cl_frequence[cl] += 1
                energy_fun = self.beta * cl_frequence
                self.matrix_prior[row, col, :] = numpy.exp(energy_fun) / numpy.sum(numpy.exp(energy_fun))
        return self.matrix_prior

    #MRF only
    def find_neighborhood(self, row, col):
        upleft = (row == 0 and col == 0)
        upright = (row == 0 and col == self.nCol - 1)
        downright = (row == self.nRow - 1 and col == 0)
        downleft = (row == self.nRow - 1 and col == self.nCol - 1)
        upline = (row == 0)
        downline = (row == self.nRow - 1)
        leftline = (col == 0)
        rightline = (col == self.nCol - 1)
        if (self.nNeighbors == 4 or self.nNeighbors ==8):
            if (upleft or upright or downright or downleft):
                neighborhood = numpy.zeros(2)
                if upleft:
                    neighborhood[0] = self.seg_img[0, 1]
                    neighborhood[1] = self.seg_img[1, 0]
                elif upright:
                    neighborhood[0] = self.seg_img[0, self.nCol - 2]
                    neighborhood[1] = self.seg_img[1, self.nCol - 1]
                elif downleft:
                    neighborhood[0] = self.seg_img[self.nRow - 2, 0]
                    neighborhood[1] = self.seg_img[self.nRow - 1, 1]
                elif downright:
                    neighborhood[0] = self.seg_img[self.nRow - 2, self.nCol - 1]
                    neighborhood[1] = self.seg_img[self.nRow - 1, self.nCol - 2]
            elif (upline, downline, leftline, rightline):
                neighborhood = numpy.zeros(3)
                if upline:
                    neighborhood[0] = self.seg_img[0, col - 1]
                    neighborhood[1] = self.seg_img[0, col + 1]
                    neighborhood[2] = self.seg_img[1, col]
                elif downline:
                    neighborhood[0] = self.seg_img[self.nRow - 1, col - 1]
                    neighborhood[1] = self.seg_img[self.nRow - 1, col + 1]
                    neighborhood[2] = self.seg_img[self.nRow - 2, col]
                elif leftline:
                    neighborhood[0] = self.seg_img[row + 1, 0]
                    neighborhood[1] = self.seg_img[row - 1, 0]
                    neighborhood[2] = self.seg_img[row, 1]
                elif rightline:
                    neighborhood[0] = self.seg_img[row + 1, self.nCol - 1]
                    neighborhood[1] = self.seg_img[row - 1, self.nCol - 1]
                    neighborhood[2] = self.seg_img[row, self.nCol - 2]
            else:
                neighborhood = numpy.zeros(4)
                neighborhood[0] = self.seg_img(row, col + 1)
                neighborhood[1] = self.seg_img(row, col - 1)
                neighborhood[2] = self.seg_img(row + 1, col)
                neighborhood[3] = self.seg_img(row - 1, col)

            if self.nNeighbors == 8:
                if upleft:
                    extra = numpy.array(self.seg_img[1, 1])
                elif upright:
                    extra = numpy.array(self.seg_img[1, self.nCol - 2])
                elif downleft:
                    extra = numpy.array(self.seg_img[self.nRow - 2, 1])
                elif downright:
                    extra = numpy.array(self.seg_img[self.nRow - 2, self.nCol - 2])
                elif upline:
                    extra = numpy.array([self.seg_img[1, col + 1], self.seg_img[1, col - 1]])
                elif downline:
                    extra = numpy.array([self.seg_img[self.nRow - 2, col + 1], self.seg_img[self.nRow - 2, col - 1]])
                elif leftline:
                    extra = numpy.array([self.seg_img[row + 1, 1], self.seg_img[row + 1, - 1]])
                elif rightline:
                    extra = numpy.array([self.seg_img[row + 1, self.nCol - 2], self.seg_img[row - 1, self.nCol - 2]])
                else:
                    extra = numpy.array([self.seg_img[row + 1, col + 1], self.seg_img[row + 1, col - 1], self.seg_img[row - 1, col + 1], self.seg_img[row - 1, col -1]])
                neighborhood = numpy.append(extra, neighborhood)

        return neighborhood


