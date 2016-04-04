import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np

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
    parametersFormLayout = qt.QVBoxLayout(parametersCollapsibleButton)

    sphereWidget = qt.QGroupBox()
    sphereLayout = qt.QFormLayout(sphereWidget)
    parametersFormLayout.addWidget(sphereWidget)

    self.sphereCenterValue = np.zeros((3),'float32')
    self.sphereCenter = qt.QTextEdit()
    self.sphereCenter.setFixedHeight(25)
    self.sphereCenter.setText(str(self.sphereCenterValue))
    self.sphereCenter.setToolTip("Insert Center Coordinates")
    sphereLayout.addRow('Sphere Center: ', self.sphereCenter)
    self.sphereCenter.connect('textChanged()', self.onSphereCenterChanged)

    self.sphereRadiusValue = np.array((0),'float32')
    self.sphereRadius = qt.QTextEdit()
    self.sphereRadius.setFixedHeight(25)
    self.sphereRadius.setText(str(self.sphereRadiusValue))
    self.sphereRadius.setToolTip("Insert Radius Dimensions")
    sphereLayout.addRow('Sphere Radius: ', self.sphereRadius)
    self.sphereRadius.connect('textChanged()', self.onSphereRadiusChanged)

    # Apply Button
    self.applySphereButton = qt.QPushButton("Create Sphere")
    self.applySphereButton.toolTip = "Generate Model"
    self.applySphereButton.enabled = True
    sphereLayout.addRow(self.applySphereButton)
    self.applySphereButton.connect('clicked(bool)', self.onApplySphereButtonClicked)

    cubeWidget = qt.QGroupBox()
    cubeLayout = qt.QFormLayout(cubeWidget)
    parametersFormLayout.addWidget(cubeWidget)

    self.cubeCenterValue = np.zeros((3), 'float32')
    self.cubeCenter = qt.QTextEdit()
    self.cubeCenter.setFixedHeight(25)
    self.cubeCenter.setText(str(self.cubeCenterValue))
    self.cubeCenter.setToolTip("Insert Center Coordinates")
    cubeLayout.addRow('Cube Center: ', self.cubeCenter)
    self.cubeCenter.connect('textChanged()', self.onCubeCenterChanged)

    self.cubeDimensionsValue = np.zeros((3), 'float32')
    self.cubeDimensions = qt.QTextEdit()
    self.cubeDimensions.setFixedHeight(25)
    self.cubeDimensions.setText(str(self.cubeDimensionsValue))
    self.cubeDimensions.setToolTip("Insert Center Coordinates")
    cubeLayout.addRow('Cube Dimensions (X,Y,Z): ', self.cubeDimensions)
    self.cubeDimensions.connect('textChanged()', self.onCubeDimensionsChanged)

    # Apply Button
    self.applyCubeButton = qt.QPushButton("Create Cube")
    self.applyCubeButton.toolTip = "Generate Model"
    self.applyCubeButton.enabled = True
    cubeLayout.addRow(self.applyCubeButton)
    self.applyCubeButton.connect('clicked(bool)', self.onApplyCubeButtonClicked)

    # Add vertical spacer
    self.layout.addStretch(1)

  def cleanup(self):
    pass

  def onSphereCenterChanged(self):

    center = self.sphereCenter.toPlainText().split(',')

    self.sphereCenterValue = np.zeros((3),'float32')

    for i in range(len(center)):
      if center[i] == '':
        continue
      self.sphereCenterValue[i] = float(center[i])

    print self.sphereCenterValue
    return True


  def onSphereRadiusChanged(self):

    radius = self.sphereRadius.toPlainText()

    m = radius.split(',')
    self.sphereRadiusValue = float(m[0])


    print self.sphereRadiusValue
    return True

  def onApplySphereButtonClicked(self):

    logic = vtkModelLogic()
    center = self.sphereCenterValue
    print 'center: ', center
    radius = self.sphereRadiusValue
    logic.createSphereModel(center, radius)

  def onCubeCenterChanged(self):

    center = self.cubeCenter.toPlainText().split(',')

    self.cubeCenterValue = np.zeros((3), 'float32')

    for i in range(len(center)):
      if center[i] == '':
        continue
      self.cubeCenterValue[i] = float(center[i])

    print self.cubeCenterValue
    return True

  def onCubeDimensionsChanged(self):

    dim = self.cubeDimensions.toPlainText().split(',')

    self.cubeDimensionsValue = np.zeros((3), 'float32')

    for i in range(len(dim)):
      if dim[i] == '':
        continue
      self.cubeDimensionsValue[i] = float(dim[i])

    print self.cubeDimensionsValue
    return True

  def onApplyCubeButtonClicked(self):

    logic = vtkModelLogic()
    center = self.cubeCenterValue
    print 'center: ', center
    dim = {'X': self.cubeDimensionsValue[0], 'Y': self.cubeDimensionsValue[1], 'Z': self.cubeDimensionsValue[2]}
    logic.create3DRectangle(center, dim)

# vtkModelLogic
class vtkModelLogic(ScriptedLoadableModuleLogic):
  """
  """

  def createSphereModel(self, center, radius):
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
    sphere.SetCenter(center)
    sphere.SetRadius(radius)
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

  def create3DRectangle(self, center, dimensions):

    # generating Nodes for displaying a new model
    modelNode = slicer.vtkMRMLModelNode()
    dispNode = slicer.vtkMRMLModelDisplayNode()
    transform = slicer.vtkMRMLLinearTransformNode()

    # Display node characteristics
    dispNode.SetVisibility(True)
    dispNode.SetSliceIntersectionVisibility(True)
    dispNode.SetOpacity(0.1)
    dispNode.SetColor(1, 0, 0)
    dispNode.SetScene(slicer.mrmlScene)

    # generate sphere data
    cube = vtk.vtkCubeSource()
    cube.SetCenter(center)
    cube.SetXLength(dimensions['X'])
    cube.SetYLength(dimensions['Y'])
    cube.SetZLength(dimensions['Z'])
    cube.Update()

    # adding necessary nodes to the Scene
    slicer.mrmlScene.AddNode(dispNode)
    slicer.mrmlScene.AddNode(transform)
    slicer.mrmlScene.AddNode(modelNode)

    # model node name and associations!
    modelNode.SetName("CubeModelNode")
    modelNode.SetScene(slicer.mrmlScene)
    modelNode.SetAndObserveTransformNodeID(transform.GetID())
    modelNode.SetAndObserveDisplayNodeID(dispNode.GetID())

    apd = vtk.vtkAppendPolyData()
    apd.AddInputData(cube.GetOutput())
    apd.Update()

    # adding model node poly data! Here there are  sphere's data!!!
    modelNode.SetAndObservePolyData(apd.GetOutput())
    # update the scene
    slicer.mrmlScene.Modified()

class vtkModelTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_vtkModel1()

  def test_vtkModel1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
      ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
    )

    for url, name, loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = vtkModelLogic()
    self.assertTrue(logic.hasImageData(volumeNode))
    self.delayDisplay('Test passed!')