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
    This module will be used to verify the correctness of the trajectories implemented. .
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

    #
    # fiducials
    #
    self.Fiducials = slicer.qMRMLNodeComboBox()
    self.Fiducials.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.Fiducials.selectNodeUponCreation = True
    self.Fiducials.addEnabled = False
    self.Fiducials.removeEnabled = False
    self.Fiducials.noneEnabled = False
    self.Fiducials.showHidden = False
    self.Fiducials.showChildNodeTypes = False
    self.Fiducials.setMRMLScene( slicer.mrmlScene )
    self.Fiducials.setToolTip( "Pick the fiducials." )
    parametersFormLayout.addRow("Fiducials: ", self.Fiducials)

    #
    # Compute Distance Button
    #
    # self.distanceButton = qt.QPushButton("Compute Distance")
    # self.distanceButton.toolTip = "Compute Distance between F1 and F2"
    # self.distanceButton.enabled = True
    # parametersFormLayout.addRow(self.distanceButton)


    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # Results
    self.distanceValue = qt.QLabel()
    parametersFormLayout.addRow("Distance between F1 and F2 ", self.distanceValue)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    #self.distanceButton.connect('clicked(bool)', self.ondistanceButton)
    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.Fiducials.currentNode()
    #self.distanceButton.enabled = self.Fiducials.currentNode()

  def onApplyButton(self):
    self.logic = surgeonEyeViewLogic()
    self.logic.CalcDistance(self.Fiducials.currentNode())
    self.logic.run() #L'errore e' qui in quanto non so come passare F calcolati con CalcDistance alla funzione run
    self.distanceValue.setText('%.3f' % self.logic.distance)

  #def ondistanceButton(self):
  #logic2 = surgeonEyeViewLogic()
  #logic2.CalcDistance(self.Fiducials.currentNode())



#
# surgeonEyeViewLogic
#

class surgeonEyeViewLogic(ScriptedLoadableModuleLogic):





  def CalcDistance(self, Fiducials):

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
    self.distance = numpy.sqrt((x + y + z))
    print self.distance  #qua il codice arriva!!!

    return F


  def run(self, F ):
    """
    Run the actual algorithm
    """

    logging.info('Processing started')


    P0 = F[0]
    P1 = F[1]
    P2 = numpy.random.uniform(-60,70 , 3)
    # Calcolo i 3 assi con Origine del sistema in P0
    Az = P1 - P0 #asse 1
    Az = Az / numpy.linalg.norm(Az)
    Ax = numpy.cross(Az, P2 - P0) #asse 3
    Ax = Ax / numpy.linalg.norm(Ax)
    Ay = numpy.cross(Ax, Az) #asse 2
    #Creo la Matrice
    MRT = numpy.zeros((4, 4))
    MRT[0, :3] = Az
    MRT[1, :3] = Ay
    MRT[2, :3] = Ax
    MRT[3, :3] = 0
    MRT[:3, 3] = P0
    MRT[3, 3] = 1
    det = numpy.linalg.det(MRT) #E' uguale a 1!!!
    print det
    # MRT[:3, 0] = Az
    # MRT[:3, 1] = Ay
    # MRT[:3, 2] = Ax
    # MRT[:3, 3] = P0
    # MRT[3, :3] = 0
    # MRT[3, 3] = 1


    print Az
    print Ay
    print Ax
    print P0
    print P1
    print P2
    print MRT


    #Abbozzo l'ultimo punto (Sicuramente sbagliato....)



    camera = vtk.vtkCamera()
    cameraPositionMultiplier = 5
    sliceNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLSliceNode')
    axialSliceNode = sliceNodes.GetItemAsObject(0)
    m = axialSliceNode.GetSliceToRAS()
    rSliceToRAS = numpy.matrix([[m.GetElement(0, 0), m.GetElement(0, 1), m.GetElement(0, 2)],
                          [m.GetElement(1, 0), m.GetElement(1, 1), m.GetElement(1, 2)],
                          [m.GetElement(2, 0), m.GetElement(2, 1), m.GetElement(2, 2)]])

    det = numpy.linalg.det(rSliceToRAS)
    if det > 0:  # right hand
      y = numpy.array([0, 0, -cameraPositionMultiplier])
    elif det < 0:  # left hand
      y = numpy.array([0, 0, cameraPositionMultiplier])

    x = numpy.matrix([[m.GetElement(0, 0), m.GetElement(0, 1), m.GetElement(0, 2)],
                   [m.GetElement(1, 0), m.GetElement(1, 1), m.GetElement(1, 2)],
                   [m.GetElement(2, 0), m.GetElement(2, 1), m.GetElement(2, 2)]])

    # Calculating position
    position = numpy.inner(x, y)
    camera.SetPosition(-position[0, 0], -position[0, 1], -position[0, 2])

    print axialSliceNode
    print m
    print rSliceToRAS



    logging.info('Processing completed')

    return True