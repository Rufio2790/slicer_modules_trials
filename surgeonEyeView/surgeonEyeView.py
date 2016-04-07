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
    self.parent.contributors = ["Giuseppe De Luca (POLIMI)"] # replace with "Firstname Lastname (Organization)"
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
    # self.distanceButton.connect('clicked(bool)', self.ondistanceButton)
    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.Fiducials.currentNode()
    # self.distanceButton.enabled = self.Fiducials.currentNode()

  def onApplyButton(self):
    self.logic = surgeonEyeViewLogic()
    self.logic.run(self.Fiducials.currentNode())
    self.distanceValue.setText('%.3f'%self.logic.distance)

  # def ondistanceButton(self):
  #   logic = surgeonEyeViewLogic()
  #   logic.distance(self.Fiducials.currentNode())
#
# surgeonEyeViewLogic
#

class surgeonEyeViewLogic(ScriptedLoadableModuleLogic):


  def run(self, Fiducials):
    """
    Run the actual algorithm
    """

    logging.info('Processing started')

    F = numpy.zeros((2, 3))
    FidCoords = numpy.empty(3)
    fiducialCount = Fiducials.GetNumberOfFiducials()
    for x in range(fiducialCount):
      Fiducials.GetNthFiducialPosition(x, FidCoords)
      F[x] = FidCoords
    x = numpy.power(F[0, 0] - F[1, 0], 2)
    y = numpy.power(F[0, 1] - F[1, 1], 2)
    z = numpy.power(F[0, 2] - F[1, 2], 2)
    self.distance = numpy.sqrt((x + y + z))
    #calcolo equazione parametrica retta passante per i due punti
    t = numpy.linspace(-10,10,50) #vettore che rappresenta il dominio della retta
    x = F[0,0] + (F[1,0]-F[0,0])*t
    y = F[0,1] + (F[1,1]-F[0,1])*t
    z = F[0,2] + (F[1,2]-F[0,2])*t
    #calcolo coordinate punto medio per usarlo come origine del sistema di riferimento
    xm = (F[0,0] + F[1,0]) / 2
    ym = (F[0,1] + F[1,1]) / 2
    zm = (F[0,2] + F[1,2]) / 2
    PM = [xm,ym,zm]
    PM = numpy.asarray(PM)
    P0 = F[0]
    P1 = F[1]
    # Calcolo i 3 assi con Origine del sistema in P0
    Az = PM - P0 #asse 1
    Az = Az / numpy.linalg.norm(Az)
    Ax = numpy.cross(Az, P1 - P0) #asse 3
    Ax = Ax / numpy.linalg.norm(Ax)
    Ay = numpy.cross(Ax, Az) #asse 2
    #Creo la Matrice
    MRT = numpy.zeros((4, 4))
    # MRT[0, :3] = A1
    # MRT[1, :3] = A2
    # MRT[2, :3] = A3
    # MRT[3, :3] = 0
    # MRT[:3, 3] = P0
    # MRT[3, 3] = 1
    MRT[:3, 0] = Az
    MRT[:3, 1] = Ay
    MRT[:3, 2] = Ax
    MRT[:3, 3] = P0
    MRT[3, :3] = 0
    MRT[3, 3] = 1

    #Calcolo l'inversa di MRT
    IMRT = numpy.transpose(MRT [0:3, 0:3])
    P00 = P0.reshape(3,1)
    NewCenter = numpy.dot(IMRT, P00)
    #test punto

    print Az
    print Ay
    print Ax
    print P0
    print P1
    print PM
    print MRT
    print NewCenter

    logging.info('Processing completed')

    return True