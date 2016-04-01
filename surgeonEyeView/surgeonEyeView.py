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
    # Inserire Tasto Compute Distance che restituisce a schermo la distanza
    # Questo dovrà essere collegato alla funzione ComputeDistance
    #......

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.Fiducials.currentNode()

  def onApplyButton(self):
    logic = surgeonEyeViewLogic()
    logic.run(self.Fiducials.currentNode())

#
# surgeonEyeViewLogic
#

class surgeonEyeViewLogic(ScriptedLoadableModuleLogic):


  def ComputeDistance(self, Fiducials):

    F = numpy.zeros((2, 3))
    FidCoords = numpy.empty(3)
    fiducialCount = Fiducials.GetNumberOfFiducials()
    for x in range(fiducialCount):
      Fiducials.GetNthFiducialPosition(x, FidCoords)
      F[x] = FidCoords

    x = numpy.power(F[0,0] - F[1,0], 2)
    y = numpy.power(F[0,1] - F[1,1], 2)
    z = numpy.power(F[0,2] - F[1,2], 2)
    distance=numpy.sqrt((x+y+z))
    #inserire distance nel label della GUI

  def ReferenceSystem(self, F): #Passare a questa funzione F che sono le coordinate già pronte calcolate prima

    # calcolo equazione parametrica retta passante per i due punti (non serve, da cancellare)
    t = numpy.linspace(-10, 10, 50)  # vettore che rappresenta il dominio della retta
    x = F[0, 0] + (F[1, 0] - F[0, 0]) * t
    y = F[0, 1] + (F[1, 1] - F[0, 1]) * t
    z = F[0, 2] + (F[1, 2] - F[0, 2]) * t
    # calcolo coordinate punto medio per usarlo come origine del sistema di riferimento
    xm = (F[0, 0] + F[1, 0]) / 2
    ym = (F[0, 1] + F[1, 1]) / 2
    zm = (F[0, 2] + F[1, 2]) / 2
    PM = [xm, ym, zm]
    PM = numpy.asarray(PM)
    P0 = F[0]
    P1 = F[1]
    # Calcolo i 3 assi con Origine del sistema in P0
    A1 = PM - P0  # asse 1
    A1 = A1 / numpy.linalg.norm(A1)
    A3 = numpy.cross(A1, P1 - P0)  # asse 3
    A3 = A3 / numpy.linalg.norm(A3)
    A2 = numpy.cross(A3, A1)  # asse 2
    # Creo la Matrice
    MT = numpy.zeros((3, 4))
    for i in range(0, 3):
      MT[0, i] = A1[i]
      MT[1, i] = A2[i]
      MT[2, i] = A3[i]
      MT[i, 3] = P0[i]

    # I 3 punti nel nuovo sistema di riferimento
    P = P0
    AT = [numpy.dot(P - P0, A1), numpy.dot(P - P0, A2), numpy.dot(P - P0, A3)]  # prodotto scalare tra P-P0 e A1
    P = PM
    BT = [numpy.dot(P - P0, A1), numpy.dot(P - P0, A2), numpy.dot(P - P0, A3)]
    P = P1
    CT = [numpy.dot(P - P0, A1), numpy.dot(P - P0, A2), numpy.dot(P - P0, A3)]
    print A1
    print A2
    print A3
    print P0
    print P1
    print PM
    print AT
    print BT
    print CT
    print MT

  def run(self): #Passare a questa funzione F (che sono le coordinate già pronte calcolate prima) e non Fiducials
    """
    Run the actual algorithm
    """

    logging.info('Processing started')

    #In questa parte avevo pensato di inserire il terzo punto mancante del modulo

    logging.info('Processing completed')

    return True
