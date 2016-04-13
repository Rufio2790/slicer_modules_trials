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

    #IMPORTANZA DEI NOMI: metti nomi che abbiano un significato, che siano capibili solo leggendoli! Rende molto piu'
    #facile, soprattutto quando il codice e' grande, capire cosa si sta facendo
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

    # Position Widget
    self.positionWidget = ctk.ctkSliderWidget()
    self.positionWidget.singleStep = 0.05 #da calcolare
    self.positionWidget.minimum = -1 #da calcolare
    self.positionWidget.maximum = 1 #da calcolare
    self.positionWidget.value = 0.0
    self.positionWidget.setToolTip("scroll slices through the new axis")
    parametersFormLayout.addRow("Position", self.positionWidget)
    self.positionWidget.connect('valueChanged(double)', self.onPositionWidget)

    # Compute Distance Button
    self.distanceButton = qt.QPushButton("Compute Distance")
    self.distanceButton.toolTip = "Compute Distance between F1 and F2"
    self.distanceButton.enabled = False
    parametersFormLayout.addRow(self.distanceButton)
    #personalmente preferisco mettere la connection subito sotto l'elemento della GUI che e' collegato alla funzione
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
    # IN QUESTO CASO self.logic NON MI SEMBRA CORRETTO. O MEGLIO, NON VEDO UTILE SALVARE LOGIC COME MEMBRO DELLA CLASSE GUI

    # self.logic = surgeonEyeViewLogic()
    # self.logic.CalcDistance(self.Fiducials.currentNode())
    # self.logic.run() #L'errore e' qui in quanto non so come passare F calcolati con CalcDistance alla funzione run
    # self.distanceValue.setText('%.3f' % self.logic.distance)

    logic = surgeonEyeViewLogic()
    self.distanceValue, self.F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText(self.distanceValue)
    logic.run(self.F)


  def ondistanceButton(self):

    logic = surgeonEyeViewLogic()
    self.distanceValue, F = logic.calcDistance(self.FiducialsSelector.currentNode())
    self.distanceValueLabel.setText(self.distanceValue)

  def onPositionWidget(self):
    logic = surgeonEyeViewLogic()
    # logic.moveSliceToNewReferenceFrame()


#
# surgeonEyeViewLogic
#

class surgeonEyeViewLogic(ScriptedLoadableModuleLogic):


  #definisci le funzioni con una convenzione di nomi: una funzione la facciamo iniziare sempre con una lettera minuscola
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

    #perché salvarla come membro dell'oggetto logic?? meglio fare un return della distance e al massimo salvarla come
    # membro della classe GUI
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
    P2 = numpy.random.uniform(-60,70 , 3)

    # Calcolo i 3 assi con Origine del sistema in P0
    Az = P1 - P0 #asse 1
    Az = Az / numpy.linalg.norm(Az)

    Ax = numpy.cross(Az, P2 - P0) #asse 3
    Ax = Ax / numpy.linalg.norm(Ax)

    Ay = numpy.cross(Az, Ax) #asse 2

    #Creo la Matrice
    MRT = numpy.zeros((4, 4))

    # MRT[0, :3] = Az
    # MRT[1, :3] = Ay
    # MRT[2, :3] = Ax
    # MRT[3, :3] = 0
    # MRT[:3, 3] = P0
    # MRT[3, 3] = 1

    #Questo e' il corretto!!!
    MRT[:3, 0] = Ax
    MRT[:3, 1] = Ay
    MRT[:3, 2] = Az
    MRT[:3, 3] = P0
    MRT[3, :3] = 0
    MRT[3, 3] = 1


    print "Asse Z:\n", Az
    print "Asse Y:\n", Ay
    print "Asse X:\n", Ax
    print "P0:\n", P0
    print "P1:\n", P1
    print "P2:\n", P2
    print "Transformation Matrix:\n", MRT

    det = numpy.linalg.det(MRT)
    print "Determinant: ", det

    self.createNewLinearTransform(MRT)
    return MRT
    #TODO Check why is not working correctly
    # if det == 1:
    #   logging.info('Processing completed')
    #   return True
    # else:
    #   logging.info("Error computing the rotation")
    #   slicer.util.delayDisplay("Error computing the rotation")
    #   return False

  def createNewLinearTransform(self, numpyMatrix):

    linearTransformNode = slicer.vtkMRMLLinearTransformNode()
    transformMatrix = vtk.vtkMatrix4x4()

    InvMat = numpy.linalg.inv(numpyMatrix)
    print "Inverse Transformation Matrix:\n", InvMat
    det2 = numpy.linalg.det(InvMat)
    print "DeterminantInvMat: ", det2

    # F0 = numpy.ones(4)
    # F1 = numpy.ones(4)
    # F0[:3] = F[0]
    # F0 = F0.reshape(4, 1)
    # F1[:3] = F[1]
    # F1 = F1.reshape(4, 1)
    # F0 = InvMat.dot(F0)
    # F1 = InvMat.dot(F1)
    # print "F0:\n", F0
    # print "F1:\n", F1
    InvMat = InvMat.ravel()
    InvMat = InvMat.squeeze()
    print InvMat
    transformMatrix.DeepCopy([InvMat[0], InvMat[1], InvMat[2], InvMat[3],
                              InvMat[4], InvMat[5], InvMat[6], InvMat[7],
                              InvMat[8], InvMat[9], InvMat[10], InvMat[11],
                              InvMat[12], InvMat[13], InvMat[14], InvMat[15]])


    slicer.mrmlScene.AddNode(linearTransformNode)
    linearTransformNode.SetAndObserveMatrixTransformToParent(transformMatrix)
    linearTransformNode.SetName('TransformToLCS')
    self.moveSliceToNewReferenceFrame(transformMatrix)

  def moveSliceToNewReferenceFrame(self, transformMatrix):


    #Recupero il nodo
    sliceNodes = slicer.mrmlScene.GetNodesByClass('vtkMRMLSliceNode')
    axialSliceNode = sliceNodes.GetItemAsObject(0)
    scene = slicer.mrmlScene
    inode = scene.GetNodeByID('vtkMRMLScalarVolumeNode1')
    event = vtk.vtkIntArray()
    event.InsertNextValue(slicer.vtkMRMLTransformableNode.TransformModifiedEvent)
    inode.SetAndObserveNodeReferenceID('transform', 'vtkMRMLLinearTransformNode4', event)
    #axialSliceNode.SetSliceOrigin(transformMatrix.GetElement[0, 3], transformMatrix.GetElement[1, 3], transformMatrix.GetElement[2, 3])



    #per il widget bisogna layoutManager
    # m = slicer.app.layoutManager()
    # rw = lm.sliceWidget('Red')
    # sl = rw.sliceLogic()


    #Probabilmente la vtkCamera non ti serve, o meglio, non ti serve se lavoriamo sulle slice.
    #Cerca di guardare la documentazione sul reformat, recupera il nodo della slice rossa dalla scena e modifica
    #posizione e orientamento della slice in modo che faccia quello che vogliamo noi.

    #Fallo in un volume, la MRHead dei sample data di slicer va benissimo, altrimenti non so se ti permette di scorrere
    # la slice. Ricorda che l'elettrodo si suppone debba stare al centro della slice, quindi é possibile che tu debba
    #applicare un offset o magari c'é giá una funzione, controlla!!


    #Infine, cerca di lavorare in modo da creare funzioni differenti, la cosa migliore é che una funzione faccia una sola cosa,
    # prenda in input una serie di argomenti e ritorni il risultato. In questo modo risparmi anche codice, nel senso che
    # l'idea di scrivere funzioni, oltre a rendere piu' semplice il debug, fa sí anche che quello che scrivi in una funzione
    # sia possibile riutilizzarlo in altre parti del codice

    #Abbozzo l'ultimo punto (Sicuramente sbagliato....)



    # camera = vtk.vtkCamera()
    # cameraPositionMultiplier = 5
    #
    #
    # m = axialSliceNode.GetSliceToRAS()
    # rSliceToRAS = numpy.matrix([[m.GetElement(0, 0), m.GetElement(0, 1), m.GetElement(0, 2)],
    #                       [m.GetElement(1, 0), m.GetElement(1, 1), m.GetElement(1, 2)],
    #                       [m.GetElement(2, 0), m.GetElement(2, 1), m.GetElement(2, 2)]])
    #
    # det = numpy.linalg.det(rSliceToRAS)
    # if det > 0:  # right hand
    #   y = numpy.array([0, 0, -cameraPositionMultiplier])
    # elif det < 0:  # left hand
    #   y = numpy.array([0, 0, cameraPositionMultiplier])
    #
    # x = numpy.matrix([[m.GetElement(0, 0), m.GetElement(0, 1), m.GetElement(0, 2)],
    #                [m.GetElement(1, 0), m.GetElement(1, 1), m.GetElement(1, 2)],
    #                [m.GetElement(2, 0), m.GetElement(2, 1), m.GetElement(2, 2)]])
    #
    # # Calculating position
    # position = numpy.inner(x, y)
    # camera.SetPosition(-position[0, 0], -position[0, 1], -position[0, 2])
    #
    # print axialSliceNode
    # print m
    # print rSliceToRAS
    #
    # return True

# Widgets for setting translation/rotation part of a linear transform
    # Position widget
    #
    # self.positionSliderWidget = slicer.qMRMLTransformSliders()
    # self.positionSliderWidget.Title = 'Position'
    # self.positionSliderWidget.TypeOfTransform = slicer.qMRMLTransformSliders.TRANSLATION
    # self.positionSliderWidget.CoordinateReference = slicer.qMRMLTransformSliders.LOCAL
    # self.positionSliderWidget.setMRMLScene(slicer.mrmlScene)
    # self.positionSliderWidget.setMRMLTransformNode(self.getPivotToRasTransformNode())
    # planesFormLayout.addRow("Translation", self.positionSliderWidget)
    #
    # # Orientation widget
    #
    # self.orientationSliderWidget = slicer.qMRMLTransformSliders()
    # self.orientationSliderWidget.Title = 'Orientation'
    # self.orientationSliderWidget.setMRMLScene(slicer.mrmlScene)
    #
    # # Setting of qMRMLTransformSliders.TypeOfTransform is not robust: it has to be set after setMRMLScene and
    # # has to be set twice (with setting the type to something else in between).
    # # Therefore the following 3 lines are needed, and they are needed here:
    #
    # self.orientationSliderWidget.TypeOfTransform = slicer.qMRMLTransformSliders.ROTATION
    # self.orientationSliderWidget.TypeOfTransform = slicer.qMRMLTransformSliders.TRANSLATION
    # self.orientationSliderWidget.TypeOfTransform = slicer.qMRMLTransformSliders.ROTATION
    # self.orientationSliderWidget.CoordinateReference = slicer.qMRMLTransformSliders.LOCAL
    # self.orientationSliderWidget.minMaxVisible = False
    # self.orientationSliderWidget.setMRMLTransformNode(self.getPivotToRasTransformNode())
    # planesFormLayout.addRow("Orientation", self.orientationSliderWidget)