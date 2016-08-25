import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import math
import numpy as np

#
# AnglePoints
#

class AnglePoints(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "angle4points" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Lisa Plaino"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    Calculates the angle between two VECTORS rapresented as:
    a) 2 pairs of fiducial markups (starting and ending points)
    b) 3 points (the second point is the common origin)
    c) 2 points (as components of vector starting from the origin)
    """
    self.parent.acknowledgementText = """
    Organization, grant and thanks
"""
    self.parent=parent
#
# angleWidget
#

class AnglePointsWidget(ScriptedLoadableModuleWidget):

  def setup(self):

    ScriptedLoadableModuleWidget.setup(self)

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )

    # Collapsible button
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)


    # Markup selector
    self.markupSelectorLabel = qt.QLabel()
    self.markupSelectorLabel.setText( "Markup list: " )
    self.markupSelector = slicer.qMRMLNodeComboBox()
    self.markupSelector.nodeTypes = ( "vtkMRMLMarkupsFiducialNode", "" )
    self.markupSelector.noneEnabled = False
    self.markupSelector.selectNodeUponCreation = True
    self.markupSelector.setMRMLScene( slicer.mrmlScene )
    self.markupSelector.setToolTip( "Pick a list that contains 2, 3 or 4 markups" )
    parametersFormLayout.addRow(self.markupSelectorLabel, self.markupSelector)


    # Angle label
    self.angleValueLabel = qt.QLabel()
    self.angleValueLabel.setText( " " )


    # Apply button
    self.computeButton = qt.QPushButton( "Compute" )
    self.computeButton.toolTip = "Compute information for the selected markup"
    parametersFormLayout.addWidget( self.computeButton )
    self.UpdatecomputeButtonState()
    parametersFormLayout.addRow("Angle: ", self.angleValueLabel)


    # connections
    self.computeButton.connect('clicked()', self.onCompute)
    self.markupSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.onMarkupSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

  def UpdatecomputeButtonState(self):
    if not self.markupSelector.currentNode() :
      self.computeButton.enabled = False
    else:
      self.computeButton.enabled = True

  def onMarkupSelect(self, node):
    self.UpdatecomputeButtonState()

  def onCompute(self):
    #slicer.app.processEvents()
    self.logic = AnglePointsLogic()
    self.logic.computeAngle(self.markupSelector.currentNode())
    self.angleValueLabel.setText(str(self.logic.alphadeg)+" deg")

#
# angleLogic
#

class AnglePointsLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def __init__(self):
    self.alphadeg = 0


  def computeAngle(self, markupNode):
    #calls the function that suits the number of points inserted
    n = markupNode.GetNumberOfFiducials()
    abcd = [] #[A, B, C, D]
    mypoints=[] #markups ras list

    for i in range(n):
      ras = [0, 0, 0]
      markupNode.GetNthFiducialPosition(i, ras)
      mypoints.append(ras)

    abcd = findABCD(n, mypoints)
    A = abcd[0]
    B = abcd[1]
    C = abcd[2]
    D = abcd[3]

    if n > 1 and n < 5:
      v1 = vect_from_2points(A, B)
      v2 = vect_from_2points(C, D)
      alpharad = angle_between_2vect(v1, v2)
      self.alphadeg = math.degrees(alpharad)
      printingResults(n, abcd, v1, v2, self.alphadeg)
      return self.alphadeg
    else:
      print'Number of points inserted incorrect'
      slicer.util.delayDisplay('Incorrect input')


class Vect(object):
  def __init__(self, a, b, c):
    self.a = a
    self.b = b
    self.c = c
  def module(self):
    m = math.sqrt(self.a ** 2 + self.b ** 2 + self.c ** 2)
    return m

def vect_from_2points(p1, p2):
  # p1, p2 are 2 Point objects
  # the function computes the vect from p1 to p2
  return Vect(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])

def scalar(v1, v2):
  # vector * vector --> float
  return (v1.a * v2.a) + (v1.b * v2.b) + (v1.c * v2.c)

def angle_between_2vect(v1, v2):
  m1 = v1.module()
  m2 = v2.module()
  s = (v1.a * v2.a) + (v1.b * v2.b) + (v1.c * v2.c)
  # scalar = m1 * m2 * cos( alpha )
  # cos(alpha) = scalar / ( m1 * m2 )
  alpha = math.acos(s / (m1 * m2))
  return alpha


def findABCD(n, mypoints):
  # n = number of markups in the list
  # assignment of the coordinates of points A, B, C, D
  # n = 4 --> A = mypoint[0], B = mypoint[1], C = mypoint[2], D = mypoint[3]
  # n = 3 --> A = mypoint[0], B = mypoint[1], C = B, D = mypoint[2]
  # n = 2 --> A = 0, B = mypoint[0], C = 0, D = mypoint[1]

  def Apoint(n):
    A = [0, 0, 0]
    if n > 2:
      A = np.array(mypoints[0])
    return A

  def Bpoint(n):
    B = [0, 0, 0]
    if n == 3 or n == 4:
      B = np.array(mypoints[1])
    elif n == 2:
      B = np.array(mypoints[0])
    return B

  def Cpoint(n):
    C = [0, 0, 0]
    if n == 4:
      C = np.array(mypoints[2])
    elif n == 3:
      C = np.array(mypoints[1])
    return C

  def Dpoint(n):
    D = [0, 0, 0]
    if n == 4:
      D = np.array(mypoints[3])
    elif n == 3:
      D = np.array(mypoints[2])
    elif n == 2:
      D = np.array(mypoints[1])
    return D

  return [Apoint(n), Bpoint(n), Cpoint(n), Dpoint(n)]

def printingResults(n, ABCD, v1, v2, alpha):
  # prints info:
  # n = number of markups inserted
  # ABCD = list of points A, B, C, D
  # v1 = vector from A to B
  # v2 = vector from C to D
  # alpha = computed angle between v1 and v2
  print n, "points inserted"
  names = ['A) ', 'B) ', 'C) ', 'D) ']
  for i in range (4):
    print names[i], ABCD[i]
    print "\n"
  print "v1: ", v1.a, v1.b, v1.c
  print "v2: ", v2.a, v2.b, v2.c
  print "The angle between v1 and v2 is: ", alpha , "degrees"
  return()
