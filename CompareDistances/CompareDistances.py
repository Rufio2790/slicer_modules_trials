import os
import unittest
import numpy
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *

#
# CompareDistances
#

class CompareDistances(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CompareDistances" # TODO make this more human readable by adding spaces
    self.parent.categories = ["TestModules"]
    self.parent.contributors = ["Davide Scorza (Vicomtech)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    Module that compare two distance list and gives back some precision indexes like RMSE, standard deviation, variance, medium value of each dataset.
    Useful to compare the precision of a fiducial list obtained from different externals tracker.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# CompareDistancesWidget
#

class CompareDistancesWidget(ScriptedLoadableModuleWidget):
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
    # input List 1
    #
    self.inputList = slicer.qMRMLNodeComboBox()
    self.inputList.nodeTypes = ( ("vtkMRMLAnnotationHierarchyNode"), "" )
    self.inputList.selectNodeUponCreation = True
    self.inputList.addEnabled = True
    self.inputList.removeEnabled = False
    self.inputList.noneEnabled = False
    self.inputList.showHidden = False
    self.inputList.showChildNodeTypes = False
    self.inputList.setMRMLScene( slicer.mrmlScene )
    self.inputList.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input List 1: ", self.inputList)

    #
    # input List 2
    #
    self.inputList2 = slicer.qMRMLNodeComboBox()
    self.inputList2.nodeTypes = ( ("vtkMRMLAnnotationHierarchyNode"), "" )
    self.inputList2.selectNodeUponCreation = False
    self.inputList2.addEnabled = True
    self.inputList2.removeEnabled = True
    self.inputList2.noneEnabled = False
    self.inputList2.showHidden = False
    self.inputList2.showChildNodeTypes = False
    self.inputList2.setMRMLScene( slicer.mrmlScene )
    self.inputList2.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Input List 2: ", self.inputList2)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Compare")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)
    
    
     # Results
    self.StatisticalValuesLabel = qt.QLabel()
    self.StatisticalValuesLabel.setText( "Metrics related to the precision of the datasets: " )
    parametersFormLayout.addRow(self.StatisticalValuesLabel)
    
    # Add vertical spacer
    self.layout.addStretch(3)
    
    self.MeanValue = qt.QLabel()
    parametersFormLayout.addRow("Error Mean Value: ", self.MeanValue)  
    self.RMSErrorValue = qt.QLabel()
    parametersFormLayout.addRow("Root Mean Square Error: ", self.RMSErrorValue)
    self.VarianceValue = qt.QLabel()
    parametersFormLayout.addRow("Error variance: ", self.VarianceValue)
    self.stdValue = qt.QLabel()
    parametersFormLayout.addRow("Error standard deviation: ", self.stdValue)    
    
    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputList.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputList2.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)
    

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputList.currentNode() and self.inputList2.currentNode()

  def onApplyButton(self):
    #slicer.app.processEvents()
    self.logic = CompareDistancesLogic()
    print("Run the algorithm")
    self.logic.run(self.inputList.currentNode(), self.inputList2.currentNode())
    self.MeanValue.setText('%.4f'%self.logic.info['MeanError'])
    self.RMSErrorValue.setText('%.4f'%self.logic.info['RMSE'])
    self.VarianceValue.setText('%.4f'%self.logic.info['Variance'])
    self.stdValue.setText('%.4f'%self.logic.info['std'])

#
# CompareDistancesLogic
#
class CompareDistancesLogic(ScriptedLoadableModuleLogic):


  def run(self,inputList,inputList2):
    """
    Run the actual algorithm
    """
    self.info={}
    
    #get the number of data 
    rulerNumber = inputList.GetNumberOfChildrenNodes() 
    
    #initialize the array of the dimensions
    r1 = numpy.empty(rulerNumber)
    r2 = numpy.empty(rulerNumber)

    #get the distances from the ruler lists and put into numpy array
    for x in xrange(rulerNumber):
        #node is supposed to contain the gold standard distances (Points provided by CMM)
        node = inputList.GetNthChildNode(x)
        ruler = node.GetAssociatedNode()
        r1[x] = ruler.GetDistanceMeasurement()
        
        #node2 is supposed to contain the distances adquired
        node2 = inputList2.GetNthChildNode(x)
        ruler2 = node2.GetAssociatedNode()
        r2[x] = ruler2.GetDistanceMeasurement()
#        r2[x] = inputList2.GetNthChildNode(x).GetAssociatedNode().GetDistanceMeasurement()
   
    savePath1 = "C:\Users\dscorza\Documents\ORXXI\Process\Slicer4\Distance\ " + inputList2.GetName()
    n = open(savePath1 + ".txt", "w")
    for fidIndex in xrange(rulerNumber):
        n.write(str(r2[fidIndex]) + "\n")
        
    n.close()
    
    savePathCommon = "C:\Users\dscorza\Documents\ORXXI\Process\Slicer4\Distance\ALLdistances.txt"
    l = open(savePathCommon, "a")
    
    for fidIndex in xrange(rulerNumber):
        if fidIndex == 0:
            l.write(inputList2.GetName() + "\t" + str(r2[fidIndex]))
        elif fidIndex == rulerNumber - 1:
            l.write("\t" + str(r2[fidIndex]) + "\n")
        else:
            l.write("\t" + str(r2[fidIndex]) )
        
    l.close()
        
        
    error = numpy.abs(r1 - r2)                                   #array of the error between the two signals
    
    meanError = numpy.sum(error)/rulerNumber                     #mean value of the error
    RMSError = numpy.sqrt(numpy.sum((error)**2)/rulerNumber)     #Root Mean Square Error
    variance = numpy.var(error)                                  # Variance of the error
    std = numpy.std(error)                                       #standard deviation of the error
    
    self.info['MeanError'] = meanError
    self.info['RMSE'] = RMSError
    self.info['Variance'] = variance
    self.info['std'] = std
    
    #save these metrics in a txt file
    savePath = "C:\Users\dscorza\Documents\ORXXI\Process\Slicer4\DataSet_Precision\ " + inputList2.GetName()
    
    f = open(savePath + ".txt", "w")
    f.write("Mean Error \t RMSE \t Variance \t std \n")
    f.write(str(meanError) + "\t" + str(RMSError) + "\t" + str(variance) + "\t" + str(std) + "\n")
    #f.write("Mean Error: " + str(meanError) + "\n")
    #f.write("RMSE: " + str(RMSError) + "\n")
    #f.write("Variance: " + str(variance) + "\n")
    #f.write("std: " + str(std) + "\n")
    
    f.close()
    
    savePathCommon2 = "C:\Users\dscorza\Documents\ORXXI\Process\Slicer4\DataSet_Precision\ALLStat.txt"
    h = open(savePathCommon2, "a")
    
   # h.write("Dataset \t Mean Error \t RMSE \t Variance \t std \n")
    h.write(inputList2.GetName() + "\t" + str(meanError) + "\t" + str(RMSError) + "\t" + str(variance) + "\t" + str(std) + "\n")
        
    h.close()
    
    return True

