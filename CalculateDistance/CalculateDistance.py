from __main__ import qt, slicer

#
# MarkupsInfo module
#

class CalculateDistance:
  def __init__(self, parent):
    import string
    parent.title = "CalculateDistance"
    parent.categories = ["TestModules"]
    parent.contributors = ["Davide Scorza (Vicomtech)"]
    parent.helpText = string.Template("""
Compute the distances between the first point and the next points of a markup list and Put the results in a Ruler List. Be careful to create a new List
in "annotation" and select the new one every time a new markup distance list want to be adquired, otherwise the distances will be put in the same ruler List
    """).substitute({ 'a':parent.slicerWikiUrl, 'b':slicer.app.majorVersion, 'c':slicer.app.minorVersion })
    parent.acknowledgementText = """
    Supported by SparKit and the Slicer Community. See http://www.slicerrt.org for details.
    """
    self.parent = parent

#
# Widget
#

class CalculateDistanceWidget:

  def __init__(self, parent=None):
    self.parent = parent
    self.logic = None
    
  def setup(self):

    frame = qt.QFrame()
    layout = qt.QFormLayout()
    frame.setLayout( layout )
    self.parent.layout().addWidget( frame )
    
    # Markup selector
    self.markupSelectorLabel = qt.QLabel()
    self.markupSelectorLabel.setText( "Markup list: " )
    self.markupSelector = slicer.qMRMLNodeComboBox()
    self.markupSelector.nodeTypes = ( "vtkMRMLMarkupsFiducialNode", "" )
    self.markupSelector.noneEnabled = False
    self.markupSelector.selectNodeUponCreation = True
    self.markupSelector.setMRMLScene( slicer.mrmlScene )
    self.markupSelector.setToolTip( "Pick the markup list to be filled" )
    layout.addRow(self.markupSelectorLabel, self.markupSelector)    
        
    # Apply button
    self.computeButton = qt.QPushButton("Compute")
    self.computeButton.toolTip = "Compute information for the selected markup"    
    layout.addWidget(self.computeButton)
    self.UpdatecomputeButtonState()
     

    # connections
    self.computeButton.connect('clicked()', self.onCompute)
    self.markupSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.onMarkupSelect)    

  def UpdatecomputeButtonState(self):
    if not self.markupSelector.currentNode() :
      self.computeButton.enabled = False
    else:
      self.computeButton.enabled = True      
    
  def onMarkupSelect(self, node):
    self.UpdatecomputeButtonState()

  def onCompute(self):
    slicer.app.processEvents()
    self.logic = CalculateDistanceLogic(self.markupSelector.currentNode())
    self.DistanceList = slicer.vtkMRMLAnnotationRulerNode()
    
#
# Logic
#
    
class CalculateDistanceLogic:
  """Implement the logic to compute markup info
  Nodes are passed in as arguments.
  Results are stored as 'info' instance variable.
  """
  
  def __init__(self, markupNode):    
    
    self.info={}
    
    #put the first fiducial coordinates into startPtCoords
    startPtCoords = [0.0, 0.0, 0.0]
    markupNode.GetNthFiducialPosition(0,startPtCoords) 
    
    endPtCoords = [0.0, 0.0, 0.0]
    
    #create rulers node in a ruler List the calculate the distances between fiducials
    #in this way it's easy to manage the distance data in other future modules
    for fidIndex in xrange(markupNode.GetNumberOfFiducials()): 
        ruler = slicer.mrmlScene.AddNode(slicer.vtkMRMLAnnotationRulerNode())
        markupNode.GetNthFiducialPosition(fidIndex,endPtCoords)
        ruler.SetPosition1(startPtCoords[0],startPtCoords[1],startPtCoords[2])
        ruler.SetPosition2(endPtCoords[0],endPtCoords[1],endPtCoords[2])

  