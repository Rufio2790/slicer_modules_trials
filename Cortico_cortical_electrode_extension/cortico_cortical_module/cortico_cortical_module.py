import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import math
import numpy as np
import copy
from decimal import Decimal
#
# cortico_cortical_module
#

class cortico_cortical_module(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Cortical electrodes"  # TODO make this more human readable by adding spaces
    self.parent.categories = ["Epilepsy"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is a loadable module that allow to generate a cortico-cotical electrode by means of fiducial points selection
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Tiziano Vallisa, Polimi.
and Davide Scorza, PoliMI and Vicomtech-IK4.
""" # replace with organization, grant and thanks.

#
# cortico_cortical_moduleWidget
#

class cortico_cortical_moduleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    self.logic = cortico_cortical_moduleLogic()
    self.layout = self.parent.layout()
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)
    #probabilmente sono da togliere
    # self.tagSourceNode = None
    # self.tagDestinationNode = None

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # Input volume selector
    #
    self.inputSelector1 = slicer.qMRMLNodeComboBox()
    self.inputSelector1.nodeTypes = ["vtkMRMLModelNode"]
    self.inputSelector1.selectNodeUponCreation = True
    self.inputSelector1.addEnabled = True
    self.inputSelector1.removeEnabled = False
    self.inputSelector1.noneEnabled = False
    self.inputSelector1.showHidden = False
    self.inputSelector1.showChildNodeTypes = False
    self.inputSelector1.setMRMLScene(slicer.mrmlScene)
    self.inputSelector1.setToolTip("Select the input volume (lh_pial or rh_pial) with the appropriate scalar curvature" )
    parametersFormLayout.addRow("surface selector: ", self.inputSelector1)
    #
    # Markup fiducial selector
    #
    self.inputSelector2 = slicer.qMRMLNodeComboBox()
    self.inputSelector2.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.inputSelector2.selectNodeUponCreation = True
    self.inputSelector2.addEnabled = True
    self.inputSelector2.removeEnabled = False
    self.inputSelector2.noneEnabled = False
    self.inputSelector2.showHidden = False
    self.inputSelector2.showChildNodeTypes = False
    self.inputSelector2.setMRMLScene(slicer.mrmlScene)
    self.inputSelector2.setToolTip("Select a list of Markup fiducials")
    parametersFormLayout.addRow("Markup selector: ", self.inputSelector2)
    #
    # Strip length selector
    #
    self.StripLengthLayout = qt.QHBoxLayout()
    self.StripLength55 = qt.QRadioButton("4-contact strip")
    self.StripLength75 = qt.QRadioButton("6-contact strip")
    self.StripLengthLayout.addWidget(self.StripLength55)
    self.StripLengthLayout.addWidget(self.StripLength75)
    self.InterpolationGroup = qt.QButtonGroup()
    self.InterpolationGroup.addButton(self.StripLength55)
    self.InterpolationGroup.addButton(self.StripLength75)
    parametersFormLayout.addRow("Strip type selector: ", self.StripLengthLayout)
    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.StripLength55.connect('clicked(bool)', self.onSetLenTo55)
    self.StripLength75.connect('clicked(bool)', self.onSetLenTo75)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector1.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputSelector2.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    self.StripLength55.setChecked(True)
    self.onSetLenTo55(True)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector1.currentNode() and self.inputSelector2.currentNode()

  def onGenerateCurve(self):  # UTILE FORSE
      self.logic.generateCurveOnce()

  def onSetLenTo55(self, s):
      self.logic.strip_option(55.0, 4)

  def onSetLenTo75(self, s):
      self.logic.strip_option(75.0, 6)

  def onApplyButton(self):

    #logic = cortico_cortical_moduleLogic()
    print("run the algorithm")
    model_node = self.inputSelector1.currentNode()
    slicer.modules.markups.logic()
    self.logic.run(self.inputSelector2.currentNode(), model_node.GetPolyData())

#
# cortico_cortical_moduleLogic
#

class cortico_cortical_moduleLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):       #aggiunto
        self.SourceNode = None
        self.DestinationNode = None
        self.TubeRadius = 4.25
        self.StripLen = 0.0
        self.ContactsNumber = 0

        # self.AutomaticUpdate = False
        self.NumberOfIntermediatePoints = 10
        self.ModelColor = [1.0, 1.0, 0.6]

        self.CurvePoly = None
        self.interpResolution = 20

        self.InterpolationMethod = 1

        self.curvature = vtk.vtkCurvatures()
        self.curvatureValues = vtk.vtkDoubleArray()

    def print_coordinate(self, fiducials):
        """
        This function take as input the Markup fiducial list and compute the coordinate in RAS for each fiducial present in the list
        :param fiducials:
        :return: coordinate_list
        """

        coordinate_list = []
        for i in range(fiducials.GetNumberOfFiducials()):
            position = [0.0, 0.0, 0.0]
            fiducials.GetNthFiducialPosition(i, position)
            print(position)
            coordinate_list.append(position)

        return coordinate_list

    def strip_option(self, length, N_of_contacts):
        self.StripLen = length
        self.ContactsNumber = N_of_contacts

    def oriented_buonding_box_locator(self, polydata):
        """
        This function take as input the polydata of a surface and build the object locator to return
        :param polydata:
        :return: obb_locator
        """

        obb_locator = vtk.vtkOBBTree()
        obb_locator.SetDataSet(polydata)
        obb_locator.BuildLocator()

        return obb_locator

    def point_locator(self, polydata):
        """
        This function take as input the polydata of a surface and build the point locator to return
        :param polydata:
        :return: p_locator
        """

        p_locator = vtk.vtkPointLocator()
        p_locator.SetDataSet(polydata)
        p_locator.BuildLocator()

        return p_locator

    def find_closest_points(self, fiducial_coordinate, locator, polydata):
        """
        This function take as input the position of the fiducial/s selected and the point locator. It return the closest
        points belonging to the surface wrt the fiducial/s position
        :param fiducial_coordinate:
        :param locator:
        :return: coord_on_surface
        """
        fiducials_id = []
        coord_on_surface = []

        try:
            len(fiducial_coordinate[0])
            for point in fiducial_coordinate:
                id = locator.FindClosestPoint(point[0], point[1], point[2])
                if id != -1:
                    if polydata:
                        coord = polydata.GetPoint(id)
                        fiducials_id.append(id)
                        coord_on_surface.append(coord)
                    else:
                        fiducials_id.append(id)
                else:
                    raise (Exception, 'Cannot find the closest point on point_locator')
        except:
            id = locator.FindClosestPoint(fiducial_coordinate[0], fiducial_coordinate[1], fiducial_coordinate[2])
            if id != -1:
                if polydata:
                    coord = polydata.GetPoint(id)
                    fiducials_id.append(id)
                    coord_on_surface.append(coord)
                else:
                    fiducials_id.append(id)
            else:
                raise (Exception, 'Cannot find the closest point on point_locator')

        if polydata:
            return fiducials_id, coord_on_surface
        else:
            return fiducials_id

    def find_vertex(self, vertex_id, polydata, array):
        """
        This function take as input the two vertexes (starting and ending position) and return the coordinates of the
        vertexes those connect the initial and final position on a mesh, computing the geodetic distance between the two
        input vertexes
        :param vertex_id:
        :param polydata:
        :return: vertex_list
        """
        vertex_list = []
        P0 = vertex_id[0]
        P1 = vertex_id[1]
        curv_value = 0.0
        temp = []
        #import math
        p0 = [0] * 3
        p1 = [0] * 3
        dist = 0.0
        pp = vtk.vtkPoints()
        vertex_list.append(polydata.GetPoint(P0))
        geodesic_path_filter = vtk.vtkDijkstraGraphGeodesicPath()
        if isinstance(polydata, vtk.vtkPolyData):
            geodesic_path_filter.SetInputData(polydata)

        else:
            geodesic_path_filter.SetInputConnection(polydata.GetOutputPort())

        geodesic_path_filter.SetStartVertex(P0)
        geodesic_path_filter.SetEndVertex(P1)
        geodesic_path_filter.Update()

        pts = geodesic_path_filter.GetOutput().GetPoints()

        for ptId in range(pts.GetNumberOfPoints() - 1):
            pts.GetPoint(ptId, p0)
            pts.GetPoint(ptId + 1, p1)
            dist += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0, p1))
            curv_value = array.GetValue(ptId)
            if curv_value < -0.2:
                temp = copy.copy(p0)
                vertex_list.append(temp)



        return dist, vertex_list  #vertex_list

    def resampling(self, geo_vertex, step, vertex_init_pos):
        """
        This function take as input all the vertex found using the geodetic distance and return a re-sampled list of these
        points by 1:5
        :param geo_vertex:
        :return:
        """
        vertex_resampled = []


        N = len(geo_vertex)
        for i in range(0, N, step):
            vertex_resampled.append(geo_vertex[i])

        vertex_resampled.append(vertex_init_pos[0])

        return vertex_resampled



    def generateCurveOnce(self):
        """
        This function call updateCurve in order to generate the electrode's model in the slicer scene
        :return:
        """

        prevAutomaticUpdate = self.AutomaticUpdate
        self.AutomaticUpdate = True
        self.updateCurve()
        self.AutomaticUpdate = prevAutomaticUpdate

    # def controlPointsUpdated(self,caller,event):
    #  if caller.IsA('vtkMRMLMarkupsFiducialNode') and event == 'ModifiedEvent':
    #    self.updateCurve()

    def nodeToPolyCardinalSpline(self, outputPoly, nOfControlPoints, fiducial_position, loc, poly, array):
        """
        This function take the fiducial control points selected in slicer and generate intermediate points whose will be
        used to interpolate the curve. Before interpolation, the function extract the curvature of each fiducial in order
        to avoid the selection of points lying on sulcis and then pass them to the vtkPoints 'points' for interpolation.
        It returns the point used for interpolation and the length of the curve.
        :param outputPoly: vtkPolyData that will be passed to the tube filter
        :param nOfControlPoints: number of fiducial control points
        :param fiducial_position: coordinates of the fiducial control points
        :param loc: point locator necessary to extract the id of the points
        :param poly: Polydata of the model
        :param array: curvature values for each point of the model passed
        :return:
        """

        new_coo = []


        # One spline for each direction.
        aSplineX = vtk.vtkCardinalSpline()
        aSplineY = vtk.vtkCardinalSpline()
        aSplineZ = vtk.vtkCardinalSpline()

        aSplineX.ClosedOff()
        aSplineY.ClosedOff()
        aSplineZ.ClosedOff()

        for i in range(0, nOfControlPoints):
            pos = fiducial_position[i]
            aSplineX.AddPoint(i, pos[0])
            aSplineY.AddPoint(i, pos[1])
            aSplineZ.AddPoint(i, pos[2])

        # Interpolate x, y and z by using the three spline filters and
        # create new points
        nInterpolatedPoints = (self.interpResolution + 2) * (
                nOfControlPoints - 1)  # One section is devided into self.interpResolution segments
        points = vtk.vtkPoints()
        r = [0.0, 0.0]
        aSplineX.GetParametricRange(r)
        t = r[0]
        tStep = (nOfControlPoints - 1.0) / (nInterpolatedPoints - 1.0)
        nOutputPoints = 0
        length = 0.0
        p = 0
        step = 1.0
        dist = 0.0
        for h in range(nOfControlPoints - 1):
            P_0 = fiducial_position[h]
            points.InsertPoint(p, P_0[0], P_0[1], P_0[2])
            new_coo.append([P_0[0], P_0[1], P_0[2]])
            P_1 = fiducial_position[h+1]
            p = p + 1
            Initial_point = P_0
            # step = 1.0 mm
            N_step = int(math.sqrt(vtk.vtkMath.Distance2BetweenPoints(P_0, P_1)))
            for j in range(N_step):
                dist = step*(j+1)
                useless, new_point = self.set_electrode_length(P_0, P_1, dist)
                ptid, new_close_pnt = self.find_closest_points(new_point, loc, poly)
                curv_val = array.GetValue(ptid[0])
                pnt = new_close_pnt[0]
                if curv_val < 0.1:
                    points.InsertPoint(p, pnt[0], pnt[1], pnt[2])
                    new_coo.append([pnt[0], pnt[1], pnt[2]])
                    length += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(Initial_point, pnt))
                    Initial_point = copy.copy(pnt)
                    p = p + 1
                    if length >= self.StripLen-1:
                        last = fiducial_position[nOfControlPoints-1]
                        points.InsertPoint(p, last[0], last[1], last[2])
                        new_coo.append([last[0], last[1], last[2]])
                        p = p + 1
                        nOutputPoints = p
                        lines = vtk.vtkCellArray()
                        lines.InsertNextCell(nOutputPoints)
                        for i in range(0, nOutputPoints):
                            lines.InsertCellPoint(i)

                        outputPoly.SetPoints(points)
                        outputPoly.SetLines(lines)
                        return new_coo, length
        last = fiducial_position[nOfControlPoints - 1]
        if last not in new_coo:
            points.InsertPoint(p, last[0], last[1], last[2])
            new_coo.append([last[0], last[1], last[2]])
            p = p + 1

        # while t < r[1]:
        #     cord = [aSplineX.Evaluate(t), aSplineY.Evaluate(t), aSplineZ.Evaluate(t)]
        #     ptid, closest_p = self.find_closest_points(cord, loc, poly)
        #     x = closest_p[0]
        #     curv_val = array.GetValue(ptid[0])
        #     if curv_val < 0.0:
        #         points.InsertPoint(p, x[0], x[1], x[2])
        #         new_coo.append([x[0], x[1], x[2]])
        #         length += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(Initial_point, x))
        #         Initial_point = copy.copy(x)
        #         p = p + 1
        #     t = t + tStep
        #     if length >= self.StripLen:
        #         t = r[1]
        nOutputPoints = p
        if length < self.StripLen:
            z = len(new_coo)
            missing_path = self.StripLen - length
            if missing_path != 0.0:
                if new_coo[z - 2] == new_coo[z - 1]:
                    new_last_points, last_point = self.set_electrode_length(new_coo[z - 3], new_coo[z - 1], missing_path)
                    end = len(new_last_points)
                else:
                    new_last_points, last_point = self.set_electrode_length(new_coo[z - 2], new_coo[z - 1], missing_path)
                    end = len(new_last_points)

                for w in range(end-1):
                    temp = new_last_points[w]
                    pid, closest_p = self.find_closest_points(temp, loc, poly)
                    curva = array.GetValue(pid[0])
                    x = closest_p[0]
                    if curva < -0.1:
                        points.InsertPoint(p, x[0], x[1], x[2])
                        new_coo.append([x[0], x[1], x[2]])
                        p = p + 1
            last_id, final_point = self.find_closest_points(last_point, loc, poly)
            final_point = final_point[0]
            points.InsertPoint(p, final_point[0], final_point[1], final_point[2])
            new_coo.append([final_point[0], final_point[1], final_point[2]])
            p = p + 1
            nOutputPoints = p


        lines = vtk.vtkCellArray()
        lines.InsertNextCell(nOutputPoints)
        for i in range(0, nOutputPoints):
            lines.InsertCellPoint(i)

        outputPoly.SetPoints(points)
        outputPoly.SetLines(lines)
        return new_coo, length

    def updateCurve(self, nOfControlPoints, fiducial_position, loc, polydata, arrayValues):
        """
        Principal function that call the PolyCardinalSpline function, to interpolate the points, and the vtkTubeFilter,
        that generates the curve model in slicer.
        :param nOfControlPoints:
        :param fiducial_position:
        :param loc:
        :param polydata:
        :param arrayValues:
        :return:
        """
        self.CurvePoly = vtk.vtkPolyData()
        pot, leng = self.nodeToPolyCardinalSpline(self.CurvePoly, nOfControlPoints, fiducial_position, loc, polydata, arrayValues)

        modelDisplayNode = slicer.vtkMRMLModelDisplayNode()
        modelDisplayNode.SetColor(self.ModelColor)
        modelDisplayNode.SetOpacity(0.4)
        slicer.mrmlScene.AddNode(modelDisplayNode)
        self.DestinationNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLModelNode')
        self.DestinationNode.SetAndObserveDisplayNodeID(modelDisplayNode.GetID())
        tubeFilter = vtk.vtkTubeFilter()


        tubeFilter.SetInputData(self.CurvePoly)
        tubeFilter.SetRadius(self.TubeRadius)
        tubeFilter.SetNumberOfSides(10)
        tubeFilter.CappingOn()
        tubeFilter.Update()


        self.DestinationNode.SetAndObservePolyData(tubeFilter.GetOutput())
        self.DestinationNode.Modified()

        if self.DestinationNode.GetScene() == None:
            slicer.mrmlScene.AddNode(self.DestinationNode)

        displayNode = self.DestinationNode.GetDisplayNode()
        if displayNode:
            displayNode.SetActiveScalarName('')

        return pot, leng

    def display_electrodes(self, pot):
        """
        This function takes points used to generate the curve, computes the distance between them and s the elcetrode
        contacts on the model, based on the type of stip model selected
        :param pot:
        :return:
        """
        displayNode = slicer.vtkMRMLMarkupsDisplayNode()
        displayNode.SetOpacity(1)
        displayNode.SetGlyphScale(5.0)
        slicer.mrmlScene.AddNode(displayNode)
        fidNode = slicer.vtkMRMLMarkupsFiducialNode()
        slicer.mrmlScene.AddNode(fidNode)
        fidNode.SetAndObserveDisplayNodeID(displayNode.GetID())
        cont = 0
        dist_1 = 0.0
        dist_2 = 0.0
        initial_point = pot[3]
        fidNode.AddFiducial(initial_point[0], initial_point[1], initial_point[2])
        fidNode.SetNthFiducialLabel(cont, '{:d}'.format(1))
        cont = cont + 1
        while cont < self.ContactsNumber - 1:
            for i in range(4, len(pot)-1):
                p_1 = copy.copy(pot[i])
                p_2 = copy.copy(pot[i+1])
                dist_1 = copy.copy(dist_2)
                dist_2 += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p_1, p_2))
                if 9.98 < dist_2 < 10.11:
                    fidNode.AddFiducial(p_2[0], p_2[1], p_2[2])
                    fidNode.SetNthFiducialLabel(cont, '{:d}'.format(cont+1))
                    dist_2 = 0.0
                    cont += 1
                    if cont == self.ContactsNumber:
                        return True
                else:
                    if dist_2 > 10.11:
                        miss = 10.0 - dist_1
                        not_useful, electrode_pos = self.set_electrode_length(p_1, p_2, miss)
                        fidNode.AddFiducial(electrode_pos[0], electrode_pos[1], electrode_pos[2])
                        fidNode.SetNthFiducialLabel(cont, '{:d}'.format(cont+1))
                        dist_2 = 0.0
                        cont += 1
                        if cont == self.ContactsNumber:
                            return True

        return True

    def set_electrode_length(self, P_0, P_1, dist):
        """
        This function is called only if the length of the strip model is shorter than the one
        :param P_0:
        :param P_1:
        :param dist:
        :return:
        """
        new_fiducial_pos = []

        distance = [P_1[0]-P_0[0], P_1[1]-P_0[1], P_1[2]-P_0[2]]
        # norm = np.linalg.norm(distance) #math.sqrt(distance[0]**2 + distance[1]**2 + distance[2]**2)
        norm = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(P_1, P_0))
        print('NORM: ')
        print(norm)
        print("--------")
        if norm == 0.0:
            # new_fiducial_pos.append(P_0)
            new_fiducial_pos.append(P_1)
            return new_fiducial_pos, P_1
        direction = [distance[0] / norm, distance[1] / norm, distance[2] / norm]
        # The problem is probably here when compute the new end point

        # new_fiducial_pos.append(P_0)
        new_end_point = [P_0[0] + direction[0]*dist, P_0[1] + direction[1]*dist, P_0[2] + direction[2]*dist]
        offset = 0.10
        dist = float(dist)
        for i in range(10):
            space = float(offset*dist)
            point = [P_0[0] + direction[0]*space, P_0[1] + direction[1]*space, P_0[2] + direction[2]*space]
            new_fiducial_pos.append(point)
            offset = offset + 0.1
        new_fiducial_pos.append(new_end_point)

        return new_fiducial_pos, new_end_point

    def run(self, fiducial, polydata):
        """
        run the algorithm
        """
        arr_name = []
        fiducial_position = self.print_coordinate(fiducial)
        print("--------")
        fiducial.GetDisplayNode().SetVisibility(False)
        for i in range(polydata.GetPointData().GetNumberOfArrays()):
            arr_name.append(polydata.GetPointData().GetArrayName(i))

        if len(arr_name) > 1 and 'curv' in arr_name[1]:
            arrayValues = polydata.GetPointData().GetArray(arr_name[1])
        else:
            raise NameError('You must import the scalar curvature values as well')
        # arr_name = polydata.GetPointData().GetArrayName(1)
        # if arr_name is None or 'curv' not in arr_name:

        # arrayValues = polydata.GetPointData().GetArray(arr_name)
        # new_fiducials = self.set_electrode_length(fiducial_position[0], fiducial_position[1], self.StripLen)
        # improved_fiducial = self.add_mid_point(new_fiducials)
        p_locator = self.point_locator(polydata)
        point_id, point_on_surface = self.find_closest_points(fiducial_position, p_locator, polydata)
        # print("point on the surface corresponding to the fiducials selected: ")
        # self.compute_curvature(polydata)
        # distance, trajectory_vertex = self.find_vertex(point_id, polydata, arrayValues)
        print('SURFACE POINTS')
        print(point_on_surface)



        numb_of_control_points = len(point_on_surface)
        print(numb_of_control_points)
        # print(trajectory_vertex)
        # temp_coord, length = self.compute_first_interpolation(numb_of_control_points, point_on_surface, p_locator, polydata, True)
        # new_coord = self.resampling(temp_coord, 20, point_on_surface)
        # new_vertex_list = self.resampling(new_coord, 3, point_on_surface)
        # N_control = len(temp_coord)

        pot, length = self.updateCurve(numb_of_control_points, point_on_surface, p_locator, polydata, arrayValues)
        self.display_electrodes(pot)
        print('-------')
        print(length)
        print('-------')
        print(pot)


        return True


