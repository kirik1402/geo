import arcpy
import numpy
import os
import dbscan as dbs
import math


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "DBSCANpoints"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [AggregateClusters, CalculateStats]


class AggregateClusters(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aggregate Clusters"
        self.description = "Aggregates cluster points into polygons"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        fc = arcpy.Parameter(
            displayName="Input Clustered Point Features",
            name="fc",
            datatype="GPLayer",
            parameterType="Required",
            direction="Input")

        fld = arcpy.Parameter(
            displayName="Cluster Field",
            name="fld",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        fld.filter.list = ['Short', 'Long']
        fld.parameterDependencies = [fc.name]

        mindist = arcpy.Parameter(
            displayName="Cluster distance",
            name="mindist",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        mindist.value = 10000

        out_polygons = arcpy.Parameter(
            displayName="Output Polygons",
            name="output_polygons",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
        out_polygons.value = 'C:/Users/kcher/Desktop/sem6/GK/Ex02/Ex1.gdb/Aggr_polygons'

        out_points = arcpy.Parameter(
            displayName="Output Points",
            name="out_points",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
        out_points.value = 'C:/Users/kcher/Desktop/sem6/GK/Ex02/Ex1.gdb/Aggr_centroids'

        # params = [fc, minpts, mindist, fld]

        params = [fc, fld, mindist, out_polygons, out_points]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        fc = parameters[0].valueAsText
        fld = str(parameters[1].valueAsText)
        mindist = float(parameters[2].valueAsText)
        out_polygons = parameters[3].valueAsText
        out_points = parameters[4].valueAsText

        arcpy.CreateFeatureclass_management(os.path.dirname(out_polygons),
                                            os.path.basename(out_polygons) + str(int(mindist)), 'POLYGON')
        arcpy.AddField_management(out_polygons + str(int(mindist)), fld, 'LONG')  # номер кластера
        arcpy.AddField_management(out_polygons + str(int(mindist)), 'COUNT', 'LONG')  # количество точек в кластере

        ids = set([row[0] for row in arcpy.da.SearchCursor(fc, fld)])

        ptslyr = 'ptslyr'
        arcpy.MakeFeatureLayer_management(fc, ptslyr)

        for id in ids:
            if id != -1:
                arcpy.SelectLayerByAttribute_management(ptslyr, 'NEW_SELECTION', fld + ' = ' + str(id))

                aggr = 'in_memory/aggr'
                arcpy.AggregatePoints_cartography(ptslyr, aggr, mindist)

                aggr_dis = 'in_memory/aggr_dis'
                arcpy.Dissolve_management(aggr, aggr_dis)

                arcpy.AddField_management(aggr_dis, fld, 'LONG')
                arcpy.AddField_management(aggr_dis, 'COUNT', 'LONG')

                npts = int(arcpy.GetCount_management(ptslyr).getOutput(0))
                arcpy.CalculateField_management(aggr_dis, 'COUNT', npts)
                arcpy.CalculateField_management(aggr_dis, fld, id)

                arcpy.Append_management(aggr_dis, out_polygons + str(int(mindist)))

        arcpy.FeatureToPoint_management(out_polygons + str(int(mindist)), out_points + str(int(mindist)), 'INSIDE')
        return


class CalculateStats(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Stats"
        self.description = "Calculates stats for clustering and regionization"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        clustered_pts = arcpy.Parameter(
            displayName="Input Clustered Point Features",
            name="clustered_pts",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")

        fld = arcpy.Parameter(
            displayName="Cluster Field",
            name="fld",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        fld.filter.list = ['Short', 'Long']
        fld.parameterDependencies = [clustered_pts.name]
        params = [clustered_pts, fld]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        clustered_pts = parameters[0].valueAsText
        fld = parameters[1].valueAsText

        # считываем координаты всех точек в список
        arcpy.AddMessage('Reading coordinates')
        pts = arcpy.da.SearchCursor(clustered_pts, [fld, 'SHAPE@XY'])
        pnt_list = []
        pnt_id = 0
        for pnt in pts:
            if pnt[0] != -1:
                pnt_id += 1
                cluster_id = pnt[0]
                x, y = pnt[1]
                pnt_list.append([pnt_id, cluster_id, x, y])

        # считаем расстояния между парами точек и в зависимости от принадлежности к кластерам добавляем их к
        # суммарному внутрикластерному или суммарному межкластерному расстоянию, также подсчитываем количество
        # внутрикластерных и межкластерных пар
        arcpy.AddMessage('Calculating distances...')
        sum_incluster_dist = 0
        incluster_pair_count = 0
        sum_intercluster_dist = 0
        intercluster_pair_count = 0
        for fst_pnt in pnt_list:
            fst_x = fst_pnt[2]
            fst_y = fst_pnt[3]
            for sec_pnt in pnt_list[fst_pnt[0]:len(pnt_list)]:
                sec_x = sec_pnt[2]
                sec_y = sec_pnt[3]
                dist = math.sqrt((fst_x - sec_x) * (fst_x - sec_x) + (fst_y - sec_y) * (fst_y - sec_y))
                if fst_pnt[1] == sec_pnt[1]:
                    sum_incluster_dist += dist
                    incluster_pair_count += 1
                else:
                    sum_intercluster_dist += dist
                    intercluster_pair_count += 1

        # рассчитываем среднее внутрикластерное и межкластерное расстояния, а также их отношение
        mean_incluster_dist = sum_incluster_dist / incluster_pair_count
        mean_intercluster_dist = sum_intercluster_dist / intercluster_pair_count
        ratio = mean_incluster_dist / mean_intercluster_dist

        # # создаем множество номеров кластеров
        # ids = set([row[0] for row in arcpy.da.SearchCursor(clustered_pts, fld)])
        #
        # cluster_stats_list = []
        # # перебираем все номера кластеров и рассчитываем для нее
        # for id in ids:
        #     if id != -1:
        #         arcpy.AddMessage('Calculating stats for cluster ' + str(id))
        #
        #         # создаем класс в оперативной памяти и заполняем его точками с текущим номером кластера
        #         pts_in_cluster = 'in_memory/pts_in_cluster'
        #         arcpy.Select_analysis(clustered_pts, pts_in_cluster, fld + ' = ' + str(id))
        #
        #         # перебираем точки по одной и вычисляем расстояния до всех других точек этого кластера
        #         links_count = 0
        #         sum_incluster_dists = 0
        #         points = arcpy.da.SearchCursor(pts_in_cluster, ['SHAPE@XY'])
        #         for fst_pnt in points:
        #             for sec_pnt in points:
        #                 fst_x, fst_y = fst_pnt[0]
        #                 sec_x, sec_y = sec_pnt[0]
        #                 dist = math.sqrt((fst_x-sec_x)*(fst_x-sec_x)+(fst_y-sec_y)*(fst_y-sec_y))
        #                 sum_incluster_dists += dist
        #                 links_count += 1
        #
        #         cluster_stats_list.append([links_count, sum_incluster_dists])
        #
        # # рассчитываем среднее внутрикластерное расстояние
        # all_cluster_dists = 0
        # all_links_count = 0
        # for lst in cluster_stats_list:
        #     all_links_count += lst[0]
        #     all_cluster_dists += lst[1]
        # mean_incluster_dist = all_cluster_dists/all_links_count
        #
        # # рассчитываем среднее межкластерное расстояние
        # sum_intercluster_dists = 0
        # centroid_links_count = 0
        # centroid_points = arcpy.da.SearchCursor(centroids, ['SHAPE@XY'])
        # for fst_pnt in centroid_points:
        #     for sec_pnt in centroid_points:
        #         fst_x, fst_y = fst_pnt[0]
        #         sec_x, sec_y = sec_pnt[0]
        #         centroids_dist = math.sqrt((fst_x - sec_x) * (fst_x - sec_x) + (fst_y - sec_y) * (fst_y - sec_y))
        #         sum_intercluster_dists += centroids_dist
        #         centroid_links_count += 1
        # mean_intercluster_dist =

        results_file = open('C:/Users/Rcpil/Desktop/Rez.txt', 'a')
        results_file.write(
            '\n' + str(fld) + '\t' + str(mean_incluster_dist) + '\t' + str(mean_intercluster_dist) + '\t' + str(ratio))
        return
