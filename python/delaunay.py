# MIT License; Copied from:
# http://code.activestate.com/recipes/579021-delaunay-triangulation/
#
# (I've had trouble with scipy.spatial,
# so going with pure Python even if it's slower.)

# Ignore long lines:
# flake8: noqa: E501

import numpy
import math
import copy


class DictDelaunay2d:
    '''
    Wrapper for Delaunay2d: Instead of supplying a list and getting indices, this takes a dict and
    returns keys for the dict.

    Start with these points:
      D   F
      B
        C E
    A

    >>> points = {'A':[0,0], 'B':[1,2], 'C':[2,1], 'D':[1,3], 'E':[3,1], 'F':[3,3]}
    >>> delaunay = DictDelaunay2d(points)
    >>> delaunay.getTriangles()
    [['B', 'C', 'F'], ['B', 'F', 'D'], ['C', 'F', 'E'], ['B', 'A', 'C'], ['C', 'A', 'E'], ['B', 'D', 'A']]

    '''
    def __init__(self, pairs_dict):
        items = pairs_dict.items()
        self.keys = [item[0] for item in items]
        pairs = [item[1] for item in items]
        self.delaunay = Delaunay2d(pairs)

    def getTriangles(self):
        triangles = self.delaunay.getTriangles()
        return [[self.keys[point] for point in tri] for tri in triangles]


class Delaunay2d:
    '''
    The algorithm uses the S-hull method by D A Sinclair (http://www.s-hull.org/paper/s_hull.pdf).
    The method involves ordering the points in increasing distance from the cloud's center of
    gravity, creating a triangle with the first three points, and adding the remaining points while
    contructing triangles between the point and boundary edges - only triangles with a definite sign
    need to be added (the edge must be visible). Finally, an edge flipping step ensures that the
    triangles are well formed, i.e. the sum of opposite angles to an edge is < 180 degree (the
    so-called Delaunay criterion). It takes about 30 seconds to triangulate 1000 points, most of the
    time is spent in two methods (flipOneEdge and getArea).

    Further improvements could include adding support for holes and automatically adding Steiner
    points. Generally speaking, triangles which have an edge on the boundary (convex hull), tend to
    be squashed; this problem can be alleviated by adding more points on the boundary. Holes can be
    treaded by removing triangles, edges, and points that are contained inside a given closed path.

    Start with these points:
      3   5
      1
        2 4
    0

    >>> points = [[0,0], [1,2], [2,1], [1,3], [3,1], [3,3]]
    >>> delaunay = Delaunay2d(points)
    >>> [list(point) for point in delaunay.points]
    [[1, 2], [2, 1], [1, 3], [3, 1], [3, 3], [0, 0]]

    Points have been reordered:
      2   4
      0
        1 3
    5

    >>> delaunay.triangles
    [[0, 1, 4], [0, 4, 2], [1, 4, 3], [0, 5, 1], [1, 5, 3], [0, 2, 5]]

    Internal order!

    >>> delaunay.getTriangles()
    [[1, 2, 5], [1, 5, 3], [2, 5, 4], [1, 0, 2], [2, 0, 4], [1, 3, 0]]

    Original order!
    '''

    EPS = 1.23456789e-14

    def __init__(self, pairs):

        # data structures
        self.points = [numpy.array(pair) for pair in pairs]
        self.triangles = []  # cells
        self.edge2Triangles = {}  # edge to triangle(s) map
        self.boundaryEdges = set()
        self.appliedBoundaryEdges = None
        self.holes = None

        # compute center of gravity
        cg = numpy.zeros((2,), numpy.float64)
        for pt in self.points:
            cg += pt
        cg /= len(self.points)

        # sort
        def distanceSquare(pt):
            d = pt - cg
            return numpy.dot(d, d)
        dSqFromCenter = numpy.array([distanceSquare(pt) for pt in self.points])
        self.order = list(dSqFromCenter.argsort())
        self.points = [self.points[i] for i in self.order]

        # create first triangle, make sure we're getting a non-zero area otherwise
        # drop the points
        area = 0.0
        index = 0
        stop = False
        while not stop and index + 2 < len(self.points):
            area = self.getArea(index, index + 1, index + 2)
            if abs(area) < self.EPS:
                del self.points[index]
            else:
                stop = True
        if index <= len(self.points) - 3:
            tri = [index, index + 1, index + 2]
            self.makeCounterClockwise(tri)
            self.triangles.append(tri)

            # boundary edges
            e01 = (tri[0], tri[1])
            self.boundaryEdges.add(e01)
            e12 = (tri[1], tri[2])
            self.boundaryEdges.add(e12)
            e20 = (tri[2], tri[0])
            self.boundaryEdges.add(e20)

            e01 = self.makeKey(e01[0], e01[1])
            self.edge2Triangles[e01] = [0, ]

            e12 = self.makeKey(e12[0], e12[1])
            self.edge2Triangles[e12] = [0, ]

            e20 = self.makeKey(e20[0], e20[1])
            self.edge2Triangles[e20] = [0, ]

        else:
            # all the points fall on a line
            return

        # add additional points
        for i in range(3, len(self.points)):
            self.addPoint(i)

        # remove all triangles inside holes
        # TO DO

    def getTriangles(self):
        """
        @return triangles
        """
        return [[self.order[point] for point in tri] for tri in self.triangles]

    def getEdges(self):
        """
        @return egdes
        """
        return self.edge2Triangles.keys()

    def getArea(self, ip0, ip1, ip2):
        """
        Compute the parallelipiped area
        @param ip0 index of first vertex
        @param ip1 index of second vertex
        @param ip2 index of third vertex
        """
        d1 = self.points[ip1] - self.points[ip0]
        d2 = self.points[ip2] - self.points[ip0]
        return (d1[0]*d2[1] - d1[1]*d2[0])

    def isEdgeVisible(self, ip, edge):
        """
        Return true iff the point lies to its right when the edge points down
        @param ip point index
        @param edge (2 point indices with orientation)
        @return True if visible
        """
        area = self.getArea(ip, edge[0], edge[1])
        if area < self.EPS:
            return True
        return False

    def makeCounterClockwise(self, ips):
        """
        Re-order nodes to ensure positive area (in-place operation)
        """
        area = self.getArea(ips[0], ips[1], ips[2])
        if area < -self.EPS:
            ip1, ip2 = ips[1], ips[2]
            # swap
            ips[1], ips[2] = ip2, ip1

    def flipOneEdge(self, edge):
        """
        Flip one edge then update the data structures
        @return set of edges to interate over at next iteration
        """

        # start with empty set
        res = set()

        # assume edge is sorted
        tris = self.edge2Triangles.get(edge, [])
        if len(tris) < 2:
                # nothing to do, just return
            return res

        iTri1, iTri2 = tris
        tri1 = self.triangles[iTri1]
        tri2 = self.triangles[iTri2]

        # find the opposite vertices, not part of the edge
        iOpposite1 = -1
        iOpposite2 = -1
        for i in range(3):
            if not tri1[i] in edge:
                iOpposite1 = tri1[i]
            if not tri2[i] in edge:
                iOpposite2 = tri2[i]

        # compute the 2 angles at the opposite vertices
        da1 = self.points[edge[0]] - self.points[iOpposite1]
        db1 = self.points[edge[1]] - self.points[iOpposite1]
        da2 = self.points[edge[0]] - self.points[iOpposite2]
        db2 = self.points[edge[1]] - self.points[iOpposite2]
        crossProd1 = self.getArea(iOpposite1, edge[0], edge[1])
        crossProd2 = self.getArea(iOpposite2, edge[1], edge[0])
        dotProd1 = numpy.dot(da1, db1)
        dotProd2 = numpy.dot(da2, db2)
        angle1 = abs(math.atan2(crossProd1, dotProd1))
        angle2 = abs(math.atan2(crossProd2, dotProd2))

        # Delaunay's test
        if angle1 + angle2 > math.pi*(1.0 + self.EPS):

            # flip the triangles
            #                         / ^ \                                        / b \
            # iOpposite1 + a|b + iOpposite2    =>     + - > +
            #                         \     /                                        \ a /

            newTri1 = [iOpposite1, edge[0], iOpposite2]  # triangle a
            newTri2 = [iOpposite1, iOpposite2, edge[1]]  # triangle b

            # update the triangle data structure
            self.triangles[iTri1] = newTri1
            self.triangles[iTri2] = newTri2

            # now handle the topolgy of the edges

            # remove this edge
            del self.edge2Triangles[edge]

            # add new edge
            e = self.makeKey(iOpposite1, iOpposite2)
            self.edge2Triangles[e] = [iTri1, iTri2]

            # modify two edge entries which now connect to
            # a different triangle
            e = self.makeKey(iOpposite1, edge[1])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri1:
                    v[i] = iTri2
            res.add(e)

            e = self.makeKey(iOpposite2, edge[0])
            v = self.edge2Triangles[e]
            for i in range(len(v)):
                if v[i] == iTri2:
                    v[i] = iTri1
            res.add(e)

            # these two edges might need to be flipped at the
            # next iteration
            res.add(self.makeKey(iOpposite1, edge[0]))
            res.add(self.makeKey(iOpposite2, edge[1]))

        return res

    def flipEdges(self):
        """
        Flip edges to statisfy Delaunay's criterion
        """

        # start with all the edges
        edgeSet = set(self.edge2Triangles.keys())

        continueFlipping = True

        while continueFlipping:

            #
            # iterate until there are no more edges to flip
            #

            # collect the edges to flip
            newEdgeSet = set()
            for edge in edgeSet:
                # union
                newEdgeSet |= self.flipOneEdge(edge)

            edgeSet = copy.copy(newEdgeSet)
            continueFlipping = (len(edgeSet) > 0)

    def addPoint(self, ip):
        """
        Add point
        @param ip point index
        """

        # collection for later updates
        boundaryEdgesToRemove = set()
        boundaryEdgesToAdd = set()

        for edge in self.boundaryEdges:

            if self.isEdgeVisible(ip, edge):

                # create new triangle
                newTri = [edge[0], edge[1], ip]
                newTri.sort()
                self.makeCounterClockwise(newTri)
                self.triangles.append(newTri)

                # update the edge to triangle map
                e = list(edge[:])
                e.sort()
                iTri = len(self.triangles) - 1
                self.edge2Triangles[tuple(e)].append(iTri)

                # add the two boundary edges
                e1 = [ip, edge[0]]
                e1.sort()
                e1 = tuple(e1)
                e2 = [edge[1], ip]
                e2.sort()
                e2 = tuple(e2)
                v1 = self.edge2Triangles.get(e1, [])
                v1.append(iTri)
                v2 = self.edge2Triangles.get(e2, [])
                v2.append(iTri)
                self.edge2Triangles[e1] = v1
                self.edge2Triangles[e2] = v2

                # keep track of the boundary edges to update
                boundaryEdgesToRemove.add(edge)
                boundaryEdgesToAdd.add((edge[0], ip))
                boundaryEdgesToAdd.add((ip, edge[1]))

        # update the boundary edges
        for bedge in boundaryEdgesToRemove:
            self.boundaryEdges.remove(bedge)
        for bedge in boundaryEdgesToAdd:
            bEdgeSorted = list(bedge)
            bEdgeSorted.sort()
            bEdgeSorted = tuple(bEdgeSorted)
            if len(self.edge2Triangles[bEdgeSorted]) == 1:
                # only add boundary edge if it does not appear
                # twice in different order
                self.boundaryEdges.add(bedge)

        # recursively flip edges
        flipped = True
        while flipped:
            flipped = self.flipEdges()

    def makeKey(self, i1, i2):
        """
        Make a tuple key such at i1 < i2
        """
        if i1 < i2:
            return (i1, i2)
        return (i2, i1)

    def show(self, width=500, height=500, showVertices=False, showCells=False, showContour=[]):
        import tkinter

        xmin = min([p[0] for p in self.points])
        ymin = min([p[1] for p in self.points])
        xmax = max([p[0] for p in self.points])
        ymax = max([p[1] for p in self.points])
        padding = 5
        w = width - 2*padding
        h = height - 2*padding

        master = tkinter.Tk()
        c = tkinter.Canvas(master, width=width, height=height)
        c.pack()
        for e in self.edge2Triangles:
          i1, i2 = e
          xp1 = padding + int(w*(self.points[i1][0] - xmin)/(xmax - xmin))
          yp1 = padding + int(h*(ymax - self.points[i1][1])/(ymax - ymin))
          xp2 = padding + int(w*(self.points[i2][0] - xmin)/(xmax - xmin))
          yp2 = padding + int(h*(ymax - self.points[i2][1])/(ymax - ymin))
          c.create_line(xp1, yp1, xp2, yp2)

        if showVertices:
          for i in range(len(self.points)):
            xp = padding + int(w*(self.points[i][0] - xmin)/(xmax - xmin))
            yp = padding + int(h*(ymax - self.points[i][1])/(ymax - ymin))
            c.create_text(xp, yp, text=str(i))

        if showCells:
          for tId, tVals in self.triangles.items():
            cg = reduce(operator.add, [self.points[i] for i in tVals])/float(len(tVals))
            xp = padding + int(w*(cg[0] - xmin)/(xmax - xmin))
            yp = padding + int(h*(ymax - cg[1])/(ymax - ymin))
            c.create_text(xp, yp, text=str(tId))

        if len(showContour) > 0:
          for i in range(len(showContour) - 1):
            xp1 = padding + int(w*(showContour[i][0] - xmin)/(xmax - xmin))
            yp1 = padding + int(h*(ymax - showContour[i][1])/(ymax - ymin))
            xp2 = padding + int(w*(showContour[i+1][0] - xmin)/(xmax - xmin))
            yp2 = padding + int(h*(ymax - showContour[i+1][1])/(ymax - ymin))
            c.create_line(xp1, yp1, xp2, yp2, fill='red')


        tkinter.mainloop()
