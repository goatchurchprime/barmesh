# Introduction #

BarMesh is a self-contained set of Python scripts that implement geometric modelling features useful for CAM and 3D printing.  

This includes functions for slicing STL files that are not watertight for a given tool shape to a required tolerance.

This code is intended as the basis for optimization experiments using technologies such as Pypy, PyOpenCL and Cython that will hopefully out-perform the native C++ applications that are commercially available, but in an open source environment where innovation is more likely.


# Getting started #

Go to command line and type:

* `python main.py --stl=stlsamples/barobox.stl -v -r2.5 -n12`

Outputs a G-code tape file barobox.tap (suggest a better format), and works with Python2.7, Python3 and Pypy.

See [freesteel blogpost](http://www.freesteel.co.uk/wpblog/2015/01/05/we-have-some-slicing/) for images and current state of play.


# License #

This code is released under [GPL3](http://www.gnu.org/).


# Internal structure #

The point structure in basicgeo.py is `P3(namedtuple('P3', ['x', 'y', 'z']))`.  Each Node has an index i and point P3 p.

Each directed Bar has a nodeFore`NF`, nodeBack`NB`, barForeRight`BFR`, barBackLeft`BBL` which allows full iteration around the individual cells with an arbitrary number of edges.

```
        _____NF
            /|
           ^ |
          /  |BFR
    |    /-> |
    | <-/    |
 BBL|  /     |
    | / 
    |/___
   NB
```

# Visualization #

Use [twistcodewiki](https://bitbucket.org/goatchurch/twistcodewiki) which allows the geometry from intermediate stages of the algorithm to be viewed in WebGL via the sendactivity() function.

# Todo checklist #

* check if the default slicer still works [DONE]
* get the voronoi making toolpaths again  [DONE]
* shift more vorpicklework to postprocess.py and make the transfers of structures more interesting [DONE]
* enable this to spiral down on a toothpick on a steel thing to pluck

Voronoi

* Drive geometry corner points to be allocated to vbbarmesh CellMark and unambiguiously split on every cell subdivision (use izone number)
* vornoi diagram to have its own barmesh as well
* use above allocation to inform DetectMissingThreeway of nodes (which are threeways) not exhibited on contour
* get voronoi to run in self-contained command line for profiling
* the InsertColPointNodes point splitting should do splitting on those lines that cross the contour!!!
* detect 3-way bisection points signifying a missing node which is on the inside of the cell [DONE]
* fix excessive subdivision along the constant Y line at the corner node!!!
* handle too short subdivision case properly
* outside and inside pointzone bit (could record with the sign of r, accounting for external quadrant in case of corner point
* Dist on arc to give +- distance 
* it will be between rtip-z0 and rdepthz-zd forming 3-way contours everywhere, but oriented by plunging and climb milling 
* voronoi creation to show better error if incomplete when called
* handle self-intersecting and infinitely short line segments, able to conceive of lower bound solution of problems
* aim to get this to work no matter what the data


* Find what proportion of time is spent in the MakePointZoneVBS function to see if it needs optimization [DONE]
* convert all drive values to P2 [DONE]
* voronoi polygon sides sequenced and properly implemented as conics within these cells [DONE]
* use pbdir on CalcCutBisect2Obj to take aim between a pair of divergent parallels to split the 3-way nodes [DONE]
* extract the voronoi structure into a V-groove toolpath to try out immediately accounting for inside and outside [DONE]
* this is vgroove down to a depth and track along the contours plunging in all directions, so find midpoints too [DONE]
* get the stack to retain some kind of order if possible [DONE]
* detect and implement special case where midpoint3 is a boundary node [DONE]
* all cases of disabled fcl4
* multiple contours in same diagram [DONE]
* complete the leathersatchel case [DONE]
* handle the skipping too short case in subdivisions of the conts[2] case
* single line subdivide demonstration [DONE]
* single line subdivide demonstration with adjancent lines predicted from line sides [DONE]
* example with arc requiring LineArcBisectCut [DONE]
* voronoi to handle arcs (make the vbnodes actually into zero rad arcs so that all componenets are tangential with dividing perp lines, then get rid of the k2 business) [DONE]
* ensure all transections are found of a line by repeated subdividing on the dividing node to find that it is inside a 3rd vbpoly [DONE]
* implement subdividing of cells that have 4 or more way vbpolys not only when 2 share the same [DONE]
* implement VBbarmeshcell whole thing possibly using cellmark as signals [DONE]
* implement LineArc0BisectCut [DONE]
* implement 3 way bisection points [DONE]
* implement 3 way bisection points from tangent vert bisectors [DONE]
* implement Arc0Arc0BisectCut [DONE]
* voronoi cell subdivisions around cells with multiple cuts [DONE]
* the voronoi experiment should use pointzones instead of vddicts [DONE]

Slicer
* slicing to work on huge pointclouds [DONE]
* silhouette boundaries 2D barmesh, see pointcloudssliceroffset for application
* slices of pointclouds to join up into STL contours
* slicing to work on huge connected 2D contours to set inside and outside
* deal with degenerate cases of nodes and edges aligned with barmesh
* use I1 range in triangle and barmesh structures
* assign vertical connectivity between adjacent contour pockets in z
* function to detect and remove contours that are inside the cell by vertical connectivity and lack of connectivity to the external world
* make new TriangleGroupFuncs that re-states the triangle mesh but with normals
* upgrade splitbarsdirectionchangesR to find new gaps where they are both not PZ_BEYOND_R to find narrow channels
* extend slice to work offsets of polygons
* demo of bumpy textures and spongey linear structures applied to the slicing (take to ChrisS)
* 2D voronoi diagrams recalculated
* attempt the drop and wrap fabric feature
* raster passes with subdivisions implemented and cylinder projection down z
* toolpath output in G-code with rotation applied
* adaptive scallop
* toroidal tool
* investigate https://developers.google.com/protocol-buffers/
* try profiling the slicer to find which function to Cython low level
* investigate OpenCL GPU calculations done on tiles of points MakePointZoneRFS and CutbarRFS from multiple levels simulateously at a time so we can iterate through all the triangles in one section over a whole set of close-positioned points at a time

* slicing to sort out Dbintolerance issues for organizing to map back in the cellmark areas we need to go to again [DONE]
* horizon slice bar cuts implemented [DONE]
* nesting of the contour areas into pockets and islands [DONE]
* move triangle set into its own barmesh structure?  [DONE]
* investigate cython [DONE] (too C++-y)



