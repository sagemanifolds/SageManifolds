r"""
Points on a manifold

The class :class:`Point` implements the concept of point on a manifold, in a 
coordinate independent manner: a :class:`Point` object can have coordinates in 
various charts defined on the manifold. Two points are declared equal if they 
have the same coordinates in the same chart. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

EXAMPLES: 

    Defining a point on `\RR^3` by its spherical coordinates::
    
        sage: m = Manifold(3, 'R3', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
        sage: p = Point(m, (1, pi/2, 0), 'spher')
        sage: p1 = Point(m, (1, pi/2, 0))  # using the default chart

    Points can be compared::
    
        sage: p == p1
        True

    Computing the coordinates of the point in a new chart::
    
        sage: c_cart = Chart(m, 'x y z', 'cart')        
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: p.change_coord('cart') # evaluate p's Cartesian coord.
        sage: p.coord('cart')
        (1, 0, 0)
    
    Listing all the coordinates of a point in different charts::
    
        sage: p.coordinates
        {'spher': (1, 1/2*pi, 0), 'cart': (1, 0, 0)}

"""

#*****************************************************************************
#       Copyright (C) 2013 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from manifold import Manifold
from chart import Chart
     
class Point(SageObject):
    r"""
    Class for points on a manifold.

    INPUT:
    
    - ``manifold`` -- the manifold to which the point belongs
    - ``coords`` -- (default: None) the point coordinates (as a tuple or a list)
    - ``chartname`` -- (default: None) string defining the chart in which the 
      coordinates are given; if none is provided, the coordinates are assumed 
      to refer to the manifold's default chart
    
    """
    def __init__(self, manifold, coords=None, chartname=None): 
        self.manifold = manifold
        if coords is None: self.coordinates = {}
        else:
            if len(coords) != self.manifold.dim: 
                raise ValueError("The number of coordinates must be equal" +
                                 " to the manifold dimension.")
            if chartname is None: chartname1 = self.manifold.def_chart.name
            else: 
                if chartname in self.manifold.atlas: chartname1 = chartname
                else:
                    raise ValueError("The chart " + chartname +
                     " has not been defined on the manifold " + str(self.manifold))
            self.coordinates = { chartname1: coords }
            
   
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "point on " + str(self.manifold)
        return description

    def coord(self, chartname=None):
        r"""
        Returns the point coordinates in the specified chart.

        INPUT:
    
        - ``chartname`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are 
          assumed to refer to the manifold's default chart

        EXAMPLES: 

        Coordinates of a point on `\RR^3`::

            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
            sage: p = Point(m, (1, pi/2, 0), 'spher')
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, 0)
            sage: p.coord('spher') # with the chart specified
            (1, 1/2*pi, 0)
        """
        if chartname is None: chartname = self.manifold.def_chart.name
        return self.coordinates[chartname]
        
    def set_coord(self, coords, chartname=None, delete_others=True):
        r"""
        Sets the point coordinates in the specified chart.

        INPUT:
        
        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chartname`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are
          assumed to refer to the manifold's default chart
        - ``delete_others`` -- (default:True) if True, the coordinates defined
          in other charts are  deleted, otherwise their values are kept. In the
          later case there is no check that the coordinates in various charts 
          are compatible.
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart = Chart(m, 'x y z', 'cart')        
            sage: p = Point(m)
            sage: p.set_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
            sage: q = Point(m)
            sage: q.set_coord((1,2,3), 'spher')
            sage: cart_from_spher = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            sage: q.set_coord( cart_from_spher(*q.coord('spher')), 'cart', delete_others=False)
            sage: q.coordinates
            {'spher': (1, 2, 3), 'cart': (sin(2)*cos(3), sin(2)*sin(3), cos(2))}
            sage: p.coordinates
            {'cart': (1, 2, 3)}
            sage: p == q  # p and q should differ because the coordinates (1,2,3) are on different charts
            False     
            
        By default, coordinates set in other charts are deleted::
        
            sage: p = Point(m, (1,2,3), 'spher')
            sage: p.coordinates
            {'spher': (1, 2, 3)}
            sage: p.set_coord((4,5,6), 'cart')
            sage: p.coordinates
            {'cart': (4, 5, 6)}
            sage: p.set_coord((7,8,9), 'spher', delete_others=False)
            sage: p.coordinates
            {'spher': (7, 8, 9), 'cart': (4, 5, 6)}
            
        """
        if len(coords) != self.manifold.dim: 
            raise ValueError("The number of coordinates must be equal " + 
                             "to the manifold dimension.")
        if chartname is None: chartname = self.manifold.def_chart.name
        else: 
            if chartname not in self.manifold.atlas:
                raise ValueError("The chart " + chartname +
                    " has not been defined on the manifold " + str(self.manifold))
        if delete_others: self.coordinates.clear()
        self.coordinates[chartname] = coords

    def change_coord(self, chartname_new, chartname_old=None):
        r"""
        Computes the point coordinates in the specified chart from the 
        coordinates in a previous chart.
        
        INPUT:
        
        - ``chartname_new`` -- string defining the chart in which the new 
          coordinates are to be computed 
        - ``chartname_old`` -- (default: None) string defining the chart in 
          which the "old" coordinates are given; if none is provided, the 
          manifold's default chart will be used
            

        EXAMPLES: 

        Computing the Cartesian coordinates of a point of `\RR^3` from its 
        spherical ones::

            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
            sage: c_cart = Chart(m, 'x y z', 'cart')        
            sage: CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z) 
            sage: p = Point(m, (1,pi/2,0)) # p is defined by its spherical coord. (default chart)
            sage: p.change_coord('cart')  # from the default chart ('spher')
            sage: p.coord('spher')
            (1, 1/2*pi, 0)
            sage: p.coord('cart')
            (1, 0, 0)
            
        """
        if chartname_old is None: chartname_old = self.manifold.def_chart.name
        chcoord = self.manifold.coord_changes[(chartname_old, chartname_new)]
        self.coordinates[chartname_new] = \
                                    chcoord(*self.coordinates[chartname_old])


    def __eq__(self, another_point):
        r"""
        Compares the current point with another one.
        """
        notfound = True
        i1 = 0 
        while notfound and i1 < len(self.coordinates):
            chartname = self.coordinates.keys()[i1]
            notfound = chartname not in another_point.coordinates
            i1 += 1 
        if notfound: 
            raise ValueError("No common chart found for the coordinates of " +
                             "the two points.")
        else:
            return self.coordinates[chartname] == \
                another_point.coordinates[chartname]
        


