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
        sage: c_spher = Chart(m, r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        sage: p = Point(m, (1, pi/2, 0), name='P') # coordinates in the manifold's default chart
        sage: p
        point 'P' on 3-dimensional manifold 'R3'
        sage: latex(p) 
        P

    Computing the coordinates of the point in a new chart::
    
        sage: c_cart = Chart(m, 'x y z', 'cart')        
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: p.coord('cart') # evaluate P's Cartesian coordinates
        (1, 0, 0)
    
    Points can be compared::
    
        sage: p1 = Point(m, (1, pi/2, 0))
        sage: p == p1
        True
        sage: q = Point(m, (1,2,3), 'cart', name='Q') # point defined by its Cartesian coordinates
        sage: p == q
        False

    Listing all the coordinates of a point in different charts::
    
        sage: p.coordinates
        {'spher': (1, 1/2*pi, 0), 'cart': (1, 0, 0)}
        sage: q.coordinates
        {'cart': (1, 2, 3)}
        
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
    - ``chart_name`` -- (default: None) string defining the chart in which the 
      coordinates are given; if none is provided, the coordinates are assumed 
      to refer to the manifold's default chart
    - ``name`` -- (default: None) name given to the point
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the point; if 
      none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    A point on a 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy = Chart(m, 'x y', 'xy')
        sage: (a, b) = var('a b') # generic coordinates for the point
        sage: p = Point(m, (a, b), name='P') ; p
        point 'P' on 2-dimensional manifold 'M'
        sage: p.coord()  # coordinates of P in the manifold's default chart
        (a, b)
        
    By default, the LaTeX symbol of the point is deduced from its name::
    
        sage: latex(p)
        P
        
    But it can be set to any value::
    
        sage: p = Point(m, (a, b), name='P', latex_name=r'\mathcal{P}')
        sage: latex(p)
        \mathcal{P}
    
    """
    def __init__(self, manifold, coords=None, chart_name=None, name=None, 
                 latex_name=None): 
        self.manifold = manifold
        if coords is None: 
            self.coordinates = {}
        else:
            if len(coords) != self.manifold.dim: 
                raise ValueError("The number of coordinates must be equal" +
                                 " to the manifold dimension.")
            if chart_name is None: 
                chart_name = self.manifold.def_chart.name
            else: 
                if chart_name not in self.manifold.atlas: 
                    raise ValueError("The chart " + chart_name +
                                " has not been defined on the " + 
                                str(self.manifold))
            self.coordinates = { chart_name: coords }
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "point"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on " + str(self.manifold)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{no symbol}'
        else:
           return self.latex_name

    def coord(self, chart_name=None, old_chart=None):
        r"""
        Return the point coordinates in the specified chart.

        If these coordinates are not already known, they are computed from 
        known ones by means of change-of-chart formulas. 

        INPUT:
    
        - ``chart_name`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are 
          assumed to refer to the manifold's default chart
        - ``old_chart`` -- (default: None) string defining the chart from 
          which the coordinates in ``chart_name`` are to be computed. If None, 
          a chart in which the point's coordinates are already known will be
          picked, priveleging the manifold's default chart.

        EXAMPLES: 

        Spherical coordinates of a point on `\RR^3`::

            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher = Chart(m, r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher') # spherical coordinates
            sage: p = Point(m, (1, pi/2, 0)) 
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, 0)
            sage: p.coord('spher') # with the chart 'spher' specified (same result as above since this is the default chart)
            (1, 1/2*pi, 0)

        Computing the Cartesian coordinates from the spherical ones::

            sage: c_cart = Chart(m, 'x y z', 'cart')  # Cartesian coordinates   
            sage: CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z) 
            sage: p.coord('cart')  # the computation is performed by means of the above change of coordinates
            (1, 0, 0)

        Coordinates of a point on a 2-dimensional manifold::
    
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: (a, b) = var('a b') # generic coordinates for the point
            sage: p = Point(m, (a, b), name='P')
            sage: p.coord()  # coordinates of P in the manifold's default chart
            (a, b)
            
        Coordinates of P in a new chart::
        
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: ch_xy_uv = CoordChange(c_xy, c_uv, x-y, x+y)
            sage: p.coord('uv')
            (a - b, a + b)

        Coordinates of P in a third chart::
        
            sage: c_wz = Chart(m, 'w z', 'wz')
            sage: ch_uv_wz = CoordChange(c_uv, c_wz, u^3, v^3)   
            sage: p.coord('wz', old_chart='uv')
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)

        Actually, in the present case, it is not necessary to specify 
        old_chart='uv'::
        
            sage: p.set_coord((a-b, a+b), 'uv') # erases all the coordinates except those in the chart 'uv'
            sage: p.coordinates                
            {'uv': (a - b, a + b)}
            sage: p.coord('wz')                
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)
            sage: p.coordinates
            {'uv': (a - b, a + b),
             'wz': (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)}

        """
        manif = self.manifold
        def_chart = manif.def_chart.name
        if chart_name is None: 
            chart_name = def_chart
        if chart_name not in self.coordinates:
            # some change of coordinates must be performed
            if old_chart is None:
                # a chart must be find as a starting point of the computation
                # the manifold's default chart is privileged:
                if def_chart in self.coordinates \
                        and (def_chart, chart_name) in manif.coord_changes:
                    old_chart = def_chart
                else:
                    for ochart in self.coordinates:
                        if (ochart, chart_name) in manif.coord_changes:
                            old_chart = ochart
                            break 
            if old_chart is not None:
                chcoord = manif.coord_changes[(old_chart, chart_name)]
                self.coordinates[chart_name] = \
                                    chcoord(*self.coordinates[old_chart])
            else:
                raise ValueError("The coordinates of " + str(self) + \
                    " in the chart " + chart_name + " cannot be computed" + \
                    " by means of known changes of charts.")
        return self.coordinates[chart_name]
        
    def set_coord(self, coords, chart_name=None, delete_others=True):
        r"""
        Sets the point coordinates in the specified chart.

        INPUT:
        
        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart_name`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are
          assumed to refer to the manifold's default chart
        - ``delete_others`` -- (default: True) if True, the coordinates defined
          in other charts are deleted, otherwise their values are kept. 
          
        .. WARNING::
        
            Setting ``delete_others`` to False is at the responsability of the 
            user, who must make sure that the various sets of coordinates are
            consistent with each other. No consistency check is performed by 
            the method :meth:`set_coord`.
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: m = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart = Chart(m, 'x y z', 'cart')
            sage: p = Point(m)
            sage: p.set_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher = Chart(m, r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: q = Point(m)
            sage: q.set_coord((1,2,3), 'spher')
            sage: cart_from_spher = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            sage: q.set_coord( cart_from_spher(*q.coord('spher')), 'cart', delete_others=False)
            sage: q.coordinates
            {'spher': (1, 2, 3), 'cart': (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
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
        if chart_name is None: 
            chart_name = self.manifold.def_chart.name
        else: 
            if chart_name not in self.manifold.atlas:
                raise ValueError("The chart " + chart_name +
                    " has not been defined on the " + str(self.manifold))
        if delete_others:
            self.coordinates.clear()
        self.coordinates[chart_name] = coords

    def __eq__(self, other):
        r"""
        Compares the current point with another one.
        """
        # Search for a common chart to compare the coordinates
        common_chart = None
        # the manifold's default chart is privileged:
        def_chart = self.manifold.def_chart.name
        if def_chart in self.coordinates and def_chart in other.coordinates:
            common_chart = def_chart
        else:
            for chart in self.coordinates:
                if chart in other.coordinates:
                    common_chart = chart
                    break
        if common_chart is None:
            raise ValueError("No common chart has been found to compare " +
                             str(self) + " and " + str(other))
        return self.coordinates[common_chart] == \
                                              other.coordinates[common_chart]
        


