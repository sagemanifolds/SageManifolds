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
    
        sage: M = Manifold(3, 'R3', r'\mathcal{M}') 
        sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        sage: p = M.point((1, pi/2, 0), name='P') # coordinates in the manifold's default chart
        sage: p
        point 'P' on 3-dimensional manifold 'R3'
        sage: latex(p) 
        P

    Computing the coordinates of the point in a new chart::
    
        sage: c_cart.<x,y,z> = M.chart('x y z', 'cart')        
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: p.coord('cart') # evaluate P's Cartesian coordinates
        (1, 0, 0)
    
    Points can be compared::
    
        sage: p1 = M.point((1, pi/2, 0))
        sage: p == p1
        True
        sage: q = M.point((1,2,3), 'cart', name='Q') # point defined by its Cartesian coordinates
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

#from sage.structure.sage_object import SageObject
from sage.structure.element import Element   

class Point(Element):
    r"""
    Class for points on a manifold.

    INPUT:
    
    - ``domain`` -- the manifold domain to which the point belongs (can be 
      the entire manifold)
    - ``coords`` -- (default: None) the point coordinates (as a tuple or a list)
    - ``chart_name`` -- (default: None) string defining the chart in which the 
      coordinates are given; if none is provided, the coordinates are assumed 
      to refer to the domain's default chart
    - ``name`` -- (default: None) name given to the point
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the point; if 
      none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    A point on a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'M')
        sage: c_xy = M.chart('x y', 'xy')
        sage: (a, b) = var('a b') # generic coordinates for the point
        sage: p = M.point((a, b), name='P') ; p
        point 'P' on 2-dimensional manifold 'M'
        sage: p.coord()  # coordinates of P in the domain's default chart
        (a, b)
        
    A point is an element of the manifold on which it has been defined::
    
        sage: p in M
        True
        sage: p.parent()
        2-dimensional manifold 'M'
        
    By default, the LaTeX symbol of the point is deduced from its name::
    
        sage: latex(p)
        P
        
    But it can be set to any value::
    
        sage: p = M.point((a, b), name='P', latex_name=r'\mathcal{P}')
        sage: latex(p)
        \mathcal{P}
    
    """
    def __init__(self, domain, coords=None, chart_name=None, name=None, 
                 latex_name=None): 
        Element.__init__(self, domain.manifold)
        self.manifold = domain.manifold
        self.domain = domain
        self.coordinates = {}
        if coords is not None: 
            if len(coords) != self.manifold.dim: 
                raise ValueError("The number of coordinates must be equal" +
                                 " to the manifold dimension.")
            if chart_name is None: 
                chart_name = self.domain.def_chart.name
            else: 
                if chart_name not in self.domain.atlas: 
                    raise ValueError("The chart " + chart_name +
                                " has not been defined on the " + 
                                str(self.domain))
            chart = self.domain.atlas[chart_name]
            #!# The following check is not performed for it would fail with 
            # symbolic coordinates:
            # if not chart.valid_coordinates(*coords):
            #    raise ValueError("The coordinates " + str(coords) + 
            #                     " are not valid on the " +
            #                     str(self.domain.atlas[chart_name]))
            for schart in chart.supercharts:
                self.coordinates[schart.name] = tuple(coords) 
            for schart in chart.subcharts:
                if schart != chart:
                    if schart.valid_coordinates(*coords):
                        self.coordinates[schart.name] = tuple(coords)
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
          assumed to refer to the domain's default chart
        - ``old_chart`` -- (default: None) string defining the chart from 
          which the coordinates in ``chart_name`` are to be computed. If None, 
          a chart in which the point's coordinates are already known will be
          picked, priveleging the domain's default chart.

        EXAMPLES: 

        Spherical coordinates of a point on `\RR^3`::

            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher') # spherical coordinates
            sage: p = M.point((1, pi/2, 0)) 
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, 0)
            sage: p.coord('spher') # with the chart 'spher' specified (same result as above since this is the default chart)
            (1, 1/2*pi, 0)

        Computing the Cartesian coordinates from the spherical ones::

            sage: c_cart.<x,y,z> = M.chart('x y z', 'cart')  # Cartesian coordinates   
            sage: CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            coordinate change from chart 'spher' (R3, (r, th, ph)) to chart 'cart' (R3, (x, y, z))
            sage: p.coord('cart')  # the computation is performed by means of the above change of coordinates
            (1, 0, 0)

        Coordinates of a point on a 2-dimensional manifold::
    
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart('x y', 'xy')
            sage: (a, b) = var('a b') # generic coordinates for the point
            sage: p = M.point((a, b), name='P')
            sage: p.coord()  # coordinates of P in the manifold's default chart
            (a, b)
            
        Coordinates of P in a new chart::
        
            sage: c_uv.<u,v> = M.chart('u v', 'uv')
            sage: ch_xy_uv = CoordChange(c_xy, c_uv, x-y, x+y)
            sage: p.coord('uv')
            (a - b, a + b)

        Coordinates of P in a third chart::
        
            sage: c_wz.<w,z> = M.chart('w z', 'wz')
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
        atlas = self.manifold.atlas
        if chart_name is None:
            dom = self.domain 
            chart_name = dom.def_chart.name
            def_chart = chart_name
        else:
            dom = atlas[chart_name].domain
            def_chart = dom.def_chart.name
            if self not in dom:
                raise ValueError("The point does not belong to the domain " + 
                                 "of chart " + chart_name)
        if chart_name not in self.coordinates:
            # Check whether chart_name corresponds to a superchart of a chart 
            # in which the coordinates are known:
            chart = atlas[chart_name]
            for ochart in self.coordinates:
                if chart in atlas[ochart].supercharts:
                    self.coordinates[chart_name] = self.coordinates[ochart]
                    return self.coordinates[chart_name]
            # If this point is reached, some change of coordinates must be 
            # performed
            if old_chart is not None:
                s_old_chart = old_chart
            else:
                # a chart must be find as a starting point of the computation
                # the domain's default chart is privileged:
                if def_chart in self.coordinates \
                        and (def_chart, chart_name) in dom.coord_changes:
                    old_chart = def_chart
                    s_old_chart = def_chart
                else:
                    for ochart in self.coordinates:
                        for subchart in atlas[ochart].subcharts:
                            if (subchart.name, chart_name) in dom.coord_changes:
                                old_chart = ochart
                                s_old_chart = subchart.name
                                break
                        if old_chart is not None:
                            break
            if old_chart is not None:
                chcoord = dom.coord_changes[(s_old_chart, chart_name)]
                self.coordinates[chart_name] = \
                                    chcoord(*self.coordinates[old_chart])
            else:
                raise ValueError("The coordinates of " + str(self) + \
                    " in the chart " + chart_name + " cannot be computed" + \
                    " by means of known changes of charts.")
        return self.coordinates[chart_name]
        
    def set_coord(self, coords, chart_name=None):
        r"""
        Sets the point coordinates in the specified chart.
        
        Coordinates with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use the method :meth:`add_coord` 
        instead.

        INPUT:
        
        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart_name`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are
          assumed to refer to the domain's default chart
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart('x y z', 'cart')
            sage: p = M.point()
            sage: p.set_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: q = M.point()
            sage: q.set_coord((1,2,3), 'spher')
            sage: cart_from_spher = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            
        If we set the coordinates of q in the chart 'cart', those in the chart 'spher'
        are lost::
            sage: q.set_coord( cart_from_spher(*q.coord('spher')), 'cart')
            sage: q.coordinates
            {'cart': (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p.coordinates
            {'cart': (1, 2, 3)}
                        
        """
        self.coordinates.clear()
        self.add_coord(coords, chart_name)

    def add_coord(self, coords, chart_name=None):
        r"""
        Adds some coordinates in the specified chart.

        The previous coordinates with respect to other charts are kept. To
        clear them, use :meth:`set_coord` instead. 

        INPUT:
        
        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart_name`` -- (default: None) string defining the chart in which 
          the coordinates are given; if none is provided, the coordinates are
          assumed to refer to the domain's default chart
          
        .. WARNING::
        
            If the point has already coordinates in other charts, it 
            is the user's responsability to make sure that the coordinates
            to be added are consistent with them. 
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart('x y z', 'cart')
            sage: p = M.point()
            sage: p.add_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: q = M.point()
            sage: q.add_coord((1,2,3), 'spher')
            sage: cart_from_spher = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            sage: q.add_coord( cart_from_spher(*q.coord('spher')), 'cart' )
            sage: q.coordinates
            {'spher': (1, 2, 3), 'cart': (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p.coordinates
            {'cart': (1, 2, 3)}
            sage: p == q  # p and q should differ because the coordinates (1,2,3) are on different charts
            False     
            
        Contrary to :meth:`set_coord`, the method :meth:`add_coord` does not 
        the coordinates in other charts::
        
            sage: p = M.point((1,2,3), 'spher')
            sage: p.coordinates
            {'spher': (1, 2, 3)}
            sage: p.set_coord((4,5,6), 'cart')
            sage: p.coordinates
            {'cart': (4, 5, 6)}
            sage: p.add_coord((7,8,9), 'spher')
            sage: p.coordinates
            {'spher': (7, 8, 9), 'cart': (4, 5, 6)}
            
        """
        if len(coords) != self.manifold.dim: 
            raise ValueError("The number of coordinates must be equal " + 
                             "to the manifold dimension.")
        if chart_name is None: 
            chart_name = self.domain.def_chart.name
        else: 
            if chart_name not in self.domain.atlas:
                raise ValueError("The chart " + chart_name +
                    " has not been defined on the " + str(self.domain))
        self.coordinates[chart_name] = coords

    def __eq__(self, other):
        r"""
        Compares the current point with another one.
        """
        if not isinstance(other, Point):
            return False
        # Search for a common chart to compare the coordinates
        common_chart = None
        # the domain's default chart is privileged:
        def_chart = self.domain.def_chart.name
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
        


