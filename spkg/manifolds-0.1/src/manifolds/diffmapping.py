r"""
Differentiable mappings on manifolds

The class :class:`DiffMapping` implements differentiable mappings from a 
differentiable manifold (class Manifold) to another one. 

The special case of diffeomorphisms is implemented through the class 
:class:`Diffeomorphism`, which inherits from :class:`DiffMapping`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES: 

    Defining a mapping between the sphere `S^2` and `\RR^3`::
    
        sage: m = Manifold(2, 'sphere', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'th:positive:\theta, ph:\phi', 'spher')
        sage: n = Manifold(3, 'R3', r'\mathcal{N}')
        sage: c_cart = Chart(n, 'x y z', 'cart')
        sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)) )
        
    The mapping acting on a point::
       
        sage: p = Point(m, (0,0)) # the North pole (th,ph) = (0,0)
        sage: q = Phi(p) ; q
        point on 3-dimensional manifold 'R3'
        sage: q.coordinates
        {'cart': (0, 0, 1)}
        sage: var('u v')
        (u, v)
        sage: a = Point(m, (u,v)) # a point defined by a symbolic expression
        sage: Phi(a).coord()
        (sin(u)*cos(v), sin(u)*sin(v), cos(u))
        
    An example of diffeomorphism: a rotation in some plane::
            
        sage: m = Manifold(2, "plane")
        sage: c_cart = Chart(m, 'x y', 'cart')
        sage: # A pi/3 rotation around the origin:
        sage: rot = Diffeomorphism(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2))
        sage: p = Point(m,(1,2))
        sage: q = rot(p)
        sage: q.coord()
        (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
        
    The inverse diffeormorphism::
    
        sage: irot = rot.inverse()
        sage: p1 = irot(q)
        sage: p1 == p
        True

        
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
from chart import Chart, MultiFunctionChart, CoordChange
from point import Point

     
class DiffMapping(SageObject):
    r"""
    Class for differentiable mappings between manifolds.

    INPUT:
    
    - ``manifold1`` -- start manifold 
    - ``manifold2`` -- arrival manifold 
    - ``coord_functions`` -- (default: None) the coordinate symbolic expression 
      of the mapping: list (or tuple) of the coordinates of the image expressed 
      in terms of the coordinates of the considered point; if the dimension of 
      ``manifold2`` is 1, a single expression is expected (not a list with a 
      single element)
    - ``chartname1`` -- (default: None) string defining the chart in which the 
      coordinates are given on manifold1; if none is provided, the coordinates 
      are assumed to refer to the manifold's default chart
    - ``chartname2`` -- (default: None) string defining the chart in which the 
      coordinates are given on manifold2; if none is provided, the coordinates 
      are assumed to refer to the manifold's default chart
    
    EXAMPLES:
    
    A mapping between the sphere `S^2` and `\RR^3`::
    
        sage: m = Manifold(2, 'sphere', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'th:positive:\theta, ph:\phi', 'spher')
        sage: n = Manifold(3, 'R3', r'\mathcal{N}')
        sage: c_cart = Chart(n, 'x y z', 'cart')
        sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)) )
        
    The mapping acting on a point::
       
        sage: p = Point(m, (0,0)) # the North pole (th,ph) = (0,0)
        sage: q = Phi(p) ; q
        point on 3-dimensional manifold 'R3'
        sage: q.coordinates
        {'cart': (0, 0, 1)}
        sage: var('u v')
        (u, v)
        sage: a = Point(m, (u,v)) # a point defined by a symbolic expression
        sage: Phi(a).coord()
        (sin(u)*cos(v), sin(u)*sin(v), cos(u))

    If the arrival manifold is 1-dimensional, the mapping is defined by a
    single symbolic expression and not a list with a single element::

        sage: n = Manifold(1)
        sage: chart_n = Chart(n, 'x', 'coord')
        sage: Phi = DiffMapping(m, n, sin(th)*cos(ph)) # and not ...,(sin(th)*cos(ph),))
        
    If the arrival manifold is the field of real numbers `\RR` (the Sage object
    :data:`RealLine`), the action on a point returns a real number, i.e. the 
    canonical coordinate of the image point, and not the image point itself::

        sage: Phi = DiffMapping(m, RealLine, sin(th)*cos(ph))
        sage: p = Point(m, (pi/2,0))
        sage: Phi(p)      
        1

    """
    def __init__(self, manifold1, manifold2, coord_functions=None, 
                 chartname1=None, chartname2=None): 
        if not isinstance(manifold1, Manifold):
            raise TypeError("The argument manifold1 must be a manifold.")
        if not isinstance(manifold2, Manifold):
            raise TypeError("The argument manifold2 must be a manifold.")
        self.manifold1 = manifold1
        self.manifold2 = manifold2
        if coord_functions is not None:
            if chartname1 is None: chartname1 = manifold1.def_chart.name
            if chartname2 is None: chartname2 = manifold2.def_chart.name
            if chartname1 not in self.manifold1.atlas:
                raise ValueError("The chart " + chartname1 +
                             " has not been defined on the manifold " + 
                             str(self.manifold1))
            if chartname2 not in self.manifold2.atlas:
                raise ValueError("The chart " + chartname2 +
                             " has not been defined on the manifold " + 
                             str(self.manifold2))
            chart1 = self.manifold1.atlas[chartname1]
            n2 = self.manifold2.dim
            if n2 > 1:
                if len(coord_functions) != n2:
                    raise ValueError(str(n2) + 
                                     " coordinate function must be provided.")
                self.coord_expression = {(chartname1, chartname2): 
                                        MultiFunctionChart(chart1, *coord_functions)}
            else:
                self.coord_expression = {(chartname1, chartname2): 
                                        MultiFunctionChart(chart1, coord_functions)}
        else: # case coord_functions is None:
            self.coord_expression = {}

    def new_coord_representation(self, chartname1, chartname2, coord_functions): 
        r"""
        Sets a new coordinate representation of the mapping.

        INPUT:
    
        - ``chartname1`` -- string defining the chart in which the coordinates
          are considered on the start manifold
        - ``chartname2`` -- string defining the chart in which the coordinates
          are considered on the arrival manifold
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
        
        .. WARNING::
        
            No check is performed regarding the consistency with any previously
            defined coordinate representation.
        
        EXAMPLES:
        
            Polar representation of a planar rotation initally defined in 
            Cartesian coordinates::
                
                sage: m = Manifold(2, "plane")
                sage: c_cart = Chart(m, 'x y', 'cart') 
                sage: c_spher = Chart(m, r'r:positive, ph:\phi', 'spher')
                sage: ch = CoordChange(c_cart, c_spher, sqrt(x*x+y*y), atan2(y,x))
                sage: # A pi/3 rotation around the origin defined in Cartesian coordinates:
                sage: rot = DiffMapping(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2))
                sage: rot.new_coord_representation('spher', 'spher', (r, ph+pi/3))
                sage: rot.coord_expression
                {('spher', 'spher'): functions (r, 1/3*pi + ph) on the chart 'spher' (r, ph), ('cart', 'cart'): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart 'cart' (x, y)}
                
            The rotation can then be applied by means of either coordinate 
            system::
            
                sage: p = Point(m, (1,2))  #  p defined by its Cartesian coord.
                sage: q = rot(p)  # q is computed by means of Cartesian coord.
                sage: p.change_coord('spher') # the spherical coord. of p are evaluated
                sage: q1 = rot(p, 'spher', 'spher') # q1 is computed by means of spherical coord.
                sage: q.change_coord('spher') ; # the spherical coord. of q are evaluated
                sage: q1 == q
                True
                
        """
        if chartname1 not in self.manifold1.atlas:
            raise ValueError("The chart " + chartname1 +
                  " has not been defined on the manifold " + str(self.manifold1))
        if chartname2 not in self.manifold2.atlas:
            raise ValueError("The chart " + chartname2 +
                  " has not been defined on the manifold " + str(self.manifold2))
        chart1 = self.manifold1.atlas[chartname1]
        n2 = self.manifold2.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chartname1, chartname2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chartname1, chartname2)] = \
                                   MultiFunctionChart(chart1, coord_functions)

            
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "differentiable mapping from " + str(self.manifold1) + \
                      " to " + str(self.manifold2)
        return description
        

    def __call__(self, p, chartname1=None, chartname2=None):
        r"""
        Computes the image of a point.

        INPUT:
    
        - ``p`` -- point on the manifold self.manifold1 (type: :class:`Point`)
        - ``chartname1`` -- (default: None) string defining the chart in which 
          the coordinates of p are to be considered; if none is provided, the 
          default chart of self.manifold1 will be used
        - ``chartname2`` -- (default: None) string defining the chart in which
          the coordinates of the image of p will be computed; if none is 
          provided, the default chart of self.manifold2 will be used
        
        OUTPUT:

        - image of the point by the mapping (type: :class:`Point`)

        EXAMPLES:
        
            Planar rotation acting on a point::
            
                sage: m = Manifold(2, "plane")
                sage: c_cart = Chart(m, 'x y', 'cart') 
                sage: # A pi/3 rotation around the origin defined in Cartesian coordinates:
                sage: rot = DiffMapping(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2))
                sage: p = Point(m, (1,2))  
                sage: q = rot(p)
                sage: q.coord()
                (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
                
            Image computed by means of coordinates different from the default 
            ones::
            
                sage: # Spherical coord. on the plane:
                sage: c_spher = Chart(m, r'r:positive, ph:\phi', 'spher')
                sage: ch = CoordChange(c_cart, c_spher, sqrt(x*x+y*y), atan2(y,x))
                sage: rot.new_coord_representation('spher', 'spher', (r, ph+pi/3))
                sage: p.change_coord('spher') # the spherical coord. of p are evaluated
                sage: q1 = rot(p, 'spher', 'spher') # q1 is computed by means of spherical coord.
                sage: q.change_coord('spher') ; # the spherical coord. of q are evaluated
                sage: q1 == q
                True
    
        """
        from manifold import RealLine
        if p.manifold != self.manifold1: 
            raise ValueError("The point " + str(p) +
                  " does not belong to the manifold " + str(self.manifold1))
            
        if chartname1 is None: chartname1 = self.manifold1.def_chart.name
        if chartname2 is None: chartname2 = self.manifold2.def_chart.name

        coord_map = self.coord_expression[(chartname1, chartname2)]
        x = p.coordinates[chartname1]
        y = coord_map(*x) 
        
        if self.manifold2 is RealLine:   # special case of a mapping to R
            return y[0]
        else: 
            return Point(self.manifold2, y, chartname2)

#*****************************************************************************

class Diffeomorphism(DiffMapping):
    r"""
    Class for manifold diffeomorphisms.

    INPUT:
    
    - ``manifold1`` -- start manifold 
    - ``manifold2`` -- arrival manifold 
    - ``coord_functions`` -- the coordinate symbolic expression of the mapping: 
      list (or tuple) of the coordinates of the image expressed in terms of the
      coordinates of the considered point
    - ``chartname1`` -- (default: None) string defining the chart in which the
      coordinates are given on manifold1; if none is provided, the coordinates are
      assumed to refer to the manifold's default chart
    - ``chartname2`` -- (default: None) string defining the chart in which the
      coordinates are given on manifold2; if none is provided, the coordinates are
      assumed to refer to the manifold's default chart
    
    """
    def __init__(self, manifold1, manifold2, coord_functions=None, chartname1=None, 
                 chartname2=None): 
        DiffMapping.__init__(self, manifold1, manifold2, coord_functions, 
                             chartname1, chartname2)
        if manifold1.dim != manifold2.dim:
            raise ValueError("The manifolds " + str(self.manifold1) + " and " +
                             str(self.manifold2) + 
                             " do not have the same dimension.")
        self._inverse = None
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        if self.manifold1 == self.manifold2:
            return "diffeomorphism on the " + str(self.manifold1)
        else:
            return "diffeomorphism between the " + str(self.manifold1) + \
                   " and the " + str(self.manifold2)

    def inverse(self, chartname1=None, chartname2=None): 
        r"""
        Returns the inverse diffeomorphism. 
        
        INPUT:
    
        - ``chartname1`` -- (default: None) string defining the chart in which
          the computation of the inverse is performed; if none is provided, the
          default chart of self.manifold1 will be used
        - ``chartname2`` -- (default: None) string defining the chart in which
          the computation of the inverse is performed; if none is provided, the
          default chart of self.manifold2 will be used
        
        OUTPUT:
        
        - the inverse diffeomorphism
        
        EXAMPLES:
        
            The inverse of a rotation in the plane::
            
                sage: m = Manifold(2, "plane")
                sage: c_cart = Chart(m, 'x y', 'cart')
                sage: # A pi/3 rotation around the origin:
                sage: rot = Diffeomorphism(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2))
                sage: p = Point(m,(1,2))
                sage: q = rot(p)
                sage: irot = rot.inverse()
                sage: p1 = irot(q)
                sage: p1 == p
                True
        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from utilities import simplify_chain
        if self._inverse is not None:
            return self._inverse
            
        if chartname1 is None: chartname1 = self.manifold1.def_chart.name
        if chartname2 is None: chartname2 = self.manifold2.def_chart.name

        coord_map = self.coord_expression[(chartname1, chartname2)]
        chart1 = self.manifold1.atlas[chartname1]
        chart2 = self.manifold2.atlas[chartname2]
        
        n1 = len(chart1.xx)
        n2 = len(chart2.xx)
        
        # New symbolic variables (different from chart2.xx to allow for a 
        #  correct solution even when chart2 = chart1):
        x2 = [ SR.var('xxxx' + str(i)) for i in range(n2) ]
        equations = [x2[i] == coord_map.functions[i] for i in range(n2) ]
        solutions = solve(equations, chart1.xx, solution_dict=True)
        if len(solutions) == 0: 
            raise ValueError("No solution found")
        if len(solutions) > 1: 
            raise ValueError("Non-unique solution found")
            
        #!# This should be the Python 2.7 form: 
        # substitutions = {x2[i]: chart2.xx[i] for i in range(n2)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(x2[i], chart2.xx[i]) for i in range(n2)])
       
        inv_functions = [solutions[0][chart1.xx[i]].subs(substitutions) 
                           for i in range(n1)]
        for i in range(n1):
            x = inv_functions[i]
            try:
                inv_functions[i] = simplify_chain(x)
            except AttributeError:
                pass        
        self._inverse = Diffeomorphism(self.manifold2, self.manifold1, 
                                       inv_functions, chartname2, chartname1)
        return self._inverse
