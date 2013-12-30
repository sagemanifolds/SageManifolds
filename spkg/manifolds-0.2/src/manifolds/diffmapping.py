r"""
Differentiable mappings on manifolds

The class :class:`DiffMapping` implements differentiable mappings from a 
differentiable manifold `\mathcal{M}` to another one, `\mathcal{N}` say: 

.. MATH::

    \Phi: \mathcal{M} \longrightarrow \mathcal{N}
    
In what follows, `\mathcal{M}` is called the *start manifold* and `\mathcal{N}`
the *arrival manifold*. The case `\mathcal{N}=\mathcal{M}` is allowed. 

The special case of *diffeomorphisms*, i.e. of invertible mappings such that
both `\Phi` and `\Phi^{-1}` are differentiable, is implemented through the 
class :class:`Diffeomorphism`, which inherits from :class:`DiffMapping`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES: 

    Defining a mapping between the sphere `S^2` and `\RR^3`::
    
        sage: m = Manifold(2, 'S^2')
        sage: c_spher = Chart(m, r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher') # spherical coord. on S^2
        sage: n = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart = Chart(n, 'x y z', 'cart') # Cartesian coord. on R^3
        sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi.show()
        Phi: S^2 --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))

        
    The mapping acting on a point::
       
        sage: p = Point(m, (0,0), name='p') # the North pole (th,ph) = (0,0)
        sage: q = Phi(p) ; q
        point 'Phi(p)' on 3-dimensional manifold 'R^3'
        sage: q.coord() # Cartesian coord. of q
        (0, 0, 1)
        sage: (u, v) = var('u v')
        sage: a = Point(m, (u,v)) # a (unnamed) point defined by a symbolic expression
        sage: Phi(a)
        point on 3-dimensional manifold 'R^3'
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))
        
    An example of diffeomorphism: a rotation in the Euclidean plane::

        sage: m = Manifold(2, 'R^2', r'\RR^2')
        sage: c_cart = Chart(m, 'x y', 'cart') # Cartesian coordinates
        sage: # A pi/3 rotation around the origin:
        sage: rot = Diffeomorphism(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
        sage: p = Point(m, (1,2), name='p')
        sage: q = rot(p) ; q
        point 'R(p)' on 2-dimensional manifold 'R^2'
        sage: q.coord()
        (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
 
    The inverse diffeormorphism::

        sage: rot.inverse() 
        diffeomorphism 'R^(-1)' on the 2-dimensional manifold 'R^2'
        sage: rot.inverse().show()
        R^(-1): R^2 --> R^2, (x, y) |--> (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)
        sage: p1 = rot.inverse()(q)
        sage: p1
        point 'R^(-1)(R(p))' on 2-dimensional manifold 'R^2'
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
from chart import Chart, FunctionChart, MultiFunctionChart, CoordChange
from point import Point
from component import Components, CompWithSym, CompFullySym, CompFullyAntiSym
from tensorfield import TensorField, tensor_field_from_comp

     
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
    - ``chart1_name`` -- (default: None) string defining the chart in which the 
      coordinates are given on manifold1; if none is provided, the coordinates 
      are assumed to refer to the manifold's default chart
    - ``chart2_name`` -- (default: None) string defining the chart in which the 
      coordinates are given on manifold2; if none is provided, the coordinates 
      are assumed to refer to the manifold's default chart
    - ``name`` -- (default: None) name given to the differentiable mapping
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      differentiable mapping; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    EXAMPLES:
    
    A mapping between the sphere `S^2` and `\RR^3`::
    
        sage: m = Manifold(2, 'S^2')
        sage: c_spher = Chart(m, r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher') # spherical coord. on S^2
        sage: n = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart = Chart(n, 'x y z', 'cart') # Cartesian coord. on R^3
        sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi.show()
        Phi: S^2 --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))

    The mapping acting on a point::
       
        sage: p = Point(m, (0,0)) # the North pole (th,ph) = (0,0)
        sage: q = Phi(p) ; q
        point on 3-dimensional manifold 'R^3'
        sage: q.coordinates
        {'cart': (0, 0, 1)}
        sage: var('u v')
        (u, v)
        sage: a = Point(m, (u,v)) # a point defined by a symbolic expression
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))

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
                 chart1_name=None, chart2_name=None, name=None, 
                 latex_name=None): 
        if not isinstance(manifold1, Manifold):
            raise TypeError("The argument manifold1 must be a manifold.")
        if not isinstance(manifold2, Manifold):
            raise TypeError("The argument manifold2 must be a manifold.")
        self.manifold1 = manifold1
        self.manifold2 = manifold2
        if coord_functions is not None:
            if chart1_name is None: chart1_name = manifold1.def_chart.name
            if chart2_name is None: chart2_name = manifold2.def_chart.name
            if chart1_name not in self.manifold1.atlas:
                raise ValueError("The chart " + chart1_name +
                             " has not been defined on the manifold " + 
                             str(self.manifold1))
            if chart2_name not in self.manifold2.atlas:
                raise ValueError("The chart " + chart2_name +
                             " has not been defined on the manifold " + 
                             str(self.manifold2))
            chart1 = self.manifold1.atlas[chart1_name]
            n2 = self.manifold2.dim
            if n2 > 1:
                if len(coord_functions) != n2:
                    raise ValueError(str(n2) + 
                                     " coordinate function must be provided.")
                self.coord_expression = {(chart1_name, chart2_name): 
                                        MultiFunctionChart(chart1, *coord_functions)}
            else:
                self.coord_expression = {(chart1_name, chart2_name): 
                                        MultiFunctionChart(chart1, coord_functions)}
        else: # case coord_functions is None:
            self.coord_expression = {}
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        # Initialization of derived quantities:
        DiffMapping._init_derived(self)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "differentiable mapping"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " from " + str(self.manifold1) + " to " + \
                       str(self.manifold2)
        return description
        
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{no symbol}'
        else:
           return self.latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        pass # no derived quantity yet

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        pass # no derived quantity yet


    def show(self, chart1_name=None, chart2_name=None):
        r""" 
        Display the expression of the differentiable mapping in a given 
        pair of charts. 
        
        If the expression is not known already, it is computed from some
        expression in other charts by means of change-of-chart formulas.
        
        INPUT:
        
        - ``chart1_name`` -- (default: None) name of the chart on the start
          manifold; if none, the start manifold's default chart will be used
        - ``chart2_name`` -- (default: None) name of the chart on the arrival
          manifold; if none, the arrival manifold's default chart will be used
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Standard embedding of the sphere `S^2` in `\RR^3`::
    
            sage: m = Manifold(2, 'S^2')
            sage: c_spher = Chart(m, r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: n = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart = Chart(n, 'x y z', 'cart')
            sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: Phi.show()
            Phi: S^2 --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: latex(Phi.show())
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3 \\ & \left(\theta, \phi\right) & \longmapsto & \left(x, y, z\right) = \left(\cos\left(\phi\right) \sin\left(\theta\right), \sin\left(\phi\right) \sin\left(\theta\right), \cos\left(\theta\right)\right) \end{array}

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if chart1_name is None:
            chart1_name = self.manifold1.def_chart.name
        if chart2_name is None:
            chart2_name = self.manifold2.def_chart.name
        expression = self.expr(chart1_name, chart2_name)
        chart1 = self.manifold1.atlas[chart1_name]
        chart2 = self.manifold2.atlas[chart2_name]
        if self.name is None:
            symbol = ""
        else:
            symbol = self.name + ": "
        result.txt = symbol + self.manifold1.name + " --> " + \
                     self.manifold2.name + ", " + repr(chart1()) + " |--> " 
        if chart2 == chart1:
            result.txt += repr(expression)
        else:
            result.txt += repr(chart2()) + " = " + repr(expression)
        if self.latex_name is None:
            symbol = ""
        else:
            symbol = self.latex_name + ":"
        result.latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self.manifold1) + r"& \longrightarrow & " + \
                       latex(self.manifold2) + r"\\ &" + latex(chart1()) + \
                       r"& \longmapsto & " 
        if chart2 == chart1:
            result.latex += latex(expression) + r"\end{array}"
        else:
            result.latex += latex(chart2()) + " = " + latex(expression) + \
                            r"\end{array}"
        return result


    def multi_function_chart(self, chart1_name=None, chart2_name=None):
        r""" 
        Return the functions of the coordinates representing the differentiable
        mapping in a given pair of charts.
        
        If these functions are not already known, they are computed from known 
        ones by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1_name`` -- (default: None) name of the chart on the start
          manifold; if None, the start manifold's default chart is assumed
        - ``chart2_name`` -- (default: None) name of the chart on the arrival
          manifold; if None, the arrival manifold's default chart is assumed

        OUTPUT:
        
        - instance of :class:`MultiFunctionChart` representing the 
          differentiable mapping in the above two charts

        EXAMPLES:

        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
        
            sage: m = Manifold(2, 'M')
            sage: n = Manifold(3, 'N')
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: c_xyz = Chart(n, 'x y z', 'xyz')
            sage: Phi = DiffMapping(m, n, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.show()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.multi_function_chart('uv', 'xyz')
            functions (u*v, u/v, u + v) on the chart 'uv' (u, v)
            sage: Phi.multi_function_chart() # equivalent to above since 'uv' and 'xyz' are default charts
            functions (u*v, u/v, u + v) on the chart 'uv' (u, v)
            sage: type(Phi.multi_function_chart())
            <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>

        Representation in other charts::
        
            sage: c_UV = Chart(m, 'U V', 'UV')  # new chart on M
            sage: ch_uv_UV = CoordChange(c_uv, c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ = Chart(n, 'X Y Z', 'XYZ') # new chart on N
            sage: ch_xyz_XYZ = CoordChange(c_xyz, c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.multi_function_chart('UV', 'xyz')
            functions (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V) on the chart 'UV' (U, V)
            sage: Phi.multi_function_chart('uv', 'XYZ')
            functions (((2*u + 1)*v^2 + u*v - 3*u)/v, -((u - 1)*v^2 - u*v - u)/v, -((u + 1)*v^2 + u*v - 2*u)/v) on the chart 'uv' (u, v)
            sage: Phi.multi_function_chart('UV', 'XYZ')
            functions (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V)) on the chart 'UV' (U, V)

        """
        manif1 = self.manifold1; manif2 = self.manifold2
        def_chart1 = manif1.def_chart.name; def_chart2 = manif2.def_chart.name
        if chart1_name is None:
            chart1_name = def_chart1
        if chart2_name is None:
            chart2_name = def_chart2
        if (chart1_name, chart2_name) not in self.coord_expression:
            # some change of coordinates must be performed
            chart1 = manif1.atlas[chart1_name]
            change_start = [] ; change_arrival = []
            for (ochart1, ochart2) in self.coord_expression:
                if chart1_name == ochart1:
                    change_arrival.append(ochart2)
                if chart2_name == ochart2:
                    change_start.append(ochart1)
            # 1/ Trying to make a change of chart only on the arrival manifold:
            # the arrival default chart is privileged:
            sel_chart2 = None # selected chart2
            if def_chart2 in change_arrival \
                    and (def_chart2, chart2_name) in manif2.coord_changes:
                sel_chart2 = def_chart2
            else:
                for ochart2 in change_arrival:
                    if (ochart2, chart2_name) in manif2.coord_changes:
                        sel_chart2 = ochart2
                        break 
            if sel_chart2 is not None:
                oexpr = self.coord_expression[(chart1_name, sel_chart2)]
                chg2 = manif2.coord_changes[(sel_chart2, chart2_name)]
                self.coord_expression[(chart1_name, chart2_name)] = \
                    MultiFunctionChart(chart1, *(chg2(*(oexpr.expr()))) )
                return self.coord_expression[(chart1_name, chart2_name)]

            # 2/ Trying to make a change of chart only on the start manifold:
            # the start default chart is privileged:
            sel_chart1 = None # selected chart1
            if def_chart1 in change_start \
                    and (chart1_name, def_chart1) in manif1.coord_changes:
                sel_chart1 = def_chart1
            else:
                for ochart1 in change_start:
                    if (chart1_name, ochart1) in manif1.coord_changes:
                        sel_chart1 = ochart1
                        break
            if sel_chart1 is not None:
                oexpr = self.coord_expression[(sel_chart1, chart2_name)]
                chg1 = manif1.coord_changes[(chart1_name, sel_chart1)]
                self.coord_expression[(chart1_name, chart2_name)] = \
                    MultiFunctionChart(chart1, 
                                       *(oexpr( *(chg1.transf.expr()) )) )
                return self.coord_expression[(chart1_name, chart2_name)]
                    
            # 3/ If this point is reached, it is necessary to perform some 
            # coordinate change both on the start manifold and the arrival one
            # the default charts are privileged:
            if (def_chart1, def_chart2) in self.coord_expression \
                    and (chart1_name, def_chart1) in manif1.coord_changes \
                    and (def_chart2, chart2_name) in manif2.coord_changes:
                sel_chart1 = def_chart1
                sel_chart2 = def_chart2
            else:
                for (ochart1, ochart2) in self.coord_expression:
                    if (chart1_name, ochart1) in manif1.coord_changes \
                        and (ochart2, chart2_name) in manif2.coord_changes:
                        sel_chart1 = ochart1
                        sel_chart2 = ochart2
                        break
            if (sel_chart1 is not None) and (sel_chart2 is not None):
                oexpr = self.coord_expression[(sel_chart1, sel_chart2)]
                chg1 = manif1.coord_changes[(chart1_name, sel_chart1)]
                chg2 = manif2.coord_changes[(sel_chart2, chart2_name)]
                self.coord_expression[(chart1_name, chart2_name)] = \
                     MultiFunctionChart(chart1, 
                                *(chg2( *(oexpr(*(chg1.transf.expr()))) )) )
                return self.coord_expression[(chart1_name, chart2_name)]
                
            # 4/ If this point is reached, the demanded value cannot be
            # computed 
            raise ValueError("The expression of the mapping in the pair of " +
                "charts (" + chart1_name + ", " + chart2_name + ") cannot " + 
                "be computed by means of known changes of charts.")
                
        return self.coord_expression[(chart1_name, chart2_name)]
            

    def expr(self, chart1_name=None, chart2_name=None):
        r""" 
        Return the expression of the differentiable mapping in terms of
        specified coordinates.
        
        If the expression is not already known, it is computed from some known 
        expression by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1_name`` -- (default: None) name of the chart on the start
          manifold; if None, the start manifold's default chart is assumed
        - ``chart2_name`` -- (default: None) name of the chart on the arrival
          manifold; if None, the arrival manifold's default chart is assumed

        OUTPUT:
        
        - symbolic expression representing the differentiable mapping in the 
          above two charts

        EXAMPLES:
        
        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
        
            sage: m = Manifold(2, 'M')
            sage: n = Manifold(3, 'N')
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: c_xyz = Chart(n, 'x y z', 'xyz')
            sage: Phi = DiffMapping(m, n, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.show()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.expr('uv', 'xyz')
            (u*v, u/v, u + v)
            sage: Phi.expr()  # equivalent to above since 'uv' and 'xyz' are default charts
            (u*v, u/v, u + v)
            sage: type(Phi.expr()[0])
            <type 'sage.symbolic.expression.Expression'>

        Expressions in other charts::
        
            sage: c_UV = Chart(m, 'U V', 'UV')  # new chart on M
            sage: ch_uv_UV = CoordChange(c_uv, c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ = Chart(n, 'X Y Z', 'XYZ') # new chart on N
            sage: ch_xyz_XYZ = CoordChange(c_xyz, c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.expr('UV', 'xyz')
            (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V)
            sage: Phi.expr('uv', 'XYZ')
            (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v,
             -((u + 1)*v^2 + u*v - 2*u)/v)
            sage: Phi.expr('UV', 'XYZ')
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))

        A rotation in some Euclidean plane::
        
            sage: m = Manifold(2, 'M') # the plane
            sage: c_spher = Chart(m, r'r:[0,+oo) ph:[0,2*pi):\phi', 'spher') # spherical coordinates on the plane
            sage: rot = DiffMapping(m, m, (r, ph+pi/3), name='R') # pi/3 rotation around r=0
            sage: rot.expr()
            (r, 1/3*pi + ph)

        Expression of the rotation in terms of Cartesian coordinates::
        
            sage: c_cart = Chart(m, 'x y', 'cart') # Declaration of Cartesian coordinates
            sage: ch_spher_cart = CoordChange(c_spher, c_cart, r*cos(ph), r*sin(ph)) # relation to spherical coordinates
            sage: ch_spher_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))              
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: rot.expr('cart', 'cart')                            
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        """
        return self.multi_function_chart(chart1_name, chart2_name).expr()
            
    def set_expr(self, chart1_name, chart2_name, coord_functions, 
                 delete_others=True): 
        r"""
        Set a new coordinate representation of the mapping.

        INPUT:
    
        - ``chart1_name`` -- string defining the chart for the coordinates
          on the start manifold
        - ``chart2_name`` -- string defining the chart for the coordinates
          on the arrival manifold
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
        - ``delete_others`` -- (default: True) determines whether the 
          coordinate expressions with respect to charts different from that 
          specified by ``chart1_name`` and ``chart2_name`` are deleted or not. 
          
        .. WARNING::
        
            Setting ``delete_others`` to False is at the responsability of the 
            user, who must make sure that the various expressions are
            consistent with each other. 
        
        
        EXAMPLES:
        
        Polar representation of a planar rotation initally defined in 
        Cartesian coordinates::
            
            sage: m = Manifold(2, 'R^2', r'\RR^2')   # Euclidean plane
            sage: c_cart = Chart(m, 'x y', 'cart')  # Cartesian coordinates
            sage: c_spher = Chart(m, r'r:[0,+oo) ph:[0,2*pi):\phi', 'spher') # spherical coordinates
            sage: # Links between spherical coordinates and Cartesian ones:
            sage: ch_cart_spher = CoordChange(c_cart, c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
               x == x
               y == y
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
            sage: # A pi/3 rotation around the origin defined in terms of Cartesian coordinates:
            sage: rot = DiffMapping(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.expr()
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        If we try to make Sage calculate the expression in terms of spherical
        coordinates, via the method :meth:`expr`, we notice some difficulties
        in arctan2 simplifications::
        
            sage: rot.expr('spher', 'spher') # correct output but could be simplified !
            (r,
             arctan2(1/2*(sqrt(3)*cos(ph) + sin(ph))*r, -1/2*(sqrt(3)*sin(ph) - cos(ph))*r))
        
        Therefore, we use the method :meth:`set_expr` to set the 
        spherical-coordinate expression by hand::

            sage: rot.set_expr('spher', 'spher', (r, ph+pi/3), delete_others=False)
            sage: rot.expr('spher', 'spher')  # the output is now satisfactory
            (r, 1/3*pi + ph)
        
        Because ``delete_others`` was set to False, the expression in Cartesian
        coordinates has been kept in the dictionary :attr:`coord_expression` 
        that stores the various representations of the differentiable mapping::
             
            sage: rot.coord_expression
            {('spher', 'spher'): functions (r, 1/3*pi + ph) on the chart 'spher' (r, ph), ('cart', 'cart'): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart 'cart' (x, y)}
             
        If, on the contrary, we use :meth:`set_expr` with ``delete_others`` set
        to True (the default), the expression in Cartesian coordinates is 
        lost::
        
            sage: rot.set_expr('spher', 'spher', (r, ph+pi/3))
            sage: rot.coord_expression                        
            {('spher', 'spher'): functions (r, 1/3*pi + ph) on the chart 'spher' (r, ph)}
 
        It is recovered by a call to :meth:`expr`::
        
            sage: rot.expr('cart','cart')
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot.coord_expression
            {('spher', 'spher'): functions (r, 1/3*pi + ph) on the chart 'spher' (r, ph), ('cart', 'cart'): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart 'cart' (x, y)}
            
        The rotation can be applied to a point by means of either coordinate 
        system::
            
            sage: p = Point(m, (1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p.coord('spher') # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, 'spher', 'spher') # q1 is computed by means of spherical coord.
            sage: q.coord('spher') ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
                
        """
        if chart1_name not in self.manifold1.atlas:
            raise ValueError("The chart " + chart1_name +
               " has not been defined on the " + str(self.manifold1))
        if chart2_name not in self.manifold2.atlas:
            raise ValueError("The chart " + chart2_name +
              " has not been defined on the " + str(self.manifold2))
        chart1 = self.manifold1.atlas[chart1_name]
        if delete_others:
            self.coord_expression.clear()
            self._del_derived()
        n2 = self.manifold2.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chart1_name, chart2_name)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chart1_name, chart2_name)] = \
                                   MultiFunctionChart(chart1, coord_functions)


    def __call__(self, p, chart1_name=None, chart2_name=None):
        r"""
        Compute the image of a point.

        INPUT:
    
        - ``p`` -- point on the start manifold (type: :class:`Point`)
        - ``chart1_name`` -- (default: None) string defining the chart in which 
          the coordinates of p are to be considered; if none is provided, a
          chart in which both p's coordinates and the expression of ``self``
          are known is searched, starting from the 
          default chart of self.manifold1 will be used
        - ``chart2_name`` -- (default: None) string defining the chart in which
          the coordinates of the image of p will be computed; if none is 
          provided, the default chart of self.manifold2 is assumed.
        
        OUTPUT:

        - image of the point by the mapping (type: :class:`Point`)

        EXAMPLES:
        
            Planar rotation acting on a point::
            
                sage: m = Manifold(2, 'R^2', r'\RR^2') # Euclidean plane
                sage: c_cart = Chart(m, 'x y', 'cart') # Cartesian coordinates
                sage: # A pi/3 rotation around the origin defined in Cartesian coordinates:
                sage: rot = DiffMapping(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
                sage: p = Point(m, (1,2), name='p')
                sage: q = rot(p) ; q
                point 'R(p)' on 2-dimensional manifold 'R^2'
                sage: q.coord()
                (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
                
            Image computed by means of coordinates different from the default 
            ones::
            
                sage: # Spherical coord. on the plane:
                sage: c_spher = Chart(m, r'r:[0,+oo) ph:[0,2*pi):\phi', 'spher')
                sage: ch = CoordChange(c_cart, c_spher, sqrt(x*x+y*y), atan2(y,x))
                sage: rot.set_expr('spher', 'spher', (r, ph+pi/3), delete_others=False)
                sage: p.coord('spher') # the spherical coord. of p are evaluated
                (sqrt(5), arctan(2))
                sage: q1 = rot(p, 'spher', 'spher') # q1 is computed by means of spherical coord.
                sage: q.coord('spher') ; # the spherical coord. of q are evaluated
                (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
                sage: q1 == q
                True
    
        """
        from manifold import RealLine
        if p.manifold != self.manifold1: 
            raise ValueError("The point " + str(p) +
                  " does not belong to the manifold " + str(self.manifold1))
        if chart2_name is None: 
            chart2_name = self.manifold2.def_chart.name
        if chart1_name is None: 
            def_chart1 = self.manifold1.def_chart.name
            if def_chart1 in p.coordinates and \
                        (def_chart1, chart2_name) in self.coord_expression:
                chart1_name = def_chart1
            else:
                for chart in p.coordinates:
                    if (chart, chart2_name) in self.coord_expression:
                        chart1_name = chart
                        break
        if chart1_name is None:
            raise ValueError("No common chart has been found to evaluate " \
                "the action of " + str(self) + " on the " + str(p) + ".")

        coord_map = self.coord_expression[(chart1_name, chart2_name)]
        y = coord_map(*(p.coordinates[chart1_name])) 
        
        if self.manifold2 is RealLine:   # special case of a mapping to R
            return y[0]
        else:
            if p.name is None or self.name is None:
                res_name = None
            else:
                res_name = self.name + '(' + p.name + ')'
            if p.latex_name is None or self.latex_name is None:
                res_latex_name = None
            else:
                res_latex_name = self.latex_name + r'\left(' + p.latex_name + \
                                 r'\right)'
            
            return Point(self.manifold2, y, chart2_name, name=res_name, 
                         latex_name=res_latex_name)

    def pullback(self, tensor):
        r""" 
        Pullback operator associated with the differentiable mapping. 
        
        INPUT:
        
        - ``tensor`` -- instance of :class:`TensorField` representing a fully 
          covariant tensor field `T` on the *arrival* manifold, i.e. a tensor 
          field of type (0,p), with p a positive or zero integer. The case p=0 
          corresponds to a scalar field.
          
        OUTPUT:
        
        - instance of :class:`TensorField` representing a fully 
          covariant tensor field on the *start* manifold that is the 
          pullback of `T` given by ``self``. 
          
        EXAMPLES:
        
        Pullback on `S^2` of a scalar field defined on `R^3`::
        
            sage: m = Manifold(2, 'S^2', start_index=1)
            sage: c_spher = Chart(m, r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher') # spherical coord. on S^2
            sage: n = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart = Chart(n, 'x y z', 'cart') # Cartesian coord. on R^3
            sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: f = ScalarField(n, x*y*z, name='f') ; f
            scalar field 'f' on the 3-dimensional manifold 'R^3'
            sage: f.show()
            f: (x, y, z) |--> x*y*z
            sage: pf = Phi.pullback(f) ; pf
            scalar field 'Phi_*(f)' on the 2-dimensional manifold 'S^2'
            sage: pf.show()
            Phi_*(f): (th, ph) |--> cos(ph)*cos(th)*sin(ph)*sin(th)^2
            
        Pullback on `S^2` of the standard Euclidean metric on `R^3`::
                
            sage: g = SymBilinFormField(n, 'g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.show()
            g = dx*dx + dy*dy + dz*dz
            sage: pg = Phi.pullback(g) ; pg
            field of symmetric bilinear forms 'Phi_*(g)' on the 2-dimensional manifold 'S^2'
            sage: pg.show()
            Phi_*(g) = dth*dth + sin(th)^2 dph*dph

        Pullback on `S^2` of a 3-form on `R^3`::
                
            sage: a = DiffForm(n, 3, 'A')
            sage: a[1,2,3] = f 
            sage: a.show()
            A = x*y*z dx/\dy/\dz
            sage: pa = Phi.pullback(a) ; pa
            3-form 'Phi_*(A)' on the 2-dimensional manifold 'S^2'
            sage: pa.show() # should be zero (as any 3-form on a 2-dimensional manifold)
            Phi_*(A) = 0

        """
        from scalarfield import ScalarField
        if not isinstance(tensor, TensorField):
            raise TypeError("The argument 'tensor' must be a tensor field.")
        manif1 = self.manifold1
        manif2 = self.manifold2
        if tensor.manifold != manif2:
            raise TypeError("The tensor field is not defined on the mapping " +
                            "arrival manifold.")
        (ncon, ncov) = tensor.tensor_type
        if ncon != 0:
            raise TypeError("The pullback cannot be taken on a tensor " + 
                            "with some contravariant part.")
        resu_name = None ; resu_latex_name = None
        if self.name is not None and tensor.name is not None:
            resu_name = self.name + '_*(' + tensor.name + ')'
        if self.latex_name is not None and tensor.latex_name is not None:
            resu_latex_name = self.latex_name + '_*' + tensor.latex_name                
        chart1_name = None; chart2_name = None
        def_chart1 = manif1.def_chart ; def_chart1_name = def_chart1.name
        def_chart2 = manif2.def_chart ; def_chart2_name = def_chart2.name
        if ncov == 0:
            # Case of a scalar field
            # ----------------------
            # A pair of chart (chart1_name, chart2_name) where the computation
            # is feasable is searched, privileging the default chart of the 
            # start manifold for chart1_name
            if def_chart2_name in tensor.express and \
                   (def_chart1_name, def_chart2_name) in self.coord_expression:
                chart1_name = def_chart1_name
                chart2_name = def_chart2_name
            else:
                for (chart1n, chart2n) in self.coord_expression:
                    if chart1n == def_chart1_name and chart2n in tensor.express:
                        chart1_name = def_chart1_name
                        chart2_name = chart2n
                        break
            if chart1_name is None:
                # It is not possible to have def_chart1 as chart for 
                # expressing the result; any other chart is then looked for:
                for (chart1n, chart2n) in self.coord_expression:
                    if chart2n in tensor.express:
                        chart1_name = chart1n
                        chart2_name = chart2n
                        break
            if chart1_name is None:
                raise ValueError("No common chart could be find to compute " +
                    "the pullback of the scalar field.")
            phi = self.coord_expression[(chart1_name, chart2_name)]
            coord1 = manif1.atlas[chart1_name].xx
            ff = tensor.express[chart2_name]
            return ScalarField(manif1, ff(*(phi(*coord1))), chart1_name, 
                               resu_name, resu_latex_name)
        else:
            # Case of tensor field of rank >= 1
            # ---------------------------------
            # A pair of chart (chart1_name, chart2_name) where the computation
            # is feasable is searched, privileging the default chart of the 
            # start manifold for chart1_name
            if def_chart2.frame.name in tensor.components and \
                   (def_chart1_name, def_chart2_name) in self.coord_expression:
                chart1_name = def_chart1_name
                chart2_name = def_chart2_name
            else:
                for (chart1n, chart2n) in self.coord_expression:
                    if chart1n == def_chart1_name and \
                         manif2.atlas[chart2n].frame.name in tensor.components:
                        chart1_name = def_chart1_name
                        chart2_name = chart2n
                        break
            if chart1_name is None:
                # It is not possible to have def_chart1 as chart for 
                # expressing the result; any other chart is then looked for:
                for (chart1n, chart2n) in self.coord_expression:
                    if manif2.atlas[chart2n].frame.name in tensor.components:
                        chart1_name = chart1n
                        chart2_name = chart2n
                        break
            if chart1_name is None:
                raise ValueError("No common chart could be find to compute " +
                    "the pullback of the tensor field.")
            chart1 = manif1.atlas[chart1_name]
            chart2 = manif2.atlas[chart2_name]
            frame1_name = chart1.frame.name
            frame2_name = chart2.frame.name
            # Computation at the component level:
            tcomp = tensor.components[frame2_name]
            if isinstance(tcomp, CompFullySym):
                ptcomp = CompFullySym(manif1, ncov, frame1_name)
            elif isinstance(tcomp, CompFullyAntiSym):
                ptcomp = CompFullyAntiSym(manif1, ncov, frame1_name)
            elif isinstance(tcomp, CompWithSym):
                ptcomp = CompWithSym(manif1, ncov, frame1_name, sym=tcomp.sym, 
                                  antisym=tcomp.antisym)
            else:
                ptcomp = Components(manif1, ncov, frame1_name)
            phi = self.coord_expression[(chart1_name, chart2_name)]
            jacob = phi.jacobian()
            # X2 coordinates expressed in terms of X1 ones via the mapping:
            coord2_1 = phi(*(chart1.xx)) 
            si1 = manif1.sindex
            si2 = manif2.sindex
            for ind_new in ptcomp.non_redundant_index_generator(): 
                res = 0 
                for ind_old in manif2.index_generator(ncov): 
                    ff = tcomp[[ind_old]].function_chart(chart2_name)
                    t = FunctionChart(chart1, ff(*coord2_1))
                    for i in range(ncov):
                        t *= jacob[ind_old[i]-si2][ind_new[i]-si1]
                    res += t
                ptcomp[ind_new] = res
            resu = tensor_field_from_comp(ptcomp, (0, ncov))
            resu.set_name(resu_name, resu_latex_name)
            return resu

        
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
    - ``chart1_name`` -- (default: None) string defining the chart in which the
      coordinates are given on manifold1; if none is provided, the coordinates are
      assumed to refer to the manifold's default chart
    - ``chart2_name`` -- (default: None) string defining the chart in which the
      coordinates are given on manifold2; if none is provided, the coordinates are
      assumed to refer to the manifold's default chart
    - ``name`` -- (default: None) name given to the differentiable mapping
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      differentiable mapping; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    """
    def __init__(self, manifold1, manifold2, coord_functions=None, 
                 chart1_name=None, chart2_name=None, name=None, 
                 latex_name=None): 
        DiffMapping.__init__(self, manifold1, manifold2, coord_functions, 
                             chart1_name, chart2_name, name, latex_name)
        if manifold1.dim != manifold2.dim:
            raise ValueError("The manifolds " + str(self.manifold1) + " and " +
                             str(self.manifold2) + 
                             " do not have the same dimension.")
        # Initialization of derived quantities:
        Diffeomorphism._init_derived(self)
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "diffeomorphism"
        if self.name is not None:
            description += " '%s'" % self.name
        if self.manifold1 == self.manifold2:
            description += " on the " + str(self.manifold1)
        else:
            description += " between the " + str(self.manifold1) + \
                           " and the " + str(self.manifold2)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        DiffMapping._init_derived(self) # derived quantities of the mother class
        self._inverse = None

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        DiffMapping._del_derived(self) # derived quantities of the mother class
        self._inverse = None


    def inverse(self, chart1_name=None, chart2_name=None): 
        r"""
        Return the inverse diffeomorphism. 
        
        INPUT:
    
        - ``chart1_name`` -- (default: None) string defining the chart in which
          the computation of the inverse is performed if necessary; if none 
          is provided, the default chart of the start manifold will be used
        - ``chart2_name`` -- (default: None) string defining the chart in which
          the computation of the inverse is performed if necessary; if none 
          is provided, the default chart of the arrival manifold will be used
        
        OUTPUT:
        
        - the inverse diffeomorphism
        
        EXAMPLES:
        
            The inverse of a rotation in the Euclidean plane::
            
                sage: m = Manifold(2, 'R^2', r'\RR^2')
                sage: c_cart = Chart(m, 'x y', 'cart')
                sage: # A pi/3 rotation around the origin:
                sage: rot = Diffeomorphism(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
                sage: rot.inverse() 
                diffeomorphism 'R^(-1)' on the 2-dimensional manifold 'R^2'
                sage: rot.inverse().show()
                R^(-1): R^2 --> R^2, (x, y) |--> (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)

            Checking that applying successively the diffeomorphism and its 
            inverse results in the identity::
            
                sage: (a, b) = var('a b')
                sage: p = Point(m, (a,b)) # a generic point on M
                sage: q = rot(p)
                sage: p1 = rot.inverse()(q)
                sage: p1 == p 
                True

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from utilities import simplify_chain
        if self._inverse is not None:
            return self._inverse
            
        if chart1_name is None: chart1_name = self.manifold1.def_chart.name
        if chart2_name is None: chart2_name = self.manifold2.def_chart.name

        coord_map = self.coord_expression[(chart1_name, chart2_name)]
        chart1 = self.manifold1.atlas[chart1_name]
        chart2 = self.manifold2.atlas[chart2_name]
        
        n1 = len(chart1.xx)
        n2 = len(chart2.xx)
        
        # New symbolic variables (different from chart2.xx to allow for a 
        #  correct solution even when chart2 = chart1):
        x2 = [ SR.var('xxxx' + str(i)) for i in range(n2) ]
        equations = [ x2[i] == coord_map.functions[i].express 
                      for i in range(n2) ]
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
        if self.name is None:
            name = None
        else:
            name = self.name + '^(-1)'
        
        if self.latex_name is None:
            latex_name = None
        else:
            latex_name = self.latex_name + r'^{-1}'
        self._inverse = Diffeomorphism(self.manifold2, self.manifold1, 
                                       inv_functions, chart2_name, chart1_name,
                                       name=name, latex_name=latex_name)
        return self._inverse
