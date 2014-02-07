r"""
Differentiable mappings between manifolds

The class :class:`DiffMapping` implements differentiable mappings from an open
domain `U` of a differentiable manifold `\mathcal{M}` to a differentiable
manifold `\mathcal{N}`: 

.. MATH::

    \Phi: U\subset \mathcal{M} \longrightarrow \mathcal{N}
    
In what follows, `\mathcal{M}` is called the *start manifold* and `\mathcal{N}`
the *arrival manifold*. The case `\mathcal{N}=\mathcal{M}` is allowed. 

The special case of *diffeomorphisms*, i.e. of invertible mappings such that
both `\Phi` and `\Phi^{-1}` are differentiable, is implemented through the 
class :class:`Diffeomorphism`, which inherits from :class:`DiffMapping`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES: 

    A mapping between the sphere `S^2` and `\RR^3`::

        sage: m = Manifold(2, 'S^2')
        sage: U = m.open_domain('U') # the subdomain of S^2 covered by regular spherical coordinates
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: n = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart.<x,y,z> = n.chart('x y z')  # Cartesian coord. on R^3
        sage: Phi = DiffMapping(U, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi.view()
        Phi: U --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
        
    The mapping acting on a point::
       
        sage: p = U.point((pi/2,pi/2), name='P')
        sage: q = Phi(p) ; q
        point 'Phi(P)' on 3-dimensional manifold 'R^3'
        sage: q.coord()
        (0, 1, 0)
        sage: (u, v) = var('u v')
        sage: a = U.point((u,v)) # a (unnamed) point defined by a symbolic expression
        sage: Phi(a)
        point on 3-dimensional manifold 'R^3'
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))
        
    An example of diffeomorphism: a rotation in the Euclidean plane::

        sage: m = Manifold(2, 'R^2', r'\RR^2')
        sage: c_cart.<x,y> = m.chart('x y') # Cartesian coordinates
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
        sage: rot.inverse().view()
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
from domain import Domain
from chart import Chart, FunctionChart, MultiFunctionChart, CoordChange
from point import Point
from component import Components, CompWithSym, CompFullySym, CompFullyAntiSym
from tensorfield import TensorField

     
class DiffMapping(SageObject):
    r"""
    Class for differentiable mappings between manifolds.

    INPUT:
    
    - ``domain1`` -- domain on the start manifold 
    - ``domain2`` -- domain on the arrival manifold 
    - ``coord_functions`` -- (default: None) the coordinate symbolic expression 
      of the mapping: list (or tuple) of the coordinates of the image expressed 
      in terms of the coordinates of the considered point; if the dimension of 
      the arrival manifold is 1, a single expression is expected 
      (not a list with a single element)
    - ``chart1`` -- (default: None) hart in which the 
      coordinates are given on domain1; if none is provided, the coordinates 
      are assumed to refer to domain's default chart
    - ``chart2`` -- (default: None) chart in which the 
      coordinates are given on domain2; if none is provided, the coordinates 
      are assumed to refer to the domain's default chart
    - ``name`` -- (default: None) name given to the differentiable mapping
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      differentiable mapping; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    EXAMPLES:
    
    A mapping between the sphere `S^2` and `\RR^3`::

        sage: m = Manifold(2, 'S^2')
        sage: U = m.open_domain('U') # the subdomain of S^2 covered by regular spherical coordinates
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: n = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart.<x,y,z> = n.chart('x y z')  # Cartesian coord. on R^3
        sage: Phi = DiffMapping(U, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi.view()
        Phi: U --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
        
    The mapping acting on a point::
       
        sage: p = U.point((pi/2,pi/2), name='P')
        sage: q = Phi(p) ; q
        point 'Phi(P)' on 3-dimensional manifold 'R^3'
        sage: q.coord()
        (0, 1, 0)
        sage: (u, v) = var('u v')
        sage: a = U.point((u,v)) # a point defined by a symbolic expression
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))

    If the arrival manifold is 1-dimensional, the mapping is defined by a
    single symbolic expression and not a list with a single element::

        sage: n = Manifold(1, 'N')
        sage: chart_n = n.chart('x')
        sage: Phi = DiffMapping(m, n, sin(th)*cos(ph)) # and not ...,(sin(th)*cos(ph),))
        
    If the arrival manifold is the field of real numbers `\RR` (the Sage object
    :data:`RealLine`), the action on a point returns a real number, i.e. the 
    canonical coordinate of the image point, and not the image point itself::

        sage: Phi = DiffMapping(m, RealLine, sin(th)*cos(ph))
        sage: p = U.point((pi/2,pi))
        sage: Phi(p)      
        -1

    """
    def __init__(self, domain1, domain2, coord_functions=None, chart1=None, 
                 chart2=None, name=None, latex_name=None): 
        if not isinstance(domain1, Domain):
            raise TypeError("The argument domain1 must be a domain.")
        if not isinstance(domain2, Domain):
            raise TypeError("The argument domain2 must be a domain.")
        self.domain1 = domain1
        self.domain2 = domain2
        if coord_functions is not None:
            if chart1 is None: chart1 = domain1.def_chart
            if chart2 is None: chart2 = domain2.def_chart
            if chart1 not in self.domain1.atlas:
                raise ValueError("The " + str(chart1) +
                                    " has not been defined on the " + 
                                    str(self.domain1))
            if chart2 not in self.domain2.atlas:
                raise ValueError("The " + str(chart2) +
                                    " has not been defined on the " + 
                                    str(self.domain2))
            n2 = self.domain2.manifold.dim
            if n2 > 1:
                if len(coord_functions) != n2:
                    raise ValueError(str(n2) + 
                                     " coordinate function must be provided.")
                self.coord_expression = {(chart1, chart2): 
                                        MultiFunctionChart(chart1, *coord_functions)}
            else:
                self.coord_expression = {(chart1, chart2): 
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
        description += " from " + str(self.domain1) + " to " + \
                       str(self.domain2)
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


    def view(self, chart1=None, chart2=None):
        r""" 
        Display the expression of the differentiable mapping in a given 
        pair of charts. 
        
        If the expression is not known already, it is computed from some
        expression in other charts by means of change-of-chart formulas.
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the start domain; if None, 
          the start domain's default chart will be used
        - ``chart2`` -- (default: None) chart on the arrival domain; if None,
          the arrival domain's default chart will be used
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Standard embedding of the sphere `S^2` in `\RR^3`::
    
            sage: m = Manifold(2, 'S^2')
            sage: c_spher.<th,ph> = m.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: n = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = n.chart('x y z')
            sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: S^2 --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: latex(Phi.view())
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3 \\ & \left(\theta, \phi\right) & \longmapsto & \left(x, y, z\right) = \left(\cos\left(\phi\right) \sin\left(\theta\right), \sin\left(\phi\right) \sin\left(\theta\right), \cos\left(\theta\right)\right) \end{array}

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if chart1 is None:
            chart1 = self.domain1.def_chart
        if chart2 is None:
            chart2 = self.domain2.def_chart
        expression = self.expr(chart1, chart2)
        if self.name is None:
            symbol = ""
        else:
            symbol = self.name + ": "
        result.txt = symbol + self.domain1.name + " --> " + \
                     self.domain2.name + ", " + repr(chart1[:]) + " |--> " 
        if chart2 == chart1:
            result.txt += repr(expression)
        else:
            result.txt += repr(chart2[:]) + " = " + repr(expression)
        if self.latex_name is None:
            symbol = ""
        else:
            symbol = self.latex_name + ":"
        result.latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self.domain1) + r"& \longrightarrow & " + \
                       latex(self.domain2) + r"\\ &" + latex(chart1[:]) + \
                       r"& \longmapsto & " 
        if chart2 == chart1:
            result.latex += latex(expression) + r"\end{array}"
        else:
            result.latex += latex(chart2[:]) + " = " + latex(expression) + \
                            r"\end{array}"
        return result


    def multi_function_chart(self, chart1=None, chart2=None):
        r""" 
        Return the functions of the coordinates representing the differentiable
        mapping in a given pair of charts.
        
        If these functions are not already known, they are computed from known 
        ones by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the start domain; if None, 
          the start domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the arrival domain; if None, 
          the arrival domain's default chart is assumed

        OUTPUT:
        
        - instance of :class:`MultiFunctionChart` representing the 
          differentiable mapping in the above two charts

        EXAMPLES:

        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
        
            sage: m = Manifold(2, 'M')
            sage: n = Manifold(3, 'N')
            sage: c_uv.<u,v> = m.chart('u v')
            sage: c_xyz.<x,y,z> = n.chart('x y z')
            sage: Phi = DiffMapping(m, n, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.multi_function_chart(c_uv, c_xyz)
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart() # equivalent to above since 'uv' and 'xyz' are default charts
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: type(Phi.multi_function_chart())
            <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>

        Representation in other charts::
        
            sage: c_UV.<U,V> = m.chart('U V')  # new chart on M
            sage: ch_uv_UV = CoordChange(c_uv, c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = n.chart('X Y Z') # new chart on N
            sage: ch_xyz_XYZ = CoordChange(c_xyz, c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.multi_function_chart(c_UV, c_xyz)
            functions (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V) on the chart (M, (U, V))
            sage: Phi.multi_function_chart(c_uv, c_XYZ)
            functions (((2*u + 1)*v^2 + u*v - 3*u)/v, -((u - 1)*v^2 - u*v - u)/v, -((u + 1)*v^2 + u*v - 2*u)/v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart(c_UV, c_XYZ)
            functions (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V)) on the chart (M, (U, V))

        """
        dom1 = self.domain1; dom2 = self.domain2
        def_chart1 = dom1.def_chart; def_chart2 = dom2.def_chart
        if chart1 is None:
            chart1 = def_chart1
        if chart2 is None:
            chart2 = def_chart2
        if (chart1, chart2) not in self.coord_expression:
            # some change of coordinates must be performed
            change_start = [] ; change_arrival = []
            for (ochart1, ochart2) in self.coord_expression:
                if chart1 == ochart1:
                    change_arrival.append(ochart2)
                if chart2 == ochart2:
                    change_start.append(ochart1)
            # 1/ Trying to make a change of chart only on the arrival domain:
            # the arrival default chart is privileged:
            sel_chart2 = None # selected chart2
            if def_chart2 in change_arrival \
                    and (def_chart2, chart2) in dom2.coord_changes:
                sel_chart2 = def_chart2
            else:
                for ochart2 in change_arrival:
                    if (ochart2, chart2) in dom2.coord_changes:
                        sel_chart2 = ochart2
                        break 
            if sel_chart2 is not None:
                oexpr = self.coord_expression[(chart1, sel_chart2)]
                chg2 = dom2.coord_changes[(sel_chart2, chart2)]
                self.coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1, *(chg2(*(oexpr.expr()))) )
                return self.coord_expression[(chart1, chart2)]

            # 2/ Trying to make a change of chart only on the start domain:
            # the start default chart is privileged:
            sel_chart1 = None # selected chart1
            if def_chart1 in change_start \
                    and (chart1, def_chart1) in dom1.coord_changes:
                sel_chart1 = def_chart1
            else:
                for ochart1 in change_start:
                    if (chart1, ochart1) in dom1.coord_changes:
                        sel_chart1 = ochart1
                        break
            if sel_chart1 is not None:
                oexpr = self.coord_expression[(sel_chart1, chart2)]
                chg1 = dom1.coord_changes[(chart1, sel_chart1)]
                self.coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1, 
                                       *(oexpr( *(chg1.transf.expr()) )) )
                return self.coord_expression[(chart1, chart2)]
                    
            # 3/ If this point is reached, it is necessary to perform some 
            # coordinate change both on the start domain and the arrival one
            # the default charts are privileged:
            if (def_chart1, def_chart2) in self.coord_expression \
                    and (chart1, def_chart1) in dom1.coord_changes \
                    and (def_chart2, chart2) in dom2.coord_changes:
                sel_chart1 = def_chart1
                sel_chart2 = def_chart2
            else:
                for (ochart1, ochart2) in self.coord_expression:
                    if (chart1, ochart1) in dom1.coord_changes \
                        and (ochart2, chart2) in dom2.coord_changes:
                        sel_chart1 = ochart1
                        sel_chart2 = ochart2
                        break
            if (sel_chart1 is not None) and (sel_chart2 is not None):
                oexpr = self.coord_expression[(sel_chart1, sel_chart2)]
                chg1 = dom1.coord_changes[(chart1, sel_chart1)]
                chg2 = dom2.coord_changes[(sel_chart2, chart2)]
                self.coord_expression[(chart1, chart2)] = \
                     MultiFunctionChart(chart1, 
                                *(chg2( *(oexpr(*(chg1.transf.expr()))) )) )
                return self.coord_expression[(chart1, chart2)]
                
            # 4/ If this point is reached, the demanded value cannot be
            # computed 
            raise ValueError("The expression of the mapping in the pair of " +
                "charts (" + str(chart1) + ", " + str(chart2) + ") cannot " + 
                "be computed by means of known changes of charts.")
                
        return self.coord_expression[(chart1, chart2)]
            

    def expr(self, chart1=None, chart2=None):
        r""" 
        Return the expression of the differentiable mapping in terms of
        specified coordinates.
        
        If the expression is not already known, it is computed from some known 
        expression by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the start domain; if None, 
          the start domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the arrival domain; if None, 
          the arrival domain's default chart is assumed

        OUTPUT:
        
        - symbolic expression representing the differentiable mapping in the 
          above two charts

        EXAMPLES:
        
        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
        
            sage: m = Manifold(2, 'M')
            sage: n = Manifold(3, 'N')
            sage: c_uv.<u,v> = m.chart('u v')
            sage: c_xyz.<x,y,z> = n.chart('x y z')
            sage: Phi = DiffMapping(m, n, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.expr(c_uv, c_xyz)
            (u*v, u/v, u + v)
            sage: Phi.expr()  # equivalent to above since 'uv' and 'xyz' are default charts
            (u*v, u/v, u + v)
            sage: type(Phi.expr()[0])
            <type 'sage.symbolic.expression.Expression'>

        Expressions in other charts::
        
            sage: c_UV.<U,V> = m.chart('U V')  # new chart on M
            sage: ch_uv_UV = CoordChange(c_uv, c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = n.chart('X Y Z') # new chart on N
            sage: ch_xyz_XYZ = CoordChange(c_xyz, c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.expr(c_UV, c_xyz)
            (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V)
            sage: Phi.expr(c_uv, c_XYZ)
            (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v,
             -((u + 1)*v^2 + u*v - 2*u)/v)
            sage: Phi.expr(c_UV, c_XYZ)
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))

        A rotation in some Euclidean plane::
        
            sage: m = Manifold(2, 'M') # the plane
            sage: c_spher.<r,ph> = m.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on the plane
            sage: rot = DiffMapping(m, m, (r, ph+pi/3), name='R') # pi/3 rotation around r=0
            sage: rot.expr()
            (r, 1/3*pi + ph)

        Expression of the rotation in terms of Cartesian coordinates::
        
            sage: c_cart.<x,y> = m.chart('x y') # Declaration of Cartesian coordinates
            sage: ch_spher_cart = CoordChange(c_spher, c_cart, r*cos(ph), r*sin(ph)) # relation to spherical coordinates
            sage: ch_spher_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))              
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: rot.expr(c_cart, c_cart)                            
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        """
        return self.multi_function_chart(chart1, chart2).expr()
            
    def set_expr(self, chart1, chart2, coord_functions): 
        r"""
        Set a new coordinate representation of the mapping.

        The expressions with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:
    
        - ``chart1`` -- chart for the coordinates on the start domain
        - ``chart2`` -- chart for the coordinates on the arrival domain
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
        
        EXAMPLES:
        
        Polar representation of a planar rotation initally defined in 
        Cartesian coordinates::
            
            sage: m = Manifold(2, 'R^2', r'\RR^2')   # Euclidean plane
            sage: c_cart.<x,y> = m.chart('x y')  # Cartesian coordinates
            sage: c_spher.<r,ph> = m.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates
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
        
        Let us use the method :meth:`set_expr` to set the 
        spherical-coordinate expression by hand::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.expr(c_spher, c_spher) 
            (r, 1/3*pi + ph)
        
        The expression in Cartesian coordinates has been lost in the dictionary :attr:`coord_expression` that stores the various representations of 
        the differentiable mapping::
             
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
            
        It is recovered by a call to :meth:`expr`::
        
            sage: rot.expr(c_cart, c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph)), (chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y))}
            
        The rotation can be applied to a point by means of either coordinate 
        system::
            
            sage: p = Point(m, (1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p.coord(c_spher) # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
            sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
                
        """
        if chart1 not in self.domain1.atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self.domain1))
        if chart2 not in self.domain2.atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self.domain2))
        self.coord_expression.clear()
        self._del_derived()
        n2 = self.domain2.manifold.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)

    def add_expr(self, chart1, chart2, coord_functions): 
        r"""
        Set a new coordinate representation of the mapping.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead. 

        INPUT:
    
        - ``chart1`` -- chart for the coordinates on the start domain
        - ``chart2`` -- chart for the coordinates on the arrival domain
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
          
        .. WARNING::
        
            If the mapping has already expressions in other charts, it 
            is the user's responsability to make sure that the expression
            to be added is consistent with them.         
        
        EXAMPLES:
        
        Polar representation of a planar rotation initally defined in 
        Cartesian coordinates::
            
            sage: m = Manifold(2, 'R^2', r'\RR^2')   # Euclidean plane
            sage: c_cart.<x,y> = m.chart('x y')  # Cartesian coordinates
            sage: c_spher.<r,ph> = m.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates
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
        
            sage: rot.expr(c_spher, c_spher) # correct output but could be simplified !
            (r,
             arctan2(1/2*(sqrt(3)*cos(ph) + sin(ph))*r, -1/2*(sqrt(3)*sin(ph) - cos(ph))*r))
        
        Therefore, we use the method :meth:`add_expr` to set the 
        spherical-coordinate expression by hand::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.expr(c_spher, c_spher)  # the output is now satisfactory
            (r, 1/3*pi + ph)
        
        The expression in Cartesian coordinates has been kept in the 
        dictionary :attr:`coord_expression` 
        that stores the various representations of the differentiable mapping::
             
            sage: rot.coord_expression
            {(chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y)), (chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
             
        If, on the contrary, we use :meth:`set_expr`, the expression in 
        Cartesian coordinates is lost::
        
            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
 
        It is recovered by a call to :meth:`expr`::
        
            sage: rot.expr(c_cart,c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph)), (chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y))}

        The rotation can be applied to a point by means of either coordinate 
        system::
            
            sage: p = Point(m, (1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p.coord(c_spher) # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
            sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
                
        """
        if chart1 not in self.domain1.atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self.domain1))
        if chart2 not in self.domain2.atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self.domain2))
        n2 = self.domain2.manifold.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)


    def __call__(self, p, chart1=None, chart2=None):
        r"""
        Compute the image of a point.

        INPUT:
    
        - ``p`` -- point on the start domain (type: :class:`Point`)
        - ``chart1`` -- (default: None) chart in which the coordinates of p 
          are to be considered; if none is provided, a chart in which both p's 
          coordinates and the expression of ``self`` are known is searched, 
          starting from the default chart of self.domain1 will be used
        - ``chart2`` -- (default: None) chart in which the coordinates of the 
          image of p will be computed; if none is provided, the default chart 
          of self.domain2 is assumed.
        
        OUTPUT:

        - image of the point by the mapping (type: :class:`Point`)

        EXAMPLES:
        
            Planar rotation acting on a point::
            
                sage: m = Manifold(2, 'R^2', r'\RR^2') # Euclidean plane
                sage: c_cart.<x,y> = m.chart('x y') # Cartesian coordinates
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
                sage: c_spher.<r,ph> = m.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
                sage: ch = CoordChange(c_cart, c_spher, sqrt(x*x+y*y), atan2(y,x))
                sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
                sage: p.coord(c_spher) # the spherical coord. of p are evaluated
                (sqrt(5), arctan(2))
                sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
                sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
                (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
                sage: q1 == q
                True
    
        """
        from manifold import RealLine
        if p not in self.domain1.manifold: 
            raise ValueError("The point " + str(p) +
                  " does not belong to the " + str(self.domain1.manifold))
        if chart2 is None: 
            chart2 = self.domain2.def_chart
        if chart1 is None: 
            def_chart1 = self.domain1.def_chart
            if def_chart1 in p.coordinates and \
                        (def_chart1, chart2) in self.coord_expression:
                chart1 = def_chart1
            else:
                for chart in p.coordinates:
                    if (chart, chart2) in self.coord_expression:
                        chart1 = chart
                        break
        if chart1 is None:
            raise ValueError("No common chart has been found to evaluate " \
                "the action of " + str(self) + " on the " + str(p) + ".")

        coord_map = self.coord_expression[(chart1, chart2)]
        y = coord_map(*(p.coordinates[chart1])) 
        
        if self.domain2.manifold is RealLine:   # special case of a mapping to R
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
            
            return Point(self.domain2.manifold, y, chart2, name=res_name, 
                         latex_name=res_latex_name)  #!# check
        
    def pullback(self, tensor):
        r""" 
        Pullback operator associated with the differentiable mapping. 
        
        INPUT:
        
        - ``tensor`` -- instance of :class:`TensorField` representing a fully 
          covariant tensor field `T` on the *arrival* domain, i.e. a tensor 
          field of type (0,p), with p a positive or zero integer. The case p=0 
          corresponds to a scalar field.
          
        OUTPUT:
        
        - instance of :class:`TensorField` representing a fully 
          covariant tensor field on the *start* domain that is the 
          pullback of `T` given by ``self``. 
          
        EXAMPLES:
        
        Pullback on `S^2` of a scalar field defined on `R^3`::
        
            sage: m = Manifold(2, 'S^2', start_index=1)
            sage: c_spher.<th,ph> = m.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coord. on S^2
            sage: n = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = n.chart('x y z') # Cartesian coord. on R^3
            sage: Phi = DiffMapping(m, n, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: f = ScalarField(n, x*y*z, name='f') ; f
            scalar field 'f' on the 3-dimensional manifold 'R^3'
            sage: f.view()
            f: (x, y, z) |--> x*y*z
            sage: pf = Phi.pullback(f) ; pf
            scalar field 'Phi_*(f)' on the 2-dimensional manifold 'S^2'
            sage: pf.view()
            Phi_*(f): (th, ph) |--> cos(ph)*cos(th)*sin(ph)*sin(th)^2
            
        Pullback on `S^2` of the standard Euclidean metric on `R^3`::
                
            sage: g = SymBilinFormField(n, 'g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.view()
            g = dx*dx + dy*dy + dz*dz
            sage: pg = Phi.pullback(g) ; pg
            field of symmetric bilinear forms 'Phi_*(g)' on the 2-dimensional manifold 'S^2'
            sage: pg.view()
            Phi_*(g) = dth*dth + sin(th)^2 dph*dph

        Pullback on `S^2` of a 3-form on `R^3`::
                
            sage: a = DiffForm(n, 3, 'A')
            sage: a[1,2,3] = f 
            sage: a.view()
            A = x*y*z dx/\dy/\dz
            sage: pa = Phi.pullback(a) ; pa
            3-form 'Phi_*(A)' on the 2-dimensional manifold 'S^2'
            sage: pa.view() # should be zero (as any 3-form on a 2-dimensional manifold)
            Phi_*(A) = 0

        """
        from scalarfield import ScalarField
        from vectorframe import CoordFrame
        from rank2field import SymBilinFormField
        from diffform import DiffForm, OneForm

        if not isinstance(tensor, TensorField):
            raise TypeError("The argument 'tensor' must be a tensor field.")
        dom1 = self.domain1
        dom2 = self.domain2
        if not tensor.domain.is_subdomain(dom2):
            raise TypeError("The tensor field is not defined on the mapping " +
                            "arrival domain.")
        (ncon, ncov) = tensor.tensor_type
        if ncon != 0:
            raise TypeError("The pullback cannot be taken on a tensor " + 
                            "with some contravariant part.")
        resu_name = None ; resu_latex_name = None
        if self.name is not None and tensor.name is not None:
            resu_name = self.name + '_*(' + tensor.name + ')'
        if self.latex_name is not None and tensor.latex_name is not None:
            resu_latex_name = self.latex_name + '_*' + tensor.latex_name                
        if ncov == 0:
            # Case of a scalar field
            # ----------------------
            resu = ScalarField(dom1, name=resu_name, 
                               latex_name=resu_latex_name)
            for chart2 in tensor.express:
                for chart1 in dom1.atlas:
                    if (chart1, chart2) in self.coord_expression:
                        phi = self.coord_expression[(chart1, chart2)]
                        coord1 = chart1.xx
                        ff = tensor.express[chart2]
                        resu.add_expr( ff(*(phi(*coord1))), chart1)
            return resu
        else:
            # Case of tensor field of rank >= 1
            # ---------------------------------
            if isinstance(tensor, OneForm):
                resu = OneForm(dom1, name=resu_name, 
                               latex_name=resu_latex_name)
            elif isinstance(tensor, DiffForm):
                resu = DiffForm(dom1, ncov, name=resu_name, 
                                latex_name=resu_latex_name)
            elif isinstance(tensor, SymBilinFormField):
                resu = SymBilinFormField(dom1, name=resu_name, 
                                         latex_name=resu_latex_name)                
            else:
                resu = TensorField(dom1, 0, ncov, name=resu_name, 
                                   latex_name=resu_latex_name, sym=tensor.sym,
                                   antisym=tensor.antisym)
            for frame2 in tensor.components:
                if isinstance(frame2, CoordFrame):
                    chart2 = frame2.chart
                    for chart1 in dom1.atlas:
                        if (chart1, chart2) in self.coord_expression:
                            # Computation at the component level:
                            frame1 = chart1.frame
                            tcomp = tensor.components[frame2]
                            if isinstance(tcomp, CompFullySym):
                                ptcomp = CompFullySym(frame1, ncov)
                            elif isinstance(tcomp, CompFullyAntiSym):
                                ptcomp = CompFullyAntiSym(frame1, ncov)
                            elif isinstance(tcomp, CompWithSym):
                                ptcomp = CompWithSym(frame1, ncov, sym=tcomp.sym, 
                                                     antisym=tcomp.antisym)
                            else:
                                ptcomp = Components(frame1, ncov)
                            phi = self.coord_expression[(chart1, chart2)]
                            jacob = phi.jacobian()
                            # X2 coordinates expressed in terms of X1 ones via the mapping:
                            coord2_1 = phi(*(chart1.xx)) 
                            si1 = dom1.manifold.sindex
                            si2 = dom2.manifold.sindex
                            for ind_new in ptcomp.non_redundant_index_generator(): 
                                res = 0 
                                for ind_old in dom2.manifold.index_generator(ncov): 
                                    ff = tcomp[[ind_old]].function_chart(chart2)
                                    t = FunctionChart(chart1, ff(*coord2_1))
                                    for i in range(ncov):
                                        t *= jacob[ind_old[i]-si2][ind_new[i]-si1]
                                    res += t
                                ptcomp[ind_new] = res
                            resu.components[frame1] = ptcomp
            return resu

        
#*****************************************************************************

class Diffeomorphism(DiffMapping):
    r"""
    Class for manifold diffeomorphisms.

    INPUT:
    
    - ``domain1`` -- domain on the start manifold 
    - ``domain2`` -- domain on the arrival manifold 
    - ``coord_functions`` -- the coordinate symbolic expression of the mapping: 
      list (or tuple) of the coordinates of the image expressed in terms of the
      coordinates of the considered point
    - ``chart1`` -- (default: None) chart in which the coordinates are given on
      domain1; if none is provided, the coordinates are assumed to refer to the
      domain's default chart
    - ``chart2`` -- (default: None) chart in which the coordinates are given on
      domain2; if none is provided, the coordinates are assumed to refer to the
      domain's default chart
    - ``name`` -- (default: None) name given to the differentiable mapping
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      differentiable mapping; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    """
    def __init__(self, domain1, domain2, coord_functions=None, 
                 chart1=None, chart2=None, name=None, latex_name=None): 
        DiffMapping.__init__(self, domain1, domain2, coord_functions, chart1, 
                             chart2, name, latex_name)
        if self.domain1.manifold.dim != self.domain2.manifold.dim:
            raise ValueError("The manifolds " + str(self.domain1.manifold) + 
                             " and " + str(self.domain2.manifold) + 
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
        if self.domain1 == self.domain2:
            description += " on the " + str(self.domain1)
        else:
            description += " between the " + str(self.domain1) + \
                           " and the " + str(self.domain2)
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


    def inverse(self, chart1=None, chart2=None): 
        r"""
        Return the inverse diffeomorphism. 
        
        INPUT:
    
        - ``chart1`` -- (default: None) chart in which the computation of the 
          inverse is performed if necessary; if none is provided, the default 
          chart of the start domain will be used
        - ``chart2`` -- (default: None) chart in which the computation of the 
          inverse is performed if necessary; if none is provided, the default 
          chart of the arrival domain will be used
        
        OUTPUT:
        
        - the inverse diffeomorphism
        
        EXAMPLES:
        
            The inverse of a rotation in the Euclidean plane::
            
                sage: m = Manifold(2, 'R^2', r'\RR^2')
                sage: c_cart.<x,y> = m.chart('x y')
                sage: # A pi/3 rotation around the origin:
                sage: rot = Diffeomorphism(m, m, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
                sage: rot.inverse() 
                diffeomorphism 'R^(-1)' on the 2-dimensional manifold 'R^2'
                sage: rot.inverse().view()
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
            
        if chart1 is None: chart1 = self.domain1.def_chart
        if chart2 is None: chart2 = self.domain2.def_chart
        coord_map = self.coord_expression[(chart1, chart2)]
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
        self._inverse = Diffeomorphism(self.domain2, self.domain1, 
                                       inv_functions, chart2, chart1,
                                       name=name, latex_name=latex_name)
        return self._inverse
