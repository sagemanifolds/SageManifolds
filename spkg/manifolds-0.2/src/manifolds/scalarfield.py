r"""
Scalar fields

The class :class:`ScalarField` implements scalar fields on differentiable 
manifolds over `\RR`. 

It inherits from the classes :class:`DiffMapping` (a scalar
field being a differentiable mapping to `\RR`) and
:class:`DiffForm` (a scalar field being a differential form of degree 0). 

The subclass :class:`ZeroScalarField` deals with null scalar fields. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

"""

#******************************************************************************
#       Copyright (C) 2013 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from manifold import Manifold
from chart import FunctionChart, ZeroFunctionChart, MultiFunctionChart
from diffmapping import DiffMapping
from diffform import DiffForm, OneForm

class ScalarField(DiffMapping, DiffForm):
    r"""
    Class for scalar fields on a differentiable manifold.
    
    INPUT:
    
    - ``manifold`` -- the manifold on which the scalar field is defined
    - ``coord_expression`` -- (default: None) coordinate expression of the 
      scalar field
    - ``chart_name`` -- (default:None) name of the chart defining the 
      coordinates used in ``coord_expression``; if none is provided and a
      coordinate expression is given, the manifold default chart is assumed.
    - ``name`` -- (default: None) name given to the scalar field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the scalar field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A scalar field on the 2-sphere::
    
        sage: m = Manifold(2, 'S^2')
        sage: f = ScalarField(m) ; f
        scalar field on the 2-dimensional manifold 'S^2'
        
    Named scalar field::
    
        sage: f = ScalarField(m, name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: latex(f)
        f

    Named scalar field with LaTeX symbol specified::
    
        sage: f = ScalarField(m, name='f', latex_name=r'\mathcal{F}') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: latex(f)
        \mathcal{F}
        
    Scalar field defined by its coordinate expression::
    
        sage: c_spher = Chart(m, r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        sage: f = ScalarField(m, sin(th)*cos(ph), name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'

    The coordinate expression of a scalar field can be read by means of 
    :meth:`expr` and set by means of :meth:`set_expr`; both methods can 
    have a chart name as argument (if not, the manifold's default chart is 
    assumed)::
    
        sage: f.expr()
        cos(ph)*sin(th)
        sage: f.expr('spher') # equivalent to above since 'spher' is the default chart
        cos(ph)*sin(th)
        sage: f.set_expr(cos(th))  # changing the value of f
        sage: f.expr()
        cos(th)
        sage: f.set_expr(sin(th)*cos(ph)) # restoring the original value
        
    The function :meth:`show` displays the coordinate expression of the scalar
    field::
    
        sage: f.show()
        f: (th, ph) |--> cos(ph)*sin(th)
        sage: f.show('spher') # equivalent to above since 'spher' is the default chart
        f: (th, ph) |-->  cos(ph)*sin(th)
        sage: latex(f.show()) # nice LaTeX formatting for the notebook
        f :\ \left(\theta, \phi\right) \mapsto \cos\left(\phi\right) \sin\left(\theta\right)

    A scalar field can also be defined by an unspecified function of the 
    coordinates::
    
        sage: g = ScalarField(m, function('G', th, ph), name='g') ; g
        scalar field 'g' on the 2-dimensional manifold 'S^2'
        sage: g.expr()
        G(th, ph)
        sage: s = f+g ; s.expr()                               
        cos(ph)*sin(th) + G(th, ph)
        
    In each chart, the scalar field is represented by a function of the 
    coordinates, which is a an instance of the class :class:`FunctionChart` 
    and can be accessed by the method :meth:`function_chart`::
    
        sage: f.function_chart('spher')
        cos(ph)*sin(th)
        sage: f.function_chart() # equivalent to above since 'spher' is the default chart
        cos(ph)*sin(th)
        sage: print type(f.function_chart())
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        
    The value returned by the method :meth:`expr` is actually the coordinate
    expression of the function::

        sage: f.expr() is f.function_chart().expr()
        True

    A scalar field is a tensor field of rank 0::
    
        sage: isinstance(f, TensorField)
        True
        sage: f.rank
        0
        
    As such, it is a differential form of degree 0::
    
        sage: isinstance(f, DiffForm)
        True
        
    It is also a differential mapping from the manifold to the field of
    real numbers (modeled by the unique instance :data:`RealLine` of the class
    :class:`RealLineManifold`)::
    
        sage: isinstance(f, DiffMapping)
        True
        sage: f.manifold1 # the start manifold
        2-dimensional manifold 'S^2'
        sage: f.manifold2 # the target manifold
        field R of real numbers
        sage: print type(f.manifold2)
        <class 'sage.geometry.manifolds.manifold.RealLineManifold'>
        sage: f.manifold2 is RealLine
        True
        sage: f.coord_expression
        {('spher', 'canonical'): functions (cos(ph)*sin(th),) on the chart 'spher' (th, ph)}
        sage: f.expr()  # expression with respect to the manifold's default chart
        cos(ph)*sin(th)

    As such, it acts on the manifold's points::
    
        sage: p = Point(m, (pi/2, pi))
        sage: f(p)
        -1

    A scalar field can be compared to another scalar field::
    
        sage: g = ScalarField(m, sin(th)*cos(ph), name='g')
        sage: f == g
        True
        sage: g.set_expr(cos(th))
        sage: f == g
        False

    ...to an instance of :class:`FunctionChart`::
    
        sage: h = FunctionChart(c_spher, sin(th)*cos(ph)) ; h
        cos(ph)*sin(th)
        sage: f == h
        True
        sage: h = FunctionChart(c_spher, cos(th))
        sage: f == h
        False
        
    ...to a symbolic expression::
    
        sage: f == sin(th)*cos(ph)
        True
        sage: f == ph + th^2
        False
        
    ...to a number::
        
        sage: f == 2
        False

    ...to zero::
    
        sage: f == 0
        False
        sage: f.set_expr(0)
        sage: f == 0
        True

    ...to anything else::
    
        sage: f == m
        False

    Scalar fields can be added::
    
        sage: f.set_expr(sin(th)*cos(ph))
        sage: g.set_expr(cos(th))
        sage: s = f + g ; s
        scalar field 'f+g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + cos(th)
        sage: s = f + cos(th) ; s # direct addition with a symbolic expression is allowed
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + cos(th)
        sage: s = 1 + f ; s 
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) + 1
        sage: s = +f ; s  # the unary plus operator
        scalar field '+f' on the 2-dimensional manifold 'S^2'
        sage: s == f
        True
      
    Scalar fields can be subtracted::
    
        sage: s = f - g ; s
        scalar field 'f-g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) - cos(th)
        sage: s = f - cos(th) ; s  # direct subtraction of a symbolic expression is allowed
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*sin(th) - cos(th)
        sage: s = 1 - f ; s
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        -cos(ph)*sin(th) + 1
        sage: s = f - g + (g - f)
        sage: s == 0
        True
        sage: s = f + (-f) # check of the unary minus operator
        sage: s == 0
        True
     
    Scalar fields can be multiplied and divided::
     
        sage: s = f*g ; s
        scalar field 'f*g' on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*cos(th)*sin(th)
        sage: h = s / g ; h 
        scalar field 'f*g/g' on the 2-dimensional manifold 'S^2'
        sage: h.expr()
        cos(ph)*sin(th)
        sage: h == f
        True
        sage: h1 = s / f ; h1 
        scalar field 'f*g/f' on the 2-dimensional manifold 'S^2'
        sage: h1.expr()
        cos(th)
        sage: h1 == g
        True
            
    The multiplication and division can be performed by a symbolic expression::
    
        sage: s = f*cos(th) ; s
        scalar field on the 2-dimensional manifold 'S^2'
        sage: s.expr()
        cos(ph)*cos(th)*sin(th)
        sage: h = s/cos(th) ; h
        scalar field on the 2-dimensional manifold 'S^2'
        sage: h.expr()
        cos(ph)*sin(th)
        sage: h == f
        True        

    The in-place operators +=, -=, \*= and /= are implemented::
    
        sage: f.expr()
        cos(ph)*sin(th)
        sage: f += cos(th)
        sage: f.expr()    
        cos(ph)*sin(th) + cos(th)
        sage: f -= cos(th)
        sage: f.expr()    
        cos(ph)*sin(th)
        sage: f *= cos(th)
        sage: f.expr()
        cos(ph)*cos(th)*sin(th)
        sage: f /= cos(th)
        sage: f.expr()    
        cos(ph)*sin(th)

    """
    def __init__(self, manifold, coord_expression=None, chart_name=None, 
                 name=None, latex_name=None):
        from manifold import RealLine
        DiffMapping.__init__(self, manifold, RealLine, coord_expression, 
                             chart_name)
        DiffForm.__init__(self, manifold, 0, name, latex_name)
        if coord_expression is not None:
            if chart_name is None:
                chart = self.manifold.def_chart
            else:
                chart = self.manifold.atlas[chart_name]
            if coord_expression == 0:
                self.express = {chart.name: chart.zero_function}
            else:
                self.express = {chart.name: FunctionChart(chart, 
                                                          coord_expression)}
        else:
            self.express = {}


    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "scalar field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.manifold)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`ScalarField` instance on the same manifold.
        
        """
        return ScalarField(self.manifold)        

    # There is no function _init_derived because ScalarField has no derived 
    # quantity per se
    # The function _del_derived is mandatory because ScalarField has two mother
    # classes:
    
    def _del_derived(self):
        r"""
        Delete the derived quantities.
        """
        DiffMapping._del_derived(self) # derived quantities of the 1st mother class
        DiffForm._del_derived(self) # derived quantities of the 2nd mother class

    def copy(self):
        r"""
        Return an exact copy of ``self``.
        
        EXAMPLES:
        
        Copy on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')  
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: f = ScalarField(m, x*y^2)
            sage: g = f.copy()
            sage: print type(g)
            <class 'sage.geometry.manifolds.scalarfield.ScalarField'>
            sage: g.expr()
            x*y^2
            sage: g == f
            True
            sage: g is f
            False
        
        """
        result = ScalarField(self.manifold)  #!# what about the name ?
        for chart_name, funct in self.express.items():
            result.express[chart_name] = funct.copy()
        for key, mfunct in self.coord_expression.items():
            result.coord_expression[key] = mfunct.copy()
        return result

    def comp(self, frame_name=None, from_frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.comp`. 
        
        This method should not be called, since a scalar field has no 
        component w.r.t. a vector frame ! 
        It therefore returns None. 
        
        """
        return None

    def set_comp(self, frame_name=None, delete_others=True):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.set_comp`.
        
        This method should not be called, since a scalar field has no 
        component w.r.t. a vector frame ! 
        It therefore returns None. 
                
        """
        return None

    def function_chart(self, chart_name=None, from_chart=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart; if none, the 
          manifold's default chart will be used
        - ``from_chart`` -- (default: None) name of the chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart_name``; if none, a chart is picked in ``self.express``

        OUTPUT:
        
        - instance of :class:`FunctionChart` representing the coordinate 
          function of the scalar field in the given chart.

        EXAMPLES:
        
        Coordinate function on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')            
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: f = ScalarField(m, x*y^2)
            sage: f.function_chart()
            x*y^2
            sage: f.function_chart('coord_xy')  # equivalent form (since 'coord_xy' is the default chart)
            x*y^2
            sage: print type(f.function_chart())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>

        Expression via a change of coordinates::
        
            sage: c_uv = Chart(m, 'u v', 'coord_uv')
            sage: CoordChange(c_uv, c_xy, u+v, u-v)
            coordinate change from chart 'coord_uv' (u, v) to chart 'coord_xy' (x, y)
            sage: f.express # at this stage, f is expressed only in terms of (x,y) coordinates
            {'coord_xy': x*y^2}
            sage: f.function_chart('coord_uv') # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: f.function_chart('coord_uv') == (u+v)*(u-v)^2  # check
            True
            sage: f.express  # f has now 2 coordinate expressions:
            {'coord_uv': u^3 - u^2*v - u*v^2 + v^3, 'coord_xy': x*y^2}

        Usage in a physical context (simple Lorentz transformation - boost in 
        x direction, with relative velocity v between o1 and o2 frames)::

            sage: m = Manifold(2, 'M')
            sage: o1 = Chart(m, 't x', 'o1')
            sage: o2 = Chart(m, 'T X', 'o2')
            sage: f = ScalarField(m, x^2 - t^2)
            sage: f.function_chart('o1')
            -t^2 + x^2
            sage: v = var('v'); gam = 1/sqrt(1-v^2)
            sage: CoordChange(o2, o1, gam*(T - v*X), gam*(X - v*T))
            coordinate change from chart 'o2' (T, X) to chart 'o1' (t, x)
            sage: f.function_chart('o2')
            -T^2 + X^2

        """
        if chart_name is None:
            chart_name = self.manifold.def_chart.name
        if chart_name not in self.express:
            # the expression must be computed from that in the chart from_chart
            if from_chart is None:
                from_chart = self.pick_a_chart()
            #!# instead of picking a chart, one should select one from which 
            # the coordinate change is possible
            change = self.manifold.coord_changes[(chart_name, from_chart)]
            # old coordinates expressed in terms of the new ones:
            coords = [ change.transf.functions[i].express 
                       for i in range(self.manifold.dim) ]
            new_expr = self.express[from_chart](*coords)
            new_chart = self.manifold.atlas[chart_name]
            self.express[chart_name] = FunctionChart(new_chart, new_expr)
            self.coord_expression[(chart_name,'canonical')] = \
                                        MultiFunctionChart(new_chart, new_expr)
            self._del_derived()
        return self.express[chart_name]


    def expr(self, chart_name=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart; if none, the 
          manifold's default chart will be used
        - ``from_chart`` -- (default: None) name of the chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart_name``; if none, a chart is picked in ``self.express``
          
        OUTPUT:
        
        - symbolic expression representing the coordinate 
          expression of the scalar field in the given chart.
        
        EXAMPLES:
        
        Expression of a scalar field on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')            
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: f = ScalarField(m, x*y^2)
            sage: f.expr()
            x*y^2
            sage: f.expr('coord_xy')  # equivalent form (since 'coord_xy' is the default chart)
            x*y^2
            sage: print type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        Expression via a change of coordinates::
        
            sage: c_uv = Chart(m, 'u v', 'coord_uv')
            sage: CoordChange(c_uv, c_xy, u+v, u-v)
            coordinate change from chart 'coord_uv' (u, v) to chart 'coord_xy' (x, y)
            sage: f.express # at this stage, f is expressed only in terms of (x,y) coordinates
            {'coord_xy': x*y^2}
            sage: f.expr('coord_uv') # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: bool( f.expr('coord_uv') == (u+v)*(u-v)^2 ) # check
            True
            sage: f.express  # f has now 2 coordinate expressions:
            {'coord_uv': u^3 - u^2*v - u*v^2 + v^3, 'coord_xy': x*y^2}

        """
        return self.function_chart(chart_name, from_chart).express
        
    def set_expr(self, coord_expression, chart_name=None, delete_others=True):
        r"""
        Set some coordinate expression of the scalar field.
        
        NB: the default value of ``delete_others`` ensures that the expressions 
        with respect to other charts are deleted, in order to avoid any 
        inconsistency. 
        
        INPUT:
        
        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart_name`` -- (default: None) name of the chart in which 
          ``coord_expression`` is defined; if none, the manifold's default 
          chart is assumed
        - ``delete_others`` -- (default: True) determines whether the 
          coordinate expressions with respect to charts different from that 
          specified by ``chart_name`` are deleted or not. 
          
        .. WARNING::
        
            Setting ``delete_others`` to False is at the responsability of the 
            user, who must make sure that the various expressions are
            consistent with each other. 
        
        EXAMPLES:
        
        Setting scalar field expressions on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: f = ScalarField(m, x^2 + 2*x*y +1)
            sage: f.express                         
            {'coord_xy': x^2 + 2*x*y + 1}
            sage: f.set_expr(3*y)
            sage: f.express  # the (x,y) expression has been changed:
            {'coord_xy': 3*y}
            sage: c_uv = Chart(m, 'u v', 'coord_uv')
            sage: f.set_expr(cos(u)-sin(v), 'coord_uv')  
            sage: f.express # the (x,y) expression has been lost:
            {'coord_uv': cos(u) - sin(v)}
            sage: f.set_expr(3*y, delete_others=False)    
            sage: f.express # the (u,v) expression has been kept:                    
            {'coord_uv': cos(u) - sin(v), 'coord_xy': 3*y}

        """
        if chart_name is None:
            chart_name = self.manifold.def_chart.name
        chart = self.manifold.atlas[chart_name]
        if delete_others:
            self.express.clear()
            self.coord_expression.clear()
        self.express[chart_name] = FunctionChart(chart, coord_expression)
        self.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(chart, coord_expression)
        self._del_derived()

    def show(self, chart_name=None):
        r""" 
        Displays the expression of the scalar field in a given chart. 
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart; if none, the 
          manifold's default chart will be used
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Various displays::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: f = ScalarField(m, sqrt(x+1), name='f')
            sage: f.show()
            f: (x, y) |--> sqrt(x + 1)
            sage: latex(f.show())
            f :\ \left(x, y\right) \mapsto \sqrt{x + 1}
            sage: g = ScalarField(m, function('G', x, y), name='g')
            sage: g.show() 
            g: (x, y) |--> G(x, y)
            sage: latex(g.show())
            g :\ \left(x, y\right) \mapsto G\left(x, y\right)

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if chart_name is None:
            chart_name = self.manifold.def_chart.name
        expression = self.expr(chart_name)
        chart = self.manifold.atlas[chart_name]
        if self.name is None:
            result.txt = repr(chart()) + " |--> " + repr(expression)
        else:
            result.txt = self.name + ": " + repr(chart()) + " |--> " + \
                         repr(expression)
        if self.latex_name is None:
            result.latex = latex(chart()) + r"\mapsto" + latex(expression)
        else:
            result.latex = latex(self) + ":\ " + latex(chart()) + r"\mapsto" + \
                           latex(expression)
        return result


    def pick_a_chart(self):
        r"""
        Return a chart for which the scalar field has an expression. 
        
        The manifold's default chart is privileged. 
        
        OUPUT:
        
        - name of the chart

        EXAMPLES:
        
        A very simple example::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: c_uv = Chart(m, 'u v', 'uv-coord')
            sage: f =  ScalarField(m, x*y^2)
            sage: f.set_expr(u-v, 'uv-coord', delete_others=False)
            sage: f.express # f has expressions on two charts:
            {'xy-coord': x*y^2, 'uv-coord': u - v}
            sage: m.default_chart()
            chart 'xy-coord' (x, y)
            sage: f.pick_a_chart()  # the manifold's default chart (xy-coord) is privileged:
            'xy-coord'
            sage: g = ScalarField(m, u+v, 'uv-coord')
            sage: g.express  # g has no expression on the manifold's default chart:
            {'uv-coord': u + v}
            sage: g.pick_a_chart()
            'uv-coord'
        
        """
        if self.manifold.def_chart.name in self.express:
            return self.manifold.def_chart.name
        return self.express.items()[0][0]  
            

    def common_chart(self, other):
        r"""
        Find a common chart for the expressions of ``self`` and ``other``. 
        
        In case of multiple charts, the manifold's default chart is 
        privileged. 
        If the current expressions of ``self`` and ``other`` are all relative to
        different charts, a common chart is searched by performing an expression
        transformation, via the transformations listed in 
        ``self.manifold.coord_changes``, still privileging transformations to 
        the manifold's default chart.
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUPUT:
        
        - name of the common chart; if a common chart is not found, None is 
          returned. 

        EXAMPLES:
        
        Search for a common chart on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: c_uv = Chart(m, 'u v', 'coord_uv')
            sage: f = ScalarField(m, x^2 + 2*x*y +1)
            sage: g = ScalarField(m, x+y)
            sage: f.common_chart(g)
            'coord_xy'
            sage: h = ScalarField(m, u+v, 'coord_uv')
            sage: f.common_chart(h)  # none is found
            sage: f.set_expr(u*v, 'coord_uv', delete_others=False)
            sage: f.common_chart(h)  # now, there is a common chart: 
            'coord_uv'
            sage: f.express  # indeed: 
            {'coord_uv': u*v, 'coord_xy': x^2 + 2*x*y + 1}
            sage: h.set_expr(x*y, delete_others=False) # setting h on (x,y) chart 
            sage: h.express  # as f, h is now defined on both charts:
            {'coord_uv': u + v, 'coord_xy': x*y}
            sage: f.common_chart(h) # the method common_chart privileges the manifold's default frame:  
            'coord_xy'

        """
        # Compatibility checks:
        if not isinstance(other, ScalarField):
            raise TypeError("The second argument must be a scalar field.")
        manif = self.manifold
        if other.manifold != manif:
            raise TypeError("The two scalar fields are not defined on the " +
                            "same manifold.")
        #
        # NB: in all that follows, the denomination "chart" stands actually for
        # a chart name. 
        def_chart = manif.def_chart.name
        #
        # 1/ Search for a common chart among the existing expressions, i.e. 
        #    without performing any expression transformation. 
        #    -------------------------------------------------------------
        if def_chart in self.express and def_chart in other.express:
            return def_chart # the manifold's default chart is privileged
        for schart in self.express:
            if schart in other.express:
                return schart
        #
        # 2/ Search for a common chart via one expression transformation
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least 
        # one expression transformation to get a common chart
        if def_chart in self.express:
            for ochart in other.express:
                if (ochart, def_chart) in manif.coord_changes:
                    other.function_chart(def_chart, from_chart=ochart)
                    return def_chart
        if def_chart in other.express:
            for schart in self.express:
                if (schart, def_chart) in manif.coord_changes:
                    self.function_chart(def_chart, from_chart=schart)
                    return def_chart
        # If this point is reached, then def_chart cannot be a common chart
        # via a single expression transformation
        for schart in self.express:
            for ochart in other.express:
                if (ochart, schart) in manif.coord_changes:
                    other.function_chart(schart, from_chart=ochart)
                    return schart
                if (schart, ochart) in manif.coord_changes:
                    self.function_chart(ochart, from_chart=schart)
                    return ochart
        #
        # 3/ Search for a common chart via two expression transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # expression transformation to get a common chart
        for schart in self.express:
            for ochart in other.express:
                if (schart, def_chart) in manif.coord_changes and \
                   (ochart, def_chart) in manif.coord_changes:
                    self.function_chart(def_chart, from_chart=schart)
                    other.function_chart(def_chart, from_chart=ochart)
                    return def_chart
                for chart in manif.atlas:
                    if (schart, chart) in manif.coord_changes and \
                       (ochart, chart) in manif.coord_changes:
                        self.function_chart(chart, from_chart=schart)
                        other.function_chart(chart, from_chart=ochart)
                        return chart
        #
        # If this point is reached, no common chart could be found, even at 
        # the price of expression transformations:
        return None
   

    def is_zero(self):
        r""" 
        Return True if the scalar field is zero and False otherwise.

        EXAMPLES:
        
        Tests on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'coord_xy')
            sage: f = ScalarField(m, x*y)
            sage: f.is_zero()
            False
            sage: f.set_expr(0)
            sage: f.is_zero()
            True
            sage: g = ScalarField(m, 0)
            sage: g.is_zero()
            True

        """
        chart_name = self.pick_a_chart()
        return self.express[chart_name].is_zero()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a scalar field or something else
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        from sage.symbolic.expression import Expression
        from sage.symbolic.ring import SR
        if isinstance(other, ScalarField):
            if other.manifold != self.manifold:
                return False
            chart_name = self.common_chart(other)
            if chart_name is None:
                raise ValueError("No common chart for the comparison.")
            return bool(self.express[chart_name] == other.express[chart_name])
        elif isinstance(other, FunctionChart):
            chart_name = other.chart.name
            if chart_name  not in self.express:
                return False
            else:
                return bool(self.express[chart_name] == other.express)
        elif isinstance(other, Expression):
            var_other = set(other.variables())
            for val in self.express.values():
                if var_other.issubset(set(val.express.variables())):
                    return bool(val == other)
        else:
            try:
                xother = SR(other)
                result = False
                for val in self.express.values():
                    if val == xother:
                        result = True
                return result
            except TypeError:
                return False

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- a scalar field or something else
        
        OUTPUT:
        
        - True if ``self`` is different from ``other``,  or False otherwise
        
        """
        return not self.__eq__(other)

    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - an exact copy of ``self``
    
        """
        result = self._new_instance()
        for chart_name in self.express:
            res = + self.express[chart_name]
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
        if self.name is not None:
            result.name = '+' + self.name 
        if self.latex_name is not None:
            result.latex_name = '+' + self.latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the negative of ``self``
    
        """
        result = self._new_instance()
        for chart_name in self.express:
            res = - self.express[chart_name]
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
        if self.name is not None:
            result.name = '-' + self.name 
        if self.latex_name is not None:
            result.latex_name = '-' + self.latex_name
        return result

    def __add__(self, other):
        r"""
        Scalar field addition. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the addition of ``self`` and 
          ``other``
        
        """
        if isinstance(other, ScalarField):
            if isinstance(other, ZeroScalarField):
                return self.copy()
            chart_name = self.common_chart(other)
            if chart_name is None:
                raise ValueError("No common chart for the addition.")
            # FunctionChart addition: 
            res = self.express[chart_name] + other.express[chart_name]
            if res.is_zero():
                return self.manifold.zero_scalar_field
            result = self._new_instance()
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
            if self.name is not None and other.name is not None:
                result.name = self.name + '+' + other.name
            if self.latex_name is not None and other.latex_name is not None:
                result.latex_name = self.latex_name + '+' + other.latex_name
            return result
        else: # other is not a ScalarField:
            return self + scalar_field(other, self.manifold)

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        """
        return self.__add__(other)

    def __iadd__(self, other):
        r"""
        In-place addition operator.  
                        
        """
        return self.__add__(other)

    def __sub__(self, other):
        r"""
        Scalar field subtraction. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the subtraction of ``other`` from 
          ``self``

        """
        if isinstance(other, ScalarField):
            if isinstance(other, ZeroScalarField):
                return self.copy()
            chart_name = self.common_chart(other)
            if chart_name is None:
                raise ValueError("No common chart for the subtraction.")
            # FunctionChart subtraction: 
            res = self.express[chart_name] - other.express[chart_name]
            if res.is_zero():
                return self.manifold.zero_scalar_field
            result = self._new_instance()
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
            if self.name is not None and other.name is not None:
                result.name = self.name + '-' + other.name
            if self.latex_name is not None and other.latex_name is not None:
                result.latex_name = self.latex_name + '-' + other.latex_name
            return result
        else: # other is not a ScalarField:
            return self - scalar_field(other, self.manifold)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 
        
        """
        return (-self).__add__(other)

    def __isub__(self, other):
        r"""
        In-place subtraction operator.  
                        
        """
        return self.__sub__(other)
        

    def __mul__(self, other):
        r"""
        Scalar field multiplication. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the multiplication of ``self`` by 
          ``other``
        
        """
        from component import Components
        from tensorfield import TensorField
        from utilities import format_mul_txt, format_mul_latex        
        if isinstance(other, ScalarField):
            if isinstance(other, ZeroScalarField):
                return other
            chart_name = self.common_chart(other)
            if chart_name is None:
                raise ValueError("No common chart for the multiplication.")
            # FunctionChart multiplication: 
            res = self.express[chart_name] * other.express[chart_name]
            if res.is_zero():
                return self.manifold.zero_scalar_field
            result = self._new_instance()
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, ' ', 
                                                 other.latex_name)
        
        elif isinstance(other, TensorField):
            result = other._new_instance()
            for frame_name in other.components:
                result.components[frame_name] = self * \
                                                other.components[frame_name]
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, ' ', 
                                                 other.latex_name)

        elif isinstance(other, Components):
            result = other._new_instance()
            for ind in other._comp:
                result._comp[ind] = self * other._comp[ind]

        else: # other is not a tensor field:
            return self * scalar_field(other, self.manifold)

        return result

    def __rmul__(self, other):
        r"""
        Multiplication on the left with ``other``. 
        
        """
        # Since the multiplication with a scalar field must be commutative:
        return self.__mul__(other)  

    def __imul__(self, other):
        r"""
        In-place multiplication operator.  
                        
        """
        return self.__mul__(other)

    def __div__(self, other):
        r"""
        Scalar field division. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the division of ``self`` by 
          ``other``
        
        """
        from utilities import format_mul_txt, format_mul_latex
        if isinstance(other, ScalarField):
            if isinstance(other, ZeroScalarField):
                raise ZeroDivisionError("Division of a scalar field by zero.")
            chart_name = self.common_chart(other)
            if chart_name is None:
                raise ValueError("No common chart for the division.")
            # FunctionChart division: 
            res = self.express[chart_name] / other.express[chart_name]
            if res.is_zero():
                return self.manifold.zero_scalar_field
            result = self._new_instance()
            result.express[chart_name] = res
            result.coord_expression[(chart_name,'canonical')] = \
                                    MultiFunctionChart(res.chart, res.express)
            result.name = format_mul_txt(self.name, '/', other.name)
            result.latex_name = format_mul_latex(self.latex_name, '/', 
                                                 other.latex_name)
            return result
        else: # other is not a ScalarField:
            return self / scalar_field(other, self.manifold)


    def __rdiv__(self, other):
        r"""
        Division of ``other`` by ``self``. 
        
        """
        raise NotImplementedError("Operator ScalarField.__rdiv__ not " + 
                                  "implemented.")
                                  
    def __idiv__(self, other):
        r"""
        In-place division operator.  
                        
        """
        return self.__div__(other)

    def exterior_der(self, chart_name=None):
        r"""
        Return the exterior derivative of the scalar field. 
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart used for the
          computation; if none is provided, a chart on which the scalar field
          is expressed is picked, privileging the manifold's defaut chart. 
        
        OUTPUT:
        
        - the 1-form exterior derivative of ``self``. 
        
        EXAMPLES:
        
        Exterior derivative on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = Chart(m, 'x y z', 'coord_xyz')
            sage: f = ScalarField(m, cos(x)*z^3 + exp(y)*z^2, name='f')
            sage: df = f.exterior_der() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.show()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        Exterior derivative computed on a chart that is not the default one::
        
            sage: c_uvw = Chart(m, 'u v w', 'coord_uvw')
            sage: g = ScalarField(m, u*v^2*w^3, 'coord_uvw', name='g')
            sage: dg = g.exterior_der() ; dg
            1-form 'dg' on the 3-dimensional manifold 'M'
            sage: dg.components
            {'coord_uvw_b': 1-index components w.r.t. the coordinate basis 'coord_uvw_b' (d/du,d/dv,d/dw)}
            sage: dg.comp('coord_uvw_b')[:, 'coord_uvw']
            [v^2*w^3, 2*u*v*w^3, 3*u*v^2*w^2]
            sage: dg.show('coord_uvw_b', 'coord_uvw')
            dg = v^2*w^3 du + 2*u*v*w^3 dv + 3*u*v^2*w^2 dw
            
        The exterior derivative is nilpotent::
        
            sage: ddf = df.exterior_der() ; ddf
            2-form 'ddf' on the 3-dimensional manifold 'M'
            sage: ddf == 0
            True
            sage: ddf[:]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: ddg = dg.exterior_der() ; ddg
            2-form 'ddg' on the 3-dimensional manifold 'M'
            sage: ddg == 0
            True

        """
        from component import Components
        from utilities import format_unop_txt, format_unop_latex
        if self._exterior_derivative is None:
           # A new computation is necessary:
            if chart_name is None:
                chart_name = self.pick_a_chart()
            chart = self.manifold.atlas[chart_name]
            frame_name = chart.frame.name
            
            n = self.manifold.dim
            si = self.manifold.sindex
            f = self.express[chart.name]
            dc = Components(self.manifold, 1, frame_name)
            for i in range(n):
                df = f.diff(chart.xx[i])
                dc[i+si] = df
            dc._del_zeros()
            # Name and LaTeX name of the result (rname and rlname):
            rname = format_unop_txt('d', self.name)
            rlname = format_unop_latex(r'\mathrm{d}', self.latex_name)
            self._exterior_derivative = OneForm(self.manifold, rname, rlname)
            self._exterior_derivative.components[frame_name] = dc
        return self._exterior_derivative

    def gradient(self, chart_name=None):
        r"""
        Return the gradient 1-form of the scalar field. 
        
        This method simply calls the method :meth:`exterior_der`.  
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart used for the
          computation; if none is provided, the manifold's default chart is 
          used.
        
        OUTPUT:
        
        - the 1-form gradient of ``self``. 
        
        EXAMPLES:
        
        Gradient on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz = Chart(m, 'x y z', 'coord_xyz')
            sage: f = ScalarField(m, cos(x)*z**3 + exp(y)*z**2, name='f')
            sage: df = f.gradient() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.show()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        """
        return self.exterior_der(chart_name)

    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.
        
        The Lie derivative is stored in the dictionary 
        :attr:`_lie_derivatives`, so that there is no need to 
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile. 
        
        In the present case (scalar field), the Lie derivative is equal to
        the scalar field resulting from the action of the vector field on 
        ``self``. 
        
        INPUT:
        
        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken
          
        OUTPUT:
        
        - the scalar field that is the Lie derivative of ``self`` with 
          respect to ``vector``
          
        EXAMPLES:
        
        Lie derivative on a 2-dimensional manifold::

            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: f = ScalarField(m, x^2*cos(y))
            sage: v = VectorField(m, 'v')
            sage: v[:] = (-y, x)
            sage: f.lie_der(v)
            scalar field on the 2-dimensional manifold 'M'
            sage: f.lie_der(v).expr()
            -x^3*sin(y) - 2*x*y*cos(y)

        Alternative expressions of the Lie derivative of a scalar field::
        
            sage: f.lie_der(v) == v(f)  # the vector acting on f
            True
            sage: f.lie_der(v) == f.gradient()(v)  # the gradient of f acting on the vector
            True

        A vanishing Lie derivative::
        
            sage: f.set_expr(x^2 + y^2)
            sage: f.lie_der(v)
            zero scalar field on the 2-dimensional manifold 'M'

        """
        from vectorfield import VectorField
        from vectorframe import CoordBasis
        if not isinstance(vector, VectorField):
            raise TypeError("The argument must be a vector field.")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            res = vector(self)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]         

#******************************************************************************

class ZeroScalarField(ScalarField):
    r"""
    Null scalar field on a differentiable manifold.
    
    INPUT:
    
    - ``manifold`` -- the manifold on which the scalar field is defined
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    Zero scalar field on a 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')                  
        sage: c_xy = Chart(m, 'x y', 'coord_xy')
        sage: f = ZeroScalarField(m) ; f
        zero scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        0
        sage: f.is_zero()
        True
        sage: p = Point(m, (1,2))
        sage: f(p)
        0

    Each manifold has a predefined zero scalar field::
    
        sage: m.zero_scalar_field
        zero scalar field on the 2-dimensional manifold 'M'
        sage: m.zero_scalar_field(p)
        0
        sage: f == m.zero_scalar_field
        True

    Arithmetics with another instance of :class:`ZeroScalarField`::
    
        sage: h = ZeroScalarField(m)
        sage: s = f+h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f-h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f*h ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f/h ; s ; s.expr()
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a scalar field by zero.
        
    Arithmetics with a non-zero instance of :class:`ScalarField`::
    
        sage: g = ScalarField(m, x+y)
        sage: s = f+g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = g+f ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f-g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        -x - y
        sage: s = g-f ; s ; s.expr()                     
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f*g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g*f ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f/g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g/f ; s ; s.expr()
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a scalar field by zero.

    Arithmetics with a non-zero instance of :class:`FunctionChart`::

        sage: g = FunctionChart(c_xy, x+y)
        sage: s = f+g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = g+f ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f-g ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        -x - y
        sage: s = g-f ; s ; s.expr()
        scalar field on the 2-dimensional manifold 'M'
        x + y
        sage: s = f*g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g*f ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = f/g ; s ; s.expr()
        zero scalar field on the 2-dimensional manifold 'M'
        0
        sage: s = g/f ; s ; s.expr()
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a FunctionChart by zero.
        
    Arithmetics with an instance of :class:`ZeroFunctionChart`::

        sage: g = ZeroFunctionChart(c_xy)
        sage: s = f+g ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = g+f ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = f-g ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = g-f ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = f*g ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = g*f ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: s = f/g ; s
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a scalar field by zero.

    """
    def __init__(self, manifold, name=None, latex_name=None):
        from manifold import RealLine
        DiffMapping.__init__(self, manifold, RealLine)
        DiffForm.__init__(self, manifold, 0, name, latex_name)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "zero scalar field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.manifold)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`ZeroScalarField` instance on the same manifold.
        
        """
        return ZeroScalarField(self.manifold)        

    def copy(self):
        r"""
        Return an exact copy of ``self``.
        """
        return ZeroScalarField(self.manifold)
        
    def function_chart(self, chart_name=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart; if none, the 
          manifold's default chart will be used
          
        OUTPUT:
        
        - instance of :class:`ZeroFunctionChart` defined in the specified chart.
        
        """
        if chart_name is None:
            chart = self.manifold.def_chart
        else:
            chart = self.manifold.atlas[chart_name]
        return chart.zero_function
 
    def expr(self, chart_name=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart_name`` -- (default: None) unused here
        - ``from_chart`` -- (default: None) unused here
                  
        OUTPUT:
        
        - number zero
        
        """
        return 0
 
    def set_expr(self, coord_expression, chart_name=None, delete_others=True):
        r"""
        Set some coordinate expression of the scalar field.
        
        Not valid for a :class:`ZeroScalarField` object. 
        """
        raise TypeError("set_expr has no meaning for a zero scalar field.")

    def show(self, chart_name=None):
        r""" 
        Displays the expression of the scalar field. 
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart for the  
        coordinate expression of the scalar field; unused here,
        since the scalar field is identically zero. 
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if self.name is None:
            result.txt = "0"
        else:
            result.txt = self.name + " = 0"
        if self.latex_name is None:
            result.latex = "0"
        else:
            result.latex = latex(self) + " = 0"
        return result

    def __call__(self, p):
        r"""
        Computes the image of a point.

        INPUT:
    
        - ``p`` -- point on the manifold self.manifold1 (type: :class:`Point`)
        
        OUTPUT:

        - the number zero. 
        
        """
        from point import Point
        if not isinstance(p, Point):
            return TypeError("The argument must be a point.")
        return 0

    def common_chart(self, other):
        r"""
        Find a common chart for the expressions of ``self`` and ``other``. 
        
        This is a redefinition of :meth:`ScalarField.common_chart` which 
        simply returns one of the charts of ``other``, privileging the 
        manifold's default chart if present. 
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUPUT:
        
        - see above. 

        """
        # Compatibility checks:
        if not isinstance(other, ScalarField):
            raise TypeError("The second argument must be a scalar field.")
        if other.manifold != self.manifold:
            raise TypeError("The two scalar fields are not defined on the " +
                            "same manifold.")
        def_chart_name = self.manifold.def_chart.name
        if def_chart_name in other.express:
            return def_chart_name
        else:
            return other.express.keys()[0]
   
    def is_zero(self):
        r""" 
        Return True if the scalar field is zero and False otherwise.

        """
        return True
            
    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a scalar field or something else
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if isinstance(other, ScalarField):
            return other.is_zero()
        else:
            return bool(isinstance(other, (int, Integer)) and other==0)
        #!# what about __req___ ?
 
    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - ``self``
    
        """
        return self

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - ``self`` (since ``self`` is zero)
    
        """
        return self

    def __add__(self, other):
        r"""
        Scalar field addition. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the addition of ``self`` and 
          ``other``
        
        """
        if isinstance(other, ScalarField):
            if other.manifold != self.manifold:
                raise TypeError("Scalar fields defined on different " + 
                                "manifolds cannot be added.")
            return other.copy()    
        else:
            return scalar_field(other, self.manifold)
            
    def __sub__(self, other):
        r"""
        Scalar field subtraction. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the subtraction of ``other`` from 
          ``self``          
        
        """
        if isinstance(other, ScalarField):
            if other.manifold != self.manifold:
                raise TypeError("Scalar fields defined on different " + 
                                "manifolds cannot be added.")
            return -other    
        else:
            return scalar_field(-other, self.manifold)

    def __mul__(self, other):
        r"""
        Scalar field multiplication. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the multiplication of ``self`` by 
          ``other``
        
        """
        from component import Components
        from tensorfield import TensorField
        if isinstance(other, ScalarField):
            return self
        elif isinstance(other, TensorField):
            result = other._new_instance()
            for frame_name in other.components:
                result.components[frame_name] = self * other.components[frame_name]
            return result
        elif isinstance(other, Components):
            return other._new_instance()  # a just created Components is zero
        else: # other is not a tensor field:
            return self

    def __div__(self, other):
        r"""
        Scalar field division. 
        
        INPUT:
        
        - ``other`` -- a scalar field or an expression
        
        OUPUT:
        
        - the scalar field resulting from the division of ``self`` by 
          ``other``
        
        """
        if other == 0:
            raise ZeroDivisionError("Division of a scalar field by zero.")
        else:
            return self

    def exterior_der(self, chart_name=None):
        r"""
        Return the exterior derivative of the scalar field, which is zero in 
        the present case. 
        
        INPUT:
        
        - ``chart_name`` -- (default: None) name of the chart used for the
          computation.
        
        OUTPUT:
        
        - the 1-form exterior derivative of ``self``. 
                
        """
        from component import Components
        if self._exterior_derivative is None:
            # A new computation is necessary:
            if chart_name is None:
                chart_name = self.pick_a_chart()
            chart = self.manifold.atlas[chart_name]
            frame_name = chart.frame.name
            
            n = self.manifold.dim
            si = self.manifold.sindex
            dc = Components(self.manifold, 1, frame_name)
            for i in range(n):      #!# Not necessary if the Components are  
                dc[i+si] = 0        #!# initialized to zero
            
            if self.name is None:
                name_r = None
            else:
                name_r = 'd' + self.name
            if self.latex_name is None:
                latex_name_r = None
            else:
                latex_name_r = r'\mathrm{d}' + self.latex_name
            self._exterior_derivative = OneForm(self.manifold, name_r, latex_name_r)
            self._exterior_derivative.components[frame_name] = dc
        return self._exterior_derivative

 
#******************************************************************************
    
def scalar_field(other, manifold=None, name=None, latex_name=None):
    r"""
    Construct a scalar field from another type of object.
    
    This is a handmade conversion function. 
    
    INPUT:
    
    - ``other`` -- an instance of :class:`FunctionChart` or a symbolic 
      expression 
    - ``manifold`` -- (default: None) a differentiable manifold on which the 
      scalar is to be defined; if none is provided, this information should be
      readable on ``other``
    - ``name`` -- (default: None) name given to the scalar field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the scalar field; 
      if none is provided, the LaTeX symbol is set to ``name``
    
    NB: if ``other`` is a symbolic expression, it would be preferable to invoke
    the :class:`ScalarField` constructor instead.

    OUTPUT:
    
    - an instance of :class:`ScalarField` or of :class:`ZeroScalarField` if
      ``other==0``. 
    
    EXAMPLES:
    
    Conversion on a 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')                  
        sage: c_xy = Chart(m, 'x y', 'coord_xy')
        sage: fc = FunctionChart(c_xy, x+2*y^3)
        sage: f = scalar_field(fc) ; f
        scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        2*y^3 + x
        sage: g = scalar_field(x*y, m) ; g
        scalar field on the 2-dimensional manifold 'M'
        sage: g.expr()
        x*y

    Conversion on a chart that is not the manifold's default one::
    
        sage: c_uv = Chart(m, 'u v', 'coord_uv') # a new chart on M
        sage: m.default_chart()  # 'coord_uv' is not the manifold's default chart:
        chart 'coord_xy' (x, y)
        sage: fc1 = FunctionChart(c_uv, u+v)
        sage: f1 = scalar_field(fc1) ; f1
        scalar field on the 2-dimensional manifold 'M'
        sage: f1.expr('coord_uv')
        u + v

    Conversion of a null object to an instance of :class:`ZeroScalarField`::
        
        sage: f = scalar_field(0, m) ; f
        zero scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        0
        sage: fc = ZeroFunctionChart(c_xy)
        sage: f = scalar_field(fc) ; f
        zero scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        0
        sage: fc = FunctionChart(c_xy, 0)
        sage: f = scalar_field(fc) ; f
        zero scalar field on the 2-dimensional manifold 'M'
        sage: f.expr()
        0
        
    """
    if isinstance(other, FunctionChart):
        if isinstance(other, ZeroFunctionChart) or other.is_zero():
            return other.chart.manifold.zero_scalar_field
        else:
            result = ScalarField(other.chart.manifold, name=name, 
                                 latex_name=latex_name)
            result.express = {other.chart.name: other}
            result.coord_expression = {(other.chart.name,'canonical'):
                                MultiFunctionChart(other.chart, other.express)}
            return result
    else:
        if other == 0:
            return manifold.zero_scalar_field
        else:
            return ScalarField(manifold, other, manifold.def_chart.name, 
                               name=name, latex_name=latex_name)
