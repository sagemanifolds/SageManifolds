r"""
Scalar fields

The class :class:`ScalarField` implements scalar fields on differentiable 
manifolds over `\RR`, i.e. differentiable mappings of the form

.. MATH::

    f: U\subset \mathcal{M} \longrightarrow \RR
    
where `U` is a domain of the differentiable manifold `\mathcal{M}`.

The class :class:`ScalarField`  inherits from the classes 
:class:`DiffMapping` (a scalar field being a differentiable mapping to `\RR`) 
and :class:`DiffForm` (a scalar field being a differential form of degree 0). 

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
from domain import Domain
from chart import FunctionChart, ZeroFunctionChart, MultiFunctionChart
from diffmapping import DiffMapping
from diffform import DiffForm, OneForm

class ScalarField(DiffMapping, DiffForm):
    r"""
    Class for scalar fields on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain `U` on which the scalar field is 
      defined (must be an instance of class :class:`Domain`)
    - ``coord_expression`` -- (default: None) coordinate expression of the 
      scalar field
    - ``chart`` -- (default:None) chart defining the coordinates used in 
      ``coord_expression``; if none is provided and a
      coordinate expression is given, the domain default chart is assumed.
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
    
        sage: c_spher.<th,ph> = m.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: f = ScalarField(m, sin(th)*cos(ph), name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'S^2'

    The coordinate expression of a scalar field can be read by means of 
    :meth:`expr` and set by means of :meth:`set_expr`; both methods can 
    have a chart as argument (if not, the manifold's default chart is 
    assumed)::
    
        sage: f.expr()
        cos(ph)*sin(th)
        sage: f.expr(c_spher) # equivalent to above since c_spher is the default chart
        cos(ph)*sin(th)
        sage: f.set_expr(cos(th))  # changing the value of f
        sage: f.expr()
        cos(th)
        sage: f.set_expr(sin(th)*cos(ph)) # restoring the original value
        
    The function :meth:`view` displays the coordinate expression of the scalar
    field::
    
        sage: f.view()
        f: (th, ph) |--> cos(ph)*sin(th)
        sage: f.view(c_spher) # equivalent to above since c_spher is the default chart
        f: (th, ph) |-->  cos(ph)*sin(th)
        sage: latex(f.view()) # nice LaTeX formatting for the notebook
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
    
        sage: f.function_chart(c_spher)
        cos(ph)*sin(th)
        sage: f.function_chart() # equivalent to above since c_spher is the default chart
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
        sage: f.domain1 # the domain on the start manifold
        2-dimensional manifold 'S^2'
        sage: f.domain2 # the target domain
        field R of real numbers
        sage: print type(f.domain2)
        <class 'sage.geometry.manifolds.manifold.RealLineManifold_with_category'>
        sage: f.domain2 is RealLine
        True
        sage: f.coord_expression
        {(chart (S^2, (th, ph)), chart (field R, (x_realline,))): functions (cos(ph)*sin(th),) on the chart (S^2, (th, ph))}
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
        
    Test of the arithmetics of scalar fields defined on multiple domains::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U')
        sage: c_xy.<x,y> = U.chart('x y')
        sage: V = M.open_domain('V')
        sage: c_uv.<u,v> = V.chart('u v')
        sage: f = ScalarField(M, x^2)
        sage: f.add_expr(u, c_uv)
        sage: g = ScalarField(M, 2*v, c_uv)
        sage: g.add_expr(y, c_xy)
        sage: f.express
        {chart (U, (x, y)): x^2, chart (V, (u, v)): u}
        sage: g.express
        {chart (V, (u, v)): 2*v, chart (U, (x, y)): y}
        sage: s = f + g ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s.express
        {chart (U, (x, y)): x^2 + y, chart (V, (u, v)): u + 2*v}
        sage: g.set_expr(3*x, c_xy)
        sage: g.express
        {chart (U, (x, y)): 3*x}
        sage: s = f + g ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s.express
        {chart (U, (x, y)): x^2 + 3*x}
        sage: g = ScalarField(U, 3*x)
        sage: g.express
        {chart (U, (x, y)): 3*x}
        sage: s = f + g ; s
        scalar field on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: s.express
        {chart (U, (x, y)): x^2 + 3*x}
    
    Vanishing result::
    
        sage: g = ScalarField(M, -x^2)
        sage: g.add_expr(-u, c_uv)
        sage: s = f + g ; s
        zero scalar field on the 2-dimensional manifold 'M'
        sage: print type(s)
        <class 'sage.geometry.manifolds.scalarfield.ZeroScalarField'>

    """
    def __init__(self, domain, coord_expression=None, chart=None, name=None, 
                 latex_name=None):
        from manifold import RealLine
        DiffMapping.__init__(self, domain, RealLine, coord_expression, chart)
        DiffForm.__init__(self, domain, 0, name, latex_name)
        if coord_expression is not None:
            if chart is None:
                chart = self.domain.def_chart
            if coord_expression == 0:
                self.express = {chart: chart.zero_function}
            else:
                self.express = {chart: FunctionChart(chart, coord_expression)}
        else:
            self.express = {}


    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "scalar field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`ScalarField` instance on the same domain.
        
        """
        return ScalarField(self.domain)        

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
            sage: c_xy.<x,y> = m.chart('x y')
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
        result = ScalarField(self.domain)  #!# what about the name ?
        for chart, funct in self.express.items():
            result.express[chart] = funct.copy()
        for key, mfunct in self.coord_expression.items():
            result.coord_expression[key] = mfunct.copy()
        return result

    def comp(self, frame=None, from_frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.comp`. 
        
        This method should not be called, since a scalar field has no 
        component w.r.t. a vector frame! 
        It therefore returns None. 
        
        """
        return None

    def set_comp(self, frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.set_comp`.
        
        This method should not be called, since a scalar field has no 
        component w.r.t. a vector frame! 
        It therefore returns None. 
                
        """
        return None

    def add_comp(self, frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.add_comp`.
        
        This method should not be called, since a scalar field has no 
        component w.r.t. a vector frame! 
        It therefore returns None. 
                
        """
        return None

    def function_chart(self, chart=None, from_chart=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the
          coordinate expression is to be returned; if None, the 
          domain's default chart will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart``; if None, a chart is picked in ``self.express``

        OUTPUT:
        
        - instance of :class:`FunctionChart` representing the coordinate 
          function of the scalar field in the given chart.

        EXAMPLES:
        
        Coordinate function on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')            
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = ScalarField(m, x*y^2)
            sage: f.function_chart()
            x*y^2
            sage: f.function_chart(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: print type(f.function_chart())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>

        Expression via a change of coordinates::
        
            sage: c_uv.<u,v> = m.chart('u v')
            sage: CoordChange(c_uv, c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f.express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.function_chart(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: f.function_chart(c_uv) == (u+v)*(u-v)^2  # check
            True
            sage: f.express  # f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        Usage in a physical context (simple Lorentz transformation - boost in 
        x direction, with relative velocity v between o1 and o2 frames)::

            sage: m = Manifold(2, 'M')
            sage: o1.<t,x> = m.chart('t x')
            sage: o2.<T,X> = m.chart('T X')
            sage: f = ScalarField(m, x^2 - t^2)
            sage: f.function_chart(o1)
            -t^2 + x^2
            sage: v = var('v'); gam = 1/sqrt(1-v^2)
            sage: CoordChange(o2, o1, gam*(T - v*X), gam*(X - v*T))
            coordinate change from chart (M, (T, X)) to chart (M, (t, x))
            sage: f.function_chart(o2)
            -T^2 + X^2

        """
        if chart is None:
            chart = self.domain.def_chart
        else:
            if chart not in self.domain.atlas:
                raise TypeError("The " + str(chart) + " has not " + \
                      " been defined on the domain " + str(self.domain))
        if chart not in self.express:
            # Check whether chart corresponds to a subchart of a chart
            # where the expression of self is known:
            for known_chart in self.express:
                if chart in known_chart.subcharts:
                    new_expr = self.express[known_chart].expr()
                    self.express[chart] = FunctionChart(chart, new_expr)
                    self.coord_expression[(chart, self.domain2.def_chart)] = \
                                            MultiFunctionChart(chart, new_expr)
                    self._del_derived()
                    return self.express[chart]
            # If this point is reached, the expression must be computed 
            # from that in the chart from_chart, by means of a 
            # change-of-coordinates formula:
            if from_chart is None:
                for known_chart in self.express:
                    if (chart, known_chart) in self.domain.coord_changes:
                        from_chart = known_chart
                        break
                if from_chart is None:
                    raise ValueError("No starting chart could be found to " + 
                           "compute the expression in the " + str(chart))
            change = self.domain.coord_changes[(chart, from_chart)]
            # old coordinates expressed in terms of the new ones:
            coords = [ change.transf.functions[i].express 
                       for i in range(self.manifold.dim) ]
            new_expr = self.express[from_chart](*coords)
            self.express[chart] = FunctionChart(chart, new_expr)
            self.coord_expression[(chart, self.domain2.def_chart)] = \
                                        MultiFunctionChart(chart, new_expr)
            self._del_derived()
        return self.express[chart]


    def expr(self, chart=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the 
          coordinate expression is required; if None, the domain's default 
          chart will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the 
          chart ``chart``; if None, a chart is picked in ``self.express``
          
        OUTPUT:
        
        - symbolic expression representing the coordinate 
          expression of the scalar field in the given chart.
        
        EXAMPLES:
        
        Expression of a scalar field on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')            
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = ScalarField(m, x*y^2)
            sage: f.expr()
            x*y^2
            sage: f.expr(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: print type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        Expression via a change of coordinates::
        
            sage: c_uv.<u,v> = m.chart('u v')
            sage: CoordChange(c_uv, c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f.express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.expr(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: bool( f.expr(c_uv) == (u+v)*(u-v)^2 ) # check
            True
            sage: f.express  # f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        """
        return self.function_chart(chart, from_chart).express
        
    def set_expr(self, coord_expression, chart=None):
        r"""
        Set the coordinate expression of the scalar field.
        
        The expressions with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.
        
        INPUT:
        
        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression`` is 
          defined; if None, the domain's default chart is assumed
        
        EXAMPLES:
        
        Setting scalar field expressions on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = ScalarField(m, x^2 + 2*x*y +1)
            sage: f.express                         
            {chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.set_expr(3*y)
            sage: f.express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = m.chart('u v')
            sage: f.set_expr(cos(u)-sin(v), c_uv)  
            sage: f.express # the (x,y) expression has been lost:
            {chart (M, (u, v)): cos(u) - sin(v)}
            sage: f.set_expr(3*y)    
            sage: f.express # the (u,v) expression has been lost:                    
            {chart (M, (x, y)): 3*y}

        """
        if chart is None:
            chart = self.domain.def_chart
        self.express.clear()
        self.coord_expression.clear()
        self.express[chart] = FunctionChart(chart, coord_expression)
        self.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(chart, coord_expression)
        self._del_derived()

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.
        
        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead. 
        
        INPUT:
        
        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression``
          is defined; if None, the domain's default chart is assumed
          
        .. WARNING::
        
            If the scalar field has already expressions in other charts, it 
            is the user's responsability to make sure that the expression
            to be added is consistent with them. 
        
        EXAMPLES:
        
        Adding scalar field expressions on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = ScalarField(m, x^2 + 2*x*y +1)
            sage: f.express
            {chart (M, (x, y)): x^2 + 2*x*y + 1}             
            sage: f.add_expr(3*y)
            sage: f.express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = m.chart('u v')
            sage: f.add_expr(cos(u)-sin(v), c_uv)  
            sage: f.express # f has now 2 expressions:
            {chart (M, (x, y)): 3*y, chart (M, (u, v)): cos(u) - sin(v)}

        """
        if chart is None:
            chart = self.domain.def_chart
        self.express[chart] = FunctionChart(chart, coord_expression)
        self.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(chart, coord_expression)
        self._del_derived()


    def view(self, chart=None):
        r""" 
        Display the expression of the scalar field in a given chart. 
        
        INPUT:
        
        - ``chart`` -- (default: None) chart with respect to which the 
          coordinate expression is to be displayed; if None, the domain's 
          default chart will be used
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Various displays::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = ScalarField(m, sqrt(x+1), name='f')
            sage: f.view()
            f: (x, y) |--> sqrt(x + 1)
            sage: latex(f.view())
            f :\ \left(x, y\right) \mapsto \sqrt{x + 1}
            sage: g = ScalarField(m, function('G', x, y), name='g')
            sage: g.view() 
            g: (x, y) |--> G(x, y)
            sage: latex(g.view())
            g :\ \left(x, y\right) \mapsto G\left(x, y\right)

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if chart is None:
            chart = self.domain.def_chart
        expression = self.expr(chart)
        if self.name is None:
            result.txt = repr(chart[:]) + " |--> " + repr(expression)
        else:
            result.txt = self.name + ": " + repr(chart[:]) + " |--> " + \
                         repr(expression)
        if self.latex_name is None:
            result.latex = latex(chart[:]) + r"\mapsto" + latex(expression)
        else:
            result.latex = latex(self) + ":\ " + latex(chart[:]) + r"\mapsto" + \
                           latex(expression)
        return result


    def pick_a_chart(self):
        r"""
        Return a chart for which the scalar field has an expression. 
        
        The domain's default chart is privileged. 
        
        OUPUT:
        
        - name of the chart

        EXAMPLES:
        
        A very simple example::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: c_uv.<u,v> = m.chart('u v')
            sage: f =  ScalarField(m, x*y^2)
            sage: f.add_expr(u-v, c_uv)
            sage: f.express # f has expressions on two charts:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u - v}
            sage: m.default_chart()
            chart (M, (x, y))
            sage: f.pick_a_chart()  # the domain's default chart (xy-coord) is privileged:
            chart (M, (x, y))
            sage: g = ScalarField(m, u+v, c_uv)
            sage: g.express  # g has no expression on the domain's default chart:
            {chart (M, (u, v)): u + v}
            sage: g.pick_a_chart()
            chart (M, (u, v))
        
        """
        if self.domain.def_chart in self.express:
            return self.domain.def_chart
        return self.express.items()[0][0]  
            

    def common_charts(self, other):
        r"""
        Find common charts for the expressions of ``self`` and ``other``. 
        
        INPUT:
        
        - ``other`` -- a scalar field
        
        OUPUT:
        
        - list of common charts; if no common chart is found, None is 
          returned (instead of an empty list). 

        EXAMPLES:
        
        Search for common charts on a 2-dimensional manifold with 2 
        overlapping domains::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U')
            sage: c_xy.<x,y> = U.chart('x y')
            sage: V = M.open_domain('V')
            sage: c_uv.<u,v> = V.chart('u v')
            sage: f = ScalarField(U, x^2)
            sage: g = ScalarField(M, x+y)
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: g.add_expr(u, c_uv)
            sage: f.express
            {chart (U, (x, y)): x^2}
            sage: g.express
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y))]

        Common charts found as subcharts: the subcharts are introduced via
        a transition map between charts c_xy and c_uv on the intersecting domain
        `W = U\cap V`::
        
            sage: trans = c_xy.transition_map(c_uv, (x+y, x-y), 'W', x<0, u>0)
            sage: c_xy_W = trans.chart1
            sage: c_uv_W = trans.chart2
            sage: trans.inverse()
            coordinate change from chart (W, (u, v)) to chart (W, (x, y))
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: f.expr(c_xy_W)  
            x^2
            sage: f.express
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2}
            sage: g.express
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: g.common_charts(f)  # c_xy_W is not returned because it is subchart of 'xy'
            [chart (U, (x, y))]
            sage: f.expr(c_uv_W)
            1/4*u^2 + 1/2*u*v + 1/4*v^2
            sage: f.express 
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2, chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2}
            sage: g.express
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y)), chart (W, (u, v))]
            sage: g.express # the expressions have been updated on the subcharts
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u, chart (W, (u, v)): u}

        Common charts found by computing some coordinate changes::
        
            sage: W = M.domains['W']
            sage: f = ScalarField(W, x^2, c_xy_W)
            sage: g = ScalarField(W, u+1, c_uv_W)
            sage: f.express
            {chart (W, (x, y)): x^2}
            sage: g.express
            {chart (W, (u, v)): u + 1}
            sage: f.common_charts(g)
            [chart (W, (u, v)), chart (W, (x, y))]
            sage: f.express
            {chart (W, (x, y)): x^2, chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2}
            sage: g.express
            {chart (W, (u, v)): u + 1, chart (W, (x, y)): x + y + 1}

        """
        if not isinstance(other, ScalarField):
            raise TypeError("The second argument must be a scalar field.")
        dom1 = self.domain
        dom2 = other.domain
        coord_changes = self.manifold.coord_changes
        resu = []
        #
        # 1/ Search for common chart among the existing expressions, i.e. 
        #    without performing any expression transformation. 
        #    -------------------------------------------------------------
        for chart1 in self.express:
            if chart1 in other.express:
                resu.append(chart1)
        # Search for a subchart:
        known_expr1 = self.express.copy()  
        known_expr2 = other.express.copy()
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if chart2 in chart1.subcharts:
                            self.expr(chart2)
                            resu.append(chart2)
                        if chart1 in chart2.subcharts:
                            other.expr(chart1)
                            resu.append(chart1)
        #
        # 2/ Search for a common chart via one expression transformation
        #    -----------------------------------------------------------
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if (chart1, chart2) in coord_changes:
                            self.function_chart(chart2, from_chart=chart1)
                            resu.append(chart2)
                        if (chart2, chart1) in coord_changes:
                            other.function_chart(chart1, from_chart=chart2)
                            resu.append(chart1)
        if resu == []:
            return None
        else:
            return resu


    def is_zero(self):
        r""" 
        Return True if the scalar field is zero and False otherwise.

        EXAMPLES:
        
        Tests on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
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
        res = True
        for funct in self.express.values():
            res = res and funct.is_zero()
        return res

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
            if other.domain != self.domain:
                return False
            com_charts = self.common_charts(other)
            if com_charts is None:
                raise ValueError("No common chart for the comparison.")
            resu = True
            for chart in com_charts:
                resu = resu and (self.express[chart] == other.express[chart])
            return resu
        elif isinstance(other, FunctionChart):
            chart = other.chart
            if chart not in self.express:
                return False
            else:
                return bool(self.express[chart] == other.express)
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
        for chart in self.express:
            res = + self.express[chart]
            result.express[chart] = res
            result.coord_expression[(chart, self.domain2.def_chart)] = \
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
        for chart in self.express:
            res = - self.express[chart]
            result.express[chart] = res
            result.coord_expression[(chart, self.domain2.def_chart)] = \
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
            com_charts = self.common_charts(other)
            if com_charts is None:
                raise ValueError("No common chart for the addition.")
            dom_result = self.domain.intersection(other.domain)
            result = ScalarField(dom_result)
            for chart in com_charts:
                # FunctionChart addition:
                res = self.express[chart] + other.express[chart]
                result.express[chart] = res
                result.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(res.chart, res.express)
            if result.is_zero():
                return dom_result.zero_scalar_field
            if self.name is not None and other.name is not None:
                result.name = self.name + '+' + other.name
            if self.latex_name is not None and other.latex_name is not None:
                result.latex_name = self.latex_name + '+' + other.latex_name
            return result
        elif isinstance(other, FunctionChart):
            return self + other.scalar_field()
        elif other == 0:
            return self.copy()
        else:
            return self + ScalarField(self.domain, coord_expression=other)

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
            com_charts = self.common_charts(other)
            if com_charts is None:
                raise ValueError("No common chart for the subtraction.")
            dom_result = self.domain.intersection(other.domain)
            result = ScalarField(dom_result)
            for chart in com_charts:
                # FunctionChart subtraction:
                res = self.express[chart] - other.express[chart]
                result.express[chart] = res
                result.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(res.chart, res.express)
            if result.is_zero():
                return dom_result.zero_scalar_field
            if self.name is not None and other.name is not None:
                result.name = self.name + '-' + other.name
            if self.latex_name is not None and other.latex_name is not None:
                result.latex_name = self.latex_name + '-' + other.latex_name
            return result
        elif isinstance(other, FunctionChart):
            return self - other.scalar_field()
        elif other == 0:
            return self.copy()
        else:
            return self - ScalarField(self.domain, coord_expression=other)

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
            com_charts = self.common_charts(other)
            if com_charts is None:
                raise ValueError("No common chart for the multiplication.")
            dom_result = self.domain.intersection(other.domain)
            result = ScalarField(dom_result)
            for chart in com_charts:
                # FunctionChart multiplication:
                res = self.express[chart] * other.express[chart]
                result.express[chart] = res
                result.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(res.chart, res.express)
            if result.is_zero():
                return dom_result.zero_scalar_field
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, ' ', 
                                                 other.latex_name)        

        elif isinstance(other, TensorField):
            #!# subdomains almost treated here 
            result = other._new_instance()
            dom1 = self.domain
            for frame in other.components:
                comp = other.components[frame]
                dom2 = comp.domain
                if dom1 in dom2.subdomains or dom2 in dom1.subdomains:
                    try:
                        result.components[frame] = self * comp
                    except ValueError:
                        pass
            result.name = format_mul_txt(self.name, '*', other.name)
            result.latex_name = format_mul_latex(self.latex_name, ' ', 
                                                 other.latex_name)

        elif isinstance(other, Components):
            #!# subdomains not treated here 
            result = other._new_instance()
            for ind in other._comp:
                result._comp[ind] = self * other._comp[ind]

        elif isinstance(other, FunctionChart):
            return self * other.scalar_field()
        elif other == 0:
            return self.domain.zero_scalar_field
        else:
            return self * ScalarField(self.domain, coord_expression=other)

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
            com_charts = self.common_charts(other)
            if com_charts is None:
                raise ValueError("No common chart for the division.")
            dom_result = self.domain.intersection(other.domain)
            result = ScalarField(dom_result)
            for chart in com_charts:
                # FunctionChart division:
                res = self.express[chart] / other.express[chart]
                result.express[chart] = res
                result.coord_expression[(chart, self.domain2.def_chart)] = \
                                    MultiFunctionChart(res.chart, res.express)
            #!# the following 2 lines could be skipped:
            if result.is_zero():
                return dom_result.zero_scalar_field
            result.name = format_mul_txt(self.name, '/', other.name)
            result.latex_name = format_mul_latex(self.latex_name, '/', 
                                                 other.latex_name)
            return result

        elif isinstance(other, FunctionChart):
            return self / other.scalar_field()
        elif other == 0:
            raise ZeroDivisionError("Division of a scalar field by zero.")
        else:
            return self / ScalarField(self.domain, coord_expression=other)


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

    def exterior_der(self, chart=None):
        r"""
        Return the exterior derivative of the scalar field. 
        
        INPUT:
        
        - ``chart`` -- (default: None) chart used for the computation; if none 
          is provided, a chart on which the scalar field is expressed is 
          picked, privileging the domain's defaut chart. 
        
        OUTPUT:
        
        - the 1-form exterior derivative of ``self``. 
        
        EXAMPLES:
        
        Exterior derivative on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = m.chart('x y z')
            sage: f = ScalarField(m, cos(x)*z^3 + exp(y)*z^2, name='f')
            sage: df = f.exterior_der() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.view()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        Exterior derivative computed on a chart that is not the default one::
        
            sage: c_uvw.<u,v,w> = m.chart('u v w')
            sage: g = ScalarField(m, u*v^2*w^3, c_uvw, name='g')
            sage: dg = g.exterior_der() ; dg
            1-form 'dg' on the 3-dimensional manifold 'M'
            sage: dg.components
            {coordinate frame (M, (d/du,d/dv,d/dw)): 1-index components w.r.t. the coordinate frame (M, (d/du,d/dv,d/dw))}
            sage: dg.comp(c_uvw.frame)[:, c_uvw]
            [v^2*w^3, 2*u*v*w^3, 3*u*v^2*w^2]
            sage: dg.view(c_uvw.frame, c_uvw)
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
            if chart is None:
                chart = self.pick_a_chart()
            n = self.manifold.dim
            si = self.manifold.sindex
            f = self.express[chart]
            dc = Components(chart.frame, 1)
            for i in range(n):
                df = f.diff(chart.xx[i])
                dc[i+si] = df
            dc._del_zeros()
            # Name and LaTeX name of the result (rname and rlname):
            rname = format_unop_txt('d', self.name)
            rlname = format_unop_latex(r'\mathrm{d}', self.latex_name)
            self._exterior_derivative = OneForm(self.domain, rname, rlname)
            self._exterior_derivative.components[chart.frame] = dc
        return self._exterior_derivative

    def differential(self, chart=None):
        r"""
        Return the differential of the scalar field. 
        
        This method simply calls the method :meth:`exterior_der`.  
        
        INPUT:
        
        - ``chart`` -- (default: None) chart used for the computation; if none 
          is provided, the domain's default chart is used.
        
        OUTPUT:
        
        - the 1-form that is the differential of ``self``. 
        
        EXAMPLES:
        
        Differential 1-form on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = m.chart('x y z')
            sage: f = ScalarField(m, cos(x)*z**3 + exp(y)*z**2, name='f')
            sage: df = f.differential() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.view()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            
        """
        return self.exterior_der(chart)

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
            sage: c_xy.<x,y> = m.chart('x y')
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
            sage: f.lie_der(v) == f.differential()(v)  # the differential of f acting on the vector
            True

        A vanishing Lie derivative::
        
            sage: f.set_expr(x^2 + y^2)
            sage: f.lie_der(v)
            zero scalar field on the 2-dimensional manifold 'M'

        """
        from vectorfield import VectorField
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
    
    - ``domain`` -- the manifold domain on which the scalar field is defined
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    Zero scalar field on a 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')                  
        sage: c_xy.<x,y> = m.chart('x y')
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
    def __init__(self, domain, name=None, latex_name=None):
        from manifold import RealLine
        DiffMapping.__init__(self, domain, RealLine)
        DiffForm.__init__(self, domain, 0, name, latex_name)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "zero scalar field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`ZeroScalarField` instance on the same domain.
        
        """
        return ZeroScalarField(self.domain)        

    def copy(self):
        r"""
        Return an exact copy of ``self``.
        """
        return ZeroScalarField(self.domain)
        
    def function_chart(self, chart=None):
        r""" 
        Return the function of the coordinates representing the scalar field 
        in a given chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) chart; if None, the domain's default 
          chart will be used
          
        OUTPUT:
        
        - instance of :class:`ZeroFunctionChart` defined in the specified chart.
        
        """
        if chart is None:
            chart = self.domain.def_chart
        return chart.zero_function
 
    def expr(self, chart=None, from_chart=None):
        r""" 
        Return the coordinate expression of the scalar field in a given 
        chart.
        
        INPUT:
        
        - ``chart`` -- (default: None) unused here
        - ``from_chart`` -- (default: None) unused here
                  
        OUTPUT:
        
        - number zero
        
        """
        return 0
 
    def set_expr(self, coord_expression, chart=None):
        r"""
        Set some coordinate expression of the scalar field.
        
        Not valid for a :class:`ZeroScalarField` object. 
        """
        raise TypeError("set_expr() has no meaning for a zero scalar field.")

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.
        
        Not valid for a :class:`ZeroScalarField` object. 
        """
        raise TypeError("add_expr() has no meaning for a zero scalar field.")

    def view(self, chart=None):
        r""" 
        Display the expression of the scalar field. 
        
        INPUT:
        
        - ``chart`` -- (default: None) chart for the  coordinate expression 
          of the scalar field; unused here, since the scalar field is 
          identically zero. 
          
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
    
        - ``p`` -- point on the manifold (type: :class:`Point`)
        
        OUTPUT:

        - the number zero. 
        
        """
        from point import Point
        if not isinstance(p, Point):
            return TypeError("The argument must be a point.")
        return 0
   
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
        elif isinstance(other, FunctionChart):
            return other.scalar_field()
        elif other == 0:
            return self.domain.zero_scalar_field
        else:
            return ScalarField(self.domain, coord_expression=other)
            
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
                                "manifolds cannot be subtracted.")
            return -other    
        elif isinstance(other, FunctionChart):
            return (-other).scalar_field()
        elif other == 0:
            return self.domain.zero_scalar_field
        else:
            return ScalarField(self.domain, coord_expression=-other)

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
            for frame in other.components:
                result.components[frame] = self * other.components[frame]
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

    def exterior_der(self, chart=None):
        r"""
        Return the exterior derivative of the scalar field, which is zero in 
        the present case. 
        
        INPUT:
        
        - ``chart`` -- (default: None) name of the chart used for the
          computation.
        
        OUTPUT:
        
        - the 1-form exterior derivative of ``self``. 
                
        """
        from component import Components
        if self._exterior_derivative is None:
            # A new computation is necessary:
            if chart is None:
                chart = self.pick_a_chart()
            n = self.manifold.dim
            si = self.manifold.sindex
            dc = Components(chart.frame, 1)
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
            self._exterior_derivative = OneForm(self.domain, name_r, latex_name_r)
            self._exterior_derivative.components[chart.frame] = dc
        return self._exterior_derivative

