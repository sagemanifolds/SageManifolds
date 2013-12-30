r"""
Charts on a manifold

Five classes are defined to deal with coordinates on a differentiable manifold
over `\RR`:

* :class:`Chart` for charts on a manifold
* :class:`FunctionChart` for real-valued functions of the coordinates of a given 
  chart
* :class:`ZeroFunctionChart` for the null function of the coordinates of a given 
  chart
* :class:`MultiFunctionChart` for sets of real-valued functions of coordinates 
  of a given chart 
* :class:`CoordChange` for transition maps between charts

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

EXAMPLES:
            
    Defining the chart for spherical coordinates on `\RR^3`::

        sage: m = Manifold(3, 'R3', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
        sage: c_spher
        chart 'spher' (r, th, ph)
        sage: m.get_atlas()
        {'spher': chart 'spher' (r, th, ph)}
        
    The defined coordinates are directly available at the user level::
    
        sage: th
        th
        sage: print type(th)
        <type 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        \theta
        
    Setting a second chart on the manifold, for instance Cartesian coordinates::
    
        sage: c_cart = Chart(m, 'x y z', 'cart')
        sage: c_cart
        chart 'cart' (x, y, z)
        sage: m.get_atlas()
        {'spher': chart 'spher' (r, th, ph), 'cart': chart 'cart' (x, y, z)}
       
    A manifold has a default chart, which, unless otherwise specified, is the 
    first defined chart::
    
        sage: m.default_chart()
        chart 'spher' (r, th, ph)
        
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
from sage.structure.element import RingElement
from sage.rings.integer import Integer
from manifold import Manifold
from utilities import simplify_chain

class Chart(SageObject):
    r"""
    Class for charts on a manifold.

    INPUT:
    
    - ``manifold`` -- the manifold on which the chart is defined
    - ``coord_symbols`` -- string defining the coordinate symbols: the 
      coordinates are separated by ',' and each coordinate has at most three 
      fields, separated by ':': 
        
        1. the coordinate symbol (a letter or a few letters)
        2. (optional) either the keyword 'positive' for a coordinate restricted 
           to `\RR^+` or the LaTeX spelling of the coordinate
        3. (optional) the LaTeX spelling of the coordinate if the second field 
           is 'positive'
        
      If it contains LaTeX expressions, the string coord_symbols must have the 
      prefix 'r' (see examples below). 
        
    - ``name`` -- name given to the chart (should be rather short)
    
    EXAMPLES: 
    
    Spherical coordinates on `\RR^3`::
    
        sage: m = Manifold(3, 'R3', r'\mathcal{M}', start_index=1)
        sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
        sage: c_spher
        chart 'spher' (r, th, ph)

    Cartesian coordinates on `\RR^3`::
    
        sage: c_cart = Chart(m, 'x y z', 'cart')
        sage: c_cart
        chart 'cart' (x, y, z)

    A defined coordinate is directly accessible by its name::

        sage: th
        th
        sage: print type(th)
        <type 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        \theta
        
    It is also accessible by its index::
    
        sage: x1 = c_spher(1); x2 = c_spher(2); x3 = c_spher(3)
        sage: print x1, x2, x3
        r th ph
        sage: (x1, x2, x3) == (r, th, ph)
        True

    If no index is provided, the () function returns all the coordinates::
    
        sage: c_spher()
        (r, th, ph)
        sage: c_cart()
        (x, y, z)
        
    Each constructed chart is automatically added to the manifold's atlas::
    
        sage: m.get_atlas()
        {'spher': chart 'spher' (r, th, ph), 'cart': chart 'cart' (x, y, z)}
        
    Each constructed chart has its zero function, mapping the coordinates to 0;
    this zero function is an instance of :class:`ZeroFunctionChart`::
    
        sage: c_spher.zero_function
        0
        sage: print type(c_spher.zero_function)
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        sage: c_cart.zero_function
        0
        sage: c_cart.zero_function == c_spher.zero_function
        False
        sage: # the result is False for the zero functions are not defined on the same chart
    
    """
    def __init__(self, manifold, coord_symbols, name): 
        from sage.symbolic.ring import SR
        from vectorframe import CoordBasis
        if not isinstance(manifold, Manifold):
            raise TypeError("The first argument must be a manifold.")
        else:
            self.manifold = manifold
        if ',' in coord_symbols: 
            coord_list = coord_symbols.split(',')
        elif ' ' in coord_symbols:
            coord_list = coord_symbols.split()
        else:
            coord_list = [coord_symbols]
        n = manifold.dim 
        if len(coord_list) != n:
            raise ValueError("The list of coordinates must contain " + \
                              str(n) + " elements.")
        # The user-level namespace (because the code is in Cython):
        glob = globals()  
        xx_list = []
        for i in range(n):
            single_coord = coord_list[i].split(':')
            coord_symb = single_coord[0].strip() # the coordinate symbol
            coord_domain = 'real' # default value, possibly redefined below 
            coord_latex = None    # default value, possibly redefined below 
            length = len(single_coord)
            if length == 2:
                if single_coord[1].strip() == 'positive':
                    coord_domain = 'positive'
                else:
                    coord_latex = single_coord[1].strip()
            if length == 3: 
                if single_coord[1].strip() != 'positive':
                    raise ValueError("The second item in coordinate " + 
                                     "descriptor should be 'positive'.")
                coord_domain = 'positive'
                coord_latex = single_coord[2].strip()
            if length > 3: 
                raise ValueError("A coordinate descriptor can have at most " + 
                                 "3 items.")
            coord_var = SR.var(coord_symb, domain=coord_domain, 
                           latex_name=coord_latex)
            xx_list.append(coord_var)
            # the coordinate is added to the dictionary of global variables:
            glob[coord_symb] = coord_var 
        self.xx = tuple(xx_list)
        self.name = name
        manifold.atlas[self.name] = self
        # The fist defined chart is considered as the default chart on the 
        # manifold:
        if manifold.def_chart is None: 
            manifold.def_chart = self
        if manifold.name != 'field R':      
            #!# to avoid circular import of RealLine
            self.frame = CoordBasis(self)
        self.zero_function = ZeroFunctionChart(self)
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "chart " 
        if self.name is not None:
            description += "'%s' " % self.name
        description += str(self.xx)
        return description
    
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        from sage.misc.latex import latex
        description = "("
        n = len(self.xx)
        for i in range(n-1):
            description += latex(self.xx[i]).strip() + ', '
        description += latex(self.xx[n-1]).strip() + ')'
        return description

    def __call__(self, i=None):
        r"""
        Access to the coordinates.
        
        INPUT:
        
        - ``i`` -- index of the coordinate; if None all the coordinates 
            are returned
            
        OUTPUT: 
        
        - the coordinate of index ``i`` or all the coordinates (as a tuple) if 
          ``i`` is None
        """
        if i is None: 
            return self.xx
        else: 
            return self.xx[i-self.manifold.sindex]

#*****************************************************************************

class FunctionChart(SageObject):
    r"""
    Real-valued function of coordinates belonging to a chart on a manifold. 
    
    Given a chart `\varphi` on a manifold `M` of dimension `n`, an instance of 
    the class :class:`FunctionChart` is a function

    .. MATH::

        \begin{array}{llcl}
        f:& U \subset\RR^n & \longrightarrow & \RR \\
          & (x^1,\ldots,x^n) & \longmapsto & f(x^1,\ldots,x^n)
        \end{array}
    
    where `U` is the domain of `\RR^n` covered by the chart `\varphi`. 
    
    INPUT:
    
    - ``chart`` -- the chart defining the coordinates
    - ``expression`` -- the coordinate expression of the function

    EXAMPLES:
    
    Function defined on a 2-dimensional chart::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy = Chart(m, 'x y', 'c_xy')
        sage: f = FunctionChart(c_xy, x^2+3*y+1)
        sage: f.chart
        chart 'c_xy' (x, y)
        sage: f.show()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    The symbolic expression is also returned when asking the direct display of
    the function::
    
        sage: f
        x^2 + 3*y + 1
        sage: latex(f)
        x^{2} + 3 \, y + 1

    or via the method :meth:`expr`::
    
        sage: f.expr()
        x^2 + 3*y + 1

    The value of the function at specified coordinates is obtained by means
    of the standard parentheses notation::
        
        sage: f(2,-1)
        2
        sage: var('a b')
        (a, b)
        sage: f(a,b)
        a^2 + 3*b + 1

    An unspecified function on a chart::
    
        sage: g = FunctionChart(c_xy, function('G', x, y))
        sage: g
        G(x, y)
        sage: g.show()
        (x, y) |--> G(x, y)
        sage: g.expr()
        G(x, y)
        sage: g(2,3)
        G(2, 3)

    Chart functions can be compared to other values::
    
        sage: f = FunctionChart(c_xy, x^2+3*y+1)
        sage: f == 2
        False
        sage: f == x^2 + 3*y + 1
        True
        sage: g = FunctionChart(c_xy, x*y) 
        sage: f == g
        False
        sage: h = FunctionChart(c_xy, x^2+3*y+1) 
        sage: f == h
        True

    Usage in a physical context (simple Lorentz transformation - boost in 
    x direction, with relative velocity v between o1 and o2 frames)::

        sage: m = Manifold(2, "M")
        sage: o1 = Chart(m, 't, x', 'o1')
        sage: o2 = Chart(m, 'T, X', 'o2')
        sage: f = ScalarField(m, x^2 - t^2); f.express
        {'o1': -t^2 + x^2}
        sage: v = var('v'); gam = 1/sqrt(1-v^2)
        sage: CoordChange(o2, o1, gam*(T - v*X), gam*(X - v*T))
        coordinate change from chart 'o2' (T, X) to chart 'o1' (t, x)
        sage: f.function_chart('o2')
        -T^2 + X^2
        sage: g = f.function_chart('o2')/(X - T)^2 # registering a new FunctionChart
        sage: g.chart
        chart 'o2' (T, X)
        sage: g
        -(T + X)/(T - X)

    """
    def __init__(self, chart, expression): 
        from sage.symbolic.ring import SR
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self.chart = chart
        self.express = SR(expression)
        self.nc = len(self.chart.xx)    # number of coordinates
        # Derived quantities:
        self._der = None  # partial derivatives

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return str(self.express)

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        from sage.misc.latex import latex
        return latex(self.express)

    def expr(self):
        r"""
        Return the expression of the image of the function.
        
        This method actually provides the access to the attribute 
        :attr:`express` that stores the coordinate expression of the function.
        
        OUTPUT:
        
        - symbolic expression, involving the chart coordinates.
        
        EXAMPLES:
        
        Function on some chart of a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: f
            x^2 + 3*y + 1
            sage: f.expr()
            x^2 + 3*y + 1
            sage: print type(f.expr())
            <type 'sage.symbolic.expression.Expression'>
            sage: f.expr() is f.express
            True

        The method :meth:`expr` is useful for accessing to all the 
        symbolic expression functionalities in Sage; for instance::
        
            sage: a = var('a')
            sage: f = FunctionChart(c_xy, a*x*y)
            sage: f.expr()
            a*x*y
            sage: f.expr().subs(a=2)
            2*x*y
        
        Note that for substituting the value of a coordinate, the function call
        can be used as well::
        
            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True

        """
        return self.express
        
    def show(self):
        r"""
        Displays the function in arrow notation.
        
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLE:
        
        Function on a 2-dimensional chart::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: f.show()
            (x, y) |--> x^2 + 3*y + 1
            sage: latex(f.show())
            (x, y) \mapsto x^{2} + 3 \, y + 1

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        result.txt = repr((self.chart)()) + " |--> " + repr(self.express)
        result.latex = latex(self.chart) + r"\mapsto" + latex(self.express)
        return result

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._der = None

    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The derived quantities are not copied, because they can be reconstructed
        if necessary.

        EXAMPLES:
        
        Copy on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: g = f.copy()
            sage: print type(g)
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: g
            x^2 + 3*y + 1
            sage: g == f    # g is mathematically equal to f:
            True
            sage: g is f    # but differs in computer memory:
            False
        
        """
        return FunctionChart(self.chart, self.express)
        
    def __call__(self, *coords):
        r"""
        Computes the value of the function at specified coordinates.
        
        INPUT:
        
        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the 
          function `f` is to be evaluated 
        
        OUTPUT:
        
        - the value `f(x^1,...,x^n)`  
         
        """
        #!# This should be the Python 2.7 form: 
        # substitutions = {self.chart.xx[j]: coords[j] for j in range(self.nc)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(self.chart.xx[j], coords[j]) for j in 
                                                               range(self.nc)])
        resu = self.express.subs(substitutions)
        return simplify_chain(resu)
                     
    def diff(self, coord):
        r""" 
        Partial derivative with respect to a coordinate.
    
        INPUT:
        
        - ``coord`` -- either the coordinate `x^i` with respect 
          to which the derivative of the function `f` is to be taken, or the 
          index `i` labelling this coordinate
          
        OUTPUT:
        
        - the partial derivative `\frac{\partial f}{\partial x^i}`, as an
          instance of :class:`FunctionChart`
          
        EXAMPLES:
        
        Partial derivatives of a function defined on a 2-dimensional chart::

            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1) ; f
            x^2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        The partial derivatives are instances of the class 
        :class:`FunctionChart`::
        
            sage: print type(f.diff(x))
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
        
        An index can be used instead of the coordinate symbol::
        
            sage: f.diff(0)
            2*x
            sage: f.diff(0) is f.diff(x)
            True
            
        The index range depends on the convention used on the manifold::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: f.diff(1)
            2*x
            sage: f.diff(1) is f.diff(x)
            True
            
        """
        from sage.calculus.functional import diff
        if self._der is None:
            # the partial derivatives have to be updated
            self._der = [FunctionChart(self.chart,
                         simplify_chain(diff(self.express, self.chart.xx[j])))
                                                    for j in range(self.nc) ]
        if isinstance(coord, (int, Integer)):
            return self._der[coord - self.chart.manifold.sindex]
        else:
            return self._der[self.chart.xx.index(coord)]

    def is_zero(self):
        r""" 
        Return True if the function is zero and False otherwise.
        
        EXAMPLES:
        
        Functions on a 2-dimensional chart::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'c_xy')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: f.is_zero()
            False
            sage: g = FunctionChart(c_xy, 0)
            sage: g.is_zero()
            True

        """
        return self.express.is_zero()
        
    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                return False
            else:
                return bool(other.express == self.express)
        else:
            return bool(self.express == other)

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
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
        return FunctionChart(self.chart, self.express)

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the opposite of the function ``self``
    
        """
        return FunctionChart(self.chart, simplify_chain(-self.express))

    def __add__(self, other):
        r"""
        Addition operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the addition of ``self`` and ``other``
        
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                raise TypeError("Two functions not defined on the same " + 
                                "chart cannot be added.")
            if isinstance(other, ZeroFunctionChart):
                return self.copy()
            res = simplify_chain(self.express + other.express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self.express + other)
        else:
            return other.__radd__(self)
        if res == 0:
            return self.chart.zero_function
        else:
            return FunctionChart(self.chart, res)

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        """
        return self.__add__(other)
        
    def __iadd__(self, other):
        r"""
        In-place addition operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
                
        """
        return self.__add__(other)

    def __sub__(self, other):
        r"""
        Subtraction operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the subtraction of ``other`` from 
          ``self``
        
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                raise TypeError("Two functions not defined on the same " + 
                                "chart cannot be subtracted.")
            if isinstance(other, ZeroFunctionChart):
                return self.copy()
            res = simplify_chain(self.express - other.express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self.express - other)
        else:
            return other.__rsub__(self)
        if res == 0:
            return self.chart.zero_function
        else:
            return FunctionChart(self.chart, res)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 
        
        """
        return (-self).__add__(other) 

    def __isub__(self, other):
        r"""
        In-place subtraction operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
                
        """
        return self.__sub__(other)


    def __mul__(self, other):
        r"""
        Multiplication operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the multiplication of ``self`` and 
          ``other`
                
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                raise TypeError("Two functions not defined on the same " + 
                                "chart cannot be multiplied.")
            if isinstance(other, ZeroFunctionChart):
                return self.chart.zero_function
            res = simplify_chain(self.express * other.express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self.express * other)
        else:
            return other.__rmul__(self)
        if res == 0:
            return self.chart.zero_function
        else:
            return FunctionChart(self.chart, res)

    def __rmul__(self, other):
        r"""
        Multiplication on the left by ``other``. 
        
        """
        return self.__mul__(other)

    def __imul__(self, other):
        r"""
        In-place multiplication operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
                
        """
        return self.__mul__(other)


    def __div__(self, other):
        r"""
        Division operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the division of ``self`` by 
          ``other`
                
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                raise TypeError("Two functions not defined on the same " + 
                                "chart cannot be divided.")
            if isinstance(other, ZeroFunctionChart):
                raise ZeroDivisionError("Division of a FunctionChart by zero.")
            res = simplify_chain(self.express / other.express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self.express / other)
        else:
            if other == 0:
                raise ZeroDivisionError("Division of a FunctionChart by zero.")
            return other.__rdiv__(self)
        if res == 0:
            return self.chart.zero_function
        else:
            return FunctionChart(self.chart, res)

    def __rdiv__(self, other):
        r"""
        Division of ``other`` by ``self``. 
        
        """
        #!# to be improved
        res = simplify_chain(other / self.express)
        return FunctionChart(self.chart, res)


    def __idiv__(self, other):
        r"""
        In-place division operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
                
        """
        return self.__div__(other)


    def factor(self):
        r"""
        Factorize the coordinate expression. 
        
        OUTPUT:
        
        - ``self``, with ``self.express`` factorized
        
        EXAMPLES:
        
        Factorization on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: X = Chart(m, 'x y', 'xy')
            sage: f = FunctionChart(X, x^2 + 2*x*y + y^2)
            sage: f
            x^2 + 2*x*y + y^2
            sage: f.factor()
            (x + y)^2
            sage: f # the method factor() has changed f:
            (x + y)^2

        """
        self.express = self.express.factor()
        self._del_derived()
        self._del_derived()
        return self
        
 
#*****************************************************************************

class ZeroFunctionChart(FunctionChart):
    r"""
    Null function of coordinates belonging to a chart on a manifold. 

    INPUT:
    
    - ``chart`` -- the chart on which the null function is defined

    EXAMPLES:
    
    Null function defined on a 2-dimensional chart::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy = Chart(m, 'x y', 'c_xy')
        sage: f = ZeroFunctionChart(c_xy) ; f
        0
        sage: f.show()
        (x, y) |--> 0
        sage: f.expr()
        0
        sage: f.is_zero()            
        True
        sage: f(1,2)
        0

    Each chart has its zero function::

        sage: c_xy.zero_function
        0
        sage: print type(c_xy.zero_function)
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        sage: f == c_xy.zero_function
        True

    Arithmetics between instances of :class:`ZeroFunctionChart`::
    
        sage: g = ZeroFunctionChart(c_xy)    
        sage: s = f+g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f-g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f*g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/g ; print type(s) ; s
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a ZeroFunctionChart by zero.

    Arithmetics with a nonzero instance of :class:`FunctionChart`::

        sage: g = FunctionChart(c_xy, x+y)
        sage: s = f+g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = g+f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = f-g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        -x - y
        sage: s = g-f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = f*g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = g*f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/g ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = g/f ; print type(s) ; s
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a FunctionChart by zero.

    Arithmetics with a symbolic expression::

        sage: s = f+x ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = x+f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = f-x ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        -x
        sage: s = x-f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = f*x ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = x*f ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/x ; print type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0

    """
    def __init__(self, chart): 
        FunctionChart.__init__(self, chart, 0)

    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The derived quantities are not copied, because they can be reconstructed
        if necessary.

        """
        return ZeroFunctionChart(self.chart)
        
    def __call__(self, *coords):
        r"""
        Computes the value of the function at specified coordinates.
        
        INPUT:
        
        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the 
          function `f` is to be evaluated 
        
        OUTPUT:
        
        - the value `f(x^1,...,x^n)`  
         
        """
        return 0    #!# SR(0) instead ? 
                     
    def diff(self, coord):
        r""" 
        Partial derivative with respect to a coordinate.
    
        INPUT:
        
        - ``coord`` -- the coordinate `x^i` with respect 
          to which the derivative of the function `f` is to be taken, or the 
          index `i` labelling this coordinate
          
        OUTPUT:
        
        - the partial derivative `\frac{\partial f}{\partial x^i}`, as an
          instance of :class:`ZeroFunctionChart`
          
        EXAMPLES:
        
        """
        if self._der is None:
            self._der = [self.chart.zero_function for j in range(self.nc)]
        return self._der[0]

    def is_zero(self):
        r""" 
        Return True if the function is zero and False otherwise.
        
        """
        return True
        
    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other``,  or False otherwise
        
        """
        if isinstance(other, FunctionChart):
            if other.chart != self.chart:
                return False
            else:
                return other.is_zero()
        else:
            return bool(isinstance(other, (int, Integer)) and other==0)

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
        
        - `self`` (since ``self`` is zero)
    
        """
        return self

    def __add__(self, other):
        r"""
        Addition operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the addition of ``self`` and ``other``
        
        """
        if isinstance(other, FunctionChart):
            if other.chart.name != self.chart.name:
                raise TypeError("Two functions not defined on the same chart " + 
                                "cannot be added.")
            return other.copy()
        elif isinstance(other, (int, RingElement)):  #!# check
            if other == 0:
                return self
            else:
                return FunctionChart(self.chart, other)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        r"""
        Subtraction operator.  
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the subtraction of ``other`` from 
          ``self``
        
        """
        if isinstance(other, FunctionChart):
            if other.chart.name != self.chart.name:
                raise TypeError("Two functions not defined on the same chart " + 
                                "cannot be subtracted.")
            return -other    
        elif isinstance(other, (int, RingElement)):  #!# check
            if other == 0:
                return self
            else:
                return FunctionChart(self.chart, -other)
        else:
            return other.__rsub__(self)
 
    def __mul__(self, other):
        r"""
        Multiplication operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the multiplication of ``self`` and 
          ``other`
                
        """
        if isinstance(other, (int, RingElement, FunctionChart)):  #!# check
            return self
        else:
            return other.__rmul__(self)
        
    def __div__(self, other):
        r"""
        Division operator. 
        
        INPUT:
        
        - ``other`` -- another instance of :class:`FunctionChart` or a value
        
        OUTPUT:
        
        - chart function resulting from the division of ``self`` by 
          ``other`
                
        """
        if isinstance(other, (int, RingElement, FunctionChart)):  #!# check
            if other == 0:
                raise ZeroDivisionError("Division of a ZeroFunctionChart by " + 
                                        "zero.")
            else:
                return self
        else:
            return other.__rdiv__(self)


#*****************************************************************************

class MultiFunctionChart(SageObject):
    r"""
    Class for handling a set of `m` real-valued functions  of
    the coordinates of a given chart. 

    Given an integer `m \geq 1` and a chart `\varphi` on a manifold `M` of 
    dimension `n`, an instance of the class :class:`MultiFunctionChart` is a 
    function

    .. MATH::

        \begin{array}{llcl}
        f:& U \subset\RR^n & \longrightarrow & \RR^m \\
          & (x^1,\ldots,x^n) & \longmapsto & (f_1(x^1,\ldots,x^n),\ldots, 
            f_m(x^1,\ldots,x^n))
        \end{array}
    
    where `U` is the domain of `\RR^n` covered by the chart `\varphi`. 

    INPUT:
    
    - ``chart`` -- the chart defining the coordinates
    - ``*expressions`` -- the list of the coordinate expressions of the `m` 
      functions (`m\geq 1`)
    
    EXAMPLES: 
    
    A set of 3 functions of 2 coordinates::
    
        sage: m = Manifold(2)
        sage: c_xy  = Chart(m, 'x y', 'xy-coord') 
        sage: f = MultiFunctionChart(c_xy, x-y, x*y, cos(x)*exp(y))
        sage: f
        functions (x - y, x*y, e^y*cos(x)) on the chart 'xy-coord' (x, y)
        sage: f.functions
        (x - y, x*y, e^y*cos(x))
        sage: latex(f)
        \left(x - y, x y, e^{y} \cos\left(x\right)\right)
        
    A single function (just as an example; one should rather employ the class
    :class:`FunctionChart` in this case)::
    
        sage: g = MultiFunctionChart(c_xy, x*y^2)
        sage: g.functions
        (x*y^2,)
    
    Evaluating the functions at specified coordinates::
 
        sage: f(1,2)
        (-1, 2, e^2*cos(1))
        sage: g(1,2)
        (4,)
        
    The Jacobian matrix::
    
        sage: f.jacobian()
        [          1          -1]
        [          y           x]
        [-e^y*sin(x)  e^y*cos(x)]
        sage: g.jacobian()
        [  y^2 2*x*y]
    
    If the number of functions equals the number of coordinates, the Jacobian
    determinant can be evaluated::
    
        sage: h = MultiFunctionChart(c_xy, x-y, x*y)
        sage: h.jacobian_det()
        x + y
        
    """
    def __init__(self, chart, *expressions): 
        from sage.symbolic.ring import SR
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self.chart = chart
        self.nc = len(self.chart.xx)    # number of coordinates
        self.nf = len(expressions)      # number of functions
        self.functions = tuple([SR(expressions[i]) for i in range(self.nf)])
        self._jacob = None
        self._jacob_det = None
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "functions " + str(self.functions) + " on the " + \
                      str(self.chart) 
        return description
        
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        from sage.misc.latex import latex
        return latex(self.functions)
    
    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The derived quantities (Jacobian matrix) are not copied, because they 
        can be reconstructed if necessary.
        
        """
        return MultiFunctionChart(self.chart, *(self.functions))

    def jacobian(self):
        r"""
        Returns the Jacobian matrix of the system of functions.
        
        ``jacobian()`` is a matrix of size `m\times n` where `m` is the number 
        of functions and `n` the number of coordinates, the generic element 
        being `J_{ij} = \frac{\partial f_i}{\partial x^j}` with `1\leq i \leq m`
        (row index) and `1\leq j \leq n` (column index). 
        """
        from sage.matrix.constructor import matrix
        from sage.calculus.functional import diff
        if self._jacob is None:
            self._jacob = matrix( [[ 
                    simplify_chain(diff(self.functions[i], self.chart.xx[j]))
                    for j in range(self.nc) ] for i in range(self.nf) ] )
        return self._jacob
        
    def jacobian_det(self):
        r"""
        Returns the Jacobian determinant of the system of functions.
        
        The number `m` of functions must equal the number `n` of 
        coordinates.
        
        """
        from utilities import simple_determinant
        if self._jacob_det is None: 
            if (self.nf != self.nc):
                raise ValueError("The Jacobian matrix is not square.")
            #!# the following is a workaround for a bug in Sage (cf. trac ticket #14403)
            self._jacob_det = simplify_chain(simple_determinant(self.jacobian()))
            # the proper writing should be this:
            # self._jacob_det = simplify_chain(self.jacobian().det())
        return self._jacob_det
  
    def __call__(self, *coords):
        r"""
        Computes the values of the functions at specified coordinates.
        
        INPUT:
        
        - ``*coords`` -- list of coordinates where the functions are to be
          evaluated 
        
        OUTPUT:
        
        - the values of the `m` functions.   
         
        """
        #!# This should be the Python 2.7 form: 
        # substitutions = {self.chart.xx[j]: coords[j] for j in range(self.nc)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(self.chart.xx[j], coords[j]) for j in 
                                                               range(self.nc)])
        return tuple(simplify_chain(self.functions[i].subs(substitutions))
                     for i in range(self.nf) )
       
           

#*****************************************************************************


class CoordChange(SageObject):
    r"""
    Class for changes of coordinates (transition maps between charts).

    The two charts may belong to different manifolds. 
    
    INPUT:
    
    - ``chart1`` -- initial chart
    - ``chart2`` -- final chart 
    - ``transformations`` -- the coordinate transformations expressed as a list 
      of the expressions of the "new" coordinates in terms of the "old" ones  
    
    EXAMPLES: 

    Change from spherical to Cartesian coordinates on `\RR^3`::
    
        sage: m = Manifold(3, 'R3', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
        sage: c_cart = Chart(m, 'x y z', 'cart')        
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: ch
        coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z)
        sage: latex(ch)
        (r, \theta, \phi) \mapsto (x, y, z)

    Each created coordinate change is automatically added to the manifold's 
    dictionary :attr:`coord_changes`; this dictionary is accessed via the method :meth:`Manifold.coord_change`::    

        sage: m.coord_change('spher', 'cart')
        coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z)
    
    It also generates a new entry in the manifold's dictionary 
    :attr:`frame_changes`, containing the relevant change-of-basis matrix; 
    this dictionary is accessed via the method :meth:`Manifold.frame_change`::

        sage: m.frame_change('cart_b', 'spher_b')
        field of tangent-space automorphisms on the 3-dimensional manifold 'R3'
        sage: m.frame_change('cart_b', 'spher_b')[:]
        [   sin(th)*cos(ph)  r*cos(ph)*cos(th) -r*sin(ph)*sin(th)]
        [   sin(ph)*sin(th)  r*sin(ph)*cos(th)  r*sin(th)*cos(ph)]
        [           cos(th)         -r*sin(th)                  0]
    
    The coordinate change can be called directly on a set of "old" coordinates 
    to get the "new" ones::
    
        sage: ch(1,pi/2,0)
        (1, 0, 0)
        
    The Jacobian matrix of the coordinate change::
    
        sage: ch.jacobian
        [   sin(th)*cos(ph)  r*cos(ph)*cos(th) -r*sin(ph)*sin(th)]
        [   sin(ph)*sin(th)  r*sin(ph)*cos(th)  r*sin(th)*cos(ph)]
        [           cos(th)         -r*sin(th)                  0]
        sage: ch.jacobian_det  # Jacobian determinant
        r^2*sin(th)
        
    """
    def __init__(self, chart1, chart2, *transformations): 
        from sage.matrix.constructor import matrix
        from sage.calculus.functional import diff
        from rank2field import AutomorphismField
        n1 = len(chart1.xx)
        n2 = len(chart2.xx)
        if len(transformations) != n2:
            raise ValueError(str(n2) + 
                             " coordinate transformations must be provided.")
        self.chart1 = chart1
        self.chart2 = chart2
        self.transf = MultiFunctionChart(chart1, *transformations)
        self._inverse = None
        # Jacobian matrix: 
        self.jacobian  = self.transf.jacobian()        
        # Jacobian determinant: 
        if n1 == n2: 
            self.jacobian_det = self.transf.jacobian_det()
        # If the two charts are on the same manifold, the coordinate change is 
        # added to the manifold dictionary and the Jacobian matrix is added to 
        # the dictionary of changes of frame:
        if chart1.manifold == chart2.manifold:
            manif = chart1.manifold
            manif.coord_changes[(chart1.name, chart2.name)] = self
            frame1 = chart1.frame.name
            frame2 = chart2.frame.name
            ch_basis = AutomorphismField(manif)
            ch_basis.set_comp(frame1)[:, chart1.name] = self.jacobian
            ch_basis.set_comp(frame2, delete_others=False)[:, chart1.name] = \
                                                                  self.jacobian
            manif.frame_changes[(chart2.frame.name, 
                                 chart1.frame.name)] = ch_basis
            if (chart1.frame.name, chart2.frame.name) not in manif.frame_changes:
                manif.frame_changes[(chart1.frame.name, 
                                     chart2.frame.name)] = ch_basis.inverse()

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "coordinate change from " + str(self.chart1) + " to " + \
                      str(self.chart2)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        from sage.misc.latex import latex
        return latex(self.chart1) + "\mapsto" + latex(self.chart2)

    def __call__(self, *old_coords):
        r"""
        Computes the new coordinates from old ones.
        """
        return self.transf(*old_coords)

    def inverse(self):
        r""" 
        Computes the inverse coordinate transformation, when the latter is
        invertible. 
        
        OUTPUT:
        
        - an instance of :class:`CoordChange` representing the inverse of
          ``self``. 
          
        EXAMPLES:
        
        Inverse of a coordinate transformation corresponding to a pi/3-rotation
        in the plane::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy-coord')
            sage: c_uv = Chart(m, 'u v', 'uv-coord')
            sage: ch_to_uv = CoordChange(c_xy, c_uv, (x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2)
            sage: m.coord_changes                                                         
            {('xy-coord',
              'uv-coord'): coordinate change from chart 'xy-coord' (x, y) to chart 'uv-coord' (u, v)}
            sage: ch_to_xy = ch_to_uv.inverse() ; ch_to_xy
            coordinate change from chart 'uv-coord' (u, v) to chart 'xy-coord' (x, y)
            sage: ch_to_xy.transf                                                         
            functions (1/2*sqrt(3)*v + 1/2*u, -1/2*sqrt(3)*u + 1/2*v) on the chart 'uv-coord' (u, v)
            sage: m.coord_changes
            {('xy-coord', 'uv-coord'): coordinate change from chart 'xy-coord' (x, y) to chart 'uv-coord' (u, v), ('uv-coord', 'xy-coord'): coordinate change from chart 'uv-coord' (u, v) to chart 'xy-coord' (x, y)}
   
        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        if self._inverse is not None:
            return self._inverse
            
        n1 = len(self.chart1.xx)
        n2 = len(self.chart2.xx)
        if n1 != n2:
            raise TypeError("The change of coordinates is not invertible.")
        
        # New symbolic variables (different from chart2.xx to allow for a 
        #  correct solution even when chart2 = chart1):
        coord_domain = ['real' for i in range(n2)]
        for i in range(n2):
            if self.chart2.xx[i].is_positive():
                coord_domain[i] = 'positive'
        x2 = [ SR.var('xxxx' + str(i), domain=coord_domain[i]) for i in 
               range(n2) ]
        equations = [x2[i] == self.transf.functions[i] for i in range(n2) ]
        solutions = solve(equations, self.chart1.xx, solution_dict=True)
        if len(solutions) == 0: 
            raise ValueError("No solution found")
        if len(solutions) > 1: 
            raise ValueError("Non-unique solution found")
            
        #!# This should be the Python 2.7 form: 
        # substitutions = {x2[i]: self.chart2.xx[i] for i in range(n2)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(x2[i], self.chart2.xx[i]) for i in range(n2)])
       
        inv_transf = [solutions[0][self.chart1.xx[i]].subs(substitutions) 
                           for i in range(n1)]
        for i in range(n1):
            try:
                inv_transf[i] = simplify_chain(inv_transf[i])
            except AttributeError:
                pass        
        self._inverse = CoordChange(self.chart2, self.chart1, *inv_transf)
        
        # Update of chart expressions of the frame changes:
        if self.chart1.manifold == self.chart2.manifold:
            manif = self.chart1.manifold
            frame1 = self.chart1.frame.name
            frame2 = self.chart2.frame.name
            fr_change12 = manif.frame_changes[(frame1,frame2)]
            fr_change21 = manif.frame_changes[(frame2,frame1)]
            for comp in fr_change12.components[frame1]._comp.values():
                comp.function_chart(self.chart1.name, 
                                    from_chart=self.chart2.name)
            for comp in fr_change12.components[frame2]._comp.values():
                comp.function_chart(self.chart1.name, 
                                    from_chart=self.chart2.name)
            for comp in fr_change21.components[frame1]._comp.values():
                comp.function_chart(self.chart2.name, 
                                    from_chart=self.chart1.name)
            for comp in fr_change21.components[frame2]._comp.values():
                comp.function_chart(self.chart2.name, 
                                    from_chart=self.chart1.name)

        return self._inverse


    def set_inverse(self, *transformations, **kwds):
        r"""
        Sets the inverse of the coordinate transformation. 
        
        This is usefull when the automatic computation via :meth:`inverse()`
        fails. 
        
        INPUT:
        
        - ``transformations`` -- the inverse transformations expressed as a list 
          of the expressions of the "old" coordinates in terms of the "new" ones
        - ``kwds`` -- keyword arguments: only ``check=True`` (default) or
          ``check=False`` is meaningfull; it determines whether the provided transformations are checked to be indeed the inverse coordinate
          transformations. 
          
        EXAMPLES:
         
        From Cartesian to spherical coordinates in the plane::
          
            sage: m = Manifold(2, 'M')
            sage: c_cart = Chart(m, 'x y', 'cart')
            sage: c_spher = Chart(m, 'r:positive ph:\phi', 'spher')
            sage: spher_to_cart = CoordChange(c_spher, c_cart, r*cos(ph), r*sin(ph))
            sage: spher_to_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))              
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: spher_to_cart.inverse()                             
            coordinate change from chart 'cart' (x, y) to chart 'spher' (r, ph)
            sage: m.coord_changes
            {('cart',
              'spher'): coordinate change from chart 'cart' (x, y) to chart 'spher' (r, ph),
             ('spher',
              'cart'): coordinate change from chart 'spher' (r, ph) to chart 'cart' (x, y)}
              
        Introducing a wrong inverse transformation is revealed by the check::
                
            sage: spher_to_cart.set_inverse(sqrt(x^3+y^2), atan2(y,x)) # note the x^3 typo
            Check of the inverse coordinate transformation:
               r == sqrt(r*cos(ph)^3 + sin(ph)^2)*r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == sqrt(x^3 + y^2)*x/sqrt(x^2 + y^2)
               y == sqrt(x^3 + y^2)*y/sqrt(x^2 + y^2)
            sage: # the check clearly fails

        """ 
        if 'check' in kwds:
            check = kwds['check']
        else:
            check = True
        self._inverse = CoordChange(self.chart2, self.chart1, *transformations)
        if check:
            print "Check of the inverse coordinate transformation:"
            x1 = self.chart1.xx
            x2 = self.chart2.xx
            n1 = len(x1)
            for i in range(n1):
                print "  ", x1[i], '==' , self._inverse(*(self(*x1)))[i]
            for i in range(n1):
                print "  ", x2[i], '==', self(*(self._inverse(*x2)))[i]
        # Update of chart expressions of the frame changes:
        if self.chart1.manifold == self.chart2.manifold:
            manif = self.chart1.manifold
            frame1 = self.chart1.frame.name
            frame2 = self.chart2.frame.name
            fr_change12 = manif.frame_changes[(frame1,frame2)]
            fr_change21 = manif.frame_changes[(frame2,frame1)]
            for comp in fr_change12.components[frame1]._comp.values():
                comp.function_chart(self.chart1.name, 
                                    from_chart=self.chart2.name)
            for comp in fr_change12.components[frame2]._comp.values():
                comp.function_chart(self.chart1.name, 
                                    from_chart=self.chart2.name)
            for comp in fr_change21.components[frame1]._comp.values():
                comp.function_chart(self.chart2.name, 
                                    from_chart=self.chart1.name)
            for comp in fr_change21.components[frame2]._comp.values():
                comp.function_chart(self.chart2.name, 
                                    from_chart=self.chart1.name)
    
