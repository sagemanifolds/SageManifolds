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
from domain import OpenDomain
from utilities import simplify_chain

class Chart(SageObject):
    r"""
    Class for charts on a manifold.
    
    Given a manifold `M` of dimension `n`, a *chart* is a pair `(U,\varphi)`, 
    where `U` is an open domain of `M` and 
    `\varphi: U \rightarrow V \subset \RR^n` is a homeomorphism from `U` to 
    an open domain `V` of `\RR^n`. 
    
    The components `(x^1,\ldots,x^n)` of `\varphi`, defined by 
    `\varphi(p) = (x^1(p),\ldots,x^n(p))`, are called the *coordinates* of the
    chart `(U,\varphi)`.

    INPUT:
    
    - ``domain`` -- open domain `U` on which the chart is defined (must be 
      an instance of :class:`OpenDomain`)
    - ``coordinates`` -- single string defining the coordinate symbols and 
      ranges: the coordinates are separated by ' ' (space) and each coordinate 
      has at most three fields, separated by ':': 
        
        1. The coordinate symbol (a letter or a few letters)
        2. (optional) The interval `I` defining the coordinate range: if not
           provided, the coordinate is assumed to span all `\RR`; otherwise 
           `I` must be provided in the form (a,b) (or equivalently ]a,b[)
           The bounds a and b can be +/-Infinity, Inf, infinity, inf or oo.
           For *singular* coordinates, non-open intervals such as [a,b] and 
           (a,b] (or equivalently ]a,b]) are allowed. 
           Note that the interval declaration must not contain any space 
           character.
        3. (optional) The LaTeX spelling of the coordinate; if not provided the
           coordinate symbol given in the first field will be used.
      
      The order of the fields 2 and 3 does not matter and each of them can be
      omitted.
      If it contains any LaTeX expression, the string ``coordinates`` must be
      declared with the prefix 'r' (for "raw") to allow for a proper treatment 
      of the backslash character (see examples below). 
    - ``names`` -- (default: None) unused argument (present only to enable
      the use of the shortcut operator <,>). 
    
    EXAMPLES: 
    
    Cartesian coordinates on `\RR^3`::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart = Chart(M, 'x y z') ; c_cart
        chart (R^3, (x, y, z))

    To have the coordinates accessible as global variables, one shall set::
    
        sage: (x,y,z) = c_cart[:]
        
    However, a shortcut is to use the declarator ``<x,y,z>`` during the Chart
    construction::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = Chart(M, 'x y z') ; c_cart
        chart (R^3, (x, y, z))
    
    The coordinates are then immediately accessible::
    
        sage: y
        y
        sage: y is c_cart[2]
        True
    
    Note that x, y, z declared in ``<x,y,z>`` are mere Python variable names 
    and do not have to coincide with the coordinate symbols; for instance, 
    one may write::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x1,y1,z1> = Chart(M, 'x y z') ; c_cart
        chart (R^3, (x, y, z))
    
    Then y is not known as a global variable and the corresponding coordinate
    is accessible only through the global variable y1::
    
        sage: y1
        y
        sage: y1 is c_cart[2]
        True
    
    However, having the name of the Python variable coincide with the 
    coordinate symbol is quite convenient; so it is recommended to declare::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = Chart(M, 'x y z')
    
    Moreover, the method :meth:`OpenDomain.chart` can be invoked instead
    of the direct call to the Chart constructor::
    
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart('x y z')
     
    Spherical coordinates on the subdomain `U` of `\RR^3` that is the 
    complement of the half-plane `\{y=0, x\geq 0\}`::
    
        sage: U = M.open_domain('U')
        sage: c_spher.<r,th,ph> = Chart(U, r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi') ; c_spher 
        chart (U, (r, th, ph))

    Note the prefix 'r' for the string defining the coordinates in the arguments of ``Chart``. 
    
    Coordinates are some Sage symbolic variables::

        sage: print type(th)
        <type 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        \theta
        sage: assumptions(th)
        [th is real, th > 0, th < pi]

    Coordinate are also accessible by their indices::
    
        sage: x1 = c_spher[1]; x2 = c_spher[2]; x3 = c_spher[3]
        sage: print x1, x2, x3
        r th ph
        sage: (x1, x2, x3) == (r, th, ph)
        True

    The full set of coordinates is obtained by means of the operator [:]::
    
        sage: c_cart[:]
        (x, y, z)
        sage: c_spher[:]
        (r, th, ph)
    
    Let us check that the declared coordinate ranges have been taken into 
    account::
        
        sage: bool(th>0 and th<pi)
        True
        sage: assumptions()  # list all current symbolic assumptions
        [x is real, y is real, z is real, r is real, r > 0, th is real, 
         th > 0, th < pi, ph is real, ph > 0, ph < 2*pi]
        
    The coordinate ranges are used for simplifications::
    
        sage: simplify(abs(r)) # r has been declared to lie in the interval (0,+oo)
        r
        sage: simplify(abs(x)) # no positive range has been declared for x
        abs(x)
        sage: from sage.geometry.manifolds.utilities import simplify_abs_trig
        sage: simplify_abs_trig(abs(sin(th)))  # sin(th) is always positive
        sin(th)
        sage: simplify_abs_trig(abs(sin(ph)))  # sin(ph) is not
        abs(sin(ph))

    Each constructed chart is automatically added to the manifold's atlas::
    
        sage: M.atlas
        [chart (R^3, (x, y, z)), chart (U, (r, th, ph))]

    and to the atlas of the domain in which it has been defined::
    
        sage: U.atlas
        [chart (U, (r, th, ph))]

    Each domain has a default chart, which, unless changed via the method
    :meth:`Domain.set_default_chart`, is the first defined chart on that 
    domain (or on a subdomain of it)::
    
        sage: M.default_chart()
        chart (R^3, (x, y, z))
        sage: U.default_chart()
        chart (U, (r, th, ph))
    
    The chart map `\varphi` acting on a point is obtained by means of the
    call operator, i.e. the operator ``()``::
    
        sage: p = M.point((1,0,-2)) ; p
        point on 3-dimensional manifold 'R^3'
        sage: c_cart(p)
        (1, 0, -2)
        sage: c_cart(p) == p.coord(c_cart)
        True
        sage: q = M.point((2,pi/2,pi/3), c_spher) # point defined by its spherical coordinates
        sage: c_spher(q)
        (2, 1/2*pi, 1/3*pi)
        sage: c_spher(q) == q.coord(c_spher)
        True
        sage: a = U.point((1,pi/2,pi)) # the default coordinates on U are the spherical ones
        sage: c_spher(a)
        (1, 1/2*pi, pi)
        sage: c_spher(a) == a.coord(c_spher)
        True

    Cartesian coordinates on U as an example of chart construction with 
    coordinate restrictions: since U is the complement of the half-plane 
    `\{y=0, x\geq 0\}`, we must have `y\not=0` or `x<0` on U. Accordingly, 
    we set::
    
        sage: c_cartU.<x,y,z> = Chart(U, 'x y z') 
        sage: c_cartU.add_restrictions((y!=0, x<0)) # the tuple (y!=0, x<0) means y!=0 or x<0
        sage: # c_cartU.add_restrictions([y!=0, x<0]) would have meant y!=0 AND x<0
        sage: U.atlas
        [chart (U, (r, th, ph)), chart (U, (x, y, z))]
        sage: M.atlas
        [chart (R^3, (x, y, z)), chart (U, (r, th, ph)), chart (U, (x, y, z))]
        sage: c_cartU.valid_coordinates(-1,0,2)
        True
        sage: c_cartU.valid_coordinates(1,0,2)
        False
        sage: c_cart.valid_coordinates(1,0,2)
        True

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
    def __init__(self, domain, coordinates, names=None): 
        from sage.symbolic.ring import SR
        from sage.symbolic.assumptions import assume
        from sage.rings.infinity import Infinity
        from vectorframe import CoordFrame
        if not isinstance(domain, OpenDomain):
            raise TypeError("The first argument must be an open domain.")
        self.manifold = domain.manifold
        self.domain = domain        
        # Treatment of the coordinates:
        if ' ' in coordinates:
            coord_list = coordinates.split()
        else:
            coord_list = [coordinates]
        n = self.manifold.dim 
        if len(coord_list) != n:
            raise ValueError("The list of coordinates must contain " + \
                              str(n) + " elements.")
        xx_list = [] # will contain the coordinates as Sage symbolic variables
        bounds_list = [] # will contain the coordinate bounds
        for coord_field in coord_list: 
            coord_properties = coord_field.split(':')
            coord_symb = coord_properties[0].strip() # the coordinate symbol
            # default values, possibly redefined below:
            coord_latex = None 
            xmin = -Infinity ; xmin_included = False   
            xmax = +Infinity ; xmax_included = False
            # scan of the properties other than the symbol: 
            for prop in coord_properties[1:]:
                prop1 = prop.strip()
                delim_min = prop1[0]
                if delim_min in ['[', ']', '(']:
                    # prop1 is the coordinate's range
                    xmin_str, xmax_str = prop1[1:len(prop1)-1].split(',')
                    if xmin_str not in ['-inf', '-Inf', '-infinity', 
                                        '-Infinity', '-oo']:
                        xmin = SR(xmin_str)
                        xmin_included = ( delim_min == '[' )
                    if xmax_str not in ['inf', '+inf', 'Inf', '+Inf', 
                                        'infinity', '+infinity', 'Infinity',
                                        '+Infinity', 'oo', '+oo']:
                        xmax = SR(xmax_str)
                        xmax_included = ( prop1[-1] == ']' )
                else:
                    # prop1 is the coordinate's LaTeX symbol
                    coord_latex = prop1
            # Construction of the coordinate as a Sage's symbolic variable:
            coord_var = SR.var(coord_symb, domain='real', 
                               latex_name=coord_latex)

            #!# To avoid the warning message
            #       Resolving lazy import maxima during startup
            #       Calling stack:
            #       ...
            # at Sage's startup, one checks that the manifold is not
            # the field R (object RealLine constructed at Sage's startup)
            # before invoking the function assume
            # cf. http://trac.sagemath.org/ticket/14187
            # This is triggered by the lazy import
            #   from sage.calculus.calculus import maxima
            # in line 112 of the file /sage/symbolic/assumptions.py
            if self.domain.name != 'field R':
                assume(coord_var, 'real')

            if xmin != -Infinity:
                if xmin_included:
                    assume(coord_var >= xmin)
                else:
                    assume(coord_var > xmin)
            if xmax != Infinity:
                if xmax_included:
                    assume(coord_var <= xmax)
                else:
                    assume(coord_var < xmax)
            xx_list.append(coord_var)
            bounds_list.append(((xmin, xmin_included), (xmax, xmax_included)))
        self.xx = tuple(xx_list)
        self.bounds = tuple(bounds_list)
        # End of the treatment of the coordinates
        
        # Additional restrictions on the coordinates
        self.restrictions = []  # to be set with method add_restrictions()

        # The chart is added to the domain's atlas, as well as to all the 
        # superdomains' atlases; moreover the fist defined chart is considered 
        # as the default chart
        for sd in self.domain.superdomains:
            sd.atlas.append(self)
            if sd.def_chart is None: 
                sd.def_chart = self
        # Construction of the coordinate frame associated to the chart:
        if self.domain.name != 'field R':      
            #!# to avoid circular import of RealLine
            self.frame = CoordFrame(self)
            self.coframe = self.frame.coframe
        # The null function of the coordinates:
        self.zero_function = ZeroFunctionChart(self)
        # Initialization of the set of charts that are restrictions of the
        # current chart to subdomains of the chart domain:
        self.subcharts = set([self]) 
        # Initialization of the set of charts which the current chart is a 
        # restriction of:
        self.supercharts = set([self]) 
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = 'chart ' + \
                      '(' + self.domain.name + ', ' + str(self.xx) + ')'
        return description
    
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        from sage.misc.latex import latex
        description = '(' + latex(self.domain).strip() + ',('
        n = len(self.xx)
        for i in range(n-1):
            description += latex(self.xx[i]).strip() + ', '
        description += latex(self.xx[n-1]).strip() + '))'
        return description

    def _latex_coordinates(self):
        r"""
        Return a LaTeX representation of the coordinates only. 
        """
        from sage.misc.latex import latex
        description = '(' 
        n = len(self.xx)
        for i in range(n-1):
            description += latex(self.xx[i]).strip() + ', '
        description += latex(self.xx[n-1]).strip() + ')'
        return description

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates.
        """
        return self[:]

    def __hash__(self):
        r"""
        Hash function (since instances of :class:`Chart` are used as
        dictionary keys).
        """
        return id(self)

    def __eq__(self, other):
        r"""
        Comparison operator
        """
        return other is self

    def __getitem__(self, i):
        r"""
        Access to the coordinates.
        
        INPUT:
        
        - ``i`` -- index of the coordinate; if [:] all the coordinates 
            are returned
            
        OUTPUT: 
        
        - the coordinate of index ``i`` or all the coordinates (as a tuple) if 
          ``i`` is [:]
        """
        if isinstance(i, slice): 
            return self.xx
        else: 
            return self.xx[i-self.manifold.sindex]

    def __call__(self, point):
        r"""
        Return the coordinates of a given point. 
        """
        return point.coord(self)

    def coord_bounds(self, i=None):
        r"""
        Return the coordinate lower and upper bounds.
        
        INPUT:
        
        - ``i`` -- index of the coordinate; if None, the bounds of all the 
            coordinates are returned
            
        OUTPUT: 
        
        - the coordinate bounds as the tuple 
          ((xmin, min_included), (xmax, max_included))
          where 
          
          - xmin is the coordinate lower bound 
          - min_included is a Boolean, indicating whether the coordinate can 
            take the value xmin, i.e. xmin is a strict lower bound iff
            min_included is False.
          - xmin is the coordinate upper bound 
          - max_included is a Boolean, indicating whether the coordinate can 
            take the value xmax, i.e. xmax is a strict upper bound iff
            max_included is False.
        
        EXAMPLES:
        
        Some coordinate bounds on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = Chart(m, 'x y:[0,1)')
            sage: c_xy.coord_bounds(0)  # x in (-oo,+oo) (the default)
            ((-Infinity, False), (+Infinity, False))
            sage: c_xy.coord_bounds(1)  # y in [0,1)
            ((0, True), (1, False))
            sage: c_xy.coord_bounds() 
            (((-Infinity, False), (+Infinity, False)), ((0, True), (1, False)))
            sage: c_xy.coord_bounds() == (c_xy.coord_bounds(0), c_xy.coord_bounds(1))
            True
    
        The coordinate bounds can also be recovered via Sage's function 
        :func:`sage.symbolic.assumptions.assumptions`::
        
            sage: assumptions(x)
            [x is real]
            sage: assumptions(y)
            [y is real, y >= 0, y < 1]
            
        """
        if i is None: 
            return self.bounds
        else: 
            return self.bounds[i-self.manifold.sindex]

    def add_restrictions(self, restrictions):
        r"""
        Add some restrictions on the coordinates.
        
        INPUT:
        
        - ``restrictions`` -- list of restrictions on the 
          coordinates, in addition to the ranges declared by the intervals 
          specified in the chart constructor. 
          A restriction can be any symbolic equality or inequality involving the 
          coordinates, such as x>y or x^2+y^2 != 0. The items of the list
          ``restrictions`` are combined with the ``and`` operator; if some
          restrictions are to be combined with the ``or`` operator instead, they 
          have to be passed as a tuple in some single item of the list 
          ``restrictions``. For example, ``restrictions`` = [x>y, (x!=0, y!=0), 
          z^2<x] means (x>y) and ((x!=0) or (y!=0)) and (z^2<x). If the list
          ``restrictions`` contains only one item, this item can be passed as 
          such, i.e. writing x>y instead of the single element list [x>y]. 
    
        EXAMPLES:

        Cartesian coordinates on the open unit disc in $\RR^2$::
        
            sage: M = Manifold(2, 'M') # the open unit disc
            sage: X.<x,y> = Chart(M, 'x y')
            sage: X.add_restrictions(x^2+y^2<1)
            sage: X.valid_coordinates(0,2)
            False
            sage: X.valid_coordinates(0,1/3)
            True

        The restrictions are transmitted to subcharts::
        
            sage: A = M.open_domain('A') # annulus 1/2 < r < 1
            sage: X_A = X.subchart(A, x^2+y^2 > 1/4)
            sage: X_A.restrictions
            [x^2 + y^2 < 1, x^2 + y^2 > (1/4)]
            sage: X_A.valid_coordinates(0,1/3)
            False
            sage: X_A.valid_coordinates(2/3,1/3)
            True

        """
        if not isinstance(restrictions, list): 
            # case of a single condition or conditions to be combined by "or"
            restrictions = [restrictions]
        self.restrictions.extend(restrictions)


    def subchart(self, domain, restrictions):
        r"""
        Construct a subchart.
        
        If ``self`` is the chart `(U,\varphi)`, a subchart is a chart `(V,\psi)`
        such that `V\subset U` and `\psi = \varphi |_V`. 

        The coordinates of the subchart bare the same names as the coordinates
        of the mother chart. 
        
        INPUT:
        
        - ``domain`` -- open subdomain `V` of the chart domain `U` 
        - ``restrictions`` -- list of coordinate restrictions defining 
          the subdomain `V`. A restriction can be any symbolic equality or 
          inequality involving the coordinates, such as x>y or x^2+y^2 != 0. 
          The items of the list ``restrictions`` are combined with the ``and`` 
          operator; if some restrictions are to be combined with the ``or`` 
          operator instead, they have to be passed as a tuple in some single 
          item of the list ``restrictions``. For example, ``restrictions`` 
          being [x>y, (x!=0, y!=0), z^2<x] means (x>y) and ((x!=0) or (y!=0)) 
          and (z^2<x). If the list ``restrictions`` contains only one item, 
          this item can be passed as such, i.e. writing x>y instead of the 
          single element list [x>y]. 

        OUTPUT:
        
        - chart `(V,\psi)`, as an instance of :class:`Chart`. 
        
        EXAMPLES:
        
        Cartesian coordinates on the unit open disc in `\RR^2` as a subchart 
        of the global Cartesian coordinates::
        
            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = Chart(M, 'x y') # Cartesian coordinates on R^2
            sage: D = M.open_domain('D') # the unit open disc
            sage: c_cart_D = c_cart.subchart(D, x^2+y^2<1) 
            sage: p = M.point((1/2, 0))
            sage: p in D
            True
            sage: q = M.point((1, 2))
            sage: q in D
            False

        Cartesian coordinates on the annulus `1<\sqrt{x^2+y^2}<2`::
        
            sage: A = M.open_domain('A')
            sage: c_cart_A = c_cart.subchart(A, [x^2+y^2>1, x^2+y^2<4])
            sage: p in A, q in A
            (False, False)
            sage: a = M.point((3/2,0))
            sage: a in A
            True

        """
        if not domain.is_subdomain(self.domain):
            raise TypeError("The argument 'domain' must be a subdomain of " + 
                            " the chart domain.")
        coordinates = ""
        for coord in self.xx:
            coordinates += repr(coord) + ' '
        res = Chart(domain, coordinates)
        res.bounds = self.bounds
        res.restrictions.extend(self.restrictions)
        res.add_restrictions(restrictions)
        # Update of supercharts and subcharts:
        res.supercharts.update(self.supercharts)
        for schart in self.supercharts:
            schart.subcharts.add(res)
        # Update of superframes and subframes:
        res.frame.superframes.update(self.frame.superframes)
        for sframe in self.frame.superframes:
            sframe.subframes.add(res.frame)
        return res
        
    def valid_coordinates(self, *coordinates):
        r""" 
        Check whether a tuple of coordinates can be the coordinates of a 
        point in the chart domain.

        INPUT:
        
        - ``*coordinates`` -- coordinate values

        OUTPUT:
        
        - True if the coordinate values are admissible in the chart domain. 

        """
        n = len(coordinates)
        if n != self.manifold.dim:
            return False
        # Check of the coordinate ranges:
        for x, bounds in zip(coordinates, self.bounds):
            xmin = bounds[0][0]
            min_included = bounds[0][1]
            xmax = bounds[1][0]
            max_included = bounds[1][1]
            if min_included:
                if x < xmin:
                    return False
            else:
                if x <= xmin:
                    return False
            if max_included:
                if x > xmax:
                    return False
            else:
                if x >= xmax:
                    return False
        # Check of additional restrictions:
        if self.restrictions != []:
            substitutions = dict([(self.xx[j], coordinates[j]) for j in 
                                                                    range(n)])
            for restrict in self.restrictions:
                if isinstance(restrict, tuple): # case of or conditions
                    combine = False
                    for expr in restrict:
                        combine = combine or bool(expr.subs(substitutions))
                    if not combine:
                        return False
                else:
                    if not restrict.subs(substitutions):
                        return False
        # All tests have been passed:
        return True

    def transition_map(self, other, transformations, intersection_name=None, 
                       restrictions1=None, restrictions2=None):
        r""" 
        Construct the transition map between the current chart, 
        `(U,\varphi)` say, and another one, `(V,\psi)` say. 
        
        If `n` is the manifold's dimension, the *transition map* is the
        map 
                
        .. MATH::
        
            \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset \RR^n 
            \rightarrow \psi(U\cap V) \subset \RR^n
        
        In other words, the 
        transition map expresses the coordinates `(y^1,\ldots,y^n)` of 
        `(V,\psi)` in terms of the coordinates `(x^1,\ldots,x^n)` of 
        `(U,\varphi)` on the domain where the two charts intersect, i.e. on 
        `U\cap V`.

        INPUT:
        
        - ``other`` -- the chart `(V,\psi)`
        - ``transformations`` -- tuple (Y_1,...,Y_2), where Y_i is a symbolic
          expression expressing the coordinate `y^i` in terms of the 
          coordinates `(x^1,\ldots,x^n)`
        - ``intersection_name`` -- (default: None) name to be given to the 
          domain `U\cap V` if the latter differs from `U` or `V`
        - ``restrictions1`` -- (default: None) list of conditions on the 
          coordinates of the current chart that define `U\cap V` if the 
          latter differs from `U`. ``restrictions1`` must be a list of 
          of symbolic equalities or inequalities involving the 
          coordinates, such as x>y or x^2+y^2 != 0. The items of the list
          ``restrictions1`` are combined with the ``and`` operator; if some
          restrictions are to be combined with the ``or`` operator instead, 
          they have to be passed as a tuple in some single item of the list 
          ``restrictions1``. For example, ``restrictions1`` = [x>y, 
          (x!=0, y!=0), z^2<x] means (x>y) and ((x!=0) or (y!=0)) and (z^2<x).
          If the list ``restrictions1`` contains only one item, this item can 
          be passed as such, i.e. writing x>y instead of the single element 
          list [x>y]. 
        - ``restrictions2`` -- (default: None) list of conditions on the 
          coordinates of the other chart that define `U\cap V` if the latter 
          differs from `V` (see ``restrictions1`` for the syntax)

        OUTPUT:
        
        - The transition map `\psi\circ\varphi^{-1}` defined on `U\cap V`, as an
          instance of :class:`CoordChange`. 
          
        EXAMPLES:
        
        Transition map between two stereographic charts on the circle `S^1`::
        
            sage: M = Manifold(1, 'S^1')
            sage: U = M.open_domain('U') # Complement of the North pole
            sage: cU.<x> = Chart(U, 'x') # Stereographic chart from the North pole
            sage: V = M.open_domain('V') # Complement of the South pole
            sage: cV.<y> = Chart(V, 'y') # Stereographic chart from the South pole
            sage: trans = cU.transition_map(cV, 1/x, 'W', x!=0, y!=0)
            sage: trans
            coordinate change from chart (W, (x,)) to chart (W, (y,))
            sage: M.domains # the domain W, intersection of U and V, has been created by transition_map()
            {'V': open domain 'V' on the 1-dimensional manifold 'S^1', 
             'U': open domain 'U' on the 1-dimensional manifold 'S^1', 
             'W': open domain 'W' on the 1-dimensional manifold 'S^1', 
             'S^1': 1-dimensional manifold 'S^1'}
            sage: W = M.domains['W']
            sage: W is U.intersection(V)
            True
            sage: M.atlas
            [chart (U, (x,)), chart (V, (y,)), chart (W, (x,)), chart (W, (y,))]

        Transition map between spherical chart and Cartesian chart on `\RR^2`::
        
            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = Chart(M, 'x y')
            sage: U = M.open_domain('U') # the complement of the segment {y=0, x >= 0}
            sage: c_spher.<r,phi> = Chart(U, r'r:(0,+oo) phi:(0,2*pi):\phi')
            sage: trans = c_spher.transition_map(c_cart, (r*cos(phi), r*sin(phi)), \
                                                 restrictions2=(y!=0, x<0))
            sage: trans
            coordinate change from chart (U, (r, phi)) to chart (U, (x, y))
            sage: M.domains # in this case, no new domain has been created since U inter M = U
            {'R^2': 2-dimensional manifold 'R^2',
             'U': open domain 'U' on the 2-dimensional manifold 'R^2'}
            sage: M.atlas # ...but a new chart has been created: (U, (x, y))
            [chart (R^2, (x, y)), chart (U, (r, phi)), chart (U, (x, y))]
        
        """
        dom1 = self.domain
        dom2 = other.domain
        dom = dom1.intersection(dom2, name=intersection_name)
        if dom is dom1:
            chart1 = self
        else:
            chart1 = self.subchart(dom, restrictions1)
        if dom is dom2:
            chart2 = other
        else:
            chart2 = other.subchart(dom, restrictions2)
        if not isinstance(transformations, (tuple, list)):
                transformations = [transformations]
        return CoordChange(chart1, chart2, *transformations)
            

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
        sage: c_xy.<x,y> = m.chart('x y')
        sage: f = FunctionChart(c_xy, x^2+3*y+1)
        sage: f.chart
        chart (M, (x, y))
        sage: f.view()
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
        sage: g.view()
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
            sage: c_xy.<x,y> = m.chart('x y')
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
        
    def view(self):
        r"""
        Displays the function in arrow notation.
        
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLE:
        
        Function on a 2-dimensional chart::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = FunctionChart(c_xy, x^2+3*y+1)
            sage: f.view()
            (x, y) |--> x^2 + 3*y + 1
            sage: latex(f.view())
            (x, y) \mapsto x^{2} + 3 \, y + 1

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        result.txt = repr((self.chart)[:]) + ' |--> ' + repr(self.express)
        result.latex = self.chart._latex_coordinates() + r' \mapsto' + latex(self.express)
        return result

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._der = None

    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The derived quantities are not copied, because they can be 
        reconstructed if necessary.

        EXAMPLES:
        
        Copy on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
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
            sage: c_xy.<x,y> = m.chart('x y')
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
            sage: c_xy.<x,y> = m.chart('x y')
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
            sage: c_xy.<x,y> = m.chart('x y')
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
            sage: X.<x,y> = m.chart('x y')
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
        return self
        
    def scalar_field(self, name=None, latex_name=None):
        r""" 
        Construct the scalar field that has ``self`` as coordinate expression. 
        
        The domain of the scalar field is the domain covered by the chart on 
        which ``self`` is defined.
        
        INPUT: 
        
        - ``name`` -- (default: None) name given to the scalar field
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the scalar 
          field; if none is provided, the LaTeX symbol is set to ``name``
        
        OUTPUT:
        
        - instance of :class:`ScalarField`
        
        EXAMPLES:

        Construction of a scalar field on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')                  
            sage: c_xy.<x,y> = M.chart('x y')
            sage: fc = FunctionChart(c_xy, x+2*y^3)
            sage: f = fc.scalar_field() ; f
            scalar field on the 2-dimensional manifold 'M'
            sage: f.expr()
            2*y^3 + x

        """
        from scalarfield import ScalarField
        result = ScalarField(self.chart.domain, name=name, 
                             latex_name=latex_name)
        result.express = {self.chart: self}
        result.coord_expression = {(self.chart, result.domain2.def_chart):
                                MultiFunctionChart(self.chart, self.express)}
        return result

 
#*****************************************************************************

class ZeroFunctionChart(FunctionChart):
    r"""
    Null function of coordinates belonging to a chart on a manifold. 

    INPUT:
    
    - ``chart`` -- the chart on which the null function is defined

    EXAMPLES:
    
    Null function defined on a 2-dimensional chart::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy.<x,y> = m.chart('x y')
        sage: f = ZeroFunctionChart(c_xy) ; f
        0
        sage: f.view()
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
            if other.chart != self.chart:
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
            if other.chart != self.chart:
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


    def scalar_field(self, name=None, latex_name=None):
        r""" 
        Return the zero scalar field on the domain covered by the chart on 
        which ``self`` is defined.
        
        INPUT: 
        
        - ``name`` -- (default: None) unused 
        - ``latex_name`` -- (default: None) unsued 
        
        OUTPUT:
        
        - instance of :class:`ZeroScalarField`
        
        EXAMPLES:

        Construction of a zero scalar field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')                  
            sage: c_xy.<x,y> = M.chart('x y')
            sage: fc = ZeroFunctionChart(c_xy)
            sage: f = fc.scalar_field() ; f
            zero scalar field on the 2-dimensional manifold 'M'
            sage: f.expr()
            0
        """
        return self.chart.domain.zero_scalar_field

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
    
    Each function `f_i` is stored as an instance of :class:`FunctionChart`.

    INPUT:
    
    - ``chart`` -- the chart defining the coordinates
    - ``*expressions`` -- the list of the coordinate expressions of the `m` 
      functions (`m\geq 1`)
    
    EXAMPLES: 
    
    A set of 3 functions of 2 coordinates::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy.<x,y>  = m.chart('x y') 
        sage: f = MultiFunctionChart(c_xy, x-y, x*y, cos(x)*exp(y))
        sage: f
        functions (x - y, x*y, cos(x)*e^y) on the chart (M, (x, y))
        sage: f.functions
        (x - y, x*y, cos(x)*e^y)
        sage: f(x,y)
        (x - y, x*y, cos(x)*e^y)
        sage: latex(f)
        \left(x - y, x y, \cos\left(x\right) e^{y}\right)        
    
    Each real-valued function `f_i` (`1\leq i \leq m`) composing `f` can be 
    accessed via the square-bracket operator, by providing `i-1` as an 
    argument::
    
        sage: f[0]
        x - y
        sage: f[1]
        x*y
        sage: f[2]
        cos(x)*e^y

    Each f[i-1] is an instance of :class:`FunctionChart`::
    
        sage: type(f[0])
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: f[0].view()
        (x, y) |--> x - y
        
    A MultiFunctionChart can contain a single function, although one should 
    rather employ the class :class:`FunctionChart` for this purpose::
    
        sage: g = MultiFunctionChart(c_xy, x*y^2)
        sage: g.functions
        (x*y^2,)
    
    Evaluating the functions at specified coordinates::
 
        sage: f(1,2)
        (-1, 2, cos(1)*e^2)
        sage: (a, b) = var('a b')
        sage: f(a,b)
        (a - b, a*b, cos(a)*e^b)
        sage: g(1,2)
        (4,)
        
    The Jacobian matrix::
    
        sage: f.jacobian()
        [[1, -1], [y, x], [-e^y*sin(x), cos(x)*e^y]]
        sage: g.jacobian()
        [[y^2, 2*x*y]]
    
    If the number of functions equals the number of coordinates, the Jacobian
    determinant can be evaluated::
    
        sage: h = MultiFunctionChart(c_xy, x-y, x*y)
        sage: h.jacobian_det()
        x + y
        
    """
    def __init__(self, chart, *expressions): 
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self.chart = chart
        self.nc = len(self.chart.xx)    # number of coordinates
        self.nf = len(expressions)      # number of functions
        self.functions = tuple(FunctionChart(chart, expressions[i]) for i in range(self.nf))
        self._jacob = None
        self._jacob_matrix = None
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
        
    def expr(self):
        r"""
        Return the symbolic expression of the image of the `m` functions, as
        
        .. MATH::
    
            (f_1(x^1,\ldots,x^n),\ldots, f_m(x^1,\ldots,x^n))
                
        OUTPUT:
        
        - tuple of symbolic expressions corresponding to the above formula
        
        EXAMPLES:
        
        A set of 3 functions of 2 coordinates::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y') 
            sage: f = MultiFunctionChart(c_xy, x-y, x*y, cos(x)*exp(y))
            sage: f.expr()
            (x - y, x*y, cos(x)*e^y)
            sage: type(f.expr()[0]) 
            <type 'sage.symbolic.expression.Expression'>
            sage: f.expr() == f(x,y)
            True

        """
        return tuple( self.functions[i].express for i in range(self.nf) )
        
    def copy(self):
        r"""
        Return an exact copy of ``self``.
        
        The derived quantities (Jacobian matrix) are not copied, because they 
        can be reconstructed if necessary.
        
        EXAMPLE:
        
        Copy of a set of 3 functions of 2 coordinates::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y') 
            sage: f = MultiFunctionChart(c_xy, x-y, x*y, cos(x)*exp(y))
            sage: g = f.copy() ; g
            functions (x - y, x*y, cos(x)*e^y) on the chart (M, (x, y))

        """
        return MultiFunctionChart(self.chart, *(self.expr()))

    def __getitem__(self, index):
        r""" 
        Return a specified function of the set represented by ``self``.
        
        INPUT:
        
        -- ``index`` -- index `i` of the function (`0\leq i \leq m-1`)
        
        OUTPUT
        
        -- instance of :class:`FunctionChart` representing the function
            
        """
        return self.functions[index]
        
    def __call__(self, *coords):
        r"""
        Compute the values of the functions at specified coordinates.
        
        INPUT:
        
        - ``*coords`` -- list of coordinates where the functions are to be
          evaluated 
        
        OUTPUT:
        
        - the values of the `m` functions.   
         
        """
        return tuple( self.functions[i](*coords) for i in range(self.nf) )

    def jacobian(self):
        r"""
        Return the Jacobian matrix of the system of functions.
        
        ``jacobian()`` is a 2-dimensional array of size `m\times n` 
        where `m` is the number of functions and `n` the number of coordinates, 
        the generic element being `J_{ij} = \frac{\partial f_i}{\partial x^j}` 
        with `1\leq i \leq m` (row index) and `1\leq j \leq n` (column index).
        
        Each `J_{ij}` is an instance of :class:`FunctionChart`.
        
        OUTPUT:
        
        - Jacobian matrix as a 2-dimensional array J of FunctionChart's, 
          J[i-1][j-1] being `J_{ij} = \frac{\partial f_i}{\partial x^j}`
          for `1\leq i \leq m` and `1\leq j \leq n`
          
        EXAMPLES:

        Jacobian of a set of 3 functions of 2 coordinates::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = MultiFunctionChart(c_xy, x-y, x*y, cos(x)*exp(y))
            sage: f.jacobian()
            [[1, -1], [y, x], [-e^y*sin(x), cos(x)*e^y]]
            sage: f.jacobian()[0][1]
            -1
            sage: type(f.jacobian()[0][1])
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: f.jacobian()[0][1].view()
            (x, y) |--> -1

        """
        from sage.matrix.constructor import matrix
        from sage.calculus.functional import diff
        if self._jacob is None:
            self._jacob = [[ FunctionChart(self.chart, 
                            simplify_chain(diff(self.functions[i].express, 
                                                self.chart.xx[j])) )
                    for j in range(self.nc) ] for i in range(self.nf) ]
            self._jacob_matrix = matrix( [[ self._jacob[i][j].express 
                    for j in range(self.nc) ] for i in range(self.nf) ] )
        return self._jacob
        
    def jacobian_det(self):
        r"""
        Return the Jacobian determinant of the system of functions.
        
        The number `m` of functions must equal the number `n` of 
        coordinates.
        
        OUTPUT:
        
        - instance of :class:`FunctionChart` representing the determinant
        
        EXAMPLE:
        
        Jacobian determinant of a set of 2 functions of 2 coordinates::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy.<x,y> = m.chart('x y')
            sage: f = MultiFunctionChart(c_xy, x-y, x*y)
            sage: f.jacobian_det()
            x + y
            
        The output of :meth:`jacobian_det` is an instance of 
        :class:`FunctionChart` and can therefore be called on specific values
        of the coordinates, e.g. (x,y)=(1,2)::
        
            sage: type(f.jacobian_det())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: f.jacobian_det()(1,2) 
            3
            
        """
        from utilities import simple_determinant
        if self._jacob_det is None: 
            if (self.nf != self.nc):
                raise ValueError("The Jacobian matrix is not square.")
            self.jacobian() # to force the computation of self._jacob_matrix
            #!# the following is a workaround for a bug in Sage (cf. trac ticket #14403)
            self._jacob_det = FunctionChart(self.chart, 
                       simplify_chain(simple_determinant(self._jacob_matrix)) )
            # the proper writing should be this:
            # self._jacob_det = FunctionChart(self.chart, simplify_chain(self._jacob_matrix.det()) )
        return self._jacob_det


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
        sage: c_spher.<r,th,ph> = m.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: c_cart.<x,y,z> = m.chart('x y z')        
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: ch
        coordinate change from chart (R3, (r, th, ph)) to chart (R3, (x, y, z))
        sage: latex(ch)
        (r, \theta, \phi) \mapsto (x, y, z)

    Each created coordinate change is automatically added to the manifold's 
    dictionary :attr:`coord_changes`; this dictionary is accessed via the method :meth:`Domain.coord_change`::    

        sage: m.coord_change(c_spher, c_cart)
        coordinate change from chart (R3, (r, th, ph)) to chart (R3, (x, y, z))
    
    It also generates a new entry in the manifold's dictionary 
    :attr:`frame_changes`, containing the relevant change-of-basis matrix; 
    this dictionary is accessed via the method :meth:`Manifold.frame_change`::

        sage: m.frame_change(c_cart.frame, c_spher.frame)
        field of tangent-space automorphisms on the 3-dimensional manifold 'R3'
        sage: m.frame_change(c_cart.frame, c_spher.frame)[:]
        [   cos(ph)*sin(th)  r*cos(ph)*cos(th) -r*sin(ph)*sin(th)]
        [   sin(ph)*sin(th)  r*cos(th)*sin(ph)  r*cos(ph)*sin(th)]
        [           cos(th)         -r*sin(th)                  0]
    
    The coordinate change can be called directly on a set of "old" coordinates 
    to get the "new" ones::
    
        sage: ch(1,pi/2,0)
        (1, 0, 0)
        
    The Jacobian matrix of the coordinate change::
    
        sage: ch.jacobian
        [[cos(ph)*sin(th), r*cos(ph)*cos(th), -r*sin(ph)*sin(th)], [sin(ph)*sin(th), r*cos(th)*sin(ph), r*cos(ph)*sin(th)], [cos(th), -r*sin(th), 0]]
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
        # If the two charts are on the same domain, the coordinate change is 
        # added to the domain (and superdomains) dictionary and the 
        # Jacobian matrix is added to the dictionary of changes of frame:
        if chart1.domain == chart2.domain:
            domain = chart1.domain
            for sdom in domain.superdomains:
                sdom.coord_changes[(chart1, chart2)] = self
            frame1 = chart1.frame
            frame2 = chart2.frame
            ch_basis = AutomorphismField(chart1.domain) 
            ch_basis.add_comp(frame1)[:, chart1] = self.jacobian
            ch_basis.add_comp(frame2)[:, chart1] = self.jacobian
            for sdom in domain.superdomains:
                sdom.frame_changes[(frame2, frame1)] = ch_basis
            if (frame1, frame2) not in domain.frame_changes:
                for sdom in domain.superdomains:
                    sdom.frame_changes[(frame1, frame2)] = ch_basis.inverse()

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
        return self.chart1._latex_coordinates() + r' \mapsto ' + \
                self.chart2._latex_coordinates()
    
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
            sage: c_xy.<x,y> = m.chart('x y')
            sage: c_uv.<u,v> = m.chart('u v')
            sage: ch_to_uv = CoordChange(c_xy, c_uv, (x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2)
            sage: m.coord_changes 
            {(chart (M, (x, y)), chart (M, (u, v))): coordinate change from chart (M, (x, y)) to chart (M, (u, v))}
            sage: ch_to_xy = ch_to_uv.inverse() ; ch_to_xy
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: ch_to_xy.transf                                                         
            functions (1/2*sqrt(3)*v + 1/2*u, -1/2*sqrt(3)*u + 1/2*v) on the chart (M, (u, v))
            sage: m.coord_changes # optional - dictionary_output
            {(chart (M, (u, v)), chart (M, (x, y))): coordinate change from chart (M, (u, v)) to chart (M, (x, y)), 
            (chart (M, (x, y)), chart (M, (u, v))): coordinate change from chart (M, (x, y)) to chart (M, (u, v))}
   
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
        x2 = [ SR.var('xxxx' + str(i), domain=coord_domain[i]) 
               for i in range(n2) ]
        equations = [ x2[i] == self.transf.functions[i].express 
                      for i in range(n2) ]
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
        if self.chart1.domain == self.chart2.domain:
            domain = self.chart1.domain
            frame1 = self.chart1.frame
            frame2 = self.chart2.frame
            fr_change12 = domain.frame_changes[(frame1,frame2)]
            fr_change21 = domain.frame_changes[(frame2,frame1)]
            for comp in fr_change12.components[frame1]._comp.values():
                comp.function_chart(self.chart1, from_chart=self.chart2)
            for comp in fr_change12.components[frame2]._comp.values():
                comp.function_chart(self.chart1, from_chart=self.chart2)
            for comp in fr_change21.components[frame1]._comp.values():
                comp.function_chart(self.chart2, from_chart=self.chart1)
            for comp in fr_change21.components[frame2]._comp.values():
                comp.function_chart(self.chart2, from_chart=self.chart1)

        return self._inverse


    def set_inverse(self, *transformations, **kwds):
        r"""
        Sets the inverse of the coordinate transformation. 
        
        This is usefull when the automatic computation via :meth:`inverse()`
        fails. 
        
        INPUT:
        
        - ``transformations`` -- the inverse transformations expressed as a 
          list of the expressions of the "old" coordinates in terms of the 
          "new" ones
        - ``kwds`` -- keyword arguments: only ``check=True`` (default) or
          ``check=False`` are meaningfull; it determines whether the provided 
          transformations are checked to be indeed the inverse coordinate
          transformations. 
          
        EXAMPLES:
         
        From Cartesian to spherical coordinates in the plane::
          
            sage: m = Manifold(2, 'M')
            sage: c_cart.<x,y> = m.chart('x y')
            sage: c_spher.<r,ph> = m.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = CoordChange(c_spher, c_cart, r*cos(ph), r*sin(ph))
            sage: spher_to_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))              
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: spher_to_cart.inverse()
            coordinate change from chart (M, (x, y)) to chart (M, (r, ph))
            sage: m.coord_changes  # random output order
            {(chart (M, (x, y)), chart (M, (r, ph))): coordinate change from chart (M, (x, y)) to chart (M, (r, ph)), (chart (M, (r, ph)), chart (M, (x, y))): coordinate change from chart (M, (r, ph)) to chart (M, (x, y))}
              
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
        if self.chart1.domain == self.chart2.domain:
            domain = self.chart1.domain
            frame1 = self.chart1.frame
            frame2 = self.chart2.frame
            fr_change12 = domain.frame_changes[(frame1,frame2)]
            fr_change21 = domain.frame_changes[(frame2,frame1)]
            for comp in fr_change12.components[frame1]._comp.values():
                comp.function_chart(self.chart1, from_chart=self.chart2)
            for comp in fr_change12.components[frame2]._comp.values():
                comp.function_chart(self.chart1, from_chart=self.chart2)
            for comp in fr_change21.components[frame1]._comp.values():
                comp.function_chart(self.chart2, from_chart=self.chart1)
            for comp in fr_change21.components[frame2]._comp.values():
                comp.function_chart(self.chart2, from_chart=self.chart1)
    
