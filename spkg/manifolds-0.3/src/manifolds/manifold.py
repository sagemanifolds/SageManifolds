r"""
Differentiable manifolds

The class :class:`Manifold` implements differentiable manifolds over `\RR`. 

Ideally this class should inherit from a class describing topological 
manifolds or at least topological spaces. Since such classes do not
exist in Sage yet, the class :class:`Manifold` inherits from the 
generic Sage class :class:`Parent` and is declared to belong to the 
category of sets (Sage category :class:`Sets`). 

The derived class :class:`RealLineManifold`  implements the real line `\RR` 
as a manifold of dimension one. The unique instance of it is :data:`RealLine`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES:
    
    The sphere `S^2` as a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'S^2')
        sage: M
        2-dimensional manifold 'S^2'
        
    Let us consider the complement of the North pole; it is an open domain
    of `S^2`, which we call U::
        
        sage: U = M.open_domain('U') ; U
        open domain 'U' on the 2-dimensional manifold 'S^2'
        
    A standard chart on U is provided by the stereographic projection from the
    North pole to the equatorial plane::
    
        sage: cU.<x,y> = U.chart('x y', 'stereo_N') ; cU
        chart 'stereo_N' (U, (x, y))
        
    Thanks to the operator <x,y> on the left-hand side, the coordinates 
    declared in a chart (here x and y), are accessible by their names; they are
    Sage's symbolic variables::
    
        sage: y
        y
        sage: type(y)
        <type 'sage.symbolic.expression.Expression'>

    The South pole is the point of coordinates `(x,y)=(0,0)` in the above
    chart::
    
        sage: S = U.point((0,0), name='S') ; S
        point 'S' on 2-dimensional manifold 'S^2'

    Let us call V the domain that is the complement of the South pole and let
    us introduce on it the chart induced by the stereographic projection from
    the South pole to the equatorial plane::
    
        sage: V = M.open_domain('V') ; V
        open domain 'V' on the 2-dimensional manifold 'S^2'
        sage: cV.<u,v> = V.chart('u v', 'stereo_S') ; cV
        chart 'stereo_S' (V, (u, v))

    The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::
    
        sage: N = V.point((0,0), name='N') ; N
        point 'N' on 2-dimensional manifold 'S^2'

    At this stage, the manifold's atlas has two charts::
    
        sage: M.atlas
        {'stereo_S': chart 'stereo_S' (V, (u, v)), 
         'stereo_N': chart 'stereo_N' (U, (x, y))}

    To finalize things, we must declare the transition map between these two
    charts: calling W the intersection of U and V, 'stereo_N_W' the 
    restriction of the chart 'stereo_N' to W (W being the subdomain of U 
    defined by `x^2+y^2\not=0`), 'stereo_S_W' the restriction of the chart 
    'stereo_S' to W (W being the subdomain of V defined by `u^2+v^2\not=0`), 
    we set::
    
        sage: trans = cU.transition_map(cV, (x/(x^2+y^2), y/(x^2+y^2)), 'W', \
                                        'stereo_N_W', x^2+y^2!=0, 'stereo_S_W', u^2+v^2!=0)
        sage: trans
        coordinate change from chart 'stereo_N_W' (W, (x, y)) to chart 'stereo_S_W' (W, (u, v))
        
    We also give the inverse of the transition map::
    
        sage: trans.set_inverse(u/(u^2+v^2), v/(u^2+v^2))
        Check of the inverse coordinate transformation:
           x == x
           y == y
           u == u
           v == v
           
    At this stage, we have four open domains on `S^2`::
        
        sage: M.domains
        {'U': open domain 'U' on the 2-dimensional manifold 'S^2', 
         'S^2': 2-dimensional manifold 'S^2', 
         'W': open domain 'W' on the 2-dimensional manifold 'S^2', 
         'V': open domain 'V' on the 2-dimensional manifold 'S^2'}

    W is the open domain that is the complement of the two poles::
    
        sage: W = M.domains['W'] ; W
        open domain 'W' on the 2-dimensional manifold 'S^2'
        sage: W is U.intersection(V)
        True
        sage: N in W
        False
        sage: S in W
        False

    The North pole lies in `V\setminus U` and the South pole in `U\setminus V`::
    
        sage: N in V, N in U
        (True, False)
        sage: S in U, S in V
        (True, False)

    Four charts have been defined on the manifold::
    
        sage: M.atlas
        {'stereo_S': chart 'stereo_S' (V, (u, v)), 
         'stereo_N_W': chart 'stereo_N_W' (W, (x, y)), 
         'stereo_S_W': chart 'stereo_S_W' (W, (u, v)), 
         'stereo_N': chart 'stereo_N' (U, (x, y))}
         
    The first defined chart is considered as the default chart on the 
    manifold (unless it is changed by the method 
    :meth:`Domain.set_default_chart`)::
    
        sage: M.default_chart()
        chart 'stereo_N' (U, (x, y))

    This means that the chart can be omitted when specifying some point coordinates::
    
        sage: p = M.point((1,2))  # a point is created with coordinates (1,2)
        sage: p.coordinates # these coordinates refer to the default chart (and its subcharts if relevant):
        {'stereo_N_W': (1, 2), 'stereo_N': (1, 2)}
        sage: p.coord() # if the chart is not specified, the default chart coordinates are returned:
        (1, 2)
        sage: p.coord('stereo_S_W') # the coordinates in the chart 'stereo_S_W' are computed by means of the transition map:
        (1/5, 2/5)
        
    Manifolds are 'Parent' Sage objects, whose elements are the points::
    
        sage: p.parent()
        2-dimensional manifold 'S^2'
        sage: p in M
        True
        sage: p == M((1,2))
        True
        
    A manifold has a default vector frame, which, unless otherwise specified, 
    is the coordinate basis associated with the first defined chart; it bears
    the name of the chart with the suffix '_b'::
    
        sage: M.default_frame()
        coordinate basis 'stereo_N_b' (d/dx,d/dy)
        sage: latex(M.default_frame())
        \left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)
        
    A manifold has a predefined zero scalar field, mapping all the points to 0; 
    it is an instance of :class:`ZeroScalarField`::
    
        sage: M.zero_scalar_field
        zero scalar field on the 2-dimensional manifold 'S^2'
        sage: M.zero_scalar_field(p)
        0
        sage: M.zero_scalar_field(N)
        0
        sage: M.zero_scalar_field(S)
        0
        
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

from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.structure.unique_representation import UniqueRepresentation
from domain import OpenDomain
from point import Point

class Manifold(OpenDomain, Parent):
    r"""  
    Base class for differentiable manifolds.
    
    This class implements differentiable manifolds over `\RR`. Ideally it 
    should inherit from a class describing topological manifolds, or at 
    least, topological spaces (not existing yet in Sage!). 
    
    INPUT:
    
    - ``n`` -- dimension of the manifold
    - ``name`` -- name given to the manifold 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the manifold; if 
      none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) lower bound of the range of indices on the
      manifold
    
    EXAMPLES:

    A 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M', r'\mathcal{M}')
        sage: m
        2-dimensional manifold 'M'
        sage: latex(m)
        \mathcal{M}
                
    The input parameter ``start_index`` becomes the attribute :attr:`sindex`
    of the manifold::
    
        sage: m = Manifold(4, 'M')  # default value of start_index is 0
        sage: m.sindex
        0
        sage: m = Manifold(4, 'M', start_index=1)
        sage: m.sindex
        1
        
    It defines the range of indices on the manifold::
    
        sage: m = Manifold(4, 'M')
        sage: list(m.irange())
        [0, 1, 2, 3]
        sage: m = Manifold(4, 'M', start_index=2)
        sage: list(m.irange())
        [2, 3, 4, 5]

    """
    
    Element = Point
    
    def __init__(self, n, name, latex_name=None, start_index=0):
        from sage.rings.integer import Integer
        if not isinstance(n, (int, Integer)):
            raise TypeError("The manifold dimension must be an integer.")
        if n<1:
            raise ValueError("The manifold dimension must be strictly " + 
                             "positive.")
        self.dim = n
        OpenDomain.__init__(self, self, name, latex_name)
        Parent.__init__(self, category=Sets())
        self.sindex = start_index
        self.domains = {self.name: self}
        
    def _element_constructor_(self, coords=None, chart_name=None, name=None, 
                 latex_name=None):
        r"""
        Construct a point on the manifold. 
        """
        return Point(self, coords, chart_name, name, latex_name)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return str(self.dim) + "-dimensional manifold '%s'" % self.name
    
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return self.latex_name

    def dimension(self):
        r"""
        Return the dimension of the manifold.
        
        EXAMPLE::
        
            sage: m = Manifold(2, 'M')
            sage: m.dimension()
            2

        """
        return self.dim


    def irange(self, start=None):
        r"""
        Single index generator.
                
        INPUT:
        
        - ``start`` -- (default: None) initial value of the index; if none is 
          provided, ``self.sindex`` is assumed

        OUTPUT:
        
        - an iterable index, starting from ``start`` and ending at
          ``self.sindex + self.dim -1``

        EXAMPLES:
        
        Index range on a 4-dimensional manifold::
        
            sage: m = Manifold(4, 'M')
            sage: for i in m.irange():
            ...       print i,
            ...     
            0 1 2 3
            sage: for i in m.irange(2):
            ...       print i,
            ...     
            2 3
            sage: list(m.irange())
            [0, 1, 2, 3]
    
        Index range on a 4-dimensional manifold with starting index=1::
        
            sage: m = Manifold(4, 'M', start_index=1)
            sage: for i in m.irange():              
            ...       print i,
            ...     
            1 2 3 4
            sage: for i in m.irange(2):             
            ...      print i,
            ...    
            2 3 4
        
        """
        si = self.sindex
        imax = self.dim + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1


    def index_generator(self, nb_indices):
        r"""
        Generator of index series.
        
        INPUT:

        - ``nb_indices`` -- number of indices in a series
        
        OUTPUT:
        
        - an iterable index series for a generic component with the specified
          number of indices

        EXAMPLES:
        
        Indices on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M', start_index=1)
            sage: for ind in m.index_generator(2):
            ...       print ind
            ...
            (1, 1)
            (1, 2)
            (2, 1)
            (2, 2)

        Loops can be nested::
        
            sage: for ind1 in m.index_generator(2):
            ...       print ind1, " : ",
            ...       for ind2 in m.index_generator(2):
            ...           print ind2,
            ...       print ""
            ...
            (1, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2) 
            (1, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2) 
            (2, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2) 
            (2, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2) 

        """
        si = self.sindex
        imax = self.dim - 1 + si
        ind = [si for k in range(nb_indices)]
        ind_end = [si for k in range(nb_indices)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(nb_indices-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

       
#******************************************************************************

class RealLineManifold(Manifold, UniqueRepresentation):
    r"""  
    Field of real numbers, as a manifold of dimension 1.
    
    This class is a of singleton type. 
    
    INPUT: 
    
    - None
    
    EXAMPLES:
                
    The pre-defined instance is :data:`RealLine`::
    
        sage: RealLine
        field R of real numbers
        sage: latex(RealLine)
        \RR
        sage: type(RealLine)
        <class 'sage.geometry.manifolds.manifold.RealLineManifold_with_category'>

    It is a manifold endoved with a canonical chart::
    
        sage: isinstance(RealLine, Manifold)
        True
        sage: RealLine.dimension()
        1
        sage: RealLine.atlas
        {'canonical': chart 'canonical' (field R, (x_realline,))}
    

    The instance is unique (singleton pattern)::
    
        sage: myRealLine = RealLineManifold()
        sage: myRealLine == RealLine
        True
        sage: myRealLine is RealLine
        True
        
    
    """
    def __init__(self):
        from chart import Chart
        Manifold.__init__(self, 1, name="field R", 
                          latex_name=r"\RR") 
        Chart(self, 'x_realline', 'canonical')

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return "field R of real numbers"

r"""
.. :data:: 
        The field of real numbers as the single instance of 
        :class:`RealLineManifold`.
"""
RealLine = RealLineManifold() 
