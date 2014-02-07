r"""
Differentiable manifolds

The class :class:`Manifold` implements differentiable manifolds over `\RR`. 

Ideally this class should inherit from a class describing topological 
manifolds or at least topological spaces. Since such classes do not
exist in Sage yet, the class :class:`Manifold` inherits from the 
class :class:`OpenDomain`. Via the latter, the class :class:`Manifold` inherits 
from the generic Sage class :class:`Parent` 
and is declared to belong to the category of sets (Sage category :class:`Sets`).
The corresponding Sage :class:`Element`'s are implemented via the class
:class:`Point`. 

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
    
        sage: stereoN.<x,y> = U.chart('x y') ; stereoN
        chart (U, (x, y))
        
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
        sage: stereoS.<u,v> = V.chart('u v') ; stereoS
        chart (V, (u, v))

    The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::
    
        sage: N = V.point((0,0), name='N') ; N
        point 'N' on 2-dimensional manifold 'S^2'

    At this stage, the manifold's atlas contains two charts::
    
        sage: M.atlas
        [chart (U, (x, y)), chart (V, (u, v))]

    To finalize things, we must declare the transition map between these two
    charts: calling W the intersection of U and V, (W is the subdomain of U 
    defined by `x^2+y^2\not=0`, as well as the subdomain of V defined by 
    `u^2+v^2\not=0`), we set::
    
        sage: trans = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)), 'W', \
                                             x^2+y^2!=0, u^2+v^2!=0)
        sage: trans
        coordinate change from chart (W, (x, y)) to chart (W, (u, v))
        sage: stereoN_W = trans.chart1
        sage: stereoS_W = trans.chart2
        
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
        [chart (U, (x, y)), chart (V, (u, v)), chart (W, (x, y)), chart (W, (u, v))]
         
    The first defined chart is considered as the default chart on the 
    manifold (unless it is changed by the method 
    :meth:`Domain.set_default_chart`)::
    
        sage: M.default_chart()
        chart (U, (x, y))

    This means that the chart can be omitted when specifying some point coordinates::
    
        sage: p = M.point((1,2))  # a point is created with coordinates (1,2)
        sage: p.coordinates # these coordinates refer to the default chart (and its subcharts if relevant):
        {chart (U, (x, y)): (1, 2), chart (W, (x, y)): (1, 2)}
        sage: p.coord() # if the chart is not specified, the default chart coordinates are returned:
        (1, 2)
        sage: p.coord(stereoS_W) # the coordinates in the chart stereoS_W are computed by means of the transition map:
        (1/5, 2/5)
        
    Manifolds are 'Parent' Sage objects, whose elements are the points::
    
        sage: p.parent()
        2-dimensional manifold 'S^2'
        sage: p in M
        True
        sage: p == M((1,2))
        True
        
    A manifold has a default vector frame, which, unless otherwise specified, 
    is the coordinate frame associated with the first defined chart::
    
        sage: M.default_frame()
        coordinate frame (U, (d/dx,d/dy))
        sage: latex(M.default_frame())
        \left(U ,\left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)\right)
        
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

from sage.structure.unique_representation import UniqueRepresentation
from domain import OpenDomain

class Manifold(OpenDomain):
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
    def __init__(self, n, name, latex_name=None, start_index=0):
        from sage.rings.integer import Integer
        if not isinstance(n, (int, Integer)):
            raise TypeError("The manifold dimension must be an integer.")
        if n<1:
            raise ValueError("The manifold dimension must be strictly " + 
                             "positive.")
        self.dim = n
        OpenDomain.__init__(self, self, name, latex_name)
        self.sindex = start_index
        self.domains = {self.name: self}
        
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
        [chart (field R, (x_realline,))]

    The instance is unique (singleton pattern)::
    
        sage: myRealLine = RealLineManifold()
        sage: myRealLine == RealLine
        True
        sage: myRealLine is RealLine
        True
        
    
    """
    def __init__(self):
        from chart import Chart
        Manifold.__init__(self, 1, name="field R", latex_name=r"\RR") 
        Chart(self, 'x_realline')

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
