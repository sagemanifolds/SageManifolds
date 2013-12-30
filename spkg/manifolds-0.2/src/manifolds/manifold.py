r"""
Differentiable manifolds

The class :class:`Manifold` implements differentiable manifolds over `\RR`. 
Ideally it should inherit from a class describing topological spaces (not 
existing yet in Sage!). 

The derived class :class:`RealLineManifold` is a singleton class that 
implements the real line `\RR` as a manifold of dimension one. The unique 
instance of it is :data:`RealLine`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES:
    
    A 2-dimensional manifold::
    
        sage: m = Manifold(2, 'M')
        sage: m
        2-dimensional manifold 'M'

    By default, the LaTeX symbol representing the manifold is deduced from the 
    manifold given name, but it can be specified explicitely::
 
        sage: latex(m)
        M
        sage: m = Manifold(2, 'M', r'\mathcal{M}')
        sage: latex(m)
        \mathcal{M}

    An unnamed 4-dimensional manifold::

        sage: m = Manifold(4)
        sage: m
        4-dimensional manifold
        
    Setting a chart on a manifold, for instance spherical coordinates on `\RR^3`::

        sage: m = Manifold(3, 'R3', r'\mathcal{M}')
        sage: c_spher = Chart(m, r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        sage: c_spher
        chart 'spher' (r, th, ph)
        sage: m.get_atlas()
        {'spher': chart 'spher' (r, th, ph)}
        
    A defined coordinate is directly accessible by its name::

        sage: th
        th
        sage: type(th)
        <type 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        \theta
        
    Setting a second chart on the manifold, for instance Cartesian coordinates::
    
        sage: c_cart = Chart(m, 'x y z', 'cart')
        sage: c_cart
        chart 'cart' (x, y, z)
        sage: m.get_atlas()
        {'spher': chart 'spher' (r, th, ph), 'cart': chart 'cart' (x, y, z)}
       
    A manifold has a default chart, which, unless changed by the method
    :meth:`set_default_chart`, is the first defined chart::
    
        sage: m.default_chart()
        chart 'spher' (r, th, ph) 
        
    A manifold has a default vector frame, which, unless otherwise specified, 
    is the coordinate basis associated with the first defined chart; it bears
    the name of the chart with the suffix '_b'::
    
        sage: m.default_frame()
        coordinate basis 'spher_b' (d/dr,d/dth,d/dph)
        sage: latex(m.default_frame())
        \left(\frac{\partial}{\partial r },\frac{\partial}{\partial \theta },\frac{\partial}{\partial \phi }\right)
     
    Defining changes of coordinates::
    
        sage: ch = CoordChange(c_spher, c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: ch
        coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z)
        sage: latex(ch)
        (r, \theta, \phi) \mapsto (x, y, z)
        sage: m.coord_changes
        {('spher', 'cart'): coordinate change from chart 'spher' (r, th, ph) to chart 'cart' (x, y, z)}
        
    Defining points on the manifold::
        
        sage: p = Point(m, (1, pi/2, 0), 'spher')
        sage: p1 = Point(m, (1, pi/2, 0))  # using the default chart
        sage: p == p1
        True
        
    A manifold has a predefined zero scalar field, mapping all the points to 0; 
    it is an instance of :class:`ZeroScalarField`::
    
        sage: m.zero_scalar_field
        zero scalar field on the 3-dimensional manifold 'R3'
        sage: m.zero_scalar_field(p)
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

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

class Manifold(SageObject):
    r"""  
    Base class for differentiable manifolds.
    
    This class implements differentiable manifolds over `\RR`. Ideally it 
    should inherit from a class describing topological spaces (not existing yet 
    in Sage!). 
    
    INPUT:
    
    - ``n`` -- dimension of the manifold
    - ``name`` -- (default: None) name given to the manifold 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the manifold; if 
      none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) lower bound of the range of indices on the
      manifold
    
    EXAMPLES:

    A 2-dimensional manifold::
    
        sage: m = Manifold(2, 'test manifold', r'\mathcal{M}')
        sage: m
        2-dimensional manifold 'test manifold'
        sage: latex(m)
        \mathcal{M}
        
    An unnamed 4-dimensional manifold::

        sage: m = Manifold(4)
        sage: m
        4-dimensional manifold
        
    The input parameter ``start_index`` becomes the attribute :attr:`sindex`
    of the manifold::
    
        sage: m = Manifold(4)  # default value of start_index is 0
        sage: m.sindex
        0
        sage: m = Manifold(4, start_index=1)
        sage: m.sindex
        1
        
    It defines the range of indices on the manifold::
    
        sage: m = Manifold(4)
        sage: list(m.irange())
        [0, 1, 2, 3]
        sage: m = Manifold(4, start_index=2)
        sage: list(m.irange())
        [2, 3, 4, 5]
        
    """
    def __init__(self, n, name=None, latex_name=None, start_index=0):
        from sage.rings.integer import Integer
        from scalarfield import ZeroScalarField
        if not isinstance(n, (int, Integer)):
            raise TypeError("The manifold dimension must be an integer.")
        if n<1:
            raise ValueError("The manifold dimension must be strictly " + 
                             "positive.")
        self.dim = n
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.sindex = start_index
        self.atlas = {}     # the dictionary of charts defined on the manifold
        self.def_chart = None
        self.coord_changes = {}
        self.frames = {}  # the dict. of vector frames defined on the manifold
        self.def_frame = None
        self.frame_changes = {}
        self.coframes = {}  # the dict. of coframes defined on the manifold
        # The zero scalar field is constructed: 
        if self.name != 'field R':  
            #!# to avoid circular import of RealLine
            self.zero_scalar_field = ZeroScalarField(self)

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{no symbol}'
        else:
           return self.latex_name

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = str(self.dim) + "-dimensional manifold"
        if self.name is not None:
            description += " '%s'" % self.name
        return description
    
    def dimension(self):
        r"""
        Return the dimension of the manifold.
        
        EXAMPLE::
        
            sage: m = Manifold(2, 'M')
            sage: m.dimension()
            2

        """
        return self.dim
    
    def chart(self, chart_name):
        r"""
        Return a chart defined on the manifold.
        
        INPUT:
        
        - ``chart_name`` -- name of the chart
        
        OUTPUT:
        
        - instance of :class:`Chart` representing the chart of the given name
        
        EXAMPLES:
        
        Charts on a 2-dimensional manifold::
            
            sage: m = Manifold(2, 'M')
            sage: Chart(m, 'x y', 'xy')
            chart 'xy' (x, y)
            sage: Chart(m, 'u v', 'uv')
            chart 'uv' (u, v)
            sage: m.chart('xy')
            chart 'xy' (x, y)
            sage: m.chart('uv')
            chart 'uv' (u, v)

        The first defined chart is the default one::
                    
            sage: m.default_chart() is m.chart('xy')
            True

        """
        if chart_name not in self.atlas:
            raise TypeError("The chart '" + chart_name + "' has not been " + 
                            "defined on the " + repr(self))
        return self.atlas[chart_name]
        
    def default_chart(self):
        r"""
        Return the default chart defined on the manifold. 
        
        Unless changed via :meth:`set_default_chart`, the default chart is the 
        first one defined on the manifold. 
        
        OUTPUT:
        
        - instance of :class:`Chart` representing the default chart.
        
        EXAMPLES:
                    
        Default chart on a 2-dimensional manifold::
            
            sage: m = Manifold(2, 'M')
            sage: Chart(m, 'x y', 'xy')
            chart 'xy' (x, y)
            sage: Chart(m, 'u v', 'uv')
            chart 'uv' (u, v)
            sage: m.default_chart()
            chart 'xy' (x, y)

        """
        return self.def_chart
        
    def set_default_chart(self, chart):
        r"""
        Changing the default chart on the manifold.
        
        INPUT:
    
        - ``chart`` -- a chart (must be defined on the current manifold)

        EXAMPLES:
                    
        Charts on a 2-dimensional manifold::
            
            sage: m = Manifold(2, 'M')
            sage: Chart(m, 'x y', 'xy')
            chart 'xy' (x, y)
            sage: Chart(m, 'u v', 'uv')
            chart 'uv' (u, v)
            sage: m.default_chart()
            chart 'xy' (x, y)
            sage: m.set_default_chart(m.chart('uv'))
            sage: m.default_chart()
            chart 'uv' (u, v)

        """
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError(str(chart) + " is not a chart.")
        if chart.manifold != self:
            raise ValueError("The chart must be defined on the manifold.")
        self.def_chart = chart

    def get_atlas(self):
        r"""
        Return the manifold's atlas.
        
        OUTPUT:
        
        - dictionary of the various charts defined on the manifold
        
        EXAMPLE:
        
        Atlas of a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: m.get_atlas()
            {'xy': chart 'xy' (x, y), 'uv': chart 'uv' (u, v)}
            sage: print type(m.get_atlas())
            <type 'dict'>
            sage: m.chart('uv') is m.get_atlas()['uv']
            True
            
        """
        return self.atlas

    def frame(self, frame_name):
        r"""
        Return a vector frame defined on the manifold.
        
        By 'vector frame' it is meant a field on the manifold that provides, 
        at each point p, a vector basis of the tangent space at p.
        
        INPUT:
        
        - ``frame_name`` -- name of the vector frame
        
        OUTPUT:
        
        - instance of :class:`VectorFrame` representing the frame of the given 
          name
        
        EXAMPLES::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: e = VectorFrame(m, 'e')
            sage: m.frame('e')
            vector frame 'e' on the 2-dimensional manifold 'M'

        """
        if frame_name not in self.frames:
            raise TypeError("The vector frame '" + frame_name + "' has not " + 
                            "been defined on the " + repr(self))
        return self.frames[frame_name]

    def default_frame(self):
        r"""
        Return the default vector frame defined on the manifold. 
        
        By 'vector frame' it is meant a field on the manifold that provides, 
        at each point p, a vector basis of the tangent space at p.

        Unless changed via :meth:`set_default_frame`, the default frame is the 
        first one defined on the manifold, usually implicitely as the coordinate
        basis associated with the first chart defined on the manifold. 
        
        OUTPUT:
        
        - instance of :class:`VectorFrame`  representing the default vector 
          frame.
        
        EXAMPLES:
                    
        The default vector frame is often the coordinate basis associated
        with the first chart defined on the manifold; it bears the name of
        the chart with the suffix '_b'::
            
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: m.default_frame()
            coordinate basis 'xy_b' (d/dx,d/dy)
            sage: m.default_frame() is m.frame('xy_b')
            True

        """
        return self.def_frame

    def set_default_frame(self, frame):
        r"""
        Changing the default vector frame on the manifold.
        
        INPUT:
    
        - ``frame`` -- a vector frame (instance of :class:`VectorFrame`) defined 
          on the current manifold
          
        EXAMPLE::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: e = VectorFrame(m, 'e')
            sage: m.default_frame()                   
            coordinate basis 'xy_b' (d/dx,d/dy)
            sage: m.set_default_frame(e)
            sage: m.default_frame()     
            vector frame 'e' on the 2-dimensional manifold 'M'

        """
        from vectorframe import VectorFrame
        if not isinstance(frame, VectorFrame):
            raise TypeError(str(frame) + " is not a vector frame.")
        if frame.manifold != self:
            raise TypeError("The frame must be defined on the manifold.")
        self.def_frame = frame

    def coframe(self, coframe_name):
        r"""
        Return a coframe defined on the manifold.
        
        By 'coframe', it is meant a n-tuple of 1-forms on a manifold that 
        provides, at each point p, a basis of the space dual to the tangent 
        space at p.
        
        INPUT:
        
        - ``coframe_name`` -- name of the coframe (the same name as the dual 
          vector frame)
        
        OUTPUT:
        
        - instance of :class:`CoFrame` representing the coframe of the given 
          name
        
        EXAMPLES:
        
        Coframes on a 2-dimensional manifold::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: e = VectorFrame(m, 'e')
            sage: m.coframe('e')
            coframe 'e' on the 2-dimensional manifold 'M'
            sage: m.coframe('xy_b')
            coordinate coframe 'xy_b' (dx,dy)
            
        Check that the coframe 'e' is the dual of the vector frame 'e'::
        
            sage: f = m.coframe('e')
            sage: f(0)
            1-form 'e^0' on the 2-dimensional manifold 'M'
            sage: e(0)      
            vector field 'e_0' on the 2-dimensional manifold 'M'
            sage: f(0)(e(0)).expr()
            1
            sage: f(0)(e(1)).expr()
            0
            sage: f(1)(e(0)).expr()
            0
            sage: f(1)(e(1)).expr()
            1

        """
        if coframe_name not in self.coframes:
            raise TypeError("The coframe '" + frame_name + "' has not " + 
                            "been defined on the " + repr(self))
        return self.coframes[coframe_name]

    def coord_change(self, chart_name1, chart_name2):
        r"""
        Return a change of coordinates (transition map) defined on the manifold.
        
        INPUT:
        
        - ``chart_name1`` -- name of chart 1
        - ``chart_name2`` -- name of chart 2
        
        OUTPUT:
        
        - instance of :class:`CoordChange` representing the transition map from 
          chart 1 to chart 2 
        
        EXAMPLES:
        
        Change of coordinates on a 2-dimensional manifold::
            
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: CoordChange(c_xy, c_uv, x+y, x-y)
            coordinate change from chart 'xy' (x, y) to chart 'uv' (u, v)
            sage: m.coord_change('xy', 'uv')
            coordinate change from chart 'xy' (x, y) to chart 'uv' (u, v)

        """
        if (chart_name1, chart_name2) not in self.coord_changes:
            raise TypeError("The change of coordinates from '" + chart_name1 + 
                            "' to '" + chart_name2 + "' has not been " + 
                            "defined on the " + repr(self))
        return self.coord_changes[(chart_name1, chart_name2)]


    def frame_change(self, frame_name1, frame_name2):
        r"""
        Return a change of vector frames defined on the manifold.
                
        INPUT:
        
        - ``frame_name1`` -- name of vector frame 1
        - ``frame_name2`` -- name of vector frame 2
        
        OUTPUT:
        
        - instance of :class:`AutomorphismField` representing, at each point,
          the vector space automorphism `P` that relates frame 1, `(e_i)` say, 
          to frame 2, `(n_i)` say, according to `n_i = P(e_i)`
        
        EXAMPLES:
        
        Change of vector frames induced by a change of coordinates::
        
            sage: m = Manifold(2, 'M')
            sage: c_xy = Chart(m, 'x y', 'xy')
            sage: c_uv = Chart(m, 'u v', 'uv')
            sage: CoordChange(c_xy, c_uv, x+y, x-y)
            coordinate change from chart 'xy' (x, y) to chart 'uv' (u, v)
            sage: m.frame_change('xy_b', 'uv_b')
            field of tangent-space automorphisms on the 2-dimensional manifold 'M'
            sage: m.frame_change('xy_b', 'uv_b')[:]
            [ 1/2  1/2]
            [ 1/2 -1/2]
            sage: m.frame_change('uv_b', 'xy_b')
            field of tangent-space automorphisms on the 2-dimensional manifold 'M'
            sage: m.frame_change('uv_b', 'xy_b')[:]
            [ 1  1]
            [ 1 -1]
            sage: m.frame_change('uv_b', 'xy_b') ==  m.frame_change('xy_b', 'uv_b').inverse()
            True            

        """
        if (frame_name1, frame_name2) not in self.frame_changes:
            raise TypeError("The change of frame from '" + frame_name1 + 
                            "' to '" + frame_name2 + "' has not been " + 
                            "defined on the " + repr(self))
        return self.frame_changes[(frame_name1, frame_name2)]



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
        
            sage: m = Manifold(4)
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
        
            sage: m = Manifold(4, start_index=1)
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
        
            sage: m = Manifold(2, start_index=1)
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

class RealLineManifold(Manifold,UniqueRepresentation):
    r"""  
    Field of real numbers, as a manifold. 
    
    This class, which inherits from :class:`UniqueRepresentation`, is a of
    singleton type. 
    
    INPUT: 
    
    - None
    
    EXAMPLES:
                
    The pre-defined instance is :data:`RealLine`::
    
        sage: RealLine
        field R of real numbers
        sage: latex(RealLine)
        \RR
        sage: type(RealLine)
        <class 'sage.geometry.manifolds.manifold.RealLineManifold'>

    It is a manifold endoved with a canonical chart::
    
        sage: isinstance(RealLine, Manifold)
        True
        sage: RealLine.dim
        1
        sage: RealLine.atlas
        {'canonical': chart 'canonical' (x_realline,)}
    

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
