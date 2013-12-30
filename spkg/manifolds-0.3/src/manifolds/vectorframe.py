r"""
Vector Frames

The class :class:`VectorFrame` implements vector frames on differentiable 
manifolds over `\RR`. By 'vector frame' it is meant a field on a manifold M that
provides, at each point p in M, a vector basis of the tangent space at p. 

A derived class of :class:`VectorFrame` is :class:`CoordBasis`; it regards the 
vector frames associated with a chart, i.e. the so-called coordinate bases. 

The vector frame duals, i.e. the coframes, are implemented via the class
:class:`CoFrame`. The derived class :class:`CoordCoFrame` is devoted to 
coframes deriving from a chart. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013): initial version

EXAMPLES:
    
    Setting a vector frame on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z', 'xyz-coord')
        sage: e = VectorFrame(m, 'e') ; e
        vector frame 'e' on the 3-dimensional manifold 'M'
        sage: latex(e)
        \left(e\right)

    The first frame defined on a manifold is its default frame; in the present
    case it is the coordinate basis defined when introducing the chart c_xyz::
    
        sage: m.default_frame()
        coordinate basis 'xyz-coord_b' (d/dx,d/dy,d/dz)
        
    The default frame can be changed via the method
    :meth:`Manifold.set_default_frame`::
    
        sage: m.set_default_frame(e)
        sage: m.default_frame()
        vector frame 'e' on the 3-dimensional manifold 'M'

    The elements of a vector frame are vector fields on the manifold::
    
        sage: e.vec
        (vector field 'e_0' on the 3-dimensional manifold 'M', vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M')   
        
    Each element can be accessed by its index::
    
        sage: e(0)
        vector field 'e_0' on the 3-dimensional manifold 'M'
            
    The index range depends on the starting index defined on the manifold::

        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z', 'xyz-coord')
        sage: e = VectorFrame(m, 'e')              
        sage: e.vec
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')
        sage: e(1), e(2), e(3)
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')
    
    Let us check that the vector fields e(i) are the frame vectors from
    their components w.r.t. to the frame 'e'::
    
        sage: e(1).comp('e')[:]
        [1, 0, 0]
        sage: e(2).comp('e')[:]
        [0, 1, 0]
        sage: e(3).comp('e')[:]
        [0, 0, 1]
    
    Defining a vector frame on a manifold automatically creates the dual 
    coframe, which bares the same name (here 'e')::
    
        sage: m.coframes
         {'xyz-coord_b': coordinate coframe 'xyz-coord_b' (dx,dy,dz), 'e': coframe 'e' on the 3-dimensional manifold 'M'}
        sage: f = m.coframes['e'] ; f
        coframe 'e' on the 3-dimensional manifold 'M'
   
    Each element of the coframe is a 1-form::
   
        sage: f(1), f(2), f(3)
        (1-form 'e^1' on the 3-dimensional manifold 'M',
        1-form 'e^2' on the 3-dimensional manifold 'M',
        1-form 'e^3' on the 3-dimensional manifold 'M')
        sage: latex(f(1)), latex(f(2)), latex(f(3))
        (e^1, e^2, e^3)

    Let us check that the coframe (e^i) is indeed the dual of the vector 
    frame (e_i)::
    
        sage: f(1)(e(1)) # the 1-form e^1 applied to the vector field e_1
        scalar field 'e^1(e_1)' on the 3-dimensional manifold 'M'
        sage: f(1)(e(1)).expr() # the explicit expression of e^1(e_1)
        1
        sage: f(1)(e(1)).expr(), f(1)(e(2)).expr(), f(1)(e(3)).expr()
        (1, 0, 0)
        sage: f(2)(e(1)).expr(), f(2)(e(2)).expr(), f(2)(e(3)).expr()
        (0, 1, 0)
        sage: f(3)(e(1)).expr(), f(3)(e(2)).expr(), f(3)(e(3)).expr()
        (0, 0, 1)
    
    The coordinate basis associated to spherical coordinates of the 
    sphere `S^2`::
    
        sage: m = Manifold(2, 'S^2', start_index=1)
        sage: c_spher.<th,ph> = m.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        sage: b = m.default_frame() ; b
        coordinate basis 'spher_b' (d/dth,d/dph)
        sage: b(1)
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b(2)
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'

    The orthonormal frame constructed from the coordinate basis::
    
        sage: change_basis = AutomorphismField(m)
        sage: change_basis[:] = [[1,0], [0, 1/sin(th)]]
        sage: e = b.new_frame(change_basis, 'e') ; e 
        vector frame 'e' on the 2-dimensional manifold 'S^2'
        sage: e(1)[:]
        [1, 0]
        sage: e(2)[:]
        [0, 1/sin(th)]
        
    The change-of-basis matrices::
    
        sage: m.frame_change('spher_b', 'e')        
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: m.frame_change('spher_b', 'e')[:]
        [        1         0]
        [        0 1/sin(th)]
        sage: m.frame_change('e', 'spher_b')
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: m.frame_change('e', 'spher_b')[:]
        [      1       0]
        [      0 sin(th)]

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
from domain import Domain

class VectorFrame(SageObject):
    r"""
    Class for vector frames on a differentiable manifold over `\RR`. 
    
    By 'vector frame', it is meant a field on a manifold M that provides, at 
    each point p in M, a vector basis of the tangent space at p. 
    
    For each instanciation of a vector frame, a coframe is automatically 
    created, as an instance of the class :class:`CoFrame`. 
    
    INPUT:
    
    - ``domain`` -- manifold domain on which the vector frame is defined
    - ``name`` -- name given to the vector frame (should be rather short)
    - ``latex_name`` -- (default: None) symbol to denote the vector frame; if 
      None, the value of ``name`` is used. The actual LaTeX name will be this 
      symbol enclosed into parentheses

    EXAMPLES:

    Setting a vector frame on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z', 'xyz-coord')
        sage: e = VectorFrame(m, 'e')
        sage: e
        vector frame 'e' on the 3-dimensional manifold 'M'
        sage: latex(e)
        \left(e\right)

    The LaTeX symbol can be specified::
    
        sage: e = VectorFrame(m, 'E', r"\epsilon")
        sage: latex(e)
        \left(\epsilon\right)

    
    """
    def __init__(self, domain, name, latex_name=None):
        from vectorfield import VectorField
        if not isinstance(domain, Domain):
            raise TypeError("The first argument must be a manifold domain.")
        self.manifold = domain.manifold
        self.domain = domain
        if name not in self.manifold.frames:
            self.name = name
        else:
            raise ValueError("The name '" + name + "' is already used for " + 
                             "another frame on " + str(self.manifold))
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        # The frame is added to the domain's dict. of frames, as well as to all 
        # the superdomains' dict. of frames; moreover the fist defined frame is 
        # considered as the default one
        for sd in self.domain.superdomains:
            sd.frames[self.name] = self
            if sd.def_frame is None: 
                sd.def_frame = self
        # Dual coframe:
        if self.name not in self.manifold.coframes:
            CoFrame(self.domain, self.name, self.latex_name)
        self.coframe = self.domain.coframes[self.name]

        vl = list()
        si = self.manifold.sindex
        n = self.manifold.dim
        for i in range(n):
            v_name = self.name + "_" + str(i+si)
            v_symb = self.latex_name + "_" + str(i+si)
            v = VectorField(self.domain, v_name, v_symb)
            for j in range(n):
                v.set_comp(self.name)[j+si] = 0
            v.set_comp(self.name)[i+si] = 1
            vl.append(v)
        self.vec = tuple(vl)
        # Derived quantities:
        self._structure_coef = None
        # Initialization of the set of frames that are restrictions of the
        # current frame to subdomains of the frame domain:
        self.subframes = set([self]) 
        # Initialization of the set of frames which the current frame is a 
        # restriction of:
        self.superframes = set([self]) 

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "vector frame"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)

        return description


    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return r"\left(" + self.latex_name + r"\right)"

    def __call__(self, index):
        r"""
        Returns the frame vector field corresponding to the index.
        
        INPUT:
        
        - ``index`` -- the index of the vector 

        """
        n = self.manifold.dim
        si = self.manifold.sindex
        i = index - si
        if i < 0 or i > n-1:
            raise ValueError("Index out of range: " +
                              str(i+si) + " not in [" + str(si) + "," +
                              str(n-1+si) + "]")
        return self.vec[i]

        
    def new_frame(self, change_of_basis, name, latex_name=None):
        r"""
        Define a new vector frame from the current one. 
        
        The new vector frame is defined on the same domain as ``self`` from
        a field of automorphisms. 
        
        INPUT:
        
        - ``change_of_basis`` -- instance of :class:`AutomorphismField`
          describing the automorphism `P` that relates the current frame 
          `(e_i)` (described by ``self``) to the new frame `(n_i)` according
          to `n_i = P(e_i)`
        - ``name`` -- name given to the new vector frame (should be rather 
          short)
        - ``latex_name`` -- (default: None) symbol to denote the new vector 
          frame; if None, the value of ``name`` is used. The actual LaTeX name 
          will be this symbol enclosed into parentheses
          
        OUTPUT:
        
        - the new frame `(n_i)`, as an instance of :class:`VectorFrame`
        
        EXAMPLES:
        
        Frame resulting from a pi/3-rotation in the Euclidean plane::
        
            sage: m = Manifold(2,'R^2')
            sage: c_xy = m.chart('x y', 'xy-coord')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: m.frame_changes
            {}
            sage: rot = AutomorphismField(m)
            sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: n = e.new_frame(rot, 'n')
            sage: n(0)[:]
            [1/2*sqrt(3), 1/2]
            sage: n(1)[:]
            [-1/2, 1/2*sqrt(3)]
            sage: m.frame_changes
            {('n', 'e'): field of tangent-space automorphisms on the 2-dimensional manifold 'R^2',
            ('e', 'n'): field of tangent-space automorphisms on the 2-dimensional manifold 'R^2'}
            sage: a =  m.frame_change('e','n')
            sage: a[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a == rot
            True
            sage: a is rot
            False
            sage: a.components
            {'e': 2-indices components w.r.t. the vector frame 'e',
            'n': 2-indices components w.r.t. the vector frame 'n'}
            sage: a.comp('n')[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a1 = m.frame_change('n','e')
            sage: a1[:]
            [1/2*sqrt(3)         1/2]
            [       -1/2 1/2*sqrt(3)]
            sage: a1 == rot.inverse()
            True
            sage: a1 is rot.inverse()
            False
            sage: e(0).comp('n')[:]
            [1/2*sqrt(3), -1/2]
            sage: e(1).comp('n')[:]
            [1/2, 1/2*sqrt(3)]
  
        """
        from rank2field import AutomorphismField
        if not isinstance(change_of_basis, AutomorphismField):
            raise TypeError("The argument change_of_basis must be an " +
                            "instance of AutomorphismField.")
        the_new_frame = VectorFrame(self.domain, name, latex_name)
        transf = change_of_basis.copy()
        inv_transf = change_of_basis.inverse().copy()
        si = self.manifold.sindex
        nsi = self.manifold.dim + si
        # Components of the new frame vectors in the old frame: 
        for i in range(si,nsi):
            for j in range(si,nsi):
                the_new_frame.vec[i-si].add_comp(self.name)[[j]] = \
                                                  transf.comp(self.name)[[j,i]]
        # Components of the new coframe 1-forms in the old coframe: 
        for i in range(si,nsi):
            for j in range(si,nsi):
                the_new_frame.coframe.form[i-si].add_comp(self.name)[[j]] = \
                                              inv_transf.comp(self.name)[[i,j]]
        # The components of the transformation and its inverse are the same in 
        # the two frames:
        for i in range(si,nsi):
            for j in range(si,nsi):
                transf.add_comp(name)[[i,j]] = transf.comp(self.name)[[i,j]]
                inv_transf.add_comp(name)[[i,j]] = \
                                              inv_transf.comp(self.name)[[i,j]]
        # Components of the old frame vectors in the new frame: 
        for i in range(si,nsi):
            for j in range(si,nsi):
                self.vec[i-si].add_comp(name)[[j]] = \
                                                   inv_transf.comp(name)[[j,i]]
        # Components of the old coframe 1-forms in the new coframe: 
        for i in range(si,nsi):
            for j in range(si,nsi):
                self.coframe.form[i-si].add_comp(name)[[j]] = \
                                                       transf.comp(name)[[i,j]]
        for sdom in self.domain.superdomains:
            sdom.frame_changes[(self.name, name)] = transf
            sdom.frame_changes[(name, self.name)] = inv_transf
        return the_new_frame
        
    def new_subframe(self, domain, name, latex_name=None):
        r"""
        Construct a subframe.
        
        If ``self`` is a vector frame defined on the domain U, a subframe
        is the restriction of ``self`` to a subdomain V of U.
        
        INPUT:
        
        - ``domain`` -- subdomain `V` of the current frame domain `U` 
        - ``name`` -- string containing the name given to the subframe (should 
          be rather short and must be unique among all frames defined on the 
          manifold)
        - ``latex_name`` -- (default: None) symbol to denote the vector frame; if 
          None, the value of ``name`` is used. The actual LaTeX name will be this 
          symbol enclosed into parentheses
        
        OUTPUT:
        
        - the subframe, as an instance of :class:`VectoFrame`. 

        """
        if not domain.is_subdomain(self.domain):
            raise TypeError("The argument 'domain' must be a subdomain of " + 
                            " the frame domain.")
        res = VectorFrame(domain, name, latex_name)
        # Update of superframes and subframes:
        res.superframes.update(self.superframes)
        for sframe in self.superframes:
            sframe.subframes.add(res)
        return res
    
    def structure_coef(self):
        r"""
        Evaluate the structure coefficients associated to the vector frame. 
        
        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}` 
        defined by 
        
        .. MATH::
            
            [e_i, e_j] = C^k_{\ \, ij} e_k
            
        OUPUT:
        
        - the structure coefficients `C^k_{\ \, ij}`, as an instance of 
          :class:`CompWithSym` with 3 indices ordered as `(k,i,j)`. 
          
        EXAMPLE:
        
        Structure coefficients of the orthonormal frame associated to
        spherical coordinates in the Euclidean space `R^3`::
        
            sage: m = Manifold(3, 'R^3', '\RR^3', start_index=1)
            sage: c_spher.<r,th,ph> = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: ch_basis = AutomorphismField(m) 
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: m.frames               
            {'spher_b': coordinate basis 'spher_b' (d/dr,d/dth,d/dph)}
            sage: e = m.frames['spher_b'].new_frame(ch_basis, 'e')
            sage: e(1)[:]  # components of e_1 in the manifold's default frame (d/dr, d/dth, d/dth)
            [1, 0, 0]
            sage: e(2)[:]
            [0, 1/r, 0]
            sage: e(3)[:]
            [0, 0, 1/(r*sin(th))]
            sage: c = e.structure_coef() ; c
            3-indices components w.r.t. the vector frame 'e', with antisymmetry on the index positions (1, 2)
            sage: c[:]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, -1/r, 0], [1/r, 0, 0], [0, 0, 0]],
             [[0, 0, -1/r], [0, 0, -cos(th)/(r*sin(th))], [1/r, cos(th)/(r*sin(th)), 0]]]
            sage: c[2,1,2]  # C^2_{12}
            -1/r
            sage: c[3,1,3]  # C^3_{13}
            -1/r
            sage: c[3,2,3]  # C^3_{23}
            -cos(th)/(r*sin(th))


        """
        from component import CompWithSym
        if self._structure_coef is None:
            self._structure_coef = CompWithSym(self, 3, antisym=(1,2))
            si = self.manifold.sindex
            nsi = si + self.manifold.dim
            for k in range(si,nsi):
                ce_k = self.coframe.form[k-si]
                for i in range(si, nsi):
                    e_i = self.vec[i-si]
                    for j in range(i+1, nsi):
                        e_j = self.vec[j-si]
                        self._structure_coef[[k,i,j]] = ce_k(e_j.lie_der(e_i))
        return self._structure_coef
            
#******************************************************************************

class CoordBasis(VectorFrame):
    r"""
    Class for coordinate bases on a differentiable manifold over `\RR`. 
    
    By 'coordinate basis', it is meant a vector frame on a manifold M that 
    is associated to a coordinate system (chart) on M. The name of the 
    coordinate basis is that of the chart with the extension "_b"
    
    INPUT:
    
    - ``chart`` -- the chart defining the coordinates

    EXAMPLES:

    The coordinate basis associated to spherical coordinates of the 
    sphere `S^2`::
    
        sage: m = Manifold(2, 'S^2', start_index=1)
        sage: m.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
        chart 'spher' (S^2, (th, ph))
        sage: b = m.default_frame()
        sage: b
        coordinate basis 'spher_b' (d/dth,d/dph)
        sage: b(1)
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b(2)
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'
        sage: latex(b)
        \left(\frac{\partial}{\partial \theta },\frac{\partial}{\partial \phi }\right)
 
    """
    def __init__(self, chart):
        from sage.misc.latex import latex
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self.chart = chart
        name = chart.name + "_b"
        VectorFrame.__init__(self, chart.domain, name)
        self.manifold.coframes[self.name] = CoordCoFrame(self.chart, 
                                                         self.name)
        self.coframe = self.domain.coframes[self.name]
        n = self.manifold.dim
        for i in range(n):
            self.vec[i].name = "d/d" + str(self.chart.xx[i])
            self.vec[i].latex_name = r"\frac{\partial}{\partial" + \
                                     latex(self.chart.xx[i]) + r"}"
        self.latex_name = ",".join([self.vec[i].latex_name for i in range(n)])


    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "coordinate basis '" + self.name + "' (" + \
          ",".join([self.vec[i].name for i in range(self.manifold.dim)]) + ")"
        return description

    def structure_coef(self):
        r"""
        Returns the structure coefficients associated to the vector frame. 
        
        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}` 
        defined by 
        
        .. MATH::
            
            [e_i, e_j] = C^k_{\ \, ij} e_k
        
        In the present case, where `(e_i)` is a coordinate basis, 
        `C^k_{\ \, ij}=0`. 
        
        OUPUT:
        
        - the structure coefficients `C^k_{\ \, ij}`, as a vanishing instance 
          of :class:`CompWithSym` with 3 indices ordered as `(k,i,j)`
          
        EXAMPLE:
        
        Structure coefficients of the coordinate basis associated to
        spherical coordinates in the Euclidean space `R^3`::
        
            sage: m = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_spher = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi', 'spher')
            sage: b = m.default_frame() ; b
            coordinate basis 'spher_b' (d/dr,d/dth,d/dph)
            sage: c = b.structure_coef() ; c
            3-indices components w.r.t. the coordinate basis 'spher_b' (d/dr,d/dth,d/dph), with antisymmetry on the index positions (1, 2)
            sage: c == 0
            True
        
        """
        from component import CompWithSym
        if self._structure_coef is None:
            self._structure_coef = CompWithSym(self, 3, antisym=(1,2))
            # A just created CompWithSym is zero
        return self._structure_coef
        


#******************************************************************************

class CoFrame(SageObject):
    r"""
    Class for coframes on a differentiable manifold over `\RR`. 
    
    By 'coframe', it is meant a n-tuple of 1-forms on a manifold M that 
    provides, at each point p in M, a basis of the space dual to the tangent 
    space at p. 
    
    INPUT:
    
    - ``domain`` -- the open domain on which the coframe is defined
    - ``name`` -- name given to the coframe; it must coincide with that of 
      the dual vector frame; if the latter has not been defined yet, it is 
      created by the CoFrame constructor. 
    - ``latex_name`` -- (default: None) symbol to denote the coframe; if 
      None, the value of ``name`` is used. The actual LaTeX name will be this 
      symbol enclosed into parentheses

    EXAMPLES:

    Coframe on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z', 'xyz')
        sage: e = CoFrame(m, 'e') ; e
        coframe 'e' on the 3-dimensional manifold 'M'
    
    The 1-forms composing the coframe are obtained via the () operator::
        
        sage: e(1), e(2), e(3)
        (1-form 'e^1' on the 3-dimensional manifold 'M',
         1-form 'e^2' on the 3-dimensional manifold 'M',
         1-form 'e^3' on the 3-dimensional manifold 'M')

    Defining a coframe automatically creates the dual vector frame, which bares
    the same name::
    
        sage: m.frames
        {'xyz_b': coordinate basis 'xyz_b' (d/dx,d/dy,d/dz), 'e': vector frame 'e' on the 3-dimensional manifold 'M'}
        sage: ev = m.frames['e'] ; ev
        vector frame 'e' on the 3-dimensional manifold 'M'
        sage: ev(1), ev(2), ev(3)
        (vector field 'e_1' on the 3-dimensional manifold 'M',
         vector field 'e_2' on the 3-dimensional manifold 'M',
         vector field 'e_3' on the 3-dimensional manifold 'M')

    Checking that e is the dual of ev::
    
        sage: e(1)(ev(1)).expr(), e(1)(ev(2)).expr(), e(1)(ev(3)).expr()
        (1, 0, 0)
        sage: e(2)(ev(1)).expr(), e(2)(ev(2)).expr(), e(2)(ev(3)).expr()
        (0, 1, 0)
        sage: e(3)(ev(1)).expr(), e(3)(ev(2)).expr(), e(3)(ev(3)).expr()
        (0, 0, 1)

    """
    def __init__(self, domain, name, latex_name=None):
        from diffform import OneForm
        if not isinstance(domain, Domain):
            raise TypeError("The first argument must be an open domain.")
        self.manifold = domain.manifold
        self.domain = domain
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        # The coframe is added to the domain's dict. of coframes, as well as to 
        # all the superdomains' dict. of frames
        for sd in self.domain.superdomains:
            sd.coframes[self.name] = self
        # Dual vector frame:
        if self.name not in self.manifold.frames:
            VectorFrame(self.domain, self.name, self.latex_name)
        self.frame = self.domain.frames[self.name]
        
        vl = list()
        si = self.manifold.sindex
        n = self.manifold.dim
        for i in range(n):
            v_name = self.name + "^" + str(i+si)
            v_symb = self.latex_name + "^" + str(i+si)
            v = OneForm(self.domain, v_name, v_symb)
            for j in range(n):
                v.set_comp(self.name)[j+si] = 0
            v.set_comp(self.name)[i+si] = 1
            vl.append(v)
        self.form = tuple(vl)
        
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "coframe"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)

        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return r"\left(" + self.latex_name + r"\right)"

    def __call__(self, index):
        r"""
        Returns the coframe 1-form corresponding to the index.
        
        INPUT:
        
        - ``index`` -- the index of the 1-form

        """
        n = self.manifold.dim
        si = self.manifold.sindex
        i = index - si
        if i < 0 or i > n-1:
            raise ValueError("Index out of range: " +
                              str(i+si) + " not in [" + str(si) + "," +
                              str(n-1+si) + "]")
        return self.form[i]


#******************************************************************************

class CoordCoFrame(CoFrame):
    r"""
    Class for coordinate coframes on a differentiable manifold over `\RR`. 
    
    By 'coordinate coframe', it is meant the n-tuple of the differentials of 
    the coordinates belonging to a chart on a manifold.
    
    INPUT:
    
    - ``chart`` -- the chart defining the coordinates
    - ``name`` -- name given to the coframe; it must coincide with that of 
      the dual vector frame; if the latter has not been defined yet, it is 
      created by the CoFrame constructor. 
    - ``latex_name`` -- (default: None) symbol to denote the coframe; if 
      None, the value of ``name`` is used. The actual LaTeX name will be this 
      symbol enclosed into parentheses

    EXAMPLES:
    
    Coordinate coframe on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z', 'xyz')
        sage: m.frames 
        {'xyz_b': coordinate basis 'xyz_b' (d/dx,d/dy,d/dz)}
        sage: m.coframes
        {'xyz_b': coordinate coframe 'xyz_b' (dx,dy,dz)}
        sage: dX = m.coframes['xyz_b'] ; dX
        coordinate coframe 'xyz_b' (dx,dy,dz)

    The 1-forms composing the coframe are obtained via the () operator::
    
        sage: dX(1)
        1-form 'dx' on the 3-dimensional manifold 'M'
        sage: dX(2)
        1-form 'dy' on the 3-dimensional manifold 'M'
        sage: dX(3)
        1-form 'dz' on the 3-dimensional manifold 'M'
        sage: dX(1)[:]
        [1, 0, 0]
        sage: dX(2)[:]
        [0, 1, 0]
        sage: dX(3)[:]
        [0, 0, 1]

    The coframe is the dual of the coordinate basis::
    
        sage: e = m.frames['xyz_b'] ; e
        coordinate basis 'xyz_b' (d/dx,d/dy,d/dz)
        sage: dX(1)(e(1)).expr(), dX(1)(e(2)).expr(), dX(1)(e(3)).expr()
        (1, 0, 0)
        sage: dX(2)(e(1)).expr(), dX(2)(e(2)).expr(), dX(2)(e(3)).expr()
        (0, 1, 0)
        sage: dX(3)(e(1)).expr(), dX(3)(e(2)).expr(), dX(3)(e(3)).expr()
        (0, 0, 1)
    
    Each 1-form of a coordinate coframe is closed::
    
        sage: dX(1).exterior_der()
        2-form 'ddx' on the 3-dimensional manifold 'M'
        sage: dX(1).exterior_der() == 0
        True

    """
    def __init__(self, chart, name=None, latex_name=None):
        from sage.misc.latex import latex
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        if name is None:
            name = chart.name + "_b"
        CoFrame.__init__(self, chart.domain, name, latex_name)

        self.chart = chart
        n = self.manifold.dim
        for i in range(n):
            self.form[i].name = "d" + str(self.chart.xx[i])
            self.form[i].latex_name = r"\mathrm{d}" + \
                                     latex(self.chart.xx[i])
        if latex_name is None:
            self.latex_name = \
                ",".join([self.form[i].latex_name for i in range(n)])

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "coordinate coframe '" + self.name + "' (" + \
          ",".join([self.form[i].name for i in range(self.manifold.dim)]) + ")"
        return description
