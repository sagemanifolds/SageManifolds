r"""
Rank-2 tensor fields

Three derived classes of :class:`TensorField` devoted to rank-2 tensor
are implemented:


* :class:`EndomorphismField` for fields of endomorphisms 
  (type (1,1) tensor fields)
* :class:`AutomorphismField` for fields of invertible endomorphisms
* :class:`IdentityMap` for the identity map on tangent spaces
* :class:`SymBilinFormField` for fields of symmetric bilinear forms 
  (type (0,2) tensor fields)

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

from tensorfield import TensorField
from vectorfield import VectorField

class SymBilinFormField(TensorField):
    r"""
    Class for fields of symmetric bilinear forms on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the field of symmetric bilinear 
      forms  is defined (must be an instance of class :class:`Domain`)
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A field of symmetric bilinear forms on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z')
        sage: t = SymBilinFormField(m, 'T') ; t
        field of symmetric bilinear forms 'T' on the 3-dimensional manifold 'M'
    
    Such a object is a tensor field of rank 2 and type (0,2)::
    
        sage: t.rank
        2
        sage: t.tensor_type
        (0, 2)

    The LaTeX symbol is deduced from the name or can be specified when creating
    the object::

        sage: latex(t)
        T
        sage: om = SymBilinFormField(m, 'Omega', r'\Omega')
        sage: latex(om)
        \Omega

    Components with respect to some vector frame::
        
        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: t.set_comp()
        fully symmetric 2-indices components w.r.t. the vector frame (M, (e_0,e_1,e_2))
        sage: type(t.comp())
        <class 'sage.geometry.manifolds.component.CompFullySym'>
           
    For the domain's default frame, the components are accessed with the square brackets::

        sage: t[0,0], t[0,1], t[0,2] = (1, 2, 3)
        sage: t[1,1], t[1,2] = (4, 5)
        sage: t[2,2] = 6
            
    The other components are deduced by symmetry::
        
        sage: t[1,0], t[2,0], t[2,1]
        (2, 3, 5)
        sage: t[:]
        [1 2 3]
        [2 4 5]
        [3 5 6]
       
    A symmetric bilinear form acts on vector pairs::
    
        sage: m = Manifold(2, 'M')
        sage: c_xy.<x,y> = m.chart('x y')
        sage: t = SymBilinFormField(m, 'T')
        sage: t[0,0], t[0,1], t[1,1] = (-1, x, y*x)
        sage: v1 = VectorField(m, 'V_1')
        sage: v1[:] = (y,x)  
        sage: v2 = VectorField(m, 'V_2')
        sage: v2[:] = (x+y,2)
        sage: s = t(v1,v2) ; s
        scalar field 'T(V_1,V_2)' on the 2-dimensional manifold 'M'
        sage: s.expr()
        x^3 + (3*x^2 + x)*y - y^2
        sage: s.expr() - t[0,0]*v1[0]*v2[0] - t[0,1]*(v1[0]*v2[1]+v1[1]*v2[0]) - t[1,1]*v1[1]*v2[1]
        0
        sage: latex(s)
        T\left(V_1,V_2\right)
    
    Adding two symmetric bilinear forms results in another symmetric bilinear
    form::

        sage: a = SymBilinFormField(m)          
        sage: a[0,0], a[0,1], a[1,1] = (1,2,3)  
        sage: b = SymBilinFormField(m)          
        sage: b[0,0], b[0,1], b[1,1] = (-1,4,5)
        sage: s = a + b ; s
        field of symmetric bilinear forms on the 2-dimensional manifold 'M'
        sage: s[:]
        [0 6]
        [6 8]

    But adding a symmetric bilinear from with a non-symmetric bilinear form 
    results in a generic type (0,2) tensor::
    
        sage: c = TensorField(m,0,2)            
        sage: c[:] = [[-2, -3], [1,7]]            
        sage: s1 = a + c ; s1
        tensor field of type (0,2) on the 2-dimensional manifold 'M'
        sage: s1[:]
        [-1 -1]
        [ 3 10]
        sage: s2 = c + a ; s2
        tensor field of type (0,2) on the 2-dimensional manifold 'M'
        sage: s2[:] 
        [-1 -1]
        [ 3 10]
        
    """
    def __init__(self, domain, name=None, latex_name=None):
        TensorField.__init__(self, domain, 0, 2, name, latex_name, sym=(0,1))

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "field of symmetric bilinear forms"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description
        
    def _new_comp(self, frame): 
        r"""
        Create some components in the given frame. 
                
        """
        from component import CompFullySym
        return CompFullySym(frame, 2)

    def _new_instance(self):
        r"""
        Create a :class:`SymBilinFormField` instance. 
        
        """
        return SymBilinFormField(self.domain)


#******************************************************************************

class EndomorphismField(TensorField):
    r"""
    Class for fields of endomorphisms (i.e. linear operators in each tangent 
    space) on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the field of endomorphisms is 
      defined (must be an instance of class :class:`Domain`)
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A field of endomorphisms on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z')
        sage: t = EndomorphismField(m, 'T') ; t
        field of endomorphisms 'T' on the 3-dimensional manifold 'M'
        
    A field of endomorphisms is a tensor field of rank 2 and of type (1,1)::
    
        sage: t.rank
        2
        sage: t.tensor_type
        (1, 1)
    
    Components with respect to a given frame::

        sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
        sage: t[:] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: t[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: change_basis = AutomorphismField(m)
        sage: change_basis[:] = [[-1,0,1], [0,1,2], [-2,1,2]]
        sage: f = e.new_frame(change_basis, 'f')
        sage: t.comp(f) # computation of t components in the frame f, from those in the frame e
        2-indices components w.r.t. the vector frame (M, (f_1,f_2,f_3))
        sage: t.comp(f)[:]
        [  9/2    -3 -15/2]
        [  -11     7    19]
        [ -5/2     2   7/2]
        sage: t.comp(f)[1,1]
        9/2

    An endomorphism maps a vector to a vector::
    
        sage: v = VectorField(m, 'v')
        sage: v[:] = (1,2,3)
        sage: w = t(v) ; w
        vector field 'T(v)' on the 3-dimensional manifold 'M'
        sage: w[:]
        [14, 32, 50]
        sage: t[:] * matrix([v[i].expr() for i in range(1,4)]).transpose()  # check:
        [14]
        [32]
        [50]
        sage: latex(t(v))
        T\left(v\right)

    """
    def __init__(self, domain, name=None, latex_name=None):
        TensorField.__init__(self, domain, 1, 1, name, latex_name)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "field of endomorphisms"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description
        
    def _new_instance(self):
        r"""
        Create a :class:`EndomorphismField` instance. 
        
        """
        return EndomorphismField(self.domain)
        
    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`TensorField.__call__` to allow for a single 
        vector argument. 
        """
        if len(arg) > 1:
            # the endomorphism acting as a type (1,1) tensor on a pair 
            # (1-form, vector), returning a scalar:
            return TensorField.__call__(self, *arg) 
        # the endomorphism acting as such, on a vector, returning a vector:
        vector = arg[0]
        if not isinstance(vector, VectorField):
            raise TypeError("The argument must be a vector field.")
        frame = self.common_frame(vector)
        t = self.components[frame]
        v = vector.components[frame]
        manif = self.manifold
        result = VectorField(self.domain)
        n = manif.dim
        si = manif.sindex
        for i in range(si, si+n):
            res = 0
            for j in range(si, si+n):
                res += t[[i,j]]*v[[j]]
            result.set_comp(frame)[i] = res
        # Name of the output:
        result.name = None
        if self.name is not None and vector.name is not None:
            result.name = self.name + "(" + vector.name + ")"
        # LaTeX symbol for the output:
        result.latex_name = None
        if self.latex_name is not None and vector.latex_name is not None:
            result.latex_name = self.latex_name + r"\left(" + \
                              vector.latex_name + r"\right)"
        return result


#******************************************************************************

class AutomorphismField(EndomorphismField):
    r"""
    Class for fields of linear automorphisms (i.e. invertible linear operators 
    in each tangent space) on a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the field of automorphisms is 
      defined (must be an instance of class :class:`Domain`)
    - ``name`` -- (default: None) name given to the field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the field; 
      if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A pi/3-rotation in the Euclidean 2-plane::
    
        sage: m = Manifold(2,'R^2')
        sage: c_xy = m.chart('x y')
        sage: rot = AutomorphismField(m, 'R') ; rot
        field of tangent-space automorphisms 'R' on the 2-dimensional manifold 'R^2'
        sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
        
    The inverse automorphism is obtained via the method :meth:`inverse`::
    
        sage: inv = rot.inverse() ; inv
        field of tangent-space automorphisms 'inv-R' on the 2-dimensional manifold 'R^2'
        sage: latex(inv)
        R^{-1}
        sage: inv[:]
        [1/2*sqrt(3)         1/2]
        [       -1/2 1/2*sqrt(3)]
        sage: rot[:]
        [1/2*sqrt(3)        -1/2]
        [        1/2 1/2*sqrt(3)]
        sage: inv[:] * rot[:]  # check
        [1 0]
        [0 1]

    """
    def __init__(self, domain, name=None, latex_name=None):
        EndomorphismField.__init__(self, domain, name, latex_name)
        self._inverse = None    # inverse automorphism not set yet

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "field of tangent-space automorphisms"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.domain)
        return description
        
    def _del_derived(self):
        r"""
        Delete the derived quantities
        
        """
        # First delete the derived quantities pertaining to the mother class:
        EndomorphismField._del_derived(self)
        # Then deletes the inverse automorphism:
        self._inverse = None
        
    def _new_instance(self):
        r"""
        Create a :class:`AutomorphismField` instance. 
        
        """
        return AutomorphismField(self.domain)

    def inverse(self):
        r"""
        Return the inverse automorphism.
        
        """
        from sage.matrix.constructor import matrix
        from component import Components
        from vectorframe import CoordFrame
        from utilities import simplify_chain
        if self._inverse is None:
            if self.name is None:
                inv_name = None
            else:
                inv_name = 'inv-' + self.name
            if self.latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self.latex_name + r'^{-1}'
            manif = self.manifold
            dom = self.domain
            si = manif.sindex ; nsi = manif.dim + si
            self._inverse = AutomorphismField(dom, inv_name, inv_latex_name)
            for frame in self.components:
                if isinstance(frame, CoordFrame):
                    chart = frame.chart
                else:
                    chart = dom.def_chart  #!# to be improved
                try:    
                    mat_self = matrix(
                              [[self.comp(frame)[i, j, chart].express 
                              for j in range(si, nsi)] for i in range(si, nsi)])
                except (KeyError, ValueError):
                    continue
                mat_inv = mat_self.inverse()
                cinv = Components(frame, 2)
                for i in range(si, nsi):
                    for j in range(si, nsi):   
                        cinv[i, j, chart] = simplify_chain(mat_inv[i-si,j-si])
                self._inverse.components[frame] = cinv
        return self._inverse


#******************************************************************************

class IdentityMap(AutomorphismField):
    r"""
    Identity map in the tangent spaces of a differentiable manifold.
    
    INPUT:
    
    - ``domain`` -- the manifold domain on which the identity map is 
      defined (must be an instance of class :class:`Domain`)
    - ``name`` -- (default: None) name given to the identity map; if none
      is provided, the value 'Id' is set. 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the identity
      map; if none is provided, the LaTeX symbol is set to `\mathrm{Id}`

    EXAMPLES:

    Identity map on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = m.chart('x y z')
        sage: a = IdentityMap(m) ; a
        Identity map 'Id' in the tangent spaces of the 3-dimensional manifold 'M'
        sage: latex(a)
        \mathrm{Id}
        sage: a[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: a.comp()
        Kronecker delta of size 3x3
        
    The components are automatically defined in any frame::
    
        sage: e = VectorFrame(m, 'e')
        sage: a.comp(e) 
        Kronecker delta of size 3x3
        sage: a.comp(e)[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The components can be read, but cannot be set::
    
        sage: a[1,1]
        1
        sage: a[1,1] = 2
        Traceback (most recent call last):
        ...
        NotImplementedError: The components of the identity map cannot be changed.

    The identity map applied to a vector field::
    
        sage: v = VectorField(m)
        sage: v[:] = (2*x, -3, y+z)
        sage: w = a(v) ; w
        vector field on the 3-dimensional manifold 'M'
        sage: w[:]
        [2*x, -3, y + z]
        sage: w is v  # the output is actually the vector v itself
        True

    The identity map acting as a type (1,1) tensor on a pair (1-form, vector)::
    
        sage: om = OneForm(m)
        sage: om[:] = (0, x*y, 2)
        sage: s = a(om, v) ; s
        scalar field on the 3-dimensional manifold 'M'
        sage: s == om(v)
        True
        
    The identity map is its own inverse::
    
        sage: a.inverse() == a
        True
        sage: a.inverse() is a
        True
        
    """
    def __init__(self, domain, name=None, latex_name=None):
        AutomorphismField.__init__(self, domain)
        if name is None:
            self.name = 'Id'
        else:    
            self.name = name
        if latex_name is None:
            self.latex_name = r'\mathrm{Id}'
        else:
            self.latex_name = latex_name
        self._inverse = self    # the identity is its own inverse
        self.comp() # Initializing the components in the domain's default frame

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "Identity map"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " in the tangent spaces of the " + str(self.domain)
        return description
        
    def _del_derived(self):
        r"""
        Delete the derived quantities
        
        """
        EndomorphismField._del_derived(self)
        
    def _new_comp(self, frame): 
        r"""
        Create some components in the given frame. 
                
        """
        from component import KroneckerDelta
        return KroneckerDelta(frame)

    def comp(self, frame=None):
        r"""
        Return the components in a given frame, as a Kronecker delta.
        
        INPUT:
        
        - ``frame`` -- (default: None) vector frame in which the components are
          required; if none is provided, the components are assumed to refer to
          the manifold's default frame
 
        OUTPUT: 
        
        - components in the frame ``frame``, as an instance of the 
          class :class:`KroneckerDelta` 
        
        EXAMPLES:

        Components of the identity map on a 3-dimensional manifold::
    
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = m.chart('x y z')
            sage: a = IdentityMap(m)
            sage: a.comp()
            Kronecker delta of size 3x3
        
        The components are automatically defined in any frame::
    
            sage: e = VectorFrame(m, 'e')
            sage: a.comp(e) 
            Kronecker delta of size 3x3
            sage: a.comp(e)[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        from component import KroneckerDelta
        if frame is None: 
            frame = self.domain.def_frame
        if frame not in self.components:
            self.components[frame] = KroneckerDelta(frame)
        return self.components[frame]

    def set_comp(self, frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.set_comp`: 
        should not be called. 
        """
        raise NotImplementedError("The components of the identity map " + 
                                  "cannot be changed.")

    def add_comp(self, frame=None):
        r"""
        Redefinition of the TensorField method :meth:`TensorField.add_comp`: 
        should not be called. 
        """
        raise NotImplementedError("The components of the identity map " + 
                                  "cannot be changed.")

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`EndomorphismField.__call__`.
        """
        from diffform import OneForm
        if len(arg) == 1:
            # the identity map acting as such, on a vector, returning a vector:
            vector = arg[0]
            if not isinstance(vector, VectorField):
                raise TypeError("The argument must be a vector field.")
            return vector
            #!# should it be return vector.copy() instead ? 
        elif len(arg) == 2:
            # the identity map acting as a type (1,1) tensor on a pair 
            # (1-form, vector), returning a scalar:
            oneform = arg[0]
            if not isinstance(oneform, OneForm):
                raise TypeError("The first argument must be a 1-form.")
            vector = arg[1]
            if not isinstance(vector, VectorField):
                raise TypeError("The second argument must be a vector field.")
            return oneform(vector)
        else:
            raise TypeError("Bad number of arguments.")
