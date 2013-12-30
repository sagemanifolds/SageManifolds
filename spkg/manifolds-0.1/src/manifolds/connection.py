r"""
Affine connections

The class :class:`AffConnection` implements affine connections on 
differentiable manifolds over `\RR`. 

A subclass of :class:`AffConnection` is :class:`LeviCivitaConnection` for
Levi-Civita connections associated to pseudo-Riemannian metrics. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

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
from manifold import Manifold
from component import Components, CompWithSym, CompFullySym
from tensorfield import TensorField, tensor_field_from_comp
from diffform import DiffForm, OneForm

class AffConnection(SageObject):
    r"""
    Base class for affine connections on a differentiable manifold.

    INPUT:
    
    - ``manifold`` -- the manifold on which the connection is defined
    - ``name`` -- name given to the affine connection
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the affine 
      connection

    EXAMPLES:
    
    Affine connection on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = Chart(m, 'x y z', 'xyz')
        sage: nab = AffConnection(m, 'nabla', r'\nabla') ; nab
        affine connection 'nabla' on the 3-dimensional manifold 'M'
        
    A just-created connection has no connection coefficients::
    
        sage: nab.coefficients
        {}

    The connection coefficients relative to the manifold's default frame 
    [here `(\partial/\partial x, \partial/\partial y, \partial/\partial z)`],
    are created by providing the relevant indices inside square brackets::
    
        sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
        sage: nab.coefficients
        {'xyz_b': 3-indices components w.r.t. the coordinate basis 'xyz_b' (d/dx,d/dy,d/dz)}
        
    Unset components are initialized to zero::
    
        sage: nab[:] # list of coefficients relative to the manifold's default vector frame
        [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]
      
    The treatment of connection coefficients in vector frames different 
    from the manifold's default one is similar to that of tensor components; 
    see therefore the class :class:`TensorField` for the documentation. 
    
    Action on a scalar field::
    
        sage: f = ScalarField(m, x^2 - y^2, name='f')
        sage: Df = nab(f) ; Df
        1-form 'df' on the 3-dimensional manifold 'M'
        sage: Df[:]
        [2*x, -2*y, 0]
        sage: Df == f.gradient()  # a necessary condition for any affine connection
        True

    A generic affine connection has some torsion::
    
        sage: DDf = nab(Df) ; DDf
        tensor field 'nabla df' of type (0,2) on the 3-dimensional manifold 'M'
        sage: DDf.antisymmetrize()[:] # nabla does not commute on scalar fields:
        [   0 -x^3    0]
        [ x^3    0    0]
        [   0    0    0]
        
    Let us check the standard formula
    
    .. MATH::
            
        \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f = T^k_{\ \, ij} \nabla_k \, f , 
        
    where the `T^k_{\ \, ij}`'s are the components of the connection's 
    torsion tensor::
            
        sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,Df,0)
        True

    The connection acting on a vector field::

        sage: v = VectorField(m, 'v')
        sage: v[:] = (y*z, x*z, x*y)
        sage: Dv = nab(v) ; Dv
        field of endomorphisms 'nabla v' on the 3-dimensional manifold 'M'
        sage: Dv[:]
        [            0 (x^2*y + 1)*z             y]
        [            z             0             x]
        [            y             x       x*y*z^2]
        
    """
    def __init__(self, manifold, name, latex_name=None):
        if not isinstance(manifold, Manifold):
            raise TypeError("The first argument must be a manifold.")
        self.manifold = manifold
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.coefficients = {}    # coefficients not set yet
        # Initialization of derived quantities:
        self._init_derived() 

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "affine connection"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " on the " + str(self.manifold)
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
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms = {}
        self._torsion_forms = {}
        self._curvature_forms = {}

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._torsion = None
        self._riemann = None
        self._ricci = None
        self._connection_forms.clear()
        self._torsion_forms.clear()
        self._curvature_forms.clear()

    def _new_coef(self, frame_name): 
        r"""
        Create the connection coefficients w.r.t. the given frame. 
        
        This method, to be called by :meth:`coef`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`Components`.
        
        """
        return Components(self.manifold, 3, frame_name)    

    def coef(self, frame_name=None):
        r"""
        Return the connection coefficients relative to the given frame.
        
        `n` being the manifold's dimension, the connection coefficients 
        relative to the vector frame `(e_i)` are the `n^3` scalar fields 
        `\Gamma^k_{\ \, ij}` defined by 
        
        .. MATH::
            
            \nabla_{e_j} e_i = \Gamma^k_{\ \, ij} e_k
        
        
        If the connection coefficients are not known already, they are computed
        from the above formula. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame relative to which the connection coefficients are 
          required; if none is provided, the manifold's default frame is
          assumed
 
        OUTPUT: 
        
        - connection coefficients relative to the frame ``frame_name``, as an 
          instance of the class :class:`Components` with 3 indices ordered 
          as `(k,i,j)`
        
        EXAMPLES:      
 
        Connection coefficient of an affine connection on a 3-dimensional
        manifold::
    
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: nab.coef()  
            3-indices components w.r.t. the coordinate basis 'xyz_b' (d/dx,d/dy,d/dz)
            sage: print type(nab.coef())
            <class 'sage.geometry.manifolds.component.Components'>
            sage: m.default_frame()
            coordinate basis 'xyz_b' (d/dx,d/dy,d/dz)
            sage: nab.coef() is nab.coef('xyz_b')
            True
            sage: nab.coef()[:]  # full list of coefficients:
            [[[0, x^2, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, y*z], [0, 0, 0]]]

        """
        if frame_name is None: 
            frame_name = self.manifold.def_frame.name
        if frame_name not in self.coefficients:
            # the coefficients must be computed
            manif = self.manifold
            ev = manif.frames[frame_name]  # the vector frame
            ef = ev.coframe # the dual frame
            gam = self._new_coef(frame_name)
            for k in manif.irange():
                for i in manif.irange():
                    for j in manif.irange():
                        gam[[k,i,j]] = self(ev(i))(ef(k),ev(j))
            self.coefficients[frame_name] = gam
        return self.coefficients[frame_name]
        

    def set_coef(self, frame_name=None, delete_others=True):
        r"""
        Return the connection coefficients in a given frame for assignment.
        
        See method :meth:`coef` for details about the definition of the 
        connection coefficents. 
        
        NB: the default value of ``delete_others`` ensures that the 
        coefficients with respect to other frames are deleted, in order to 
        avoid any inconsistency. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the connection coefficients are defined; if 
          none is provided, the manifold's default frame is assumed.
        - ``delete_others`` -- (default: True) determines whether the 
          connection coefficients relative to vector frames different from that 
          specified by ``frame_name`` are deleted or not. 
          
        .. WARNING::
        
            Setting ``delete_others`` to False is at the responsability of the 
            user, who must make sure that the various sets of connection 
            coefficients are consistent between each other. 
         
        OUTPUT: 
        
        - connection coefficients in the given frame, as an instance of the 
          class :class:`Components`; if such connection coefficients did not 
          exist previously, they are created.  See method :meth:`coef` for the
          storage convention of the connection coefficents. 
        

        """
        if frame_name is None: frame_name = self.manifold.def_frame.name
        if frame_name not in self.coefficients:
            if frame_name not in self.manifold.frames:
                raise ValueError("The vector frame " + frame_name +
                                 " has not been defined on the " + 
                                 str(self.manifold))
            self.coefficients[frame_name] = self._new_coef(frame_name)
        self._del_derived() # deletes the derived quantities
        if delete_others:
            self.del_other_coef(frame_name)
        return self.coefficients[frame_name]


    def del_other_coef(self, frame_name=None):
        r"""
        Delete all the coefficients but those corresponding to ``frame_name``.
        
        """
        if frame_name is None: frame_name = self.manifold.def_frame.name
        if frame_name not in self.coefficients:
            raise ValueError("The coefficients w.r.t. the vector frame " + 
                             frame_name + " have not been defined.")
        to_be_deleted = []
        for other_frame_name in self.coefficients:
            if other_frame_name != frame_name:
                to_be_deleted.append(other_frame_name)
        for other_frame_name in to_be_deleted:
            del self.coefficients[other_frame_name]

    def __getitem__(self, indices):
        r"""
        Return the connection coefficient w.r.t. the manifold default frame 
        corresponding to the given indices. 

        INPUT:
        
        - ``indices`` -- list of indices
    
        """
        return self.coef()[indices]
        
    def __setitem__(self, indices, value):
        r"""
        Set the connection coefficient w.r.t. the manifold default frame 
        corresponding to the given indices.

        INPUT:
        
        - ``indices`` -- list of indices
    
        """        
        self.set_coef()[indices] = value
        
    def common_frame(self, other):
        r"""
        Find a common vector frame for the components of ``self`` and 
        ``other``. 
        
        In case of multiple common frames, the manifold's default frame is 
        privileged. 
        
        INPUT:
        
        - ``other`` -- a tensor field
        
        OUPUT:
        
        - name of the common frame; if a common frame is not found, None is 
          returned. 
        
        """
        from tensorfield import TensorField
        # Compatibility checks:
        if not isinstance(other, TensorField):
            raise TypeError("The argument must be a tensor field.")
        # Does each object have components on the manifold's default frame ? 
        def_frame_name = self.manifold.def_frame.name
        if def_frame_name in self.coefficients and \
           def_frame_name in other.components:
            frame_name = def_frame_name
        else:
            # Search for a common frame
            frame_name = None
            for frame_name0 in self.coefficients:
                if frame_name0 in other.components:
                    frame_name = frame_name0
                    break
        return frame_name

    def __call__(self, tensor):
        r"""
        Action of the connection on a tensor field.
        
        INPUT:
        
        - ``tensor`` -- a tensor field `T`, of type `(k,\ell)`
        
        OUTPUT:
        
        - tensor field `\nabla T`. 
          
        """
        from tensorfield import TensorField, tensor_field_from_comp
        from scalarfield import ScalarField
        from utilities import format_unop_txt, format_unop_latex
        if not isinstance(tensor, TensorField):
            raise TypeError("The argument must be a tensor field.")
        manif = self.manifold
        if tensor.manifold != manif:
            raise TypeError("The tensor field is not defined on the same " + 
                            "manifold as the connection.")
        if isinstance(tensor, ScalarField):
            return tensor.gradient()
        frame_name = self.common_frame(tensor)
        if frame_name is None:
            raise ValueError("No common frame found for the computation.")
        # Component computation in the common frame:
        frame = manif.frames[frame_name]
        tc = tensor.components[frame_name]
        gam = self.coefficients[frame_name]
        if tensor.sym == [] and tensor.antisym == []:
            resc = Components(manif, tensor.rank+1, frame_name)
        else:
            resc = CompWithSym(manif, tensor.rank+1, frame_name, sym=tensor.sym,
                              antisym=tensor.antisym)
        n_con = tensor.tensor_type[0]
        n_cov = tensor.tensor_type[1]
        for ind in resc.non_redundant_index_generator():
            p = ind[-1]  # derivation index
            ind0 = ind[:-1]
            rsum = frame(p)(tc[[ind0]])
            # loop on contravariant indices:
            for k in range(n_con): 
                for i in manif.irange():
                    indk = list(ind0)
                    indk[k] = i  
                    rsum += gam[[ind0[k], i, p]] * tc[[indk]]
            # loop on covariant indices:
            for k in range(n_con, tensor.rank): 
                for i in manif.irange():
                    indk = list(ind0)
                    indk[k] = i  
                    rsum -= gam[[i, ind0[k], p]] * tc[[indk]]
            resc[[ind]] = rsum
        # Resulting tensor
        result = tensor_field_from_comp(resc, (n_con, n_cov+1))
        result.name = format_unop_txt(self.name + ' ', tensor.name)
        result.latex_name = format_unop_latex(self.latex_name + ' ', 
                                              tensor.latex_name)
        return result

    def torsion(self, frame_name=None):
        r""" 
        Return the connection's torsion tensor.
        
        The torsion tensor is the tensor field `T` of type (1,2) defined by

        .. MATH::
            
            T(\omega, u, v) = \left\langle \omega, \nabla_u v - \nabla_v u
                - [u, v] \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u` and `v`. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          connection coefficients are known, privileging the manifold's default 
          frame.
          
        OUTPUT:
        
        - the torsion tensor `T`, as an instance of :class:`TensorField`
        
        EXAMPLES:
        
        Torsion of an affine connection on a 3-dimensional manifold::
    
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla')
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz
            sage: t = nab.torsion() ; t
            tensor field of type (1,2) on the 3-dimensional manifold 'M'
            sage: t.symmetries()
            no symmetry;  antisymmetry: (1, 2)
            sage: t[:]
            [[[0, -x^2, 0], [x^2, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, -y*z], [0, y*z, 0]]]

        The torsion expresses the lack of commutativity of two successive 
        derivatives of a scalar field::
        
            sage: f = ScalarField(m, x*z^2 + y^2 - z^2, name='f')
            sage: DDf = nab(nab(f)) ; DDf
            tensor field 'nabla df' of type (0,2) on the 3-dimensional manifold 'M'
            sage: DDf.antisymmetrize()[:]  # two successive derivatives do not commute:
            [             0   -1/2*x^2*z^2              0]
            [   1/2*x^2*z^2              0 -(x - 1)*y*z^2]
            [             0  (x - 1)*y*z^2              0]
            sage: 2*DDf.antisymmetrize() == nab.torsion().contract(0,nab(f),0)
            True

        The above identity is the standard formula
         
        .. MATH::
            
            \nabla_j \nabla_i \, f - \nabla_i \nabla_j \, f = T^k_{\ \, ij} \nabla_k \, f , 
        
        where the `T^k_{\ \, ij}`'s are the components of the torsion tensor. 
  
        """
        if self._torsion is None:
            manif = self.manifold
            if frame_name is None:
                if manif.def_frame.name in self.coefficients:
                    frame_name = manif.def_frame.name
                else: # a random frame is picked
                    frame_name = self.coefficients.items()[0][0]
            frame = manif.frames[frame_name]
            gam = self.coef(frame_name)
            sc = frame.structure_coef()
            self._torsion = TensorField(manif, 1, 2, antisym=(1,2))   
            res = self._torsion.set_comp(frame_name)
            for k in manif.irange():
                for i in manif.irange():
                     for j in manif.irange(start=i+1):
                         res[[k,i,j]] = gam[[k,j,i]] - gam[[k,i,j]] - \
                                        sc[[k,i,j]]
        return self._torsion 

    def riemann(self, frame_name=None):
        r""" 
        Return the connection's Riemann curvature tensor.

        The Riemann curvature tensor is the tensor field `R` of type (1,3) 
        defined by

        .. MATH::
            
            R(\omega, u, v, w) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`. 

        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          connection coefficients are known, privileging the manifold's default 
          frame.
          
        OUTPUT:
        
        - the Riemann curvature tensor `R`, as an instance of 
          :class:`TensorField`
        
        EXAMPLES:
        
        Curvature of an affine connection on a 3-dimensional manifold::
            
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla') ; nab
            affine connection 'nabla' on the 3-dimensional manifold 'M'
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
            sage: r = nab.riemann() ; r
            tensor field of type (1,3) on the 3-dimensional manifold 'M'
            sage: r.symmetries()
            no symmetry;  antisymmetry: (2, 3)
            sage: r[:]
            [[[[0, 2*x, 0], [-2*x, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, z], [0, -z, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]]]
       
        """
        if self._riemann is None:
            manif = self.manifold
            if frame_name is None:
                if manif.def_frame.name in self.coefficients:
                    frame_name = manif.def_frame.name
                else: # a random frame is picked
                    frame_name = self.coefficients.items()[0][0]
            ev = manif.frames[frame_name]
            gam = self.coef(frame_name)
            sc = ev.structure_coef()
            gam_gam = gam.contract(1, gam, 0)
            gam_sc = gam.contract(2, sc, 0)
            self._riemann = TensorField(manif, 1, 3, antisym=(2,3))   
            res = self._riemann.set_comp(frame_name)
            for i in manif.irange():
                for j in manif.irange():
                    for k in manif.irange():
                        # antisymmetry of the Riemann tensor taken into account 
                        # by l>k: 
                        for l in manif.irange(start=k+1):
                            res[i,j,k,l] = ev(k)(gam[[i,j,l]]) - \
                                           ev(l)(gam[[i,j,k]]) + \
                                           gam_gam[[i,k,j,l]] -  \
                                           gam_gam[[i,l,j,k]] -  \
                                           gam_sc[[i,j,k,l]]
        return self._riemann 
        

    def ricci(self, frame_name=None):
        r""" 
        Return the connection's Ricci tensor.
        
        The Ricci tensor is the tensor field `Ric` of type (0,2) 
        defined from the Riemann curvature tensor `R` by 

        .. MATH::
            
            Ric(u, v) = R(e^i, u, e_i, v)
        
        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe. 
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          connection coefficients are known, privileging the manifold's default 
          frame.
          
        OUTPUT:
        
        - the Ricci  tensor `Ric`, as an instance of :class:`TensorField`
        
        EXAMPLES:
        
        Ricci tensor of an affine connection on a 3-dimensional manifold::
            
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla') ; nab
            affine connection 'nabla' on the 3-dimensional manifold 'M'
            sage: nab[1,1,2], nab[3,2,3] = x^2, y*z  # Gamma^1_{12} = x^2, Gamma^3_{23} = yz 
            sage: r = nab.ricci() ; r
            tensor field of type (0,2) on the 3-dimensional manifold 'M'
            sage: r[:]
            [  0 2*x   0]
            [  0  -z   0]
            [  0   0   0]
       
        """
        if self._ricci is None:
            self._ricci = self.riemann(frame_name).self_contract(0,2)
        return self._ricci 
        
    def connection_form(self, i, j, frame_name=None):
        r"""
        Return the connection 1-form corresponding to the given index and
        vector frame.
        
        The connection 1-forms with respect to the frame `(e_i)` are the 
        `n^2` 1-forms `\omega^i_{\ \, j}` defined by 
        
        .. MATH::
        
            \nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle
                \, e_i
                
        for any vector `v`. 
        
        The components of `\omega^i_{\ \, j}` in the coframe `(e^i)` dual to 
        `(e_i)` are nothing but the connection coefficients `\Gamma^i_{\ \, jk}`
        relative to the frame `(e_i)`:

        .. MATH::
        
            \omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k
        
        
        INPUT:
        
        - ``i``, ``j`` -- indices identifying the 1-form `\omega^i_{\ \, j}`
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame relative to which the connection 1-forms are defined; if 
          none is provided, the manifold's default frame is assumed. 
          
        OUTPUT:
        
        - the 1-form `\omega^i_{\ \, j}`, as an instance of :class:`OneForm`
        
        EXAMPLES:
        
        Connection 1-forms on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.connection_form(1,1)  # connection 1-form (i,j)=(1,1) w.r.t. M's default frame
            1-form 'nabla connection 1-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.connection_form(1,1)[:]
            [x*y*z, x^2, -y*z]
            
        Connection 1-forms w.r.t. a non-holonomic frame::
        
            sage: ch_basis = AutomorphismField(m)
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = m.default_frame().new_frame(ch_basis, 'e')
            sage: e(1)[:], e(2)[:], e(3)[:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: nab.connection_form(1,1,'e')
            1-form 'nabla connection 1-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.connection_form(1,1,'e').comp('e')[:]
            [x*y^2*z, (x^2*y + 1)*z/y, -x*y*z]
            
        Check of the formula `\omega^i_{\ \, j} = \Gamma^i_{\ \, jk} e^k`::
        
            sage: #... on the manifold's default frame (d/dx, d/dy, d:dz)
            sage: dx = m.default_frame().coframe ; dx
            coordinate coframe 'xyz_b' (dx,dy,dz)
            sage: check = []
            sage: for i in m.irange():
            ...       for j in m.irange():
            ...           check.append( nab.connection_form(i,j) == sum( nab[i,j,k]*dx(k) for k in m.irange() ) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]
            sage: #... on the frame 'e'
            sage: ef = e.coframe ; ef
            coframe 'e' on the 3-dimensional manifold 'M'
            sage: check = []
            sage: for i in m.irange():
            ...       for j in m.irange():
            ...           s = nab.connection_form(i,j,'e').comp('xyz_b', from_frame='e') 
            ...           check.append( nab.connection_form(i,j,'e') == sum( nab.coef('e')[i,j,k]*ef(k) for k in m.irange() ) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]

        Check of the formula `\nabla_v e_j = \langle \omega^i_{\ \, j}, v \rangle e_i`::
        
            sage: v = VectorField(m)
            sage: v[:] = (x*y, z^2-3*x, z+2*y)
            sage: b = m.default_frame()
            sage: for j in m.irange():  # check on M's default frame
            ...       nab(b(j)).contract(1,v,0) == sum( nab.connection_form(i,j)(v)*b(i) for i in m.irange())
            True
            True
            True
            sage: for j in m.irange():  # check on frame 'e'
            ...       nab(e(j)).contract(1,v,0) == sum( nab.connection_form(i,j,'e')(v)*e(i) for i in m.irange())
            True
            True
            True

        """
        if frame_name is None:
            frame_name = self.manifold.def_frame.name
        if frame_name not in self._connection_forms:
            forms = {}
            for i1 in self.manifold.irange():
                for j1 in self.manifold.irange():
                    name = self.name + " connection 1-form (" + str(i1) + \
                           "," + str(j1) + ")"
                    latex_name = r"\omega^" + str(i1) + r"_{\ \, " + str(j1) + \
                                 "}"
                    omega = OneForm(self.manifold, name, latex_name)
                    comega = omega.set_comp(frame_name)
                    for k in self.manifold.irange():
                        comega[k] = self.coef(frame_name)[[i1,j1,k]]
                    forms[(i1,j1)] = omega
            self._connection_forms[frame_name] = forms
        return  self._connection_forms[frame_name][(i,j)] 
                    
    def torsion_form(self, i, frame_name=None):
        r"""
        Return the torsion 2-form corresponding to the given index and
        vector frame.
        
        The torsion 2-forms with respect to the frame `(e_i)` are the 
        `n` 2-forms `\theta^i` defined by 
        
        .. MATH::
        
            \theta^i(u,v) = T(e^i, u, v)
            
        where `T` is the connection's torsion tensor (cf. :meth:`torsion`),
        `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is a generic pair of 
        vectors.
        
        INPUT:
        
        - ``i`` -- index identifying the 2-form `\theta^i`
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame relative to which the torsion 2-forms are defined; if 
          none is provided, the manifold's default frame is assumed. 
          
        OUTPUT:
        
        - the 2-form `\theta^i`, as an instance of :class:`DiffForm`
        
        EXAMPLES:
        
        Torsion 2-forms on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.torsion_form(1)
            2-form 'nabla torsion 2-form (1)' on the 3-dimensional manifold 'M'
            sage: nab.torsion_form(1)[:]                               
            [               0             -x^2      (y^2 + y)*z]
            [             x^2                0  x^3 - x^2 + y^2]
            [    -(y^2 + y)*z -x^3 + x^2 - y^2                0]
            
        Torsion 2-forms w.r.t. a non-holonomic frame::
            
            sage: ch_basis = AutomorphismField(m)                      
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = m.default_frame().new_frame(ch_basis, 'e')
            sage: e(1)[:], e(2)[:], e(3)[:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: ef = e.coframe 
            sage: ef(1)[:], ef(2)[:], ef(3)[:]
            ([1/y, 0, 0], [0, 1/z, 0], [0, 0, 1/x])
            sage: nab.torsion_form(1, 'e')
            2-form 'nabla torsion 2-form (1)' on the 3-dimensional manifold 'M'
            sage: nab.torsion_form(1, 'e').comp('e')[:]
            [                       0                   -x^2*z          (x*y^2 + x*y)*z]
            [                   x^2*z                        0  (x^4 - x^3 + x*y^2)*z/y]
            [        -(x*y^2 + x*y)*z -(x^4 - x^3 + x*y^2)*z/y                        0]
            
        Cartan's first structure equation is
        
        .. MATH::
        
            \theta^i = \mathrm{d} e^i + \omega^i_{\ \, j} \wedge e^j
            
        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf. 
        :meth:`connection_form`). Let us check it on the frame 'e'::
        
            sage: for i in m.irange():
            ...       nab.torsion_form(i, 'e') == xder(ef(i)) + sum(nab.connection_form(i,k,'e').wedge(ef(k)) for k in m.irange())
            ...
            True
            True
            True

        """
        if frame_name is None:
            frame_name = self.manifold.def_frame.name
        if frame_name not in self._torsion_forms:
            forms = {}
            for i1 in self.manifold.irange():
                name = self.name + " torsion 2-form (" + str(i1) + ")"
                latex_name = r"\theta^" + str(i1)
                theta = DiffForm(self.manifold, 2, name, latex_name)
                ctheta = theta.set_comp(frame_name)
                for k in self.manifold.irange():
                    for l in self.manifold.irange(start=k+1):
                        ctheta[k,l] = \
                            self.torsion(frame_name).comp(frame_name)[[i1,k,l]]
                forms[i1] = theta
            self._torsion_forms[frame_name] = forms
        return  self._torsion_forms[frame_name][i] 

                    
    def curvature_form(self, i, j, frame_name=None):
        r"""
        Return the curvature 2-form corresponding to the given index and
        vector frame.
        
        The curvature 2-forms with respect to the frame `(e_i)` are the 
        `n^2` 2-forms `\Omega^i_{\ \, j}` defined by 
        
        .. MATH::
        
            \Omega^i_{\ \, j}(u,v) = R(e^i, u, v, e_j)
            
        where `R` is the connection's Riemann curvature tensor (cf. 
        :meth:`riemann`), `(e^i)` is the coframe dual to `(e_i)` and `(u,v)` is a 
        generic pair of vectors.
        
        INPUT:
        
        - ``i``, ``j`` -- indices identifying the 2-form `\Omega^i_{\ \, j}`
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame relative to which the curvature 2-forms are defined; if 
          none is provided, the manifold's default frame is assumed. 
          
        OUTPUT:
        
        - the 2-form `\Omega^i_{\ \, j}`, as an instance of :class:`DiffForm`
        
        EXAMPLES:
        
        Curvature 2-forms on a 3-dimensional manifold::
        
            sage: m = Manifold(3, 'M', start_index=1)
            sage: c_xyz = Chart(m, 'x y z', 'xyz')
            sage: nab = AffConnection(m, 'nabla', r'\nabla')
            sage: nab[1,1,1], nab[1,1,2], nab[1,1,3] = x*y*z, x^2, -y*z
            sage: nab[1,2,3], nab[1,3,1], nab[1,3,2] = -x^3, y^2*z, y^2-x^2
            sage: nab[2,1,1], nab[2,1,2], nab[2,2,1] = z^2, x*y*z^2, -x^2
            sage: nab[2,3,1], nab[2,3,3], nab[3,1,2] = x^2+y^2+z^2, y^2-z^2, x*y+z^2
            sage: nab[3,2,1], nab[3,2,2], nab[3,3,3] = x*y+z, z^3 -y^2, x*z^2 - z*y^2
            sage: nab.curvature_form(1,1)
            2-form 'nabla curvature 2-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.curvature_form(1,1)[:]
            [                             0  y^2*z^3 + (x*y^3 - x)*z + 2*x                  x^3*z^2 - x*y]
            [-y^2*z^3 - (x*y^3 - x)*z - 2*x                              0                  x^4*y*z^2 - z]
            [                -x^3*z^2 + x*y                 -x^4*y*z^2 + z                              0]
            
        Curvature 2-forms w.r.t. a non-holonomic frame::
            
            sage: ch_basis = AutomorphismField(m)                      
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = y, z, x
            sage: e = m.default_frame().new_frame(ch_basis, 'e')
            sage: e(1)[:], e(2)[:], e(3)[:]
            ([y, 0, 0], [0, z, 0], [0, 0, x])
            sage: ef = e.coframe 
            sage: ef(1)[:], ef(2)[:], ef(3)[:]
            ([1/y, 0, 0], [0, 1/z, 0], [0, 0, 1/x])
            sage: nab.curvature_form(1,1,'e')
            2-form 'nabla curvature 2-form (1,1)' on the 3-dimensional manifold 'M'
            sage: nab.curvature_form(1,1,'e').comp('e')[:]
            [                                     0  y^3*z^4 + (x*y^4 - x*y)*z^2 + 2*x*y*z                    x^4*y*z^2 - x^2*y^2]
            [-y^3*z^4 - (x*y^4 - x*y)*z^2 - 2*x*y*z                                      0                      x^5*y*z^3 - x*z^2]
            [                  -x^4*y*z^2 + x^2*y^2                     -x^5*y*z^3 + x*z^2                                      0]
            
        Cartan's second structure equation is
         
        .. MATH::
        
            \Omega^i_{\ \, j} = \mathrm{d} \omega^i_{\ \, j} + \omega^i_{\ \, k} \wedge \omega^k_{\ \, j}
            
        where the `\omega^i_{\ \, j}`'s are the connection 1-forms (cf. 
        :meth:`connection_form`). Let us check it on the frame 'e'::
       
            sage: omega = nab.connection_form
            sage: check = []
            sage: for i in m.irange():
            ...       for j in m.irange():
            ...           check.append( nab.curvature_form(i,j,'e') == xder(omega(i,j,'e')) + \
            ...           sum( omega(i,k,'e').wedge(omega(k,j,'e')) for k in m.irange()) )
            ...
            sage: check
            [True, True, True, True, True, True, True, True, True]
            
        """
        if frame_name is None:
            frame_name = self.manifold.def_frame.name
        if frame_name not in self._curvature_forms:
            forms = {}
            for i1 in self.manifold.irange():
                for j1 in self.manifold.irange():
                    name = self.name + " curvature 2-form (" + str(i1) + \
                           "," + str(j1) + ")"
                    latex_name = r"\Omega^" + str(i1) + r"_{\ \, " + str(j1) + \
                                 "}"
                    omega = DiffForm(self.manifold, 2, name, latex_name)
                    comega = omega.set_comp(frame_name)
                    for k in self.manifold.irange():
                        for l in self.manifold.irange(start=k+1):
                            comega[k,l] = \
                          self.riemann(frame_name).comp(frame_name)[[i1,j1,k,l]]
                    forms[(i1,j1)] = omega
            self._curvature_forms[frame_name] = forms
        return  self._curvature_forms[frame_name][(i,j)] 

            
#******************************************************************************
            
class LeviCivitaConnection(AffConnection):
    r"""
    Levi-Civita connection on a differentiable manifold.

    INPUT:
    
    - ``metric`` -- the metric defining the Levi-Civita connection, as an
      instance of :class:`Metric`
    - ``name`` -- name given to the connection
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the  
      connection

    EXAMPLES:
    
    Levi-Civita connection associated with the Euclidean metric on `\RR^3`
    expressed in spherical coordinates::
    
        sage: m = Manifold(3, 'R^3', start_index=1)
        sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
        sage: g = Metric(m, 'g')
        sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2
        sage: g.show()
        g = dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
        sage: nab = LeviCivitaConnection(g, 'nabla', r'\nabla') ; nab
        Levi-Civita connection 'nabla' associated with the pseudo-Riemannian metric 'g' on the 3-dimensional manifold 'R^3'
    
    Let us check that the connection is compatible with the metric::
    
        sage: Dg = nab(g) ; Dg
        tensor field 'nabla g' of type (0,3) on the 3-dimensional manifold 'R^3'
        sage: Dg == 0
        True

    and that it is torsionless::
    
        sage: nab.torsion() == 0
        True
        sage: AffConnection.torsion(nab) == 0  # forces the computation of the torsion
        True
 
    The connection coefficients in the manifold's default frame are Christoffel 
    symbols, since the default frame is a coordinate basis::
    
        sage: m.default_frame()
        coordinate basis 'spher_b' (d/dr,d/dth,d/dph)
        sage: nab.coef()  
        3-indices components w.r.t. the coordinate basis 'spher_b' (d/dr,d/dth,d/dph), with symmetry on the index positions (1, 2)
        sage: # note that the Christoffel symbols are symmetric with respect to their last two indices (positions (1,2))
        sage: nab.coef()[:]
        [[[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(th)^2]],
        [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -sin(th)*cos(th)]],
        [[0, 0, 1/r], [0, 0, cos(th)/sin(th)], [1/r, cos(th)/sin(th), 0]]]


    """
    def __init__(self, metric, name, latex_name=None):
        AffConnection.__init__(self, metric.manifold, name, latex_name)
        self.metric = metric
        # initialization of the Christoffel symbols in the manifold's default chart:
        self.coef(self.manifold.def_chart.frame.name) 
        
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "Levi-Civita connection"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " associated with the " + str(self.metric)
        return description

    def _init_derived(self):
        r"""
        Initialise the derived quantities
        """
        AffConnection._init_derived(self)
        self._ricci_scalar = None

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        AffConnection._del_derived(self)
        self._ricci_scalar = None

    def coef(self, frame_name=None):
        r"""
        Return the connection coefficients relative to the given frame.
        
        `n` being the manifold's dimension, the connection coefficients 
        relative to the vector frame `(e_i)` are the `n^3` scalar fields 
        `\Gamma^k_{\ \, ij}` defined by 
        
        .. MATH::
            
            \nabla_{e_j} e_i = \Gamma^k_{\ \, ij} e_k
        
        If the connection coefficients are not known already, they are computed

         * as Christoffel symbols if the frame `(e_i)` is a coordinate basis
         * frome the above formula otherwise 
                
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame relative to which the connection coefficients are 
          required; if none is provided, the manifold's default frame is
          assumed
 
        OUTPUT: 
        
        - connection coefficients relative to the frame ``frame_name``, as an 
          instance of the class :class:`Components` with 3 indices ordered 
          as `(k,i,j)`; for Christoffel symbols, an instance of the subclass
          :class:`CompWithSym` is returned. 
        
        EXAMPLES:
        
        Christoffel symbols of the Levi-Civita connection associated to 
        the Euclidean metric on `\RR^3` expressed in spherical coordinates::
        
            sage: m = Manifold(3, 'R^3', start_index=1)
            sage: c_spher = Chart(m, r'r:positive, th:positive:\theta, ph:\phi', 'spher')
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[2,2], g[3,3] = 1, r^2 , (r*sin(th))^2
            sage: g.show()
            g = dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
            sage: nab = LeviCivitaConnection(g, 'nabla', r'\nabla')
            sage: gam = nab.coef() ; gam
            3-indices components w.r.t. the coordinate basis 'spher_b' (d/dr,d/dth,d/dph), with symmetry on the index positions (1, 2)
            sage: gam[:]
            [[[0, 0, 0], [0, -r, 0], [0, 0, -r*sin(th)^2]],
            [[0, 1/r, 0], [1/r, 0, 0], [0, 0, -sin(th)*cos(th)]],
            [[0, 0, 1/r], [0, 0, cos(th)/sin(th)], [1/r, cos(th)/sin(th), 0]]]
            sage: # The only non-zero Christoffel symbols:
            sage: gam[1,2,2], gam[1,3,3]
            (-r, -r*sin(th)^2)
            sage: gam[2,1,2], gam[2,3,3]
            (1/r, -sin(th)*cos(th))
            sage: gam[3,1,3], gam[3,2,3]
            (1/r, cos(th)/sin(th))
            
        Connection coefficients of the same connection with respect to the 
        orthonormal frame associated to spherical coordinates::
        
            sage: ch_basis = AutomorphismField(m) 
            sage: ch_basis[1,1], ch_basis[2,2], ch_basis[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: e = m.frame('spher_b').new_frame(ch_basis, 'e')
            sage: gam_e = nab.coef('e') ; gam_e
            3-indices components w.r.t. the vector frame 'e'
            sage: gam_e[:]
            [[[0, 0, 0], [0, -1/r, 0], [0, 0, -1/r]],
            [[0, 1/r, 0], [0, 0, 0], [0, 0, -cos(th)/(r*sin(th))]],
            [[0, 0, 1/r], [0, 0, cos(th)/(r*sin(th))], [0, 0, 0]]]
            sage: # The only non-zero connection coefficients:
            sage: gam_e[1,2,2], gam_e[2,1,2]
            (-1/r, 1/r)
            sage: gam_e[1,3,3], gam_e[3,1,3]
            (-1/r, 1/r)
            sage: gam_e[2,3,3], gam_e[3,2,3]
            (-cos(th)/(r*sin(th)), cos(th)/(r*sin(th)))

        """
        from vectorframe import CoordBasis
        if frame_name is None: 
            frame_name = self.manifold.def_frame.name
        if frame_name not in self.coefficients:
            # the coefficients must be computed
            manif = self.manifold
            frame = manif.frames[frame_name]
            if isinstance(frame, CoordBasis):
                # Christoffel symbols 
                chart_name = frame_name[:-2]
                gam = CompWithSym(manif, 3, frame_name, sym=(1,2))
                gg = self.metric.comp(frame_name)
                ginv = self.metric.inverse().comp(frame_name)
                for ind in gam.non_redundant_index_generator():
                    i, j, k = ind
                    # The computation is performed at the FunctionChart level:
                    rsum = 0
                    for s in manif.irange():
                        rsum += ginv[i,s, chart_name] * ( 
                                            gg[s,k, chart_name].diff(j)
                                          + gg[j,s, chart_name].diff(k)
                                          - gg[j,k, chart_name].diff(s) )
                    gam[i,j,k, chart_name] = rsum / 2
                    self.coefficients[frame_name] = gam
            else:
                # Computation from the formula defining the connection coef.
                return AffConnection.coef(self, frame_name)
        return self.coefficients[frame_name]

    def torsion(self, frame_name=None):
        r""" 
        Return the connection's torsion tensor (identically zero for a 
        Levi-Civita connection). 
        
        See :meth:`AffConnection.torsion` for the general definition of the 
        torsion tensor. 
        
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the torsion must be initialized; if none is 
          provided, the manifold's default frame is assumed.
          
        OUTPUT:
        
        - the torsion tensor `T`, as a vanishing instance of :class:`TensorField`
          
        """
        if self._torsion is None:
            manif = self.manifold
            if frame_name is None:
                frame_name = manif.def_frame.name
            self._torsion = TensorField(manif, 1, 2, antisym=(1,2))
            # Initialization of the frame_name components to zero: 
            self._torsion.set_comp(frame_name) 
        return self._torsion 

    def riemann(self, frame_name=None, name=None, latex_name=None):
        r""" 
        Return the Riemann curvature tensor associated with the metric.

        This method redefines :meth:`AffConnection.riemann` to set some name
        and the latex_name to the output.
        
        The Riemann curvature tensor is the tensor field `R` of type (1,3) 
        defined by

        .. MATH::
            
            R(\omega, u, v, w) = \left\langle \omega, \nabla_u \nabla_v w
                - \nabla_v \nabla_u w - \nabla_{[u, v]} w \right\rangle
        
        for any 1-form  `\omega`  and any vector fields `u`, `v` and `w`. 

        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          metric components are known, privileging the manifold's default 
          frame.
        - ``name`` -- (default: None) name given to the Riemann tensor; 
          if none, it is set to "Riem(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
          Riemann tensor; if none, it is set to "\\mathrm{Riem}(g)", where "g" 
          is the metric's name

        OUTPUT:
        
        - the Riemann curvature tensor `R`, as an instance of 
          :class:`TensorField`
          
        """
        if self._riemann is None:
            AffConnection.riemann(self, frame_name)
            if name is None:
                self._riemann.name = "Riem(" + self.metric.name + ")"
            else:
                self._riemann.name = name
            if latex_name is None:
                self._riemann.latex_name = r"\mathrm{Riem}\left(" + self.metric.latex_name + r"\right)"
            else:
                self._riemann.latex_name = latex_name
        return self._riemann
            


    def ricci(self, frame_name=None, name=None, latex_name=None):
        r""" 
        Return the connection's Ricci tensor.
        
        This method redefines :meth:`AffConnection.ricci` to take into account
        the symmetry of the Ricci tensor for a Levi-Civita connection. 

        The Ricci tensor is the tensor field `Ric` of type (0,2) 
        defined from the Riemann curvature tensor `R` by 

        .. MATH::
            
            Ric(u, v) = R(e^i, u, e_i, v)
        
        for any vector fields `u` and `v`, `(e_i)` being any vector frame and
        `(e^i)` the dual coframe. 
                
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          connection coefficients are known, privileging the manifold's default 
          frame.
        - ``name`` -- (default: None) name given to the Ricci tensor; 
          if none, it is set to "Ric(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
          Ricci tensor; if none, it is set to "\\mathrm{Ric}(g)", where "g" 
          is the metric's name
          
        OUTPUT:
        
        - the Ricci tensor `Ric`, as an instance of :class:`SymBilinFormField`
        
        EXAMPLES:
        
        Ricci tensor of the standard connection on the 2-dimensional sphere::
        
            sage: m = Manifold(2, 'S^2', start_index=1)
            sage: c_spher = Chart(m, r'th:\theta, ph:\phi', 'spher')
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[2,2] = 1, sin(th)^2
            sage: g.show() # standard metric on S^2:
            g = dth*dth + sin(th)^2 dph*dph
            sage: nab = g.connection() ; nab
            Levi-Civita connection 'nabla_g' associated with the pseudo-Riemannian metric 'g' on the 2-dimensional manifold 'S^2'
            sage: ric = nab.ricci() ; ric             
            field of symmetric bilinear forms 'Ric(g)' on the 2-dimensional manifold 'S^2'
            sage: ric.show()
            Ric(g) = dth*dth + sin(th)^2 dph*dph            
        
        Checking that the Ricci tensor of the Levi-Civita connection associated
        to Schwarzschild metric is identically zero::
        
            sage: m = Manifold(4, 'M')
            sage: c_BL = Chart(m, r't, r:positive, th:positive:\theta, ph:\phi', 'BL') # Boyer-Linquist coordinates
            sage: var('M')
            M
            sage: g = Metric(m, 'g')
            sage: g[0,0], g[1,1] = -(1-2*M/r), 1/(1-2*M/r)
            sage: g[2,2], g[3,3] = r^2, (r*sin(th))^2
            sage: g.show()
            g = (2*M/r - 1) dt*dt - 1/(2*M/r - 1) dr*dr + r^2 dth*dth + r^2*sin(th)^2 dph*dph
            sage: nab = g.connection() ; nab
            Levi-Civita connection 'nabla_g' associated with the pseudo-Riemannian metric 'g' on the 4-dimensional manifold 'M'
            sage: ric = nab.ricci() ; ric
            field of symmetric bilinear forms 'Ric(g)' on the 4-dimensional manifold 'M'
            sage: ric == 0
            True

        """
        if self._ricci is None:
            manif = self.manifold
            riem = self.riemann(frame_name)
            if frame_name is None:
                if manif.def_frame.name in riem.components:
                    frame_name = manif.def_frame.name
                else: # a random frame is picked
                    frame_name = riem.components.items()[0][0]
            criem = riem.components[frame_name]
            cric = CompFullySym(manif, 2, frame_name)
            si = manif.sindex
            for i in manif.irange():
                # symmetry of the Ricci tensor taken into account by j>=i: 
                for j in manif.irange(start=i):  
                    rsum = criem[[si,i,si,j]].copy()
                    for k in manif.irange(start=si+1):
                        rsum += criem[[k,i,k,j]]
                    cric[i,j] = rsum
            self._ricci = tensor_field_from_comp(cric, (0,2)) 
            if name is None:
                self._ricci.name = "Ric(" + self.metric.name + ")"
            else:
                self._ricci.name = name
            if latex_name is None:
                self._ricci.latex_name = r"\mathrm{Ric}\left(" + self.metric.latex_name + r"\right)"
            else:
                self._ricci.latex_name = latex_name
        return self._ricci 

    def ricci_scalar(self, frame_name=None, name=None, latex_name=None):
        r""" 
        Return the connection's Ricci scalar.
        
        The Ricci scalar is the scalar field `r` defined from the Ricci tensor 
        `Ric` and the metric tensor `g` by 

        .. MATH::
            
            r = g^{ij} Ric_{ij}
        
        INPUT:
        
        - ``frame_name`` -- (default: None) string containing the name of the 
          vector frame in which the computation must be performed; if none is 
          provided, the computation is performed in a frame for which the 
          connection coefficients are known, privileging the manifold's default 
          frame.
        - ``name`` -- (default: None) name given to the Ricci scalar; 
          if none, it is set to "r(g)", where "g" is the metric's name
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
          Ricci scalar; if none, it is set to "\\mathrm{r}(g)", where "g" 
          is the metric's name
          
        OUTPUT:
        
        - the Ricci scalar `r`, as an instance of :class:`ScalarField`
        
        EXAMPLES:
        
        Ricci scalar of the standard connection on the 2-dimensional sphere::
        
            sage: m = Manifold(2, 'S^2', start_index=1)
            sage: c_spher = Chart(m, r'th:\theta, ph:\phi', 'spher')
            sage: a = var('a') # the sphere radius
            sage: g = Metric(m, 'g')
            sage: g[1,1], g[2,2] = a^2, a^2*sin(th)^2
            sage: g.show() # standard metric on S^2 with radius a
            g = a^2 dth*dth + a^2*sin(th)^2 dph*dph
            sage: nab = g.connection() ; nab
            Levi-Civita connection 'nabla_g' associated with the pseudo-Riemannian metric 'g' on the 2-dimensional manifold 'S^2'
            sage: r = nab.ricci_scalar() ; r
            scalar field 'r(g)' on the 2-dimensional manifold 'S^2'
            sage: r.expr()
            2/a^2        

        """
        if self._ricci_scalar is None:            
            manif = self.manifold
            ric = self.ricci(frame_name)
            ig = self.metric.inverse()
            frame_name = ig.common_frame(ric)
            cric = ric.components[frame_name]
            cig = ig.components[frame_name]
            si = manif.sindex
            rsum1 = cig[[si,si]] * cric[[si,si]] 
            for i in manif.irange(start=si+1):
                rsum1 += cig[[i,i]] * cric[[i,i]]
            rsum2 = 0
            for i in manif.irange():
                for j in manif.irange(start=i+1):
                    rsum2 += cig[[i,j]] * cric[[i,j]]
            self._ricci_scalar = rsum1 + 2*rsum2
            if name is None:
                self._ricci_scalar.name = "r(" + self.metric.name + ")"
            else:
                self._ricci_scalar.name = name
            if latex_name is None:
                self._ricci_scalar.latex_name = r"\mathrm{r}\left(" + self.metric.latex_name + r"\right)"
            else:
                self._ricci_scalar.latex_name = latex_name
        return self._ricci_scalar 

