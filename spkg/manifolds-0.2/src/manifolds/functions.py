r"""
SageManifolds functions. 

This module defines functions that are not class methods, but rather shortcuts
to them. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

"""

def xder(form):
    r"""
    Exterior derivative of a differential form
    """
    return form.exterior_der()

def Lie(vector, tensor):
    r"""
    Lie derivative of a tensor field with respect to a vector field
    """
    return tensor.lie_der(vector)
    
def ctr(tensor1, pos1, tensor2, pos2):
    r"""
    Contraction of two tensors
    """
    return tensor1.contract(pos1, tensor2, pos2)

    

    
    
    
    

