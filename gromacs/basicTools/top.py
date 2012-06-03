from math import sqrt, acos, pi

def dot(v1, v2):
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z

def angle(v1, v2, v3=None):
    if v3 is not None:
        v1 = v1-v2
        v2 = v3-v2
    r = dot(v1, v2) / sqrt(dot(v1, v1) * dot(v2, v2))
    if r <= -1.0:
        return pi
    elif r >= 1:
        return 0
    else:
        return acos(r)
        
def torsion(v1, v2, v3, v4):
    axis = v3-v2
    v1 = (v1 - v2) * axis
    v4 = (v3 - v4) * axis
    return angle(v1, v4)