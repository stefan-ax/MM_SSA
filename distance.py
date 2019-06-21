'''
The following function return the manhattan distance

__author__ = 'Stefan Alexandru Obada'
__date__ = '09/06/2019'
'''
from scipy.spatial import distance

def manhattan_distance(arr1, arr2):
    length = min(arr1.shape[0], arr2.shape[0])
    _distance = 0
    for i in range(length):
        #print(arr1[i] - arr2[i])
        _distance = _distance + distance.minkowski(arr1[i], arr2[i], p = 2)#np.linalg.norm(arr1[i] - arr2[i]) #BUILD MANHATTAN !!!
        #print(_distance)
    #print(_distance)
    return _distance