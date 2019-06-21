'''
The following function return the mean and the variance of the product of
2 different Gaussian multivariate distributions

Note : It does not include the normalizing Constant. Please check the doc above
for it.

Input : dict1, dict2 - dictionaries with {'mean':mean, 'cov':cov}
Output : The tuple (mean, covariance). 

__author__ = 'Stefan Alexandru Obada'
__date__ = '09/06/2019'

Please check here for the full mathematical path:
    http://www.tina-vision.net/docs/memos/2003-003.pdf
    Thanks to P.A. Bromiley from Uni. of Manchester
'''
import numpy as np

def Prod_Multivar_Normal(dict1, dict2):
    mean1, mean2 = dict1['mean'], dict2['mean']
    cov1, cov2 = dict1['cov'], dict2['cov']
    cov = np.linalg.inv( np.linalg.inv(cov1) + np.linalg.inv(cov2) )
    mean = cov.dot( np.linalg.inv(cov1).dot(mean1) + np.linalg.inv(cov2).dot(mean2) )
    return (mean, cov)