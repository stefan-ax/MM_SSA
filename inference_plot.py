from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def p(x, y, z, V, miu):
	X = np.array([x, y, z])
	f = (X-miu).T
	f2 = X-miu
	V_inv = np.linalg.inv(V)
	return (1/ (np.pi**(3/2) * np.sqrt(np.linalg.det(V))) * np.exp( -1/2 * f.dot(V_inv).dot(f2)  ))

true_theta = [1, 0.1, 0.05]
miu = [0.985, 0.15, 0.04]

x = np.linspace(0.485, 1.485,num = 70)
y = np.linspace(0.05, 0.25,num = 70)
z = np.linspace(-0.06, 0.14,num = 70)

V = np.array([[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

T = []
colors = []

for xx in x:
	for yy in y:
		for zz in z:
			val = p(xx, yy, zz, V, miu)
			if(val > 0.1793):
				T.append([xx, yy, zz])
				colors.append('darkgreen')
			elif(val > 0.1791):
				T.append([xx, yy, zz])
				colors.append('olivedrab')
#			else:
#				T.append([xx, yy, zz])
#				colors.append('lightseagreen')
T = np.array(T)
ax.scatter(true_theta[0], true_theta[1], true_theta[2], c = 'red', s = 500, marker = '*', zorder = 10, label = 'True parameters')
ax.scatter(T[:, 0], T[:, 1], T[:, 2], c = colors, zorder = 0)


ax.set_xlabel('k1')
ax.set_ylabel('k2')
ax.set_zlabel('k3')
plt.title('American fotball ball')

plt.show()

