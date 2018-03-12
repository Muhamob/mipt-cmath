import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d, correlate2d
from scipy.ndimage import correlate1d
import matplotlib.cm as cm
from tqdm import tqdm


class Nodes():
	def __init__(self, depth=1, hx=.01, hy=.01, xmin=0., xmax=1., ymin=0., ymax=1., f=None):
		"""
			depth : number of states to memorize
			f     : zero-distribution
		"""
		self.hx = hx
		self.hy = hy
		self.x = np.arange(xmin, xmax, hx)
		self.y = np.arange(ymin, ymax, hy)
		if f is None:
			self.f = np.zeros((self.x.shape[0], self.y.shape[0]))
		else:
			if f.shape != (self.x.shape[0], self.y.shape[0]):
				raise Exception('shapes aren\'t matching')
			else:
				self.f = f
		self.depth = int(depth)
		if depth != 1:
			self.states = [0]*(int(depth)-1) + [self.f.copy()]

	def correlate(self, kernel):
		"""
			mode:
				-> valid :  uses when there are border conditions
							return (shape[0]-kernel.shape[0]//2, shape[1]-kernel.shape[1]//2)
				-> same  :  uses when bording conditions are neglected
							return same shape array
			kernel: 
		"""
		self.f = correlate2d(self.f, kernel, mode='same')
		if self.depth != 1:
			self.states.pop(0)
			self.states.append(self.f)

	def draw(self, im_cmap=cm.gray, contour_cmap='jet', colorbar=True, contours=True, figsize=(8, 12), interpolate=True, phix=None, phiy=None):
		plt.figure(figsize=figsize)
		interpolation = 'gaussian' if interpolate else 'none'
		im = plt.imshow(self.f, interpolation=interpolation, cmap=im_cmap)
		if colorbar:
			CBI = plt.colorbar(im, orientation='horizontal', shrink=.9)
		if contours:
			xx, yy = np.meshgrid(self.x, self.y)
			cnt = plt.contour(self.f, 12, cmap=contour_cmap)
			plt.clabel(cnt, inline=1, fontsize=8, fmt='%1.1f')
			plt.colorbar(cnt, shrink=.8, extend='both')
		# if phix is not None and phiy is not None:
		# 	xx, yy = np.meshgrid(self.x, self.y)
		# 	plt.quiver(-phix[::3, ::3], -phiy[::3, ::3],
		# 		units='x', pivot='tip', edgecolor='k', alpha=0.5, color='red')
		# plt.imshow(self.f, interpolation=interpolation, cmap=im_cmap, origin='lower',
		#	extent=[self.y.min(), self.y.max(), self.x.min(), self.x.max()])
		plt.imshow(self.f, interpolation=interpolation, cmap=im_cmap, origin='lower')



class Solver():
	def __init__(self, depth=1, hx=.01, hy=.01, xmin=0., xmax=1., ymin=0., ymax=1., ro=None):
		self.phi = Nodes(depth, hx, hy, xmin, xmax, ymin, ymax)
		self.kernel = np.array([[0, 0.5*hx**2/(hx**2+hy**2), 0], 
			[0.5*hy**2/(hx**2+hy**2), 0, 0.5*hy**2/(hx**2+hy**2)],
			[0, 0.5*hx**2/(hx**2+hy**2), 0]])
		if ro is not None:
			self.ro = Nodes(1, hx, hy, xmin, xmax, ymin, ymax, ro)
		else:
			self.ro = None

	def setBorderCond(self, boundary_mask, c=-999999999999):
		"""
			mask	: numpy.ndarray like self.phi.f where mask[i,j]=boundaryCond if
					  there is a border else -999999999999
					  NOTE !!! boundary conditions should not be equal to 0. If it so,
					  then add some constant and finally substract this constant
			c 		: constant to recognize borders
		"""
		if self.phi.f.shape != boundary_mask.shape:
			raise Exception('shapes aren\'t matching')

		self.maskinner = np.where(boundary_mask == c, 1., 0.)
		self.mask = np.where(boundary_mask != c, boundary_mask, 0)
		mean = self.mask.mean()
		self.phi.f = (mean*np.ones_like(self.mask))+self.mask

	def step(self):
		"""
			computes next step of iteration
		"""
		self.phi.correlate(self.kernel)
		if self.ro is None:
			self.phi.f = self.phi.f*self.maskinner + self.mask
		else:
			self.phi.f = (self.phi.f-(1/(1/self.ro.hx**2 + 1/self.ro.hy**2))*self.ro.f)*self.maskinner + self.mask

	def fit(self, iterations=None, e=1e-6):
		"""
			iterations  : number of iterations use if depth >= 1
			e           : precision, use if depth >= 2

		"""
		self.step()
		if iterations is not None:
			for _ in tqdm(range(iterations)):
				self.step()
		else:
			if self.phi.depth == 1:
				raise Exception('depth should be more than 1')
			else:
				err = np.linalg.norm(self.phi.states[-1]-self.phi.states[-2]) / np.linalg.norm(self.phi.states[-1])
				pbar = tqdm(total=50000, desc='e={}'.format(err))
				while err > e:
					self.step()
					pbar.update()
					pbar.set_description('e={:.5%}'.format(err))
					err = np.linalg.norm(self.phi.states[-1]-self.phi.states[-2]) / np.linalg.norm(self.phi.states[-1])
				pbar.close()
		return self.phi

	def diff(self):
		phix = correlate1d(self.phi.f, np.array([-0.5, 0, 0.5])/self.phi.hx, mode='constant', cval=0., axis=1)
		phiy = correlate1d(self.phi.f, np.array([-0.5, 0, 0.5])/self.phi.hy, mode='constant', cval=0., axis=0)
		phix[:, 0] = (self.phi.f[:, 1] - self.phi.f[:, 0])/self.phi.hx
		phix[:, -1] = (self.phi.f[:, -1] - self.phi.f[:, -2])/self.phi.hx
		phiy[0, :] = (self.phi.f[1, :] - self.phi.f[0, :])/self.phi.hy
		phiy[-1, :] = (self.phi.f[-1, :] - self.phi.f[-2, :])/self.phi.hy
		return phix, phiy



if __name__ == '__main__':
	solver = s.Solver(2, .1, .1, 0, 10, 0, 10)
	x = np.arange(0, 10, .1)
	y = np.arange(0, 10, .1)
	c = -999999999999
	bmask = c*np.ones((x.shape[0], y.shape[0]))
	bmask[:,0] = 1*np.ones_like(y)
	bmask[:,-1] = 1*np.ones_like(y)
	bmask[0,:] = 4*np.ones_like(x)
	bmask[-1,:] = 4*np.ones_like(x)
	cind2 = np.array([np.floor(x.shape[0]/2), np.floor(y.shape[0]/2)], dtype=np.uint8)
	bmask[cind2[0]-2:cind2[0]+2, cind2[1]-10:cind2[1]+10] = 10
	solver.setBorderCond(bmask)
	solution = solver.fit(10000)
	solution.draw()