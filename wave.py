import numpy
from scipy.special import hankel1

ZHAT = numpy.array([0., 0., 1.])
PI = numpy.pi
TWOPI = 2. * PI
FOURPI = 2. * TWOPI


def incident(kvec, points):
	"""
	Incident wave

	@param kvec incident wave vector
	@param point target points, of size n * 3
	@return complex number
	"""
	n = points.shape[0]
	return numpy.exp([1j*kvec.dot(points[i, :]) for i in range(n)])


def gradIncident(kvec, nDotK, points):
	"""
	Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)

	@param kvec incident wave vector
	@param nDotK nvec . kvec array of size n
	@param points (source) point, array of size n * 3
	@return complex number
	"""
	return 1j*nDotK*incident(kvec, points)


def computeScatteredWave(kvec, xc, yc, point):
	"""
	Total scattered wave response, summing up 
	contributions from each segment

	@param kvec incident wave vector
	@param xc list of x coordinates representing the contour, must close
	@param yc list of y coordinates representing the contour, must close
	@param point observer point
	@return complex value
	"""

	kmod = numpy.sqrt(kvec.dot(kvec))

	n = len(xc)
	nm1 = n - 1
	pc = numpy.array([(xc[i], yc[i], 0.) for i in range(n)])
	xdot = pc[1:, :] - pc[:-1, :]
	pmid = 0.5*(pc[1:, :] + pc[:-1, :])
	dsdt = numpy.sqrt(xdot[:, 0]*xdot[:, 0] + xdot[:, 1]*xdot[:, 1])
	nvec = numpy.array([numpy.cross(ZHAT, xdot[i, :]) for i in range(nm1)])
	rvec = numpy.array([point - pmid[i, :] for i in range(nm1)])
	r = numpy.sqrt(rvec[:, 0]*rvec[:, 0] + rvec[:, 1]*rvec[:, 1])
	kr = kmod * r
	nDotR = numpy.array([nvec[i, :].dot(rvec[i, :]) for i in range(nm1)])
	nDotK = numpy.array([nvec[i, :].dot(kvec) for i in range(nm1)])

	g = (1j/4.) * hankel1(0, kr)
	dgdn = (-1j/4.) * hankel1(1, kr) * kmod * nDotR / r

	# contribution from the gradient of the incident wave on the surface
	# of the obstacle. The normal derivative of the scattered wave is 
	# - normal derivative of the incident wave.
	scattered_wave = - dsdt * g * gradIncident(kvec, nDotK, pmid)

	# shadow side: total wave is nearly zero 
	#              => scattered wave amplitude = -incident wave ampl.
	#
	# illuminated side:
	#              => scattered wave amplitude = +incident wave ampl.
	shadow = (nDotK > 0.) - 0.5
	scattered_wave += shadow * dsdt * dgdn * incident(kvec, pmid)

	return scattered_wave.sum()

