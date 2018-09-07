#include <cmath.h>
#include <complex>

#ifndef WAVE_H
#define WAVE_H

#define TWOPI 2. * M_PI
#define FOURPI 2. * TWOPI
#define J1 std::complex<double>(0., 1.)

/**
 * Incident wave
 * 
 * @param kvec incident wave vector
 * @param point target point
 * @return amplitude
 */
std::complex<double> 
incident(const double kvec[], const double point[]);

/**
 * Normal gradient of the incident wave, assumes incident wave is exp(1j * kvec.x)
 *
 * @param nvec normal vector pointing inwards
 * @param kvec incident wave vector
 * @param point (source) point
 * @return amplitude
 */
std::complex<double> 
gradIncident(const double nvec[], const double kvec[], 
	         const double point[]);

/**
 * Scattered wave contribution from a single segment
 *
 * @param kvec incident wave vector
 * @param p0 starting point of the segment
 * @param p1 end point of the segment
 * @param point observer point
 * @return wave contribution
 */
std::complex<double> 
computeScatteredWaveElement(const double kvec[], const double p0[], 
	                        const double p1[], const double point[]);

/**
 * Total scattered wave response, summing up contributions from each segment
 *
 * @param kvec incident wave vector
 * @param nc number of contour points
 * @param xc list of x coordinates representing the contour, must close
 * @param yc list of y coordinates representing the contour, must close
 * @param point observer point
 * @return wave response
*/
std::complex<double> 
computeScatteredWave(const double kvec[], int nc, const double xc[], const double yc[], 
	                 const double point[]);

#endif // WAVE_H
