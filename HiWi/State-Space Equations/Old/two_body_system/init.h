#ifndef __INIT_H_
#define __INIT_H_

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the 
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 */
int read_parameters_1(
  const char *inFileName1,
  double *stiffCoeff1,
  double *dampCoeff1,
  double *mass1,
  double *tEnd,
  double *tStart,
  double *stepHeight,
  double *dt,
  double *g,
  double *springFreeLength1);

int read_parameters_2(
  const char *inFileName2,
  double *stiffCoeff2,
  double *dampCoeff2,
  double *mass2,
  double *tEnd,
  double *tStart,
  double *stepHeight,
  double *dt,
  double *g,
  double *springFreeLength2);


int write_output( char* outFileName, std::vector<double> time, std::vector<std::vector<double> > disp,
					std::vector<std::vector<double> > vel, std::vector<std::vector<double> > acc, std::vector<std::vector<double> > force,
					std::vector<double> massIniPosition);


#endif

