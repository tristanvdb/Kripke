#include <Kripke/Directions.h>
#include <Kripke/Grid.h>
#include <Kripke/Input_Variables.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <algorithm>

namespace {
  /*
    GaussLegendre returns the n point Gauss-Legendre quadrature rule for
    the integral between x1 and x2.
  */
  void GaussLegendre(double x1, double x2, std::vector<double> &x,
      std::vector<double> &w, double eps)
  {
    int n = x.size();
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for(i=1; i<=m; i++){
      z=cos(M_PI*(i-0.25)/(n+0.5));
      do {
        p1=1.0;
        p2=0.0;
        for(j=1; j<=n; j++){
          p3=p2;
          p2=p1;
          p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        }
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
      } while(fabs(z-z1) > eps);
      x[i-1]=xm-xl*z;
      x[n-i]=xm+xl*z;
      w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);

      w[n-i]=w[i-1];
    }
  }


  bool dirSortFcn(Directions const &a, Directions const &b){
    return b.octant > a.octant;
  }
}

/**
 * Initializes the quadrature set information for a Grid_Data object.
 * This guarantees that each <GS,DS> pair have a single originating octant.
 */
void InitDirections(Grid_Data *grid_data, Input_Variables *input_vars)
{
  std::vector<Directions> &directions = grid_data->directions;

  // Get set description from user
  int num_directions_per_octant = input_vars->num_dirs_per_dirset *
                                  input_vars->num_dirsets_per_octant;
  int num_directions = 8*num_directions_per_octant;

  // allocate storage
  directions.resize(num_directions);

  // Are we running a REAL quadrature set?
  int num_polar = input_vars->quad_num_polar;
  int num_azimuth = input_vars->quad_num_azimuthal;

  std::vector<double> polar_cos;
  std::vector<double> polar_weight;
  if(num_polar > 0){
    // make sure the user specified the correct number of quadrature points
    if(num_polar % 4 != 0){
      printf("Must have number of polar angles be a multiple of 4\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(num_azimuth % 2 != 0){
      printf("Must have number of azimuthal angles be a multiple of 2\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(num_polar*num_azimuth != num_directions){
      printf("You need to specify %d total directions, not %d\n",
          num_polar*num_azimuth, num_directions);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Compute gauss legendre weights
    polar_cos.resize(num_polar);
    polar_weight.resize(num_polar);
    GaussLegendre(-1.0, 1.0, polar_cos, polar_weight, DBL_EPSILON);

    // compute azmuhtal angles and weights
    std::vector<double> az_angle(num_azimuth);
    std::vector<double> az_weight(num_azimuth);
    double dangle = M_PI*2.0/((double) num_azimuth);
    az_angle[0] = dangle/2.0;
    for(int i=1; i<num_azimuth; i++){
      az_angle[i] = az_angle[i-1] + dangle;
    }
    for(int i=0; i<num_azimuth; i++){
      az_weight[i] = dangle;
    }

    // Loop over polar 'octants
    int d = 0;
    double total_w = 0;
    for(int i=0; i< num_polar; i++){
      for(int j=0; j< num_azimuth; j++){
        directions[d].xcos = sqrt(1.0-polar_cos[i]*polar_cos[i]) * cos(az_angle[j]);
        directions[d].ycos = sqrt(1.0-polar_cos[i]*polar_cos[i]) * sin(az_angle[j]);
        directions[d].zcos = polar_cos[i];
        directions[d].w = polar_weight[i]*az_weight[j];
        total_w += directions[d].w;

        directions[d].id = (directions[d].xcos > 0.) ? 1 : -1;
        directions[d].jd = (directions[d].ycos > 0.) ? 1 : -1;
        directions[d].kd = (directions[d].zcos > 0.) ? 1 : -1;

        directions[d].octant = 0;
        if(directions[d].id == -1){
          directions[d].octant += 1;
        }
        if(directions[d].jd == -1){
          directions[d].octant += 2;
        }
        if(directions[d].kd == -1){
          directions[d].octant += 4;
        }

        ++ d;
      }
    }
    // normalize weights
    for(int d = 0;d < num_directions;++ d){
      directions[d].w /= total_w;
    }

    // Sort by octant.. so each set has same directions
    std::sort(directions.begin(), directions.end(), dirSortFcn);
  }
  else{
    // Do (essentialy) an S2 quadrature.. but with multiple directions
    // co-located
    int d = 0;
    for(int octant = 0;octant < 8;++ octant){
      double omegas[3];
      omegas[0] = octant & 0x1;
      omegas[1] = (octant>>1) & 0x1;
      omegas[2] = (octant>>2) & 0x1;

      for(int sd=0; sd<num_directions_per_octant; sd++, d++){
        // Store which logical direction of travel we have
        directions[d].id = (omegas[0] > 0.) ? 1 : -1;
        directions[d].jd = (omegas[1] > 0.) ? 1 : -1;
        directions[d].kd = (omegas[2] > 0.) ? 1 : -1;

        // Store quadrature point's weight
        directions[d].w = 1.0 / (double)num_directions;

        // Get the direction of this quadrature point
        double theta = M_PI/4;
        double omega = M_PI/4;

        // Compute x,y,z cosine values
        double mu  = cos(theta);
        double eta = sqrt(1-mu*mu) * cos(omega);
        double xi  = sqrt(1-mu*mu) * sin(omega);
        directions[d].xcos = mu;
        directions[d].ycos = eta;
        directions[d].zcos = xi;
      }
    }
  }
}




