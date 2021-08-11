#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GAMMA (5.0 / 3.0)
#define PI 3.14159
#define min2(a, b) (a) < (b) ? (a) : (b)
#define max2(a, b) (a) > (b) ? (a) : (b)

// Convert primitive quantities (density, velocity, and pressure) to conserved
// quantities (density, momentum, and energy).
void primitive_to_conserved(double *prim, double *cons) {
  double density  = prim[0];
  double velocity = prim[1];
  double pressure = prim[2];

  double momentum = velocity * density;
  double energy   = (0.5 * density * velocity * velocity) +
                    (pressure / (GAMMA - 1));

  cons[0] = density;
  cons[1] = momentum;
  cons[2] = energy;
}

// Convert conserved quantities (density, momentum, and energy) to primitive
// quantities (density, velocity, and pressure).
void conserved_to_primitive(double *cons, double *prim) {
  double density  = cons[0];
  double momentum = cons[1];
  double energy   = cons[2];

  double velocity = momentum / density;
  double pressure = (energy - (0.5 * density * velocity * velocity)) *
                    (GAMMA - 1);

  prim[0] = density;
  prim[1] = velocity;
  prim[2] = pressure;
}

// Compute the sound speed of the medium.
double cs(double *prim) {
  double density  = prim[0];
  double pressure = prim[2];

  return sqrt(GAMMA*pressure/density);
}

// Compute outer wave speeds.
void speed(double *ul, double *ur, double *s) {
	double pl[3];
	double pr[3];

	conserved_to_primitive(ul, pl);
	conserved_to_primitive(ur, pr);

	double vl = pl[1];
	double vr = pr[1];

	s[0] = min2((vl - cs(pl)), (vr - cs(pr)));
  s[1] = max2((vl + cs(pl)), (vr + cs(pr)));
}

// Compute wave speed in star region.
double speed_star(double *ul, double *ur, double sl, double sr) {
  double pl[3];
	double pr[3];

	conserved_to_primitive(ul, pl);
	conserved_to_primitive(ur, pr);

	double rhol = pl[0];
	double vl   = pl[1];
	double prel = pl[2];
	double rhor = pr[0];
	double vr   = pr[1];
	double prer = pr[2];

	double num = prer - prel + rhol * vl * (sl - vl) - rhor * vr * (sr - vr);
	double den = rhol * (sl - vl) - rhor * (sr - vr);

	return num / den;
}

// Compute cell flux vector.
void flux_vector(double *cons, double *flux) {
	double prim[3];
	conserved_to_primitive(cons, prim);

  double rho = prim[0];
	double v = prim[1];
	double pre = prim[2];

	flux[0] = rho * v;
	flux[1] = rho * v * v + pre;
	flux[2] = v * (pre + 0.5 * rho * v * v + pre / (GAMMA - 1));
}

// Compute source vector.
void source_vector(double *cons, double *source, double r) {
	double prim[3];
	conserved_to_primitive(cons, prim);

	double pre = prim[2];

	source[0] = 0;
	source[1] = 2 * pre / r;
	source[2] = 0;
}

// Compute flux vector in star region.
void flux_star_vector(double *cons, double *d_star, double *flux, double s_star,
                      double s, double *flux_star) {
  double prim[3];
	conserved_to_primitive(cons, prim);

  double rho = prim[0];
	double v   = prim[1];
	double pre = prim[2];

  for (int i = 0; i < 3; ++i) {
    flux_star[i] = (s_star * (s * cons[i] - flux[i]) + s * (pre + rho*(s - v) *
                   (s_star - v)) * d_star[i]) / (s - s_star);
  }
}

// Return the sign of a value.
double sign(double a) {
  if (a < 0) {
    return -1;
  } else if (a > 0) {
    return 1;
  } else {
    return 0;
  }
}

// Select the flattest slope from a set of three.
double minmod(double a, double b, double c) {
  return fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) * (min2(fabs(a), (min2(fabs(b), fabs(c))))) / 4;
}

// Compute HLL interface flux.
void hll_flux(double *ul2, double *ul1, double *ur1, double *ur2, double *flux_half, double plm_theta, double *s) {
  double pl2[3];
  double pl1[3];
  double pr1[3];
  double pr2[3];
  conserved_to_primitive(ul2, pl2);
  conserved_to_primitive(ul1, pl1);
  conserved_to_primitive(ur1, pr1);
  conserved_to_primitive(ur2, pr2);
  double pl[3];
  double pr[3];
  double ul[3];
  double ur[3];

  for (int i = 0; i < 3; ++i) {
    pl[i] = pl1[i] + 0.5 * minmod(plm_theta * (pl1[i] - pl2[i]),
                                  0.5 * (pr1[i] - pl2[i]),
                                  plm_theta * (pr1[i] - pl1[i]));
    pr[i] = pr1[i] - 0.5 * minmod(plm_theta * (pr1[i] - pl1[i]),
                                  0.5 * (pr2[i] - pl1[i]),
                                  plm_theta * (pr2[i] - pr1[i]));
  }

  primitive_to_conserved(pl, ul);
  primitive_to_conserved(pr, ur);

  speed(ul, ur, s);
  double s_l       = s[0];
  double s_r       = s[1];

  double fl[3];
  double fr[3];
  flux_vector(ul, fl);
	flux_vector(ur, fr);

  if (0 < s_l) {
	  flux_half[0] = fl[0];
    flux_half[1] = fl[1];
    flux_half[2] = fl[2];
  } else if ((s_l <= 0) && (0 < s_r)) {
	  flux_half[0] = (s_r * fl[0] - s_l * fr[0] + s_l * s_r * (ur[0] - ul[0])) / (s_r - s_l);
    flux_half[1] = (s_r * fl[1] - s_l * fr[1] + s_l * s_r * (ur[1] - ul[1])) / (s_r - s_l);
    flux_half[2] = (s_r * fl[2] - s_l * fr[2] + s_l * s_r * (ur[2] - ul[2])) / (s_r - s_l);
  } else {
	  flux_half[0] = fr[0];
    flux_half[1] = fr[1];
    flux_half[2] = fr[2];
  }
}

// Compute HLLC interface flux.
void hllc_flux(double *ul2, double *ul1, double *ur1, double *ur2, double *flux_half, double plm_theta, double *s) {
  double pl2[3];
  double pl1[3];
  double pr1[3];
  double pr2[3];
  conserved_to_primitive(ul2, pl2);
  conserved_to_primitive(ul1, pl1);
  conserved_to_primitive(ur1, pr1);
  conserved_to_primitive(ur2, pr2);
  double pl[3];
  double pr[3];
  double ul[3];
  double ur[3];

  for (int i = 0; i < 3; ++i) {
    pl[i] = pl1[i] + 0.5 * minmod(plm_theta * (pl1[i] - pl2[i]),
                                  0.5 * (pr1[i] - pl2[i]),
                                  plm_theta * (pr1[i] - pl1[i]));
    pr[i] = pr1[i] - 0.5 * minmod(plm_theta * (pr1[i] - pl1[i]),
                                  0.5 * (pr2[i] - pl1[i]),
                                  plm_theta * (pr2[i] - pr1[i]));
  }

  primitive_to_conserved(pl, ul);
  primitive_to_conserved(pr, ur);

  speed(ul, ur, s);
	double s_l       = s[0];
  double s_r       = s[1];
	double s_star    = speed_star(ul, ur, s_l, s_r);
	double d_star[3] = {0, 1, s_star};
  double fl[3];
  double fr[3];
  double fl_star[3];
  double fr_star[3];

	flux_vector(ul, fl);
	flux_vector(ur, fr);
	flux_star_vector(ul, d_star, fl, s_star, s_l, fl_star);
	flux_star_vector(ur, d_star, fr, s_star, s_r, fr_star);

	if (0 <= s_l) {
	  flux_half[0] = fl[0];
    flux_half[1] = fl[1];
    flux_half[2] = fl[2];
  } else if ((s_l <= 0) && (0 < s_star)) {
	  flux_half[0] = fl_star[0];
    flux_half[1] = fl_star[1];
    flux_half[2] = fl_star[2];
  } else if ((s_star <= 0) && (0 < s_r)) {
	  flux_half[0] = fr_star[0];
    flux_half[1] = fr_star[1];
    flux_half[2] = fr_star[2];
  } else {
	  flux_half[0] = fr[0];
    flux_half[1] = fr[1];
    flux_half[2] = fr[2];
  }
}

double gaussian(double x) {
  double mu = 0.5;
  double sigma = 0.2;
  return (1/sigma/sqrt(2*PI)) * exp(-(x-mu)*(x-mu)/2/sigma/sigma);
}

// Establish initial conditions for primitive quantities.
void initialize_primitive(double *primitive, double dx, int n, double x0) {
  for (int i = 0; i < n; ++i) {
    double x = x0 + (i + 0.5) * dx;
    double *prim = &primitive[3*i];
    prim[0] = 10 / x / x;
    prim[1] = 2;
    prim[2] = 0.1 / x / x;
  }
}

int main() {
  const double tmax  = 1;              // Final time
  const int n        = 1000;           // Number of spatial cells
  const double rl    = 1;              // Left spatial boundary
  const double rr    = 10;             // Right spatial boundary
  const double dr    = (rr - rl) / n;  // Size of each cell
  const double chkpt = 0.001;          // Time intervals over which data is collected
  double dt          = 0.00001;        // Initial timestep size (varied later)
  double t           = 0;              // Initial time
  double plm_theta   = 1.5;            // Slope limiter parameter
  int j              = 0;              // Initial value for checkpoint counter
  double a           = 0;              // Dummy variable for fastest wave speed
  double cfl_number  = 0.4;            // Convergence condition parameter
  double s[2];                         // Outermost wave speed array

  // Initialize array of primitives
  double primitive[3*n];

  // Initialize 4 separate arrays of conserved quantities for RK3 in time
  double conserved[3*n];
  double conserved1[3*n];
  double conserved2[3*n];
  double conserved3[3*n];

  // Implement initial conditions
  initialize_primitive(primitive, dr, n, rl);
  for (int i = 0; i < n; ++i) {
    double *prim = &primitive[3*i];
    double *cons = &conserved[3*i];
    double *cons1 = &conserved1[3*i];
    double *cons2 = &conserved2[3*i];
    primitive_to_conserved(prim, cons);
    primitive_to_conserved(prim, cons1);
    primitive_to_conserved(prim, cons2);
  }

  // Evolve the simulation in time (third order via RK3).
  while (t < tmax) {
    // Dummy variables for fastest wave speeds in each RK3 cycle
    double a1 = 0;
    double a2 = 0;
    double a3 = 0;

    // Update the simulation in space (higher order via slope limiting).
    for (int i = 2; i < (n-2); ++i) {
      double *cons_im2 = &conserved[3*(i-2)];
			double *cons_im1 = &conserved[3*(i-1)];
    	double *cons_i00 = &conserved[3*(i+0)];
    	double *cons_ip1 = &conserved[3*(i+1)];
      double *cons_ip2 = &conserved[3*(i+2)];
      double f_iph[3];
      double f_imh[3];
			hllc_flux(cons_im2, cons_im1, cons_i00, cons_ip1, f_imh, plm_theta, s);
      hllc_flux(cons_im1, cons_i00, cons_ip1, cons_ip2, f_iph, plm_theta, s);

      double r_i = rl + (i + 0.5) * dr;
      double r_imh = r_i - (0.5 * dr);
      double r_iph = r_i + (0.5 * dr);

      double source[3];
      source_vector(cons_i00, source, r_i);

			conserved1[3*i+0] = conserved[3*i+0] - (r_iph * r_iph * f_iph[0] - r_imh * r_imh * f_imh[0]) * dt / dr / r_i / r_i + dt * source[0];
      conserved1[3*i+1] = conserved[3*i+1] - (r_iph * r_iph * f_iph[1] - r_imh * r_imh * f_imh[1]) * dt / dr / r_i / r_i + dt * source[1];
      conserved1[3*i+2] = conserved[3*i+2] - (r_iph * r_iph * f_iph[2] - r_imh * r_imh * f_imh[2]) * dt / dr / r_i / r_i + dt * source[2];

      a  = max2(abs(s[0]), abs(s[1]));
      a1 = max2(a, abs(a1));
    }

    // Inflow boundary conditions (repeated for each RK3 cycle)
    for (int i = 0; i < 2; ++i) {
      conserved1[3*i]   = conserved[3*i];
      conserved1[3*i+1] = conserved[3*i+1];
      conserved1[3*i+2] = conserved[3*i+2];
    }

    // Outflow boundary conditions (repeated for each RK3 cycle)
    for (int i = 0; i < 2; ++i) {
      conserved1[3*(n-i)-1] = conserved1[3*n-7];
      conserved1[3*(n-i)-2] = conserved1[3*n-8];
      conserved1[3*(n-i)-3] = conserved1[3*n-9];
    }

    // Periodic boundary conditions (repeated for each RK3 cycle)
    /* for (int i = 0; i < 2; ++i) {
      conserved1[3*i+0] = conserved1[3*(n-3+i)-3];
      conserved1[3*i+1] = conserved1[3*(n-3+i)-2];
      conserved1[3*i+2] = conserved1[3*(n-3+i)-1];

      conserved1[3*(n-1+i)-3] = conserved1[3*(i+2)+0];
      conserved1[3*(n-1+i)-2] = conserved1[3*(i+2)+1];
      conserved1[3*(n-1+i)-1] = conserved1[3*(i+2)+2];
    } */

    for (int i = 2; i < (n-2); ++i) {
      double *cons_im2_1 = &conserved1[3*(i-2)];
      double *cons_im1_1 = &conserved1[3*(i-1)];
      double *cons_i00_1 = &conserved1[3*(i+0)];
      double *cons_ip1_1 = &conserved1[3*(i+1)];
      double *cons_ip2_1 = &conserved1[3*(i+2)];
      double f_iph1[3];
      double f_imh1[3];
      hllc_flux(cons_im2_1, cons_im1_1, cons_i00_1, cons_ip1_1, f_imh1, plm_theta, s);
      hllc_flux(cons_im1_1, cons_i00_1, cons_ip1_1, cons_ip2_1, f_iph1, plm_theta, s);

      double r_i = rl + (i + 0.5) * dr;
      double r_imh = r_i - (0.5 * dr);
      double r_iph = r_i + (0.5 * dr);

      double source1[3];
      source_vector(cons_i00_1, source1, r_i);

      conserved2[3*i+0] = 3 * conserved[3*i+0] / 4 + conserved1[3*i+0] / 4 -
                          (r_iph * r_iph * f_iph1[0] - r_imh * r_imh * f_imh1[0]) * dt / dr / r_i / r_i / 4 + dt * source1[0] / 4;
      conserved2[3*i+1] = 3 * conserved[3*i+1] / 4 + conserved1[3*i+1] / 4 -
                          (r_iph * r_iph * f_iph1[1] - r_imh * r_imh * f_imh1[1]) * dt / dr / r_i / r_i / 4 + dt * source1[1] / 4;
      conserved2[3*i+2] = 3 * conserved[3*i+2] / 4 + conserved1[3*i+2] / 4 -
                          (r_iph * r_iph * f_iph1[2] - r_imh * r_imh * f_imh1[2]) * dt / dr / r_i / r_i / 4 + dt * source1[2] / 4;

      a = max2(abs(s[0]), abs(s[1]));
      a2 = max2(a, abs(a2));
    }

    for (int i = 0; i < 2; ++i) {
      conserved2[3*i]   = conserved[3*i];
      conserved2[3*i+1] = conserved[3*i+1];
      conserved2[3*i+2] = conserved[3*i+2];
    }

    for (int i = 0; i < 2; ++i) {
      conserved2[3*(n-i)-1] = conserved2[3*n-7];
      conserved2[3*(n-i)-2] = conserved2[3*n-8];
      conserved2[3*(n-i)-3] = conserved2[3*n-9];
    }

    /* for (int i = 0; i < 2; ++i) {
      conserved2[3*i+0] = conserved2[3*(n-3+i)-3];
      conserved2[3*i+1] = conserved2[3*(n-3+i)-2];
      conserved2[3*i+2] = conserved2[3*(n-3+i)-1];

      conserved2[3*(n-1+i)-3] = conserved2[3*(i+2)+0];
      conserved2[3*(n-1+i)-2] = conserved2[3*(i+2)+1];
      conserved2[3*(n-1+i)-1] = conserved2[3*(i+2)+2];
    } */

    for (int i = 2; i < (n-2); ++i) {
      double *cons_im2_2 = &conserved2[3*(i-2)];
      double *cons_im1_2 = &conserved2[3*(i-1)];
      double *cons_i00_2 = &conserved2[3*(i+0)];
      double *cons_ip1_2 = &conserved2[3*(i+1)];
      double *cons_ip2_2 = &conserved2[3*(i+2)];
      double f_iph2[3];
      double f_imh2[3];
      hllc_flux(cons_im2_2, cons_im1_2, cons_i00_2, cons_ip1_2, f_imh2, plm_theta, s);
      hllc_flux(cons_im1_2, cons_i00_2, cons_ip1_2, cons_ip2_2, f_iph2, plm_theta, s);

      double r_i = rl + (i + 0.5) * dr;
      double r_imh = r_i - (0.5 * dr);
      double r_iph = r_i + (0.5 * dr);

      double source2[3];
      source_vector(cons_i00_2, source2, r_i);

      conserved3[3*i+0] = conserved[3*i+0] / 3 + 2 * conserved2[3*i+0] / 3 -
                          2 * (r_iph * r_iph * f_iph2[0] - r_imh * r_imh * f_imh2[0]) * dt / dr / r_i / r_i / 3 + dt * 2 * source2[0] / 3;
      conserved3[3*i+1] = conserved[3*i+1] / 3 + 2 * conserved2[3*i+1] / 3 -
                          2 * (r_iph * r_iph * f_iph2[1] - r_imh * r_imh * f_imh2[1]) * dt / dr / r_i / r_i / 3 + dt * 2 * source2[1] / 3;
      conserved3[3*i+2] = conserved[3*i+2] / 3 + 2 * conserved2[3*i+2] / 3 -
                          2 * (r_iph * r_iph * f_iph2[2] - r_imh * r_imh * f_imh2[2]) * dt / dr / r_i / r_i / 3 + dt * 2 * source2[2] / 3;

      a = max2(abs(s[0]), abs(s[1]));
      a3 = max2(a, abs(a3));
    }

    for (int i = 0; i < 2; ++i) {
      conserved3[3*i]   = conserved[3*i];
      conserved3[3*i+1] = conserved[3*i+1];
      conserved3[3*i+2] = conserved[3*i+2];
    }

    for (int i = 0; i < 2; ++i) {
      conserved3[3*(n-i)-1] = conserved3[3*n-7];
      conserved3[3*(n-i)-2] = conserved3[3*n-8];
      conserved3[3*(n-i)-3] = conserved3[3*n-9];
    }

    /* for (int i = 0; i < 2; ++i) {
      conserved3[3*i+0] = conserved3[3*(n-3+i)-3];
      conserved3[3*i+1] = conserved3[3*(n-3+i)-2];
      conserved3[3*i+2] = conserved3[3*(n-3+i)-1];

      conserved3[3*(n-1+i)-3] = conserved3[3*(i+2)+0];
      conserved3[3*(n-1+i)-2] = conserved3[3*(i+2)+1];
      conserved3[3*(n-1+i)-1] = conserved3[3*(i+2)+2];
    } */

    // Pick out fastest wave speed on grid at current time
    double ws = max2(a1, a2);
    double a_max = max2(ws, a3);

    // Update timestep such that it satisfies necessary stability criteria
    dt = cfl_number * dr / a_max;

    // Save conserved vectors to text files in checkpoint intervals.
    if (t >= (chkpt * j)) {
      FILE *fp;
      char filepath[256];
      snprintf (filepath, sizeof(filepath), "output/data%d.txt", j);
      fp = fopen(filepath, "w");
      for (int k = 0; k < n; ++k) {
        fprintf(fp, "%f %f %f %f\n", rl + (k + 0.5) * dr, conserved[3*k],
                conserved[3*k+1], conserved[3*k+2]);
      }
      fclose(fp);
      j += 1;
    } else {
      ;
    }

    // Update array of conserved quantities for next timestep
    for (int i = 0; i < 3*n; ++i) {
      conserved[i] = conserved3[i];
    }

    t += dt;
  }

  return 0;
}
