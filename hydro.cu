#include <stdio.h>
#include <math.h>

#define GAMMA (5.0 / 3.0)
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

// Compute HLLC interface flux.
void hllc_flux(double *ul, double *ur, double *flux_half){
  double s[2];
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

// Establish initial conditions for primitive quantities.
void initialize_primitive(double *primitive, double dx, int n, double x0) {
  for (int i = 0; i < n; ++i) {
    double x = x0 + (i + 0.5) * dx;
    double *prim = &primitive[3*i];
    if (x < 0.5) {
      prim[0] = 1;
      prim[1] = 0;
      prim[2] = 1;
    } else {
      prim[0] = 0.1;
      prim[1] = 0;
      prim[2] = 0.125;
    }
  }
}

int main() {
  const double tmax  = 0.5;
  const int n        = 1000;
  const double xl    = 0;
  const double xr    = 1;
  const double dx    = (xr - xl) / n;
  const double dt    = 0.00025;
  const double chkpt = 0.0025;

  double primitive[3*n];
  double conserved[3*n];
  double conserved1[3*n];

  initialize_primitive(primitive, dx, n, xl);
  for (int i = 0; i < n; ++i) {
    double *prim = &primitive[3*i];
    double *cons = &conserved[3*i];
    primitive_to_conserved(prim, cons);
  }

  double t = 0;
  int j = 0;

  // Evolve the simulation in time.
  while (t < tmax) {
    // Update the simulation in space.
    for (int i = 1; i < (n-1); ++i) {
			double *cons_im1 = &conserved[3*(i-1)];
    	double *cons_i   = &conserved[3*(i)];
    	double *cons_ip1 = &conserved[3*(i+1)];

      double f_iph[3];
      double f_imh[3];
			hllc_flux(cons_im1, cons_i, f_imh);
      hllc_flux(cons_i, cons_ip1, f_iph);

			conserved1[3*i]   = conserved[3*i] - (f_iph[0] - f_imh[0]) * dt / dx;
      conserved1[3*i+1] = conserved[3*i+1] - (f_iph[1] - f_imh[1]) * dt / dx;
      conserved1[3*i+2] = conserved[3*i+2] - (f_iph[2] - f_imh[2]) * dt / dx;
    }

    for (int i = 3; i < (3*(n-1)); ++i) {
			conserved[i] = conserved1[i];
    }

    // Save conserved vectors to text files in checkpoint intervals.
    if (t >= (chkpt * j)) {
      FILE *fp;
      char filepath[256];
      snprintf (filepath, sizeof(filepath), "output/data%d.txt", j);
      fp = fopen(filepath, "w");
      for (int k = 0; k < n; ++k) {
        fprintf(fp, "%f %f %f %f\n", (k + 0.5) * dx, conserved[3*k],
                conserved[3*k+1], conserved[3*k+2]);
      }
      fclose(fp);
      j += 1;
    } else {
      ;
    }

    t += dt;
  }

  return 0;
}
