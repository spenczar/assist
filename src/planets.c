#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spk.h"
#include "planets.h"
#include "assist.h"

/*
 *  assist_jpl_work
 *
 *  Interpolate the appropriate Chebyshev polynomial coefficients.
 *
 *      ncf - number of coefficients per component
 *      ncm - number of components (ie: 3 for most)
 *      niv - number of intervals / sets of coefficients
 *
 */

void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w)
{
        double T[24], S[24];
        double U[24];
        double t, c;
        int p, m, n, b;

        // adjust to correct interval
        t = t0 * (double)niv;
        t0 = 2.0 * fmod(t, 1.0) - 1.0;
        c = (double)(niv * 2) / t1 / 86400.0;

        b = (int)t;

        // set up Chebyshev polynomials and derivatives
        T[0] = 1.0; T[1] = t0;
        S[0] = 0.0; S[1] = 1.0;
        U[0] = 0.0; U[1] = 0.0;	U[2] = 4.0;

        for (p = 2; p < ncf; p++) {
                T[p] = 2.0 * t0 * T[p-1] - T[p-2];
                S[p] = 2.0 * t0 * S[p-1] + 2.0 * T[p-1] - S[p-2];
        }
        for (p = 3; p < ncf; p++) {
                U[p] = 2.0 * t0 * U[p-1] + 4.0 * S[p-1] - U[p-2];
        }

        // compute the position/velocity
        for (m = 0; m < ncm; m++) {
                u[m] = v[m] = w[m] = 0.0;
                n = ncf * (m + b * ncm);

                for (p = 0; p < ncf; p++) {
                        u[m] += T[p] * P[n+p];
                        v[m] += S[p] * P[n+p] * c;
                        w[m] += U[p] * P[n+p] * c * c;
                }
        }
}

/*
 *	assist_jpl_position
 *
 *	Interpolate the appropriate Chebyshev polynomial coefficients to get position only,
 *	without velocity or acceleration.
 */
void assist_jpl_position(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u)
{
  double T[24];
  double t;
  int p, m, n, b;

  // adjust to correct interval
  t = t0 * (double)niv;
  t0 = 2.0 * fmod(t, 1.0) - 1.0;

  b = (int)t;

  // set up Chebyshev polynomials and derivatives
  T[0] = 1.0; T[1] = t0;

  for (p = 2; p < ncf; p++) {
    T[p] = 2.0 * t0 * T[p-1] - T[p-2];
  }

  // compute the position
  for (m = 0; m < ncm; m++) {
    u[m] = 0.0;
    n = ncf * (m + b * ncm);

    for (p = 0; p < ncf; p++) {
      u[m] += T[p] * P[n+p];
    }
  }
}
 
/*
 *  assist_jpl_init
 *
 *  Initialise everything needed ... probaly not be compatible with a non-430 file.
 *
 */

static double getConstant(struct jpl_s* jpl, char* name){
    for (int p = 0; p < jpl->num; p++) {
        if (strncmp(name,jpl->str[p],6)==0){
            return jpl->con[p];
        }
    }
    fprintf(stderr,"WARNING: Constant [%s] not found in ephemeris file.\n",name); 
    return 0;
}

struct jpl_s * assist_jpl_init(char *str)
{
    struct stat sb;
    ssize_t ret;
    int fd;

    if ((fd = open(str, O_RDONLY)) < 0){
        return NULL;
    }

    if (fstat(fd, &sb) < 0){
        close(fd);
        fprintf(stderr, "Error while trying to determine filesize.\n");
        return NULL;
    }
    

    // skip the header and constant names for now
    if (lseek(fd, 0x0A5C, SEEK_SET) < 0){
        close(fd);
        fprintf(stderr, "Error while seeking to header.\n");
        return NULL;
    }

    struct jpl_s* jpl = calloc(1, sizeof(struct jpl_s));

    // read header
    ret  = read(fd, &jpl->beg, sizeof(double));     // Start JD
    ret += read(fd, &jpl->end, sizeof(double));     // End JD
    ret += read(fd, &jpl->inc, sizeof(double));     // Days per block
    ret += read(fd, &jpl->num, sizeof(int32_t));    // Number of constants
    ret += read(fd, &jpl->cau, sizeof(double));     // AU to km 
    ret += read(fd, &jpl->cem, sizeof(double));     // Ratio between Earth/Moon

    // number of coefficients for all components
    for (int p = 0; p < JPL_N; p++){
        jpl->ncm[p] = 3;
    }
    // exceptions:
    jpl->ncm[JPL_NUT] = 2; // nutations
    jpl->ncm[JPL_TDB] = 1; // TT-TDB

    for (int p = 0; p < 12; p++) {                      // Columns 1-12 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    ret += read(fd, &jpl->ver,     sizeof(int32_t));    // Version. e.g. 440
    ret += read(fd, &jpl->off[12], sizeof(int32_t));    // Columns 13 of Group 1050
    ret += read(fd, &jpl->ncf[12], sizeof(int32_t));
    ret += read(fd, &jpl->niv[12], sizeof(int32_t));

    // Get all the constant names
    jpl->str = calloc(jpl->num, sizeof(char *));

    // retrieve the names of the first 400 constants
    lseek(fd, 0x00FC, SEEK_SET);    
    for (int p = 0; p < 400; p++) {     // Group 1040
        jpl->str[p] = calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    // read the remaining constant names
    lseek(fd, 0x0B28, SEEK_SET);
    for (int p = 400; p < jpl->num; p++) {
        jpl->str[p] = calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    for (int p = 13; p < 15; p++) {                     // Columns 14 and 15 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    // adjust for correct indexing (ie: zero based)
    for (int p = 0; p < JPL_N; p++){
        jpl->off[p] -= 1;
    }

    // save file size, and determine 'kernel size' or 'block size' (=8144 bytes for DE440/441)
    jpl->len = sb.st_size;
    jpl->rec = sizeof(double) * 2;

    for (int p = 0; p < JPL_N; p++){
        jpl->rec += sizeof(double) * jpl->ncf[p] * jpl->niv[p] * jpl->ncm[p];
    }

    // memory map the file, which makes us thread-safe with kernel caching
    jpl->map = mmap(NULL, jpl->len, PROT_READ, MAP_SHARED, fd, 0);

    if (jpl->map == NULL){ 
        close(fd);
        free(jpl); // note constants leak
        fprintf(stderr, "Error while calling mmap().\n");
        return NULL;
    }

    // Read constants
    jpl->con = calloc(jpl->num, sizeof(double));
    lseek(fd, jpl->rec, SEEK_SET); // Starts at offset of 1 block size
    for (int p = 0; p < jpl->num; p++){
        read(fd, &jpl->con[p], sizeof(double));
        //printf("%6d  %s   %.5e\n",p,jpl->str[p],jpl->con[p]);
    }


    // Find masses
    jpl->mass[0] = getConstant(jpl, "GMS   ");  // Sun 
    jpl->mass[1] = getConstant(jpl, "GM1   ");  // Mercury
    jpl->mass[2] = getConstant(jpl, "GM2   ");
	double emrat = getConstant(jpl, "EMRAT  "); // Earth Moon Ratio
    double gmb = getConstant(jpl, "GMB   ");    // Earth Moon combined
    jpl->mass[3] = (emrat/(1.+emrat)) * gmb;    // Earth 
    jpl->mass[4] = 1./(1+emrat) * gmb;          // Moon 
    jpl->mass[5] = getConstant(jpl, "GM4   ");  // Mars
    jpl->mass[6] = getConstant(jpl, "GM5   ");  // Jupiter
    jpl->mass[7] = getConstant(jpl, "GM6   ");
    jpl->mass[8] = getConstant(jpl, "GM7   ");
    jpl->mass[9] = getConstant(jpl, "GM8   ");
    jpl->mass[10] = getConstant(jpl, "GM9   "); // Pluto


    // Other constants
    jpl->J2E = getConstant(jpl, "J2E   ");
    jpl->J3E = getConstant(jpl, "J3E   ");
    jpl->J4E = getConstant(jpl, "J4E   ");
    jpl->J2SUN = getConstant(jpl, "J2SUN ");
    jpl->AU = getConstant(jpl, "AU    ");
    jpl->RE = getConstant(jpl, "RE    ");
    jpl->CLIGHT = getConstant(jpl, "CLIGHT");
    jpl->ASUN = getConstant(jpl, "ASUN  ");

    // this file descriptor is no longer needed since we are memory mapped
    if (close(fd) < 0) { 
        fprintf(stderr, "Error while closing file.\n");
    }
#if defined(MADV_RANDOM)
    if (madvise(jpl->map, jpl->len, MADV_RANDOM) < 0){
        fprintf(stderr, "Error during madvise.\n");
    }
#endif

    return jpl;

}

void assist_jpl_free(struct jpl_s *jpl) {
    if (jpl == NULL){
        return;
    }
    if (munmap(jpl->map, jpl->len) < 0){ 
        fprintf(stderr, "Error during munmap().\n");
    }
    for (int p = 0; p < jpl->num; p++){
        free(jpl->str[p]);
    }
    free(jpl->str);
    free(jpl->con);
    free(jpl);
}

/*
 *  jpl_calc
 *
 *  Caculate the position+velocity in _equatorial_ coordinates.
 *  Assumes pos is initially zero.
 */
enum ASSIST_STATUS assist_jpl_calc(struct jpl_s *jpl, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const out_x, double* const out_y, double* const out_z,
		 double* const out_vx, double* const out_vy, double* const out_vz,
         double* const out_ax, double* const out_ay, double* const out_az){
    if (jpl == NULL || jpl->map == NULL)
        return ASSIST_ERROR_EPHEM_FILE;
    if(body<0 || body >= ASSIST_BODY_NPLANETS)
	    return(ASSIST_ERROR_NEPHEM);
    
    // check if covered by this file
    if (jd_ref + jd_rel < jpl->beg || jd_ref + jd_rel > jpl->end)
      return ASSIST_ERROR_COVERAGE;

    struct mpos_s pos;

    // Get mass, position, velocity, and mass of body i in barycentric coords.
    *GM = jpl->mass[body];

    // special case for earth and the moon, since JPL stores the barycenter only.
    if (body == ASSIST_BODY_EARTH || body == ASSIST_BODY_MOON) {
      struct mpos_s emb, lun;
      struct jpl_record rec_emb, rec_lun;
      rec_emb = assist_jpl_get_record(jpl, jd_ref, jd_rel, JPL_EMB);
      rec_lun = assist_jpl_get_record(jpl, jd_ref, jd_rel, JPL_LUN);
      assist_jpl_work(rec_emb.data, rec_emb.ncm, rec_emb.ncf, rec_emb.niv, rec_emb.t, jpl->inc, emb.u, emb.v, emb.w); // earth moon barycenter
      assist_jpl_work(rec_lun.data, rec_lun.ncm, rec_lun.ncf, rec_lun.niv, rec_lun.t, jpl->inc, lun.u, lun.v, lun.w); // moon

      vecpos_set(pos.u, emb.u);
      vecpos_set(pos.v, emb.v);
      vecpos_set(pos.w, emb.w);
      if (body == ASSIST_BODY_EARTH) {
	vecpos_off(pos.u, lun.u, -1.0 / (1.0 + jpl->cem));
	vecpos_off(pos.v, lun.v, -1.0 / (1.0 + jpl->cem));
	vecpos_off(pos.w, lun.w, -1.0 / (1.0 + jpl->cem));
      } else {
	vecpos_off(pos.u, lun.u, jpl->cem / (1.0 + jpl->cem));
	vecpos_off(pos.v, lun.v, jpl->cem / (1.0 + jpl->cem));
	vecpos_off(pos.w, lun.w, jpl->cem / (1.0 + jpl->cem));
      }
    } else {
      struct jpl_record record;
      enum JPL_COL col;
      switch (body) { // The indices in the jpl-> arrays match the JPL component index for the body
      case ASSIST_BODY_SUN:
	col = 10;
	break;
      case ASSIST_BODY_MERCURY:
	col = JPL_MER;
	break;
      case ASSIST_BODY_VENUS:
	col = JPL_VEN;
	break;
      case ASSIST_BODY_MARS:
	col = JPL_MAR;
	break;
      case ASSIST_BODY_JUPITER:
	col = JPL_JUP;
	break;
      case ASSIST_BODY_SATURN:
	col = JPL_SAT;
	break;
      case ASSIST_BODY_URANUS:
	col = JPL_URA;
      	break;
      case ASSIST_BODY_NEPTUNE:
	col = JPL_NEP;
	break;
      case ASSIST_BODY_PLUTO:
	col = JPL_PLU;
	break;
      default:
	return ASSIST_ERROR_NEPHEM;
      }
      record = assist_jpl_get_record(jpl, jd_ref, jd_rel, col);

      assist_jpl_work(record.data, record.ncm, record.ncf, record.niv, record.t, jpl->inc, pos.u, pos.v, pos.w);
    }

    // Convert to au/day and au/day^2
    vecpos_div(pos.u, jpl->cau);
    vecpos_div(pos.v, jpl->cau/86400.);
    vecpos_div(pos.w, jpl->cau/(86400.*86400.));

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];
    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];
    *out_ax = pos.w[0];
    *out_ay = pos.w[1];
    *out_az = pos.w[2];

    return(ASSIST_SUCCESS);

}

/*
 * assist_jpl_record_number
 *
 * Calculate the record block number for a given date.
 */
u_int32_t assist_jpl_record_number(struct jpl_s *jpl, const double jd_ref, const double jd_rel) {
  return (u_int32_t)((jd_ref + jd_rel - jpl->beg) / jpl->inc);
}

struct jpl_record assist_jpl_get_record(struct jpl_s *jpl, const double jd_ref, const double jd_rel, enum JPL_COL col) {
  struct jpl_record rec;
  u_int32_t blk = assist_jpl_record_number(jpl, jd_ref, jd_rel);
  rec.data = &((double*)jpl->map + (blk + 2) * jpl->rec/sizeof(double))[jpl->off[col]];
  rec.ncm = jpl->ncm[col];
  rec.ncf = jpl->ncf[col];
  rec.niv = jpl->niv[col];
  rec.t = ((jd_ref - jpl->beg - (double)blk * jpl->inc) + jd_rel) / jpl->inc;
  return rec;
}

enum ASSIST_STATUS assist_helio_to_bary(struct jpl_s *jpl, double jd_ref, double jd_rel,
					double* const x,  double* const y, double* const z) {
  /*
   * Translate from heliocentric to barycentric coordinates.
   */
  if (jpl == NULL || jpl->map == NULL) {
    return ASSIST_ERROR_EPHEM_FILE;
  }

  double jd = jd_ref + jd_rel;
  if (jd < jpl->beg || jd > jpl->end) {
    return ASSIST_ERROR_COVERAGE;
  }

  struct mpos_s pos;
  // Get position of sun in barycentric coordinates.

  u_int32_t blk = (u_int32_t)((jd - jpl->beg) / jpl->inc);
  double *zz = (double*)jpl->map + (blk + 2) * jpl->rec / sizeof(double);
  double t = (jd - jpl->beg - (double)blk * jpl->inc) / jpl->inc;

  assist_jpl_position(&zz[jpl->off[JPL_SUN]], jpl->ncm[JPL_SUN], jpl->ncf[JPL_SUN], jpl->niv[JPL_SUN], t, jpl->inc, pos.u);

  vecpos_div(pos.u, jpl->cau);

  // Result is translation of input position by sun's position.
  *x = *x + pos.u[0];
  *y = *y + pos.u[1];
  *z = *z + pos.u[2];

  return (ASSIST_SUCCESS);
}
