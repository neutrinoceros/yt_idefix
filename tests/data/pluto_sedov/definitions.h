#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
/* This comment is meant to challenge geometry parsing tests
#define  GEOMETRY                       SPHERICAL
*/
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ENRG0                          0
#define  DNST0                          1
#define  GAMMA                          2

/* [Beg] user-defined constants (do not change this line) */

#define  INITIAL_SMOOTHING              YES
// Define the code units below in cgs, some pre-defined constants can be used
#define  UNIT_LENGTH                    1000*CONST_pc
#define  UNIT_VELOCITY                  1e5
#define  UNIT_DENSITY                   CONST_mp

/* [End] user-defined constants (do not change this line) */
