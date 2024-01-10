#ifndef UNIT_H
#define UNIT_H

#define Max(a,b) ((a) > (b) ? (a) : (b))
#define Min(a,b) ((a) < (b) ? (a) : (b))
/* length unit */
#define DM2M 1.0E-1
#define CM2M 1.0E-2
#define MM2M 1.0E-3
#define UM2M 1.0E-6
#define NM2M 1.0E-9
#define AM2M 1.0E-10

/* time unit */
#define S2MS 1.0E+3
#define S2US 1.0E+6
#define S2NS 1.0E+9
#define S2PS 1.0E+12

#define M2A (1 / A)
#define DM2A (DM2M / A)
#define CM2A (CM2M / A)
#define MM2A (MM2M / A)
#define UM2A (UM2M / A)
#define NM2A (NM2M / A)
#define AM2A (AM2M / A)

#endif