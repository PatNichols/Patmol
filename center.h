#ifndef _CENTER_H_
#define _CENTER_H_
#include <stdio.h>

typedef struct {
    double q, r[3];
} Center;

#define center_position(cen) (cen).r
#define center_charge(cen) (cen).q
#define center_x_position(cen) (cen).r[0]
#define center_y_position(cen) (cen).r[1]
#define center_z_position(cen) (cen).r[2]

extern void read_center ( FILE * infile, Center * cen );
extern void write_center ( FILE * outfile, const Center * cen );
#endif
