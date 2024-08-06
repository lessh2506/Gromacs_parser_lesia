#ifndef _NEIGHBORS_H_
#define _NEIGHBORS_H_

#ifdef DOUBLE
#define REAL double
#else
#define REAL float
#endif

typedef REAL rvect[3];
typedef int  ivect[3]; 

typedef struct {
  int n;
  int maxn;
  int *idx;
  ivect *images;
  REAL *dists;
} atomlist;

typedef struct {
  rvect box[3];
  REAL  rcut;
  ivect gridsize;
  int   atomnum;
  rvect *coords;
  ivect *cellidx;
  int   *incellidx;
  atomlist *cells;
  atomlist *neighbors;
  char *nghstatus;
} grid;

extern grid *new_grid(void);
extern void free_grid(grid *g);
extern void build_grid(rvect box[3], REAL rcut, int atomnum, rvect *coords, grid *g);
extern atomlist find_neighbors(int atom, grid *g);

#endif
