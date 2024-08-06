#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "neighbors.h"

grid *new_grid(void) {
  grid *g;
  int i, j;

  g = (grid *) malloc(sizeof(grid));
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      g->box[i][j] = 0;
  g->rcut = 0;
  for(i=0;i<3;i++)
    g->gridsize[i] = 0;
  g->atomnum = 0;
  g->coords = (rvect *) NULL;
  g->cellidx = (ivect *) NULL;
  g->incellidx = (int *) NULL;
  g->cells = (atomlist *) NULL;
  g->neighbors = (atomlist *) NULL;
  g->nghstatus = (char *) NULL;

  return g;
}

void free_grid(grid *g) {
  int i;

  if(g->coords != (rvect *) NULL)
    free(g->coords);

  if(g->cellidx != (ivect *) NULL)
    free(g->cellidx);

  if(g->incellidx != (int *) NULL)
    free(g->incellidx);

  if(g->cells != (atomlist *) NULL) {
    for(i=0; i<g->gridsize[0]*g->gridsize[1]*g->gridsize[2]; i++) {
      free(g->cells[i].idx);
      free(g->cells[i].images);
    }
    free(g->cells);
  }

  if(g->neighbors != (atomlist *) NULL) {
    for(i=0; i<g->atomnum; i++) {
      free(g->neighbors[i].idx);
      free(g->neighbors[i].images);
      free(g->neighbors[i].dists);
    }
    free(g->neighbors);
  }

  if(g->nghstatus != (char *) NULL)
    free(g->nghstatus);

  free(g);
}

void build_grid(rvect box[3], REAL rcut, int atomnum, rvect *coords, grid *g) {
  int i, j, k, s;
  ivect image;
  int cx, cy, cz, cidx;
  REAL coarse;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      g->box[i][j] = box[i][j];

  g->rcut = rcut;

  s = g->gridsize[0]*g->gridsize[1]*g->gridsize[2];
  for(i=0;i<3;i++)
    g->gridsize[i] = (int) floor(box[i][i]/(rcut+1e-5));

  if(((REAL) atomnum)/((REAL) g->gridsize[0]*g->gridsize[1]*g->gridsize[2]) < 20.0) {
    coarse = 20.0/(((REAL) atomnum)/((REAL) g->gridsize[0]*g->gridsize[1]*g->gridsize[2]));
#ifdef DOUBLE
    coarse = pow(coarse, 1.0/3.0);
#else
    coarse = powf(coarse, 1.0/3.0);
#endif

    for(i=0;i<3;i++)
      g->gridsize[i] = (int) ceil(((REAL) g->gridsize[i])/coarse);
  }

  g->coords = (rvect *) realloc(g->coords, sizeof(rvect)*atomnum);
  for(i=0; i<atomnum; i++)
    for(j=0; j<3; j++)
      g->coords[i][j] = coords[i][j];

  g->cellidx = (ivect *) realloc(g->cellidx, sizeof(ivect)*atomnum);
  g->incellidx = (int *) realloc(g->incellidx, sizeof(int)*atomnum);

  for(i=g->gridsize[0]*g->gridsize[1]*g->gridsize[2]; i<s; i++) {
    free(g->cells[i].idx);
    free(g->cells[i].images);
  }
  g->cells = (atomlist *) realloc(g->cells, sizeof(atomlist)*g->gridsize[0]*g->gridsize[1]*g->gridsize[2]);
  for(i=s; i<g->gridsize[0]*g->gridsize[1]*g->gridsize[2]; i++) {
    g->cells[i].maxn = 100;
    g->cells[i].idx = (int *) malloc(sizeof(int)*g->cells[i].maxn);
    g->cells[i].images = (ivect *) malloc(sizeof(ivect)*g->cells[i].maxn);
  }
  for(i=0; i<g->gridsize[0]*g->gridsize[1]*g->gridsize[2]; i++)
    g->cells[i].n = 0;

  for(i=atomnum; i<g->atomnum; i++) {
    free(g->neighbors[i].idx);
    free(g->neighbors[i].images);
    free(g->neighbors[i].dists);
  }
  g->neighbors = (atomlist *) realloc(g->neighbors, sizeof(atomlist)*atomnum);
  for(i=g->atomnum; i<atomnum; i++) {
    g->neighbors[i].maxn = 100;
    g->neighbors[i].idx = (int *) malloc(sizeof(int)*g->neighbors[i].maxn);
    g->neighbors[i].images = (ivect *) malloc(sizeof(ivect)*g->neighbors[i].maxn);
    g->neighbors[i].dists = (REAL *) malloc(sizeof(REAL)*g->neighbors[i].maxn);
  }
  for(i=0; i<atomnum; i++)
    g->neighbors[i].n = 0;

  g->nghstatus = (char *) realloc(g->nghstatus, sizeof(char)*atomnum);
  for(i=0; i<atomnum; i++)
    g->nghstatus[i] = 0;

  g->atomnum = atomnum;

  for(i=0; i<g->atomnum; i++) {
    for(j=0; j<3; j++) {
      image[j] = 0;
      while(g->coords[i][j] < 0) {
        for(k=0; k<3; k++)
          g->coords[i][k] += g->box[j][k];
        image[j]++;
      }
      while(g->coords[i][j] >= g->box[j][j]) {
        for(k=0; k<3; k++)
          g->coords[i][k] -= g->box[j][k];
        image[j]--;
      }
    }

    cx = (int) floor(g->coords[i][0]/(g->box[0][0]/g->gridsize[0]));
    cy = (int) floor(g->coords[i][1]/(g->box[1][1]/g->gridsize[1]));
    cz = (int) floor(g->coords[i][2]/(g->box[2][2]/g->gridsize[2]));

    cidx = cz*g->gridsize[0]*g->gridsize[1]+cy*g->gridsize[0]+cx;
    g->cells[cidx].n++;
    if(g->cells[cidx].n == g->cells[cidx].maxn) {
      g->cells[cidx].maxn = (int) (((REAL) g->cells[cidx].maxn)*1.5);
      g->cells[cidx].idx = (int *) realloc(g->cells[cidx].idx, sizeof(int)*g->cells[cidx].maxn);
      g->cells[cidx].images = (ivect *) realloc(g->cells[cidx].images, sizeof(ivect)*g->cells[cidx].maxn);
    }
    g->cells[cidx].idx[g->cells[cidx].n-1] = i;
    for(j=0; j<3; j++)
       g->cells[cidx].images[g->cells[cidx].n-1][j] = image[j];

    g->cellidx[i][0] = cx;
    g->cellidx[i][1] = cy;
    g->cellidx[i][2] = cz;
    g->incellidx[i] = g->cells[cidx].n-1;
  }
}

atomlist find_neighbors(int atom, grid *g) {
  int i, j, k, l;
  int cx, cy, cz, cidx, cidx2;
  int sx, sy, sz;
  int num;
  REAL vx, vy, vz, d;


  if(g->nghstatus[atom] == 1)
    return g->neighbors[atom];

  cidx = g->cellidx[atom][2]*g->gridsize[0]*g->gridsize[1] + 
         g->cellidx[atom][1]*g->gridsize[0] +
         g->cellidx[atom][0];

  for(i=g->cellidx[atom][0]-1; i<=g->cellidx[atom][0]+1; i++) {
    if(i == -1) {
      cx=g->gridsize[0]-1;
      sx=-1;
    }
    else if(i == g->gridsize[0]) {
      cx=0;
      sx=1;
    }
    else {
      cx=i;
      sx=0;
    }

    for(j=g->cellidx[atom][1]-1; j<=g->cellidx[atom][1]+1; j++) {
      if(j == -1) {
        cy=g->gridsize[1]-1;
        sy=-1;
      }
      else if(j == g->gridsize[1]) {
       cy=0;
       sy=1;
      }
      else {
        cy=j;
        sy=0;
      }

      for(k=g->cellidx[atom][2]-1; k<=g->cellidx[atom][2]+1; k++) {
        if(k == -1) {
          cz=g->gridsize[2]-1;
          sz=-1;
        }
        else if(k == g->gridsize[2]) {
         cz=0;
         sz=1;
        }
        else {
          cz=k;
          sz=0;
        }

        vx = ((REAL) sx)*g->box[0][0] + ((REAL) sy)*g->box[1][0] + ((REAL) sz)*g->box[2][0];
        vy = ((REAL) sx)*g->box[0][1] + ((REAL) sy)*g->box[1][1] + ((REAL) sz)*g->box[2][1];
        vz = ((REAL) sx)*g->box[0][2] + ((REAL) sy)*g->box[1][2] + ((REAL) sz)*g->box[2][2];

        cidx2 = cz*g->gridsize[0]*g->gridsize[1]+cy*g->gridsize[0]+cx;
        for(l=0; l<g->cells[cidx2].n; l++) {
          num = g->cells[cidx2].idx[l];
          if(num == atom) continue;
          d = sqrt((g->coords[atom][0] - g->coords[num][0] - vx)*(g->coords[atom][0] - g->coords[num][0] - vx) +
                   (g->coords[atom][1] - g->coords[num][1] - vy)*(g->coords[atom][1] - g->coords[num][1] - vy) +
                   (g->coords[atom][2] - g->coords[num][2] - vz)*(g->coords[atom][2] - g->coords[num][2] - vz));
          if(d < g->rcut) {
            g->neighbors[atom].n++;
            if(g->neighbors[atom].n == g->neighbors[atom].maxn) {
              g->neighbors[atom].maxn = (int) (((REAL) g->neighbors[atom].maxn)*1.5);
              g->neighbors[atom].idx = (int *) realloc(g->neighbors[atom].idx, sizeof(int)*g->neighbors[atom].maxn);
              g->neighbors[atom].images = (ivect *) realloc(g->neighbors[atom].images, sizeof(ivect)*g->neighbors[atom].maxn);
              g->neighbors[atom].dists = (REAL *) realloc(g->neighbors[atom].dists, sizeof(REAL)*g->neighbors[atom].maxn);
            }
            g->neighbors[atom].idx[g->neighbors[atom].n-1] = num;
            g->neighbors[atom].images[g->neighbors[atom].n-1][0] = g->cells[cidx2].images[l][0] -
                                                                   g->cells[cidx].images[g->incellidx[atom]][0] + sx;
            g->neighbors[atom].images[g->neighbors[atom].n-1][1] = g->cells[cidx2].images[l][1] -
                                                                   g->cells[cidx].images[g->incellidx[atom]][1] + sy;
            g->neighbors[atom].images[g->neighbors[atom].n-1][2] = g->cells[cidx2].images[l][2] -
                                                                   g->cells[cidx].images[g->incellidx[atom]][2] + sz;
            g->neighbors[atom].dists[g->neighbors[atom].n-1] = d;
          }
        }
      }
    }
  }

  g->nghstatus[atom]=1;
  return g->neighbors[atom];
}

