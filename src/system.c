#include "system.h"
#include <stdlib.h>

sys* system_create(int natoms, int nbfs){
  sys *ret_sys = malloc(sizeof (sys));
  if (ret_sys == NULL)
    return NULL;

  ret_sys->natoms = natoms;
  ret_sys->nbfs = nbfs;
  ret_sys->basisfunctions = malloc(sizeof (bfn) * nbfs);
  if (ret_sys->basisfunctions == NULL){
    free(ret_sys);
    return NULL;
  }
  ret_sys->S = malloc(sizeof (double) * nbfs*nbfs);
  if (ret_sys->S == NULL){
    free(ret_sys);
    return NULL;
  }
  ret_sys->T = malloc(sizeof (double) *  nbfs*nbfs);
  if (ret_sys->T == NULL){
    free(ret_sys);
    return NULL;
  }
  ret_sys->V = malloc(sizeof (double) *  nbfs*nbfs);
  if (ret_sys->V == NULL){
    free(ret_sys);
    return NULL;
  }
  ret_sys->ERI = malloc(sizeof (double) *  nbfs*nbfs*nbfs*nbfs);
  if (ret_sys->ERI == NULL){
    free(ret_sys);
    return NULL;
  }
  return(ret_sys);
}

void system_destroy(sys* in_system){
  if (in_system != NULL) {
    free(in_system->basisfunctions);
    free(in_system->S);
    free(in_system->T);
    free(in_system->V);
    free(in_system->ERI);
    free(in_system);
  }
}

void print_system_info(sys* system){

}

void system_from_file(FILE* infile, char* basis){

}
