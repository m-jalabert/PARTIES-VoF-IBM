#ifndef COLLISION_H
#define COLLISION_H

int Collision_above_critical(Particle *p);
Particle_list *Collision_evaluate(Cart3d_bag *data_bag, Debug_trace *dtrace);
void Collision_particle(Particle *p, Particle *p_start, int stage,
		MAC_grid *grid, Parameters *params);
void Collision_wall(Particle *p, MAC_grid *grid, Parameters *params);

#endif
