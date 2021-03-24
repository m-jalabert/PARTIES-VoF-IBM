//NOTE: This is an example file included to show how to comment PARTIES functions.
//If you are adding to PARTIES, make sure you have similar comments to maintain
//clarity for new users.

/******************************************************************************/
/*
 * Update and transmit particles between processors
 *
 * Input(s):
 *     p_list - particle linked list structure to communicate
 *			requires p_list->state = LIST_STATE_LOCAL/LIST_STATE_FOREIGN
 *          (for more details see definitions.h)
 *     No arguments can be a null pointer
 *
 * Memory allocation:
 *     Does not allocate memory we have to free later; creates new memory for
 *     shifted particles and frees the old memory
 *
 * Processor Communication:
 *     Particles are sent to neighboring processor if they extend beyond that
 *     neighboring edge.
 *         - Ghost nodes are not taken into account for these considerations.
 *         - Particle coordinates ('X') are updated for periodic boundaries.
 *     Potential speedups:
 *         - Use same send buffer for both lower and upper communications
 *         - Have Np be an input parameter (then Particle_read() should output Np)
 *
 * Assumptions:
 *     - No y-periodicity
 *     - uniform grid
 *
 * NOTE-2021.03.05-Zack Maches: Is "no y-periodicity" still a true assumption,
 * or is this just left over? Y-periodicity seems to be present in the current version.
 */
/******************************************************************************/

void Particle_MPI_update(Particle_list *p_list, Cart3d_bag *data_bag,
                         Debug_trace *dtrace) {
    double T1, T2;
    T1 = MPI_Wtime();

    Parameters *params = data_bag->params;
    MAC_grid   *grid = data_bag->grid;
    Lagrangian *lag = data_bag->lag;

    // Calculate range over which particles can have an influence
    double h = grid->dx_u[1];
    double range = DELTA_FUNC_RADIUS * h;
#if defined LAG_MARKER_FLAG || defined LAG_MARKER_PRIORITY
    range += LAG_FLAG_RANGE * h;
#endif
#if defined LUBRICATION_NORMAL || defined LUBRICATION_TANGENTIAL
    range = max(range, (0.5*params->lub_range*h));
#endif
#ifdef ELECTROSTATIC_REPULSION
    range = max(range, (0.5*params->erm_range*h));
#endif
#ifdef COHESION
    range = max(range, (0.5*params->coh_range*h));
#endif

    // Grid bounds of local subdomain
    double G_xmin = grid->xu[grid->G_Is] + range;
    double G_ymin = grid->yv[grid->G_Js] + range;
    double G_zmin = grid->zw[grid->G_Ks] + range;
    double G_xmax = grid->xu[min(grid->G_Ie, grid->NX - 1)] - range;
    double G_ymax = grid->yv[min(grid->G_Je, grid->NY - 1)] - range;
    double G_zmax = grid->zw[min(grid->G_Ke, grid->NZ - 1)] - range;

#ifdef XPERIODIC
    double Lx = params->Lx;
    int xproccoord = params->xproccoord;
    int NPX = params->NPX;
#endif
#ifdef YPERIODIC
    double Ly = params->Ly;
    int yproccoord = params->yproccoord;
    int NPY = params->NPY;
#endif
#ifdef ZPERIODIC
    double Lz = params->Lz;
    int NPZ = params->NPZ;
    int zproccoord = params->zproccoord;
#endif

    Particle *p;
    int send_Np;  // Number of particles to be sent (to other processor)
    int recv_Np;  // Number of particles received (from other processor)
    Particle *send_part = lag->send_part;  // Particles to be sent
    Particle *recv_part = lag->recv_part;  // Received particles

    Collision *pc;
    int *send_Nc = lag->send_Nc;  // Number of collisions per particle to be sent
    int *recv_Nc = lag->recv_Nc;  // Number of collisions per particle received
    Collision *send_coll = lag->send_coll;  // Collision to be sent
    Collision *recv_coll = lag->recv_coll;  // Received collisions

    MPI_Status status;
    int src;   // Index ofsource procesor (where reveiving from)
    int dest;  // Index of destination processor (where sending to)
    int tag = 1;

    // Proactively set state to reflect state of p_list after function execution
    if (p_list->state == LIST_STATE_LOCAL) {
        p_list->state = LIST_STATE_BOTH;
    } else if (p_list->state == LIST_STATE_FOREIGN) {
        p_list->state = LIST_STATE_EDGE;
    } else {
        // Any other initial state will produce a useless set of particles
        char message[100], name0[50], name1[50], name2[50];
        Display_print_list_name(name0, p_list->state);
        Display_print_list_name(name1, LIST_STATE_LOCAL);
        Display_print_list_name(name2, LIST_STATE_FOREIGN);
        sprintf(message, "Invalid input linked list state '%s'\n"
                "Should be '%s' or '%s'", name0, name1, name2);
        Display_throw_error(message, params, DTRACE("Display_throw_error"));
    }

#ifdef XPERIODIC
    // Ensure we don't communicate foreign particles to ourself
    if (params->NPX != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
        // Lower x communication ie. communicating with processor to the 'left'
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            if (p->X[0] - R_EFF < G_xmin) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef XPERIODIC
                if (xproccoord == 0) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[0] += Lx;
                    send_part[send_Np].X_old[0] += Lx;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
        dest = params->npxminus;
        src  = params->npxplus;
        SENDRECV_PARTICLES_COLLISIONS();

        // Upper x communication ie. communicating with processor to the 'right'
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            // Send nearby particles
            if (p->X[0] + R_EFF > G_xmax) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef XPERIODIC
                if (xproccoord == NPX - 1) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[0] -= Lx;
                    send_part[send_Np].X_old[0] -= Lx;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
        // Add particles received during lower x communication.
        // This placement (after filling send buffer) ensures that we do not
        // send particles we just received back to the host.
        Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

        dest = params->npxplus;
        src  = params->npxminus;
        SENDRECV_PARTICLES_COLLISIONS();
        // Add particles received during upper x communication.
        Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef XPERIODIC
    }
#endif

#ifdef YPERIODIC
    // Ensure we don't communicate foreign particles to ourself
    if (params->NPY != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
        // Lower y communication
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            if (p->X[1] - R_EFF < G_ymin) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef YPERIODIC
                if (yproccoord == 0) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[1] += Ly;
                    send_part[send_Np].X_old[1] += Ly;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
        dest = params->npyminus;
        src  = params->npyplus;
        SENDRECV_PARTICLES_COLLISIONS();

        // Upper y communication
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            if (p->X[1] + R_EFF > G_ymax) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef YPERIODIC
                if (yproccoord == NPY - 1) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[1] -= Ly;
                    send_part[send_Np].X_old[1] -= Ly;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
    // Add particles received during lower y communication.
    // This placement (after filling send buffer) ensures that we do not
    // send particles we just received back to the host.
    Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

    dest = params->npyplus;
    src  = params->npyminus;
    SENDRECV_PARTICLES_COLLISIONS();
    // Add particles received during upper y communication.
    Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef YPERIODIC
    }
#endif

#ifdef ZPERIODIC
    // Ensure we don't communicate foreign particles to ourself
    if (params->NPZ != 1 || p_list->state != LIST_STATE_EDGE) {
#endif
        // Lower z communication
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            if (p->X[2] - R_EFF < G_zmin) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef ZPERIODIC
                if (zproccoord == 0) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[2] += Lz;
                    send_part[send_Np].X_old[2] += Lz;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
        dest = params->npzminus;
        src  = params->npzplus;
        SENDRECV_PARTICLES_COLLISIONS();

        // Upper z communication
        send_Np = 0;
        send_Nc[0] = 0;
        p = p_list->start;
        while (p != NULL) {
            // Send nearby particles
            if (p->X[2] + R_EFF > G_zmax) {
                ADD_TO_BUFFER(send_part, send_Np, *p, Particle);
                NULLIFY_PARTICLE_PTRS(send_part[send_Np]);
#ifdef ZPERIODIC
                if (zproccoord == NPZ - 1) {
                    // Update particle position to reflect periodicity
                    send_part[send_Np].X[2] -= Lz;
                    send_part[send_Np].X_old[2] -= Lz;
                }
#endif
                FILL_COLLISION_SEND_BUFFER();
                send_Np++;
            }
            p = p->next;
        }
        // Add particles received during lower x communication.
        // This placement (after filling send buffer) ensures that we do not
        // send particles we just received back to the host.
        Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

        dest = params->npzplus;
        src  = params->npzminus;
        SENDRECV_PARTICLES_COLLISIONS();
        // Add particles received during upper z communication.
        Particle_list_add_array(p_list, recv_part, recv_Np, recv_coll, recv_Nc, grid);

#ifdef ZPERIODIC
    }
#endif
    T2 = MPI_Wtime();
    data_bag->timer->Wtime_particle_comm += T2 - T1;
}
