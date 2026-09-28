    /* Test-only manufactured translation. This include is injected into the
     * isolated gate build, never the production solver. Shift the old sphere
     * mask one x cell, seed a smooth mixed phase, then rebuild the real mask. */
    {
        Parameters *p=db->params;
        Array_copy_withghost(vof->C_S,vof->ch_rhs_nm1,grid,p);
        double local[2]={0.0,0.0},global[2];
        for(int k=grid->G_Ks;k<min(grid->G_Ke,grid->NZ-1);++k)
        for(int j=grid->G_Js;j<min(grid->G_Je,grid->NY-1);++j)
        for(int i=grid->G_Is;i<min(grid->G_Ie,grid->NX-1);++i) {
            const double cs=vof->ch_rhs_nm1[k][j][i+1];
            vof->C_S[k][j][i]=cs;
            const double cap=cs<DIFFUSE_SOLID_MASS_CUTOFF?1.0-cs:0.0;
            const double liquid=.65+.1*cos(2.0*acos(-1.0)*grid->xc[i]/4.0);
            vof->C_L[k][j][i]=cap*liquid;
            local[0]+=cap*(1.0-liquid)*diffuse_physical_cell_measure(grid,p,i,j,k);
        }
        VOF_DIFFUSE_set_boundary_values(vof->C_S,db);
        VOF_DIFFUSE_set_boundary_values(vof->C_L,db);
        VOF_DIFFUSE_move_ice_mask(db);
        for(int k=grid->G_Ks;k<min(grid->G_Ke,grid->NZ-1);++k)
        for(int j=grid->G_Js;j<min(grid->G_Je,grid->NY-1);++j)
        for(int i=grid->G_Is;i<min(grid->G_Ie,grid->NX-1);++i)
            if(vof->C_S[k][j][i]<DIFFUSE_SOLID_MASS_CUTOFF)
                local[1]+=(1.0-vof->C_S[k][j][i]-vof->C_L[k][j][i])*diffuse_physical_cell_measure(grid,p,i,j,k);
        MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,PCW);
        char msg[180];
        snprintf(msg,sizeof(msg),"REMAP_FIXTURE before=%.17g after=%.17g delta=%.17g\n",global[0],global[1],global[1]-global[0]);
        Display_progress(p,msg);
    }
