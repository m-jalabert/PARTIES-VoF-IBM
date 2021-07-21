# y-Periodicity

It seems that the single/individual periodicity in y-direction is not implemented into the code. For example, simulating a box with no-slip walls in x- and z-directions and periodic boundaries in y-direction results in an error within the first iterations.

<pre>
***************************************************************************
!!! Error called by Display_throw_error (Temporal_int.c:1052)           !!!
!!!              at Temporal_int_all_the_equations (Temporal_int.c:473) !!!
!!!              at Temporal_int_rk3 (Cart3d.c:446)                     !!!
!!! Poisson solver did not converge                                     !!!
***************************************************************************
</pre>

In contrast, the same setup with periodic boundary conditions in all directions (TRIPLE_Periodic) runs successfully. The analysis of the code regarding the used expressions of "XPERIODIC", "ZPERIODIC" and "YPERIODIC" explains the above statement that the y-periodicity is not implemented yet. 

---

The search for "XPERIODIC" leads to 186 results in 30 files (within PARTIES).<br>
The search for "ZPERIODIC" leads to 159 results in 30 files (within PARTIES).<br>
The search for "YPERIODIC" leads to 57 results in 17 files (within PARTIES).<br>

---
	
It must be emphasized that this only applies to the single y-periodicity, as the combined one included in "TRIPLE_PERIODIC" works fine.

