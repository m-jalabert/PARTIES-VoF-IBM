#ifndef TIMER_H
 #define TIMER_H
Timer *Timer_create(MAC_grid *grid, Parameters *params);
void Timer_destroy(Timer *timer, MAC_grid *grid, Parameters *params);
#endif
