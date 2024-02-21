#ifndef MYTIME_H
#define MYTIME_H

#include <ctime>

//si on n'est pas sous windows
#ifndef _WIN32
  #include <sys/resource.h>
  #include <sys/time.h>
  #define my_time double

#else
  #define my_time double //time_t
#endif


//renvoie une mesure de temps en seconde a l'instant donne
inline my_time give_time() 
{
  #ifndef _WIN32

    struct rusage ru;
    struct timeval tim;
    getrusage(RUSAGE_SELF, &ru);
   

    tim = ru.ru_utime;
    //temps systeme en secondes
    double stime = (double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
    
    return stime;

  #else
   //return time(NULL);
   return (double)(clock()) / CLOCKS_PER_SEC; //temps d'attente
  #endif
}

#endif
