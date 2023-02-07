#include<iostream>
#include<math.h>

typedef double real;
//Pendulum eqn is mx'' + cx' + kx=0;
real SMD_func(real m,real k,real c,real x,real xdot,real &t); //SMD_func->Spring-mass damper function

real SMD_analy(real m,real k,real c,real &t);                                      //SMD_analy->Analytical solution for fixed mass 

real SMD_func(real m,real k,real c,real x,real xdot,real &t)
{
    real xddot=-(k/m)*x -(c/m)*xdot;
    return xddot;
}

real SMD_analy(real m,real k,real c,real &t)
{   //The initial conditions are at t=0,x=1 and xdot=0. Based on these ICs, we get the values of constants
    real x;
    if (pow((c/(2*m)),2) < (k/m))        //The roots are imaginary
    {
        real A=1;
        real omega=sqrt( (k/m) - pow((c/(2*m)),2));
        real B=k/(m*omega);
        real x=exp(-c*t/(2*m))*(A*cos(omega*t) + B*sin(omega*t));
        return x;
    }

    else if (pow((c/(2*m)),2) == (k/m))  //The roots are real and equal
    {
        real A=1;
        real B=c/(2*m);
        real x = (A + B*t)*exp(-c*t/(2*m));
        return x;
    }

    else if (pow((c/(2*m)),2) > (k/m))   //The roots are distinct and real
    {
        real alpha_term= -(c/(2*m)) + sqrt(-(k/m) + pow((c/(2*m)),2));
        real beta_term= -(c/(2*m)) - sqrt(-(k/m) + pow((c/(2*m)),2));
        real A= - beta_term/(alpha_term - beta_term);
        real B= alpha_term/(alpha_term - beta_term);
        real x=A*exp(alpha_term*t) + B*exp(beta_term*t);
        return x;
    }
}
