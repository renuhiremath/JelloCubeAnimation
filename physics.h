/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void calculateHookesLawForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2, double side, point &force);

point calculateStructuralForce(struct world * jello, int i, int j, int k, point &force);

point calculateShearForce(struct world * jello, int i, int j, int k, point &force);

point calculateBendForce(struct world * jello, int i, int j, int k, point &force);

point calculateExternalForce(struct world * jello, int i, int j, int k, point &force);

void calculateCollisionForceForSingleAxis(struct world * jello, point L, point dv, point &force);

void calculateCollisionForce(struct world * jello, int i, int j, int k, point &force);

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif
