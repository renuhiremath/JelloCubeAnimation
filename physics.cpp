/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

void calculateHookesLawForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2, double side, point &force)
{
  point hooksForce, dampingForce, dv;
  double length;

  point L, normL;
  double magL = 0.0;
  double temp = 0.0;

  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  pCPY(force, L);
  pCPY(force, normL);
  pCPY(force, hooksForce);
  pCPY(force, dampingForce);

  pDIFFERENCE(jello->p[i1][j1][k1], jello->p[i2][j2][k2], L);
  pDIFFERENCE(jello->v[i1][j1][k1], jello->v[i2][j2][k2], dv);
  magL = pMAGNITUDE(L);

  if (magL == 0)
  {
    return;
  }

  pCPY(L, normL);
  pNORMALIZE(normL);
  temp = pDOTPRODUCT(dv, L);

  hooksForce.x = (-1 * jello->kElastic * (magL - side) * normL.x);
  hooksForce.y = (-1 * jello->kElastic * (magL - side) * normL.y);
  hooksForce.z = (-1 * jello->kElastic * (magL - side) * normL.z);

  dampingForce.x = (-1 * jello->dElastic * (temp/magL) * normL.x);
  dampingForce.y = (-1 * jello->dElastic * (temp/magL) * normL.y);
  dampingForce.z = (-1 * jello->dElastic * (temp/magL) * normL.z);

  pSUM(dampingForce, hooksForce, force);

}

point calculateStructuralForce(struct world * jello, int i, int j, int k, point &force)
{
  double side = 1.0/7.0;

  point temp;
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;

  if (i-1 >= 0)
  {
    calculateHookesLawForce(jello, i, j, k, i-1, j, k, side, temp);
    pSUM(force, temp, force);
  }
  if (i+1 < 8)
  {
    calculateHookesLawForce(jello, i, j, k, i+1, j, k, side, temp);
    pSUM(force, temp, force);
  }

  if (j-1 >= 0)
  {
    calculateHookesLawForce(jello, i, j, k, i, j-1, k, side, temp);
    pSUM(force, temp, force);
  }
  if (j+1 < 8)
  {
    calculateHookesLawForce(jello, i, j, k, i, j+1, k, side, temp);
    pSUM(force, temp, force);
  }

  if (k-1 >= 0)
  {
    calculateHookesLawForce(jello, i, j, k, i, j, k-1, side, temp);
    pSUM(force, temp, force);
  }
  if (k+1 < 8)
  {
    calculateHookesLawForce(jello, i, j, k, i, j, k+1, side, temp);
    pSUM(force, temp, force);
  }

  return force;
}

point calculateShearForce(struct world * jello, int i, int j, int k, point &force)
{
  double side = 1.0/7.0;
  double diagonal1 = sqrt(2.0) * side;
  double diagonal2 = sqrt(3.0) * side;

  point temp;
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;

  //diagonal 3 type - 8
  {
    if (i+1<8 && j+1<8 && k+1<8)
    {
       calculateHookesLawForce(jello, i, j, k, i+1, j+1, k+1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && j-1>=0 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j-1, k+1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && j+1<8 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j+1, k-1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && j-1>=0 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j-1, k-1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j-1>=0 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j-1, k-1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j-1>=0 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j-1, k+1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j+1<8 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j+1, k-1, diagonal2, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j+1<8 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j+1, k+1, diagonal2, temp);
      pSUM(force, temp, force);
    }
  }

  //diagonal 2 type - 12
  {
    if (i-1>=0 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j, k-1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j, k+1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j, k-1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j, k+1, diagonal1, temp);
      pSUM(force, temp, force);
    }

    if (j+1<8 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i, j+1, k+1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (j+1<8&& k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i, j+1, k-1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (j-1>=0 && k+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i, j-1, k+1, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (j-1>=0 && k-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i, j-1, k-1, diagonal1, temp);
      pSUM(force, temp, force);
    }

    if (i+1<8 && j+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j+1, k, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i+1<8 && j-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i+1, j-1, k, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j+1<8)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j+1, k, diagonal1, temp);
      pSUM(force, temp, force);
    }
    if (i-1>=0 && j-1>=0)
    {
      calculateHookesLawForce(jello, i, j, k, i-1, j-1, k, diagonal1, temp);
      pSUM(force, temp, force);
    }
  }

  return force;
}

point calculateBendForce(struct world * jello, int i, int j, int k, point &force)
{
  double side = 2.0/7.0;

  point temp;
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;

  if (i-2>=0)
  {
    calculateHookesLawForce(jello, i, j, k, i-2, j, k, side, temp);
    pSUM(force, temp, force);
  }
  if (i+2<8)
  {
    calculateHookesLawForce(jello, i, j, k, i+2, j, k, side, temp);
    pSUM(force, temp, force);
  }
  if (j+2<8)
  {
    calculateHookesLawForce(jello, i, j, k, i, j+2, k, side, temp);
    pSUM(force, temp, force);
  }
  if (j-2>=0)
  {
    calculateHookesLawForce(jello, i, j, k, i, j-2, k, side, temp);
    pSUM(force, temp, force);
  }
  if (k-2>=0)
  {
    calculateHookesLawForce(jello, i, j, k, i, j, k-2, side, temp);
    pSUM(force, temp, force);
  }
  if (k+2<8)
  {
    calculateHookesLawForce(jello, i, j, k, i, j, k+2, side, temp);
    pSUM(force, temp, force);
  }

  return force;
}

point calculateExternalForce(struct world * jello, int i, int j, int k, point &force)
{
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  point p;
  p.x = jello->p[i][j][k].x + 2;
  p.y = jello->p[i][j][k].y + 2;
  p.z = jello->p[i][j][k].z + 2;

  int r = jello->resolution;
  double xl, xh, yl, yh, zl, zh;
  double x, y, z;

  x = ((r-1) * p.x) / 4.0;
  y = ((r-1) * p.y) / 4.0;
  z = ((r-1) * p.z) / 4.0;

  xl = floor(x);
  yl = floor(y);
  zl = floor(z);

  if (xl < (r-1))
  {
    if (xl < 0)
      xl = 0;
    xh = xl + 1;
  }
  else
  {
    xl = r - 1;
    xh = r - 1;
  }

  if (yl < (r-1))
  {
    if (yl < 0)
      yl = 0;
    yh = yl + 1;
  }
  else
  {
    yl = r - 1;
    yh = r - 1;
  }

  if (zl < (r-1))
  {
    if (zl < 0)
      zl = 0;
    zh = zl + 1;
  }
  else
  {
    zl = r - 1;
    zh = r - 1;
  }

  double wx = (x - xl);
  double wy = (y - yl);
  double wz = (z - zl);

  point forceField[8];
  forceField[0] = jello->forceField[int((xl)*r*r + (yl)*r + (zl))];
  forceField[1] = jello->forceField[int((xl)*r*r + (yl)*r + (zl+1))];
  forceField[2] = jello->forceField[int((xl)*r*r + (yl+1)*r + (zl))];
  forceField[3] = jello->forceField[int((xl)*r*r + (yl+1)*r + (zl+1))];
  forceField[4] = jello->forceField[int((xl+1)*r*r + (yl)*r + (zl))];
  forceField[5] = jello->forceField[int((xl+1)*r*r + (yl)*r + (zl+1))];
  forceField[6] = jello->forceField[int((xl+1)*r*r + (yl+1)*r + (zl))];
  forceField[7] = jello->forceField[int((xl+1)*r*r + (yl+1)*r + (zl+1))];

  force.x = ((1 - wx)*(1 - wy)*(1 - wz))*forceField[0].x +
            ((1 - wx)*(1 - wy)*(wz))*jello->forceField[1].x +
            ((1 - wx)*(wy)*(1 - wz))*jello->forceField[2].x +
            ((1 - wx)*(wy)*(wz))*jello->forceField[3].x +
            ((wx)*(1 - wy)*(1 - wz))*jello->forceField[4].x +
            ((wx)*(1 - wy)*(wz))*jello->forceField[5].x +
            ((wx)*(wy)*(1 - wz))*jello->forceField[6].x +
            ((wx)*(wy)*(wz))*jello->forceField[7].x;

  force.y = ((1 - wx)*(1 - wy)*(1 - wz))*forceField[0].y +
            ((1 - wx)*(1 - wy)*(wz))*jello->forceField[1].y +
            ((1 - wx)*(wy)*(1 - wz))*jello->forceField[2].y +
            ((1 - wx)*(wy)*(wz))*jello->forceField[3].y +
            ((wx)*(1 - wy)*(1 - wz))*jello->forceField[4].y +
            ((wx)*(1 - wy)*(wz))*jello->forceField[5].y +
            ((wx)*(wy)*(1 - wz))*jello->forceField[6].y +
            ((wx)*(wy)*(wz))*jello->forceField[7].y;

  force.z = ((1 - wx)*(1 - wy)*(1 - wz))*forceField[0].z +
            ((1 - wx)*(1 - wy)*(wz))*jello->forceField[1].z +
            ((1 - wx)*(wy)*(1 - wz))*jello->forceField[2].z +
            ((1 - wx)*(wy)*(wz))*jello->forceField[3].z +
            ((wx)*(1 - wy)*(1 - wz))*jello->forceField[4].z +
            ((wx)*(1 - wy)*(wz))*jello->forceField[5].z +
            ((wx)*(wy)*(1 - wz))*jello->forceField[6].z +
            ((wx)*(wy)*(wz))*jello->forceField[7].z;


  //printf("%lf\t%lf\t%lf\n", force.x, force.y, force.z);
}

void calculateCollisionForceForSingleAxis(struct world * jello, point L, point dv, point &force)
{
  point hooksForce, dampingForce;
  double length;

  point normL;
  double magL = 0.0;
  double temp = 0.0;

  magL = pMAGNITUDE(L);

  if (magL == 0)
    return;

  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  pCPY(force, normL);
  pCPY(force, hooksForce);
  pCPY(force, dampingForce);

  pCPY(L, normL);
  pNORMALIZE(normL);
  temp = pDOTPRODUCT(dv, L);

  hooksForce.x = (-1 * jello->kCollision * magL * normL.x);
  hooksForce.y = (-1 * jello->kCollision * magL * normL.y);
  hooksForce.z = (-1 * jello->kCollision * magL * normL.z);

  dampingForce.x = (1 * jello->dCollision * (temp/magL) * normL.x);
  dampingForce.y = (1 * jello->dCollision * (temp/magL) * normL.y);
  dampingForce.z = (1 * jello->dCollision * (temp/magL) * normL.z);

  pSUM(dampingForce, hooksForce, force);
  //printf("%lf\t%lf\t%lf\n", force.x, force.y, force.z);

}

void calculateCollisionForce(struct world * jello, int i, int j, int k, point &force)
{
  point temp, dx, wall;
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  pCPY(force, temp);

  if (jello->p[i][j][k].x < -2)
  {
    dx.x = jello->p[i][j][k].x - (-2);
    dx.y = 0.0;
    dx.z = 0.0;

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }
  else if (jello->p[i][j][k].x > 2)
  {
    dx.x = jello->p[i][j][k].x - (2);
    dx.y = 0.0;
    dx.z = 0.0;

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }

  if (jello->p[i][j][k].y < -2)
  {
    dx.x = 0.0;
    dx.y = jello->p[i][j][k].y - (-2);
    dx.z = 0.0;

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }
  else if (jello->p[i][j][k].y > 2)
  {
    dx.x = 0.0;
    dx.y = jello->p[i][j][k].y - (2);
    dx.z = 0.0;

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }

  if (jello->p[i][j][k].z < -2)
  {
    dx.x = 0.0;
    dx.y = 0.0;
    dx.z = jello->p[i][j][k].z - (-2);

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }
  else if (jello->p[i][j][k].z > 2)
  {
    dx.x = 0.0;
    dx.y = 0.0;
    dx.z = jello->p[i][j][k].z - (2);

    calculateCollisionForceForSingleAxis(jello, dx, jello->v[i][j][k], temp);
    pSUM(temp, force, force);
  }
}

void calculateCollisionForceWithPlane(struct world * jello, int i, int j, int k, point &force)
{
  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;

  //calculateCollisionForceForSingleAxis


}

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. It returns the acceleration for each of the 512 points. It also adds any effects of an external force field, if such a field is present of course. In general, the acceleration will of course be different for each of the 512 simulation points. To compute the acceleration, the function must take into account the forces due to

structural, shear and bend springs,
external forces (force field, if any), and
bouncing off the walls.*/
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  point hooksForce, externalForce, collisionForce, finalForce;
  point temp;

  for (int i=0; i<8; i++)
  {
    for (int j=0; j<8; j++)
    {
      for (int k=0; k<8; k++)
      {
        temp.x = 0.0;
        temp.y = 0.0;
        temp.z = 0.0;

        pCPY(temp, hooksForce);
        pCPY(temp, externalForce);
        pCPY(temp, collisionForce);
        pCPY(temp, finalForce);

        //calculating hooke's force
        //structural force
        calculateStructuralForce(jello, i, j, k, temp);
        pSUM(temp, hooksForce, hooksForce);

        //shear force
        temp.x = 0.0;
        temp.y = 0.0;
        temp.z = 0.0;
        calculateShearForce(jello, i, j, k, temp);
        pSUM(temp, hooksForce, hooksForce);

                // printf("%lf\t", hooksForce.x);
                // printf("%lf\t", hooksForce.y);
                // printf("%lf\t", hooksForce.z);
        //bind force
        temp.x = 0.0;
        temp.y = 0.0;
        temp.z = 0.0;
        calculateBendForce(jello, i, j, k, temp);
        pSUM(temp, hooksForce, hooksForce);

        //calculating collision force
        //printf("%lf\t%lf\t%lf\n", jello->p[i][j][k].x, jello->p[i][j][k].y, jello->p[i][j][k].z);
        if (jello->p[i][j][k].x < -2 || jello->p[i][j][k].x > 2 ||
            jello->p[i][j][k].y < -2 || jello->p[i][j][k].y > 2 ||
            jello->p[i][j][k].z < -2 || jello->p[i][j][k].z > 2 )
        {
          calculateCollisionForce(jello, i, j, k, collisionForce);
        }

        //calculating external force
        calculateExternalForce(jello, i, j, k, externalForce);

        //calculating the final force
        pSUM(hooksForce, externalForce, finalForce);
        pSUM(finalForce, collisionForce, finalForce);

        //calculating acceleration for the mass point
        a[i][j][k].x = finalForce.x / jello->mass;
        a[i][j][k].y = finalForce.y / jello->mass;
        a[i][j][k].z = finalForce.z / jello->mass;
      }
    }
  }

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8],
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;
}
