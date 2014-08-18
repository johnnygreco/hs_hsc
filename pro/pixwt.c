#include <stdio.h>
#include <math.h>
#include "export.h"

/* To compile, link, and use:
     make pixwt    (see associated makefile)

     linkimage,'PIXWT','/gryll/data1/buie/idl/Custom/pixwt.so',1,'pixwt', $
        max_args=5,min_args=5

     After this statement is executed, pixwt will be drawn from the custom
     routine and not the IDL pixwt.pro file.
*/

/* compute the area within an arc of a circle.  The arc is defined by
 * the two points (x,y0) and (x,y1) in the following manner:  The circle
 * is of radius r and is positioned at the origin.  The origin and each
 * individual point define a line which intersect the circle at some
 * point.  The angle between these two points on the circle measured
 * from y0 to y1 defines the sides of a wedge of the circle.  The area
 * returned is the area of this wedge.  If the area is traversed clockwise
 * the the area is negative, otherwise it is positive. */

static double arc (x, y0, y1, r)

double   x;      /* X coordinate of the two points. */
double   y0;     /* Y coordinate of the first point. */
double   y1;     /* Y coordinate of the second point. */
double   r;      /* radius of the circle. */

{
   return( 0.5 * r*r * (atan( y1/x) - atan( y0/x) ) );
}

/* compute the area of a triangle defined by the origin and two points,
 * (x,y0) and (x,y1).  This is a signed area.  If y1 > y0 then the area
 * will be positive, otherwise it will be negative.
 */

static double chord (x, y0, y1)

double   x;      /* X coordinate of the two points. */
double   y0;     /* Y coordinate of the first point. */
double   y1;     /* Y coordinate of the second point. */

{
   return( 0.5 * x * (y1-y0) );
}

/* Compute the area of intersection between a triangle and a circle.
 * The circle is centered at the origin and has a radius of r.  The
 * triangle has verticies at the origin and at (x,y0) and (x,y1).
 * This is a signed area.  The path is traversed from y0 to y1.  If
 * this path takes you clockwise the area will be negative.
 */

static double oneside (x, y0, y1, r)

double   x;      /* X coordinate of the two points. */
double   y0,y1;  /* Y coordinates of the two points. */
double   r;      /* radius of the circle */

{
double   yh;

   if (x == 0) return(0);
   else if ( fabs(x) >=  r ) return( arc(x,y0,y1,r) ); 
   else {
      yh = sqrt(r*r-x*x);
      if (y0 <= -yh) {
         if (y1 <= -yh) return( arc(x,y0,y1,r) );
         else if (y1 <= yh) return( arc(x,y0,-yh,r) + chord(x,-yh,y1) );
         else return( arc(x,y0,-yh,r) + chord(x,-yh,yh)
                                      + arc(x,yh,y1,r) ); }
      else if (y0 < yh) {
         if (y1 < -yh) return( chord(x,y0,-yh) + arc(x,-yh,y1,r) );
         else if (y1 <= yh) return( chord(x,y0,y1) );
         else return( chord(x,y0,yh) + arc(x,yh,y1,r) ); }
      else {
         if (y1<-yh) return( arc(x,y0,yh,r) + chord(x,yh,-yh)
                                            + arc(x,-yh,y1,r) );
         else if (y1 < yh) return( arc(x,y0,yh,r) + chord(x,yh,y1) );
         else return( arc(x,y0,y1,r) ); }
   }
}

/* Compute the area of overlap between a circle and a rectangle. */

double intarea (xc, yc, r, x0, x1, y0, y1)

double   xc,yc;   /* Center of the circle. */
double   r;       /* Radius of the circle. */
double   x0,y0;   /* Corner of the rectangle. */
double   x1,y1;   /* Opposite corner of the rectangle. */

{

x0 -= xc;   y0 -= yc;  /* Shift the objects so that circle is at the orgin. */
x1 -= xc;   y1 -= yc;

return( oneside(x1,y0,y1,r) + oneside(y1,-x1,-x0,r)
        + oneside(-x0,-y1,-y0,r) + oneside(-y0,x0,x1,r) );

}

/*
char *result_var(IDL_VPTR template, int type, IDL_VPTR *res)
*/
/*
 * Allocate a result variable, using the template IDL_VPTR to determine
 * the structure, and type to determine the type.
 * res is set to
 * the new variable, and a pointer to its data area is returned.
 */
/*
{
   char *data;
   IDL_VPTR lres;

   if (template->flags & IDL_V_ARR) {
      data =
         IDL_MakeTempArray(type, template->value.arr->n_dim,
         template->value.arr->dim, IDL_BARR_INI_NOP, res);
      }
   else {
      lres = *res = IDL_Gettmp();
      lres->type = type;
      data = (char *) &(lres->value.c);
      }

   return(data);

}	
*/
/* Compute the fraction of a unit pixel that is interior to a circle.
 * The circle has a radius r and is centered at (xc,yc).  The center of
 * the unit pixel (length of sides = 1) is at (x,y).
 */

IDL_VPTR pixwt(int argc, IDL_VPTR argv[])
{
   IDL_VPTR vxc,vyc,vr,vx,vy;
   IDL_VPTR vp;

   double *xc, *yc, *r, *x, *y, *p;

   int outtype;

   IDL_LONG nx,ny;

/*   fprintf(stderr,"inside pixwt\r\n"); fflush(stderr); */

/*   return(IDL_StrToSTRING("Hello World!")); */

   vxc = argv[0];
   vyc = argv[1];
   vr  = argv[2];
   vx  = argv[3];
   vy  = argv[4];

   IDL_ENSURE_SCALAR(vxc);
   IDL_ENSURE_SCALAR(vyc);
   IDL_ENSURE_SCALAR(vr);

   outtype = IDL_TYP_FLOAT;
   if (vxc->type == IDL_TYP_DOUBLE)      outtype = IDL_TYP_DOUBLE;
   else if (vyc->type == IDL_TYP_DOUBLE) outtype = IDL_TYP_DOUBLE;
   else if (vr->type  == IDL_TYP_DOUBLE) outtype = IDL_TYP_DOUBLE;
   else if (vx->type  == IDL_TYP_DOUBLE) outtype = IDL_TYP_DOUBLE;
   else if (vy->type  == IDL_TYP_DOUBLE) outtype = IDL_TYP_DOUBLE;

   vxc = IDL_BasicTypeConversion(1,&argv[0],IDL_TYP_DOUBLE);
   vyc = IDL_BasicTypeConversion(1,&argv[1],IDL_TYP_DOUBLE);
   vr  = IDL_BasicTypeConversion(1,&argv[2],IDL_TYP_DOUBLE);
   vx  = IDL_BasicTypeConversion(1,&argv[3],IDL_TYP_DOUBLE);
   vy  = IDL_BasicTypeConversion(1,&argv[4],IDL_TYP_DOUBLE);

   IDL_VarGetData(vxc,(long long *)&nx, (char **)&xc, FALSE);
   IDL_VarGetData(vyc,(long long *)&nx, (char **)&yc, FALSE);
   IDL_VarGetData(vr, (long long *)&nx, (char **)&r,  FALSE);
   IDL_VarGetData(vx, (long long *)&nx, (char **)&x,  FALSE);
   IDL_VarGetData(vy, (long long *)&ny, (char **)&y,  FALSE);

   if (nx != ny) {
      IDL_DELTMP(vxc);
      IDL_DELTMP(vyc);
      IDL_DELTMP(vr);
      IDL_DELTMP(vx);
      IDL_DELTMP(vy);
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
          "x and y varibles must be same length");
      }

   if (vx->flags & IDL_V_ARR) {
      p = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, vx->value.arr->n_dim,
              vx->value.arr->dim, IDL_BARR_INI_NOP, &vp);
      }
   else {
      vp = IDL_Gettmp();
      vp->type = IDL_TYP_DOUBLE;
      p = (double *) &(vp->value.c);
      }

/*   fprintf(stderr,"nvals=%ld\r\n",nx); fflush(stderr); */

/*   fprintf(stderr,"just before working loop\r\n"); fflush(stderr); */
   for (; nx--; p++,x++,y++) {
/*      fprintf(stderr," %ld %f %f %f %f %f\r\n",nx,*xc,*yc,*r,*x,*y); */
      *p = intarea( *xc, *yc, *r, *x-0.5, *x+0.5, *y-0.5, *y+0.5 );

      }
/*   fprintf(stderr,"end of pixwt\r\n"); fflush(stdout); */

   IDL_DELTMP(vxc);
   IDL_DELTMP(vyc);
   IDL_DELTMP(vr);
   IDL_DELTMP(vx);
   IDL_DELTMP(vy);
   return(vp);
}

/*
cc -K pic -G -I /opt/local/rsi/idl/external -c pixwt.c
cc -G -o pixwt pixwt.o

linkimage,'pixwt','pixwt',1,'pixwt',max_args=5,min_args=5
*/
