#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3580642894513962333) {
   out_3580642894513962333[0] = delta_x[0] + nom_x[0];
   out_3580642894513962333[1] = delta_x[1] + nom_x[1];
   out_3580642894513962333[2] = delta_x[2] + nom_x[2];
   out_3580642894513962333[3] = delta_x[3] + nom_x[3];
   out_3580642894513962333[4] = delta_x[4] + nom_x[4];
   out_3580642894513962333[5] = delta_x[5] + nom_x[5];
   out_3580642894513962333[6] = delta_x[6] + nom_x[6];
   out_3580642894513962333[7] = delta_x[7] + nom_x[7];
   out_3580642894513962333[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2148588102071341310) {
   out_2148588102071341310[0] = -nom_x[0] + true_x[0];
   out_2148588102071341310[1] = -nom_x[1] + true_x[1];
   out_2148588102071341310[2] = -nom_x[2] + true_x[2];
   out_2148588102071341310[3] = -nom_x[3] + true_x[3];
   out_2148588102071341310[4] = -nom_x[4] + true_x[4];
   out_2148588102071341310[5] = -nom_x[5] + true_x[5];
   out_2148588102071341310[6] = -nom_x[6] + true_x[6];
   out_2148588102071341310[7] = -nom_x[7] + true_x[7];
   out_2148588102071341310[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6563764178898461836) {
   out_6563764178898461836[0] = 1.0;
   out_6563764178898461836[1] = 0;
   out_6563764178898461836[2] = 0;
   out_6563764178898461836[3] = 0;
   out_6563764178898461836[4] = 0;
   out_6563764178898461836[5] = 0;
   out_6563764178898461836[6] = 0;
   out_6563764178898461836[7] = 0;
   out_6563764178898461836[8] = 0;
   out_6563764178898461836[9] = 0;
   out_6563764178898461836[10] = 1.0;
   out_6563764178898461836[11] = 0;
   out_6563764178898461836[12] = 0;
   out_6563764178898461836[13] = 0;
   out_6563764178898461836[14] = 0;
   out_6563764178898461836[15] = 0;
   out_6563764178898461836[16] = 0;
   out_6563764178898461836[17] = 0;
   out_6563764178898461836[18] = 0;
   out_6563764178898461836[19] = 0;
   out_6563764178898461836[20] = 1.0;
   out_6563764178898461836[21] = 0;
   out_6563764178898461836[22] = 0;
   out_6563764178898461836[23] = 0;
   out_6563764178898461836[24] = 0;
   out_6563764178898461836[25] = 0;
   out_6563764178898461836[26] = 0;
   out_6563764178898461836[27] = 0;
   out_6563764178898461836[28] = 0;
   out_6563764178898461836[29] = 0;
   out_6563764178898461836[30] = 1.0;
   out_6563764178898461836[31] = 0;
   out_6563764178898461836[32] = 0;
   out_6563764178898461836[33] = 0;
   out_6563764178898461836[34] = 0;
   out_6563764178898461836[35] = 0;
   out_6563764178898461836[36] = 0;
   out_6563764178898461836[37] = 0;
   out_6563764178898461836[38] = 0;
   out_6563764178898461836[39] = 0;
   out_6563764178898461836[40] = 1.0;
   out_6563764178898461836[41] = 0;
   out_6563764178898461836[42] = 0;
   out_6563764178898461836[43] = 0;
   out_6563764178898461836[44] = 0;
   out_6563764178898461836[45] = 0;
   out_6563764178898461836[46] = 0;
   out_6563764178898461836[47] = 0;
   out_6563764178898461836[48] = 0;
   out_6563764178898461836[49] = 0;
   out_6563764178898461836[50] = 1.0;
   out_6563764178898461836[51] = 0;
   out_6563764178898461836[52] = 0;
   out_6563764178898461836[53] = 0;
   out_6563764178898461836[54] = 0;
   out_6563764178898461836[55] = 0;
   out_6563764178898461836[56] = 0;
   out_6563764178898461836[57] = 0;
   out_6563764178898461836[58] = 0;
   out_6563764178898461836[59] = 0;
   out_6563764178898461836[60] = 1.0;
   out_6563764178898461836[61] = 0;
   out_6563764178898461836[62] = 0;
   out_6563764178898461836[63] = 0;
   out_6563764178898461836[64] = 0;
   out_6563764178898461836[65] = 0;
   out_6563764178898461836[66] = 0;
   out_6563764178898461836[67] = 0;
   out_6563764178898461836[68] = 0;
   out_6563764178898461836[69] = 0;
   out_6563764178898461836[70] = 1.0;
   out_6563764178898461836[71] = 0;
   out_6563764178898461836[72] = 0;
   out_6563764178898461836[73] = 0;
   out_6563764178898461836[74] = 0;
   out_6563764178898461836[75] = 0;
   out_6563764178898461836[76] = 0;
   out_6563764178898461836[77] = 0;
   out_6563764178898461836[78] = 0;
   out_6563764178898461836[79] = 0;
   out_6563764178898461836[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_9217663828102947708) {
   out_9217663828102947708[0] = state[0];
   out_9217663828102947708[1] = state[1];
   out_9217663828102947708[2] = state[2];
   out_9217663828102947708[3] = state[3];
   out_9217663828102947708[4] = state[4];
   out_9217663828102947708[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_9217663828102947708[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_9217663828102947708[7] = state[7];
   out_9217663828102947708[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4765097794067878070) {
   out_4765097794067878070[0] = 1;
   out_4765097794067878070[1] = 0;
   out_4765097794067878070[2] = 0;
   out_4765097794067878070[3] = 0;
   out_4765097794067878070[4] = 0;
   out_4765097794067878070[5] = 0;
   out_4765097794067878070[6] = 0;
   out_4765097794067878070[7] = 0;
   out_4765097794067878070[8] = 0;
   out_4765097794067878070[9] = 0;
   out_4765097794067878070[10] = 1;
   out_4765097794067878070[11] = 0;
   out_4765097794067878070[12] = 0;
   out_4765097794067878070[13] = 0;
   out_4765097794067878070[14] = 0;
   out_4765097794067878070[15] = 0;
   out_4765097794067878070[16] = 0;
   out_4765097794067878070[17] = 0;
   out_4765097794067878070[18] = 0;
   out_4765097794067878070[19] = 0;
   out_4765097794067878070[20] = 1;
   out_4765097794067878070[21] = 0;
   out_4765097794067878070[22] = 0;
   out_4765097794067878070[23] = 0;
   out_4765097794067878070[24] = 0;
   out_4765097794067878070[25] = 0;
   out_4765097794067878070[26] = 0;
   out_4765097794067878070[27] = 0;
   out_4765097794067878070[28] = 0;
   out_4765097794067878070[29] = 0;
   out_4765097794067878070[30] = 1;
   out_4765097794067878070[31] = 0;
   out_4765097794067878070[32] = 0;
   out_4765097794067878070[33] = 0;
   out_4765097794067878070[34] = 0;
   out_4765097794067878070[35] = 0;
   out_4765097794067878070[36] = 0;
   out_4765097794067878070[37] = 0;
   out_4765097794067878070[38] = 0;
   out_4765097794067878070[39] = 0;
   out_4765097794067878070[40] = 1;
   out_4765097794067878070[41] = 0;
   out_4765097794067878070[42] = 0;
   out_4765097794067878070[43] = 0;
   out_4765097794067878070[44] = 0;
   out_4765097794067878070[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4765097794067878070[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4765097794067878070[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4765097794067878070[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4765097794067878070[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4765097794067878070[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4765097794067878070[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4765097794067878070[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4765097794067878070[53] = -9.8000000000000007*dt;
   out_4765097794067878070[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4765097794067878070[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4765097794067878070[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4765097794067878070[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4765097794067878070[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4765097794067878070[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4765097794067878070[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4765097794067878070[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4765097794067878070[62] = 0;
   out_4765097794067878070[63] = 0;
   out_4765097794067878070[64] = 0;
   out_4765097794067878070[65] = 0;
   out_4765097794067878070[66] = 0;
   out_4765097794067878070[67] = 0;
   out_4765097794067878070[68] = 0;
   out_4765097794067878070[69] = 0;
   out_4765097794067878070[70] = 1;
   out_4765097794067878070[71] = 0;
   out_4765097794067878070[72] = 0;
   out_4765097794067878070[73] = 0;
   out_4765097794067878070[74] = 0;
   out_4765097794067878070[75] = 0;
   out_4765097794067878070[76] = 0;
   out_4765097794067878070[77] = 0;
   out_4765097794067878070[78] = 0;
   out_4765097794067878070[79] = 0;
   out_4765097794067878070[80] = 1;
}
void h_25(double *state, double *unused, double *out_4673975711821749475) {
   out_4673975711821749475[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1459547032221503854) {
   out_1459547032221503854[0] = 0;
   out_1459547032221503854[1] = 0;
   out_1459547032221503854[2] = 0;
   out_1459547032221503854[3] = 0;
   out_1459547032221503854[4] = 0;
   out_1459547032221503854[5] = 0;
   out_1459547032221503854[6] = 1;
   out_1459547032221503854[7] = 0;
   out_1459547032221503854[8] = 0;
}
void h_24(double *state, double *unused, double *out_3178447534282387031) {
   out_3178447534282387031[0] = state[4];
   out_3178447534282387031[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5012907806422176583) {
   out_5012907806422176583[0] = 0;
   out_5012907806422176583[1] = 0;
   out_5012907806422176583[2] = 0;
   out_5012907806422176583[3] = 0;
   out_5012907806422176583[4] = 1;
   out_5012907806422176583[5] = 0;
   out_5012907806422176583[6] = 0;
   out_5012907806422176583[7] = 0;
   out_5012907806422176583[8] = 0;
   out_5012907806422176583[9] = 0;
   out_5012907806422176583[10] = 0;
   out_5012907806422176583[11] = 0;
   out_5012907806422176583[12] = 0;
   out_5012907806422176583[13] = 0;
   out_5012907806422176583[14] = 1;
   out_5012907806422176583[15] = 0;
   out_5012907806422176583[16] = 0;
   out_5012907806422176583[17] = 0;
}
void h_30(double *state, double *unused, double *out_571291877412285187) {
   out_571291877412285187[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1058785926285744773) {
   out_1058785926285744773[0] = 0;
   out_1058785926285744773[1] = 0;
   out_1058785926285744773[2] = 0;
   out_1058785926285744773[3] = 0;
   out_1058785926285744773[4] = 1;
   out_1058785926285744773[5] = 0;
   out_1058785926285744773[6] = 0;
   out_1058785926285744773[7] = 0;
   out_1058785926285744773[8] = 0;
}
void h_26(double *state, double *unused, double *out_6863851557016415626) {
   out_6863851557016415626[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5201050351095560078) {
   out_5201050351095560078[0] = 0;
   out_5201050351095560078[1] = 0;
   out_5201050351095560078[2] = 0;
   out_5201050351095560078[3] = 0;
   out_5201050351095560078[4] = 0;
   out_5201050351095560078[5] = 0;
   out_5201050351095560078[6] = 0;
   out_5201050351095560078[7] = 1;
   out_5201050351095560078[8] = 0;
}
void h_27(double *state, double *unused, double *out_64170337685064455) {
   out_64170337685064455[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1115977385514680138) {
   out_1115977385514680138[0] = 0;
   out_1115977385514680138[1] = 0;
   out_1115977385514680138[2] = 0;
   out_1115977385514680138[3] = 1;
   out_1115977385514680138[4] = 0;
   out_1115977385514680138[5] = 0;
   out_1115977385514680138[6] = 0;
   out_1115977385514680138[7] = 0;
   out_1115977385514680138[8] = 0;
}
void h_29(double *state, double *unused, double *out_4430381390248560458) {
   out_4430381390248560458[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1569017270600136957) {
   out_1569017270600136957[0] = 0;
   out_1569017270600136957[1] = 1;
   out_1569017270600136957[2] = 0;
   out_1569017270600136957[3] = 0;
   out_1569017270600136957[4] = 0;
   out_1569017270600136957[5] = 0;
   out_1569017270600136957[6] = 0;
   out_1569017270600136957[7] = 0;
   out_1569017270600136957[8] = 0;
}
void h_28(double *state, double *unused, double *out_3341351713126243052) {
   out_3341351713126243052[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3513381746469393617) {
   out_3513381746469393617[0] = 1;
   out_3513381746469393617[1] = 0;
   out_3513381746469393617[2] = 0;
   out_3513381746469393617[3] = 0;
   out_3513381746469393617[4] = 0;
   out_3513381746469393617[5] = 0;
   out_3513381746469393617[6] = 0;
   out_3513381746469393617[7] = 0;
   out_3513381746469393617[8] = 0;
}
void h_31(double *state, double *unused, double *out_7001933135537451205) {
   out_7001933135537451205[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5827258453328911554) {
   out_5827258453328911554[0] = 0;
   out_5827258453328911554[1] = 0;
   out_5827258453328911554[2] = 0;
   out_5827258453328911554[3] = 0;
   out_5827258453328911554[4] = 0;
   out_5827258453328911554[5] = 0;
   out_5827258453328911554[6] = 0;
   out_5827258453328911554[7] = 0;
   out_5827258453328911554[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3580642894513962333) {
  err_fun(nom_x, delta_x, out_3580642894513962333);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2148588102071341310) {
  inv_err_fun(nom_x, true_x, out_2148588102071341310);
}
void car_H_mod_fun(double *state, double *out_6563764178898461836) {
  H_mod_fun(state, out_6563764178898461836);
}
void car_f_fun(double *state, double dt, double *out_9217663828102947708) {
  f_fun(state,  dt, out_9217663828102947708);
}
void car_F_fun(double *state, double dt, double *out_4765097794067878070) {
  F_fun(state,  dt, out_4765097794067878070);
}
void car_h_25(double *state, double *unused, double *out_4673975711821749475) {
  h_25(state, unused, out_4673975711821749475);
}
void car_H_25(double *state, double *unused, double *out_1459547032221503854) {
  H_25(state, unused, out_1459547032221503854);
}
void car_h_24(double *state, double *unused, double *out_3178447534282387031) {
  h_24(state, unused, out_3178447534282387031);
}
void car_H_24(double *state, double *unused, double *out_5012907806422176583) {
  H_24(state, unused, out_5012907806422176583);
}
void car_h_30(double *state, double *unused, double *out_571291877412285187) {
  h_30(state, unused, out_571291877412285187);
}
void car_H_30(double *state, double *unused, double *out_1058785926285744773) {
  H_30(state, unused, out_1058785926285744773);
}
void car_h_26(double *state, double *unused, double *out_6863851557016415626) {
  h_26(state, unused, out_6863851557016415626);
}
void car_H_26(double *state, double *unused, double *out_5201050351095560078) {
  H_26(state, unused, out_5201050351095560078);
}
void car_h_27(double *state, double *unused, double *out_64170337685064455) {
  h_27(state, unused, out_64170337685064455);
}
void car_H_27(double *state, double *unused, double *out_1115977385514680138) {
  H_27(state, unused, out_1115977385514680138);
}
void car_h_29(double *state, double *unused, double *out_4430381390248560458) {
  h_29(state, unused, out_4430381390248560458);
}
void car_H_29(double *state, double *unused, double *out_1569017270600136957) {
  H_29(state, unused, out_1569017270600136957);
}
void car_h_28(double *state, double *unused, double *out_3341351713126243052) {
  h_28(state, unused, out_3341351713126243052);
}
void car_H_28(double *state, double *unused, double *out_3513381746469393617) {
  H_28(state, unused, out_3513381746469393617);
}
void car_h_31(double *state, double *unused, double *out_7001933135537451205) {
  h_31(state, unused, out_7001933135537451205);
}
void car_H_31(double *state, double *unused, double *out_5827258453328911554) {
  H_31(state, unused, out_5827258453328911554);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
