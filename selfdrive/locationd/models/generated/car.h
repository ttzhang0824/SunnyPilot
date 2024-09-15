#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3580642894513962333);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2148588102071341310);
void car_H_mod_fun(double *state, double *out_6563764178898461836);
void car_f_fun(double *state, double dt, double *out_9217663828102947708);
void car_F_fun(double *state, double dt, double *out_4765097794067878070);
void car_h_25(double *state, double *unused, double *out_4673975711821749475);
void car_H_25(double *state, double *unused, double *out_1459547032221503854);
void car_h_24(double *state, double *unused, double *out_3178447534282387031);
void car_H_24(double *state, double *unused, double *out_5012907806422176583);
void car_h_30(double *state, double *unused, double *out_571291877412285187);
void car_H_30(double *state, double *unused, double *out_1058785926285744773);
void car_h_26(double *state, double *unused, double *out_6863851557016415626);
void car_H_26(double *state, double *unused, double *out_5201050351095560078);
void car_h_27(double *state, double *unused, double *out_64170337685064455);
void car_H_27(double *state, double *unused, double *out_1115977385514680138);
void car_h_29(double *state, double *unused, double *out_4430381390248560458);
void car_H_29(double *state, double *unused, double *out_1569017270600136957);
void car_h_28(double *state, double *unused, double *out_3341351713126243052);
void car_H_28(double *state, double *unused, double *out_3513381746469393617);
void car_h_31(double *state, double *unused, double *out_7001933135537451205);
void car_H_31(double *state, double *unused, double *out_5827258453328911554);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}