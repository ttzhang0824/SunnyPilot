#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1443495834072089830);
void live_err_fun(double *nom_x, double *delta_x, double *out_7269670487349612617);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_5168123923996786904);
void live_H_mod_fun(double *state, double *out_1452774069512641077);
void live_f_fun(double *state, double dt, double *out_6830993965878836661);
void live_F_fun(double *state, double dt, double *out_2745763742263175106);
void live_h_4(double *state, double *unused, double *out_8974275436426955909);
void live_H_4(double *state, double *unused, double *out_3251180476872113937);
void live_h_9(double *state, double *unused, double *out_4146100310906222837);
void live_H_9(double *state, double *unused, double *out_4036038458392333533);
void live_h_10(double *state, double *unused, double *out_112020028582643827);
void live_H_10(double *state, double *unused, double *out_1724120192745202477);
void live_h_12(double *state, double *unused, double *out_2427243725057012526);
void live_H_12(double *state, double *unused, double *out_8814305219794704683);
void live_h_35(double *state, double *unused, double *out_3374267511616249904);
void live_H_35(double *state, double *unused, double *out_115481580500493439);
void live_h_32(double *state, double *unused, double *out_5546416164689322160);
void live_H_32(double *state, double *unused, double *out_3339031481602674495);
void live_h_13(double *state, double *unused, double *out_1839926880760156660);
void live_H_13(double *state, double *unused, double *out_9164343041642234366);
void live_h_14(double *state, double *unused, double *out_4146100310906222837);
void live_H_14(double *state, double *unused, double *out_4036038458392333533);
void live_h_33(double *state, double *unused, double *out_6560067115827767053);
void live_H_33(double *state, double *unused, double *out_3266038585139351043);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}