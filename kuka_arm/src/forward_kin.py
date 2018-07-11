#!/usr/bin/env python

import numpy as np
from numpy import array
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2
from sympy.matrices import Matrix

### Create symbols for joint variables
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')    # theta_1
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')


### KUKA KR210 ###
# DH Parameters
#   alpha: twist angles
#   a:      link lengths
#   d:      link offsets
#   q:      joint variables (the thetas)
s = {alpha0:    0,      a0:     0,      d1: 0.75,
     alpha1:    -pi/2,  a1:     0.35,   d2: 0,      q2: q2-pi/2,
     alpha2:    0,      a2:     1.25,   d3: 0,
     alpha3:    -pi/2,  a3:     -0.054, d4: 1.5,
     alpha4:    pi/2,   a4:     0,      d5: 0,
     alpha5:    -pi/2,  a5:     0,      d6: 0,
     alpha6:    0,      a6:     0,      d7: 0.303,  q7: 0    }


def HomTransform(q_1, alpha_0, d_1, a_0):
    transf = Matrix([   [       cos(q_1),                -sin(q_1),               0,              a_0              ],
                        [       sin(q_1)*cos(alpha_0),    cos(q_1)*cos(alpha_0),    -sin(alpha_0),   -sin(alpha_0)*d_1 ],
                        [       sin(q_1)*sin(alpha_0),    cos(q_1)*sin(alpha_0),    cos(alpha_0),    cos(alpha_0)*d_1  ],
                        [       0,                        0,                      0,              1               ]])
    return transf


#### Homogenous Transforms between neighboring links
# base_link to link1

T0_1 = HomTransform(q1, alpha0, d1, a0).T0_1.subs(s)
T1_2 = HomTransform(q2, alpha1, d2, a1).T1_2.subs(s)
T2_3 = HomTransform(q3, alpha2, d3, a2).T2_3.subs(s)
T3_4 = HomTransform(q4, alpha3, d4, a3).T3_4.subs(s)
T4_5 = HomTransform(q5, alpha4, d5, a4).T4_5.subs(s)
T5_6 = HomTransform(q6, alpha5, d6, a5).T5_6.subs(s)
T6_G = HomTransform(q7, alpha6, d7, a6).T6_G.subs(s)

# Composition of Homogeneous Tranforms
T0_2 = simplify(T0_1 * T1_2)    # base_link to link_2
T0_3 = simplify(T0_2 * T2_3)
T0_4 = simplify(T0_3 * T3_4)
T0_5 = simplify(T0_4 * T4_5)
T0_6 = simplify(T0_5 * T5_6)
T0_G = simplify(T0_6 * T6_G)

# Correction needed to account for orienation difference between definition of
# Gripper_link in URDF versus DH convention

R_z = Matrix([  [   cos(pi),    -sin(pi),   0,          0],
                [   sin(pi),    cos(pi),    0,          0],
                [   0,          0,          1,          0],
                [   0,          0,          0,          1]])

R_y = Matrix([[     cos(-pi/2), 0,          sin(-pi/2), 0],
                [   0,          1,          0,          0],
                [   -sin(pi/2), 0,          cos(-pi/2), 0],
                [   0,          0,          0,          1]])

R_corr = simplify(R_z * R_y)


#### Numerically evaluate transforms (compare this with output of tf_echo!)
print("T0_1 = ",T0_1.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_2 = ",T0_2.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_3 = ",T0_3.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_4 = ",T0_4.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_5 = ",T0_5.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_6 = ",T0_6.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))
print("T0_G = ",T0_G.evalf(subs={q1: -1.33, q2: -0.5, q3: 0.6, q4: -2.99, q5: 1.74, q6: -2.15}))

# Total Homogenous transform between Base_link and Gripper_link with
# orientation correction applied
T_total = simplify(T0_G * R_corr)
