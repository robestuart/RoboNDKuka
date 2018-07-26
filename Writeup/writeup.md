## Project: Kinematics Pick & Place

[//]: # (Image References)

[robo_kin]: kinematics_1.jpg
[elb_up2]: kinematics_2.jpg
[elb_up3]: kinematics_3.jpg
[elb_dwn2]: kinematics_4.jpg
[elb_dwn3]: kinematics_5.jpg
[rob_pick_place]: Robot_pick_place.png
[Euler_Rot_Mat]: EulerAnglesFromRotMat.jpg

---
### Writeup / README

## Kinematic Analysis

I used the XACRO file to determine the distances between links and fill out the modified DH table according to the orientations of the joint axes that I have drawn.

![alt text][robo_kin]

Below is the modified DH table I derived from the XACRO file and my definition of the axes orientations and locations

Links | alpha(i-1) | a(i-1) | d(i-1) | theta(i)
--- | --- | --- | --- | ---
0->1 | 0 | 0 | 0.75 | q1
1->2 | -pi/2 | 0.35 | 0 | q2 - pi/2
2->3 | 0 | 1.25 | 0 | q3
3->4 |  -pi/2 | -0.054 | 1.5 | q4
4->5 | pi/2 | 0 | 0 | q5
5->6 | -pi/2 | 0 | 0 | q6
6->EE | 0 | 0 | 0.303 | 0

## Forward Kinematics

I used python to create homogeneous transforms for each joint and proceeded with creating a composite transform from the base link to the end gripper by post-multiplying each transform:
```
def HomTransform(q_1, alpha_0, d_1, a_0):
    """
    Utility method for evaluation of Homogenous Transforms
    """
    transf = Matrix([   [       cos(q_1),                -sin(q_1),               0,              a_0              ],
                        [       sin(q_1)*cos(alpha_0),    cos(q_1)*cos(alpha_0),    -sin(alpha_0),   -sin(alpha_0)*d_1 ],
                        [       sin(q_1)*sin(alpha_0),    cos(q_1)*sin(alpha_0),    cos(alpha_0),    cos(alpha_0)*d_1  ],
                        [       0,                        0,                      0,              1               ]])
    return transf


T0_1 = HomTransform(self.q1, self.alpha0, self.d1, self.a0)
T1_2 = HomTransform(self.q2, self.alpha1, self.d2, self.a1)
T2_3 = HomTransform(self.q3, self.alpha2, self.d3, self.a2)
T3_4 = HomTransform(self.q4, self.alpha3, self.d4, self.a3)
T4_5 = HomTransform(self.q5, self.alpha4, self.d5, self.a4)
T5_6 = HomTransform(self.q6, self.alpha5, self.d6, self.a5)
T6_G = HomTransform(self.q7, self.alpha6, self.d7, self.a6)


self.T0_1 = T0_1.subs(self.s)
self.T0_2 = self.T0_1*T1_2.subs(self.s)
self.T0_3 = self.T0_2*T2_3.subs(self.s)
self.T0_4 = self.T0_3*T3_4.subs(self.s)
self.T0_5 = self.T0_4*T4_5.subs(self.s)
self.T0_6 = self.T0_5*T5_6.subs(self.s)
self.T0_G = self.T0_6*T6_G.subs(self.s)
```

## Decoupled Inverse Kinematics

To cacluate the wrist center from the desired end effector orientation and location I needed to correct the end effector orientation due to the difference between the model and DH parameter axis orientation.  To do so I had to rotate about the end effector's z axis by pi and y axis by -pi/2.

```
ROT_EE = (ROT_z*ROT_y*ROT_x).subs(self.s)
Rot_Error = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))
self.base_to_EE_RMat = (ROT_EE*Rot_Error)
```

The location of the wrist center was calculated by a translation up the robot arm via the 
z-axis of the end effector and link length d7:

```
base_to_EE_RMat_val =  self.base_to_EE_RMat.subs({'r':roll, 'p':pitch, 'y':yaw})
self.dWC_g = self.s[self.d7]
z_vect_EE = base_to_EE_RMat_val[:,2]
self.wc = Matrix([[px], [py], [pz]]) - self.dWC_g * z_vect_EE
```

Once the wrist center was found the first three joint angles could be calculated using trigonometery.  The pictures below depict the trigonometry used to calculate the joint angles of the robot arm in the up pose and down poses for joints J2 and J3.  For the actual program only the "elbow up" configuration was used as there was no need to switch the pose of the robot.

![alt text][elb_up2]
![alt text][elb_up3]
![alt text][elb_dwn2]
![alt text][elb_dwn3]

To determine J4, J5, and J6 joint angles I used the fact that `R0_G = R0_3 * R3_6` to find the rotation matrix for the last three joints I took the inverse of R0_3 which is equal to its transpose because it is a unitary matrix.  This results in the following equation: `R3_6 = R0_3.transpose()*R0_G`

With a numerical solution to `R3_6` you can determine the joint angles by matching it up with the equivalent analytical euler angle rotation matrix and use some trig identities as shown below:

![alt text][Euler_Rot_Mat]


## Project Implementation



To implement the code I thought it would be easier to write a separate module to keep things organized.  I created a module called IK_Calcs.py.  After having implemented my calculations I ran the debugger and found that in some of the test cases I was not getting the correct answer and in fact it was off by 4 orders of magnitude.  I checked the walkthrough equations and determined them to be numerically equivalent to mine except for some numerical approximations the walkthrough makes as opposed to my numerically more precise equations.  After quite a lot of struggle and frustration I finally checked the slack channel and found that a student suggested the portion of the assignment directing us to use the sympy inverse with LU decomposition does not produce the correct answers and that using .transpose will fix things.  In fact they were correct and after replacing .inv("LU") with .transpose my code worked perfectly.  Using the transpose makes sense from a mathematical standpoint because rotation matricies are unitary transformations and therefore the transpose is equivalent to the inverse.  Since inverses are computationally expensive operations I imagine the problem with using the inverse is that numerical precision becomes a factor with the inverse operation whereas the transpose operation does not have this issue.  Please change this in the lesson as I wasted many hours trying to figure out what I was doing wrong.  There are also a couple of other errors in the lesson and the walkthrough which create confusion and I think should be addressed that I have brought up with my mentor.

I have some extra code that was not used which could help solve for the other possible solutions to the kinematic problem but since it was not needed in this assignment I did not bother troubleshooting them.


Below you can see the successul result from 10 trials of performing the pick and place operations:

![alt text][rob_pick_place]


