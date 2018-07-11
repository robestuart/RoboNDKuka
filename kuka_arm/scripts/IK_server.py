#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

def Rot(dir, angle):
    if dir == 'Z':
        return Matrix([ [   cos(angle),     -sin(angle),    0],
                        [   sin(angle),     cos(angle),     0],
                        [   0,              0,              1]])
    elif dir == 'Y':
        return Matrix([ [   cos(angle),     0,          sin(angle)  ], 
                        [   0,              1,          0           ],
                        [   -sin(angle),    0,          cos(angle)]])
    elif dir == 'X':
        return Matrix([ [   1,              0,          0           ],
                        [   0,              cos(angle), -sin(angle) ],
                        [   0,              sin(angle), cos(angle)  ]])

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print("No valid poses received")
        return -1
    else:

        ### Your FK code here
        # Create symbols	
        ### Create symbols for joint variables
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')    # theta_1
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')



	    # Create Modified DH parameters
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
	    # Define Modified DH Transformation matrix
        def HomTransform(q_1, alpha_0, d_1, a_0):
            transf = Matrix([   [       cos(q_1),                -sin(q_1),               0,              a_0              ],
                                [       sin(q_1)*cos(alpha_0),    cos(q_1)*cos(alpha_0),    -sin(alpha_0),   -sin(alpha_0)*d_1 ],
                                [       sin(q_1)*sin(alpha_0),    cos(q_1)*sin(alpha_0),    cos(alpha_0),    cos(alpha_0)*d_1  ],
                                [       0,                        0,                      0,              1               ]])
            return transf
	
        T0_1 = HomTransform(q1, alpha0, d1, a0).subs(s)
        T1_2 = HomTransform(q2, alpha1, d2, a1).subs(s)
        T2_3 = HomTransform(q3, alpha2, d3, a2).subs(s)
        T3_4 = HomTransform(q4, alpha3, d4, a3).subs(s)
        T4_5 = HomTransform(q5, alpha4, d5, a4).subs(s)
        T5_6 = HomTransform(q6, alpha5, d6, a5).subs(s)
        T6_G = HomTransform(q7, alpha6, d7, a6).subs(s)

	    # Create individual transformation matrices
	
        T0_2 = simplify(T0_1 * T1_2)    # base_link to link_2
        T0_3 = simplify(T0_2 * T2_3)
        T0_4 = simplify(T0_3 * T3_4)
        T0_5 = simplify(T0_4 * T4_5)
        T0_6 = simplify(T0_5 * T5_6)
        T0_G = simplify(T0_6 * T6_G)

	    # Extract rotation matrices from the transformation matrices
        R_z = Matrix([  [   cos(pi),    -sin(pi),   0,          0],
                [   sin(pi),    cos(pi),    0,          0],
                [   0,          0,          1,          0],
                [   0,          0,          0,          1]])

        R_y = Matrix([[     cos(-pi/2), 0,          sin(-pi/2), 0],
                        [   0,          1,          0,          0],
                        [   -sin(pi/2), 0,          cos(-pi/2), 0],
                        [   0,          0,          0,          1]])

        R_corr = simplify(Rot('Z', pi) * Rot('Y', -pi/2))#simplify(R_z * R_y)
        T_corr = R_corr.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]))

        ###

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here
            # Compensate for rotation discrepancy between DH parameters and Gazebo
            
            # Take the euler angles from the simulation for the end effector pose
            # use these to calculate the rotation matrix from the base_link to the end effector
            #Rrpy = Rot(Z, yaw)*Rot(Y,pitch)*Rot(X, roll)*R_corr
            
            base_to_EE_RMat = Rot('Z', yaw)*Rot('Y',pitch)*Rot('X', roll)*R_corr
            # extract the z-normal components from the Rrpy matrix
            # nx = Rrpy[0,2]
            # ny = Rrpy[1,2]
            # nz = Rrpy[2,2]
            
            dWc_g = s[d7] #s[d6] - l

            n = base_to_EE_RMat[:,2]

            wc = Matrix([[px],
                        [py],
                        [pz]])
            wc = wc-dWc_g*n

            # Calculate joint angles using Geometric IK method
            
                    ## calculating triangle angles for IK
            A = sqrt(s[a3]**2 + s[d4]**2)
            Bx = sqrt(wc[0]**2 + wc[1]**2) - s[a1]
            Bz = wc[2] - s[d1]
            B = sqrt(Bx**2 + Bz**2)
            C = s[a2]

            angle_a = acos((-A**2 + B**2 + C**2)/(2*B*C))
            angle_b = acos((-B**2 + A**2 + C**2)/(2*A*C))
            angle_c = acos((-C**2 + B**2 + A**2)/(2*B*A))

            theta1 = atan2(wc[1], wc[0])
            theta2 = pi/2 - atan2(Bz, Bx) - angle_a
            theta3 = pi/2 - angle_b + atan2(s[a3], s[d4])

            # Create rotation matrix for first 3 links using the solved for joint angles
            T0_3 = T0_1 * T1_2 * T2_3
            R0_3 = T0_3[0:3, 0:3].evalf(subs={q1:theta1, q2:theta2, q3:theta3})

            # create rotation matrix from link 3 to link 6 using the fact that 
            # base_to_EE_RMat = R0_3 * R3_6
            R3_6 = R0_3.inv("LU") * base_to_EE_RMat

            # beta = atan2(-R3_6[2,0], sqrt(R_6[0,0]**2 + R_6[1,0]**2))
            # gamma = atan2(R3_6[2,1], R3_6[2,2])
            # alpha = atan2(R3_6[1,0], R3_6[0,0])

            theta4 = atan2(R3_6[2,2], -R3_6[0,2])
            theta5 = atan2(sqrt(R3_6[0,2]**2 + R3_6[2,2]**2), R3_6[1,2])
            theta6 = atan2(-R3_6[1,1], R3_6[1,0]) 

            ###

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print("Ready to receive an IK request")
    rospy.spin()

if __name__ == "__main__":
    IK_server()
