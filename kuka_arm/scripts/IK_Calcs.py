
class IK_Calcs:
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


    def __init__(self):
         # Create individual transformation matricies
        T0_1 = HomTransform(self.q1, self.alpha0, self.d1, self.a0)
        T1_2 = HomTransform(self.q2, self.alpha1, self.d2, self.a1)
        T2_3 = HomTransform(self.q3, self.alpha2, self.d3, self.a2)
        T3_4 = HomTransform(self.q4, self.alpha3, self.d4, self.a3)
        T4_5 = HomTransform(self.q5, self.alpha4, self.d5, self.a4)
        T5_6 = HomTransform(self.q6, self.alpha5, self.d6, self.a5)
        T6_G = HomTransform(self.q7, self.alpha6, self.d7, self.a6)

        T0_G = T0_1*T1_2*T2_3*T3_4*T4_5*T5_6*T6_G.sub(s)

        # correct the end effector orientation difference between DH table and the urdf file
        # by rotating about z pi, and about y -pi/2
        # add the extra row and column to make it a hom transform
        R_corr = simplify(Rot('Z', pi) * Rot('Y', -pi/2))
        T_corr = R_corr.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]))


        # Take the euler angles from the simulation for the end effector pose
        # use these to calculate the rotation matrix from the base_link to the end effector including the correction rotations
        # extrinsic euler angles XYZ
        # base_to_EE_RMat = Rot('Z', yaw)*Rot('Y',pitch)*Rot('X', roll)*R_corr

        r, p, y = symbols('r p y')

        ROT_x = Matrix([[1,     0,      0],
                        [0,     cos(r), -sin(r)],
                        [0, sin(r), cos(r)]])
        ROT_y = Matrix([[cos(p), 0, sin(p)],
                        [0,     1,      0],
                        [-sin(p),   0,  cos(p)]])
        ROT_z = Matrix([[cos(y), -sin(y), 0],
                        [sin(y), cos(y), 0],
                        [0,     0,      1]])

        ROT_EE = ROT_z*ROT_y*ROT_x
        Rot_Error = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))

        # store symbolic representation for iterative evaluation
        self.base_to_EE_RMat = simplify(ROT_EE*Rot_Error)
        # length from wrist center to end effector gripper
        self.dWC_g = s[d7]

        # Create analytical rotation matrix for first 3 joints from the cascade of homogenous transformations
        # plug in the theta values for those joints found from the trig analysis above
        # pre-calculate the Rotation matrix from base to 3rd link and store symbolically for later
        T0_3= simplify(T0_1*T1_2*T2_3).sub(s)
        self.R0_3 = T0_3[0:3,0:3]

    def calc_thetas(self, roll, pitch, yaw):
        """
        Calculates joint angles using my calculations
        """
        # evaluate the rotation matrix from the base to end effector using roll pitch and yaw values
        base_to_EE_RMat_val =  self.base_to_EE_RMat.evalf({'r':roll, 'p':pitch, 'y':yaw})
        # use the z_vector from the EE to point to the WC
        z_vect_EE = base_to_EE_RMat_val[:,2]
        # WC is located along the EE z-axis in the negative direction
        wc = Matrix([[px], [py], [pz]]) - self.dWc_g * z_vect_EE

        # solve for first 3 thetas using the wrist center and the DH table
        theta1, theta2, theta3 = calc_theta_123(wc)

        # evaluate the rotation matrix from base to link 3 using the thetas solved above
        R0_3_val = self.R0_3.evalf(subs={q1:theta1, q2:theta2, q3:theta3})

        # solve for remaining thetas
        theta4, theta5, theta6 = calc_theta_456(base_to_EE_RMat_val, R0_3_val)


        return (theta1, theta2, theta3, theta4, theta5, theta6)
        
    def calc_thetas_walkthrough(self, roll, pitch, yaw):
        """
        Calculates joint angles using the walkthrough calculations
        """

        # evaluate the rotation matrix from the base to end effector using roll pitch and yaw values
        base_to_EE_RMat_val =  self.base_to_EE_RMat.evalf({'r':roll, 'p':pitch, 'y':yaw})
        # use the z_vector from the EE to point to the WC
        z_vect_EE = base_to_EE_RMat_val[:,2]
        # WC is located along the EE z-axis in the negative direction
        wc = Matrix([[px], [py], [pz]]) - self.dWc_g * z_vect_EE

        # solve for first 3 thetas using the wrist center and the DH table
        theta1, theta2, theta3 = calc_theta_123_walkthrough(wc)

        # evaluate the rotation matrix from base to link 3 using the thetas solved above
        R0_3_val = self.R0_3.evalf(subs={q1:theta1, q2:theta2, q3:theta3})

        # solve for remaining thetas
        theta4, theta5, theta6 = calc_theta_456(base_to_EE_RMat_val, R0_3_val)


        return (theta1, theta2, theta3, theta4, theta5, theta6)


    @staticmethod
    def Rot(dir, angle):
        """
        Utility method for numerical evaluation of rotation matricies
        """
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


    # Define Modified DH Transformation matrix
    @staticmethod
    def HomTransform(q_1, alpha_0, d_1, a_0):
        """
        Utility method for evaluation of Homogenous Transforms
        """
        transf = Matrix([   [       cos(q_1),                -sin(q_1),               0,              a_0              ],
                            [       sin(q_1)*cos(alpha_0),    cos(q_1)*cos(alpha_0),    -sin(alpha_0),   -sin(alpha_0)*d_1 ],
                            [       sin(q_1)*sin(alpha_0),    cos(q_1)*sin(alpha_0),    cos(alpha_0),    cos(alpha_0)*d_1  ],
                            [       0,                        0,                      0,              1               ]])
        return transf


    def calc_theta_123_walkthrough(self,wc):
        # look down on link1 from above, project onto x-y plane to get first joint angle
        theta1 = atan2(self.wc[1], self.wc[0])

        side_a = 1.501
        # print("side_a: " + str(A))
        # print(side_a)
        # print('\n\n\n')
        side_b = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35), 2) + pow((wc[2] - 0.75), 2))
        # print("side_B: " + str(side_b))
        side_c = 1.25

        # print(A)
        # print(side_a)
        # print(B)
        # print(side_b)
        # print(C)
        # print(side_c)

        angle_a = acos((side_b * side_b + side_c * side_c - side_a*side_a) / (2*side_b*side_c))
        angle_b = acos((side_a * side_a + side_c * side_c - side_b*side_b) / (2*side_a*side_c))
        angle_c = acos((side_a * side_a + side_b * side_b - side_c*side_c) / (2*side_a*side_b))
        
        theta2 = pi/2 - angle_a - atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)
        theta3 = pi/2 - (angle_b + 0.036)

        return (theta1, theta2, theta3)



    def calc_theta_123(self, wc):
        # look down on link1 from above, project onto x-y plane to get first joint angle
        theta1 = atan2(self.wc[1], self.wc[0])

        # projecting link2-link3-link5(WC) triangle onto a vertical plane and calculating joint angles using law of cosines
        A = sqrt(self.s[a3]*self.s[a3] + self.s[d4]*self.s[d4])
        Bx = sqrt(pow(self.wc[0],2) + pow(self.wc[1],2)) - self.s[a1]
        Bz = self.wc[2] - self.s[d1]
        B = sqrt(pow(Bx,2) + pow(Bz,2))
        # print("B: " + str(B))
        C = self.s[a2]

        # internal angles of triangle formed by link2-link3-link5
        angle_a = acos((-pow(A,2) + pow(B,2) + pow(C,2))/(2*B*C))
        angle_b = acos((-pow(B,2) + pow(A,2) + pow(C,2))/(2*A*C))
        angle_c = acos((-pow(C,2) + pow(B,2) + pow(A,2))/(2*B*A))
        
        # use geometric identities and directions found from inspecting the model calculate the joint angles
        # print("atan2(Bz,Bx) " + str(atan2(Bz,Bx)))
        theta2 = pi/2  - angle_a - atan2(Bz, Bx)
        # print("atan2(s[a3], s[d4]" + str(atan2(s[a3],s[d4])))
        theta3 = pi/2 - angle_b + atan2(self.s[a3], self.s[d4])

        return (theta1, theta2, theta3)

    def calc_theta_456(self, base_to_EE_RMat_val, R0_3_val):
        # create rotation matrix from link 3 to link 6 using the fact that 
        # base_to_EE_RMat = R0_3 * R3_6
        R3_6 = R0_3.inv("LU") * base_to_EE_RMat_val
        
        # R3_6_analytical = simplify(T3_4[0:3,0:3]*T4_5[0:3,0:3]*T5_6[0:3,0:3])
        # print the analytical matrix to solve for euler angles using trig identities and atan2
        # print(R3_6_analytical)

        # use the rotation matrix for the last 3 joints to solver for the three joint angles
        theta4 = atan2(R3_6[2,2], -R3_6[0,2])
        theta5 = atan2(sqrt(pow(R3_6[0,2],2) + pow(R3_6[2,2],2)), R3_6[1,2])
        # theta5_neg = atan2(-sqrt(pow(R3_6[0,2],2) + pow(R3_6[2,2],2)), R3_6[1,2])
        # theta5_2 = atan2(sqrt(pow(R3_6[1,0],2) + pow(R3_6[1,1],2)), R3_6[1,2]) equivalent atan2 calculation
        theta6 = atan2(-R3_6[1,1], R3_6[1,0])    ## 
        ## 
        return (theta4, theta5, theta6)
