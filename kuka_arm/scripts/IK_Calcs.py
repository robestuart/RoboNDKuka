from sympy import *
from mpmath import radians


class IK_Calcs:
    debug = False

    ### theta joint limits
    theta_min = {'theta1':-3.23, 'theta2':-0.79, 'theta3':-3.67, 'theta4':-6.11, 'theta5':-2.18, 'theta6':-6.11}
    theta_max = {'theta1':3.23, 'theta2':1.48, 'theta3':1.13, 'theta4':6.11,'theta5':2.18, 'theta6':6.11}

    ### defining symbols for sympy
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


    def __init__(self, simplify_=False):
        self.setup(simplify_)

    def forward_kin(self, theta1, theta2, theta3, theta4, theta5, theta6, print_=False):
        J1 = extract_pos(self.T0_1.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        J2 = extract_pos(self.T0_2.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        J3 = extract_pos(self.T0_3.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        J4 = extract_pos(self.T0_4.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        J5 = extract_pos(self.T0_5.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        J6 = extract_pos(self.T0_6.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))
        JG = extract_pos(self.T0_G.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6}))

        if print_:
            # print("J1:\t{}".format(J1))
            # print("J2:\t{}".format(J2))
            # print("J3:\t{}".format(J3))
            # print("J4:\t{}".format(J4))
            print("J5:\t{}".format(J5))
            # print("J6:\t{}".format(J6))
            print("JG:\t{}".format(JG))

        return (J1, J2, J3, J4, J5, J6, JG)

    def setup(self, simplify_=False):
        # Create individual transformation matricies between each joint
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

        print("simplifying T0_1")
        self.T0_1.simplify()
        print("simplifying T0_2")
        self.T0_2.simplify()
        print("simplifying T0_3")
        self.T0_3.simplify()
        print("simplifying T0_4")
        self.T0_4.simplify()
        print("simplifying T0_5")
        self.T0_5.simplify()
        print("simplifying T0_6")
        self.T0_6.simplify()
        print("simplifying T0_G")
        self.T0_G.simplify()


        # self.T0_G = (T0_1*T1_2*T2_3*T3_4*T4_5*T5_6*T6_G).subs(self.s)
        
        
        # correct the end effector orientation difference between DH table and the urdf file
        # by rotating about z pi, and about y -pi/2
        # add the extra row and column to make it a homogenous transform
        R_corr = simplify(Rot('Z', pi) * Rot('Y', -pi/2))
        T_corr = R_corr.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]))

        # Use the euler angles for the desired end effector pose
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

        # XYZ euler angles from base to EndEffector in world_fixed (extrinsic) frame
        ROT_EE = (ROT_z*ROT_y*ROT_x).subs(self.s)
        Rot_Error = ROT_z.subs(y, radians(180)) * ROT_y.subs(p, radians(-90))

        # store symbolic representation for iterative evaluation
        self.base_to_EE_RMat = (ROT_EE*Rot_Error)
        # length from wrist center to end effector gripper
        self.dWC_g = self.s[self.d7]

        # Create analytical rotation matrix for first 3 joints from the cascade of homogenous transformations
        # plug in the theta values for those joints found from the trig analysis above
        # pre-calculate the Rotation matrix from base to 3rd link and store symbolically for later
        self.R0_3 = self.T0_3[0:3,0:3]
        
        if simplify_ == True:
            print("simplifying the rotation matricies")
            # print("1/3")
            # self.T0_G.simplify()
            print("2/3")
            self.R0_3.simplify()
            print("3/3")
            self.base_to_EE_RMat.simplify()
            print("simplification completed")

    def calc_thetas(self, px, py, pz, roll, pitch, yaw, calc_func=None):
        """
        Given desired end effector position (px,py,pz) and end effector orientation (roll, pitch, yaw)
        calculate the joint angles using inverse kinematics
        """
        
        if self.debug: print("px:\t{}\npy:\t{}\npz:\t{}\nroll:\t{}\npitch:\t{}\nyaw:\t{}".format(px,py,pz,roll,pitch,yaw))
        
        if calc_func is None:
            calc_func = self.calc_theta_123

        # evaluate the rotation matrix from the base to end effector using desired roll pitch and yaw values
        base_to_EE_RMat_val =  self.base_to_EE_RMat.subs({'r':roll, 'p':pitch, 'y':yaw})
        # use the z_vector from the EE to point to the WC
        z_vect_EE = base_to_EE_RMat_val[:,2]
        # WC is located along the EE z-axis in the negative direction
        self.wc = Matrix([[px], [py], [pz]]) - self.dWC_g * z_vect_EE

        # solve for first 3 thetas using the wrist center and the DH table
        theta1, theta2, theta3 = calc_func(self.wc)

        # self.check_theta_123(theta1,theta2,theta3)

        # evaluate the rotation matrix from base to link 3 using the thetas solved above
        R0_3_val = self.R0_3.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3})

        # solve for remaining thetas
        theta4, theta5, theta6 = self.calc_theta_456(base_to_EE_RMat_val, R0_3_val, (theta1, theta2, theta3), px, py, pz)

        print("theta1:\t{}\ntheta2:\t{}\ntheta3:\t{}\ntheta4:\t{}\ntheta5:\t{}\ntheta6:\t{}".format(theta1,theta2,theta3,theta4,theta5,theta6))
        # self.check_theta_456(theta4,theta5,theta6)
        self._check_limit('theta1', theta1)
        self._check_limit('theta2', theta2)
        self._check_limit('theta3', theta3)
        self._check_limit('theta4', theta4)
        self._check_limit('theta5', theta5)
        self._check_limit('theta6', theta6)

        return (theta1, theta2, theta3, theta4, theta5, theta6)
    
    def _calc_theta2_up(self, angle_a, theta2_atan2):
        return pi/2 - angle_a - theta2_atan2

    def _calc_theta2_down(self, angle_a, theta2_atan2):
        return pi/2 - theta2_atan2 + angle_a

    def _calc_theta3_up(self, angle_b, theta3_atan2):
        return pi/2 - (angle_b + theta3_atan2)
    
    def _calc_theta3_down(self, angle_b, theta3_atan2):
        return 3*pi/2 - angle_b + theta3_atan2
        
    def calc_theta_123_walkthrough(self,wc):
        """
        Calculate theta 1,2,3 using the code presented in the project walk through
        """
        # look down on link1 from above, project onto x-y plane to get first joint angle
        theta1 = atan2(wc[1], wc[0])

        side_a = 1.501
        # print("side_a " + str(side_a))
        side_b = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35), 2) + pow((wc[2] - 0.75), 2))
        side_c = 1.25

        angle_a = acos((side_b * side_b + side_c * side_c - side_a*side_a) / (2*side_b*side_c))
        angle_b = acos((side_a * side_a + side_c * side_c - side_b*side_b) / (2*side_a*side_c))

        #### angle_c is not used
        # angle_c = acos((side_a * side_a + side_b * side_b - side_c*side_c) / (2*side_a*side_b))
        # print("angle_c: {}".format(angle_c))

        theta2_atan2_y = wc[2] - 0.75
        theta2_atan2_x = sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35

        theta2_atan2 = atan2(theta2_atan2_y, theta2_atan2_x)
        theta2 = self._calc_theta2_up(angle_a, theta2_atan2)#pi/2 - angle_a - theta2_atan2

        theta3_atan2 = 0.036
        theta3 = self._calc_theta3_up(angle_b, theta3_atan2)#pi/2 - (angle_b + theta3_atan2)
        
        if self.debug:
            print("\nWALKTHROUGH")
            print("side_a:\t{}".format(side_a))
            print("side_b:\t{}".format(side_b))
            print("side_c:\t{}".format(side_c))

            print("angle_a: {}".format(angle_a))
            print("angle_b: {}".format(angle_b))
            print("theta2 atan2:\t{}".format(theta2_atan2))
            print("theta3 atan2:\t{}".format(theta3_atan2))
            print("theta2:\t{}".format(theta2))
            print("theta3:\t{}".format(theta3))
        return (theta1, theta2, theta3)

    # def check_theta_123(self, theta1, theta2, theta3):
    #     """
    #     limits for joints 1-3:
    #     -3.23 < J1 < 3.23
    #     -0.79 < J2 < 1.48
    #     -3.67 < J3 < 1.13
    #     """
    #     self._check_theta('theta1', theta1, -3.23, 3.23)
    #     self._check_theta('theta2', theta2, -0.79, 1.48)
    #     self._check_theta('theta3', theta3, -3.67, 1.13)

        
    # def _check_theta(self,lbl, theta, t_min, t_max):
    #     if t_min < theta and theta < t_max:
    #         return true
    #     else:
    #         print(lbl + " IS OUTSIDE OF RANGE!!!!")
    #         print(theta)
    #         return false

    # def check_theta_456(self, theta4, theta5, theta6):
    #     """
    #     limits for joints 4-6:
    #     -6.11 < J4 < 6.11
    #     -2.18 < J5 < 2.18
    #     -6.11 < J6 < 6.11
    #     """
    #     self._check_theta('theta4', theta4, -6.11, 6.11)
    #     self._check_theta('theta5', theta5, -2.18, 2.18)
    #     self._check_theta('theta6', theta6, -6.11, 6.11)


    def _check_limit(self, theta_name, theta_val):
        if self.theta_min[theta_name] < theta_val and theta_val < self.theta_max[theta_name]:
            return true
        else:
            print(theta_name + " IS OUTSIDE RANGE\t" + str(theta_val))
            return false

    def calc_theta_123(self, wc):
        """
        Calculate theta 1,2,3 using the more numerically accurate equations I came up with
        """

        # look down on link1 from above, project onto x-y plane to get first joint angle
        theta1 = atan2(wc[1], wc[0])

        # projecting link2-link3-link5(WC) triangle onto a vertical plane and calculating joint angles using law of cosines
        A = sqrt(self.s[self.a3]*self.s[self.a3] + self.s[self.d4]*self.s[self.d4])
        Bx = sqrt(pow(wc[0],2) + pow(wc[1],2)) - self.s[self.a1]
        Bz = wc[2] - self.s[self.d1]
        B = sqrt(pow(Bx,2) + pow(Bz,2))
        C = self.s[self.a2]



        # internal angles of triangle formed by link2-link3-link5
        angle_a = acos((-pow(A,2) + pow(B,2) + pow(C,2))/(2*B*C))
        angle_b = acos((-pow(B,2) + pow(A,2) + pow(C,2))/(2*A*C))


        ##### angle_c is not used
        # angle_c = acos((-pow(C,2) + pow(B,2) + pow(A,2))/(2*B*A))
        # print("angle_c: {}".formt(angle_c))

        # use geometric identities and directions found from inspecting the model calculate the joint angles

        # up pose
        theta2_atan2 = atan2(Bz, Bx)
        theta2 = self._calc_theta2_up(angle_a, theta2_atan2)#pi/2  - angle_a - theta2_atan2
        theta3 = None
        if self._check_limit('theta2', theta2):
            theta3_atan2 = -1.0 * atan2(self.s[self.a3], self.s[self.d4])
            theta3 = self._calc_theta3_up(angle_b, theta3_atan2)#pi/2 - (angle_b + theta3_atan2)
        else:
            # down pose
            theta2 = self.calc_theta2_down(angle_a, theta2_atan2)
            theta3 = self.calc_theta3_down(angle_b, theta3_atan2)

        if self.debug:
            print("\nMINE")
            print("side_a:\t{}".format(A))
            print("side_b:\t{}".format(B))
            print("side_c:\t{}".format(C))

            print("angle_a: {}".format(angle_a))
            print("angle_b: {}".format(angle_b))
            print("theta2 atan2:\t{}".format(theta2_atan2))
            print("theta3_atan2\t{}".format(theta3_atan2))

            print("theta2:\t{}".format(theta2))
            print("theta3:\t{}".format(theta3))

        return (theta1, theta2, theta3)

    def calc_theta_456(self, base_to_EE_RMat_val, R0_3_val, theta123, px, py, pz, multi_sol_check=False):
        """
        Calculates theta 456 using the base to end effector rotation matrix and the rotation matrix from joint 0 to joint 3
        """
        theta1 = theta123[0]
        theta2 = theta123[1]
        theta3 = theta123[2]

        # create rotation matrix from link 3 to link 6 using the fact that 
        # base_to_EE_RMat = R0_3 * R3_6
        R3_6 = R0_3_val.transpose()*base_to_EE_RMat_val#R0_3_val.inv("LU") * base_to_EE_RMat_val
        
        # R3_6_analytical = simplify(T3_4[0:3,0:3]*T4_5[0:3,0:3]*T5_6[0:3,0:3])
        # print the analytical matrix to solve for euler angles using trig identities and atan2
        # print(R3_6_analytical)

        # use the rotation matrix for the last 3 joints to solve for the three joint angles
        theta4 = atan2(R3_6[2,2], -R3_6[0,2])
        # theta5 = atan2(sqrt(pow(R3_6[0,2],2) + pow(R3_6[2,2],2)), R3_6[1,2])
        theta5 = atan2(sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2]*R3_6[2,2]), R3_6[1,2])

        theta5_neg = atan2(-sqrt(pow(R3_6[0,2],2) + pow(R3_6[2,2],2)), R3_6[1,2])
        theta5_2 = atan2(sqrt(pow(R3_6[1,0],2) + pow(R3_6[1,1],2)), R3_6[1,2]) #equivalent atan2 calculation
        
        theta6 = atan2(-R3_6[1,1], R3_6[1,0])    ## 

        # sols = [(theta4, theta5, theta6), (theta4, theta5_neg, theta6), (theta4, theta5_2, theta6)]
        # enumerate_sols((theta4), (theta5, theta5_neg, theta5_2), (theta6))
        
        
        ### iterate through potential combinations of solutions and compare to desired end effector position
        if multi_sol_check:
            min_error = 1000.0
            min_error_idx = -1
            
            sols = enumerate_sols(self._check_2pi('theta4', [theta4.evalf()]), self._check_2pi('theta5', [theta5.evalf(), theta5_neg.evalf(), theta5_2.evalf()]), self._check_2pi('theta6',[theta6.evalf()]))
            
            for i,s in enumerate(sols):
                print("theta4:\t{}\ttheta5:\t{}\ttheta6:\t{}".format(s[0],s[1],s[2]))

                T0_G = self.calc_T0_G(theta1,theta2,theta3, s[0],s[1],s[2])
                errors = calc_error(T0_G, px, py, pz)
                print(errors[3])
                if errors[3] < min_error:
                    min_error = errors[3]
                    min_error_idx = i
            
            return sols[min_error_idx]
        else:
            return (theta4, theta5, theta6)


    def calc_T0_G(self, theta1, theta2, theta3, theta4, theta5, theta6):
        return self.T0_G.evalf(subs={self.q1:theta1, self.q2:theta2, self.q3:theta3, self.q4:theta4, self.q5:theta5, self.q6:theta6})

    def _check_2pi(self, theta_name, theta_val):
        """
        Finds all the equivalent angles within the joint ranges that would have produced the same atan2 values
        by adding or subtracting multiples of 2pi
        """
        possible = []
        if type(theta_val) is list or tuple:
            possible = list(theta_val)
        else:
            possible = [theta_val]
            theta_val = [theta_val]     # makes single value iterable
        
        mx = self.theta_max[theta_name]
        mn = self.theta_min[theta_name]
        # print(type(theta_val))
        # print(theta_val)
        
        for c, theta in enumerate(theta_val):
            # adding positive 2pi
            i = 1
            brk = False
            while(not brk):
                add_2pi = theta + 2.0*pi.evalf()*i

                if (add_2pi < mx):
                    possible.append(add_2pi)
                    i = i + 1
                else:
                    brk = True

            # adding negative 2pi
            i = 1
            brk = False
            while(not brk):
                add_2pi = theta - 2.0*pi.evalf()*i

                if (mn < add_2pi):
                    possible.append(add_2pi)
                    i+=1
                else:
                    brk = True

        return possible


def enumerate_sols(theta4, theta5, theta6):
    combos = []
    for t4 in theta4:
        for t5 in theta5:
            for t6 in theta6:
                combos.append((t4,t5,t6))

    return combos

def calc_error(T0_G, px, py, pz):
    x = abs(T0_G[0,3]-px)
    y = abs(T0_G[1,3] - py)
    z = abs(T0_G[2,3]-pz)

    offset = sqrt(pow(T0_G[0,3]-px,2) + pow(T0_G[1,3] - py,2) + pow(T0_G[2,3]-pz,2))
    return (x, y, z, offset)
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
def HomTransform(q_1, alpha_0, d_1, a_0):
    """
    Utility method for evaluation of Homogenous Transforms
    """
    transf = Matrix([   [       cos(q_1),                -sin(q_1),               0,              a_0              ],
                        [       sin(q_1)*cos(alpha_0),    cos(q_1)*cos(alpha_0),    -sin(alpha_0),   -sin(alpha_0)*d_1 ],
                        [       sin(q_1)*sin(alpha_0),    cos(q_1)*sin(alpha_0),    cos(alpha_0),    cos(alpha_0)*d_1  ],
                        [       0,                        0,                      0,              1               ]])
    return transf


def extract_pos(T):
    return (T[0,3], T[1,3], T[2,3])





if __name__ == "__main__":
    ##### forward kinematics verified!!!!
    ##### joint5 and gripper joint are the ones to check
    ik = IK_Calcs()
    theta1 = 0.72
    theta2 = 0.93
    theta3 = -1.32
    theta4 = -1.6
    theta5 = -1.41
    theta6 = 1.57
    ik.forward_kin(theta1, theta2, theta3, theta4, theta5, theta6, True)