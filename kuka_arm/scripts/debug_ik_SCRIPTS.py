from sympy import *
from time import time
from mpmath import radians
import tf
from IK_Calcs import IK_Calcs
'''
Format of test case is [ [[EE position],[EE orientation as quaternions]],[WC location],[joint angles]]
You can generate additional test cases by setting up your kuka project and running `$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract positions and orientation (in quaternion xyzw) and lastly use link 5
to find the position of the wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              4:[[[2.11731,0.0888341,1.83512],              # J4,J5,J6
                  [0.777226, 0.234364, -0.0612385,0.580727]],
                  [1.84986, 0, 1.94645],
                  [0.0,0.0,0.0,0.67,0.49,1.19]],
              5:[[[2.15286, 0, 1.94653],
                  [0.0, -0.000148353, 0.0, 1.0]],
                  [1.84986, 0.0, 1.94645],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
              6:[[[2.11871, 0.0, 1.64667],                  ### ONLY J3
                [0.0, 0.0835083, 0.0, 0.996507]],
                [1.81994, 0.0, 1.6971],
                [0.0, 0.0, 0.17, 0.0, 0.0, 0.0]],
              7:[[[1.99022, 1.15968, 1.19411],              ### J1,J2,J3
                [-0.0527774, 0.195405, 0.255352, 0.945424]],
                [1.74987, 1.01964, 1.31423],
                [0.53, 0.26, 0.15, 0.0, 0.0, 0.0]],
              8:[[[1.99022, 1.15968, 1.19411],              ### J1,J2,J3,J4
                [0.176941, 0.251253, 0.200647, 0.930217]],
                [1.74987, 1.01964, 1.31423],
                [0.53, 0.26, 0.15, 0.49, 0.0, 0.0]],
              9:[[[1.89371, 1.17584, 1.09809],
                [0.126962, 0.455179, 0.235474, 0.849261]],
                [1.74987, 1.01964, 1.31423],
                [0.53, 0.26, 0.15, 0.49, 0.46, 0.0]],
              10:[[[1.89371, 1.17584, 1.09809],
                [0.469032, 0.511891, 0.0245747, 0.719286]],
                [1.74987, 1.01964, 1.31423],
                [0.53, 0.26, 0.15, 0.49, 0.46, 0.86]],
              11:[[[2.28271,0,2.72781],
                [0,-0.231536,0,0.972826]],
                [2.0122, 0, 2.59131],
                [0.0,0.24,-0.71,0.0,0.0,0.0]],              ## H2==J2, J3
              12:[[[2.56852, 0.0, 2.73449],
                [0.0,-0.455856, 0.0, 0.890054]],
                [2.39145, 0.0, 2.48861],
                [0.0,1.11,-2.06,0.0,0.0,0.0]],              
              13:[[[2.15286, 0.0, 1.94653],                         #J6 
                [0.874662, -7.19116e-05, 0.000129759, 0.484733]],
                [1.84986, 0.0, 1.94645],
                [0.0,0.0,0.0,0.0,0.0,2.13]],
              14:[[[1.98384, 0.0, 1.67467],                         #J5
                [0.0, 0.528127, 0.0, 0.849165]],
                [1.84986, 0.0, 1.94645],
                [0.0,0.0,0.0,0.0,1.11,0.0]],
              15:[[[2.15286, 0.0, 1.94653],                         #J4
                [0.998537, 8.02276e-06, 0.000148136, -0.0540789]],
                [1.84986, 0.0, 1.94645],
                [0.0,0.0,0.0,3.25,0.0,0.0]],
              16:[[[2.15286, 0.0, 1.94653],                         #NO MOVE
                [0.0, -0.000148353, 0.0, 1]],
                [1.84986, 0.0, 1.94645],
                [0.0,0.0,0.0,0.0,0.0,0.0]],
              17:[[[1.40328, 1.63267, 1.94653],                         #J1
                [6.18989E-05, -0.000134823, 0.417241, 0.908796]],
                [1.20578, 1.40289, 1.94645],
                [0.86,0.0,0.0,0.0,0.0,0.0]],
              18:[[[2.4902, 0.0, 1.06866],                         #j2
                [0.0, 0.217172, 0.0, 0.976133]],
                [2.21578, 0.0, 1.19712],
                [0.0,0.44,0.0,0.0,0.0,0.0]]}


# 13 only joint 6
# 14 only joint 5
# 15 only joint 4

# 4 has first 3 links at 0
# 5 everything is 0
# 6 only moves J3
# 7 moves J1-J2-J3
# 8 same as 7 but moves J4
# 9 same as 8 but moves J5
# 10 same as 9 but moves J6

# 11 elbow up pose with negative J3
# 12 elbow down pose with negative J3


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


def sub_test(suppress_output, myCode, px, py, pz, roll, pitch, yaw):
    ik = IK_Calcs()

    theta1 = None
    theta2 = None
    theta3 = None
    theta4 = None
    theta5 = None
    theta6 = None

    if myCode:
        theta1, theta2, theta3, theta4, theta5, theta6 = ik.calc_thetas(px, py, pz, roll, pitch, yaw)
    else:
        theta1, theta2, theta3, theta4, theta5, theta6 = ik.calc_thetas(px, py, pz, roll, pitch, yaw, ik.calc_theta_123_walkthrough)
     ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    T_total =  ik.calc_T0_G(theta1,theta2,theta3,theta4,theta5,theta6)#ik.T0_G.evalf(subs={ik.q1:theta1, ik.q2:theta2, ik.q3:theta3, ik.q4:theta4, ik.q5:theta5, ik.q6:theta6})
    wc = ik.wc

    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = wc#[1,1,1] # <--- Load your calculated WC values in this array
    your_ee = [T_total[0,3], T_total[1,3], T_total[2,3]]#T_total #[1,1,1] # <--- Load your calculated end effector value from your forward kinematics
   ########################################################################################
    print(theta1)
    print(theta2)
    print(theta3)
    print(theta4)
    print(theta5)
    print(theta6)
    print(your_ee)

def test_code(suppress_output, myCode, test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()

    px = req.poses[x].position.x
    py = req.poses[x].position.y
    pz = req.poses[x].position.z

    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [req.poses[x].orientation.x, req.poses[x].orientation.y,
            req.poses[x].orientation.z, req.poses[x].orientation.w])


    ik = IK_Calcs(True)

    theta1 = None
    theta2 = None
    theta3 = None
    theta4 = None
    theta5 = None
    theta6 = None

    if myCode:
        theta1, theta2, theta3, theta4, theta5, theta6 = ik.calc_thetas(px, py, pz, roll, pitch, yaw)
    else:
        theta1, theta2, theta3, theta4, theta5, theta6 = ik.calc_thetas(px, py, pz, roll, pitch, yaw, ik.calc_theta_123_walkthrough)
     ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    T_total =  ik.calc_T0_G(theta1,theta2,theta3,theta4,theta5,theta6)#ik.T0_G.evalf(subs={ik.q1:theta1, ik.q2:theta2, ik.q3:theta3, ik.q4:theta4, ik.q5:theta5, ik.q6:theta6})
    wc = ik.wc

    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = wc#[1,1,1] # <--- Load your calculated WC values in this array
    your_ee = [T_total[0,3], T_total[1,3], T_total[2,3]]#T_total #[1,1,1] # <--- Load your calculated end effector value from your forward kinematics
   ########################################################################################

    ## Error analysis
    if not suppress_output: print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time)) 

    ret = {}
    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(pow(wc_x_e,2) + pow(wc_y_e,2) + pow(wc_z_e,2))
        ret['wc_x_e'] = wc_x_e
        ret['wc_y_e'] = wc_y_e
        ret['wc_z_e'] = wc_z_e
        ret['wc_offset'] = wc_offset
        if not suppress_output:
            print ("\nWrist error for x position is: %04.8f" % wc_x_e)
            print ("Wrist error for y position is: %04.8f" % wc_y_e)
            print ("Wrist error for z position is: %04.8f" % wc_z_e)
            print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])

    ret['theta1'] = theta1.evalf()
    ret['theta2'] = theta2.evalf()
    ret['theta3'] = theta3.evalf()
    ret['theta4'] = theta4.evalf()
    ret['theta5'] = theta5.evalf()
    ret['theta6'] = theta6.evalf()

    ret['t_1_e'] = t_1_e
    ret['t_2_e'] = t_2_e
    ret['t_3_e'] = t_3_e
    ret['t_4_e'] = t_4_e
    ret['t_5_e'] = t_5_e
    ret['t_6_e'] = t_6_e
    if not suppress_output:
        print ("\nTheta 1 error is: %04.8f" % t_1_e)
        print ("Theta 2 error is: %04.8f" % t_2_e)
        print ("Theta 3 error is: %04.8f" % t_3_e)
        print ("Theta 4 error is: %04.8f" % t_4_e)
        print ("Theta 5 error is: %04.8f" % t_5_e)
        print ("Theta 6 error is: %04.8f" % t_6_e)
        print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
            \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
            \nconfirm whether your code is working or not**")
        print (" ")



    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(pow(ee_x_e,2) + pow(ee_y_e,2) + pow(ee_z_e,2))
        
        ret['ee_x_e'] = ee_x_e
        ret['ee_y_e'] = ee_y_e
        ret['ee_z_e'] = ee_z_e
        ret['ee_offset'] = ee_offset

        if not suppress_output:
            print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
            print ("End effector error for y position is: %04.8f" % ee_y_e)
            print ("End effector error for z position is: %04.8f" % ee_z_e)
            print ("Overall end effector offset is: %04.8f units \n" % ee_offset)


    return ret


def test_codea(myCode, test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()
    
    ########################################################################################
    ## 
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

    # T0_2 = simplify(T0_1 * T1_2)    # base_link to link_2
    # # print('T0_2')
    # T0_3 = simplify(T0_2 * T2_3)
    # # print('T0_3')
    # T0_4 = simplify(T0_3 * T3_4)
    # # print('T0_4')
    # T0_5 = simplify(T0_4 * T4_5)
    # # print('T0_5')
    # T0_6 = simplify(T0_5 * T5_6)
    # # print('T0_6')
    # T0_G = simplify(T0_6 * T6_G)
    # print('T0_G')
    T0_G = T0_1*T1_2*T2_3*T3_4*T4_5*T5_6*T6_G

    # correct the end effector orientation difference between DH table and the urdf file
    # by rotating about z pi, and about y -pi/2
    # add the extra row and column to make it a hom transform
    R_corr = simplify(Rot('Z', pi) * Rot('Y', -pi/2))
    T_corr = R_corr.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]))

    px = req.poses[x].position.x
    py = req.poses[x].position.y
    pz = req.poses[x].position.z

    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
        [req.poses[x].orientation.x, req.poses[x].orientation.y,
            req.poses[x].orientation.z, req.poses[x].orientation.w])

    print("Roll: " + str(roll))
    print("Pitch: " + str(pitch))
    print("Yaw: " + str(yaw))
    print("px: " + str(px))
    print("py: " + str(py))
    print("pz: " + str(pz))


    # start of proper IK code

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

    base_to_EE_RMat = ROT_EE*Rot_Error
    base_to_EE_RMat = base_to_EE_RMat.subs({'r':roll, 'p':pitch, 'y':yaw})

    # length from wrist center to end effector gripper
    dWc_g = 0.303#s[d7] #s[d6] - l

    # the vector of the z-axis of the end effector
    z_vect_EE = base_to_EE_RMat[:,2]


    # wrist center is located along the end effector's z-axis in the negative direction 
    wc = Matrix([[px], [py], [pz]]) - dWc_g * z_vect_EE
    

    # look down on link1 from above, project onto x-y plane to get first joint angle
    theta1 = atan2(wc[1], wc[0])
    
    def calcWlkSidesMyThetas(s,wc):
        side_a = 1.501
        # print("side_a: " + str(A))
        # print(side_a)
        # print('\n\n\n')

        Bx = sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35
        # print("BX: " + str(Bx))
        Bz = wc[2] - 0.75
        # print("Bz: " + str(Bz))
        
        side_b = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35), 2) + pow((wc[2] - 0.75), 2))
        # print("B: " + str(side_b))
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

        print(atan2(Bz,Bx))
        # equivalently => pi/2 - angle_a - 0.121437586958607
        theta2 = pi/2 - angle_a - atan2(Bz, Bx)
        # equivalently => pi/2 - angle_b + (-0.0359844600820516)
        theta3 = pi/2 - angle_b + atan2(s[a3], s[d4])

        return (theta2, theta3)

    def calcMySidesWlkThetas(s,wc):
        ######## MY CALCS
        # projecting link2-link3-link5(WC) triangle onto a vertical plane and calculating joint angles using law of cosines
        A = sqrt(s[a3]*s[a3] + s[d4]*s[d4])
        # print("SIDE A: " + str(A))
        Bx = sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - s[a1]
        Bz = wc[2] - s[d1]
        B = sqrt(Bx*Bx + Bz*Bz)
        C = s[a2]

        # internal angles of triangle formed by link2-link3-link5
        angle_a = acos((-pow(A,2) + pow(B,2) + pow(C,2))/(2*B*C))
        angle_b = acos((-pow(B,2) + pow(A,2) + pow(C,2))/(2*A*C))
        angle_c = acos((-pow(C,2) + pow(B,2) + pow(A,2))/(2*B*A))
        ########## END OF MY CALCS

        # use the formulas for theta2 and theta3 from the walkthrough but using my side calculations
        # print(atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35))
        wlk_mysides_T2 = pi/2 - angle_a - atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)
        wlk_mysides_T3 = pi/2 - (angle_b + 0.036)

        return (wlk_mysides_T2, wlk_mysides_T3)

    def calcMyThetas(s, wc):

        # projecting link2-link3-link5(WC) triangle onto a vertical plane and calculating joint angles using law of cosines
        A = sqrt(s[a3]*s[a3] + s[d4]*s[d4])
        Bx = sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - s[a1]
        Bz = wc[2] - s[d1]
        B = sqrt(Bx*Bx + Bz*Bz)
        # print("B: " + str(B))
        C = s[a2]

        # internal angles of triangle formed by link2-link3-link5
        angle_a = acos((-pow(A,2) + pow(B,2) + pow(C,2))/(2*B*C))
        angle_b = acos((-pow(B,2) + pow(A,2) + pow(C,2))/(2*A*C))
        angle_c = acos((-pow(C,2) + pow(B,2) + pow(A,2))/(2*B*A))
        
        # use geometric identities and directions found from inspecting the model calculate the joint angles
        # print("atan2(Bz,Bx) " + str(atan2(Bz,Bx)))
        theta2 = pi/2  - angle_a - atan2(Bz, Bx)
        # print("atan2(s[a3], s[d4]" + str(atan2(s[a3],s[d4])))
        theta3 = pi/2 - angle_b + atan2(s[a3], s[d4])

        myThetas = (theta2, theta3)
        return myThetas
    
    def calcWalkThrough(s, wc):
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
        # theta2_wlk = pi/2 - angle_a_wlk - atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)
        # theta3_wlk = pi/2 - (angle_b_wlk + 0.036)

        # print("atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)" + str(atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)))
        theta2 = pi/2 - angle_a - atan2(wc[2] - 0.75, sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - 0.35)
        theta3 = pi/2 - (angle_b + 0.036)
        # theta3 = pi/2 - angle_b - 0.35

        # use geometric identities and directions found from inspecting the model calculate the joint angles
        # theta2 = pi/2 - atan2(Bz, Bx) - angle_a
        # theta3 = pi/2 - angle_b + atan2(s[a3], s[d4])
        # theta3 = pi/2 - angle_b - 0.35
        # theta3 = pi/2 - angle_b - 0.0359844600820516
        
        # print(atan2(s[a3], s[d4]))
        # print(atan2(s[a3], s[d4]))
        # print(atan(s[a3]/s[d4]))
        # print(s[a3])
        # print(s[d4])
        # return None
        # theta3 = pi/2 - angle_b + atan2(s[a3], s[d4])
        
        # theta2 = theta2_wlk
        # theta3 = theta3_wlk
        if False:
            print(angle_a)
            print(angle_a_wlk)
            print(angle_b)
            print(angle_b_wlk)
            print(angle_c)
            print(angle_c_wlk)
            print(A)
            print(side_a)
            print(B)
            print(side_b)
            print(C)
            print(side_c)
            print(theta2)
            print(theta2_wlk)
            print(theta2_wlkorg)
            print(theta3)
            print(theta3_wlk)
            print(theta3_wlkorg)
            return None
        wlkThetas = (theta2, theta3)

        return wlkThetas


    def calcTheta456(theta1, theta2, theta3, base_to_EE_RMat,T0_1, T1_2, T3_4):
        # Create analytical rotation matrix for first 3 joints from the cascade of homogenous transformations
        # plug in the theta values for those joints found from the trig analysis above
        T0_3 = T0_1 * T1_2 * T2_3
        R0_3 = T0_3[0:3, 0:3].evalf(subs={q1:theta1, q2:theta2, q3:theta3})
        
        # create rotation matrix from link 3 to link 6 using the fact that 
        # base_to_EE_RMat = R0_3 * R3_6
        R3_6 = R0_3.inv("LU") * base_to_EE_RMat
        R3_6_analytical = simplify(T3_4[0:3,0:3]*T4_5[0:3,0:3]*T5_6[0:3,0:3])
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
    

    mSwT2, mSwT3 = calcMySidesWlkThetas(s,wc)
    wSmT2, wSmT3 = calcWlkSidesMyThetas(s,wc)

    mT2, mT3 = calcMyThetas(s, wc)
    mT4, mT5, mT6 = calcTheta456(theta1, mT2, mT3, base_to_EE_RMat, T0_1, T1_2, T2_3)

    wT2, wT3 = calcWalkThrough(s, wc)
    wT4, wT5, wT6 = calcTheta456(theta1, wT2, wT3, base_to_EE_RMat, T0_1, T1_2, T2_3)
    
    print("my sides my thetas: ")
    # print((theta1, mT2, mT3, mT4, mT5, mT6))
    print((mT2, mT3))
    print("walk sides walk thetas: ")
    # print((theta1, wT2, wT3, wT4, wT5, wT6))
    print((wT2, wT3))

    print("walk sides my thetas: ")
    print((wSmT2, wSmT3))

    print("my sides walk thetas: ")
    print((mSwT2, mSwT3))
    # return

    # myCode = True


    theta2 = None
    theta3 = None

    if myCode:
        theta2 = mT2
        theta3 = mT3

    else:
        theta2 = wT2
        theta3 = wT3

    theta4, theta5, theta6 = calcTheta456(theta1, theta2, theta3, base_to_EE_RMat, T0_1, T1_2, T2_3)
    ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    T_total =  T0_G.evalf(subs={q1:theta1, q2:theta2, q3:theta3, q4:theta4, q5:theta5, q6:theta6})


    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = wc#[1,1,1] # <--- Load your calculated WC values in this array
    your_ee = [T_total[0,3], T_total[1,3], T_total[2,3]]#T_total #[1,1,1] # <--- Load your calculated end effector value from your forward kinematics
   ########################################################################################

    ## Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(pow(wc_x_e,2) + pow(wc_y_e,2) + pow(wc_z_e,2))
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(pow(ee_x_e,2) + pow(ee_y_e,2) + pow(ee_z_e,2))
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)


def printfloat(var_key, dct1, dct2):
    if dct2 is None:
        if len(var_key) <= 7:
            print('{}:\t\t{:.6f}'.format(var_key, dct1[var_key]))
        else:
            print('{}:\t{:.6f}'.format(var_key, dct1[var_key]))
    else:
        if len(var_key) <= 7:
            print('{}:\t\t{:.6f}\t{:.6f}'.format(var_key, dct1[var_key], dct2[var_key]))
        else:
            print('{}:\t{:.6f}\t{:.6f}'.format(var_key, dct1[var_key], dct2[var_key]))
        

if __name__ == "__maina__":
    px = 0.0
    py = 0.0
    pz =0.0
    
    sub_test(False,False,)


# 17 WORKS


if __name__ == "__main__":

    tests = [3]
    for test_case_number in tests:
        print("test case: " + str(test_case_number))
        # Change test case number for different scenarios
        # test_case_number = 2
        # print(test_cases[6])
        walk = test_code(True, False, test_cases[test_case_number])
        # mine = None
        mine = test_code(True, True, test_cases[test_case_number])

        print("\nend effector errors: ")
        # print("x:\t{:.6f}\t{:.6f}".format(walk['ee_x_e']))
        printfloat('wc_x_e',  walk, mine)
        printfloat('wc_y_e',  walk, mine)
        printfloat('wc_z_e',  walk, mine)
        printfloat('wc_offset', walk, mine)
        print('\n')
        printfloat('ee_x_e',  walk, mine)
        printfloat('ee_y_e',  walk, mine)
        printfloat('ee_z_e',  walk, mine)
        printfloat('ee_offset', walk, mine)
        print('\n')
        printfloat('theta1', walk, mine)
        printfloat('theta2', walk, mine)
        printfloat('theta3', walk, mine)
        printfloat('theta4', walk, mine)
        printfloat('theta5', walk, mine)
        printfloat('theta6', walk, mine)


######## when joint 4 is used everything gets fucked


# 4 has first 3 links at 0
# 5 everything is 0
# 6 only moves J3
# 7 moves J1-J2-J3
# 8 same as 7 but moves J4
# 9 same as 8 but moves J5
# 10 same as 9 but moves J6