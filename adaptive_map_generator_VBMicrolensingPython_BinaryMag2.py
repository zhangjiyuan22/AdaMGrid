import numpy as np
#from MulensModel.binarylens import BinaryLens
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpathes

import multiprocessing as mp
import time
import sys
import ctypes

#import VBBinaryLensing
import VBMicrolensing
import os

#s = 0.93#0.8
#q = 1.2e-4#0.1
#rho = 1e-2
box_size = 3.5#0.6#3.5#0.4#10.0#0.8

#bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)

#map_set_name = 'kb192991_close_densed_beltlike'
#map_set_name = 'kb192991_resonant'
#map_set_name = 'kb220440_publish_dense_close'
#map_set_name = 'kb220440_resonant_extreme_detail'
#map_set_name = 'ob161195_publish_dense'
#map_set_name = 'kb220886_dense'
#map_set_name = 'kb231190_sparse_close_singlerho_limb'
#map_set_name = 'kb231431_logs_minus0p15_to_0p15_logq_minus5p2_to_minus3p5_logrho_minus2p5'
#map_set_name = 'kb231431_logs_minus0p15_to_0p15_logq_minus5p2_to_minus3p5_logrho_minus2p7_to_minus2p6'
#map_set_name = 'kb210736_logs_0_to_0p65_dlogs_0p01_logq_minus6_to_minus2_dlogq_0p1_logrho_minus3p7_to_minus1p6_dlogrho_0p3_layer_12'
#map_set_name = 'kb210736_logs_minus0p65_to_0_dlogs_0p01_logq_minus6_to_minus2_dlogq_0p05_logrho_minus3p7_to_minus1p6_dlogrho_0p3_layer_12'
#map_set_name = 'kb221684_logs_minus1p5_to_0p65_dlogs_0p05_logq_minus6_to_minus0_dlogq_0p2_logrho_minus3p7_to_minus1p6_dlogrho_0p3_layer_15_boxsize_10'
#map_set_name = 'kb221684_logs_0_to_0p9_dlogs_0p05_logq_0_to_plus4_dlogq_0p2_logrho_minus3p7_to_minus1p6_dlogrho_0p3_layer_15_boxsize_10'
#map_set_name = 'kb221684_logs_minus1p5_to_0p65_dlogs_0p05_logq_minus6_to_minus0_dlogq_0p2_logrho_minus1p3_to_minus0p7_dlogrho_0p3_layer_15_boxsize_10'
#map_set_name = 'VBBLPython_logs_minus1p5_to_0p65_dlogs_0p05_logq_minus6_to_minus0_dlogq_0p2_logrho_minus3p7_to_minus1p9_dlogrho_0p3_layer_16_boxsize_3p5'
#map_set_name = 'VBBLPython_logs_minus1p5_to_1p3_dlogs_0p05_logq_minus6_to_minus0_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'
#map_set_name = 'VBBLPython_logs_0p05_to_1p3_dlogs_0p05_logq_0p1_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'
#map_set_name = 'VBML5p0Python_logs_minus1p5_to_1p4_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'
# map_set_name = 'VBML5p0Python_test'

# map_set_name = 'VBMicrolensing5p0Python_logs_minus1p5_to_1p5_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'
# map_set_name = 'VBMicrolensing5p0Python_logs_minus1p5_to_1p5_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus3p85_to_minus1p45_dlogrho_0p3_layer_16_boxsize_3p5'
# map_set_name = 'VBMicrolensing5p0Python_logs_minus0p2_to_0p2_dlogs_0p002_logq_minus4p5_to_minus3p0_dlogq_0p05_logrho_minus3p1_to_minus1p9_dlogrho_0p1_layer_16_boxsize_3p5'
# map_set_name = 'VBMicrolensing5p0Python_logs_minus0p2_to_0p2_dlogs_0p002_logq_minus5p5_to_minus3p5_dlogq_0p05_logrho_minus3p7_to_minus2p5_dlogrho_0p3_layer_16_boxsize_3p5'
map_set_name = 'VBMicrolensing5p0Python_logs_minus0p2_to_0p2_dlogs_0p002_logq_minus5p5_to_minus3p0_dlogq_0p1_logrho_minus3p7_to_minus2p5_dlogrho_0p1_layer_16_boxsize_3p5'


def densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, s, q, rho , shift_x , shift_y , threshold_coefficient ):
    
    current_nmesh = 2**current_layer
    current_step_length = 2.*box_size / current_nmesh
    
    mag_interpolation_center = 0.25 * ( current_corner_mag[0] + current_corner_mag[1] + current_corner_mag[2] + current_corner_mag[3] )
    mag_interpolation_up     =  0.5 * ( current_corner_mag[2] + current_corner_mag[3] )
    mag_interpolation_down   =  0.5 * ( current_corner_mag[0] + current_corner_mag[1] )
    mag_interpolation_right  =  0.5 * ( current_corner_mag[1] + current_corner_mag[3] )
    mag_interpolation_left   =  0.5 * ( current_corner_mag[0] + current_corner_mag[2] )
    
    wrong_judgement_threshold = 200000.

    if abs(mag_interpolation_center) > wrong_judgement_threshold or \
       abs(mag_interpolation_up) > wrong_judgement_threshold or \
       abs(mag_interpolation_down) > wrong_judgement_threshold or \
       abs(mag_interpolation_right) > wrong_judgement_threshold or \
       abs(mag_interpolation_left) > wrong_judgement_threshold :
        return [] , [] , False , False 
    #the second False means this map is wrong due to the mistake of VBBL
    
    arg_y = current_serial_number // current_nmesh
    arg_x = current_serial_number % current_nmesh

    y = ( arg_y + 0.5 ) * current_step_length - box_size + shift_y
    x = ( arg_x + 0.5 ) * current_step_length - box_size + shift_x
 
    next_nmesh = current_nmesh * 2
    next_step_length = 0.5 * current_step_length

    
    #mag_vbbl_center = bilens.vbbl_magnification( x , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_up     = bilens.vbbl_magnification( x , y + next_step_length , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_down   = bilens.vbbl_magnification( x , y - next_step_length , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_right  = bilens.vbbl_magnification( x + next_step_length , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_left   = bilens.vbbl_magnification( x - next_step_length , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    mag_vbbl_center = VBML.BinaryMag2( s, q, x , y ,                    rho )
    mag_vbbl_up     = VBML.BinaryMag2( s, q, x , y + next_step_length , rho )
    mag_vbbl_down   = VBML.BinaryMag2( s, q, x , y - next_step_length , rho )
    mag_vbbl_right  = VBML.BinaryMag2( s, q, x + next_step_length , y , rho )
    mag_vbbl_left   = VBML.BinaryMag2( s, q, x - next_step_length , y , rho )


    #if np.abs( mag_interpolation - mag_vbbl_center ) > threshold_coefficient * pow( mag_vbbl_center , 0.5 ) :
    if np.abs( mag_interpolation_center - mag_vbbl_center ) > threshold_coefficient * pow( mag_vbbl_center , 0.5 ) or \
       np.abs( mag_interpolation_up - mag_vbbl_up ) > threshold_coefficient * pow( mag_vbbl_up , 0.5 ) or \
       np.abs( mag_interpolation_down - mag_vbbl_down ) > threshold_coefficient * pow( mag_vbbl_down , 0.5 ) or \
       np.abs( mag_interpolation_right - mag_vbbl_right ) > threshold_coefficient * pow( mag_vbbl_right , 0.5 ) or \
       np.abs( mag_interpolation_left - mag_vbbl_left ) > threshold_coefficient * pow( mag_vbbl_left , 0.5 ) or \
       current_layer <= 5:

        next_serial_number_0 = 2 * arg_y * next_nmesh + 2 * arg_x
        next_serial_number_1 = 2 * arg_y * next_nmesh + ( 2 * arg_x + 1 )
        next_serial_number_2 = ( 2 * arg_y + 1 ) * next_nmesh + 2 * arg_x
        next_serial_number_3 = ( 2 * arg_y + 1 ) * next_nmesh + ( 2 * arg_x + 1 )

        #mag_vbbl_up    = bilens.vbbl_magnification( x + shift_x , y + next_step_length + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_down  = bilens.vbbl_magnification( x + shift_x , y - next_step_length + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_right = bilens.vbbl_magnification( x + next_step_length + shift_x , y + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_left  = bilens.vbbl_magnification( x - next_step_length + shift_x , y + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
    
        next_corner_mag_0 = [ current_corner_mag[0] , mag_vbbl_down , mag_vbbl_left , mag_vbbl_center ]
        next_corner_mag_1 = [ mag_vbbl_down , current_corner_mag[1] , mag_vbbl_center , mag_vbbl_right ]
        next_corner_mag_2 = [ mag_vbbl_left , mag_vbbl_center , current_corner_mag[2] , mag_vbbl_up ]
        next_corner_mag_3 = [ mag_vbbl_center , mag_vbbl_right , mag_vbbl_up , current_corner_mag[3] ]

        ### use the following to avoid that previous 15 layers all correct, while the last layer has some errors, 
        ### which won't be found because there won't be checks on whether last layer needed to be densed to next layer. 
        if abs(mag_vbbl_center) > wrong_judgement_threshold or \
           abs(mag_vbbl_up)     > wrong_judgement_threshold or \
           abs(mag_vbbl_down)   > wrong_judgement_threshold or \
           abs(mag_vbbl_right)  > wrong_judgement_threshold or \
           abs(mag_vbbl_left)   > wrong_judgement_threshold :
            return [] , [] , False , False 
        #the second False means this map is wrong due to the mistake of VBBL

        return [ next_serial_number_0 , next_serial_number_1 , next_serial_number_2 , next_serial_number_3 ] , next_corner_mag_0 + next_corner_mag_1 + next_corner_mag_2 + next_corner_mag_3 , True , True
    
    else:

        return [] , [] , False , True

def generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, s, q, rho , shift_x , shift_y ):
    
    next_serial_number_sum = []
    next_corner_mag_sum = []

    current_whether_densed_sum = []
    current_sequence_number_in_next_layer_file_sum = []

    current_sequence_number_in_next_layer_file = 0

    for i in range( len(current_serial_number_sum) ):
        
        current_serial_number = current_serial_number_sum[i]
        current_corner_mag = [ current_corner_mag_sum[4*i] , current_corner_mag_sum[4*i+1] , current_corner_mag_sum[4*i+2] , current_corner_mag_sum[4*i+3] ]

        next_serial_number , next_corner_mag , current_whether_densed , whether_correct_map = densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, s, q, rho , shift_x , shift_y , threshold_coefficient = 1e-2 )
        
        if whether_correct_map == False :
            return next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , False
        
        next_serial_number_sum.extend( next_serial_number )
        next_corner_mag_sum.extend( next_corner_mag )

        current_whether_densed_sum.append( current_whether_densed )
        if current_whether_densed :
            current_sequence_number_in_next_layer_file_sum.append( current_sequence_number_in_next_layer_file )
            current_sequence_number_in_next_layer_file += 1
        else :
            current_sequence_number_in_next_layer_file_sum.append( None )
        
    return next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , True

### the following is to draw the caustics###
# def getCaustic(separation,mass_ratio,npt=1000):
#     masses = np.array([1.,mass_ratio])
#     totalMass = sum(masses)
#     masses /= totalMass
#     nlens = len(masses)
#     offset = masses[1]*separation
#     zlens1 = np.complex(-offset,0.)
#     zlens2 = np.complex(separation-offset,0.)
#     zlenses = np.array([zlens1,zlens2])
#     ######
#     f0 = np.zeros(2*nlens+1)*1j
#     gc = np.zeros([nlens,2*nlens])*1j
#     fc = np.zeros([nlens,2*nlens])*1j
#     phis = np.linspace(0,2*np.pi,npt+1)[:npt]
#     zcauList,zcriList = [],[]
#     for phi in phis:
#         f0[0] = zlenses[0]**2
#         f0[1] = -2.*zlenses[0]
#         f0[2] = 1.
#         k = 1
#         for ilens in range(1,nlens):
#             k += 2
#             f0[k+1] = f0[k-1]
#             f0[k] = f0[k-2]-2.*f0[k-1]*zlenses[ilens]
#             for j in range(k-2,0,-1):
#                 f0[j+1] = f0[j-1]-2.*f0[j]*zlenses[ilens]+f0[j+1]*zlenses[ilens]**2
#             f0[1] = -2.*f0[0]*zlenses[ilens]+f0[1]*zlenses[ilens]**2
#             f0[0] = f0[0]*zlenses[ilens]**2
#         for ilens in range(nlens):
#             gc[ilens,2*nlens-1] = f0[2*nlens]
#             for j in range(2*nlens,1,-1):
#                 gc[ilens,j-2] = gc[ilens,j-1]*zlenses[ilens]+f0[j-1]
#             fc[ilens,2*nlens-2] = gc[ilens,2*nlens-1]
#             for j in range(2*nlens-1,1,-1):
#                 fc[ilens,j-2] = fc[ilens,j-1]*zlenses[ilens]+gc[ilens,j-1]
#         hc = np.zeros(2*nlens+1)*1j
#         eiphi = np.exp(1j*phi)
#         hc[2*nlens] = f0[2*nlens]*eiphi
#         hc[2*nlens-1] = f0[2*nlens-1]*eiphi
#         for order in range(2*nlens-1,0,-1):
#             hc[order-1] = f0[order-1]*eiphi
#             secondTerm = 0.
#             for ilens in range(nlens):
#                 secondTerm += masses[ilens]*fc[ilens,order-1]
#             hc[order-1] -= secondTerm
#         orders = 2*nlens+1
#         coeffs = np.zeros(orders)*1j
#         for ith in range(orders):
#             coeffs[ith] = hc[orders-ith-1]
#         zcri = np.roots(coeffs)
#         zcau = np.zeros_like(zcri)*1j
#         for eachRoot in range(len(zcri)):
#             zcau[eachRoot] = zcri[eachRoot]
#             for ilens in range(nlens):
#                 zcau[eachRoot] -= masses[ilens]/(np.conj(zcri[eachRoot])-np.conj(zlenses[ilens]))
#         zcauList.extend(zcau)
#         zcriList.extend(zcri)
#     return np.array(zcauList),np.array(zcriList)
"""
if __name__ == '__main__':

    fig = plt.figure(dpi=180, figsize=(8,8))
    ax = fig.add_subplot(111)
    
    current_layer = 0
    current_serial_number_sum = [0]
     
    current_corner_mag_0 = bilens.vbbl_magnification( -box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_1 = bilens.vbbl_magnification( box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_2 = bilens.vbbl_magnification( -box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_3 = bilens.vbbl_magnification( box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    
    current_corner_mag_sum = [ current_corner_mag_0 , current_corner_mag_1 , current_corner_mag_2 , current_corner_mag_3 ]
    
    while(1):
        
        print('the following is layer = %s'%current_layer)
        #print(current_serial_number_sum)
        print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
        print( 'number of real exists square = %s'%len(current_serial_number_sum) )
        #print(current_corner_mag_sum)

        current_nmesh = 2**current_layer
        current_step_length = 2*box_size / current_nmesh

        for i in current_serial_number_sum :
            arg_y = i // current_nmesh
            arg_x = i % current_nmesh
            
            lower_left_y = ( arg_y ) * current_step_length - box_size
            lower_left_x = ( arg_x ) * current_step_length - box_size

            lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            ax.add_patch(rect)
        
        if current_layer == 10 :
            break

        next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum )

        print( 'whether densed = %s'%current_whether_densed_sum )
        print( 'sequence number in next layer file = %s\n'%current_sequence_number_in_next_layer_file_sum )

        current_layer = current_layer + 1
        current_serial_number_sum = next_serial_number_sum
        current_corner_mag_sum = next_corner_mag_sum

    zcauList,zcriList = getCaustic( s , q )
    #print(zcauList)
    ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
    ### end ###

    plt.axis('equal')

    plt.show()


"""
def generating_a_map(arg):

    #fig = plt.figure(dpi=180, figsize=(8,8))
    #ax = fig.add_subplot(111)
    
    print(arg)
    print(parm[arg])
    logs1,logq1,logrho = parm[arg]
    
    q = 10.**logq1
    s = 10.**logs1
    rho = 10.**logrho
    
    if s <= 1 :
        shift_x = 0.
    else : #s>1
        shift_x = (-1.)*(s-1./s)*q/(1.+q)
    
    shift_y = 0.
    
    VBML = VBMicrolensing.VBMicrolensing()
    VBML.Tol=1e-3
    VBML.RelTol=1e-4
    VBML.a1=0.

    #bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
    
    all_layer_serial_number = []
    all_layer_corner_mag = []
    all_layer_whether_densed = []
    all_layer_sequence_number_in_next_layer_file = []
    layer_length = []

    current_layer = 0
    current_serial_number_sum = [0]
    
    #current_corner_mag_0 = bilens.vbbl_magnification( -box_size + shift_x , -box_size + shift_y , rho, u_limb_darkening=None, accuracy=5e-3)
    #current_corner_mag_1 = bilens.vbbl_magnification( box_size + shift_x , -box_size + shift_y , rho, u_limb_darkening=None, accuracy=5e-3)
    #current_corner_mag_2 = bilens.vbbl_magnification( -box_size + shift_x , box_size + shift_y , rho, u_limb_darkening=None, accuracy=5e-3)
    #current_corner_mag_3 = bilens.vbbl_magnification( box_size + shift_x , box_size + shift_y , rho, u_limb_darkening=None, accuracy=5e-3)
    current_corner_mag_0 = VBML.BinaryMag2( s, q, -box_size + shift_x , -box_size + shift_y , rho)
    current_corner_mag_1 = VBML.BinaryMag2( s, q, box_size + shift_x , -box_size + shift_y , rho)
    current_corner_mag_2 = VBML.BinaryMag2( s, q, -box_size + shift_x , box_size + shift_y , rho)
    current_corner_mag_3 = VBML.BinaryMag2( s, q, box_size + shift_x , box_size + shift_y , rho)
    
    current_corner_mag_sum = [ current_corner_mag_0 , current_corner_mag_1 , current_corner_mag_2 , current_corner_mag_3 ]
    
    all_layer_serial_number.extend(current_serial_number_sum)
    all_layer_corner_mag.extend(current_corner_mag_sum)
    layer_length.append(len(current_serial_number_sum))

    while(1):
        
        #print('the following is layer = %s'%current_layer)
        ##print(current_serial_number_sum)
        #print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
        #print( 'number of real exists square = %s'%len(current_serial_number_sum) )
        #print(max(current_corner_mag_sum))
        """
        current_nmesh = 2**current_layer
        current_step_length = 2*box_size / current_nmesh
        
        for i in current_serial_number_sum :
            arg_y = i // current_nmesh
            arg_x = i % current_nmesh
            
            lower_left_y = ( arg_y ) * current_step_length - box_size
            lower_left_x = ( arg_x ) * current_step_length - box_size

            lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            ax.add_patch(rect)
        """
        if current_layer == 16 : #15
            current_whether_densed_sum_final = []
            current_sequence_number_in_next_layer_file_sum_final = []
            for i in range(len(current_serial_number_sum)):
                current_whether_densed_sum_final.append(False)
                current_sequence_number_in_next_layer_file_sum_final.append(None)
            all_layer_whether_densed.extend(current_whether_densed_sum_final)
            all_layer_sequence_number_in_next_layer_file.extend(current_sequence_number_in_next_layer_file_sum_final)

            break

        next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , whether_correct_map_2 = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, s, q, rho, shift_x, shift_y )
        
        if whether_correct_map_2 == False :
            print('map %s is wrong'%parm[arg])
            return

        all_layer_serial_number.extend(next_serial_number_sum)
        all_layer_corner_mag.extend(next_corner_mag_sum)
        all_layer_whether_densed.extend(current_whether_densed_sum)
        all_layer_sequence_number_in_next_layer_file.extend(current_sequence_number_in_next_layer_file_sum)
        layer_length.append(len(next_serial_number_sum))

        #print( 'whether densed = %s'%current_whether_densed_sum )
        #print( 'sequence number in next layer file = %s\n'%current_sequence_number_in_next_layer_file_sum )
        
        current_layer = current_layer + 1
        current_serial_number_sum = next_serial_number_sum
        current_corner_mag_sum = next_corner_mag_sum
    """
    zcauList,zcriList = getCaustic( s , q )
    #print(zcauList)
    ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
    ### end ###

    plt.axis('equal')

    plt.show()
    """
    #print(all_layer_serial_number)
    #print(all_layer_corner_mag)
    #print(all_layer_whether_densed)
    #print(all_layer_sequence_number_in_next_layer_file)
    box_size_array = np.array([box_size])
    all_layer_serial_number_array = np.array(all_layer_serial_number)
    all_layer_corner_mag_array = np.array(np.float32(all_layer_corner_mag))
    all_layer_whether_densed_array = np.array(all_layer_whether_densed)
    all_layer_sequence_number_in_next_layer_file_array = np.array(all_layer_sequence_number_in_next_layer_file)
    layer_length_array = np.array(layer_length)

    np.savez(dir_path + '/%s'%(arg), all_layer_serial_number=all_layer_serial_number_array , \
                                     all_layer_corner_mag=all_layer_corner_mag_array , \
                                     all_layer_whether_densed=all_layer_whether_densed_array , \
                                     all_layer_sequence_number_in_next_layer_file=all_layer_sequence_number_in_next_layer_file_array , \
                                     layer_length=layer_length_array , \
                                     box_size=box_size_array)
    

def saveparm():
    # logsstep = 0.05
    logsstep = 0.002
    #logsstep = 0.0003
    #logsstep = 0.01
    #logsstart,logsend = -0.033,-0.033
    #logsstart,logsend = -1.5, 0.65

    # logsstart,logsend = -1.5, 1.5
    #logsstart,logsend = 0.05, 1.3

    #logsstart,logsend = 0.002, 0.25
    #logsstart,logsend = -0.25, 0.
    #logsstart,logsend = -0.015, 0.015
    #logsstart,logsend = 0.0, 0.9
    #logsstart,logsend = -0.65, 0.0
    #logsstart,logsend = -0.035,-0.03
    logsstart,logsend = -0.2,0.2
    
    logqstep = 0.1
    #logqstep = 0.2
    # logqstep = 0.05
    #logqstart,logqend = -3.842,-3.842
    
    #logqstart,logqend = -6.,0.
    #logqstart,logqend = 0.1, 4.0
    # logqstart,logqend = -6.0, 4.0

    logqstart,logqend = -5.5,-3.0
    #logqstart,logqend = -5.,-4.4
    #logqstart,logqend = -6.,-2.
    #logqstart,logqend = 0.0,+4.0

    # logrhostep = 0.3
    logrhostep = 0.1
    #logrhostart,logrhoend = -2.558,-2.558
    
    # logrhostart,logrhoend = -4.0,-1.6
    # logrhostart,logrhoend = -3.85,-1.45

    # logrhostart,logrhoend = -3.1,-1.9
    #logrhostart,logrhoend = -4.0,-2.8
    logrhostart,logrhoend = -3.7,-2.5
    #logrhostart,logrhoend = -2.077641 ,-2.077641 
    #logrhostart,logrhoend = -2.7,-2.6

    # alphas = np.linspace(0,180,19)
    allparm = []
    logs1 = logsstart
    n=0
    while(logs1<=logsend):
        logq1 = logqstart
        while(logq1<=logqend):
            logrho = logrhostart
            while(logrho<=logrhoend):
                #if 45*logs1+13*logq1 > -58.5 and 65*logs1+19*logq1 < -49.4 and abs(logrho-(-2.5))>1e-3:
                #if abs(logrho-(-2.5))>1e-3:
                #    allparm.append([logs1,logq1,logrho])
                #    n += 1

                # if logs1 <= 0.65 :
                #     if -4*abs(logs1)+logq1+8 > 0 :
                #         allparm.append([logs1,logq1,logrho])
                #         n += 1
                # else : # logs1 >= 0.70
                #     if -4*abs(logs1)+logq1+7 > 0 :
                #         allparm.append([logs1,logq1,logrho])
                #         n += 1
                #if -4*abs(logs1)-logq1+7 > 0 and (10**logs1)**2 > (1+(10**(-logq1))**(1/3))**3 / (1+10**(-logq1) ):
                if (   logq1>0   and  (10**logs1)**2 > (1+(10**(-logq1))**(1/3))**3 / (1+10**(-logq1) )   ) or\
                   (   logq1<=0  and  logs1 <= 0.65  and  (  -4*abs(logs1)+logq1+8 > 0  )  ) or \
                   (   logq1<=0  and  logs1 >= 0.70  and  (  -4*abs(logs1)+logq1+7 > 0  )  ) :

                    allparm.append([logs1,logq1,logrho])
                    n += 1

                #allparm.append([logs1,logq1,logrho])
                #n += 1
                logrho = round(logrho+logrhostep,5)
            logq1 = round(logq1+logqstep,5)
        logs1 = round(logs1+logsstep,5)
    allparm = np.array(allparm)
                
                
    #np.save('parms_for_lc_contrast.npy',allparm)
    np.save('parms_%s.npy'%map_set_name,allparm)
    return n
def loadparm(parm,n):
    #tempparm = np.load('parms_for_lc_contrast.npy')
    tempparm = np.load('parms_%s.npy'%map_set_name)
    for i in range(n):
        for j in range(3):
            parm[i][j] = float(tempparm[i,j])

if __name__ == '__main__':

    ### directory estabilish ###
    dir_path = '/work/zhangjiyuan/adaptive_map_binary_test/map_set_%s'%(map_set_name)

    whether_exist = os.path.exists(dir_path)
    if not whether_exist :
        os.makedirs(dir_path)
        print('%s estabilished successfully'%dir_path)
    else :
        print('%s already exists, overwrite the files inside'%dir_path)

    n = saveparm()
    print ('number of map = %d'%n)

    parm = mp.Array(ctypes.c_float, n*3)
    parm = np.ctypeslib.as_array(parm.get_obj())
    parm = parm.reshape(n,3)
    loadparm(parm,n)

#    start = time.clock()
    pool = mp.Pool(processes=96)
    args = range(len(parm))
    pool.map(generating_a_map,args)
    pool.close()
    pool.join()



#    print ("Total time: %f"%(time.clock()-start))

