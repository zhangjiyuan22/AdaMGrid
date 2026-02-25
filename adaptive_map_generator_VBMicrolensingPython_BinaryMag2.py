import numpy as np
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

box_size = 3.5#0.6#3.5#0.4#10.0#0.8

map_set_name = 'VBMicrolensing5p0Python_logs_minus1p5_to_1p5_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'



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



def generating_a_map(arg):
    
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
    
    all_layer_serial_number = []
    all_layer_corner_mag = []
    all_layer_whether_densed = []
    all_layer_sequence_number_in_next_layer_file = []
    layer_length = []

    current_layer = 0
    current_serial_number_sum = [0]
    
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
    
    logsstep = 0.05
    logsstart,logsend = -1.5, 1.5
    
    logqstep = 0.1
    logqstart,logqend = -6.0, 4.0

    logrhostep = 0.3
    logrhostart,logrhoend = -4.0,-1.6
    
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
    dir_path = './map_set_%s'%(map_set_name)

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

