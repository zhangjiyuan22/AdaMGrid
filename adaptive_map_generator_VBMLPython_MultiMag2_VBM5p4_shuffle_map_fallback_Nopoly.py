import numpy as np
#from MulensModel.binarylens import BinaryLens
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpathes

import multiprocessing
import multiprocessing.pool
import time
import sys
import ctypes
# A return only returns inside that child process.
# It does not return a value to the parent process.
# So you need some IPC mechanism
import queue

#import VBBinaryLensing
import VBMicrolensing
import os

import fcntl
import socket
import traceback
from datetime import datetime

### request all mem might cause hard to schedule
#SBATCH --mem 1536G

## known also by subprocess, even initiated with spawn
JOB_ID = int(os.environ.get("JOB_ID", "0"))        # current job index: 0,1,2,...,49
NUM_JOBS = int(os.environ.get("NUM_JOBS", "1"))    # total number of jobs
SHUFFLE_SEED = 20260402                            # fixed seed so all jobs use the same shuffle
RUN_MODE = os.environ.get("RUN_MODE", "maps")   # "prepare" or "maps"

runtime_log_path = None
slow_log_path = None
# active_dir_path = None

MAP_TIMEOUT_SEC = 1800         # hard timeout per map = 30 min
SLOW_MAP_THRESHOLD_SEC = 1800  # log to slow file if elapsed >= 30 min
HOSTNAME = socket.gethostname()
CHILD_START_METHOD = "spawn"   # safer isolation than fork for the per-map child

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
# map_set_name = 'VBML5p0Python_logs_minus1p5_to_1p4_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'
# map_set_name = 'VBML5p0Python_triple_test_layer_16_boxsize_3p5_apply_async'

# map_set_name = 'VBML5p4Python_triple_test_layer_16_boxsize_3p5_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map'
# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs'
# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly'

# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_test1' # original
# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_test2' # layer16
# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_test3' # reltol 1e-3
# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_test4' # threshold=0.02

# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_OB190468' # original

# map_set_name = 'VBML5p4Python_triple_test_layer_15_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_DynamicAssign_new_cluster_test_shuffle_map_no_large_logs_subprocess_TimeLimit_fallback_Nopoly_fullset_test1' # original
map_set_name = 'VBML5p4Python_triple_test_layer_16_boxsize_3p5_RelTol_5eminus4_MultiBatch_plus_imap_unordered_shuffle_map_subprocess_TimeLimit_fallback_Nopoly_fullset_final' # final

# already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array = np.loadtxt('/storage/maoshudeLab/zhangjy225001/map_set_VBML5p0Python_triple_test_MPI_Plus_Multiprocessing_NoTimeLimit_layer_16_boxsize_3p5/already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit.txt', dtype=int)
# already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_set   = set(already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array)
# # search in set is O(1) through hashing, while search in array is O(n)
# already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array_sorted = np.sort(already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array)



def safe_str(x):
    s = str(x)
    s = s.replace('\t', ' ')
    s = s.replace('\n', ' | ')
    return s


def format_walltime(ts):
    return datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')


def append_tsv_row(log_file, header, row):
    with open(log_file, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        if f.tell() == 0:
            f.write('\t'.join(header) + '\n')
        f.write('\t'.join(safe_str(x) for x in row) + '\n')
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)


# def write_active_file(active_file, info_dict):
#     with open(active_file, 'w') as f:
#         for k, v in info_dict.items():
#             f.write(f'{k}\t{safe_str(v)}\n')


def remove_file_if_exists(path):
    if path is not None and os.path.exists(path):
        os.remove(path)


def calculate_and_cache_mag(
    VBML,
    rho,
    x,
    y,
    key,
    mag_cache,
):
    # """
    # For right/up points:
    # always calculate with VBML, then cache the result.
    # """
    mag = VBML.MultiMag2(x, y, rho)
    mag_cache[key] = mag

    return mag


def get_or_calculate_mag(
    VBML,
    rho,
    x,
    y,
    key,
    mag_cache,
):
    # """
    # For left/down points:
    # use cached value if available.
    # Otherwise calculate directly, but do not cache it.
    # """

    # mag = mag_cache.get(key) # hash, O(1), 100 ns per call
    mag = mag_cache.pop(key, None) # hash, O(1), 100 ns per call

    if mag is None:
        mag = VBML.MultiMag2(x, y, rho)

    return mag


# def densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, s, q, rho , shift_x , shift_y , threshold_coefficient ):
def densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, rho , shift_x , shift_y , threshold_coefficient , mag_cache ):
    
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


    # ============================================================
    # Integer coordinates on the half-grid.
    #
    # Physical position corresponding to (gx, gy):
    #
    # x = gx * next_step_length - box_size + shift_x
    # y = gy * next_step_length - box_size + shift_y
    # ============================================================
    # (gx,gy) is the current-layer-square's center point's index in the next layer odd grid
    # for example, if the current layer is 2*2 grid, then the next layer odd grid runs 0,1,2,3,4 on the bottom edge
    # and for the current layer's left bottom cell, it's gx=1; for the current layer's right bottom cell, it's gx=3
    gy = 2 * arg_y + 1
    gx = 2 * arg_x + 1

    probe_width = 2 * current_nmesh + 1


    #mag_vbbl_center = bilens.vbbl_magnification( x , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_up     = bilens.vbbl_magnification( x , y + next_step_length , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_down   = bilens.vbbl_magnification( x , y - next_step_length , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_right  = bilens.vbbl_magnification( x + next_step_length , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    #mag_vbbl_left   = bilens.vbbl_magnification( x - next_step_length , y , rho , u_limb_darkening=None , accuracy=5e-3 )
    # mag_vbbl_center = VBML.BinaryMag2( s, q, x , y ,                    rho )
    # mag_vbbl_up     = VBML.BinaryMag2( s, q, x , y + next_step_length , rho )
    # mag_vbbl_down   = VBML.BinaryMag2( s, q, x , y - next_step_length , rho )
    # mag_vbbl_right  = VBML.BinaryMag2( s, q, x + next_step_length , y , rho )
    # mag_vbbl_left   = VBML.BinaryMag2( s, q, x - next_step_length , y , rho )

    ### center:       direct calculate, don't cache
    mag_vbbl_center = VBML.MultiMag2( x                    , y                    , rho )

    ### right and up: direct calculate, cache
    # this is the right-point's index in the next-layer-odd-grid
    key_right = gy * probe_width + (gx + 1) 
    mag_vbbl_right = calculate_and_cache_mag(VBML , rho , x + next_step_length , y                    , key_right , mag_cache)
    
    key_up    = (gy + 1) * probe_width + gx
    mag_vbbl_up    = calculate_and_cache_mag(VBML , rho , x                    , y + next_step_length , key_up    , mag_cache)

    ### left and down: check cache, if not cached, then calculate and don't cache
    key_left  = gy * probe_width + (gx - 1)
    mag_vbbl_left  =    get_or_calculate_mag(VBML , rho , x - next_step_length , y                    , key_left  , mag_cache)

    key_down  = (gy - 1) * probe_width + gx
    mag_vbbl_down  =    get_or_calculate_mag(VBML , rho , x                    , y - next_step_length , key_down  , mag_cache)

    # mag_vbbl_left   = VBML.MultiMag2( x - next_step_length , y                    , rho )
    # mag_vbbl_down   = VBML.MultiMag2( x                    , y - next_step_length , rho )
    # mag_vbbl_right  = VBML.MultiMag2( x + next_step_length , y                    , rho )
    # mag_vbbl_up     = VBML.MultiMag2( x                    , y + next_step_length , rho )

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

# def generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, s, q, rho , shift_x , shift_y ):
def generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, rho , shift_x , shift_y ):
      
    next_serial_number_sum = []
    next_corner_mag_sum = []

    current_whether_densed_sum = []
    current_sequence_number_in_next_layer_file_sum = []

    current_sequence_number_in_next_layer_file = 0

    # one cache for this layer
    mag_cache = {}

    for i in range( len(current_serial_number_sum) ):
        
        current_serial_number = current_serial_number_sum[i]
        current_corner_mag = [ current_corner_mag_sum[4*i] , current_corner_mag_sum[4*i+1] , current_corner_mag_sum[4*i+2] , current_corner_mag_sum[4*i+3] ]

        # next_serial_number , next_corner_mag , current_whether_densed , whether_correct_map = densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, s, q, rho , shift_x , shift_y , threshold_coefficient = 1e-2 )
        next_serial_number , next_corner_mag , current_whether_densed , whether_correct_map \
        = densing_a_square( current_layer , current_serial_number , current_corner_mag , VBML, rho , shift_x , shift_y , \
                            threshold_coefficient = 1e-2 , mag_cache=mag_cache ) #0.02
        
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

# """
# if __name__ == '__main__':

#     fig = plt.figure(dpi=180, figsize=(8,8))
#     ax = fig.add_subplot(111)
    
#     current_layer = 0
#     current_serial_number_sum = [0]
     
#     current_corner_mag_0 = bilens.vbbl_magnification( -box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
#     current_corner_mag_1 = bilens.vbbl_magnification( box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
#     current_corner_mag_2 = bilens.vbbl_magnification( -box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
#     current_corner_mag_3 = bilens.vbbl_magnification( box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    
#     current_corner_mag_sum = [ current_corner_mag_0 , current_corner_mag_1 , current_corner_mag_2 , current_corner_mag_3 ]
    
#     while(1):
        
#         print('the following is layer = %s'%current_layer)
#         #print(current_serial_number_sum)
#         print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
#         print( 'number of real exists square = %s'%len(current_serial_number_sum) )
#         #print(current_corner_mag_sum)

#         current_nmesh = 2**current_layer
#         current_step_length = 2*box_size / current_nmesh

#         for i in current_serial_number_sum :
#             arg_y = i // current_nmesh
#             arg_x = i % current_nmesh
            
#             lower_left_y = ( arg_y ) * current_step_length - box_size
#             lower_left_x = ( arg_x ) * current_step_length - box_size

#             lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
#             rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
#             ax.add_patch(rect)
        
#         if current_layer == 10 :
#             break

#         next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum )

#         print( 'whether densed = %s'%current_whether_densed_sum )
#         print( 'sequence number in next layer file = %s\n'%current_sequence_number_in_next_layer_file_sum )

#         current_layer = current_layer + 1
#         current_serial_number_sum = next_serial_number_sum
#         current_corner_mag_sum = next_corner_mag_sum

#     zcauList,zcriList = getCaustic( s , q )
#     #print(zcauList)
#     ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
#     ### end ###

#     plt.axis('equal')

#     plt.show()


# """




# def generating_a_map(arg):
#     with multiprocessing.Pool(processes=1) as sub_pool :
#         try:
#             # pool.apply_async: can only run one function with one set of argument; 
#             #                   non-blocking
#             #                   support multi-arguments, need to be in the form of tuple: (args1, args2, ...)
#             #
#             # pool.map:         can run one function with different arguments;      
#             #                   blocking
#             #                   only support single-argument function
#             future = sub_pool.apply_async(generating_a_map_in_sub_process, (arg,))

#             # we can run some non-blocking code here in main process
#             # ...

#             # below is blocking waiting; main process wait here till sub process return/timeout
#             future.get(timeout=2000) # if sub process doesn't return in 2000 s, then treat it as error

#         except TimeoutError:
#             print(f"Timeout Error in Subprocess in map {arg}")
#             #return None,None,None,None, 1

#         except Exception as e:
#             print(f"Critical Error in Subprocess: {e} in map {arg}")
#             #return None,None,None,None, 1



# def generating_a_map_in_sub_process(arg):
# def generating_a_map(arg):
def generating_a_map_child(arg, params, output_file, result_queue, solver_algorithm_name):
    """
    Heavy work runs here in one fresh subprocess per map.
    This function should NOT touch the shared runtime log files directly.
    The parent/controller will do the final logging.
    """

    # output_file = os.path.join(dir_path, '%s.npz' % arg)

    pid = os.getpid()
    # active_file = os.path.join(active_dir_path, f'pid_{pid}.txt')

    # params = np.array(parm[arg]).copy()
    params = np.array(params, dtype=float)
    logs2, logq2, logs3, logq3, psi, logrho = params

    start_ts = time.time()
    status = 'started'
    error_msg = ''

    # t_geom_sec = -1.0
    # t_build_sec = -1.0
    # t_save_sec = -1.0

    # write_active_file(active_file, {
    #     'job_id': JOB_ID,
    #     'pid': pid,
    #     'hostname': HOSTNAME,
    #     'map_index': arg,
    #     'logs2': logs2,
    #     'logq2': logq2,
    #     'logs3': logs3,
    #     'logq3': logq3,
    #     'psi': psi,
    #     'logrho': logrho,
    #     'status': 'running',
    #     'start_time': format_walltime(start_ts),
    #     'elapsed_sec': 0.0,
    #     'current_layer': -1,
    # })

    try: ### not designed to deal with C++ collapse, but the python collapse

        #fig = plt.figure(dpi=180, figsize=(8,8))
        #ax = fig.add_subplot(111)
        
        # print(arg, already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array_sorted[arg])
        # print(parm[arg])
        # current_map_index = already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_array_sorted[arg]
        # logs2,logq2,logs3,logq3,psi,logrho = parm[arg]

        # print(arg, flush=True)
        # print(parm[arg], flush=True)
        # logs2,logq2,logs3,logq3,psi,logrho = parm[arg]

        # print(f'[START] map={arg} pid={pid} params={params}', flush=True)
        print(f'[CHILD_START] map={arg} pid={pid} method={solver_algorithm_name} params={params}', flush=True)
        
        # s2 = 10.**logs2
        # q2 = 10.**logq2
        # s3 = 10.**logs3
        # q3 = 10.**logq3
        # psi = psi * np.pi/180
        rho = 10.**logrho
        
        ### should set the origin at the mass center of 3 lenses or to say magnification center of 3 lenses, 
        ### but the formula of magnification center of 3 lenses should be explored
        shift_x = 0.
        # if s <= 1 :
        #     shift_x = 0.
        # else : #s>1
        #     shift_x = (-1.)*(s-1./s)*q/(1.+q)
        
        shift_y = 0.
        
        # t0_geom = time.time()

        VBML = VBMicrolensing.VBMicrolensing()

        # VBML.SetMethod(VBML.Method.Multipoly)
        if solver_algorithm_name == "Multipoly":
            VBML.SetMethod(VBML.Method.Multipoly)
        elif solver_algorithm_name == "Nopoly":
            VBML.SetMethod(VBML.Method.Nopoly)
        else:
            raise ValueError(f'Unknown solver_algorithm_name: {solver_algorithm_name}')


        VBML.Tol=1e-3
        VBML.RelTol=5e-4 # 1e-3 #
        VBML.a1=0

        len1_x = 0
        len1_y = 0
        len1_mass = 1

        len2_x = 10**logs2
        len2_y = 0
        len2_mass = 10**logq2

        len3_x = 10**logs3 * np.cos(psi /180*np.pi)
        len3_y = 10**logs3 * np.sin(psi /180*np.pi)
        len3_mass = 10**logq3

        parameters = [len1_x,len1_y,len1_mass,
                    len2_x,len2_y,len2_mass,
                    len3_x,len3_y,len3_mass
                    ]

        VBML.SetLensGeometry(parameters) #Initialize the lens configuration

        #bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
        
        # t_geom_sec = time.time() - t0_geom

        # t0_build = time.time()

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
        # current_corner_mag_0 = VBML.BinaryMag2( s, q, -box_size + shift_x , -box_size + shift_y , rho)
        # current_corner_mag_1 = VBML.BinaryMag2( s, q, box_size + shift_x , -box_size + shift_y , rho)
        # current_corner_mag_2 = VBML.BinaryMag2( s, q, -box_size + shift_x , box_size + shift_y , rho)
        # current_corner_mag_3 = VBML.BinaryMag2( s, q, box_size + shift_x , box_size + shift_y , rho)
        current_corner_mag_0 = VBML.MultiMag2(-box_size + shift_x , -box_size + shift_y , rho)
        current_corner_mag_1 = VBML.MultiMag2( box_size + shift_x , -box_size + shift_y , rho)
        current_corner_mag_2 = VBML.MultiMag2(-box_size + shift_x ,  box_size + shift_y , rho)
        current_corner_mag_3 = VBML.MultiMag2( box_size + shift_x ,  box_size + shift_y , rho)
        
        current_corner_mag_sum = [ current_corner_mag_0 , current_corner_mag_1 , current_corner_mag_2 , current_corner_mag_3 ]
        
        all_layer_serial_number.extend(current_serial_number_sum)
        all_layer_corner_mag.extend(current_corner_mag_sum)
        layer_length.append(len(current_serial_number_sum))

        while(1):

            # write_active_file(active_file, {
            #     'job_id': JOB_ID,
            #     'pid': pid,
            #     'hostname': HOSTNAME,
            #     'map_index': arg,
            #     'logs2': logs2,
            #     'logq2': logq2,
            #     'logs3': logs3,
            #     'logq3': logq3,
            #     'psi': psi,
            #     'logrho': logrho,
            #     'status': 'running',
            #     'start_time': format_walltime(start_ts),
            #     'elapsed_sec': round(time.time() - start_ts, 3),
            #     'current_layer': current_layer,
            # })

            #print('the following is layer = %s'%current_layer)
            ##print(current_serial_number_sum)
            #print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
            #print( 'number of real exists square = %s'%len(current_serial_number_sum) )
            #print(max(current_corner_mag_sum))

            # """
            # current_nmesh = 2**current_layer
            # current_step_length = 2*box_size / current_nmesh
            
            # for i in current_serial_number_sum :
            #     arg_y = i // current_nmesh
            #     arg_x = i % current_nmesh
                
            #     lower_left_y = ( arg_y ) * current_step_length - box_size
            #     lower_left_x = ( arg_x ) * current_step_length - box_size

            #     lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            #     rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            #     ax.add_patch(rect)
            # """
            if current_layer == 16: #16 : #15
                current_whether_densed_sum_final = []
                current_sequence_number_in_next_layer_file_sum_final = []
                for i in range(len(current_serial_number_sum)):
                    current_whether_densed_sum_final.append(False)
                    current_sequence_number_in_next_layer_file_sum_final.append(None)
                all_layer_whether_densed.extend(current_whether_densed_sum_final)
                all_layer_sequence_number_in_next_layer_file.extend(current_sequence_number_in_next_layer_file_sum_final)

                break

            # next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , whether_correct_map_2 = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, s, q, rho, shift_x, shift_y )
            next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , whether_correct_map_2 = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , VBML, rho, shift_x, shift_y )
            
            if whether_correct_map_2 == False:
                status = 'wrong_map'
                error_msg = 'VBML returned > 200000 magnification'
                print(f'[WRONG_MAP] map={arg} method={solver_algorithm_name} params={params}', flush=True)
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
        # """
        # zcauList,zcriList = getCaustic( s , q )
        # #print(zcauList)
        # ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
        # ### end ###

        # plt.axis('equal')

        # plt.show()
        # """

        # t_build_sec = time.time() - t0_build

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

        # t0_save = time.time()

        np.savez_compressed(output_file, all_layer_serial_number=all_layer_serial_number_array , \
                                         all_layer_corner_mag=all_layer_corner_mag_array , \
                                         all_layer_whether_densed=all_layer_whether_densed_array , \
                                         all_layer_sequence_number_in_next_layer_file=all_layer_sequence_number_in_next_layer_file_array , \
                                         layer_length=layer_length_array , \
                                         box_size=box_size_array)
        # np.savez(dir_path + '/%s'%(current_map_index), all_layer_serial_number=all_layer_serial_number_array , \
        #                                 all_layer_corner_mag=all_layer_corner_mag_array , \
        #                                 all_layer_whether_densed=all_layer_whether_densed_array , \
        #                                 all_layer_sequence_number_in_next_layer_file=all_layer_sequence_number_in_next_layer_file_array , \
        #                                 layer_length=layer_length_array , \
        #                                 box_size=box_size_array)

        # t_save_sec = time.time() - t0_save

        status = 'ok'
        # print(f'[DONE] map={arg} pid={pid} elapsed={time.time()-start_ts:.2f}s', flush=True)
        print(f'[CHILD_DONE] map={arg} pid={pid} method={solver_algorithm_name} elapsed={time.time() - start_ts:.2f}s', flush=True)

    except Exception as e:
        status = 'exception'
        error_msg = f'{type(e).__name__}: {e}'
        # print(f'[FAIL] map={arg} pid={pid} error={error_msg}', flush=True)
        print(f'[CHILD_FAIL] map={arg} pid={pid} method={solver_algorithm_name} error={error_msg}', flush=True)
        traceback.print_exc()

    finally: ### either try or except will conduct finally
        end_ts = time.time()
        elapsed_sec = end_ts - start_ts
        size_bytes = os.path.getsize(output_file) if os.path.exists(output_file) else -1

        # header = [
        #     'job_id', 'pid', 'hostname', 'map_index',
        #     'logs2', 'logq2', 'logs3', 'logq3', 'psi', 'logrho',
        #     'status', 'error_msg',
        #     'start_ts', 'end_ts', 'elapsed_sec',
        #     'start_time', 'end_time',
        #     't_geom_sec', 't_build_sec', 't_save_sec',
        #     'size_bytes'
        # ]

        # row = [
        #     JOB_ID, pid, HOSTNAME, arg,
        #     logs2, logq2, logs3, logq3, psi, logrho,
        #     status, error_msg,
        #     f'{start_ts:.6f}', f'{end_ts:.6f}', f'{elapsed_sec:.6f}',
        #     format_walltime(start_ts), format_walltime(end_ts),
        #     f'{t_geom_sec:.6f}', f'{t_build_sec:.6f}', f'{t_save_sec:.6f}',
        #     size_bytes
        # ]

        # append_tsv_row(runtime_log_path, header, row)

        # if elapsed_sec >= SLOW_MAP_THRESHOLD_SEC:
        #     append_tsv_row(slow_log_path, header, row)

        # remove_file_if_exists(active_file)

        ###############################################
        ### below is changing to subprocess version ###
        ###############################################
        # result = {
        #     'job_id': JOB_ID,
        #     'pid': pid,
        #     'hostname': HOSTNAME,
        #     'map_index': int(arg),
        #     'logs2': float(logs2),
        #     'logq2': float(logq2),
        #     'logs3': float(logs3),
        #     'logq3': float(logq3),
        #     'psi': float(psi),
        #     'logrho': float(logrho),
        #     'status': status,
        #     'error_msg': error_msg,
        #     'start_ts': start_ts,
        #     'end_ts': end_ts,
        #     'elapsed_sec': elapsed_sec,
        #     'start_time': format_walltime(start_ts),
        #     'end_time': format_walltime(end_ts),
        #     't_geom_sec': t_geom_sec,
        #     't_build_sec': t_build_sec,
        #     't_save_sec': t_save_sec,
        #     'size_bytes': size_bytes,
        # }
        result = {
            'child_pid': pid,
            'status': status,       ### ok/wrong_map/exception
            'error_msg': error_msg, ### ' '/'>200000'/f'{type(e).__name__}: {e}'
            'child_start_ts': start_ts,
            'child_end_ts': end_ts,
            'child_elapsed_sec': elapsed_sec,
            'size_bytes': size_bytes,
        }

        try:
            result_queue.put_nowait(result)
        except Exception:
            pass



def run_one_map_attempt(arg, params, output_file, solver_algorithm_name):
    """
    Run one child attempt with the requested VBML algorithm.
    Returns a dict describing the outcome of this one attempt.
    """
    # attempt_start_ts = time.time()

    child_status = 'unknown'
    error_msg = ''
    child_exitcode = None
    child_pid = -1

    ### defaults only for result == None
    # t_geom_sec = -1.0
    # t_build_sec = -1.0
    # t_save_sec = -1.0
    child_start_ts = -1.0
    child_end_ts = -1.0
    child_elapsed_sec = -1.0
    size_bytes = -1

    ### low-level manual process control: 
    ###   manage everything yourself: start, wait, timeout, kill, communicate with queue/pipe
    ### while map/imap_unordered are high-level batch interface
    ctx = multiprocessing.get_context(CHILD_START_METHOD)
    # A return only returns inside that child process.
    # It does not return a value to the parent process.
    # So you need some IPC mechanism, like queue, shared mem, or temp_file
    result_queue = ctx.Queue(maxsize=1)

    proc = ctx.Process(
        target=generating_a_map_child,
        args=(int(arg), params, output_file, result_queue, solver_algorithm_name)
    )

    proc.start()
    child_pid = proc.pid if proc.pid is not None else -1

    proc.join(MAP_TIMEOUT_SEC) # wait up to MAP_TIMEOUT_SEC; but return earlier if the child already finished

    if proc.is_alive():                                        # timeout
        # hard timeout
        child_status = 'timeout'
        error_msg = f'exceeded {MAP_TIMEOUT_SEC} sec'

        proc.terminate()
        proc.join(5)

        if proc.is_alive():
            proc.kill()
            proc.join() # wait for the process to actually exit and be cleaned up

        child_exitcode = proc.exitcode

        # if killed during np.savez, partial file may exist
        remove_file_if_exists(output_file)

    else:
        child_exitcode = proc.exitcode

        result = None
        try:
            result = result_queue.get_nowait()
        except queue.Empty:
            result = None
        except Exception:
            result = None

        if result is not None:                                 # ok, wrong_map, exception(python collapse by try-except)
            # result = {
            #     'child_pid': pid,
            #     'status': status,       ### ok/wrong_map/exception
            #     'error_msg': error_msg, ### ' '/'>200000'/f'{type(e).__name__}: {e}'
            #     'child_start_ts': start_ts,
            #     'child_end_ts': end_ts,
            #     'child_elapsed_sec': elapsed_sec,
            #     'size_bytes': size_bytes,
            # }
            child_pid    = result.get('child_pid', child_pid)
            child_status = result.get('status', 'unknown')
            error_msg    = result.get('error_msg', '')
            child_start_ts    = result.get('child_start_ts', -1.0)
            child_end_ts      = result.get('child_end_ts', -1.0)
            child_elapsed_sec = result.get('child_elapsed_sec', -1.0)
            # t_geom_sec = result.get('t_geom_sec', -1.0)
            # t_build_sec = result.get('t_build_sec', -1.0)
            # t_save_sec = result.get('t_save_sec', -1.0)
            size_bytes = result.get('size_bytes', -1)
        else:
            # child died before sending result
            if child_exitcode == 0:                            # (very rare)
                child_status = 'unexpected_no_result'
                error_msg = 'child exited cleanly but returned no result'
            else:                                              # crash of C++
                child_status = 'crash'
                error_msg = f'child exitcode={child_exitcode}'

            remove_file_if_exists(output_file)

    # attempt_end_ts = time.time()
    # attempt_elapsed_sec = attempt_end_ts - attempt_start_ts

    try:
        result_queue.close()
        result_queue.join_thread()
    except Exception:
        pass

    return {
        'solver_algorithm_name': solver_algorithm_name,
        ###
        'child_status': child_status,
        'error_msg': error_msg,
        'child_exitcode': child_exitcode,
        'child_pid': child_pid,
        'child_start_ts': child_start_ts,
        'child_end_ts': child_end_ts,
        'child_elapsed_sec': child_elapsed_sec,
        'size_bytes': size_bytes,
        ###
        # 'attempt_start_ts': attempt_start_ts,
        # 'attempt_end_ts': attempt_end_ts,
        # 'attempt_elapsed_sec': attempt_elapsed_sec,
    }



def generating_a_map_parent(arg):
    """
    Lightweight controller.
    This runs inside the long-lived outer Pool worker.
    The heavy VBML work is done in one fresh subprocess per map.

    First try Multipoly.
    If the child crashes natively, retry once with Nopoly.
    """
    output_file = os.path.join(dir_path, f'{arg}.npz')
    params = tuple(float(x) for x in parm[arg])

    logs2, logq2, logs3, logq3, psi, logrho = params

    ## remove stale partial file before starting
    # remove_file_if_exists(output_file)

    parent_start_ts = time.time()

    # child_status = 'unknown'
    # error_msg = ''
    # child_exitcode = None
    # child_pid = -1

    # ### defaults only for result == None
    # # t_geom_sec = -1.0
    # # t_build_sec = -1.0
    # # t_save_sec = -1.0
    # child_start_ts = -1.0
    # child_end_ts = -1.0
    # child_elapsed_sec = -1.0
    # size_bytes = -1

    # ### low-level manual process control: 
    # ###   manage everything yourself: start, wait, timeout, kill, communicate with queue/pipe
    # ### while map/imap_unordered are high-level batch interface
    # ctx = multiprocessing.get_context(CHILD_START_METHOD)
    # # A return only returns inside that child process.
    # # It does not return a value to the parent process.
    # # So you need some IPC mechanism, like queue, shared mem, or temp_file
    # result_queue = ctx.Queue(maxsize=1)

    # proc = ctx.Process(
    #     target=generating_a_map_child,
    #     args=(int(arg), params, output_file, result_queue, solver_algorithm_name)
    # )

    # proc.start()
    # child_pid = proc.pid if proc.pid is not None else -1

    # proc.join(MAP_TIMEOUT_SEC) # wait up to MAP_TIMEOUT_SEC; but return earlier if the child already finished

    # if proc.is_alive():                                        # timeout
    #     # hard timeout
    #     child_status = 'timeout'
    #     error_msg = f'exceeded {MAP_TIMEOUT_SEC} sec'

    #     proc.terminate()
    #     proc.join(5)

    #     if proc.is_alive():
    #         proc.kill()
    #         proc.join() # wait for the process to actually exit and be cleaned up

    #     child_exitcode = proc.exitcode

    #     # if killed during np.savez, partial file may exist
    #     remove_file_if_exists(output_file)

    # else:
    #     child_exitcode = proc.exitcode

    #     result = None
    #     try:
    #         result = result_queue.get_nowait()
    #     except queue.Empty:
    #         result = None
    #     except Exception:
    #         result = None

    #     if result is not None:                                 # ok, wrong_map, exception(python collapse by try-except)
    #         # result = {
    #         #     'child_pid': pid,
    #         #     'status': status,       ### ok/wrong_map/exception
    #         #     'error_msg': error_msg, ### ' '/'>200000'/f'{type(e).__name__}: {e}'
    #         #     'child_start_ts': start_ts,
    #         #     'child_end_ts': end_ts,
    #         #     'child_elapsed_sec': elapsed_sec,
    #         #     'size_bytes': size_bytes,
    #         # }
    #         child_pid    = result.get('child_pid', child_pid)
    #         child_status = result.get('status', 'unknown')
    #         error_msg    = result.get('error_msg', '')
    #         child_start_ts    = result.get('child_start_ts', -1.0)
    #         child_end_ts      = result.get('child_end_ts', -1.0)
    #         child_elapsed_sec = result.get('child_elapsed_sec', -1.0)
    #         # t_geom_sec = result.get('t_geom_sec', -1.0)
    #         # t_build_sec = result.get('t_build_sec', -1.0)
    #         # t_save_sec = result.get('t_save_sec', -1.0)
    #         size_bytes = result.get('size_bytes', -1)
    #     else:
    #         # child died before sending result
    #         if child_exitcode == 0:                            # (very rare)
    #             child_status = 'unexpected_no_result'
    #             error_msg = 'child exited cleanly but returned no result'
    #         else:                                              # crash of C++
    #             child_status = 'crash'
    #             error_msg = f'child exitcode={child_exitcode}'

    #         remove_file_if_exists(output_file)



    # First attempt: Multipoly
    attempt1 = run_one_map_attempt(arg, params, output_file, "Multipoly")

    fallback_used = False
    first_try_status = attempt1['child_status']
    first_try_exitcode = attempt1['child_exitcode']

    # Default final result = first attempt
    final_attempt = attempt1

    # Only retry native crash with Nopoly
    if attempt1['child_status'] == 'crash':
        fallback_used = True
        print(f'[RETRY_NOPOLY] map={arg} first_exitcode={attempt1["child_exitcode"]}', flush=True)

        remove_file_if_exists(output_file)
        attempt2 = run_one_map_attempt(arg, params, output_file, "Nopoly")
        final_attempt = attempt2



    parent_end_ts = time.time()
    parent_elapsed_sec = parent_end_ts - parent_start_ts


    header = [
        'job_id', 'parent_pid', 'child_pid', 'hostname', 'map_index',
        'logs2', 'logq2', 'logs3', 'logq3', 'psi', 'logrho',
        'solver_algorithm_name', 'fallback_used', 'first_try_status', 'first_try_exitcode',
        'child_status', 'child_error_msg', 'child_exitcode',
        'child_start_ts','child_end_ts','child_elapsed_sec',
        'parent_start_ts', 'parent_end_ts', 'parent_elapsed_sec',
        'parent_start_time', 'parent_end_time',
        # 't_geom_sec', 't_build_sec', 't_save_sec',
        'size_bytes'
    ]

    # row = [
    #     JOB_ID, os.getpid(), child_pid, HOSTNAME, int(arg),                            # definitely have
    #     logs2, logq2, logs3, logq3, psi, logrho,                                       # definitely have
    #     child_status, error_msg, child_exitcode,                                       # definitely have
    #     f'{child_start_ts:.6f}', f'{child_end_ts:.6f}', f'{child_elapsed_sec:.6f}',    # only result != None have
    #     f'{parent_start_ts:.6f}', f'{parent_end_ts:.6f}', f'{parent_elapsed_sec:.6f}', # definitely have
    #     format_walltime(parent_start_ts), format_walltime(parent_end_ts),
    #     # f'{t_geom_sec:.6f}', f'{t_build_sec:.6f}', f'{t_save_sec:.6f}',
    #     size_bytes                                                                     # only result != None have
    # ]
    row = [
        JOB_ID, os.getpid(), final_attempt['child_pid'], HOSTNAME, int(arg),                              # definitely have
        logs2, logq2, logs3, logq3, psi, logrho,                                                          # definitely have
        final_attempt['solver_algorithm_name'], int(fallback_used), first_try_status, first_try_exitcode,
        final_attempt['child_status'], final_attempt['error_msg'], final_attempt['child_exitcode'],       # definitely have
        f'{final_attempt["child_start_ts"]:.6f}', f'{final_attempt["child_end_ts"]:.6f}', f'{final_attempt["child_elapsed_sec"]:.6f}', # only result != None have
        f'{parent_start_ts:.6f}', f'{parent_end_ts:.6f}', f'{parent_elapsed_sec:.6f}',                    # definitely have
        format_walltime(parent_start_ts), format_walltime(parent_end_ts),
        # f'{t_geom_sec:.6f}', f'{t_build_sec:.6f}', f'{t_save_sec:.6f}',
        final_attempt['size_bytes']                                                                       # only result != None have
    ]


    append_tsv_row(runtime_log_path, header, row)

    if parent_elapsed_sec >= SLOW_MAP_THRESHOLD_SEC or final_attempt['child_status'] in ('timeout', 'crash', 'unexpected_no_result'):
        append_tsv_row(slow_log_path, header, row)

    # print(
    #     f'[PARENT_DONE] map={arg} status={child_status} child_exitcode={child_exitcode} elapsed={parent_elapsed_sec:.2f}s',
    #     flush=True
    # )
    print(
        f"[PARENT_DONE] map={arg} method={final_attempt['solver_algorithm_name']} "
        f"status={final_attempt['child_status']} child_exitcode={final_attempt['child_exitcode']} "
        f"fallback_used={fallback_used} elapsed={parent_elapsed_sec:.2f}s",
        flush=True
    )



SQ_GRID = np.array([
    [-1.5, 0], [-1.425, 0], [-1.35, 0], [-1.275, 0], [-1.2, 0], [-1.125, 0], [-1.05, 0], [-0.975, 0], [-0.9, 0], [-0.825, 0], [-0.75, 0], [-0.675, 0], [-0.6, 0], [-0.525, 0], [-0.45, 0], [-0.375, 0], [-0.3, 0], [-0.225, 0], [-0.15, 0], [-0.075, 0], [0, 0], [0.075, 0], [0.15, 0], [0.225, 0], [0.3, 0], [0.375, 0], [0.45, 0], [0.525, 0], [0.6, 0], [0.675, 0], [0.75, 0], [0.825, 0], [0.9, 0], [0.975, 0], [1.05, 0], [1.125, 0], [1.2, 0], [1.275, 0], [1.35, 0], [1.425, 0],
    [-1.4751271108, -0.315789473684], [-1.4001271108, -0.315789473684], [-1.3251271108, -0.315789473684], [-1.2501271108, -0.315789473684], [-1.1751271108, -0.315789473684], [-1.1001271108, -0.315789473684], [-1.0251271108, -0.315789473684], [-0.950127110799, -0.315789473684], [-0.875127110799, -0.315789473684], [-0.800127110799, -0.315789473684], [-0.725127110799, -0.315789473684], [-0.650127110799, -0.315789473684], [-0.575127110799, -0.315789473684], [-0.500127110799, -0.315789473684], [-0.432338741437, -0.315789473684], [-0.36028228453, -0.315789473684], [-0.288225827624, -0.315789473684], [-0.216169370718, -0.315789473684], [-0.144112913812, -0.315789473684], [-0.0720564569061, -0.315789473684], [0, -0.315789473684], [0.0720564569061, -0.315789473684], [0.144112913812, -0.315789473684], [0.216169370718, -0.315789473684], [0.288225827624, -0.315789473684], [0.36028228453, -0.315789473684], [0.432338741437, -0.315789473684], [0.504395198343, -0.315789473684], [0.583663285887, -0.315789473684], [0.658663285887, -0.315789473684], [0.733663285887, -0.315789473684], [0.808663285887, -0.315789473684], [0.883663285887, -0.315789473684], [0.958663285887, -0.315789473684], [1.03366328589, -0.315789473684], [1.10866328589, -0.315789473684], [1.18366328589, -0.315789473684], [1.25866328589, -0.315789473684], [1.33366328589, -0.315789473684], [1.40866328589, -0.315789473684],
    [-1.44019446679, -0.631578947368], [-1.36519446679, -0.631578947368], [-1.29019446679, -0.631578947368], [-1.21519446679, -0.631578947368], [-1.14019446679, -0.631578947368], [-1.06519446679, -0.631578947368], [-0.990194466795, -0.631578947368], [-0.915194466795, -0.631578947368], [-0.840194466795, -0.631578947368], [-0.765194466795, -0.631578947368], [-0.690194466795, -0.631578947368], [-0.615194466795, -0.631578947368], [-0.540194466795, -0.631578947368], [-0.465194466795, -0.631578947368], [-0.399081739264, -0.631578947368], [-0.332568116053, -0.631578947368], [-0.266054492842, -0.631578947368], [-0.199540869632, -0.631578947368], [-0.133027246421, -0.631578947368], [-0.0665136232106, -0.631578947368], [0, -0.631578947368], [0.0665136232106, -0.631578947368], [0.133027246421, -0.631578947368], [0.199540869632, -0.631578947368], [0.266054492842, -0.631578947368], [0.332568116053, -0.631578947368], [0.399081739264, -0.631578947368], [0.465595362474, -0.631578947368], [0.540996258154, -0.631578947368], [0.615996258154, -0.631578947368], [0.690996258154, -0.631578947368], [0.765996258154, -0.631578947368], [0.840996258154, -0.631578947368], [0.915996258154, -0.631578947368], [0.990996258154, -0.631578947368], [1.06599625815, -0.631578947368], [1.14099625815, -0.631578947368], [1.21599625815, -0.631578947368], [1.29099625815, -0.631578947368], [1.36599625815, -0.631578947368],
    [-1.31867775177, -0.947368421053], [-1.24367775177, -0.947368421053], [-1.16867775177, -0.947368421053], [-1.09367775177, -0.947368421053], [-1.01867775177, -0.947368421053], [-0.943677751768, -0.947368421053], [-0.868677751768, -0.947368421053], [-0.793677751768, -0.947368421053], [-0.718677751768, -0.947368421053], [-0.643677751768, -0.947368421053], [-0.568677751768, -0.947368421053], [-0.493677751768, -0.947368421053], [-0.418677751768, -0.947368421053], [-0.354702084334, -0.947368421053], [-0.295585070278, -0.947368421053], [-0.236468056222, -0.947368421053], [-0.177351042167, -0.947368421053], [-0.118234028111, -0.947368421053], [-0.0591170140556, -0.947368421053], [0, -0.947368421053], [0.0591170140556, -0.947368421053], [0.118234028111, -0.947368421053], [0.177351042167, -0.947368421053], [0.236468056222, -0.947368421053], [0.295585070278, -0.947368421053], [0.354702084334, -0.947368421053], [0.413819098389, -0.947368421053], [0.48396044501, -0.947368421053], [0.55896044501, -0.947368421053], [0.63396044501, -0.947368421053], [0.70896044501, -0.947368421053], [0.78396044501, -0.947368421053], [0.85896044501, -0.947368421053], [0.93396044501, -0.947368421053], [1.00896044501, -0.947368421053], [1.08396044501, -0.947368421053], [1.15896044501, -0.947368421053], [1.23396044501, -0.947368421053], [1.30896044501, -0.947368421053], [1.38396044501, -0.947368421053],
    [-1.27014749174, -1.26315789474], [-1.19514749174, -1.26315789474], [-1.12014749174, -1.26315789474], [-1.04514749174, -1.26315789474], [-0.970147491744, -1.26315789474], [-0.895147491744, -1.26315789474], [-0.820147491744, -1.26315789474], [-0.745147491744, -1.26315789474], [-0.670147491744, -1.26315789474], [-0.595147491744, -1.26315789474], [-0.520147491744, -1.26315789474], [-0.445147491744, -1.26315789474], [-0.370147491744, -1.26315789474], [-0.308239340885, -1.26315789474], [-0.256866117404, -1.26315789474], [-0.205492893923, -1.26315789474], [-0.154119670443, -1.26315789474], [-0.102746446962, -1.26315789474], [-0.0513732234808, -1.26315789474], [0, -1.26315789474], [0.0513732234808, -1.26315789474], [0.102746446962, -1.26315789474], [0.154119670443, -1.26315789474], [0.205492893923, -1.26315789474], [0.256866117404, -1.26315789474], [0.308239340885, -1.26315789474], [0.359612564366, -1.26315789474], [0.424077636988, -1.26315789474], [0.499077636988, -1.26315789474], [0.574077636988, -1.26315789474], [0.649077636988, -1.26315789474], [0.724077636988, -1.26315789474], [0.799077636988, -1.26315789474], [0.874077636988, -1.26315789474], [0.949077636988, -1.26315789474], [1.02407763699, -1.26315789474], [1.09907763699, -1.26315789474], [1.17407763699, -1.26315789474], [1.24907763699, -1.26315789474], [1.32407763699, -1.26315789474],
    [-1.21709654603, -1.57894736842], [-1.14209654603, -1.57894736842], [-1.06709654603, -1.57894736842], [-0.992096546031, -1.57894736842], [-0.917096546031, -1.57894736842], [-0.842096546031, -1.57894736842], [-0.767096546031, -1.57894736842], [-0.692096546031, -1.57894736842], [-0.617096546031, -1.57894736842], [-0.542096546031, -1.57894736842], [-0.467096546031, -1.57894736842], [-0.392096546031, -1.57894736842], [-0.321410493639, -1.57894736842], [-0.265488419134, -1.57894736842], [-0.221240349278, -1.57894736842], [-0.176992279422, -1.57894736842], [-0.132744209567, -1.57894736842], [-0.0884961397112, -1.57894736842], [-0.0442480698556, -1.57894736842], [0, -1.57894736842], [0.0442480698556, -1.57894736842], [0.0884961397112, -1.57894736842], [0.132744209567, -1.57894736842], [0.176992279422, -1.57894736842], [0.221240349278, -1.57894736842], [0.265488419134, -1.57894736842], [0.309736488989, -1.57894736842], [0.364434589123, -1.57894736842], [0.435120641515, -1.57894736842], [0.510120641515, -1.57894736842], [0.585120641515, -1.57894736842], [0.660120641515, -1.57894736842], [0.735120641515, -1.57894736842], [0.810120641515, -1.57894736842], [0.885120641515, -1.57894736842], [0.960120641515, -1.57894736842], [1.03512064151, -1.57894736842], [1.11012064151, -1.57894736842], [1.18512064151, -1.57894736842], [1.26012064151, -1.57894736842],
    [-1.17060074425, -1.89473684211], [-1.09560074425, -1.89473684211], [-1.02060074425, -1.89473684211], [-0.945600744251, -1.89473684211], [-0.870600744251, -1.89473684211], [-0.795600744251, -1.89473684211], [-0.720600744251, -1.89473684211], [-0.645600744251, -1.89473684211], [-0.570600744251, -1.89473684211], [-0.495600744251, -1.89473684211], [-0.420600744251, -1.89473684211], [-0.345600744251, -1.89473684211], [-0.279449084512, -1.89473684211], [-0.229213277911, -1.89473684211], [-0.191011064926, -1.89473684211], [-0.152808851941, -1.89473684211], [-0.114606638956, -1.89473684211], [-0.0764044259704, -1.89473684211], [-0.0382022129852, -1.89473684211], [0, -1.89473684211], [0.0382022129852, -1.89473684211], [0.0764044259704, -1.89473684211], [0.114606638956, -1.89473684211], [0.152808851941, -1.89473684211], [0.191011064926, -1.89473684211], [0.229213277911, -1.89473684211], [0.267415490896, -1.89473684211], [0.312685216758, -1.89473684211], [0.378836876497, -1.89473684211], [0.453836876497, -1.89473684211], [0.528836876497, -1.89473684211], [0.603836876497, -1.89473684211], [0.678836876497, -1.89473684211], [0.753836876497, -1.89473684211], [0.828836876497, -1.89473684211], [0.903836876497, -1.89473684211], [0.978836876497, -1.89473684211], [1.0538368765, -1.89473684211], [1.1288368765, -1.89473684211], [1.2038368765, -1.89473684211],
    [-1.13366632883, -2.21052631579], [-1.05866632883, -2.21052631579], [-0.98366632883, -2.21052631579], [-0.90866632883, -2.21052631579], [-0.83366632883, -2.21052631579], [-0.75866632883, -2.21052631579], [-0.68366632883, -2.21052631579], [-0.60866632883, -2.21052631579], [-0.53366632883, -2.21052631579], [-0.45866632883, -2.21052631579], [-0.38366632883, -2.21052631579], [-0.30866632883, -2.21052631579], [-0.246143012686, -2.21052631579], [-0.200186529154, -2.21052631579], [-0.166822107628, -2.21052631579], [-0.133457686102, -2.21052631579], [-0.100093264577, -2.21052631579], [-0.0667288430512, -2.21052631579], [-0.0333644215256, -2.21052631579], [0, -2.21052631579], [0.0333644215256, -2.21052631579], [0.0667288430512, -2.21052631579], [0.100093264577, -2.21052631579], [0.133457686102, -2.21052631579], [0.166822107628, -2.21052631579], [0.200186529154, -2.21052631579], [0.233550950679, -2.21052631579], [0.27100552096, -2.21052631579], [0.333528837105, -2.21052631579], [0.408528837105, -2.21052631579], [0.483528837105, -2.21052631579], [0.558528837105, -2.21052631579], [0.633528837105, -2.21052631579], [0.708528837105, -2.21052631579], [0.783528837105, -2.21052631579], [0.858528837105, -2.21052631579], [0.933528837105, -2.21052631579], [1.0085288371, -2.21052631579], [1.0835288371, -2.21052631579], [1.1585288371, -2.21052631579],
    [-1.09751173265, -2.52631578947], [-1.02251173265, -2.52631578947], [-0.947511732646, -2.52631578947], [-0.872511732646, -2.52631578947], [-0.797511732646, -2.52631578947], [-0.722511732646, -2.52631578947], [-0.647511732646, -2.52631578947], [-0.572511732646, -2.52631578947], [-0.497511732646, -2.52631578947], [-0.422511732646, -2.52631578947], [-0.347511732646, -2.52631578947], [-0.27662693124, -2.52631578947], [-0.220985595562, -2.52631578947], [-0.178052274167, -2.52631578947], [-0.148376895139, -2.52631578947], [-0.118701516111, -2.52631578947], [-0.0890261370836, -2.52631578947], [-0.0593507580557, -2.52631578947], [-0.0296753790279, -2.52631578947], [0, -2.52631578947], [0.0296753790279, -2.52631578947], [0.0593507580557, -2.52631578947], [0.0890261370836, -2.52631578947], [0.118701516111, -2.52631578947], [0.148376895139, -2.52631578947], [0.178052274167, -2.52631578947], [0.207727653195, -2.52631578947], [0.238982779369, -2.52631578947], [0.294624115047, -2.52631578947], [0.365508916453, -2.52631578947], [0.440508916453, -2.52631578947], [0.515508916453, -2.52631578947], [0.590508916453, -2.52631578947], [0.665508916453, -2.52631578947], [0.740508916453, -2.52631578947], [0.815508916453, -2.52631578947], [0.890508916453, -2.52631578947], [0.965508916453, -2.52631578947], [1.04050891645, -2.52631578947], [1.11550891645, -2.52631578947],
    [-1.071225028, -2.84210526316], [-0.996225027997, -2.84210526316], [-0.921225027997, -2.84210526316], [-0.846225027997, -2.84210526316], [-0.771225027997, -2.84210526316], [-0.696225027997, -2.84210526316], [-0.621225027997, -2.84210526316], [-0.546225027997, -2.84210526316], [-0.471225027997, -2.84210526316], [-0.396225027997, -2.84210526316], [-0.321225027997, -2.84210526316], [-0.253373036103, -2.84210526316], [-0.202786382945, -2.84210526316], [-0.161877290104, -2.84210526316], [-0.134897741754, -2.84210526316], [-0.107918193403, -2.84210526316], [-0.0809386450522, -2.84210526316], [-0.0539590967014, -2.84210526316], [-0.0269795483507, -2.84210526316], [0, -2.84210526316], [0.0269795483507, -2.84210526316], [0.0539590967014, -2.84210526316], [0.0809386450522, -2.84210526316], [0.107918193403, -2.84210526316], [0.134897741754, -2.84210526316], [0.161877290104, -2.84210526316], [0.188856838455, -2.84210526316], [0.215396616491, -2.84210526316], [0.265983269649, -2.84210526316], [0.333835261543, -2.84210526316], [0.408835261543, -2.84210526316], [0.483835261543, -2.84210526316], [0.558835261543, -2.84210526316], [0.633835261543, -2.84210526316], [0.708835261543, -2.84210526316], [0.783835261543, -2.84210526316], [0.858835261543, -2.84210526316], [0.933835261543, -2.84210526316], [1.00883526154, -2.84210526316], [1.08383526154, -2.84210526316],
    [-1.04658508489, -3.15789473684], [-0.971585084892, -3.15789473684], [-0.896585084892, -3.15789473684], [-0.821585084892, -3.15789473684], [-0.746585084892, -3.15789473684], [-0.671585084892, -3.15789473684], [-0.596585084892, -3.15789473684], [-0.521585084892, -3.15789473684], [-0.446585084892, -3.15789473684], [-0.371585084892, -3.15789473684], [-0.296585084892, -3.15789473684], [-0.232747940767, -3.15789473684], [-0.188852700559, -3.15789473684], [-0.163875563443, -3.15789473684], [-0.140464768665, -3.15789473684], [-0.117053973888, -3.15789473684], [-0.0936431791103, -3.15789473684], [-0.0702323843327, -3.15789473684], [-0.0468215895551, -3.15789473684], [-0.0234107947776, -3.15789473684], [0, -3.15789473684], [0.0234107947776, -3.15789473684], [0.0468215895551, -3.15789473684], [0.0702323843327, -3.15789473684], [0.0936431791103, -3.15789473684], [0.117053973888, -3.15789473684], [0.140464768665, -3.15789473684], [0.163875563443, -3.15789473684], [0.197425413271, -3.15789473684], [0.241320653479, -3.15789473684], [0.305157797604, -3.15789473684], [0.380157797604, -3.15789473684], [0.455157797604, -3.15789473684], [0.530157797604, -3.15789473684], [0.605157797604, -3.15789473684], [0.680157797604, -3.15789473684], [0.755157797604, -3.15789473684], [0.830157797604, -3.15789473684], [0.905157797604, -3.15789473684], [0.980157797604, -3.15789473684],
    [-1.03441874014, -3.47368421053], [-0.959418740143, -3.47368421053], [-0.884418740143, -3.47368421053], [-0.809418740143, -3.47368421053], [-0.734418740143, -3.47368421053], [-0.659418740143, -3.47368421053], [-0.584418740143, -3.47368421053], [-0.509418740143, -3.47368421053], [-0.434418740143, -3.47368421053], [-0.359418740143, -3.47368421053], [-0.28444665869, -3.47368421053], [-0.221993189603, -3.47368421053], [-0.180357543544, -3.47368421053], [-0.155439745284, -3.47368421053], [-0.133234067387, -3.47368421053], [-0.111028389489, -3.47368421053], [-0.088822711591, -3.47368421053], [-0.0666170336933, -3.47368421053], [-0.0444113557955, -3.47368421053], [-0.0222056778978, -3.47368421053], [0, -3.47368421053], [0.0222056778978, -3.47368421053], [0.0444113557955, -3.47368421053], [0.0666170336933, -3.47368421053], [0.088822711591, -3.47368421053], [0.111028389489, -3.47368421053], [0.133234067387, -3.47368421053], [0.155439745284, -3.47368421053], [0.186036141769, -3.47368421053], [0.227671787827, -3.47368421053], [0.290125256914, -3.47368421053], [0.365097338367, -3.47368421053], [0.440097338367, -3.47368421053], [0.515097338367, -3.47368421053], [0.590097338367, -3.47368421053], [0.665097338367, -3.47368421053], [0.740097338367, -3.47368421053], [0.815097338367, -3.47368421053], [0.890097338367, -3.47368421053], [0.965097338367, -3.47368421053],
    [-1.02373428323, -3.78947368421], [-0.948734283229, -3.78947368421], [-0.873734283229, -3.78947368421], [-0.798734283229, -3.78947368421], [-0.723734283229, -3.78947368421], [-0.648734283229, -3.78947368421], [-0.573734283229, -3.78947368421], [-0.498734283229, -3.78947368421], [-0.423734283229, -3.78947368421], [-0.348734283229, -3.78947368421], [-0.275110073248, -3.78947368421], [-0.214903056612, -3.78947368421], [-0.174765045521, -3.78947368421], [-0.149848574738, -3.78947368421], [-0.12844163549, -3.78947368421], [-0.107034696242, -3.78947368421], [-0.0856277569933, -3.78947368421], [-0.064220817745, -3.78947368421], [-0.0428138784966, -3.78947368421], [-0.0214069392483, -3.78947368421], [0, -3.78947368421], [0.0214069392483, -3.78947368421], [0.0428138784966, -3.78947368421], [0.064220817745, -3.78947368421], [0.0856277569933, -3.78947368421], [0.107034696242, -3.78947368421], [0.12844163549, -3.78947368421], [0.149848574738, -3.78947368421], [0.178449452076, -3.78947368421], [0.218587463166, -3.78947368421], [0.278794479802, -3.78947368421], [0.352418689784, -3.78947368421], [0.427418689784, -3.78947368421], [0.502418689784, -3.78947368421], [0.577418689784, -3.78947368421], [0.652418689784, -3.78947368421], [0.727418689784, -3.78947368421], [0.802418689784, -3.78947368421], [0.877418689784, -3.78947368421], [0.952418689784, -3.78947368421],
    [-0.941816443372, -4.10526315789], [-0.866816443372, -4.10526315789], [-0.791816443372, -4.10526315789], [-0.716816443372, -4.10526315789], [-0.641816443372, -4.10526315789], [-0.566816443372, -4.10526315789], [-0.491816443372, -4.10526315789], [-0.416816443372, -4.10526315789], [-0.341816443372, -4.10526315789], [-0.269067222275, -4.10526315789], [-0.210318520446, -4.10526315789], [-0.171152719226, -4.10526315789], [-0.146218991219, -4.10526315789], [-0.125330563902, -4.10526315789], [-0.104442136585, -4.10526315789], [-0.0835537092679, -4.10526315789], [-0.062665281951, -4.10526315789], [-0.041776854634, -4.10526315789], [-0.020888427317, -4.10526315789], [0, -4.10526315789], [0.020888427317, -4.10526315789], [0.041776854634, -4.10526315789], [0.062665281951, -4.10526315789], [0.0835537092679, -4.10526315789], [0.104442136585, -4.10526315789], [0.125330563902, -4.10526315789], [0.146218991219, -4.10526315789], [0.173506331504, -4.10526315789], [0.212672132723, -4.10526315789], [0.271420834552, -4.10526315789], [0.34417005565, -4.10526315789], [0.41917005565, -4.10526315789], [0.49417005565, -4.10526315789], [0.56917005565, -4.10526315789], [0.64417005565, -4.10526315789], [0.71917005565, -4.10526315789], [0.79417005565, -4.10526315789], [0.86917005565, -4.10526315789], [0.94417005565, -4.10526315789],
    [-0.862402167125, -4.42105263158], [-0.787402167125, -4.42105263158], [-0.712402167125, -4.42105263158], [-0.637402167125, -4.42105263158], [-0.562402167125, -4.42105263158], [-0.487402167125, -4.42105263158], [-0.412402167125, -4.42105263158], [-0.337402167125, -4.42105263158], [-0.265212295657, -4.42105263158], [-0.207395843211, -4.42105263158], [-0.16885154158, -4.42105263158], [-0.143898726088, -4.42105263158], [-0.123341765219, -4.42105263158], [-0.102784804349, -4.42105263158], [-0.0822278434791, -4.42105263158], [-0.0616708826093, -4.42105263158], [-0.0411139217395, -4.42105263158], [-0.0205569608698, -4.42105263158], [0, -4.42105263158], [0.0205569608698, -4.42105263158], [0.0411139217395, -4.42105263158], [0.0616708826093, -4.42105263158], [0.0822278434791, -4.42105263158], [0.102784804349, -4.42105263158], [0.123341765219, -4.42105263158], [0.143898726088, -4.42105263158], [0.170338312771, -4.42105263158], [0.208882614402, -4.42105263158], [0.266699066848, -4.42105263158], [0.338888938316, -4.42105263158], [0.413888938316, -4.42105263158], [0.488888938316, -4.42105263158], [0.563888938316, -4.42105263158], [0.638888938316, -4.42105263158], [0.713888938316, -4.42105263158], [0.788888938316, -4.42105263158], [0.863888938316, -4.42105263158],
    [-0.784614556267, -4.73684210526], [-0.709614556267, -4.73684210526], [-0.634614556267, -4.73684210526], [-0.559614556267, -4.73684210526], [-0.484614556267, -4.73684210526], [-0.409614556267, -4.73684210526], [-0.334614556267, -4.73684210526], [-0.262778345057, -4.73684210526], [-0.205551326372, -4.73684210526], [-0.167399980583, -4.73684210526], [-0.142431690948, -4.73684210526], [-0.122084306527, -4.73684210526], [-0.101736922106, -4.73684210526], [-0.0813895376844, -4.73684210526], [-0.0610421532633, -4.73684210526], [-0.0406947688422, -4.73684210526], [-0.0203473844211, -4.73684210526], [0, -4.73684210526], [0.0203473844211, -4.73684210526], [0.0406947688422, -4.73684210526], [0.0610421532633, -4.73684210526], [0.0813895376844, -4.73684210526], [0.101736922106, -4.73684210526], [0.122084306527, -4.73684210526], [0.142431690948, -4.73684210526], [0.168331862366, -4.73684210526], [0.206483208155, -4.73684210526], [0.26371022684, -4.73684210526], [0.33554643805, -4.73684210526], [0.41054643805, -4.73684210526], [0.48554643805, -4.73684210526], [0.56054643805, -4.73684210526], [0.63554643805, -4.73684210526], [0.71054643805, -4.73684210526], [0.78554643805, -4.73684210526],
    [-0.707866982161, -5.05263157895], [-0.632866982161, -5.05263157895], [-0.557866982161, -5.05263157895], [-0.482866982161, -5.05263157895], [-0.407866982161, -5.05263157895], [-0.332866982161, -5.05263157895], [-0.261252659608, -5.05263157895], [-0.204395455352, -5.05263157895], [-0.166490652515, -5.05263157895], [-0.141511263925, -5.05263157895], [-0.121295369078, -5.05263157895], [-0.101079474232, -5.05263157895], [-0.0808635793855, -5.05263157895], [-0.0606476845391, -5.05263157895], [-0.0404317896928, -5.05263157895], [-0.0202158948464, -5.05263157895], [0, -5.05263157895], [0.0202158948464, -5.05263157895], [0.0404317896928, -5.05263157895], [0.0606476845391, -5.05263157895], [0.0808635793855, -5.05263157895], [0.101079474232, -5.05263157895], [0.121295369078, -5.05263157895], [0.141511263925, -5.05263157895], [0.16707161245, -5.05263157895], [0.204976415287, -5.05263157895], [0.261833619542, -5.05263157895], [0.333447942096, -5.05263157895], [0.408447942096, -5.05263157895], [0.483447942096, -5.05263157895], [0.558447942096, -5.05263157895], [0.633447942096, -5.05263157895], [0.708447942096, -5.05263157895],
    [-0.631776972888, -5.36842105263], [-0.556776972888, -5.36842105263], [-0.481776972888, -5.36842105263], [-0.406776972888, -5.36842105263], [-0.331776972888, -5.36842105263], [-0.260301118676, -5.36842105263], [-0.203674694989, -5.36842105263], [-0.165923745865, -5.36842105263], [-0.140936876731, -5.36842105263], [-0.120803037198, -5.36842105263], [-0.100669197665, -5.36842105263], [-0.080535358132, -5.36842105263], [-0.060401518599, -5.36842105263], [-0.040267679066, -5.36842105263], [-0.020133839533, -5.36842105263], [0, -5.36842105263], [0.020133839533, -5.36842105263], [0.040267679066, -5.36842105263], [0.060401518599, -5.36842105263], [0.080535358132, -5.36842105263], [0.100669197665, -5.36842105263], [0.120803037198, -5.36842105263], [0.140936876731, -5.36842105263], [0.166284606429, -5.36842105263], [0.204035555554, -5.36842105263], [0.26066197924, -5.36842105263], [0.332137833452, -5.36842105263], [0.407137833452, -5.36842105263], [0.482137833452, -5.36842105263], [0.557137833452, -5.36842105263], [0.632137833452, -5.36842105263],
    [-0.556099512083, -5.68421052632], [-0.481099512083, -5.68421052632], [-0.406099512083, -5.68421052632], [-0.331099512083, -5.68421052632], [-0.25970974612, -5.68421052632], [-0.203226802847, -5.68421052632], [-0.165571507332, -5.68421052632], [-0.140579769922, -5.68421052632], [-0.120496945647, -5.68421052632], [-0.100414121373, -5.68421052632], [-0.0803312970983, -5.68421052632], [-0.0602484728237, -5.68421052632], [-0.0401656485491, -5.68421052632], [-0.0200828242746, -5.68421052632], [0, -5.68421052632], [0.0200828242746, -5.68421052632], [0.0401656485491, -5.68421052632], [0.0602484728237, -5.68421052632], [0.0803312970983, -5.68421052632], [0.100414121373, -5.68421052632], [0.120496945647, -5.68421052632], [0.140579769922, -5.68421052632], [0.165795093198, -5.68421052632], [0.203450388713, -5.68421052632], [0.259933331985, -5.68421052632], [0.331323097948, -5.68421052632], [0.406323097948, -5.68421052632], [0.481323097948, -5.68421052632], [0.556323097948, -5.68421052632],
    [-0.4806795043, -6], [-0.4056795043, -6], [-0.3306795043, -6], [-0.259343121503, -6], [-0.202949150174, -6], [-0.165353169289, -6], [-0.140358328639, -6], [-0.120307138834, -6], [-0.100255949028, -6], [-0.0802047592224, -6], [-0.0601535694168, -6], [-0.0401023796112, -6], [-0.0200511898056, -6], [0, -6], [0.0200511898056, -6], [0.0401023796112, -6], [0.0601535694168, -6], [0.0802047592224, -6], [0.100255949028, -6], [0.120307138834, -6], [0.140358328639, -6], [0.165491462504, -6], [0.203087443389, -6], [0.259481414718, -6], [0.330817797515, -6], [0.405817797515, -6], [0.480817797515, -6],
], dtype=float)

def saveparm(parm_file):

    psi_grid = np.arange(
        0,
        180 + 1e-10,
        3
    )

    # logrho_grid = np.arange(
    #     -4.0,
    #     -2.0 + 1e-10,
    #     0.5
    # )
    logrho_grid = np.arange(
        -3.6,
        -2.0 + 1e-10,
        0.4
    )

    allparm = []
    n = 0

    # for logs2, logq2 in SQ_GRID:

    #     for logs3, logq3 in SQ_GRID:

    #         # remove q2 < q3 exchange duplication
    #         if (logq2+1e-5) < logq3:
    #             continue
    for i_unordered_pairs, (logs2, logq2) in enumerate(SQ_GRID):
        
        for logs3, logq3 in SQ_GRID[i_unordered_pairs:]:

            for psi in psi_grid:

                for logrho in logrho_grid:
                    
                    # if -0.4<logs3<0.4 and -2.5<logq3<-1.5 and -0.4<logs2<0.4 and -3.0<logq2<-2.0 and 40<psi<55 and logrho<-2.9: #OB190468
                    allparm.append([
                        logs2,
                        logq2,
                        logs3,
                        logq3,
                        psi,
                        logrho
                    ])

                    n += 1

    allparm = np.asarray(allparm)

    np.save(parm_file, allparm)

    return n



# def saveparm(parm_file):
    
    

#     # # 21 s2
#     # logs2_step = 0.15
#     # # logs2_start,logs2_end = -1.5, 1.5
#     # logs2_start,logs2_end = -1.35, 0.6
    
#     # # 11 q2
#     # logq2_step = 0.6
#     # logq2_start,logq2_end = -6.0, 0.0

#     # # 21 s3
#     # logs3_step = 0.15
#     # # logs3_start,logs3_end = -1.5, 1.5
#     # logs3_start,logs3_end = -1.35, 0.6
    
#     # # 11 q3
#     # logq3_step = 0.6
#     # logq3_start,logq3_end = -6.0, 0.0

#     # 61 psi
#     psi_step = 3
#     psi_start,psi_end = 0., 180. #degree

#     # 5 rho
#     logrho_step = 0.5
#     logrho_start,logrho_end = -4.0,-2.0

#     ################ 
#     # around 5 million maps
#     ################

#     # alphas = np.linspace(0,180,19)
#     allparm = []
#     logs2 = logs2_start
#     n=0
#     while(logs2 <= logs2_end):
#         logq2 = logq2_start
#         while(logq2 <= logq2_end):
#             logs3 = logs3_start
#             while(logs3 <= logs3_end):
#                 logq3 = logq3_start
#                 while(logq3 <= logq3_end):
#                     psi = psi_start
#                     while(psi <= psi_end):
#                         logrho = logrho_start
#                         while(logrho <= logrho_end):
#                             #if 45*logs1+13*logq1 > -58.5 and 65*logs1+19*logq1 < -49.4 and abs(logrho-(-2.5))>1e-3:
#                             #if abs(logrho-(-2.5))>1e-3:
#                             #    allparm.append([logs1,logq1,logrho])
#                             #    n += 1

#                             # if logs1 <= 0.65 :
#                             #     if -4*abs(logs1)+logq1+8 > 0 :
#                             #         allparm.append([logs1,logq1,logrho])
#                             #         n += 1
#                             # else : # logs1 >= 0.70
#                             #     if -4*abs(logs1)+logq1+7 > 0 :
#                             #         allparm.append([logs1,logq1,logrho])
#                             #         n += 1
#                             #if -4*abs(logs1)-logq1+7 > 0 and (10**logs1)**2 > (1+(10**(-logq1))**(1/3))**3 / (1+10**(-logq1) ):
                            
#                             # if (   logq1>0   and  (10**logs1)**2 > (1+(10**(-logq1))**(1/3))**3 / (1+10**(-logq1) )   ) or\
#                             # (   logq1<=0  and  logs1 <= 0.65  and  (  -4*abs(logs1)+logq1+8 > 0  )  ) or \
#                             # (   logq1<=0  and  logs1 >= 0.70  and  (  -4*abs(logs1)+logq1+7 > 0  )  ) :

#                             #     allparm.append([logs1,logq1,logrho])
#                             #     n += 1

#                             if logq2 >= logq3 and (-4*abs(logs2)+logq2+7 > 0) and (-4*abs(logs3)+logq3+7 > 0) :

#                                 # if n in already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_set: 
#                                 #     allparm.append([logs2,logq2,logs3,logq3,psi,logrho])

#                                 allparm.append([logs2,logq2,logs3,logq3,psi,logrho])
#                                 n += 1

#                             #allparm.append([logs1,logq1,logrho])
#                             #n += 1

#                             logrho = round(logrho+logrho_step,5)
#                         psi = round(psi+psi_step,5)
#                     logq3 = round(logq3+logq3_step,5)
#                 logs3 = round(logs3+logs3_step,5)
#             logq2 = round(logq2+logq2_step,5)
#         logs2 = round(logs2+logs2_step,5)
#     allparm = np.array(allparm)
                
                
#     #np.save('parms_for_lc_contrast.npy',allparm)
#     # np.save('parms_%s.npy'%map_set_name, allparm)
#     np.save(parm_file, allparm)
#     # return n, len(already_generated_map_index_MPI_Plus_Multiprocessing_NoTimeLimit_set)
#     return n

# def loadparm(parm,n, n_already_generated):
def loadparm(parm,n, parm_file):
    #tempparm = np.load('parms_for_lc_contrast.npy')
    tempparm = np.load(parm_file)
    # for i in range(n_already_generated):
    for i in range(n):
        for j in range(6):
            parm[i][j] = float(tempparm[i,j])






class NoDaemonProcess(multiprocessing.Process):
    # Make 'daemon' attribute always False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass  # Ignore attempts to set it

    daemon = property(_get_daemon, _set_daemon)

class NonDaemonPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwargs):
        proc = super().Process(*args, **kwargs)
        proc.__class__ = NoDaemonProcess
        return proc


if __name__ == '__main__':

    ### directory estabilish ###
    # dir_path = '/storage2/maoshudeLab/zhangjy225001/map_set_%s'%(map_set_name)
    # whether_exist = os.path.exists(dir_path)
    # if not whether_exist :
    #     os.makedirs(dir_path)
    #     print('%s estabilished successfully'%dir_path)
    # else :
    #     print('%s already exists, overwrite the files inside'%dir_path)
    base_dir = '/storage2/maoshudeLab/zhangjy225001/map_set_%s' % (map_set_name)
    job_dir = os.path.join(base_dir, 'job_%04d' % JOB_ID)
    assigned_dir = os.path.join(base_dir, 'assigned_indices')

    os.makedirs(base_dir, exist_ok=True)
    # os.makedirs(job_dir, exist_ok=True)
    os.makedirs(assigned_dir, exist_ok=True)

    parm_file = os.path.join(base_dir, 'parms_%s.npy' % map_set_name)

    print(f'RUN_MODE = {RUN_MODE}')
    print(f'base_dir = {base_dir}')
    # print(f'job_dir = {job_dir}')
    print(f'assigned_dir = {assigned_dir}')
    print(f'parm_file = {parm_file}')
    
    # --------------------------------------------------
    # prepare mode: only create parameter file, then exit
    # --------------------------------------------------
    if RUN_MODE == 'prepare':
        if os.path.exists(parm_file):
            tempparm = np.load(parm_file, mmap_mode='r')
            n = len(tempparm)
            print(f'parameter file already exists: {parm_file}')
            print(f'number of map = {n}')
        else:
            n = saveparm(parm_file)
            print(f'parameter file created: {parm_file}')
            print(f'number of map = {n}')
        sys.exit(0)

    # --------------------------------------------------
    # map-generation mode
    # --------------------------------------------------
    if not os.path.exists(parm_file):
        raise FileNotFoundError(
            f'parameter file does not exist: {parm_file}\n'
            f'Please run first with RUN_MODE=prepare'
        )

    os.makedirs(job_dir, exist_ok=True)

    print(f'job_dir = {job_dir}')
    dir_path = job_dir


    ###### for elapse time log #####
    runtime_log_path = os.path.join(job_dir, 'map_runtime_log.tsv')
    slow_log_path = os.path.join(job_dir, 'slow_map_runtime_log.tsv')
    # active_dir_path = os.path.join(job_dir, 'active')

    # os.makedirs(active_dir_path, exist_ok=True)

    print(f'runtime_log_path = {runtime_log_path}')
    print(f'slow_log_path = {slow_log_path}')
    # print(f'active_dir_path = {active_dir_path}')
    ###### for elapse time log #####

    # parm = np.load(parm_file, mmap_mode='r')
    # n = len(parm)

    # print(f'number of map = {n}')
    # print(f'parm dtype = {parm.dtype}')
    # print(f'parm shape = {parm.shape}')

    tempparm = np.load(parm_file, mmap_mode='r')
    n = len(tempparm)
    print(f'number of map = {n}')

    parm = multiprocessing.Array(ctypes.c_float, n * 6)
    parm = np.ctypeslib.as_array(parm.get_obj())
    parm = parm.reshape(n, 6)
    loadparm(parm, n, parm_file)

    
    
    # global shuffle, then split among independent jobs
    all_indices = np.arange(len(parm))
    rng = np.random.default_rng(SHUFFLE_SEED)
    shuffled_indices = rng.permutation(all_indices)

    # args = shuffled_indices[JOB_ID::NUM_JOBS]
    chunks = np.array_split(shuffled_indices, NUM_JOBS)
    args = chunks[JOB_ID]

    print(f"JOB_ID = {JOB_ID}")
    print(f"NUM_JOBS = {NUM_JOBS}")
    print(f"Assigned number of maps = {len(args)}")
    print(f"First 20 assigned indices = {args[:20]}")

    assigned_txt = os.path.join(
        assigned_dir,
        'job_%04d_assigned_indices.txt' % JOB_ID
    )
    np.savetxt(
        assigned_txt,
        args.reshape(-1, 1),
        fmt='%d',
        header='map_index',
        comments=''
    )
    print(f'assigned indices saved to {assigned_txt}')

    # Why the outer pool should be ThreadPool, not multiprocessing.Pool:
    ###########################################################
    ## Pool is daemonic, not allowed to create child process, while ThreadPool allow
    ## Pool are 256 processes, ThreadPool are 256 threads in the same process and don't use 256 cpu cores

    ### A process pool means:
    # multiple separate OS processes
    # each has its own Python interpreter
    # each has its own memory space
    # no shared GIL between processes

    ### A thread pool means:
    # one process
    # many threads inside that one process
    # threads share the same memory space
    # Python threads are affected by the GIL for Python bytecode 

    ### GIL note:(Global Interpreter Lock, only one thread at a time can execute Python bytecode inside one Python process)
    # - Python threads do not scale well for CPU-heavy Python code because of the GIL.
    # - But that is not a problem here, since the outer threads mostly wait, and the
    #   heavy work is in child processes, each with its own interpreter/GIL.

    ### it's not like “this thread does not have its own core, so its child process cannot get a core”
    ### The correct picture is:
    # the thread asks the OS: “please start this new process”
    # the OS starts the process
    # the OS scheduler then places that new process on an available core

    ### my previous NonDaemonPool is non-standard and fragile, and also not better in efficiency

    ### it's not weried that process is not allowed to generate subprocess, while thread is allowed
    # Python makes strong assumptions about process, including daemon-like behavior, because the pool wants tight control.
    # Letting pool workers freely create more child processes *complicates* shutdown, cleanup, and bookkeeping.

    # A thread pool worker is just a thread inside the main process.
    # Threads are not separate OS processes, so spawning a subprocess from a thread is much more normal and common.
    ###########################################################

    ###########################################################
    ### but the ThreadPool cause one process has too many 'file'(including process,queue) opened, so crash
    ### so we have to change back to NonDaemonPool
    ###########################################################


    # with multiprocessing.Pool(processes=256) as pool:
    # with multiprocessing.pool.ThreadPool(processes=256) as pool:
    with NonDaemonPool(processes=256) as pool: 
        # Outer pool = lightweight controllers.
        # Each controller launches one fresh child process per map.
        list(pool.imap_unordered(generating_a_map_parent, args, chunksize=1))
    


# #    start = time.clock()
#     # pool = mp.Pool(processes=64)
#     args = range(len(parm))
#     # pool.map(generating_a_map,args)
#     # pool.close()
#     # pool.join()
#     with multiprocessing.Pool(processes=256) as pool:
#     # with NonDaemonPool(processes=64) as pool: 
#         #args = range(len(parm))
#         #pool.map(generating_a_map, args)
        
#         #results = []
#         #for i in range(len(parm)):
#         #    result = pool.apply_async(generating_a_map, (i,))
#         #    results.append(result)
#         #
#         #output = [result.get() for result in results]

#         list(pool.imap_unordered(generating_a_map, args, chunksize=1))

# #    print ("Total time: %f"%(time.clock()-start))


