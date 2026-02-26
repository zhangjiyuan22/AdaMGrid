import numpy as np
import ctypes
from scipy.optimize import fmin
import time 
import multiprocessing as mp
import emcee
#from mpi4py import MPI
from PyAstronomy import pyasl
#Chi2Lib = ctypes.cdll.LoadLibrary('./chi2_calculator.so')
#Chi2Lib = ctypes.cdll.LoadLibrary('./chi2_calculator_pipeline.so') # C version, not C++
Chi2Lib = ctypes.cdll.LoadLibrary('./chi2_calculator_pipeline_struct.so') # C version, not C++


map_set_name = 'VBMicrolensing5p0Python_logs_minus1p5_to_1p5_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5'


#def chi2(tempparms,data, all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,use_mcmc):
def chi2(tempparms,data, one_map_all_layer_struct, layer_length_accumulated_in_front_of_each_layer,box_size,use_mcmc,s,q):
    initial_type = ctypes.c_float*4
    initial = initial_type()
    initial[0] = tempparms[0]
    initial[1] = tempparms[1]
    initial[2] = tempparms[2]
    initial[3] = tempparms[3]

    chi2 = 0.
    ichi2type = ctypes.c_float
    ichi2 = ichi2type()
    ip = ctypes.byref(ichi2)

    box_size_type = ctypes.c_float
    box_size_final = box_size_type()
    box_size_final = box_size

    s_type = ctypes.c_float
    s_final = s_type()
    s_final = s

    q_type = ctypes.c_float
    q_final = q_type()
    q_final = q
    
    
    #Chi2Lib.wrapgetchi2.argtypes = [ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_float,ctypes.POINTER(ctypes.c_float)]
    Chi2Lib.wrapgetchi2.argtypes = [ctypes.POINTER(ctypes.c_float),\
                                    ctypes.POINTER(ctypes.c_float),\
                                    ctypes.POINTER(ctypes.c_float),\
                                    ctypes.c_int,\
                                    ctypes.POINTER(ctypes.c_float),\
                                    ctypes.POINTER(whether_offset_4mag),\
                                    ctypes.POINTER(ctypes.c_int),\
                                    ctypes.c_float,\
                                    ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_float]
    #Chi2Lib.wrapprintlc.argtypes = [ctypes.POINTER(ctypes.c_float),ctypes.c_int,ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_float]

    #total_time = 0 #ns
    #print ("Total time: %f (ns)"%total_time)
    for i in range(len(data)):

        
        #start = time.time_ns()



        nhjd = len(data[i])
        #temptype = ctypes.c_float*nhjd
        hjds  = each_dataset_hjds[i]
        iflux = each_dataset_fluxs[i]
        iferr = each_dataset_ferrs[i]

        """ 
        start = time.time_ns()
        idata = data[i]
        for i in range(nhjd):
            hjds[i]=idata[i,0]
            iflux[i]=idata[i,1]
            iferr[i]=idata[i,2]
        end   = time.time_ns()
        """

        

        
        #start = time.time_ns()
        #Chi2Lib.wrapgetchi2(hjds,iflux,iferr,nhjd,initial,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size_final,ip)
        Chi2Lib.wrapgetchi2(hjds,iflux,iferr,nhjd,initial, \
                            one_map_all_layer_struct, \
                            layer_length_accumulated_in_front_of_each_layer,\
                            box_size_final, ip, s_final, q_final)
        #end   = time.time_ns()

        
        if ichi2.value <0 :
            chi2 = np.inf
            break
        chi2 += ichi2.value



        #end   = time.time_ns()
        #total_time += ( end - start )

    #print ("Total time: %e (ns)"%total_time)

    #print chi2
    if use_mcmc == True:
        return -0.5*chi2 * hot_mcmc_factor
    else:
        return chi2 * hot_mcmc_factor


#def grid(n):
def grid(p):
    #start = time.time()
    n,logs,logq,logrho = p
    s = 10**logs
    q = 10**logq
    rho = 10**logrho
    print(p)

    try:
        print (n)
        global use_mcmc
        
        map_content = np.load(mapdir+"%s.npz"%n,allow_pickle=True)
        all_layer_corner_mag_raw = map_content['all_layer_corner_mag']
        all_layer_whether_densed_raw = map_content['all_layer_whether_densed']
        all_layer_sequence_number_in_next_layer_file_raw = map_content['all_layer_sequence_number_in_next_layer_file']
        layer_length_raw = map_content['layer_length']
        box_size = map_content['box_size'][0]
        
        nlayer = len(layer_length_raw)

        layer_length_accumulated_in_front_of_each_layer_raw = []
        accumulation = 0
        for i in range(nlayer):
            layer_length_accumulated_in_front_of_each_layer_raw.append(accumulation)
            accumulation += layer_length_raw[i] 

        layer_length_raw = list(map(lambda x:int(x+0.5),layer_length_raw))
        sum_layer_length = sum(layer_length_raw)
        layer_length_accumulated_in_front_of_each_layer_raw = list(map(lambda x:int(x+0.5),layer_length_accumulated_in_front_of_each_layer_raw))

        layer_length_accumulated_type = ctypes.c_int*nlayer
        layer_length_accumulated_in_front_of_each_layer = layer_length_accumulated_type()
        for i in range(nlayer):
            layer_length_accumulated_in_front_of_each_layer[i] = layer_length_accumulated_in_front_of_each_layer_raw[i]


        ### using following too see the ctypes.c_bool is aligned to 4 bytes ###
        #print(whether_offset_4mag.whether_densed)
        #print(whether_offset_4mag.sequence_number_in_next_layer_file)
        #print(whether_offset_4mag.corner_mag_4)
        one_map_all_layer_struct_type = whether_offset_4mag * sum_layer_length
        one_map_all_layer_struct = one_map_all_layer_struct_type()

        for i in range(sum_layer_length):
            one_map_all_layer_struct[i].whether_densed = all_layer_whether_densed_raw[i]
            one_map_all_layer_struct[i].corner_mag_1 = all_layer_corner_mag_raw[ 4 * i + 0 ]
            one_map_all_layer_struct[i].corner_mag_2 = all_layer_corner_mag_raw[ 4 * i + 1 ]
            one_map_all_layer_struct[i].corner_mag_3 = all_layer_corner_mag_raw[ 4 * i + 2 ]
            one_map_all_layer_struct[i].corner_mag_4 = all_layer_corner_mag_raw[ 4 * i + 3 ]
            try:
                one_map_all_layer_struct[i].sequence_number_in_next_layer_file = int(all_layer_sequence_number_in_next_layer_file_raw[i]+0.5) 
            except: 
                one_map_all_layer_struct[i].sequence_number_in_next_layer_file = 65535
                continue

        """
        all_layer_corner_mag_type = ctypes.c_float*(4*sum_layer_length)
        all_layer_corner_mag = all_layer_corner_mag_type()
        for i in range(4*sum_layer_length):
            all_layer_corner_mag[i] = all_layer_corner_mag_raw[i]
        
        all_layer_whether_densed_type = ctypes.c_bool*(sum_layer_length)
        all_layer_whether_densed = all_layer_whether_densed_type()
        for i in range(sum_layer_length):
            all_layer_whether_densed[i] = all_layer_whether_densed_raw[i]

        #all_layer_sequence_number_in_next_layer_file_raw = list(map(lambda x:int(x+0.5),all_layer_sequence_number_in_next_layer_file_raw))
        
        all_layer_sequence_number_in_next_layer_file_type = ctypes.c_int*(sum_layer_length)
        all_layer_sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file_type()
        for i in range(sum_layer_length):
            try:
                all_layer_sequence_number_in_next_layer_file[i] = int(all_layer_sequence_number_in_next_layer_file_raw[i]+0.5)
                
            except: 
                all_layer_sequence_number_in_next_layer_file[i] = 65535
                continue
        
        magmap,magmap_sparse = RB.ReadBin(mapdir+'%d.dat'%n,filelength=sum(filelength))
    
        datatype = ctypes.c_float*filelength
        data = datatype()
        for n in range(filelength):
            data[n] = data_raw[n]
        datatype1 = ctypes.c_float*filelength1
        data1 = datatype1()
        for n in range(filelength1):
            data1[n] = data_raw[filelength+n]
        """
        
        
        nalpha = 16
        alphalist = np.linspace(0,360.0-360.0/nalpha,nalpha)
        
        
        allparm = []
        allchi2 = []
        all_iter_number = []
        all_funcalls = []
        

        for alpha in alphalist:
            
            start = time.time()

            tempparms = [t0,u0,te,alpha]
            #parmbest,chi2min,iter_number,funcalls,warnflag,allevcs = fmin(chi2,tempparms,args=(data,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,False),full_output=True,retall=True,disp=0,maxiter=500,maxfun=1000)
            parmbest,chi2min,iter_number,funcalls,warnflag,allevcs = fmin(chi2,tempparms,args=(data, one_map_all_layer_struct, layer_length_accumulated_in_front_of_each_layer,box_size,False,s,q),full_output=True,retall=True,disp=0,maxiter=500,maxfun=1000)
            if use_mcmc==True:
                nburn_in = 500
                nsample = 1000
                ndim = len(parmbest)
                nwalkers = 2*ndim
                pos = [parmbest+1e-4*np.random.randn(ndim) for i in range(nwalkers)]
                #sampler = emcee.EnsembleSampler(nwalkers,ndim,chi2,args=(data,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,True),threads=1)
                sampler = emcee.EnsembleSampler(nwalkers,ndim,chi2,args=(data, one_map_all_layer_struct, layer_length_accumulated_in_front_of_each_layer,box_size,True,s,q),threads=1)
                ## run EMCEE ##
                pos,lnprob,rstate = sampler.run_mcmc(pos,nburn_in)#,skip_initial_state_check=True)  ## burn-in
                sampler.reset()
                sampler.run_mcmc(pos,nsample)#,skip_initial_state_check=True)
                ## save EMCEE results ##
                chain = sampler.chain.reshape((-1,ndim),order='F')
                #fsfbs = np.array(sampler.blobs).reshape((-1,2*len(data)),order='F')
                chi2s = -2*sampler.lnprobability.reshape(-1,order='F')
                chi2min = min(chi2s)
                parmbest = chain[chi2s==chi2min][0]

            allparm.append(parmbest)
            allchi2.append(chi2min)
            all_iter_number.append(iter_number)
            all_funcalls.append(funcalls)

            end   = time.time()
            total_time = ( end - start )
            print ("Total time of map %s alpha=%s: %.2f (s)"%(n,alpha,total_time))

        allchi2 = np.array(allchi2)
        allparm = np.array(allparm)
        all_iter_number = np.array(all_iter_number)
        all_funcalls = np.array(all_funcalls)

        return_2d_array = np.vstack((   allchi2, \
                                        allparm[:,0], \
                                        allparm[:,1], \
                                        allparm[:,2], \
                                        allparm[:,3], \
                                        np.ones(nalpha)*parms[n][0], \
                                        np.ones(nalpha)*parms[n][1], \
                                        np.ones(nalpha)*parms[n][2] )).T

        # arg = np.argmin(allchi2)
        # parmopt = allparm[arg]
        # chi2opt = allchi2[arg]
        # #iter_number_opt = all_iter_number[arg]
        # #funcalls_opt = all_funcalls[arg]


        #end   = time.time()
        #total_time = ( end - start )
        #print ("Total time of map %s: %.1f (s)"%(n,total_time))


        #return chi2opt,parmopt[0],parmopt[1],parmopt[2],parmopt[3],parms[n][0],parms[n][1],parms[n][2]
        return return_2d_array
    except:
        print('map %s is wrong!'%n)
        return np.array( [np.inf, 0., 0., 0., 0., logs,logq,logrho] )

    #return chi2opt,parmopt[0],parmopt[1],parmopt[2],parmopt[3],parms[n][0],parms[n][1],parms[n][2],all_iter_number[0],all_iter_number[1],all_iter_number[2],all_iter_number[3],all_iter_number[4],all_iter_number[5],all_iter_number[6],all_iter_number[7],all_funcalls[0],all_funcalls[1],all_funcalls[2],all_funcalls[3],all_funcalls[4],all_funcalls[5],all_funcalls[6],all_funcalls[7]

if __name__ == '__main__': # variable defined here is global variable

    class whether_offset_4mag(ctypes.Structure):
        _fields_ = [("whether_densed", ctypes.c_bool),
                    ("sequence_number_in_next_layer_file", ctypes.c_int),
                    ("corner_mag_1", ctypes.c_float),
                    ("corner_mag_2", ctypes.c_float),
                    ("corner_mag_3", ctypes.c_float),
                    ("corner_mag_4", ctypes.c_float)]

    eventname = 'kb251616'
    source_alpha = 265.86424999999997
    source_delta = -24.88635
    datadir = './data/%s/'%eventname
    datanames = ['KMTA01_I.pysis.dat', 'KMTA42_I.pysis.dat',\
                 'KMTC01_I.pysis.dat', 'KMTC42_I.pysis.dat',\
                 'KMTS01_I.pysis.dat', 'KMTS42_I.pysis.dat']
    fluxnames = []
    # fluxnames = ['KMTA01_I.pysis.dat'] # if 'KMTA01_I.pysis.dat' is in flux while others are in magnitude, then do this
    errfac =    [1.033, 1.071, 1.000, 1.038, 1.110, 1.056]
    errsys =    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    inverse  =  ['False', 'False', 'False', 'False', 'False', 'False']
    jd_to_hjd = ['False', 'False', 'False', 'False', 'False', 'False']



    t0 = 10857.76
    u0 = 0.003
    te = 17.55


    
    tbegin = 10820.0
    tend = 10900.0
    test = False

    use_mcmc = True
    #use_mcmc = False

    hot_mcmc_factor = 0.1


    parms = np.load('parms_%s.npy'%map_set_name)
    #parms = np.load('parms_kb220371_all_q_close.npy')
#    args = range(len(parms))
    ##############################################
    ''' hongjing: pick out a fraction of parms '''
    args = []
    parm_select = []
    p = []

    def select(ilogs,ilogq,ilogrho):
        ls,lq,lrho = round(ilogs,5),round(ilogq,5),round(ilogrho,5)

        logs_range   = [-1.50,1.50]
        logq_range   = [-6.0,4.0]
        logrho_range = [-4.0,-1.6]
        
        logs_step = 0.05
        logq_step = 0.1
        logrho_step = 0.3

        logss = np.round(np.arange(logs_range[0],logs_range[1]+0.5*logs_step,logs_step),5)
        logqs = np.round(np.arange(logq_range[0],logq_range[1]+0.5*logq_step,logq_step),5)
        logrhos = np.round(np.arange(logrho_range[0],logrho_range[1]+0.5*logrho_step,logrho_step),5)
        if   (ls not in logss):
            return False
        elif (lq not in logqs):
            return False
        elif (lrho not in logrhos):
            return False

        ### exclude a sub region ###
        #elif lq <= (-3.5/0.8*(ls+0.2)-5.5):
        #    return False
        #elif lq >= (-2.0/0.7*(ls-0.0)-2.5):
        #    return False

            ### map already has a sub region ###
        # elif (4*ls + lq + 8) < 0:
        #     return False
        # elif (-4*ls + lq + 8) < 0:
        #     return False
            ### end ###

        ### end sub region ###

        else:
            return True
    """
    for i,pi in enumerate(parms):
        ilogs,ilogq,ilogrho = pi
        if select(ilogs,ilogq,ilogrho): # == True:
            args.append(i)
        else:
            continue
    """

    mapdir = "./map_set_%s/"%map_set_name
    
    for i,pi in enumerate(parms):
        ilogs,ilogq,ilogrho = pi
         
        fp = mapdir + "%s.npz"%i
        if not os.path.isfile(fp):
            print('map: %s do not exist, so skips'%pi)
            continue
        
        if select(ilogs,ilogq,ilogrho): # == True:
            args.append(i)
            parm_select.append([ilogs,ilogq,ilogrho])
            p.append([i,ilogs,ilogq,ilogrho])
        else:
            continue


        

    '''hongjing: pick out a fraction of parms '''
    '''               - END -                 '''
    #############################################
    print('Number of grid: %d'%(len(args)))

    initial_type = ctypes.c_float*4

    data = []
    total_number_of_data_points = 0
    for i in range(len(datanames)):
        iname = datanames[i]

        ### special dealing ###
        tempdata = np.loadtxt(datadir+iname,usecols=(0,1,2))
        ###       end       ###
        
        if tempdata[0,0]>2450000:
            tempdata[:,0]-=2450000
        if tempdata[0,0]>50000:
            tempdata[:,0]-=50000
        if jd_to_hjd[i] == True:
            hjd = []
            for t in tempdata[:,0]:
                hjd.append(pyasl.helio_jd(t+50000,source_alpha,source_delta)-50000)
            tempdata[:,0] = hjd

        arg = (tempdata[:,0]>tbegin)*(tempdata[:,0]<tend)
        tempdata = tempdata[arg]

        
        if len(tempdata)==0:
            print ("No data in %s satisfy the time domain"%iname)
            continue
        if inverse[i] == True:
            tempdata[:,1] *= -1.
        if iname not in fluxnames:
            tempdata[:,2] = errfac[i]*np.sqrt(tempdata[:,2]**2+errsys[i]**2)
            tempdata[:,1] = 10.**(0.4*(18.-tempdata[:,1]))
            tempdata[:,2] = tempdata[:,2]*tempdata[:,1]*np.log(10.)/2.5
        else:
            tempdata[:,2] = tempdata[:,2]*errfac[i]
        data.append(tempdata)
        print ("%s has %d data points" %(iname,len(tempdata)))
        total_number_of_data_points += len(tempdata)

    data = np.array(data)
    print('total %d data points'%total_number_of_data_points)



    each_dataset_hjds = []
    each_dataset_fluxs = []
    each_dataset_ferrs = []

    for idata in data :
        nhjd = len(idata)
        temptype = ctypes.c_float*nhjd
        hjds = temptype()
        iflux = temptype()
        iferr = temptype()

        for i in range(nhjd):
            hjds[i]=idata[i,0]
            iflux[i]=idata[i,1]
            iferr[i]=idata[i,2]
        
        each_dataset_hjds.append(hjds)
        each_dataset_fluxs.append(iflux)
        each_dataset_ferrs.append(iferr)



    print ("Start Grid Search")

    #pool = mp.Pool(processes=300)


  
    if test == True:
        print (grid(p[3]))
    else:
        pool = mp.Pool(processes=96)
        #results = pool.map(grid,args)
        results = pool.map(grid, p)
        #results = np.array(results,dtype='float')
        results = np.vstack(results)
        print (results)
        
        np.save("result/%s_pipeline_struct_test_boxsize3p5_layer16_mcmc_16_alpha.npy"%(eventname),results)
        
    exit(0)



####### map resolution ######
### 0.1 day needs 10 points, te<10 day
"""
### need small square <= 0.05/te/10 = 5e-5 < 1e-4, but still good, just a factor of two

### tanomaly = 1 day, needs 10 points, te = 50 day
### need small square <= tanomaly/te/10 = 0.002
"""

#######  map coverage  ######
"""
### full coverage because of long te

### te = 50 day, data range = 9600-9900,t0 = 9800
### thus -4te to +2te
"""
