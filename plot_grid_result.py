import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter,MaxNLocator
import time
from scipy.optimize import fsolve

plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'


eventname = 'kb240697'
dirpath = '/work/zhangjiyuan/AdaMGrid/result'
rname = ['pipeline_struct_test_boxsize3p5_layer16_mcmc_16_alpha']#name of result files

max_shown = 1000
hot_mcmc_factor = 0.1


# result_original = np.load('%s/%s_%s.npy'%(dirpath,eventname,rname[0]))
# result_append   = np.load('%s/%s_%s.npy'%(dirpath,eventname,rname[1]))
# #result = np.load('%s/%s_%s.npy'%(dirpath,eventname,rname[0]))
# result = np.vstack( (result_original,result_append) )
#print(result)
result = np.load('%s/%s_%s.npy'%(dirpath,eventname,rname[0]))


result = result[(np.argsort(result[:,0]))]
chi2list = result[:,0]
chi2list[np.isnan(chi2list)] = np.inf
result[:,0] = chi2list / hot_mcmc_factor



# choose logrho
chi2,t0,u0,te,alpha,logs,logq,logrho = result.T
#logrho_filter = (logrho == -)
#result = result[logrho_filter]

chi2,t0,u0,te,alpha,logs,logq,logrho = result.T
logs = np.round(logs,5)
logq = np.round(logq,5)
logrho = np.round(logrho,3)
chi2min = np.min(chi2)

all_s = np.unique(logs)
all_q = np.unique(logq)
all_rho = np.unique(logrho)
dlogs = all_s[1]-all_s[0]
dlogq = all_q[1]-all_q[0]

# sort by chi2
ind = np.argsort(chi2)[::-1]
logs,logq,alpha,logrho,chi2 = logs[ind],logq[ind],alpha[ind],logrho[ind],chi2[ind]
alpha = np.mod(alpha,360)

cmap = 'jet_r'
# cmap = hj_cmap

"""
for i_logrho in all_rho :
    result_logrho = []
    result_logs = []
    result_logq = []
    result_alpha = []
    result_chi2 = []
    for i in range(len(logrho)) :
        if np.abs(logrho[i] - i_logrho)<=0.0001 :
            result_logrho.append(logrho[i])
            result_logs.append(logs[i])
            result_logq.append(logq[i])
            result_alpha.append(alpha[i])
            result_chi2.append(chi2[i])
            
    result_logrho = np.array(result_logrho)
    result_logs = np.array(result_logs)
    result_logq = np.array(result_logq)
    result_alpha = np.array(result_alpha)
    result_chi2 = np.array(result_chi2)
"""

fig = plt.figure(figsize=(13,13))
gs = gridspec.GridSpec(3,3,hspace=0.1,wspace=0.15)
gscb = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0,1],width_ratios=[1,1,0.2])

ax_sq = plt.subplot(gs[2,0])
ax_aq = plt.subplot(gs[2,1],sharey=ax_sq)
ax_rq = plt.subplot(gs[2,2],sharey=ax_aq)



ax_sr = plt.subplot(gs[1,0],sharex=ax_sq)
ax_ar = plt.subplot(gs[1,1],sharex=ax_aq,sharey=ax_sr)

ax_sa = plt.subplot(gs[0,0],sharex=ax_sr)

cb = plt.subplot(gscb[0,2])


sqmap = ax_sq.scatter(logs,logq,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)#140
ax_aq.scatter(alpha,logq,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)
ax_rq.scatter(logrho,logq,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)

ax_sa.scatter(logs,alpha,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)

ax_sr.scatter(logs,logrho,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)
ax_ar.scatter(alpha,logrho,c=chi2-chi2min,s=100,marker='o',cmap = cmap,vmax=max_shown)

cb = plt.colorbar(sqmap,cax=cb)

"""
### resonance boundary
# wide-resonancce
logs = np.arange(1e-3,1.5,0.002)
s = 10**logs
q = (s**(2/3)-1)**3
logq = np.log10(q)
ax_sq.plot(logs,logq,c='#BBBBBB',linewidth=3,linestyle='--')
# close-resonance
logq,logs = np.load('logsc_vs_logq.npy').T
ax_sq.plot(logs,logq,c='#BBBBBB',linewidth=3,linestyle='--')
"""

ax_sq.set_xlabel(r'$\log{\ s}$')
ax_sq.set_ylabel(r'$\log{\ q}$')
ax_sa.set_ylabel(r'$\alpha$')
ax_aq.set_xlabel(r'$\alpha$')
ax_rq.set_xlabel(r'$\log{\rho}$')
ax_sr.set_ylabel(r'$\log{\rho}$')
cb.set_label(r'$\Delta\chi^2=\chi^2-\chi^2_{\rm min}$',fontsize =30)


ax_sq.set_xlim(all_s[0]-0.5*dlogs,all_s[-1]+0.5*dlogs)
ax_sq.set_ylim(all_q[0]-0.5*dlogq,all_q[-1]+0.5*dlogq)
for ax in [ax_sq,ax_sa,ax_aq,ax_sr,ax_ar,ax_rq,cb]:
    ax.minorticks_on()

    
    
theoretical_logq_array = np.linspace(-6,0,500)
theoretical_logs_close_array = np.zeros(500)
theoretical_logs_wide_array = np.zeros(500)

def close_resonant_boundary(s,q):
    #return s**8 - ( (1+q)**2 - s**4 )**3 /27./q
    return s**8 - (1+q)**2 /27./q * (1-s**4)**3

for arg in range(len(theoretical_logq_array)): 
    arg_q = 10**theoretical_logq_array[arg]
    
    arg_s_close = fsolve(close_resonant_boundary, 0.7, args=(arg_q))
    #print(i_s)
    theoretical_logs_close_array[arg] = np.log10(arg_s_close)
    
    arg_s_wide = np.sqrt( (1+arg_q**(1/3))**3 / (1+arg_q) )
    
    theoretical_logs_wide_array[arg] = np.log10(arg_s_wide)

    
ax_sq.plot(theoretical_logs_close_array, theoretical_logq_array, linewidth = 1.2 , color = 'r' , linestyle='--',label='Cassan 2008')
ax_sq.plot(theoretical_logs_wide_array, theoretical_logq_array, linewidth = 1.2 , color = 'r' , linestyle='--',label='Cassan 2008')



# plt.savefig('output/kb190505/kb190505_grid.pdf',dpi=300,bbox_inches='tight')
plt.show()

print(len(result))




######### change here ############
show_number = 400
######### change here ############

current_show_number = 0

for i in range(100):
#     #if result[i][5] > -0.1 and result[i][5] < 0.1 : # logs
#     #if result[i][4] > 0 and result[i][4] < 40 :     # alpha
#     if result[i][6] > 0 :                            # logq
#         print('%.1f  %.3f  %.4f  %.2f  %.4f  %.3f  %.5f  %.5f'%(tuple(result[i][:8])))
#         current_show_number += 1
#         if current_show_number > show_number :
#             break
    print('%.1f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f'%(tuple(result[i][:8])))
