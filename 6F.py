# -*- coding: utf-8 -*-
import os
import json
import warnings
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage import *
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D 


from scipy.optimize import curve_fit
T_plot = np.array([1/8, 1/12, 1/15, 1/20, 1/25, 1/30, 1/40, 1/50, 1/65, 1/75, 1/100])
chi_with_plot = np.array([0.7899, 0.8936, 0.932, 0.9665, 0.9779, 0.9805, 0.9833, 0.9779, 0.9753, 0.9691, 0.9634])
chi_without_plot = np.array([0.9414, 1.13, 1.222, 1.324, 1.394, 1.437, 1.501, 1.539, 1.573, 1.595, 1.609])

def curie_weiss(x, a, Tc):
	y = a/(x-Tc)
	return y

lower = [0, -1]
upper = [1, 0]

parameters, covariance = curve_fit(curie_weiss, T_plot, chi_without_plot, bounds=(lower, upper))
# print(covariance)
fit_a = parameters[0]
# print('a:', fit_a)
fit_Tc = parameters[1]
# print("Tc:", fit_Tc)

fit_y_plot = curie_weiss(T_plot, fit_a, fit_Tc)
# print('Error:', (np.sqrt((fit_y - chi_without)**2).sum()))


matplotlib.rcParams['axes.linewidth'] = 4
plt.rc('text', usetex=True)
plt.rc('font', family='Times New Roman')
params={
        #        'font.family':'serif',
        #        'font.serif':'Times New Roman',
                        'font.weight':'normal',
                        'font.size':'15'
                                }
                                
matplotlib.rcParams.update(params)

fig = plt.figure(figsize=(22,10), facecolor="white")


plt.subplots_adjust(wspace=0.32, hspace=0.45)

########################################################################################################

# ax1 = fig.add_subplot(231)
# #plt.subplot(221)


# x, dx2y2, dzr, dxy, ins, total = np.loadtxt('./U0/Spectral_loc_U0.out', dtype=np.float, skiprows=0 , unpack=True)

# plt.xlim(-4,10)
# plt.ylim(0,2)
# plt.vlines(0,0,2, colors = "black", linestyles = "dashed")
# pdx2y2, = plt.plot(x, dx2y2, color='red',lw=2.5)
# plt.fill_between(x, dx2y2, facecolor='salmon')
# pdzr, = plt.plot(x, dzr, color='blue',lw=2.5)
# pdxy, = plt.plot(x, dxy, color='purple',lw=2.5)
# pins, = plt.plot(x, ins, color='orange',lw=2.5)
# ptotal, = plt.plot(x, total, color='green',lw=2.5)
# ax1.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax1.set_xticks([-4,-2,0,2,4,6,8,10])
# ax1.set_xticklabels([-4,-2,0,2,4,6,8,10],fontsize=30)
# ax1.set_yticks([0,0.5,1,1.5,2])
# ax1.set_yticklabels([0,'',1,'',2],fontsize=30)
# plt.legend([pdx2y2,pdzr,pdxy,pins,ptotal],[r'Ni $d_{x^2-y^2}$',r'Nd $d_{3z^2-r^2}$',r'Nd $d_{xy}$',r'interstitial $s$','total'],loc='best',frameon=False, handlelength=0.4 ,fontsize=19)
# ax1.set_title(r'$U_{\rm Ni}=0\,$eV',fontsize =30)
# plt.ylabel(r'$A(\omega)$\,(eV$^{-1}$)', fontsize=30)
# plt.text(-5, 2.2,r'\textbf{(a)}',fontsize=35, color='black')
# plt.xlabel(r'$\omega$(eV)', fontsize=30)
# #plt.xlabel(r'$\omega$(eV)', fontsize=25)

# ########################################################################################################


# ########################################################################################################

# ax2 = fig.add_subplot(232)
# #plt.subplot(221)


# x3, dx2y23, dzr3, dxy3, ins3, total3 = np.loadtxt('./U3/Spectral_loc_U3.out', dtype=np.float, skiprows=0 , unpack=True)

# plt.xlim(-4,10)
# plt.ylim(0,2)
# plt.vlines(0,0,2, colors = "black", linestyles = "dashed")
# pdx2y23, = plt.plot(x3, dx2y23, color='red',lw=2.5)
# plt.fill_between(x3, dx2y23, facecolor='salmon')
# pdzr3, = plt.plot(x3, dzr3, color='blue',lw=2.5)
# pdxy3, = plt.plot(x3, dxy3, color='purple',lw=2.5)
# pins3, = plt.plot(x3, ins3, color='orange',lw=2.5)
# ptotal3, = plt.plot(x3, total3, color='green',lw=2.5)
# ax2.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax2.set_xticks([-4,-2,0,2,4,6,8,10])
# ax2.set_xticklabels([-4,-2,0,2,4,6,8,10],fontsize=30)
# ax2.set_yticks([0,0.5,1,1.5,2])
# ax2.set_yticklabels([0,'',1,'',2],fontsize=30)
# ax2.set_title(r'$U_{\rm Ni}=3\,$eV',fontsize =30)
# plt.text(-5, 2.2,r'\textbf{(b)}',fontsize=35, color='black')
# plt.ylabel(r'$A(\omega)$\,(eV$^{-1}$)', fontsize=30)
# plt.xlabel(r'$\omega$(eV)', fontsize=30)
# #plt.ylabel(r'$A(\omega)$(eV$^{-1}$)', fontsize=25)
# #plt.xlabel(r'$\omega$(eV)', fontsize=25)

# ########################################################################################################


# ########################################################################################################

# ax3 = fig.add_subplot(233)
# x7, dx2y27, dzr7, dxy7, ins7, total7 = np.loadtxt('./U7/Spectral_loc_U7.out', dtype=np.float, skiprows=0 , unpack=True)

# plt.xlim(-4,10)
# plt.ylim(0,2)
# plt.vlines(0,0,2, colors = "black", linestyles = "dashed")
# pdx2y27, = plt.plot(x7, dx2y27, color='red',lw=2.5)
# plt.fill_between(x7, dx2y27, facecolor='salmon')
# pdzr7, = plt.plot(x7, dzr7, color='blue',lw=2.5)
# pdxy7, = plt.plot(x7, dxy7, color='purple',lw=2.5)
# pins7, = plt.plot(x7, ins7, color='orange',lw=2.5)
# ptotal7, = plt.plot(x7, total7, color='green',lw=2.5)
# ax3.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax3.set_xticks([-4,-2,0,2,4,6,8,10])
# ax3.set_xticklabels([-4,-2,0,2,4,6,8,10],fontsize=30)
# ax3.set_yticks([0,0.5,1,1.5,2])
# ax3.set_yticklabels([0,'',1,'',2],fontsize=30)
# ax3.set_title(r'$U_{\rm Ni}=7\,$eV',fontsize =30)
# plt.ylabel(r'$A(\omega)$\,(eV$^{-1}$)', fontsize=30)
# plt.xlabel(r'$\omega$(eV)', fontsize=30)
# plt.text(-5, 2.2,r'\textbf{(c)}',fontsize=35, color='black')


# ########################################################################################################



# '''
# ax4 = fig.add_subplot(235)
# #plt.subplot(221)


# x9, dx2y29, dzr9, dxy9, ins9, total9 = np.loadtxt('./U9/Spectral_loc_U9.out', dtype=np.float, skiprows=0 , unpack=True)

# plt.xlim(-4,10)
# plt.ylim(0,2)
# plt.vlines(0,0,2, colors = "black", linestyles = "dashed")
# pdx2y29, = plt.plot(x9, dx2y29, color='red',lw=2.5)
# plt.fill_between(x9, dx2y29, facecolor='salmon')
# pdzr9, = plt.plot(x9, dzr9, color='blue',lw=2.5)
# pdxy9, = plt.plot(x9, dxy9, color='purple',lw=2.5)
# pins9, = plt.plot(x9, ins9, color='orange',lw=2.5)
# ptotal9, = plt.plot(x9, total9, color='green',lw=2.5)
# ax4.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax4.set_xticks([-4,-2,0,2,4,6,8,10])
# ax4.set_xticklabels([-4,-2,0,2,4,6,8,10],fontsize=30)
# ax4.set_yticks([0,0.5,1,1.5,2])
# ax4.set_yticklabels([0,'',1,'',2],fontsize=30)
# #plt.lege9([pdx2y2,pdzr,pdxy,pins,ptotal],[r'Ni $d_{x^2-y^2}$',r'9 $d_{3z^2-r^2}$',r'9 $d_{xy}$','interstitial s','total'],loc='best',frameon=False)
# ax4.set_title(r'$U_{\rm Ni}=9\,$eV',fontsize =30)
# #plt.ylabel(r'$A(\omega)$(eV$^{-1}$)', fontsize=30)
# plt.xlabel(r'$\omega$(eV)', fontsize=30)
# plt.text(-6, 2.2,r'\textbf{(a4)}',fontsize=35, color='black')
# '''

# ########################################################################################################

# ax5 = fig.add_subplot(234)
# xnd, dx2y2nd, dzrnd, dxynd, insnd= np.loadtxt('nd', dtype=np.float, skiprows=0 , unpack=True)

# plt.xlim(0,9)
# plt.ylim(0.8,1.0)
# pdx2y2nd, = ax5.plot(xnd, dx2y2nd, color='red',lw=2.5,marker='o',markerfacecolor='none',markersize=10, markeredgewidth=2)
# plt.axvspan(0, 2, facecolor='powderblue', alpha=0.5)
# plt.axvspan(2, 9, facecolor='coral', alpha=0.5)

# #pplusnd, = ax5.plot(xnd, dzrnd+dxynd+insnd, color='red',lw=2.5)

# #plt.fill_between(xnd, dx2y2nd, facecolor='salmon')
# #ptotalnd, = plt.plot(xnd, totalnd, color='green',lw=2.5)
# ax5.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax5.set_xticks([0,1,2,3,4,5,6,7,8,9])
# ax5.set_xticklabels([0,'','',3,'','',6,'','',9],fontsize=30)
# ax5.set_yticks([0.8,0.9,1.0])
# ax5.set_yticklabels([0.0,'',0.2],fontsize=30)
# #ax5.annotate('',xy=(1, 0.15), xytext=(2.5, 0.15), transform=ax5.transAxes, arrowprops=dict(facecolor='navy',arrowstyle='->'))


# plt.text(-1, 0.9, r'$N_{\alpha}$',rotation=90,fontsize=35,)
# plt.text(-1, 1.02,r'\textbf{(d)}',fontsize=35, color='black')
# plt.xlabel(r'$U_{\rm Ni}\,$(eV)', fontsize=30)

# #ax5.set_ylim()
# #arrow(1,1,4,0,head_width=0.2)
# #arrow2 = mpatches.Arrow(x_tail, y_tail, dx, dy)
# #ax6.add_patch(arrow2)


# ax51 = ax5.twinx()  # this is the important function
# p3nd, = ax51.plot(xnd, dzrnd,color='blue',lw=2.5, marker='o',markersize=10)
# p4nd, = ax51.plot(xnd, dxynd, color='purple',lw=2.5, marker='o',markersize=10)
# p5nd, = ax51.plot(xnd, insnd, color='orange',lw=2.5, marker='o',markersize=10)
# pplusnd, = ax51.plot(xnd, dzrnd+dxynd+insnd, color='navy',lw=2.5, marker='o',markersize=10)
# ax51.set_ylim([0, 0.2])
# ax51.set_yticks([0,0.1,0.2])
# ax51.set_yticklabels([0.8,'',1.0],fontsize=30)
# #ax51.arrow(2.5, 0.15, -1.5, 0,head_width=0.1,color='navy')
# ax51.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax51.annotate('',xy=(0.5, 0.12), xytext=(2.5, 0.12), arrowprops=dict(color='navy',width=5))
# ax51.annotate('',xy=(8.5, 0.17), xytext=(6.5, 0.17), arrowprops=dict(edgecolor='red',facecolor='none',width=6))
# plt.legend([pdx2y2nd,p3nd,p4nd,p5nd,pplusnd],[r'Ni $d_{x^2-y^2}$',r'Nd $d_{3z^2-r^2}$',r'Nd $d_{xy}$',r'interstitial $s$',r'$d_{3z^2-r^2}$'+'+'+r'$d_{xy}$'+'+'+r'$s$'],loc='center right',frameon=False, handlelength=0.8,ncol=1, fontsize=19)
# #transform=ax5.transAxes,


# ########################################################################################################

# '''
# ax6 = fig.add_subplot(236)
# plt.xlim(0,9)
# plt.ylim(0,0.1)
# p3nd, = plt.plot(xnd, dzrnd,color='blue',lw=2.5, marker='o',markersize=10)
# p4nd, = plt.plot(xnd, dxynd, color='purple',lw=2.5, marker='o',markersize=10)
# p5nd, = plt.plot(xnd, insnd, color='orange',lw=2.5, marker='o',markersize=10)
# ax6.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
# ax6.set_xticks([0,1,2,3,4,5,6,7,8,9])
# ax6.set_xticklabels([0,'','',3,'','',6,'','',9],fontsize=30)
# ax6.set_yticks([0,0.05,0.1])
# ax6.set_yticklabels([0,'',0.1],fontsize=30)
# plt.xlabel(r'$U_{\rm Ni}\,$(eV)', fontsize=30)
# plt.text(-1, 0.05, r'$N_{\alpha}$',rotation=90,fontsize=35,)
# #plt.ylabel(r'N$\alpha$', fontsize=30)

# plt.text(-2, 0.11,r'\textbf{(b2)}',fontsize=35, color='black')
# plt.legend([p3nd,p4nd,p5nd],[r'Nd $d_{3z^2-r^2}$',r'Nd $d_{xy}$',r'interstitial $s$'],loc='best',frameon=False, handlelength=0.8, fontsize=19)
# '''

########################################################################################################

########################################################################################################

ax7 = fig.add_subplot(234)
u2c1,u2c2,u2c3,= np.loadtxt('chiU2', dtype=np.float, skiprows=0 , unpack=True)
u7c1,u7c2,u7c3,= np.loadtxt('chiU7', dtype=np.float, skiprows=0 , unpack=True)

curvex7 = np.arange(0.0, 1500, 1)
curvey7 = 2771.29/(curvex7+0.0613135)
plt.xlim(0,1500)
plt.ylim(0,50)
pu2, = ax7.plot(u2c2, u2c3, color='saddlebrown',lw=2.5,marker='o',markersize=10)
pu7, = ax7.plot(u7c2, u7c3, color='none',marker='o',markersize=10,markerfacecolor='dodgerblue')
curvexy7 = ax7.plot(curvex7, curvey7, color='black',ls='dashed',lw=2.5)
# = ax7.plot(u7c3,u7c1,color='black',ls='dashed',lw=1.5)
psuedo7 = Line2D([0],[0],color='black',ls='dashed',lw=2.5)
ax7.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
ax7.set_xticks([0,250,500,750,1000,1250,1500])
ax7.set_xticklabels([0,'',500,'',1000,'',1500],fontsize=30)
ax7.set_yticks([0,10,20,30,40,50])
ax7.set_yticklabels([0,10,20,30,40,50],fontsize=30)
plt.legend([pu2,pu7,psuedo7],[r'$U_{\rm Ni}=2\,$eV',r'$U_{\rm Ni}=7\,$eV','Curie-Weiss fit'],loc='best',frameon=False, handlelength=0.8,ncol=1, fontsize=19)
#,'Curie-Weiss fit'
plt.text(-200, 56,r'\textbf{(a)}',fontsize=35, color='black')
plt.ylabel(r'$\chi^{\omega=0}_{loc}(T)\,(\mu^{2}_{B}\,\rm eV^{-1})$', fontsize=24)
plt.xlabel(r'$T$(K)', fontsize=30)

# ########################################################################################################

# ########################################################################################################


ax8 = fig.add_subplot(235)
un2c1,un2c2,un2c3,= np.loadtxt('chiU2-nohyb', dtype=np.float, skiprows=0 , unpack=True)

curvex8 = np.arange(0.0, 1500, 1)
curvey8 = 2166.34/(curvex8+1200.98)

plt.xlim(0,1500)
plt.ylim(0.6,1.8)
pu28, = ax8.plot(u2c2, u2c3, color='saddlebrown',lw=2.5,marker='o',markersize=10)
pun2, = ax8.plot(un2c2, un2c3, color='none',marker='s',markersize=10,markerfacecolor='none',markeredgecolor='saddlebrown',markeredgewidth=3)
#pun2, = ax8.plot(un2c2, un2c3, color='saddlebrown',lw=2.5,marker='s',markerfacecolor='none',markersize=10,markeredgewidth=3)
curvexy8 = ax8.plot(curvex8, curvey8, color='black',ls='dashed',lw=2.5)
psuedo8 = Line2D([0],[0],color='black',ls='dashed',lw=2.5)

ax8.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
ax8.set_xticks([0,250,500,750,1000,1250,1500])
ax8.set_xticklabels([0,'',500,'',1000,'',1500],fontsize=30)
ax8.set_yticks([0.6,1.2,1.8])
ax8.set_yticklabels([0.6,1.2,1.8],fontsize=30)
ax8.set_title(r'$U_{\rm Ni}=2\,$eV',fontsize =30)
plt.legend([pu28,pun2,psuedo8],['with hybridization','without hybridization','Curie-Weiss fit'],loc='best',frameon=False, handlelength=0.8,ncol=1, fontsize=19)

plt.text(-200, 1.95,r'\textbf{(b)}',fontsize=35, color='black')
plt.ylabel(r'$\chi^{\omega=0}_{loc}(T)\,(\mu^{2}_{B}\,\rm eV^{-1})$', fontsize=24)
plt.xlabel(r'$T$(K)', fontsize=30)

########################################################################################################

########################################################################################################

ax9 = fig.add_subplot(236)
# un2c1,un2c2,un2c3,= np.loadtxt('chiU2-nohyb', dtype=np.float, skiprows=0 , unpack=True)

# curvex8 = np.arange(0.0, 1500, 1)
# curvey8 = 2166.34/(curvex8+1200.98)

plt.xlim(0,1500)
plt.ylim(0.6,1.8)
T_plot = T_plot * 1.16 * 1e4
# pu28, = ax8.plot(u2c2, u2c3, color='saddlebrown',lw=2.5,marker='o',markersize=10)
pu29, = ax9.plot(T_plot, chi_with_plot, color='saddlebrown',lw=2.5,marker='o',markersize=10)
# pun2, = ax8.plot(un2c2, un2c3, color='none',marker='s',markersize=10,markerfacecolor='none',markeredgecolor='saddlebrown',markeredgewidth=3)
pun29, = ax9.plot(T_plot, chi_without_plot, color='none',marker='s',markersize=10,markerfacecolor='none',markeredgecolor='saddlebrown',markeredgewidth=3)
# # pun2, = ax8.plot(un2c2, un2c3, color='saddlebrown',lw=2.5,marker='s',markerfacecolor='none',markersize=10,markeredgewidth=3)
# curvexy8 = ax8.plot(curvex8, curvey8, color='black',ls='dashed',lw=2.5)
curvexy9 = ax9.plot(T_plot, fit_y_plot, color='black',ls='dashed',lw=2.5)
psuedo9 = Line2D([0],[0],color='black',ls='dashed',lw=2.5)

ax9.tick_params('both', left = True, right = True, top = True, direction='in', width=3, length=7.5,)
ax9.set_xticks([0,250,500,750,1000,1250,1500])
# ax9.set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
ax9.set_xticklabels([0,'',500,'',1000,'',1500],fontsize=30)
# ax9.set_xticklabels([0,'',0.4,'',0.8,'',1.2],fontsize=30)
ax9.set_yticks([0.6,1.2,1.8])
ax9.set_yticklabels([0.6,1.2,1.8],fontsize=30)
ax9.set_title(r'$U_{\rm Ni}=2\,$eV,$U_{\rm Nd}=5\,$eV,$J_{\rm Nd}=0.7\,$eV',fontsize =20)
plt.legend([pu29,pun29,psuedo9],['with hybridization','without hybridization','Curie-Weiss fit'],loc='best',frameon=False, handlelength=0.8,ncol=1, fontsize=19)

plt.text(-200, 1.95,r'\textbf{(c)}',fontsize=35, color='black')
plt.ylabel(r'$\chi^{\omega=0}_{loc}(T)\,(\mu^{2}_{B}\,\rm eV^{-1})$', fontsize=24)
plt.xlabel(r'$T$(K)', fontsize=30)


########################################################################################################


fig.savefig('fig_3.pdf', bbox_inches='tight', dpi=300)
plt.close('all')
