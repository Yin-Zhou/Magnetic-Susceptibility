from scipy import *
import os,sys,subprocess,shutil,glob
from time import gmtime, localtime, strftime, time
import Struct, Fileio, copy 
import scipy.interpolate
import generate_cix
import header

############### This module performs DMFT loop. ######################

def DMFT_sigi(ommesh,T,U,J,dm,TB,Nd_atom,E_dc,norb,cor_at,cor_orb,it):
   """This function calculates self energy sigi for correlated orbitals."""

   ### define cor_at_flat ###
   cor_at_flat = [item for sublist in cor_at for item in sublist]

   ### define cor_orb_flat ###
   cor_orb_flat = []
   for i in range(len(cor_orb)):
       cor_orb_flat.append([item for sublist in cor_orb[i] for item in sublist])

   ### calculate the number of correlated orbitals (HF+DMFT) ###
   ncor_orb = 0 
   for at in cor_at_flat: ncor_orb += len(TB.TB_orbs[at])

   ### calculate the number of correlated orbitals (DMFT only) ###
   n_DMFT_orb = 0
   for i,orb in enumerate(cor_orb): n_DMFT_orb += len(cor_orb[i])

   ### label each orbital: LDA, HF and DMFT ### 
   cor_idx = zeros(ncor_orb, dtype=int)

   ### produce the index for each orbital ###
   # 0: s,p. 1: d(HF). 2,3,4,...: d(DMFT for non-eqival. orbitals=> reads from Sig.out file)
   idx1 = 0
   for i,ats in enumerate(cor_at):
      for at in ats:
         for orb in TB.TB_orbs[at]:
            idx = TB.idx[at][orb]
            if (header.s_orb.count(orb)==1 or header.p_orb.count(orb)==1):
               if (it == 0):
                  OUTPUT = open('OUTPUT', 'a')
                  OUTPUT.write('%s atom %s orbital is treated by LDA.\n'%(at, orb))
                  OUTPUT.close()
               cor_idx[idx] = 0
            elif (cor_orb_flat[i].count(orb)==0 and header.d_orb.count(orb)==1):
               if (it == 0):
                  OUTPUT = open('OUTPUT', 'a')
                  OUTPUT.write('%s atom %s orbital is treated by HF.\n'%(at, orb))
                  OUTPUT.close()
               cor_idx[idx] = 1
            elif (cor_orb_flat[i].count(orb)==1 and header.d_orb.count(orb)==1):
               if (it == 0):
                  OUTPUT = open('OUTPUT', 'a')
                  OUTPUT.write('%s atom %s orbital is treated by DMFT.\n'%(at, orb))
                  OUTPUT.close()
               for j,cor in enumerate(cor_orb[i]):
                  if cor.count(orb):
                     cor_idx[idx] = 2+idx1+j; break
            else:
               OUTPUT = open('OUTPUT', 'a')
               OUTPUT.write('ERROR: %s atom %s orbital is not supported.\n'%(at, orb))
               OUTPUT.close()
               exit() 
      idx1 += len(cor_orb[i])
   
   ### sig is the DMFT self-energy calculated directly from ctqmc ### 
   sig = zeros((n_DMFT_orb,len(ommesh)), dtype=complex)

   TrSigmaG = []; nf_qmc = [];
   loc_idx = 0
   
   for i,ats in enumerate(cor_at):
      Sig_name = 'Sig_'+str(cor_at[i][0])+'.out'
      
      ### if Sig file exists, then start from that file ###
      if (os.path.exists(Sig_name)): 
         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('%s atom: self-energy file starts from the existing file %s.\n'%(cor_at[i][0], Sig_name))
         OUTPUT.close()
         (ommesh,sig_atom,TrS,nf_q,mom) = Fileio.Read_complex_Data(Sig_name)

         ### sanity check ###
         if (len(sig_atom) != len(cor_orb[i])): 
            OUTPUT = open('OUTPUT', 'a')
            OUTPUT.write('ERROR: The number of correated orbital is not same as Sig file column.\n')
            OUTPUT.close()
            exit()
         if (len(mom) != len(cor_orb[i])): 
            OUTPUT = open('OUTPUT', 'a')
            OUTPUT.write('ERROR: The number of correated orbital is not same as mom list in Sig file.\n')
            OUTPUT.close()
            exit()

         ### calculate TrSigmaG and nf_qmc ###
         TrSigmaG.append(TrS-0.5*sum([mom[jj]*sig_atom[jj][-1].real for jj in range(len(mom))]));
         nf_qmc.append(list(mom))

         ### subtract Sig0 from DMFT self-energy ###
         for j in range(len(cor_orb[i])):
            Sig0 = sig_atom[j][-1].real
            sig[loc_idx]  = copy.deepcopy(sig_atom[j])
            sig[loc_idx] -= Sig0 
            loc_idx += 1

      ### if Sig file does not exist, start from Hartree-Fock value ###
      else:            
         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('%s atom: self-energy file starts from Hartree-Fock value.\n'%(cor_at[i][0]))
         OUTPUT.close()

         TrSigmaG.append(0.0)
         nf_qmc.append([])
         for orbs in cor_orb[i]: nf_qmc[i].append(dm[TB.idx[ats[0]][orbs[0]]])

   ### initialize sigi and sig_st ###
   sigi       = zeros((len(cor_idx),len(ommesh)), dtype=complex)
   sigi_avg   = zeros((len(cor_idx),len(ommesh)), dtype=complex)
   sigi_var   = zeros(len(cor_idx), dtype=float)
   sig_st     = zeros(len(cor_at),  dtype=float)

   ### calculate sigi and sig_st ###
   for i,ats in enumerate(cor_at):

      ### initialize the equivalent orbital index ###
      idx_equi     = [[] for ats_i in range(len(ats)*len(TB.TB_orbs[ats[0]]))] #equivalent atoms must have the same number of orbitals  
      idx_equi_orb = []

      for at in ats:
         ### calculate sig_st for each correlated atom ###     
         sig_st_total = (U[i]-2.5*J[i])*sum([dm[TB.idx[at][orb]] for orb in TB.TB_orbs[at] if cor_idx[TB.idx[at][orb]]==1])
         sig_st[i] += sig_st_total/len(ats)

         ### calculate sigi ###
         for orb in TB.TB_orbs[at]:
            ### define idx and idx_equi_orb ###
            idx = TB.idx[at][orb]
            idx_equi_orb.append([idx, orb])

            if (cor_idx[idx] == 1):
               ### HF ###
               sigi[idx,:] = 0.5*U[i]*dm[idx] + (U[i]-2.5*J[i])*(Nd_atom[i]-dm[idx])-E_dc[i] # Sig0
            elif(cor_idx[idx] > 1):
               ### DMFT ###
               sigi[idx,:] = 0.5*U[i]*dm[idx] + (U[i]-2.5*J[i])*(Nd_atom[i]-dm[idx])-E_dc[i] # Sig0

               sigi[idx,:] += copy.deepcopy(sig[cor_idx[idx]-2,:])


      ### finally subtract E_dc from sig_st ###
      sig_st[i] -= E_dc[i]

      ### calculate idx_equi before averaging ###  
      for ii in range(len(idx_equi_orb)):
         idx_equi[ii].append(idx_equi_orb[ii][0])
         for jj in range(len(idx_equi_orb)):
            if (idx_equi_orb[jj][1] == idx_equi_orb[ii][1]):
               if (jj != ii):
                  idx_equi[ii].append(idx_equi_orb[jj][0])

      ### average over equivalent atoms on each orbital ### 
      for ii in range(len(idx_equi)):
         idx = idx_equi[ii][0]
         for jj in idx_equi[ii]:
            sigi_avg[idx,:] += sigi[jj,:]/len(ats)

   ### calculate the variance of sigi and write to output ###
   for i in range(len(cor_idx)):
      for j in range(len(ommesh)):
         sigi_var[i] += (abs(sigi[i,j]-sigi_avg[i,j]))**2
      sigi_var[i] = sqrt(sigi_var[i])

   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Total variance of the self energy is: %f.\n'%(sum(sigi_var)))
   OUTPUT.close()

   ### finally update sigi with sigi_avg and print out Sig_loc ###
   sigi = sigi_avg

   ### write to output files ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Self energy is calculated.\n')
   OUTPUT.close()

   return (TrSigmaG,nf_qmc,sigi,sig_st)

def DMFT_chi_static(para_com_muk,q,temperature,Nc,mu,norb,nk,TB,Hopping,sigi):
   """This function computes static chi."""
 
   ### number of R points for Wannier functions ###
   nr = len(Hopping.keys())

   ### initialize the Hamiltonian ###
   ham = zeros((nr,norb,norb), dtype=complex)

   ### initialize the tran ###
   tran = zeros((nr,3), dtype=int)

   ### read in the non-interacting Hamiltonian and tran ###
   for i,ii in enumerate(sorted(Hopping.keys())):
      ham[i,:,:] = Hopping[ii]
   for i,ii in enumerate(sorted(Hopping.keys())):
      tran[i,:] = array(ii)

   ### prepare the input file for DMFT_chi.f ###
   CHI_INPUT = open('DMFT_chi_static.in', 'w')

   numq = len(q[:])
   CHI_INPUT.write('%s\n'%numq)

   for i in range(numq):
      for j in range(3):
         CHI_INPUT.write('%s '%q[i,j])
   CHI_INPUT.write('\n')

   CHI_INPUT.write('%s\n'%temperature)
   CHI_INPUT.write('%s\n'%Nc)
   CHI_INPUT.write('%s\n'%mu)
   for i in range(3):
     CHI_INPUT.write('%s '%nk[i])
   CHI_INPUT.write('\n')
   CHI_INPUT.write('%s\n'%nr)
   CHI_INPUT.write('%s\n'%norb)
   CHI_INPUT.write('%s\n'%len(sigi))  # ncor_orb

   for i in range(nr):
       for j in range(3):
           CHI_INPUT.write('%s '%tran[i,j])
   CHI_INPUT.write('\n')

   # currently only single-orbital case #
   for i in range(nr):
       for j in range(norb):
           for l in range(norb):
               CHI_INPUT.write('( %s , %s ) '%(ham[i,j,l].real, ham[i,j,l].imag))
   CHI_INPUT.write('\n')

   for i in range(len(sigi)):
       for j in range(-Nc, Nc):
           if (j>=0):
              CHI_INPUT.write('( %s , %s ) '%(sigi[i,j].real, sigi[i,j].imag))
           else:
              CHI_INPUT.write('( %s , %s ) '%(sigi[i,-j-1].real, -sigi[i,-j-1].imag))
   CHI_INPUT.close()
 
   ### calculate the chi ###
   OUTPUT = open('OUTPUT', 'a')

   OUTPUT.write('Parameters for calculating static chi:\n')
   OUTPUT.write('Number of points: %s.\n'%numq)
   OUTPUT.write('Temperature is: %s.\n'%temperature)
   OUTPUT.write('Chemical potential is: %s.\n'%mu)
   OUTPUT.write('Kgrid: %s.\n'%nk)
   OUTPUT.write('Cutoff of Matsubara frequencies: %s.\n'%Nc)
   OUTPUT.write('Calculating static chi...')
   OUTPUT.close()

   cmd = para_com_muk+"DMFT_chi_static > DMFT_chi_static.errors 2>> DMFT_chi_static.errors || { echo 'Parallel run failed for chi!'; exit 1; }"
   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()   
   print out, err

   ### move the files to the directory ###
   shutil.move('DMFT_chi_static.in', './CHI')
   shutil.move('DMFT_chi_static.out', './CHI')
   shutil.move('DMFT_chi_static.errors', './CHI')

   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Done.\n')
   OUTPUT.close()

   return 


def DMFT_Eliashberg(para_com_muk,initial,U,temperature,mu,mix,num_iter,nk,Nc,TB,Hopping,sigi):
   """This function computes the linearized Eliashberg equation."""

   ### number of R points for Wannier functions ###
   nr = len(Hopping.keys())

   # currently only single-orbital #
   norb = 1

   ### initialize the Hamiltonian ###
   ham = zeros((nr,norb,norb), dtype=complex)

   ### initialize the tran ###
   tran = zeros((nr,3), dtype=int)

   ### read in the non-interacting Hamiltonian and tran ###
   for i,ii in enumerate(sorted(Hopping.keys())):
      ham[i,:,:] = Hopping[ii]
   for i,ii in enumerate(sorted(Hopping.keys())):
      tran[i,:] = array(ii)

   ### prepare the input file for DMFT_chi.f ###
   ELI_INPUT = open('DMFT_Eliashberg.in', 'w')

   ELI_INPUT.write('%s\n'%initial)
   ELI_INPUT.write('%s\n'%U)
   ELI_INPUT.write('%s\n'%temperature)
   ELI_INPUT.write('%s\n'%mu)
   ELI_INPUT.write('%s\n'%mix)
   ELI_INPUT.write('%s\n'%num_iter)
   for i in range(3):
     ELI_INPUT.write('%s '%nk[i])
   ELI_INPUT.write('\n')
   ELI_INPUT.write('%s\n'%Nc)

   ELI_INPUT.write('%s\n'%nr)

   for i in range(nr):
       for j in range(3):
           ELI_INPUT.write('%s '%tran[i,j])
   ELI_INPUT.write('\n')

   # currently only single-orbital case #
   for i in range(nr):
       for j in range(norb):
           for l in range(norb):
               ELI_INPUT.write('( %s , %s ) '%(ham[i,j,l].real, ham[i,j,l].imag))
   ELI_INPUT.write('\n')

   for i in range(1):
       for j in range(-Nc, Nc):
           if (j>=0):
              ELI_INPUT.write('( %s , %s ) '%(sigi[i,j].real, sigi[i,j].imag))
           else:
              ELI_INPUT.write('( %s , %s ) '%(sigi[i,-j-1].real, -sigi[i,-j-1].imag))
   ELI_INPUT.close()

   ### calculate the chi ###
   OUTPUT = open('OUTPUT', 'a')

   OUTPUT.write('Parameters for calculating linearized Eliashberg equation:\n')
   OUTPUT.write('Temperature is: %s.\n'%temperature)
   OUTPUT.write('Chemical potential is: %s.\n'%mu)
   OUTPUT.write('Kgrid: %s.\n'%nk)
   OUTPUT.write('Cutoff of Matsubara frequencies: %s.\n'%Nc)
   OUTPUT.write('Calculating static linearized Eliashberg equation...')
   OUTPUT.close()


   cmd = para_com_muk+"DMFT_Eliashberg > DMFT_Eliashberg.errors 2>> DMFT_Eliashberg.errors || { echo 'Parallel run failed for Eliashberg!'; exit 1; }"
   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
   print out, err

   ### move the files to the directory ###
   shutil.move('DMFT_Eliashberg.in', './ELI')
   shutil.move('DMFT_Eliashberg.out', './ELI')
   shutil.move('DMFT_Eliashberg.errors', './ELI')

   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Done.\n')
   OUTPUT.close()

   return




def DMFT_ksum(para_com_muk,mu_conv,mu_iter,n_tot,mix_mu,norb,q,mu,T,TB,Hopping,ommesh,sigi,it):
   """This function computes mu and dm from non-interacting Hamiltonian and Sigma."""

   ### read in q grid ###
   qx = q[0]; qy = q[1]; qz = q[2]
   
   ### number of R points for Wannier functions ###
   nr = len(Hopping.keys())

   ### initialize the Hamiltonian ###
   ham = zeros((nr,norb,norb), dtype=complex)

   ### initialize the tran ###
   tran = zeros((nr,3), dtype=int)

   ### read in the non-interacting Hamiltonian and tran ###
   for i,ii in enumerate(sorted(Hopping.keys())):
      ham[i,:,:] = Hopping[ii]
   for i,ii in enumerate(sorted(Hopping.keys())):
      tran[i,:] = array(ii)

   ### prepare the input file for DMFT_mu.f ###
   MU_INPUT = open('DMFT_mu.in', 'w')
   MU_INPUT.write('%s\n'%qx)
   MU_INPUT.write('%s\n'%qy)
   MU_INPUT.write('%s\n'%qz)
   MU_INPUT.write('%s\n'%len(ommesh))
   MU_INPUT.write('%s\n'%T)
   MU_INPUT.write('%s\n'%mu)
   MU_INPUT.write('%s\n'%norb)
   MU_INPUT.write('%s\n'%nr)
   MU_INPUT.write('%s\n'%len(sigi))
   MU_INPUT.write('%s\n'%n_tot)
   MU_INPUT.write('%s\n'%mix_mu)
   MU_INPUT.write('%s\n'%mu_iter)
   MU_INPUT.write('%s\n'%mu_conv)

   for i in range(nr):
      for j in range(3):
         MU_INPUT.write('%s '%tran[i,j])
   MU_INPUT.write('\n')

   for i in range(nr):
      for j in range(norb):
         for l in range(norb):
            MU_INPUT.write('( %s , %s ) '%(ham[i,j,l].real, ham[i,j,l].imag))
   MU_INPUT.write('\n')

   for i in range(len(ommesh)):
      MU_INPUT.write('%s '%ommesh[i])
   MU_INPUT.write('\n')

   for i in range(len(sigi)):
      for j in range(len(ommesh)):
         MU_INPUT.write('( %s , %s ) '%(sigi[i,j].real, sigi[i,j].imag))

   MU_INPUT.close()

   ### write to output ### 
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Calculating chemical potential...')
   OUTPUT.close()

   cmd = para_com_muk+"DMFT_mu > DMFT_mu.errors 2>> DMFT_mu.errors || { echo 'Parallel run failed!'; exit 1; }"
   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
   print out, err

   ### find the chemical potential mu ###
   NEW_MU = open('DMFT_mu.out','r')
   mu = float(NEW_MU.readline())
   NEW_MU.close()

   ### write chemical potential to output files ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Done.\n')
   OUTPUT.close()

   ### store the DMFT_mu input and output files ###
   if (it != -1):
      if (it != 0):
         os.remove('./MU_KSUM/DMFT_mu.in')
         os.remove('./MU_KSUM/DMFT_mu.out')
         os.remove('./MU_KSUM/DMFT_mu.errors')
      shutil.move('DMFT_mu.in', './MU_KSUM')
      shutil.move('DMFT_mu.out', './MU_KSUM')
      shutil.move('DMFT_mu.errors', './MU_KSUM')
   
   ### prepare the input file for DMFT_ksum.f ###
   KSUM_INPUT = open('DMFT_ksum.in', 'w')
   KSUM_INPUT.write('%s\n'%qx)
   KSUM_INPUT.write('%s\n'%qy)
   KSUM_INPUT.write('%s\n'%qz)
   KSUM_INPUT.write('%s\n'%len(ommesh))
   KSUM_INPUT.write('%s\n'%T)
   KSUM_INPUT.write('%s\n'%mu)
   KSUM_INPUT.write('%s\n'%norb)
   KSUM_INPUT.write('%s\n'%nr)
   KSUM_INPUT.write('%s\n'%len(sigi))

   for i in range(nr):
       for j in range(3):
           KSUM_INPUT.write('%s '%tran[i,j])
   KSUM_INPUT.write('\n')

   for i in range(nr):
       for j in range(norb):
           for l in range(norb):
               KSUM_INPUT.write('( %s , %s ) '%(ham[i,j,l].real, ham[i,j,l].imag))
   KSUM_INPUT.write('\n')

   for i in range(len(ommesh)):
       KSUM_INPUT.write('%s '%ommesh[i])
   KSUM_INPUT.write('\n')

   for i in range(len(sigi)):
       for j in range(len(ommesh)):
           KSUM_INPUT.write('( %s , %s ) '%(sigi[i,j].real, sigi[i,j].imag))

   KSUM_INPUT.close()

   ### perform the ksum in the DMFT loop ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Running ksum...')
   OUTPUT.close()

   cmd = para_com_muk+"DMFT_ksum > DMFT_ksum.errors 2>> DMFT_ksum.errors || { echo 'Parallel run failed!'; exit 1; }"
   out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
   print out, err

   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Done.\n')
   OUTPUT.close()

   ### initialize density matrix and local Green function ###
   dm   = zeros(norb,dtype=float)
   Gloc = zeros((norb,len(ommesh)),dtype=complex)
   Ekin = 0.0

   ### read in the output of ksum ###
   KSUM_OUTPUT = open('DMFT_ksum.out','r')

   ### find the kinetic energy ###
   Ekin = float(KSUM_OUTPUT.readline())

   ### find the density matrix ###
   for i in range(norb):
      dm[i] = float(KSUM_OUTPUT.readline())

   ### find the local Green function ###
   for i in range(norb):
      for j in range(len(ommesh)):
         Gl = KSUM_OUTPUT.readline()
         Gloc[i,j] = float(Gl.split()[0])+1j*float(Gl.split()[1])

   KSUM_OUTPUT.close()

   ### store the ksum input and output files ###
   if (it != -1):
      if (it != 0):
         os.remove('./MU_KSUM/DMFT_ksum.in')
         os.remove('./MU_KSUM/DMFT_ksum.out')
         os.remove('./MU_KSUM/DMFT_ksum.errors')
      shutil.move('DMFT_ksum.in', './MU_KSUM')
      shutil.move('DMFT_ksum.out', './MU_KSUM')       
      shutil.move('DMFT_ksum.errors', './MU_KSUM')

   ### deallocate tran and ham ###
   del tran,ham

   ### for spin ###
   dm = array(dm)*2 
   Ekin = Ekin*2

   ### return the results ###
   return (mu,dm,Ekin,Gloc)


def DMFT_Delta(Gloc,sigi,ommesh,cor_at,cor_orb,TB,ed,mu):
   """This function calculates the hybridization function Delta."""

   ### calculate the local Green function and local self energy ###
   Gf_loc  = []
   Sig_loc = []
   for i,ats in enumerate(cor_at):
       Gf_loc.append([])
       Sig_loc.append([])
       for j,orbs in enumerate(cor_orb[i]):
           Gf_avg  = zeros(len(ommesh), dtype=complex)
           Sig_avg = zeros(len(ommesh), dtype=complex)
           for at in ats:
               for orb in orbs:
                   idx = TB.idx[at][orb]
                   Gf_avg  += Gloc[idx]
                   Sig_avg += sigi[idx]
           Gf_avg  /= len(ats)*len(orbs)
           Sig_avg /= len(ats)*len(orbs)#;Sig_avg-=sig_st[i]
           Gf_loc[i].append(list(Gf_avg))
           Sig_loc[i].append(list(Sig_avg))

   ### compute hybridization Delta ###
   for i in range(len(Gf_loc)): # index i labels different correlated atoms
      if len(Gf_loc[i])>0:
         Delta = zeros((len(Gf_loc[i]),len(ommesh)), dtype=complex)
         for j in range(len(Gf_loc[i])): # index j labels each orbital for every correlated atom 
            Delta[j,:] = 1j*ommesh+mu-ed[i][j]-array(Sig_loc[i][j])-1.0/array(Gf_loc[i][j]) # 1j is \sqrt(-1)

         delta_name = 'Delta_'+str(cor_at[i][0])+'.in'
         Fileio.Print_complex_multilines(Delta, ommesh, delta_name)

         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('Hybridization function for %s is calculated.\n'%cor_at[i][0])
         OUTPUT.close()

   return


def DMFT_energy(Gloc,sigi,ommesh,TB,cor_at,cor_type,cor_orb,T,TrSigmaG,dm,Ekin,Uprime,J,dc_type,Nd_atom,E_dc):
   """This function calculates the various types of energies."""

   ### calculate the local Green function and local self energy for preparation ###
   Gf_loc  = []
   Sig_loc = []
   for i,ats in enumerate(cor_at):
       Gf_loc.append([])
       Sig_loc.append([])
       for j,orbs in enumerate(cor_orb[i]):
           Gf_avg  = zeros(len(ommesh), dtype=complex)
           Sig_avg = zeros(len(ommesh), dtype=complex)
           for at in ats:
               for orb in orbs:
                   idx = TB.idx[at][orb]
                   Gf_avg  += Gloc[idx]
                   Sig_avg += sigi[idx]
           Gf_avg  /= len(ats)*len(orbs)
           Sig_avg /= len(ats)*len(orbs)#;Sig_avg-=sig_st[i]
           Gf_loc[i].append(list(Gf_avg))
           Sig_loc[i].append(list(Sig_avg))

   ### calculate the potential energy ###
   ### potential energy is computed using Eq. E_pot=1/2*Tr[G_loc*Sig_loc] ###
   ### potential energy based on ctqmc sampled values is also calculated  ###
   Epot = 0.0; Estat = 0.0; Edyna = 0.0 
   Epot_ctqmc = 0.0; Edyna_ctqmc = 0.0

   for i,ats in enumerate(cor_at):
      for at in ats:
         if (cor_type[i][0] == 'full'):
            for orb in header.d_orb:
               idx = TB.idx[at][orb]
               Estat += 0.5*sigi[idx,-1].real*dm[idx]
         elif (cor_type[i][0] == 't2g'):
            for orb in header.t2g_orb:
               idx = TB.idx[at][orb]
               Estat += 0.5*sigi[idx,-1].real*dm[idx]
         elif (cor_type[i][0] == 'd_x2y2'):
            for orb in header.dx2y2_orb:
               idx = TB.idx[at][orb]
               Estat += 0.5*sigi[idx,-1].real*dm[idx]
         elif (cor_type[i][0] == 'd_test'):
            for orb in header.dtest_orb:
               idx = TB.idx[at][orb]
               Estat += 0.5*sigi[idx,-1].real*dm[idx]
         for j,orbs in enumerate(cor_orb[i]):
            Sig0 =  Sig_loc[i][j][-1].real
            Sig1 = -Sig_loc[i][j][-1].imag*ommesh[-1]

            Edyna += len(orbs)*(2.0*T*(sum(array(Gf_loc[i][j])*(array(Sig_loc[i][j])-Sig0)).real)+2.0*T*Sig1*sum(1.0/ommesh**2.0)-Sig1/T/4.0)
         if (len(cor_orb[i]) > 0):
            Edyna_ctqmc += TrSigmaG[i]

   ### sum over E_stat and E_dyna to get E_pot ###
   Epot       = Estat + Edyna
   Epot_ctqmc = Estat + Edyna_ctqmc

   ### calculate the double counting energy ###
   Eprime = 0.0
   EdcNd  = 0.0
  
   for i in range(len(cor_at)):
      ### double counting energy Eprime ###
      if (dc_type == 1):
         Eprime += len(cor_at[i])*(Uprime[i]*Nd_atom[i]*(Nd_atom[i]-1)/2.0-5.0*J[i]*Nd_atom[i]/2.0*(Nd_atom[i]/2.0-1.0))
      elif (dc_type == 2):
         Eprime += len(cor_at[i])*(Uprime[i]*Nd_atom[i]*(0.9*Nd_atom[i])/2.0-5.0*J[i]*Nd_atom[i]/2.0*(2.0*Nd_atom[i]/5.0))
      ### 0.5*V_dc*N_d ###
      EdcNd += 0.5*(E_dc[i]*Nd_atom[i])*len(cor_at[i])

   ### finally calculate the total energy ###
   Etot = Epot  + Ekin + EdcNd - Eprime

   return (Epot,Epot_ctqmc,EdcNd,Eprime,Etot)


def Output_density(mu,E_dc,TB,cor_at,cor_type,dm):
   """This function generates the DENSITY file."""

   ### write the DENSITY file ###
   DENSITY = open('DENSITY','w')

   DENSITY.write('mu=\t%.16f\n'%mu)
   DENSITY.write('Edc=')
   DENSITY.writelines('\t%.16f'%item for item in E_dc)
   DENSITY.write('\n')

   ### initialize Nd, Nd_eg and Nd_t2g ###
   Nd_eg  = zeros(len(cor_at), dtype = float)
   Nd_t2g = zeros(len(cor_at), dtype = float)
   Nd     = zeros(len(cor_at), dtype = float)

   ### write the occupancy of eg ### 
   DENSITY.write('Nd_eg=')
   for i,ats in enumerate(cor_at):
       if (cor_type[i][0] == 'full'):
          for at in ats:
              Nd_eg[i] += float(dm[TB.idx[at]['d_z2']])+float(dm[TB.idx[at]['d_x2y2']])
          Nd_eg[i] /= len(ats)
          DENSITY.write('\t%.16f'%(Nd_eg[i]))
       elif (cor_type[i][0] == 't2g'):
          DENSITY.write('\tNA')
       elif (cor_type[i][0] == 'd_x2y2'):
          for at in ats:
              Nd_eg[i] += float(dm[TB.idx[at]['d_x2y2']])
          Nd_eg[i] /= len(ats)
          DENSITY.write('\t%.16f'%(Nd_eg[i]))
       elif (cor_type[i][0] == 'd_test'):
          for at in ats:
              Nd_eg[i] += float(dm[TB.idx[at]['d_x2y2']])
          Nd_eg[i] /= len(ats)
          DENSITY.write('\t%.16f'%(Nd_eg[i]))
   DENSITY.write('\n')

   ### write the occupancy of t2g ###
   DENSITY.write('Nd_t2g=')
   for i,ats in enumerate(cor_at):
       if (cor_type[i][0] == 'd_x2y2' or cor_type[i][0] == 'd_test'):
          DENSITY.write('\tNA')
       else:
          for at in ats:
             Nd_t2g[i] += float(dm[TB.idx[at]['d_xz']])+float(dm[TB.idx[at]['d_yz']])+float(dm[TB.idx[at]['d_xy']])
          Nd_t2g[i] /= len(ats)
          DENSITY.write('\t%.16f'%(Nd_t2g[i]))
   DENSITY.write('\n')

   ### write the occupancy of Nd ###
   DENSITY.write('Nd=')
   for i,ats in enumerate(cor_at):
       Nd[i] = Nd_eg[i] + Nd_t2g[i]
       DENSITY.write('\t%.16f'%(Nd[i]))
   DENSITY.write('\n')

   ### write the occupancy of each orbital ###
   DENSITY.writelines('%.16f\t'%item for item in dm)
   DENSITY.write('\n')

   DENSITY.close()

   return (Nd_eg,Nd_t2g,Nd)


def Create_PARAMS(params_ctqmc):
   """This function creates the input file (PARAMS) for ctqmc."""

   PARAMS = open('PARAMS', 'w')
   PARAMS.write('# Input file for continuous time quantum Monte Carlo\n')
   for item in params_ctqmc:
       PARAMS.write('%s\t%s\t%s\n'%(item, params_ctqmc[item][0], params_ctqmc[item][1]))
   PARAMS.close()

   return


def Compute_DMFT(Hopping,params,params2,params_ctqmc,para_com,para_com_muk):
   """This is the main function for DMFT loop."""

   ### delete previous output file if it exists ###
   if os.path.exists('OUTPUT'):
      os.remove('OUTPUT')

   ### start timing ###
   starting = time()
   start_time = strftime("%a, %d %b %Y %H:%M:%S", localtime())

   ### open the output file ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Self-consistent-loop v%s starts at %s.\n\n'%(header.version,start_time))
   OUTPUT.close()

   ### read in parameters ###
   Niter = params['Niter'][0]
   n_tot = params['n_tot'][0] 
   nom = params['nom'][0]
   atomnames = params['atomnames'][0]
   cor_at = params['cor_at'][0]
   cor_type = params['cor_type'][0]
   cor_orb = params['cor_orb'][0]
   co_at = params['co_at'][0]
   mix_mu = params['mix_mu'][0]
   mix_sig = params['mix_sig'][0]
   Nd_f = params['Nd_f'][0]
   mu_conv = params['mu_conv'][0]
   mu_iter = params['mu_iter'][0]
   mu = params['mu'][0]
   Ethreshold = params['e_threshold'][0]
   dc_type = params['dc_type'][0]
   q = params['q'][0]
   U = params['U'][0]
   J = params['J'][0]
   Uprime = params['Uprime'][0]
   Jterm = params['Jterm'][0]
   T = 1.0/params_ctqmc['beta'][0]

   norb=int(len(Hopping[Hopping.keys()[0]]))

   ### initialize some variables ###
   dir_name = []
   Etot     = 0.0
   Etot_old = 0.0
   dm       = zeros(norb, dtype=float)

   ### create Matsubara frequency ###
   ommesh = zeros(nom, dtype=float)
   for i in range(nom):
       ommesh[i] = pi*T*(2*i+1)

   ### sanity check ###
   if (len(U) != len(cor_at) or len(J) != len(cor_at) or len(Uprime) != len(cor_at)):
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write("ERROR: set U or J or Uprime for all the inquivalent correlated atoms.\n")
      OUTPUT.close()
      exit()
   if (len(cor_type) != len(cor_at) or len(Jterm) != len(cor_at)):
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('ERROR: set cor_type or Jterm for all the inquivalent correlated atoms.\n')
      OUTPUT.close()
      exit()
   if (dc_type != 1 and dc_type != 2):
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('ERROR: there is something wrong with dc_type.\n')
      OUTPUT.close()
      exit()
   for i in range(len(Jterm)):
      if (Jterm[i] != 0 and Jterm[i] != 1):
         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('ERROR: there is something wrong with Jterm.\n')
         OUTPUT.close()
         exit()
   for i in range(len(cor_type)):
      if (cor_type[i][0] != 'full' and cor_type[i][0] != 't2g' and cor_type[i][0] != 'd_x2y2' and cor_type[i][0] != 'd_test'):
         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('ERROR: there is something wrong with cor_type.\n')
         OUTPUT.close()
         exit()


   ###### write to output file ######
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('===== Basic Information =====\n')

   ### read the file ORBITAL ###
   TB = Struct.TBstructure('ORBITAL',atomnames)
   TB.Construct_idx()

   ### read TB Hamiltonian and calculate ed ###
   H0 = copy.deepcopy(array(Hopping[(0,0,0)]))
   ed = []  # onsite energy of each correlated Wannier orbital

   ### write the orbital information to output file ###
   OUTPUT.write('Total number of orbitals is: %i.\n'%(norb))

   for i,ats in enumerate(cor_at):
       ed.append([])
       for j,orbs in enumerate(cor_orb[i]):
           ed_loc = 0
           for at in ats:
               for orb in orbs:
                   OUTPUT = open('OUTPUT', 'a')
                   OUTPUT.write('Correlated atom is: %s.\n'%(at))
                   OUTPUT.write('Correlated orbital is: %s.\n'%(orb))
                   OUTPUT.close()
                   ed_loc += H0[TB.idx[at][orb],TB.idx[at][orb]].real
           ed[i].append(ed_loc/len(ats)/len(orbs))
           OUTPUT = open('OUTPUT', 'a')
           OUTPUT.write('Onsite energy Ed[%s][%s](eV) is: %f.\n'%(i, j, ed[i][j]))

   ### write other parameters to output file ###
   OUTPUT.write("Correlated atoms are: ")
   for i in range(len(cor_at)):
       OUTPUT.write("  %s"%(cor_at[i][0]))
   OUTPUT.write(".\n")
   OUTPUT.write("Hund's term: ")
   for i in range(len(Jterm)):
       if (Jterm[i] == 0):
          OUTPUT.write("\tNo")
       elif (Jterm[i] == 1):
          OUTPUT.write("\tYes")
   OUTPUT.write(".\n")
   OUTPUT.write('U(eV)   =')
   OUTPUT.writelines('\t%f'%item for item in U)
   OUTPUT.write('.\n')
   OUTPUT.write('Up(eV)  =')
   OUTPUT.writelines('\t%f'%item for item in Uprime)
   OUTPUT.write('.\n')
   OUTPUT.write('J(eV)   =') 
   OUTPUT.writelines('\t%f'%item for item in J)
   OUTPUT.write('.\n')
   OUTPUT.write('Mix_mu  =\t%f.\n'%mix_mu)
   OUTPUT.write('Mix_sig =\t%f.\n'%mix_sig)
   OUTPUT.write('Chemical potential convergence threshold = %f eV.\n'%mu_conv)
   OUTPUT.write('Maximum number of iterations for chemical potential = %i.\n'%mu_iter)
   OUTPUT.write('Kpoint sampling = %s.\n'%q)
   OUTPUT.write('Temperature = %f eV.\n'%T)
   OUTPUT.write('Energy convergence threshold = %f eV.\n'%Ethreshold)
   OUTPUT.write('Maximum number of iterations for total energy = %i.\n'%Niter)
   OUTPUT.write('Number of Matsubara frequencies = %i.\n'%nom)
   OUTPUT.write('Dynamical mean field calculations.\n')
   OUTPUT.write('Paramagnetic calculations.\n\n')
   OUTPUT.close()

   ### write to output ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('===== Initialization =====\n')

   ### read GF_SIG/Sig_loc.final file if it exists ###
   if (os.path.exists('GF_SIG/Sig_loc.final')):

      OUTPUT.write('Self energy is read from the previous calculations.\n')
      OUTPUT.close()

      sigi     = genfromtxt('GF_SIG/Sig_loc.final')
      ### strip the column of frequency 
      sigi     = sigi[:,1:]  
      ### combine the real and imaginary parts to make a complex
      sigi_combined     = zeros((len(sigi), int(len(sigi[0])/2)), dtype=complex)

      for i in range(len(sigi)):
         for j in range(len(sigi[0])/2):
            sigi_combined[i,j]     =  sigi[i,2*j]+sigi[i,2*j+1]*1j

      sigi_combined     = transpose(sigi_combined)
      sigi              = sigi_combined
      sigi_ref          = sigi

   ### if GF_SIG/Sig_loc.final file does not exist, compute non-int DENSITY file ###
   else: 

      OUTPUT.write('Start from a zero self energy.\n')
      OUTPUT.close()

      ### define cor_at_flat ###
      cor_at_flat = [item for sublist in cor_at for item in sublist]

      ### calculate ncor_orb ###
      ncor_orb = 0
      for at in cor_at_flat: ncor_orb += len(TB.TB_orbs[at])

      ### get the self energy, which is zero for non-interacting case ###
      sigi = zeros((ncor_orb,len(ommesh)), dtype=complex)


   ### from the self energy sigi, get the dm for initialization ###
   it = -1
   (mu,dm,Ekin,Gloc) = DMFT_ksum(para_com_muk,mu_conv,mu_iter,n_tot,mix_mu,norb,q,mu,T,TB,Hopping,ommesh,sigi,it)

   ### compute Nd from dm for initialization ### 
   Nd_atom = zeros(len(cor_at), dtype=float)
   for i,ats in enumerate(cor_at):
       for at in ats:
           for orb in TB.TB_orbs[at]:
               idx = TB.idx[at][orb]
               if (header.d_orb.count(orb) == 1):
                  Nd_atom[i] += dm[idx]
       Nd_atom[i] /= len(ats)

   ### compute double counting energy for initialization ###
   E_dc = zeros(len(cor_at), dtype=float)
   for i in range(len(cor_at)):         
       if (dc_type == 1):
          if (cor_type[i][0] == 'full'): 
             E_dc[i] = (Uprime[i]-2.0*J[i])*(Nd_atom[i]-0.5)-J[i]/2.0*(Nd_atom[i]-3.0)
          elif (cor_type[i][0] == 't2g'):
             E_dc[i] = (Uprime[i]-4.0/3.0*J[i])*(Nd_atom[i]-0.5)-5.0/6.0*J[i]*(Nd_atom[i]-1.0)
          elif (cor_type[i][0] == 'd_x2y2'):
             E_dc[i] = 0.0  # for d-only calculation
          elif (cor_type[i][0] == 'd_test'): 
             E_dc[i] = (Uprime[i]-2.0*J[i])*(Nd_atom[i]-0.5)-J[i]/2.0*(Nd_atom[i]-3.0)
       elif (dc_type == 2):
          E_dc[i] = (Uprime[i]-2.0*J[i])*(0.9*Nd_atom[i])-J[i]/2.0*(2.0*Nd_atom[i]/5.0)

   ### write to output files ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Chemical potential is: %f eV.\n'%mu)
   OUTPUT.write('Density matrix is: \n')
   OUTPUT.writelines('%f\t'%item for item in dm)
   OUTPUT.write('\nTotal occupancy is: %f.\n'%sum(dm))
   OUTPUT.write('Nd is: \t')
   OUTPUT.writelines('\t%f'%item for item in Nd_atom)
   OUTPUT.write('.\n')
   OUTPUT.write('Edc is: ')
   OUTPUT.writelines('\t%f'%item for item in E_dc)
   OUTPUT.write('.\n')
   OUTPUT.close()

   ############ DMFT Self-consistent Loop ##############

   ### timing up till now ###
   uptillnow = time()

   ### write elapsed time to output file ###
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write('Time elapsed is: %f s.\n\n'%(uptillnow - starting))
   OUTPUT.write('===== Entering DMFT Loops =====')
   OUTPUT.close()

   ### loop for Niter iterations ###
   for it in range(Niter):

      ###### set up ######
      ### write to output file ###
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('\nIteration %i of DMFT loop:\n'%(it+1))
      OUTPUT.close()

      ### creat imp_X, MU_KSUM and GF_SIG directories for the first iteration ###
      if (it == 0):
         OUTPUT = open('OUTPUT', 'a')
         ### imp_X ###
         for i in range(len(cor_at)):
            dir_name.append('imp_'+str(cor_at[i][0])+'/')
            if (len(cor_orb[i]) > 0):
               if os.path.exists(dir_name[i]):
                  OUTPUT.write('The %s directory exists.\n'%dir_name[i])
               else:
                  os.makedirs(dir_name[i])
                  OUTPUT.write('Creating the directory: %s.\n'%dir_name[i])  
         ### MU_KSUM ###
         if os.path.exists('MU_KSUM'):
            shutil.rmtree('MU_KSUM')
         os.makedirs('MU_KSUM')
         OUTPUT.write('Creating the directory: MU_KSUM.\n')
         ### GF_SIG ###
         if os.path.exists('GF_SIG'):
            OUTPUT.write('The GF_SIG directory exists.\n')
         else:
            os.makedirs('GF_SIG')
            OUTPUT.write('Creating the directory: GF_SIG.\n')

         OUTPUT.close()

      ### calculate the self energy sigi ###  
      (TrSigmaG,nf_qmc,sigi,sig_st) = DMFT_sigi(ommesh,T,U,J,dm,TB,Nd_atom,E_dc,norb,cor_at,cor_orb,it)

      ### update the self energy sigi ###
      if (it == 0 and os.path.exists('GF_SIG/Sig_loc.final')):
         sigi = sigi_ref
      if (it > 0): 
         sigi = sigi_old+mix_sig*(sigi-sigi_old)

      sigi_old = copy.deepcopy(sigi)

      ### perform the ksum ###
      (mu,dm,Ekin,Gloc) = DMFT_ksum(para_com_muk,mu_conv,mu_iter,n_tot,mix_mu,norb,q,mu,T,TB,Hopping,ommesh,sigi,it)

      ### calculate hybridization function ###
      DMFT_Delta(Gloc,sigi,ommesh,cor_at,cor_orb,TB,ed,mu)

      ### update Nd_atom from dm for E_dc ###
      Nd_atom = zeros(len(cor_at), dtype=float)
      for i,ats in enumerate(cor_at):
          for at in ats:
              for orb in TB.TB_orbs[at]:
                  idx = TB.idx[at][orb]
                  if (header.d_orb.count(orb) == 1):
                     Nd_atom[i] += dm[idx]
          Nd_atom[i] /= len(ats)

      ### update double counting energy from Nd_atom ###
      E_dc = zeros(len(cor_at), dtype=float)
      for i in range(len(cor_at)):
         if (dc_type == 1):
            if (cor_type[i][0] == 'full'):
                E_dc[i] = (Uprime[i]-2.0*J[i])*(Nd_atom[i]-0.5)-J[i]/2*(Nd_atom[i]-3.0)
            elif (cor_type[i][0] == 't2g'):
                E_dc[i] = (Uprime[i]-4.0/3.0*J[i])*(Nd_atom[i]-0.5)-5.0/6.0*J[i]*(Nd_atom[i]-1.0)
            elif (cor_type[i][0] == 'd_x2y2'):
                E_dc[i] = 0.0  # for d-only calculations, no double counting
            elif (cor_type[i][0] == 'd_test'):
                E_dc[i] = (Uprime[i]-2.0*J[i])*(Nd_atom[i]-0.5)-J[i]/2*(Nd_atom[i]-3.0)
         elif (dc_type == 2):
            E_dc[i]=(Uprime[i]-2.0*J[i])*(0.9*Nd_atom[i])-J[i]/2.0*(2.0*Nd_atom[i]/5.0)

      ### calculate various types of energy based on Nd and Edc ###
      (Epot,Epot_ctqmc,EdcNd,Eprime,Etot) = DMFT_energy(Gloc,sigi,ommesh,TB,cor_at,cor_type,cor_orb,T,TrSigmaG,dm,Ekin,Uprime,J,dc_type,Nd_atom,E_dc)

      ###### solve the impurity problem for each correlated atom ######
      for i in range(len(cor_at)):
          if (len(cor_orb[i]) > 0):
             ### provide mu, U and J for PARAMS ###
             params_ctqmc['mu'] = [float(mu-ed[i][0]-sig_st[i]), "# Chemical potential"] # mu for impurity solver needs to be shifted
             params_ctqmc['U'] = [float(U[i]), "# Coulomb repulsion (F0)"]
             params_ctqmc['J'] = [float(J[i]), "# Hund's coupling"]

             ### provide Delta and Cix file names for PARAMS ###
             params_ctqmc['Delta'] = ['Delta_'+str(cor_at[i][0])+'.in', "# Input hybridization function"]
             params_ctqmc['cix']   = ['impurity_'+str(cor_at[i][0])+'.cix', "# Input file with atomic state"]
              
             ### create PARAMS file ###
             Create_PARAMS(params_ctqmc) # mu needs to be updated

             ### create Cix file ###
             generate_cix.Print_cix(cor_orb[i],cor_at[i][0],J[i],ed[i],Jterm[i],it)

             ### shuffle STATUS files before execution ###
             if ((it == 0 and os.path.isfile(dir_name[i]+'status.000')) or it > 0):
                Status_name = glob.glob(dir_name[i]+'status*')
                for file in Status_name: shutil.move(file, './')

             ### running ctqmc... ###  
             OUTPUT = open('OUTPUT', 'a')
             OUTPUT.write('Running ctqmc for atom %s...'%(cor_at[i][0]))
             OUTPUT.close()
     
             cmd = para_com+"ctqmc > ctqmc.out  2> ctqmc.errors || { echo 'Parallel run failed!'; exit 1; }"
             out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
             print out, err

             OUTPUT = open('OUTPUT', 'a')
             OUTPUT.write('Done.\n')
             OUTPUT.close()

             ### store data after each iteration ###
             Params_name = 'PARAMS_'+str(cor_at[i][0])+'.in' 
             Params_name_history = dir_name[i]+Params_name

             Delta_name = 'Delta_'+str(cor_at[i][0])+'.in'
             Delta_name_history = dir_name[i]+Delta_name+str(it+1)

             impurity_name = 'impurity_'+str(cor_at[i][0])+'.cix'
             impurity_name_history = dir_name[i]+impurity_name

             Gf_name = 'Gf_'+str(cor_at[i][0])+'.out'
             Gf_name_history = dir_name[i]+Gf_name+str(it+1)

             Sig_name = 'Sig_'+str(cor_at[i][0])+'.out'
             Sig_name_history = dir_name[i]+Sig_name+str(it+1)

             Gtau_name = 'Gtau_'+str(cor_at[i][0])+'.out'
             Gtau_name_history = dir_name[i]+Gtau_name+str(it+1)

             Prob_name = 'Probability_'+str(cor_at[i][0])+'.dat'
             Prob_name_history = dir_name[i]+Prob_name+str(it+1)

             hist_name = 'Histogram_'+str(cor_at[i][0])+'.dat'
             hist_name_history = dir_name[i]+hist_name+str(it+1)

             g_hb0_name = 'g_hb0_'+str(cor_at[i][0])+'.dat'
             g_hb0_name_history = dir_name[i]+g_hb0_name+str(it+1) 
           
             g_hb1_name = 'g_hb1_'+str(cor_at[i][0])+'.dat'
             g_hb1_name_history = dir_name[i]+g_hb1_name+str(it+1)

             s_hb1_name = 's_hb1_'+str(cor_at[i][0])+'.dat'
             s_hb1_name_history = dir_name[i]+s_hb1_name+str(it+1)

             g_qmc_name = 'g_qmc_'+str(cor_at[i][0])+'.dat'
             g_qmc_name_history = dir_name[i]+g_qmc_name+str(it+1)

             ctqmc_o_name = 'ctqmc_'+str(cor_at[i][0])+'.out'
             ctqmc_o_name_history = dir_name[i]+ctqmc_o_name+str(it+1)

             ctqmc_e_name = 'ctqmc_'+str(cor_at[i][0])+'.errors'
             ctqmc_e_name_history = dir_name[i]+ctqmc_e_name+str(it+1)

             if os.path.isfile(Params_name_history): os.remove(Params_name_history)
             shutil.move('PARAMS', Params_name_history)
             shutil.move('Delta_'+str(cor_at[i][0])+'.in', Delta_name_history)
             if os.path.isfile(impurity_name_history): os.remove(impurity_name_history)
             shutil.move('impurity_'+str(cor_at[i][0])+'.cix', impurity_name_history)
             shutil.move('Gf.out', Gf_name_history)
             shutil.copyfile('Sig.out', Sig_name_history)
             os.rename('Sig.out', Sig_name) # for the next iteration
             shutil.move('Gtau.dat', Gtau_name_history)
             shutil.move('Probability.dat', Prob_name_history)
             shutil.move('Histogram.dat', hist_name_history)
             shutil.move('g_hb0.dat', g_hb0_name_history)
             shutil.move('g_hb1.dat', g_hb1_name_history)
             shutil.move('g_qmc.dat', g_qmc_name_history)
             shutil.move('s_hb1.dat', s_hb1_name_history)
             shutil.move('ctqmc.out', ctqmc_o_name_history)
             shutil.move('ctqmc.errors', ctqmc_e_name_history)

             ### shuffle the STATUS files after execution ### 
             Status_name = glob.glob('status.*') 
             for file in Status_name:  shutil.move(file, dir_name[i]) 
     
      ###### output ######
      ### print out the local Green function and spectra function ###
      Gf_loc_name = 'Gf_loc.out'
      Gf_loc_name_history = 'GF_SIG/'+Gf_loc_name+str(it+1)
      Fileio.Print_complex_multilines(Gloc,ommesh,Gf_loc_name)
      shutil.move('Gf_loc.out', Gf_loc_name_history)

      Sig_loc_name = 'Sig_loc.out'
      Sig_loc_name_history = 'GF_SIG/'+Sig_loc_name+str(it+1)

      Fileio.Print_complex_multilines(sigi,ommesh,Sig_loc_name)
      shutil.move('Sig_loc.out', Sig_loc_name_history)

      ### write the DENSITY file ###
      (Nd_eg,Nd_t2g,Nd) = Output_density(mu,E_dc,TB,cor_at,cor_type,dm)

      ### write to output files ###
      OUTPUT = open('OUTPUT', 'a')
      ### write chemical potential to output files ###
      OUTPUT.write("Chemical potential is: %f eV.\n"%(mu))

      ### write density matrix and n_tot to output files ###       
      OUTPUT.write('Density matrix is: \n')
      OUTPUT.writelines('%f\t'%item for item in dm)
      OUTPUT.write('\nTotal occupancy is: %f.\n'%sum(dm))

      ### write Nd_atom and Edc to output file ###
      OUTPUT.write('Nd is: \t')
      OUTPUT.writelines('\t%f'%item for item in Nd_atom)
      OUTPUT.write('.\n')
      OUTPUT.write('Nd(eg) is: \t')
      OUTPUT.writelines('\t%f'%item for item in Nd_eg)
      OUTPUT.write('.\n')
      OUTPUT.write('Nd(t2g) is: \t')
      OUTPUT.writelines('\t%f'%item for item in Nd_t2g)
      OUTPUT.write('.\n')
      OUTPUT.write('Edc is: ')
      OUTPUT.writelines('\t%f'%item for item in E_dc)
      OUTPUT.write('.\n')

      ### write energies to output files ###
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('Kinetic energy is: %f eV.\n'%Ekin)
      OUTPUT.write('Potential energy (potential energy from ctqmc) is: %f (%f) eV.\n'%(Epot,Epot_ctqmc))
      OUTPUT.write('EdcNd is: %f eV.\n'%EdcNd)
      OUTPUT.write('Eprime is: %f eV.\n'%Eprime)
      OUTPUT.write('Total energy is: %f eV.\n'%Etot)
      OUTPUT.close()
                    

      ### finally update the energy and check the energy difference ###
      if (it != 0):
         Ediff = Etot - Etot_old

         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write('Energy difference between two consecutive iterations is: %f meV.\n'%(Ediff*1000))
         OUTPUT.close()

         if (abs(Ediff) < abs(Ethreshold)):
            OUTPUT = open('OUTPUT', 'a')
            OUTPUT.write('\nEnergy difference is smaller than the threshold. Exit the loop.\n')
            OUTPUT.close()
            break
      Etot_old = Etot

      ### timing up till now ###
      uptillnow = time()

      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('Time elapsed is: %.f s.\n'%(uptillnow - starting))
      OUTPUT.close()

      ### for the last iteration ###
      if (it+1 == Niter):
         OUTPUT = open('OUTPUT', 'a')
         OUTPUT.write("\nMaximum number of iteration is reached.\n")
         OUTPUT.close()

   ### loop ends here ###

   ### write the final results ###                                          
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write("\n\n")
   OUTPUT.write("===== DMFT Final Results =====\n")
   OUTPUT.write("Total energy is: %f eV.\n"%(Etot))
   OUTPUT.write("Chemical potential is: %f eV.\n"%(mu))
   OUTPUT.write("Total occupancy is: %f.\n"%(sum(dm)))
   OUTPUT.write("Correlated atoms are: ")
   for i in range(len(cor_at)):
       OUTPUT.write(" %s"%(cor_at[i]))
   OUTPUT.write("\n")

   OUTPUT.write("Occupancies of the correlated atoms are: \n")
   for i,ats in enumerate(cor_at):
       if (cor_type[i][0] == 'full'):
          OUTPUT.write("Nd(%s)     is: %f.\n"%(ats[0], Nd[i]))
          OUTPUT.write("Nd_eg(%s)  is: %f.\n"%(ats[0], Nd_eg[i]))
          OUTPUT.write("Nd_t2g(%s) is: %f.\n"%(ats[0], Nd_t2g[i]))
       elif (cor_type[i][0] == 't2g'):
          OUTPUT.write("Nd(%s)     is: %f.\n"%(ats[0], Nd[i]))
          OUTPUT.write("Nd_eg(%s)  is: NA.\n"%(ats[0]))
          OUTPUT.write("Nd_t2g(%s) is: %f.\n"%(ats[0], Nd_t2g[i]))
       elif (cor_type[i][0] == 'd_x2y2'):
          OUTPUT.write("Nd(%s)     is: %f.\n"%(ats[0], Nd[i]))
          OUTPUT.write("Nd_eg(%s)  is: %f.\n"%(ats[0], Nd_eg[i]))
          OUTPUT.write("Nd_t2g(%s) is: NA.\n"%(ats[0]))          
   OUTPUT.write('Density matrix is:\n')
   OUTPUT.writelines("%f\t"%item for item in dm)
   OUTPUT.write("\n")

   OUTPUT.write("Sig0 is: \n")

   for j in range(len(sigi)):
       OUTPUT.write("%f\t"%(sigi[j, -1].real))
   OUTPUT.write("\n")
   OUTPUT.write("Sig1 is: \n")
   for j in range(len(sigi)):
       OUTPUT.write("%f\t"%(-sigi[j, -1].imag*ommesh[-1]))
   OUTPUT.write('\n')

   uptillnow = time()
   OUTPUT.write('DMFT loop ends. \nTime elapsed is: %.f s.\n'%(uptillnow - starting))
   OUTPUT.close()

   ### doing two-particle response function ###

   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write("\n\n")
   OUTPUT.write("===== Two-particle Response Functions =====\n")
   OUTPUT.close()

   ### calculate chi ###

   do_chi_static = params2['chi_static'][0]
   qhigh = params2['q_high'][0]
   numq_cut = params2['numq_cut'][0]
   nk = params2['k_points'][0]
   Nc = params2['Nc'][0]
 
   do_Eliashberg = params2['Eliashberg'][0]
   initial = params2['initial'][0]
   U = params2['U_single'][0]
   mix = params2['mix'][0]
   num_iter = params2['num_iter'][0] 
   temperature = 1.0/params_ctqmc['beta'][0]


   numqhigh = len(qhigh)
   numq = (numqhigh-1)*numq_cut + numqhigh
   q = zeros((numq,3), dtype=float)
   qindex = 0

   for i in range(numqhigh-1):
      for j in range(numq_cut+1):
         q[qindex][0] = 2*pi*(qhigh[i][0]+(qhigh[i+1][0]-qhigh[i][0])*float(j)/float(numq_cut+1))
         q[qindex][1] = 2*pi*(qhigh[i][1]+(qhigh[i+1][1]-qhigh[i][1])*float(j)/float(numq_cut+1))
         q[qindex][2] = 2*pi*(qhigh[i][2]+(qhigh[i+1][2]-qhigh[i][2])*float(j)/float(numq_cut+1))
         qindex = qindex + 1      
   q[-1][0] = 2*pi*qhigh[-1][0]
   q[-1][1] = 2*pi*qhigh[-1][1]
   q[-1][2] = 2*pi*qhigh[-1][2] 

   if (do_chi_static == 1):

      if os.path.exists('CHI'):
         shutil.rmtree('CHI')
      os.makedirs('CHI')
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('Creating the directory: CHI.\n')
      OUTPUT.close()

      DMFT_chi_static(para_com_muk,q,temperature,Nc,mu,norb,nk,TB,Hopping,sigi)

      uptillnow = time()
      OUTPUT = open('OUTPUT', 'a')   
      OUTPUT.write('Time elapsed is: %.f s.\n'%(uptillnow - starting))
      OUTPUT.close()

   if (do_Eliashberg == 1):

      if os.path.exists('ELI'):
         shutil.rmtree('ELI')
      os.makedirs('ELI')
      OUTPUT = open('OUTPUT', 'a')
      OUTPUT.write('Creating the directory: ELI.\n')
      OUTPUT.close()

      DMFT_Eliashberg(para_com_muk,initial,U,temperature,mu,mix,num_iter,nk,Nc,TB,Hopping,sigi)

      uptillnow = time()
      OUTPUT = open('OUTPUT', 'a') 
      OUTPUT.write('Time elapsed is: %.f s.\n'%(uptillnow - starting))
      OUTPUT.close()
 
   ### timing ###

   ending = time()
   end_time = strftime("%a, %d %b %Y %H:%M:%S", localtime())
 
   OUTPUT = open('OUTPUT', 'a')
   OUTPUT.write("\n\n===== End =====\n")
   OUTPUT.write("The calculation ends at %s.\n"%end_time)
   OUTPUT.write("Total elapsed time is: %.f s.\n"%(ending-starting))
   OUTPUT.close()


   ### prepare for the post-processing
   if os.path.exists('PP'):
      shutil.move('PP', 'PP.old')
   else:
      os.makedirs('PP')
      os.makedirs('PP/maxent')
      os.makedirs('PP/spectral')
      os.makedirs('PP/maxent/continue_Gf')
      os.makedirs('PP/maxent/continue_Sig')

   ### write the input file
   INPUT = open('PP/maxent/continue_Gf/INPUT', 'a')
   INPUT.write('### Input parameters for analytical continuation ###\n')
   INPUT.write('\n')
   INPUT.write('params = {"Nspin":      [1,                  "1: paramagnetic, 2: spin-polarized"],\n')
   INPUT.write('          "Choice":     [0,                  "0: G; 1: Sig; 2: G & Sig"],\n')
   INPUT.write('          "Width":      [0.01,               "# Smearing width"], \n')
   INPUT.write('          "N_tau":      [20000,              "number of tau"], \n')
   INPUT.write('          "Max_freq":   [30,                 "max frequency"] \n')
   INPUT.write(' }\n')
   INPUT.write('\n')
   INPUT.write('#for c++ pointer issue, the dictionary value must a constant rather than a list.\n')
   INPUT.write('params_maxent = {"Fit_tail":      False, \n')
   INPUT.write('                 "Nfit":          10,  \n')
   INPUT.write('                 "Logfile":       "log.dat", \n')
   INPUT.write('                 "Finegrid":      1, \n')
   INPUT.write('                 "Maxomega":      30., \n')
   INPUT.write('                 "Domega":        0.02, \n')
   INPUT.write('                 "Enforce_norm":  True, \n')
   INPUT.write('                 "Model_select":  2, \n')
   INPUT.write('                 "Alpha":         [5], \n')
   INPUT.write('                 "Compute_logPa": True,\n')
   INPUT.write('                 "Firststep":     3000,\n')
   INPUT.write('                 "Prefix":        "", \n')
   INPUT.write('                 "Output_dir":    "results", \n')
   INPUT.write('                 "Pnorm_tol":     0.01 \n')
   INPUT.write('}\n')
   INPUT.close()


   ### write the input file
   INPUT = open('PP/maxent/continue_Sig/INPUT', 'a')
   INPUT.write('### Input parameters for analytical continuation ###\n')
   INPUT.write('\n')
   INPUT.write('params = {"Nspin":      [1,                  "1: paramagnetic, 2: spin-polarized"],\n')
   INPUT.write('          "Choice":     [1,                  "0: G; 1: Sig; 2: G & Sig"],\n')
   INPUT.write('          "Width":      [0.01,               "# Smearing width"], \n')
   INPUT.write('          "N_tau":      [20000,              "number of tau"], \n')
   INPUT.write('          "Max_freq":   [30,                 "max frequency"] \n')
   INPUT.write(' }\n')
   INPUT.write('\n')
   INPUT.write('#for c++ pointer issue, the dictionary value must a constant rather than a list.\n')
   INPUT.write('params_maxent = {"Fit_tail":      False, \n')
   INPUT.write('                 "Nfit":          10,  \n')
   INPUT.write('                 "Logfile":       "log.dat", \n')
   INPUT.write('                 "Finegrid":      1, \n')
   INPUT.write('                 "Maxomega":      30., \n')
   INPUT.write('                 "Domega":        0.02, \n')
   INPUT.write('                 "Enforce_norm":  True, \n')
   INPUT.write('                 "Model_select":  2, \n')
   INPUT.write('                 "Alpha":         [5], \n')
   INPUT.write('                 "Compute_logPa": True,\n')
   INPUT.write('                 "Firststep":     3000,\n')
   INPUT.write('                 "Prefix":        "", \n')
   INPUT.write('                 "Output_dir":    "results", \n')
   INPUT.write('                 "Pnorm_tol":     0.01 \n')
   INPUT.write('}\n')
   INPUT.close()

   ### write an input file
   INPUT = open('PP/spectral/INPUT', 'a')
   INPUT.write('### Input parameters for analytical continuation ###\n')
   INPUT.write('\n')
   INPUT.write('params = {"nspin":       [1,                  "1: nsp, 2: spin-polarized"],\n')
   INPUT.write('          "choice":      [0,                  "0: ksum; 1: kcut, 2: kfermi"],\n')
   INPUT.write('          "q":           [[12, 12, 12],       "[Nq_x, Nq_y, Nq_z]"],\n')
   INPUT.write('          "qlist":       [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.5, 0.5], [0.0, 0.0, 0.0]],  "qlist for kcut"],\n')
   INPUT.write('          "nqgrid":      [20,                 "number of q grid for each cut"],\n')
   INPUT.write('          "mu":          [3.77,               "chemical potential"],\n')
   INPUT.write('          "n_tot":       [164,                "total number of electron"],\n')
   INPUT.write('          "delta":       [0.01,               "i*delta"],\n')
   INPUT.write('          "strength":    [1.0,                "strength"],\n')
   INPUT.write('          "kz":          [0.0,                "kz"],\n')
   INPUT.write('          "beta":        [40,                 "inverse temperature"],\n')
   INPUT.write('          "omega_range": [20.0,               "omega_range"]\n')   
   INPUT.write('}')

   ### copy file
   shutil.copy('GF_SIG/Gf_loc.out'+str(Niter), 'PP/maxent/continue_Gf/Gf_loc.out')
   shutil.copy('GF_SIG/Sig_loc.out'+str(Niter), 'PP/maxent/continue_Gf/Sig_loc.out')
   shutil.copy('GF_SIG/Gf_loc.out'+str(Niter), 'PP/maxent/continue_Sig/Gf_loc.out')
   shutil.copy('GF_SIG/Sig_loc.out'+str(Niter), 'PP/maxent/continue_Sig/Sig_loc.out')
   shutil.copy('RHAM', 'PP/spectral')
   shutil.copy('ORBITAL', 'PP/spectral')
