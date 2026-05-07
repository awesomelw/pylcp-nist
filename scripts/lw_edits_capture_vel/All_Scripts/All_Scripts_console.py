import numpy as np
import matplotlib.pyplot as plt

# Edits
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import odeint
from scipy.optimize import fsolve
import scipy.constants as cts
import pylcp
import pylcp.atom as atom
from pylcp.fields import conventional3DMOTBeams
from pylcp.common import bisectFindChangeValue, progressBar
#plt.style.use('paper')
import time


# 50 gauss/cm for Sr 
# mub on database
# 0.43 cm for x0 30/1.4/50

# x is length unit
# k is wavenumber which relates to wavelength
x0 = 0.1 # cm
k = 2*np.pi/461E-7 # cm^{-1}
kbar = k*x0

# D source for position
d = np.array([0., -8.839, -8.839])

# Gamma is decay rate
# t0 is normalized time units of decary
# wb is width factor
gamma = 2*np.pi*30e6
t0 = 1e-5 # s
gammabar = gamma*t0
alpha = 2*np.pi*cts.value('Bohr magneton in Hz/T')*50*1e-4*x0*t0

mbar = 87.9056122571*cts.value('atomic mass constant')*(x0*1e-2)**2/(cts.hbar*t0)

wb = 4.7 # mm
print(x0, t0, k, kbar, gamma, gammabar, mbar, wb, alpha, gammabar/kbar)


class aPHIMOT(pylcp.laserBeams):
     def __init__(self, *args, **kwargs):
        super().__init__()
    
        beam_type = kwargs.pop('beam_type', pylcp.laserBeam)
        pol = kwargs.pop('pol', +1)
        kmag = kwargs.pop('k', 1.)
        
        self.add_laser(beam_type(kmag*(np.array([0., 0.,  17.678])/(np.linalg.norm(np.array([0., 0.,  17.678])))), +pol, *args, **kwargs))
        self.add_laser(beam_type(kmag*(np.array([0., 0.,  -17.678])/(np.linalg.norm(np.array([0., 0.,  -17.678])))), +pol, *args, **kwargs))
        self.add_laser(beam_type(kmag*(np.array([10.825, -13.258, 4.419])/(np.linalg.norm(np.array([10.825, -13.258, 4.419])))), -pol, *args, **kwargs))
        self.add_laser(beam_type(kmag*(np.array([10.825, 13.258, -4.419])/(np.linalg.norm(np.array([10.825, 13.258, -4.419])))), -pol, *args, **kwargs))
        self.add_laser(beam_type(kmag*(np.array([-10.825, -13.258, 4.419])/(np.linalg.norm(np.array([-10.825, -13.258, 4.419])))), -pol, *args, **kwargs))
        self.add_laser(beam_type(kmag*(np.array([-10.825, 13.258, -4.419])/(np.linalg.norm(np.array([-10.825, 13.258, -4.419])))), -pol, *args, **kwargs))

#beam_to_sim = pylcp.infinitePlaneWaveBeam
#beam_to_sim = pylcp.gaussianBeam
#beam_to_sim = pylcp.clippedGaussianBeam
list_beam_to_sim = [pylcp.infinitePlaneWaveBeam, pylcp.gaussianBeam, pylcp.collimatedGaussianBeam]

#MOT_to_sim = conventional3DMOTBeams
#MOT_to_sim = aPHIMOT
list_MOT_to_sim = [conventional3DMOTBeams, aPHIMOT]

#MOT_to_sim_kwargs = {'rotation_angles':[np.pi/4, 0., 0.]} # Extra arguments for conventional3DMOTBeams
#MOT_to_sim_kwargs = {} # Extra arguments for aPHIMOT
list_MOT_to_sim_kwargs = [{'rotation_angles':[np.pi/4, 0., 0.]}, {}]

#laser_kwargs = {}
#laser_kwargs = {'wb':wb} # Extra arguments for GaussianBeam
#laser_kwargs = {'wb':1000*wb, 'rs':wb} # Extra arguments for clippedGaussianBeam
list_laser_kwargs = [{}, {'wb':wb}, {'wb':1000*wb, 'rs':wb}]

fileName1 = ""
fileName2 = ""

for g in range(len(list_MOT_to_sim)):
    if g==0:
        fileName1 = "conventional"
    elif g==1:
        fileName1 = "APHI"
        
    for b in range(len(list_beam_to_sim)):
        if b==0:
            fileName2 = "infinite"
        elif b==1:
            fileName2 = "gauss"
        elif b==2:
            fileName2 = "collimatedgauss"
        
        # Set up the laser beams with their appropriate characteristics
        # det is detuning of frequency
        # alpha relates to magnetic field gradient
        # beta is normalized intensity
        det = -1.5*gammabar
        beta = 1.0
        
        laserBeams = list_MOT_to_sim[g](beta=beta, delta=det, beam_type=list_beam_to_sim[b], k=kbar, **list_laser_kwargs[b], **list_MOT_to_sim_kwargs[g]) # kbar
        
        #print(laserBeams.pol())

        if (np.sum(np.isnan(laserBeams.pol()), axis=(0, 1)) > 0):
            del laserBeams
            laserBeams = list_MOT_to_sim[g](beta=beta, delta=det, beam_type=list_beam_to_sim[b], k=kbar, **list_laser_kwargs[b], **list_MOT_to_sim_kwargs[g])

        magField = pylcp.quadrupoleMagneticField(alpha)

        equation = pylcp.heuristiceq(laserBeams, magField, gamma=gammabar, mass=mbar)

#         Hg, mugq = pylcp.hamiltonians.singleF(F=0, muB=1)
#         He, mueq = pylcp.hamiltonians.singleF(F=1, muB=1)

#         dq = pylcp.hamiltonians.dqij_two_bare_hyperfine(0, 1)

#         hamiltonian = pylcp.hamiltonian(mass=mbar)
#         hamiltonian.add_H_0_block('g', Hg)
#         hamiltonian.add_H_0_block('e', He)
#         hamiltonian.add_mu_q_block('g', mugq)
#         hamiltonian.add_mu_q_block('e', mueq)
#         hamiltonian.add_d_q_block('g','e', dq, gamma=gammabar, k=kbar)

#         equation = pylcp.rateeq(laserBeams, magField, hamiltonian, include_mag_forces=False)
#         equation.set_initial_pop(np.array([1.,0,0,0]))

        # Create the forces
        #dz = 0.025
        #dv = 0.25 # Like pixels, if you increase this it runs faster
        z = np.linspace(-10, 10, 101) # mm
        v = np.linspace(-1, 1, 101) # mm/10 microsec

        dz = np.mean(np.diff(z))
        dv = np.mean(np.diff(v))

        Z, V = np.meshgrid(z, v)

        # Rfull = np.array([np.zeros(Z.shape), np.zeros(Z.shape), Z]) Changed when X = Y = 0
        # Vfull = np.array([np.zeros(Z.shape), np.zeros(Z.shape), V])

        # Computing force along this line
        Rfull = np.array([np.zeros(Z.shape), np.zeros(Z.shape), Z])
        Vfull = np.array([np.zeros(Z.shape), np.zeros(Z.shape), V]) 

        equation.generate_force_profile(Rfull, Vfull, name='Fz', progress_bar=True)


        # Define when the program stops and also describes the trajectories. Sets the initial motion/velocity and then evolves the motion

        v0s = np.arange(0.1, 1., 0.1)

        # See solve_ivp documentation for event function discussion:
        def captured_condition(t, y, threshold=0.1):
            if np.linalg.norm(y[-6:-3])<threshold  and np.linalg.norm(y[-3:])<0.1:
                val = -1.
            else:
                val = 1.

            return val

        def lost_condition(t, y):
            if np.linalg.norm(y[-3:])>15.:
                val = -1.
            else:
                val = 1.

            return val

        def nan_found(t, y):
            if np.sum(np.isnan(y)) > 0:
                return -1.
            else:
                return 1.

        captured_condition.terminal=True
        lost_condition.terminal=True
        nan_found.terminal=True

        sols = []

        if b==2:
            for v0 in v0s:
                equation.set_initial_position_and_velocity(d, -1*v0*d/np.linalg.norm(d))
            #    equation.set_initial_position_and_velocity(np.array([0,0,z[0]]), np.array([0,0,v0]))
                #equation.set_initial_pop(np.array([1.,0,0,0]))
                equation.evolve_motion([0., 1000.], events=[captured_condition, lost_condition, nan_found], max_step=1., progress_bar=True)
                sols.append(equation.sol)
        else:
            for v0 in v0s:
                equation.set_initial_position_and_velocity(d, -1*v0*d/np.linalg.norm(d))
            #    equation.set_initial_position_and_velocity(np.array([0,0,z[0]]), np.array([0,0,v0]))
                #equation.set_initial_pop(np.array([1.,0,0,0]))
                equation.evolve_motion([0., 1000.], events=[captured_condition, lost_condition, nan_found], progress_bar=True)
                sols.append(equation.sol)


        # Phase diagrams are only useful right now in certain axises
        fig, ax = plt.subplots(1, 1, figsize=(4, 2.75))
        plt.imshow(equation.profile['Fz'].F[2]/kbar/gammabar, origin='bottom',
                   extent=(np.amin(z)-dz/2, np.amax(z)+dz/2,
                           np.amin(v)-dv/2, np.amax(v)+dv/2), aspect='auto')
        cb1 = plt.colorbar()
        cb1.set_label('$F/(\hbar k \Gamma)$')
        ax.set_xlabel('$z$ (mm)')
        ax.set_ylabel('$v$ (mm/$10$ $\mu s$)')
        fig.subplots_adjust(left=0.12,right=0.9)

        for sol in sols:
            ax.plot(sol.r[2], sol.v[2], 'w-', linewidth=0.375)

        #plt.show()
        plt.savefig(fileName1 + fileName2 + "1" + ".png")

        # Phase diagrams are only useful right now in certain axises
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # sol.r[0] is x and so on
        # You can use dir to debug
        for sol in sols:
            ax.plot(sol.r[0], sol.r[1], sol.r[2], linewidth=0.375)

        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')  
        ax.set_zlabel('$z$')

        # Set axises
        # lims = np.array([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
        # bottom = np.amin(lims[:,0])
        # top = np.amax(lims[:,1])
        bottom = -1.5
        top = 1.5      

        ax.set_xlim(bottom, top)
        ax.set_ylim(bottom, top)
        ax.set_zlim(bottom, top)

        #plt.show()
        plt.savefig(fileName1 + fileName2 + "2" + ".png")

        # Phase diagrams are only useful right now in certain axises
        fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.75))

        for sol in sols:
            ax[0].plot(sol.t, sol.r[0], linewidth=0.375)
            ax[1].plot(sol.t, sol.r[1], linewidth=0.375)
            ax[2].plot(sol.t, sol.r[2], linewidth=0.375)

        #plt.show()
        plt.savefig(fileName1 + fileName2 + "3" + ".png")

        # Phase diagrams are only useful right now in certain axises
        fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.75))

        for sol in sols:
            ax[0].plot(sol.t, sol.v[0], linewidth=0.375)
            ax[1].plot(sol.t, sol.v[1], linewidth=0.375)
            ax[2].plot(sol.t, sol.v[2], linewidth=0.375)

        #plt.show()
        plt.savefig(fileName1 + fileName2 + "4" + ".png")

        for sol in sols:
            if len(sol.t_events[0]) == 1:
                print('captured')
            elif len(sol.t_events[1]) == 1:
                print('lost')

        def iscaptured(v0, r0, eqn, captured_condition, lost_condition, nan_found, tmax=1000, **kwargs):
            eqn.set_initial_position_and_velocity(r0, -1*v0*r0/np.linalg.norm(r0))
            eqn.evolve_motion([0., tmax], events=[captured_condition, lost_condition, nan_found],
                              **kwargs)

            return len(eqn.sol.t_events[0]) == 1

#         iscaptured(0.1, d, equation, captured_condition, lost_condition, nan_found, tmax=1000)



#         start_time = time.time()

        bisectFindChangeValue(iscaptured, 0.1,
                              args=(d, equation, captured_condition, lost_condition, nan_found),
                              kwargs={'tmax':10000},
                              tol=1e-1
                             )

#         print("--- %s seconds ---" % (time.time() - start_time))



#         dets = -gammabar*np.logspace(-0.6, np.log10(3), 10)[::-1]
        dets = gammabar*np.linspace(-3, -0.25, 41) #41
        betas = np.array([0.5,1.,2.])

        DETS, BETAS = np.meshgrid(dets, betas)

        it = np.nditer([DETS, BETAS, None, None])

        progress = progressBar()
        for (det, beta, vc, iterations) in it:
        #     #beam_to_sim = pylcp.infinitePlaneWaveBeam
        #     #beam_to_sim = pylcp.gaussianBeam
        #     #beam_to_sim = pylcp.clippedGaussianBeam
        #     list_beam_to_sim = [pylcp.infinitePlaneWaveBeam, pylcp.gaussianBeam, pylcp.clippedGaussianBeam]

        #     #MOT_to_sim = conventional3DMOTBeams
        #     #MOT_to_sim = aPHIMOT
        #     list_MOT_to_sim = [conventional3DMOTBeams, aPHIMOT]

        #     #MOT_to_sim_kwargs = {'rotation_angles':[np.pi/4, 0., 0.]} # Extra arguments for conventional3DMOTBeams
        #     #MOT_to_sim_kwargs = {} # Extra arguments for aPHIMOT
        #     list_MOT_to_sim_kwargs = [{'rotation_angles':[np.pi/4, 0., 0.]}, {}]

        #     #laser_kwargs = {}
        #     #laser_kwargs = {'wb':wb} # Extra arguments for GaussianBeam
        #     #laser_kwargs = {'wb':1000*wb, 'rs':wb} # Extra arguments for clippedGaussianBeam
        #     list_laser_kwargs = [{}, {'wb':wb}, {'wb':1000*wb, 'rs':wb}]

        #     for g in range(len(list_MOT_to_sim)):
        #         for b in range(len(list_beam_to_sim)):
            laserBeams = list_MOT_to_sim[g](beta=beta, delta=det, beam_type=list_beam_to_sim[b], k=kbar, **list_laser_kwargs[b], **list_MOT_to_sim_kwargs[g]) # kbar

            #print(laserBeams.pol())

            if (np.sum(np.isnan(laserBeams.pol()), axis=(0, 1)) > 0):
                del laserBeams
                laserBeams = list_MOT_to_sim[g](beta=beta, delta=det, beam_type=list_beam_to_sim[b], k=kbar, **list_laser_kwargs[b], **list_MOT_to_sim_kwargs[g])

            equation = pylcp.heuristiceq(laserBeams, magField, gamma=gammabar, mass=mbar)

#             Hg, mugq = pylcp.hamiltonians.singleF(F=0, muB=1)
#             He, mueq = pylcp.hamiltonians.singleF(F=1, muB=1)

#             dq = pylcp.hamiltonians.dqij_two_bare_hyperfine(0, 1)

#             hamiltonian = pylcp.hamiltonian(mass=mbar)
#             hamiltonian.add_H_0_block('g', Hg)
#             hamiltonian.add_H_0_block('e', He)
#             hamiltonian.add_mu_q_block('g', mugq)
#             hamiltonian.add_mu_q_block('e', mueq)
#             hamiltonian.add_d_q_block('g','e', dq, gamma=gammabar, k=kbar)

#             equation = pylcp.rateeq(laserBeams, magField, hamiltonian, include_mag_forces=False)
#             equation.set_initial_pop(np.array([1.,0,0,0]))
            
            if b==2.:
                vc[...], iterations[...] = bisectFindChangeValue(
                iscaptured, 0.1,
                args=(d, equation, captured_condition, lost_condition, nan_found),
                kwargs={'tmax':10000, 'max_step':1.},
                tol=1e-1, debug=False
                )
            else:
                vc[...], iterations[...] = bisectFindChangeValue(
                iscaptured, 0.1,
                args=(d, equation, captured_condition, lost_condition, nan_found),
                kwargs={'tmax':10000},
                tol=1e-1, debug=False
                )

            progress.update((it.iterindex+1)/it.itersize)

        def vc_from_paper(delta, beta, mbar):
            return 1/(2*mbar)**(2./3.)*(beta**2/(1+beta)**(3./2.))**(1./3.)*(8*np.pi*delta**2/(1+beta+4*delta**2))**(1./3.)

        fig, ax = plt.subplots(1, 1)
        fig, ax = plt.subplots(1, 1)
        for (beta, vc_vs_det) in zip(betas, it.operands[2]):
            ax.plot(dets/gammabar, vc_vs_det, label='$\\beta=%.1f$' % beta)
        ax.legend(fontsize=8)
        ax.set_xlabel('$\Delta/\Gamma$')
        ax.set_ylabel('$v_c$ ($10^2$ m/s)')

        #plt.show()
        plt.savefig(fileName1 + fileName2 + "5" + ".png")
        print("------------------------------------------------------------------------------")