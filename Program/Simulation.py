import numpy as np
import time
from scipy.optimize import curve_fit


class HMC:
    def __init__(self, method, matrix, iteration_len, leapfrog_len, lat_len=3, Nt=1):
        self.iteration_len = iteration_len
        self.leapfrog_len = leapfrog_len
        self.lat_len = lat_len
        self.method = method
        self.matrix = matrix

        self.Nt = int(Nt)

        self.Ng = lat_len**2*self.Nt*2
        self.delta_range = np.append(np.arange(1/1000, self.leapfrog_len/100, 1/200), np.arange(self.leapfrog_len/100, self.leapfrog_len/50, 1/100))
        self.delta_range = np.append(self.delta_range, np.arange(self.leapfrog_len / 50, self.leapfrog_len / 2, 1 / 70))

    # Simulation
    def get_Nt(self):
        return self.Nt

    def find_delta(self):
        start = time.time()
        acceptance_rate_list = list()

        for delta in self.delta_range:
            # Set initial Configuration
            acceptance_rate = 0
            phi = np.random.normal(0, 1, size=(self.Ng, 1))
            if self.method == 'A':
                phi = self.matrix @ (phi)
            if self.method == 'A_FT':
                phi = np.fft.fft(phi).real / np.sqrt(self.Ng)

            for t in range(self.iteration_len):
                # Define phi and momentum at t
                phi_0 = phi
                momentum = np.random.normal(0, 2, size=(self.Ng, 1))

                # Calculate trajectory
                phi, momentum, tot_energy_1, tot_energy_2 = self.leapfrog(phi, momentum, delta)

                # Metropolis check
                r = float(np.random.random())
                energy_dif = abs(tot_energy_2 - tot_energy_1)
                if r > np.exp(-energy_dif):
                    phi = phi_0
                else:
                    acceptance_rate += 1

            # Produce plotdata
            acceptance_rate = acceptance_rate / self.iteration_len
            acceptance_rate_list.append(acceptance_rate)
            if acceptance_rate <= 0.01:
                acceptance_rate_list.extend([0]*60)
                print('low accep rat acceplist', acceptance_rate_list)
                return acceptance_rate_list
        print('It took {0:0.01f} seconds to find delta'.format(time.time() - start))

        return acceptance_rate_list

    def create_sample_space(self, delta):
        start = time.time()
        a = 0
        # Set initial Configuration
        phi = np.random.normal(0, 1, size=(self.Ng, 1))
        if self.method == 'A_FT':
            configurations = phi.T
            phi = np.fft.fft(phi).real/np.sqrt(self.Ng)
        else:
            configurations = np.array(phi.T)

        for t in range(self.iteration_len):
            # Save Configuration at t
            if self.method == 'A_inv':
                configurations = np.vstack((configurations, phi.T))
            elif self.method == 'A':
                u = (self.matrix @ phi).T
                configurations = np.vstack((configurations, u))
            elif self.method == 'A_FT':
                u = (np.fft.ifft(phi)).T.real/np.sqrt(self.Ng)
                configurations = np.vstack((configurations, u))

            # Define phi and momentum at t
            phi_0 = phi
            momentum = np.random.normal(0, 1, size=(self.Ng, 1))
            # Calculate trajectory
            phi, momentum, tot_energy_1, tot_energy_2 = self.leapfrog(phi, momentum, delta)

            # Metropolis check
            r = float(np.random.random())
            energy_dif = abs(tot_energy_2 - tot_energy_1)
            if r > np.exp(-energy_dif):
                phi = phi_0
            else:
                a +=1
        print('Configuration is generated', a/self.iteration_len)
        print('It took {0:0.01f} seconds to generate the configuration'.format(time.time() - start))
        return configurations

    def leapfrog(self, phi, momentum, delta):
        # Calculate the energy at t
        tot_energy_1 = self.energy(phi, momentum)
        # Calculate Trajectory
        i = delta
        while i <= self.leapfrog_len:
            phi = phi + (delta / 2) * momentum
            if self.method == 'A_inv':
                momentum = momentum - delta * (self.matrix @ phi - np.tanh(phi))
            elif self.method == 'A':
                momentum = momentum - delta * (self.matrix @ phi - self.matrix @ (np.tanh(self.matrix @ phi)))

            phi = phi + (delta / 2) * momentum
            i += delta
        # Calculate the energy at t+1
        tot_energy_2 = self.energy(phi, momentum)

        return phi, momentum, tot_energy_1, tot_energy_2

    def energy(self, phi, momentum):
        kin_energy = sum(momentum ** 2) / 2
        if self.method == 'A_inv':
            pot_energy = np.transpose(phi) @ (self.matrix @ phi) / 2 - sum(np.log(np.cosh(phi)))
        elif self.method == 'A':
            pot_energy = np.transpose(phi) @ (self.matrix @ phi) / 2 - sum(np.log(np.cosh(self.matrix @ phi)))
        tot_energy = complex(kin_energy + pot_energy).real
        return tot_energy

    def deltaplot_data(self, y2):
        # Generate curve fit
        print('acceptancerate', len(y2), y2)
        try:
            fit_para, fitter = curve_fit(self.fit_func_delta, self.delta_range[:len(y2)], y2)
            delta_fit = float(fit_para[0])
            fit = self.fit_func_delta(self.delta_range, *fit_para)
        except ValueError:
            delta_fit = 0.015
            fit = []
        except RuntimeError:
            delta_fit = 0.01
            fit = []
        print('sim lenx leny', len(self.delta_range))
        return self.delta_range, abs(delta_fit), fit, False

    def fit_func_delta(self, x, delta_0, W):
        return 0.5*np.tanh(-(x-delta_0)*W)+0.5

    def find_max_index(self, data, value, a, b):
        i=0
        for x in data:
            if self.fit_func_delta(x, a, b) <= value:
                return i
            i += 1
        return 0
