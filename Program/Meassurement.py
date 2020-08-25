import numpy as np
import time


class Correlation:
    def __init__(self, configurations, therm, iteration_len=10000, lat_len=3, delta=1, mean_len=100,  step=4, Nt=10):
        self.configurations = configurations
        self.therm = therm
        self.iteration_len = iteration_len
        self.lat_len = lat_len

        self.mean_len = mean_len
        self.step = step
        self.delta = delta

        self.Nt = Nt
        self.Ng = self.lat_len**2*self.Nt*2

    def phi2_correlation(self):
        y = []
        for j in range(1, self.iteration_len-4):
            mv_phi_sq = float(sum((self.configurations[j]) ** 2)+sum((self.configurations[j+1]) ** 2)+
                                   sum((self.configurations[j+2]) ** 2)+sum((self.configurations[j+3]) ** 2)+sum((self.configurations[j+4]) ** 2))
            y.append(mv_phi_sq/self.Ng**2/5)

        x = list(range(1, self.iteration_len-4))
        return y, x, 0

    def fit_func_phi(self, x, a, b):
        return 0.5*np.tanh((x-a)*b)+0.5

    def phi_auto_correlation(self):
        size = self.iteration_len - self.therm
        y = list()
        x = list()
        for N in range(1, int(size/2), self.step*30):
            mv_phi_cor = 0
            for j in range(size):
                mv_phi_cor += float(self.configurations[self.therm + j].T @ (self.configurations[((self.therm + j + N)%size)])) / (size * self.Ng ** 2)
            y.append(mv_phi_cor)
            x.append(N)

        return y, x

    def spin_correlation(self, axis):
        start = time.time()
        mv_spin_cor_list = list()
        a = np.arange(0, self.Ng)
        x = list()
        max = len(self.configurations)

        print('therm', self.therm)
        # Calculating MV of Spincorrelation of a given axis for each configuration
        for n in range(self.therm, max):
            spin_cor = sum(map(lambda z, y: np.tanh(self.configurations[n][z]) * np.tanh(self.configurations[n][y]), a, self.map_generator()[axis]))/self.Ng
            mv_spin_cor_list.append(spin_cor)
            x.append(n)

        # Calculate the Mean of the Spincorrelation with standard deviation
        mv_spin_cor = sum(mv_spin_cor_list) / len(mv_spin_cor_list)
        sd_spin_cor = np.sqrt(
            1 / ((max - self.therm) * self.Ng - 1) * sum((np.array(mv_spin_cor_list) - mv_spin_cor) ** 2))
        print('dim of mv spin', np.array(mv_spin_cor_list).shape)

        print('It took {0:0.01f} seconds to calculate spin correlation'.format(time.time() - start))

        return mv_spin_cor, sd_spin_cor, mv_spin_cor_list, x

    def spin_auto_cor(self, spin_cor_list, mv_spin_cor):
        auto_spin_cor_list = list()
        x = list()
        size = len(spin_cor_list)
        for N in range(1, int(size/2), self.step*30):
            auto_spin_cor = 0
            for j in range(size):
                auto_spin_cor += (spin_cor_list[j]-mv_spin_cor)*(spin_cor_list[(j+N)%size]-mv_spin_cor) / size
            auto_spin_cor_list.append(auto_spin_cor)
            x.append(N)
        try:
            autocor_len = list(map(lambda i: i < auto_spin_cor_list[0]/np.e, auto_spin_cor_list)).index(True)
        except ValueError:
            autocor_len = 0
        print('acceplist[corlen]', auto_spin_cor_list[autocor_len])
        return auto_spin_cor_list, x, autocor_len*120

    # Defining spin correlation mapping
    def map_generator(self):
        x = np.array([i + 1 if (i + 1) % self.lat_len != 0 else i + 1 - self.lat_len for i in range(self.Ng)])
        y = np.array([i + self.lat_len if ((i + self.lat_len) % (self.lat_len ** 2)) >= self.lat_len else i + self.lat_len - self.lat_len ** 2 for i in range(self.Ng)])
        tau = np.array([i+self.lat_len**2 if ((i+self.lat_len**2) % (self.lat_len**2*self.Nt)) >= self.lat_len**2 else i + self.lat_len**2 - self.lat_len**2*self.Nt for i in range(self.Ng)])
        dir = np.hstack((np.arange(int(self.Ng / 2), self.Ng), np.arange(0, int(self.Ng / 2))))
        map_list = [x, y, tau, dir]
        return map_list
