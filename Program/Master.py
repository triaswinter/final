import Simulation as SIM
import Meassurement as MES
import Matrix_generator as MTX
from multiprocessing.dummy import Pool as ThreadPool
from itertools import product
import mysql.connector

mydb = mysql.connector.connect(
    host="localhost",
    user="root",
    passwd="zoomg21u",
    database="hmcdb"
)

mycurser = mydb.cursor()

sqlAddData = "INSERT INTO hmcdb (J, Nt, dir-correlation, dc-error) VALUES (%s, %s, %s, %s)"


# Global Variables
iteration_len1 = 2000
iteration_len2 = 300
therm = 30
leapfrog_len = 2
lat_len = 3


def run(J, Nt):
    # Create Data
    print('0')
    hmc, hmc_delta = initialize('A', J, Nt)
    print('1')
    delta, delta_pd = find_delta(hmc_delta)
    print(J,Nt, '2: here i am')
    configurations = hmc.create_sample_space(delta)
    cor = MES.Correlation(configurations, therm, iteration_len1, lat_len, delta, Nt=Nt)
    print(J,Nt, 'rock you')

    # Analyse Data
    phi_pd1, phi_pd2 = phi_correlation(cor)
    cor_pd1, cor_pd2, mv_spin_correlation_dir, dc_error = spin_correlation(cor)
    print(J,Nt, 'like a')
    # Maybe Safe into SQL

    result_tupel = (J, Nt, mv_spin_correlation_dir, dc_error)
    print(J,Nt, 'Hurricane')
    return result_tupel


def who_knows():  # Run simulation and Measurements
    delta, delta_pd = find_delta()
    self.config_gen(delta)

    phi_pd1, phi_pd2 = self.phi_correlation()
    cor_pd1, cor_pd2, mv_spin_correlation_dir = self.spin_correlation()

    plot_list = [delta_pd, phi_pd1, phi_pd2, cor_pd1, cor_pd2]
    parameters = self.parameters
    self.SAVE_ARRAY.append([*plot_list, mv_spin_correlation_dir])

    # df = pd.DataFrame([[11, 21, 31], [12, 22, 32], [31, 32, 33]], index=['one', 'two', 'three'], columns=['a', 'b', 'c'])

    if self.enable_plot.get() == 1:
        self.plot(plot_list)

    self.counter += 1


def initialize(method, J, Nt):  # Set Parameters
    matrix = create_matrix(method, J, Nt)
    hmc = SIM.HMC(method, matrix, iteration_len1, leapfrog_len, lat_len, Nt)
    hmc_delta = SIM.HMC(method, matrix, iteration_len2, leapfrog_len, lat_len, Nt)
    return hmc, hmc_delta


def create_matrix(method, J, Nt):
    if method == 'A_inv':
        matrix = MTX.create_inv(lat_len, para1=J, para2=Nt)
    elif method == 'A':
        A, matrix = MTX.create_a(lat_len, para1=J, para2=Nt)
    return matrix


def find_delta(hmc_delta):
    # Find delta for the simulation with the given parameters (provided that they are proper)
    y = hmc_delta.find_delta()
    print('find delta is no problem')
    x_max = len(y) - 1
    x, delta, fit, error = hmc_delta.deltaplot_data(y)
    print('deltaplot data neigher')
    if error is True:
        print('error')


    # Generate delta plot
    labels = ['Acceptance rate scaling to ' + r'$\delta$', r'$\delta$', 'Acceptancerate']
    print('gui lenx leny', len(x[:x_max]), len(y[:x_max]))
    delta_plot_data = [x[:x_max], y[:x_max], labels, fit[:x_max]]

    return delta, delta_plot_data


'''def config_gen(delta):
    # Get Parameters
    lat_len = getattr(self.hmc, 'lat_len')
    Nt = getattr(self.hmc, 'Nt')
    configuration = self.hmc.create_sample_space(delta)

    # Create a corralation objekt
    self.cor = SIM.Correlation(configuration, therm, iteration_len1, lat_len, delta, Nt=Nt)

    self.results.insert('end', "Configuration is generated")
    self.results.see('end')'''


'''def plot(plot_list):
    for i in range(5):
        fig = SIM.add_subplot(self.fig_list[i], self.counter, self.SAVE_ARRAY[self.counter - 1][i])
        canvas = FigureCanvasTkAgg(fig, master=self.fig_canvas_list[i])
        canvas.draw()
        size = fig.get_size_inches() * fig.dpi
        self.fig_canvas_list[i].create_window((size[0] // 2, size[1] // 2), window=canvas.get_tk_widget())
        self.fig_canvas_list[i].configure(scrollregion=self.fig_canvas_list[i].bbox("all"))'''


# Measurements
def phi_correlation(cor):
    y1, x1 = cor.phi_auto_correlation()
    y2, x2, therm = cor.phi2_correlation()

    label_1 = ['Correlation <Phi_0*Phi_N>', 'N', '<Phi_0*Phi_N>']
    label_2 = ['Meanvalue of Phi² in dependence of N', 'N', '<Phi²>']
    phi_pd1 = [x1, y1, label_1, []]
    phi_pd2 = [x2, y2, label_2, []]

    return phi_pd1, phi_pd2


def spin_correlation(cor):
    mv_spin_cor, sd_spin_cor, spin_cor_list, x1 = cor.spin_correlation(3)
    spin_auto_cor_list, x2, autocor_len = cor.spin_auto_cor(spin_cor_list, mv_spin_cor)

    labels_z = ['Spin correlation mean value', 'N', '<S>']
    labels_y = ['Spin autocorrelation mean value', 'N', '<S_0*S_N>']

    print('print the correlation', spin_cor_list)

    cor_pd1 = [x1, spin_cor_list, labels_z, []]
    cor_pd2 = [x2, spin_auto_cor_list, labels_y, []]

    return cor_pd1, cor_pd2, mv_spin_cor, sd_spin_cor


J_list = (0.5,2)
Nt_list = (3,10)
input_para = product(J_list, Nt_list)
pool = ThreadPool(4)

results = pool.starmap(run, list(input_para))

print(results)