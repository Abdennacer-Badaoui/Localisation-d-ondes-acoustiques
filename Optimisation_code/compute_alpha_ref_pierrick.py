# -*- coding: utf-8 -*-


# Python packages
import matplotlib.pyplot
import numpy
import scipy
from scipy.optimize import minimize
import scipy.io
from math import exp
from math import sin


def real_to_complex(z):
    return z[0] + 1j * z[1]


def complex_to_real(z):
    return numpy.array([numpy.real(z), numpy.imag(z)])


class Memoize:
    def __init__(self, f):
        self.f = f
        self.memo = {}

    def __call__(self, *args):
        if args not in self.memo:
            self.memo[args] = self.f(*args)
        # .. todo: deepcopy here if returning objects
        return self.memo[args]


""" Définition des fonctions étudié, fft_etudie correspond à g_1, g_k_doppler correspond à g_2 et g_3 est la dernière fonction"""


def id(k, omega):
    if k == 0:
        return 1
    else:
        return 0


def compute_alpha(omega, material='Birch_LT', fft_etudier=id):
    """
    .. warning: $w = 2 \pi f$
    w is called circular frequency
    f is called frequency
    """
    """ On a entré la fonction étudié en paramètre de la fonction pour plus de reproductibilité"""
    # parameters of the material
    # * Ici On suppose que la condition de neumann y'a un 0 à côté
    phi = 0.0
    gamma_p = 0.0
    sigma = 0.0
    rho_0 = 0.0
    alpha_h = 0.0
    c_0 = 0.0

    # Birch LT

    gamma_p = 7.0 / 5.0

    rho_0 = 1.2

    c_0 = 340.0
    if material == 'Birch_LT':
        phi = 0.529  # porosity
        sigma = 151429.0  # resitivity
        alpha_h = 1.37  # tortuosity
    if material == 'beton armé':
        phi = 0.38
        alpha_h = 1
        sigma = 2206.63
    if material == 'fibre de carbone':
        phi = 0.508
        alpha_h = 2.23
        sigma = 3700
        rho_0 = 2000
    if material == 'ISOREL':
        phi = 0.7
        alpha_h = 1.15
        sigma = 142300
    if material == 'Liner 1':
        phi = 0.8
        alpha_h = 0.98
        sigma = 50000

    # parameters of the geometry
    L = 0.01

    # parameters of the mesh
    resolution = 10  # := number of elements along L

    # parameters of the material (cont.)
    mu_0 = 1.0
    ksi_0 = 1.0 / (c_0 ** 2)
    mu_1 = phi / alpha_h
    ksi_1 = phi * gamma_p / (c_0 ** 2)
    a = sigma * (phi ** 2) * gamma_p / ((c_0 ** 2) * rho_0 * alpha_h)

    ksi_volume = phi * gamma_p / (c_0 ** 2)
    a_volume = sigma * (phi ** 2) * gamma_p / ((c_0 ** 2) * rho_0 * alpha_h)
    mu_volume = phi / alpha_h
    k2_volume = (1.0 / mu_volume) * ((omega ** 2) / (c_0 ** 2)) * \
        (ksi_volume + 1j * a_volume / omega)
    # print(k2_volume)

    # parameters of the objective function
    A = 1.0
    B = 1.0
    # if fft_etudier == id:
    #print(" la fft_étudié n'est pas initialisée ---- ")

    # defining k, omega and alpha dependant parameters' functions
    @Memoize
    def lambda_0(k, omega):  # * Lambda 0
        if k ** 2 >= (omega ** 2) * ksi_0 / mu_0:
            return numpy.sqrt(k ** 2 - (omega ** 2) * ksi_0 / mu_0)
        else:
            return numpy.sqrt((omega ** 2) * ksi_0 / mu_0 - k ** 2) * 1j

    @Memoize
    def lambda_1(k, omega):  # * Lambda 1
        temp1 = (omega ** 2) * ksi_1 / mu_1
        temp2 = numpy.sqrt((k ** 2 - temp1) ** 2 + (a * omega / mu_1) ** 2)
        real = (1.0 / numpy.sqrt(2.0)) * numpy.sqrt(k ** 2 - temp1 + temp2)
        im = (-1.0 / numpy.sqrt(2.0)) * numpy.sqrt(temp1 - k ** 2 + temp2)
        return complex(real, im)

    @Memoize
    def g(y, w):  # ! SIGNAL PHYSIQUE
        # * La fonction de numpy qu'on évalue
        # cas simple
        return y*exp(-w)

    """def g_fft(array,n=None,axis=-1,norm=None):
        # array = g. evalué en k fois
        return np.fft.ftt(array,n,axis,norm) #! Devrait remplacer g_k ??
    """
    @Memoize
    def g_k(k, omega):
        # * Enfaite la seule variable qu'on va faire varier c'est k
        # ! Transformée de fourier par rapport à y,w
        # * C'est la futur fonction de fourier qu'on va devoir changer
        return fft_etudier(k, omega)  # ! Ca change que ici

    @Memoize
    def f(x, k):  # * La fonction f du truc
        return ((lambda_0(k, omega) * mu_0 - x) * numpy.exp(-lambda_0(k, omega) * L)
                + (lambda_0(k, omega) * mu_0 + x) * numpy.exp(lambda_0(k, omega) * L))

    @Memoize
    def chi(k, alpha, omega):  # * la fonction xi du truc
        return (g_k(k, omega) * ((lambda_0(k, omega) * mu_0 - lambda_1(k, omega) * mu_1)
                                 / f(lambda_1(k, omega) * mu_1, k) - (lambda_0(k, omega) * mu_0 - alpha) / f(alpha, k)))

    @Memoize
    def eta(k, alpha, omega):
        return (g_k(k, omega) * ((lambda_0(k, omega) * mu_0 + lambda_1(k, omega) * mu_1)
                                 / f(lambda_1(k, omega) * mu_1, k) - (lambda_0(k, omega) * mu_0 + alpha) / f(alpha, k)))

    @Memoize
    def e_k(k, alpha, omega):  # * Les e_k souhaités
        expm = numpy.exp(-2.0 * lambda_0(k, omega) * L)
        expp = numpy.exp(+2.0 * lambda_0(k, omega) * L)

        if k ** 2 >= (omega ** 2) * ksi_0 / mu_0:
            return ((A + B * (numpy.abs(k) ** 2))
                    * (
                (1.0 / (2.0 * lambda_0(k, omega)))
                * ((numpy.abs(chi(k, alpha, omega)) ** 2) * (1.0 - expm)
                   + (numpy.abs(eta(k, alpha, omega)) ** 2) * (expp - 1.0))
                + 2 * L * numpy.real(chi(k, alpha, omega) * numpy.conj(eta(k, alpha, omega))))
                + B * numpy.abs(lambda_0(k, omega)) / 2.0 * ((numpy.abs(chi(k, alpha, omega)) ** 2) * (1.0 - expm)
                                                             + (numpy.abs(eta(k, alpha, omega)) ** 2) * (
                    expp - 1.0))
                - 2 * B * (lambda_0(k, omega) ** 2) * L * numpy.real(
                        chi(k, alpha, omega) * numpy.conj(eta(k, alpha, omega))))
        else:
            return ((A + B * (numpy.abs(k) ** 2)) * (L
                                                     * ((numpy.abs(chi(k, alpha, omega)) ** 2) + (
                                                         numpy.abs(eta(k, alpha, omega)) ** 2))
                                                     + complex(0.0, 1.0) * (1.0 / lambda_0(k, omega)) * numpy.imag(
                                                         chi(k, alpha, omega) * numpy.conj(eta(k, alpha, omega)
                                                                                           * (1.0 - expm))))) + B * L * (
                numpy.abs(lambda_0(k, omega)) ** 2) \
                * ((numpy.abs(chi(k, alpha, omega)) ** 2) + (numpy.abs(eta(k, alpha, omega)) ** 2)) \
                + complex(0.0, 1.0) * B * lambda_0(k, omega) * numpy.imag(
                chi(k, alpha, omega) * numpy.conj(eta(k, alpha, omega)
                                                  * (1.0 - expm)))

    @Memoize
    def sum_e_k(omega):  # * Renvoie la somme des e_k
        def sum_func(alpha):
            s = 0.0
            for n in range(-resolution, resolution + 1):
                k = n * numpy.pi / L
                s += e_k(k, alpha, omega)
            return s

        return sum_func

    @Memoize
    def alpha(omega):  # * Renvoie le Alpha souhaité
        alpha_0 = numpy.array(complex(40.0, -40.0))
        temp = real_to_complex(minimize(lambda z: numpy.real(sum_e_k(omega)(
            real_to_complex(z))), complex_to_real(alpha_0), tol=1e-4).x)
        #print(temp, "------", "je suis temp")
        return temp

    @Memoize
    def error(alpha, omega):
        temp = numpy.real(sum_e_k(omega)(alpha))
        return temp

    temp_alpha = alpha(omega)
    #temp_error = error(temp_alpha, omega)

    return temp_alpha


def run_compute_alpha(material, fun=id):
    """ Compute alpha, on a ajouté aux mtx les matériaux pour plus efficaces"""
    print('Computing alpha...')
    numb_omega = 600  # 1000
    # omegas = numpy.logspace(numpy.log10(600), numpy.log10(30000), num=numb_omega)
    omegas = numpy.linspace(2.0 * numpy.pi, numpy.pi * 10000, num=numb_omega)
    temp = [compute_alpha(omega, material, fun) for omega in omegas]
    print("temp:", "------", temp)
    alphas, errors = map(list, zip(*temp))
    alphas = numpy.array(alphas)
    errors = numpy.array(errors)

    print('Writing alpha...')
    output_filename = 'dta_omega_' + str(material) + namefun + '.mtx'
    scipy.io.mmwrite(output_filename, omegas.reshape(
        alphas.shape[0], 1), field='complex', symmetry='general')
    output_filename = 'dta_alpha_' + str(material) + namefun + '.mtx'
    scipy.io.mmwrite(output_filename, alphas.reshape(
        alphas.shape[0], 1), field='complex', symmetry='general')
    output_filename = 'dta_error_' + str(material) + namefun + '.mtx'
    scipy.io.mmwrite(output_filename, errors.reshape(
        errors.shape[0], 1), field='complex', symmetry='general')

    return


def run_plot_alpha(lmaterial):
    color = 'darkblue'
    color = {}
    """ Couleurs des différents matériaux
    On entre la fonction en paramètre pour la reproductibilité"""
    color['beton'] = 'indigo'
    color['Birch_LT'] = 'darkblue'
    color['pierre'] = 'y'
    # * liste des material
    omega_liste = {}
    alphas_liste = {}
    for material in lmaterial:
        # * Récupère les infos sur les datas
        input_filename = 'dta_omega_' + str(material) + namefun + '.mtx'
        omega_liste[material] = scipy.io.mmread(input_filename)
        omega_liste[material] = omega_liste[material].reshape(
            omega_liste[material].shape[0])
    print('Reading alpha...')
    # print(" on verifie pour omega",
    # omega_liste['beton'] == omega_liste['Birch_LT'])  # ! C'est normal?
    # input_filename = 'dta_omega_' + str(material) + '.mtx'
    # omegas = scipy.io.mmread(input_filename)
    # omegas = omegas.reshape(omegas.shape[0])
    for material in lmaterial:
        # * Récupère pour alphas les datas
        input_filename = 'dta_alpha_' + str(material) + namefun + '.mtx'
        alphas_liste[material] = scipy.io.mmread(input_filename)
        alphas_liste[material] = alphas_liste[material].reshape(
            alphas_liste[material].shape[0])

    """input_filename = 'dta_alpha_' + str(material) + '.mtx'
    alphas = scipy.io.mmread(input_filename)
    alphas = alphas.reshape(alphas.shape[0])
    input_filename = 'dta_error_' + str(material) + '.mtx'"""

    errors_liste = {}
    for material in lmaterial:
        input_filename = 'dta_alpha_' + str(material) + namefun + '.mtx'
        errors_liste[material] = scipy.io.mmread(input_filename)
        errors_liste[material] = alphas_liste[material].reshape(
            alphas_liste[material].shape[0])
    print('Plotting alpha...')
    # print(" ICI ON VA VERIFIER")
    # print("c'est vrai? ", alphas_liste['beton']
    # == alphas_liste['Birch_LT'])
    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.subplot(1, 1, 1)
    for material in lmaterial:

        matplotlib.pyplot.plot(numpy.real(omega_liste[material]), numpy.real(
            alphas_liste[material]), color=color[material], label=str(material))
    matplotlib.pyplot.xlabel(r'$\omega$')
    # matplotlib.pyplot.ylim(0, 35)
    matplotlib.pyplot.ylabel(r'$\operatorname{Re}(\alpha)$')
    # matplotlib.pyplot.show()
    matplotlib.pyplot.savefig(
        'fig_alpha_real ' + namefun + '.jpg')
    matplotlib.pyplot.close(fig)
    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.subplot(1, 1, 1)
    for material in lmaterial:

        matplotlib.pyplot.plot(numpy.real(omega_liste[material]), numpy.imag(
            alphas_liste[material]), color=color[material], label=str(material))
        # matplotlib.pyplot.plot(numpy.real(omegas), numpy.imag(alphas), color=color)
    matplotlib.pyplot.xlabel(r'$\omega$')
    matplotlib.pyplot.ylabel(r'$\operatorname{Im}(\alpha)$')
    # matplotlib.pyplot.ylim(-120, 10)
    # matplotlib.pyplot.show()
    matplotlib.pyplot.savefig(
        'fig_alpha_imag_' + namefun + '.jpg')
    matplotlib.pyplot.close(fig)

    for material in lmaterial:
        fig = matplotlib.pyplot.figure()
        ax = matplotlib.pyplot.axes()

        ax.fill_between(numpy.real(omega_liste[material]), numpy.real(
            errors_liste[material]), color=color[material], label="fonction étudié:"+str(namefun) + "matériaux"+str(material))
        # matplotlib.pyplot.ylim(1.e-9, 1.e-4)
        matplotlib.pyplot.yscale('log')
        matplotlib.pyplot.xlabel(r'$\omega$')
        matplotlib.pyplot.ylabel(r'$e(\alpha)$')
        # matplotlib.pyplot.show()
        # ! NE DEPEND PAS DU MATERIAU
        matplotlib.pyplot.savefig(
            'fig_error_' + "matériaux_" + str(material) + "_fun_" + namefun + '.jpg')
        matplotlib.pyplot.close(fig)
    fig = matplotlib.pyplot.figure()  # ! Calcul de Im/Re
    matplotlib.pyplot.subplot(1, 1, 1)
    for material in lmaterial:

        matplotlib.pyplot.plot(numpy.real(omega_liste[material]), numpy.real(
            alphas_liste[material]), color=color[material],)
    matplotlib.pyplot.xlabel(r'$\omega$')
    matplotlib.pyplot.ylabel(r'$\operatorname{Re/Im}(\alpha)$')
    matplotlib.pyplot.ylim(0, 35)
    matplotlib.pyplot.legend(loc='best')
    matplotlib.pyplot.show()
    matplotlib.pyplot.savefig('fig_alpha_real_im' +
                              str(material) + namefun + '.jpg')
    matplotlib.pyplot.close(fig)
    return


def fonction_final(fun, namefun):
    lmaterial = ['Birch_LT', 'beton', 'nanofiber']
    #lmaterial = ['Birch_LT']
    for material in lmaterial:  # * Pour chaque matériaux on compute son alpha
        run_compute_alpha(material, fun)  # ! A ENLEVER QUAND FINI
    run_plot_alpha(lmaterial, fun, str(fun.__name__))


#! Ici Pierre=Sandstone
if __name__ == '__main__':
    """omega = 100
    compute_alpha(omega, 'Birch_LT', id)
    """
    fonction_final(id, "g_usuel")
