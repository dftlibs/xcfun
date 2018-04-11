import sys
from matplotlib import pyplot as plt
from numpy import *
import subprocess

evaluator_args = ['../funeval', '--quiet', '--order', '0', '--ngnntaun', 'm06l', '1.0']


def gen_grad2_tau(rho, reduced_gradient_s, iso_orbital_indicator_z):
    """Synthesize square gradients and tau from rho and fixed values of s and z"""
    out = zeros((len(rho), 3))
    out[:, 0] = rho
    kf = (3 * pi**2 * rho)**(1.0 / 3.0)
    # gnn = s^2*(2*kf*n)
    out[:, 1] = reduced_gradient_s**2 * (2 * kf * rho)**2  # |gn|^2
    tauw = 1.0 / 8.0 * out[:, 1] / rho  # Tau weizecker
    out[:, 2] = tauw / iso_orbital_indicator_z
    return out


def eval_fun(rho, s, z):
    try:
        evaluator = subprocess.Popen(evaluator_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    except OSError:
        print "Trouble running '%s'. This script must be run from the test/ directory and funeval must be present in ../." % evaluator_args[
            0]
        sys.exit(-1)
    inpd = gen_grad2_tau(rho, s, z)
    inp = '\n'.join(["%.15e %.15e %.15e\n" % tuple(inpd[i, :]) for i in range(inpd.shape[0])])
    out = evaluator.communicate(inp)
    return array(map(float, out[0].split()))


plt.suptitle(' '.join(evaluator_args[5:]))
rho = linspace(1e-6, 2e1, 100)
ss = [1e-6, 0.1, 0.5, 1.0]
zs = [1e-3, 0.5, 1.0]
for i in range(len(zs)):
    plt.subplot(len(zs), 1, len(zs) - i)
    for s in ss:
        plt.plot(rho, eval_fun(rho, s, zs[i]), '-', label='s = %.1e' % s)
    plt.ylabel('F[rho]')
    plt.title('z = %.1e' % zs[i])
    if i == 0:
        plt.xlabel('rho')
        plt.legend()
plt.show()
