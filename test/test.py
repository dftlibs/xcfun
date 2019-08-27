#!/usr/bin/env python

import pytest

import numpy
from numpy.testing import assert_allclose

import xcfun


def test_splash():

    version = xcfun.xcfun_version()
    splash = xcfun.xcfun_splash()

    print("\n Sample use of python interface to xcfun")

    print("\n\nXCFun version ", version)
    print(splash)


def test_lda_old():
    fun_xc_name = 'LDA'
    fun_xc_weight = 1.0
    fun_xc = xcfun.Functional({fun_xc_name: fun_xc_weight})

    n = numpy.array([1.0])
    result = fun_xc.eval_potential_n(n)

    energy = result[0, 0]
    potential = result[0, 1]

    print("Functional used: ", fun_xc_name, " weight ", fun_xc_weight)
    print("        density: ", n)
    print("         energy: ", energy)
    print("      potential: ", potential)

    assert_allclose(energy, -0.8101513)
    assert_allclose(potential, -1.06468341)


@pytest.fixture
def lda_fun():
    return xcfun.Functional({'LDA': 1.0})

@pytest.fixture
def pbe_fun():
    return xcfun.Functional({'PBE': 1.0})

def test_type(lda_fun, pbe_fun):
    assert lda_fun.type == 0
    assert pbe_fun.type == 1

@pytest.fixture
def dens():
    return numpy.arange(0.1, 5.0, 0.5)

@pytest.fixture
def spindens(dens):
    return 0.1*numpy.sin(dens)

@pytest.fixture
def densgrad(dens):    
    dg = numpy.zeros((dens.shape[0], 3))
    dg[:,0] = 0.2*numpy.cos(dens)
    dg[:,1] = 0.3*numpy.sin(dens)
    dg[:,2] = 0.1*numpy.sin(dens)**2
    return dg

@pytest.fixture
def denshess(dens):
    dh = numpy.zeros((dens.shape[0], 6))
    dh[:,0] = 0.5
    dh[:,1] = 0.1 * dens
    dh[:,2] = numpy.cos(dens)**2
    dh[:,3] = numpy.sin(dens)
    dh[:,4] = 0.6*numpy.cos(dens)**2
    dh[:,5] = 0.3*numpy.sin(dens)
    return dh 

@pytest.fixture
def refen_lda():
    return numpy.array([-0.03962059, -0.41418026, -0.91826807, -1.50300217, -2.14968647, 
                        -2.84774973, -3.59023064, -4.37213299, -5.18966423, -6.03982917])
@pytest.fixture
def refpot_lda():
    return numpy.array([-0.51789018, -0.90610372, -1.09730105, -1.23582171, -1.34751238, 
                        -1.44246858, -1.52581304, -1.60054497, -1.66858921, -1.73126329])

def test_lda_energy_restr(lda_fun, dens, refen_lda):
    out = lda_fun.eval_energy_n(dens)
    assert_allclose(out, refen_lda)

def test_lda_potential_restr(lda_fun, dens, refen_lda, refpot_lda):
    out = lda_fun.eval_potential_n(dens)
    assert_allclose(out[:,0], refen_lda)
    assert_allclose(out[:,1], refpot_lda)

def test_lda_energy_unrestr_closedshell(lda_fun, dens, refen_lda):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:]
    rho[:,1] = 0.5 * dens[:]

    out = lda_fun.eval_energy_ab(rho)
    assert_allclose(out, refen_lda)
 
def test_lda_potential_unrestr_closedshell(lda_fun, dens, refen_lda, refpot_lda):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:]
    rho[:,1] = 0.5 * dens[:]

    out = lda_fun.eval_potential_ab(rho)
    assert_allclose(out[:,0], refen_lda)
    assert_allclose(out[:,1], refpot_lda)
    assert_allclose(out[:,2], refpot_lda)

def test_lda_energy_unrestr_openshell(lda_fun, dens, spindens):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:] + spindens[:]
    rho[:,1] = 0.5 * dens[:] - spindens[:]

    out = lda_fun.eval_energy_ab(rho)

    refen_lda_openshell = numpy.array([-0.03985427, -0.41666061, -0.92248897, -1.50718908, -2.15231233,
                                       -2.84856692, -3.59023539, -4.3726223,  -5.19120401, -6.04193876])

    assert_allclose(out, refen_lda_openshell)

def test_lda_potential_unrestr_openshell(lda_fun, dens, spindens):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:] + spindens[:]
    rho[:,1] = 0.5 * dens[:] - spindens[:]

    out = lda_fun.eval_potential_ab(rho)

    refpot_lda_alpha = numpy.array([-0.53994997, -0.94756428, -1.14234107, -1.27610516, -1.37715554,
                                     -1.4581231, -1.52695501, -1.58940021, -1.64952872, -1.70973604])
    refpot_lda_beta  = numpy.array([-0.49296826, -0.85942489, -1.04738939, -1.1922106,  -1.31627816,
                                    -1.42641378, -1.52466911, -1.61151607, -1.68716914, -1.75220302])

    assert_allclose(out[:,1], refpot_lda_alpha)
    assert_allclose(out[:,2], refpot_lda_beta)

    outen = lda_fun.eval_energy_ab(rho)
    assert_allclose(out[:,0], outen)


@pytest.fixture
def refen_pbe():
    return numpy.array([-0.04051725, -0.41397904, -0.91782841, -1.50231227, -2.14873903,
                        -2.84654021, -3.58875595, -4.37039097, -5.18765335, -6.03754832])

@pytest.fixture
def refpot_pbe():
    return numpy.array([-0.48524858, -0.90474127, -1.09660442, -1.23524903, -1.34697602,
                        -1.44193844, -1.52527933, -1.60000832, -1.66805094, -1.73072351])

def test_pbe_energy_restr(pbe_fun, dens, densgrad, refen_pbe):
    out = pbe_fun.eval_energy_n(dens, densgrad)
    assert_allclose(out, refen_pbe)

def test_pbe_potential_restr(pbe_fun, dens, densgrad, denshess, refen_pbe, refpot_pbe):
    out = pbe_fun.eval_potential_n(dens, densgrad, denshess)
    assert_allclose(out[:,0], refen_pbe)
    assert_allclose(out[:,1], refpot_pbe)

def test_pbe_energy_unrestr_closedshell(pbe_fun, dens, densgrad, refen_pbe):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:]
    rho[:,1] = 0.5 * dens[:]

    rhograd = numpy.zeros((dens.shape[0], 3, 2))
    rhograd[:,0:3,0] = 0.5 * densgrad[:,0:3]
    rhograd[:,0:3,1] = 0.5 * densgrad[:,0:3]

    out = pbe_fun.eval_energy_ab(rho, rhograd)
    assert_allclose(out, refen_pbe)
 
def test_pbe_potential_unrestr_closedshell(pbe_fun, dens, densgrad, denshess, refen_pbe, refpot_pbe):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = 0.5 * dens[:]
    rho[:,1] = 0.5 * dens[:]

    rhograd = numpy.zeros((dens.shape[0], 3, 2))
    rhograd[:,0:3,0] = 0.5 * densgrad[:,0:3]
    rhograd[:,0:3,1] = 0.5 * densgrad[:,0:3]
    
    rhohess = numpy.zeros((dens.shape[0], 6, 2))
    rhohess[:,0:6,0] = 0.5 * denshess[:,0:6]
    rhohess[:,0:6,1] = 0.5 * denshess[:,0:6]

    out = pbe_fun.eval_potential_ab(rho, rhograd, rhohess)
    assert_allclose(out[:,0], refen_pbe)
    assert_allclose(out[:,1], refpot_pbe)
    assert_allclose(out[:,2], refpot_pbe)

def test_pbe_energy_unrestr_openshell(pbe_fun, dens, spindens, densgrad, refen_pbe):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = dens[:] + spindens[:]
    rho[:,1] = dens[:] - spindens[:]

    rhograd = numpy.zeros((dens.shape[0], 3, 2))
    rhograd[:,0:3,0] = 0.4 * densgrad[:,0:3]
    rhograd[:,0:3,1] = 0.6 * densgrad[:,0:3]

    out = pbe_fun.eval_energy_ab(rho, rhograd)

    refen_pbe_openshell = numpy.array([-0.09848243, -1.03067727, -2.28717791, -3.74476386, -5.35679716,
                                       -7.09762348, -8.95059771,-10.90338775,-12.94601448,-15.07029572])

    assert_allclose(out, refen_pbe_openshell)
 
def test_pbe_potential_unrestr_closedshell(pbe_fun, dens, spindens, densgrad, denshess, refen_pbe, refpot_pbe):
    rho = numpy.zeros((dens.shape[0],2))
    rho[:,0] = dens[:] + spindens[:]
    rho[:,1] = dens[:] - spindens[:]

    rhograd = numpy.zeros((dens.shape[0], 3, 2))
    rhograd[:,0:3,0] = 0.4 * densgrad[:,0:3]
    rhograd[:,0:3,1] = 0.6 * densgrad[:,0:3]

    rhohess = numpy.zeros((dens.shape[0], 6, 2))
    rhohess[:,0:6,0] = 0.7 * denshess[:,0:6]
    rhohess[:,0:6,1] = 0.2 * denshess[:,0:6]

    out = pbe_fun.eval_potential_ab(rho, rhograd, rhohess)

    print out[:,1]
    print out[:,2]

    refpot_pbe_alpha = numpy.array([-0.62768401, -1.15201996, -1.39499424, -1.56595725, -1.69939803,
                                    -1.80975816, -1.90530719, -1.99151265, -2.07206904, -2.1492557 ])
    refpot_pbe_beta  = numpy.array([-0.64069084, -1.10200908, -1.33808862, -1.51488446, -1.66214897,
                                    -1.79045119, -1.90425337, -2.00565175, -2.09589401, -2.17609771])
    assert_allclose(out[:,1], refpot_pbe_alpha)
    assert_allclose(out[:,2], refpot_pbe_beta)

    outen = pbe_fun.eval_energy_ab(rho, rhograd)
    assert_allclose(out[:,0], outen)

def test_raises_eval_energy_n(lda_fun, pbe_fun):
    # correct shape for all arguments (LDA)
    lda_fun.eval_energy_n(numpy.zeros((10,)))

    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_energy_n(numpy.zeros((10,2)))

    # correct shape for all arguments (PBE)
    pbe_fun.eval_energy_n(numpy.zeros((10,)), numpy.zeros((10,3)))

    with pytest.raises(xcfun.XCFunException, match="Density gradient required"):
        pbe_fun.eval_energy_n(numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_n(numpy.zeros((10,)), numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_n(numpy.zeros((10,)), numpy.zeros((10,4)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_n(numpy.zeros((10,)), numpy.zeros((20,3)))


def test_raises_eval_potential_n(lda_fun, pbe_fun):
    # correct shape for all arguments (LDA)
    lda_fun.eval_potential_n(numpy.zeros((10,)))

    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_potential_n(numpy.zeros((10,2)))

    # correct shape for all arguments (PBE)
    pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,3)), numpy.zeros((10,6)))

    with pytest.raises(xcfun.XCFunException, match="Density gradient and Hessian required"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Density gradient and Hessian required"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,3)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,)), numpy.zeros((10,6)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,4)), numpy.zeros((10,6)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((20,3)), numpy.zeros((10,6)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,3)), numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,3)), numpy.zeros((10,3)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_n(numpy.zeros((10,)), numpy.zeros((10,3)), numpy.zeros((20,6)))

def test_raises_eval_energy_ab(lda_fun, pbe_fun):
    # correct shape for all arguments (LDA)
    lda_fun.eval_energy_ab(numpy.zeros((10,2)))

    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_energy_ab(numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_energy_ab(numpy.zeros((10,4)))

    # correct shape for all arguments (PBE)
    pbe_fun.eval_energy_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)))

    with pytest.raises(xcfun.XCFunException, match="Density gradient required"):
        pbe_fun.eval_energy_ab(numpy.zeros((10,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_ab(numpy.zeros((10,2)), numpy.zeros((10,3)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_ab(numpy.zeros((10,2)), numpy.zeros((10,4,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_energy_ab(numpy.zeros((10,2)), numpy.zeros((20,3,2)))

def test_raises_eval_potential_ab(lda_fun, pbe_fun):
    # correct shape for all arguments (LDA)
    lda_fun.eval_potential_ab(numpy.zeros((10,2)))

    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_potential_ab(numpy.zeros((10,)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of density argument"):
        lda_fun.eval_potential_ab(numpy.zeros((10,4)))

    # correct shape for all arguments (PBE)
    pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)), numpy.zeros((10,6,2)))

    with pytest.raises(xcfun.XCFunException, match="Density gradient and Hessian required"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)))
    with pytest.raises(xcfun.XCFunException, match="Density gradient and Hessian required"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3)), numpy.zeros((10,6,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,4,2)), numpy.zeros((10,6,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of densgrad argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((20,3,2)), numpy.zeros((10,6,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)), numpy.zeros((10,6)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)), numpy.zeros((10,3,2)))
    with pytest.raises(xcfun.XCFunException, match="Wrong shape of denshess argument"):
        pbe_fun.eval_potential_ab(numpy.zeros((10,2)), numpy.zeros((10,3,2)), numpy.zeros((20,6,2)))

