import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


P = 1.01325 * 10**5 # 1 atm [N m**-2]

R = 8.3145 # [J / (K mol) = N m / (K mol) = kg (m/s)**2 / (K mol)]

N_avo = 6.022 * 10**23 # [1/mol]

k_B = 1.38 * 10**-23 # [J/K = N m / K = kg m/s**2 m/K = kg (m/s)**2 / K]

h = 6.626*10**-34   # [m**2 kg/s]


def mass(M_r): # molecular weight to mass per particle

    # 분자량(M_r)을 넣어주면, 이 분자 하나 당 질량을 계산해 줍니다.

    return M_r / N_avo * 10**-3 # 단위: [kg]

def Λ(m,T):

    # Thermal de Broglie wavelength를 계산해줍니다.

    # np.sqrt는 square root를 의미합니다.

    return h / np.sqrt(2*np.pi*m*k_B*T) # 단위: [m]



def q_rot_rod(Θ_rot, σ, T):

    # rod 모양의 분자의 rotational partition function을 계산해줍니다.

    # rotational partition function for rod shape molecule: N2, H2
    return 1/σ * T / Θ_rot

def q_rot_top(Θ_rot_1,Θ_rot_2,Θ_rot_3,σ, T):

    # top 모양의 분자의 rotational partition function을 계산해줍니다.

    # rotational partition function for top shape molecule: NH3
    return np.sqrt(np.pi) / σ * (T/ Θ_rot_1)**(1/2) * (T/ Θ_rot_2)**(1/2) * (T/ Θ_rot_3)**(1/2) 

def q_vib(Θ_vib,T):

    # 각 vibrational mode에 대한 partition function을 계산해줍니다.
    # 여기에서 조심해야될 점은 Ammonia의 vibrational mode는 degeneracy를 가졌다는 점입니다. 이 경우 degeneracy 횟수만큼 이를 반복해주셔야 됩니다.

    return 1 / (1 - np.exp(-Θ_vib/T))

def q_elec(D_O, T):

    # 분자의 electrical partition function을 계산해줍니다.

    # D_O : [J / mol]
    # R : [J / (K mol)]
    return np.exp(D_O * (R * T)**-1) # D_O = D_e + \sum_i exp(\theta_vib_i / 2T)

def log10_q_elec(D_O, T):

    # 분자의 electrical partition function을 계산해줍니다.

    # D_O : [J / mol]
    # R : [J / (K mol)]
    return D_O * (R * T)**-1 * np.log10(np.exp(1)) # D_O = D_e + \sum_i exp(\theta_vib_i / 2T)


def K_f(T):

    # part 1.

    m_A = 28.014 # Nitrogen (N2)
    m_B = 2.01588 # Hydrogen (H2)
    m_C = 17.031 # Ammonia (NH3)

    # part 1. result

    value_trans = (Λ(mass(m_A),T)**1.5 * Λ(mass(m_B),T)**4.5 * Λ(mass(m_C),T)**-3) * (P / (k_B * T))

    # part 2.

    # 각 분자의 rotational temperatures, vibrational temperatures, D_O 등의 자료는 중간 레포트때 받으신 McQuarrie 책 Table 18.2, 18.4, 18.5 에 나와있습니다. 

    # internal part 1. Nitrogen

    Θ_rot_A = 2.88
    Θ_vib_A = 3374
    D_O_A = 941.6 * 10**3 # J / mol
    g_e_A = 1
    σ_A = 2

    q_int_A = q_rot_rod(Θ_rot_A, σ_A, T) *  q_vib(Θ_vib_A, T) * g_e_A * q_elec(D_O_A, T)

    # internal part 2. Hydrogen

    Θ_rot_B = 85.3
    Θ_vib_B = 6332
    D_O_B = 432.1 * 10**3 # J / mol
    g_e_B = 1
    σ_B = 2

    q_int_B = q_rot_rod(Θ_rot_B, σ_B, T) *  q_vib(Θ_vib_B, T) * g_e_B * q_elec(D_O_B, T)

    # internal part 3. Ammonia

    Θ_rot_C_1 = 13.6
    Θ_rot_C_2 = 13.6
    Θ_rot_C_3 = 8.92

    Θ_vib_C_1 = 4800
    Θ_vib_C_2 = 1360
    Θ_vib_C_3 = 4880 # with degeneracy 2
    Θ_vib_C_4 = 2330 # with degeneracy 2

    D_O_C = 1158 * 10**3 # J / mol
    σ_C = 3

    q_int_C = q_rot_top(Θ_rot_C_1, Θ_rot_C_2, Θ_rot_C_3, σ_C, T) 
    q_int_C *= q_vib(Θ_vib_C_1, T) * q_vib(Θ_vib_C_2, T) * q_vib(Θ_vib_C_3, T)**2 * q_vib(Θ_vib_C_4,T)**2 
    q_int_C *= q_elec(D_O_C, T)


    # part 2. result

    value_internal = q_int_A**-0.5 * q_int_B**-1.5 * q_int_C 

    # final value

    return value_trans * value_internal

def log10_K_f(T):

    # part 1.

    m_A = 28.014 # Nitrogen (N2)
    m_B = 2.01588 # Hydrogen (H2)
    m_C = 17.031 # Ammonia (NH3)

    # part 1. result

    value_trans = (Λ(mass(m_A),T)**1.5 * Λ(mass(m_B),T)**4.5 * Λ(mass(m_C),T)**-3) * (P / (k_B * T))

    log_value_trans = np.log10(value_trans)

    # part 2.

    # 각 분자의 rotational temperatures, vibrational temperatures, D_O 등의 자료는 중간 레포트때 받으신 McQuarrie 책 Table 18.2, 18.4, 18.5 에 나와있습니다. 

    # internal part 1. Nitrogen

    Θ_rot_A = 2.88
    Θ_vib_A = 3374
    D_O_A = 941.6 * 10**3 # J / mol
    g_e_A = 1
    σ_A = 2

    # q_int_A = q_rot_rod(Θ_rot_A, σ_A, T) *  q_vib(Θ_vib_A, T) * g_e_A * q_elec(D_O_A, T)

    log_q_int_A = np.log10(q_rot_rod(Θ_rot_A, σ_A, T) *  q_vib(Θ_vib_A, T) * g_e_A) + log10_q_elec(D_O_A, T)


    # internal part 2. Hydrogen

    Θ_rot_B = 85.3
    Θ_vib_B = 6332
    D_O_B = 432.1 * 10**3 # J / mol
    g_e_B = 1
    σ_B = 2

    # q_int_B = q_rot_rod(Θ_rot_B, σ_B, T) *  q_vib(Θ_vib_B, T) * g_e_B * q_elec(D_O_B, T)

    log_q_int_B = np.log10(q_rot_rod(Θ_rot_B, σ_B, T) *  q_vib(Θ_vib_B, T) * g_e_B) + log10_q_elec(D_O_B, T)

    # internal part 3. Ammonia

    Θ_rot_C_1 = 13.6
    Θ_rot_C_2 = 13.6
    Θ_rot_C_3 = 8.92

    Θ_vib_C_1 = 4800
    Θ_vib_C_2 = 1360
    Θ_vib_C_3 = 4880 # with degeneracy 2
    Θ_vib_C_4 = 2330 # with degeneracy 2

    D_O_C = 1158 * 10**3 # J / mol
    σ_C = 3

    # q_int_C = q_rot_top(Θ_rot_C_1, Θ_rot_C_2, Θ_rot_C_3, σ_C, T) 
    # q_int_C *= q_vib(Θ_vib_C_1, T) * q_vib(Θ_vib_C_2, T) * q_vib(Θ_vib_C_3, T)**2 * q_vib(Θ_vib_C_4,T)**2 
    # q_int_C *= q_elec(D_O_C, T)

    log_q_int_C = np.log10(q_rot_top(Θ_rot_C_1, Θ_rot_C_2, Θ_rot_C_3, σ_C, T)) 
    log_q_int_C += np.log10(q_vib(Θ_vib_C_1, T) * q_vib(Θ_vib_C_2, T) * q_vib(Θ_vib_C_3, T)**2 * q_vib(Θ_vib_C_4,T)**2) 
    log_q_int_C += log10_q_elec(D_O_C, T)

    # part 2. result

    # value_internal = q_int_A**-0.5 * q_int_B**-1.5 * q_int_C 

    log_value_internal = -0.5 * log_q_int_A - 1.5 * log_q_int_B + log_q_int_C 


    # final value

    return log_value_trans + log_value_internal