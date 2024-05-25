import quaternion_math.quaternionMath as qm
from astropy import units as u

import numpy as np


def solveQLaw(orbit, Wp, rp_min, f_mag,
              Wa, We, Wi, Wraan, Wargp,
              aT, eT, iT, raanT, argpT):
    
    dxx_dtheta, dxx_dr, dxx_dh = return_oe_partials(orbit)
    dQ_da, dQ_de, dQ_di, dQ_draan, dQ_dargp = return_Q_partials(orbit, Wp, rp_min, f_mag,
                                                                Wa, We, Wi, Wraan, Wargp,
                                                                aT, eT, iT, raanT, argpT)
    
    da_dtheta, de_dtheta, di_dtheta, draan_dtheta, dargp_dtheta = dxx_dtheta[0], dxx_dtheta[1], dxx_dtheta[2], dxx_dtheta[3], dxx_dtheta[4]
    da_dr, de_dr, di_dr, draan_dr, dargp_dr = dxx_dr[0], dxx_dr[1], dxx_dr[2], dxx_dr[3], dxx_dr[4]
    da_dh, de_dh, di_dh, draan_dh, dargp_dh = dxx_dh[0], dxx_dh[1], dxx_dh[2], dxx_dh[3], dxx_dh[4]
    
    D1 = dQ_da * da_dtheta + dQ_de * de_dtheta + dQ_di * di_dtheta + dQ_dargp * dargp_dtheta + dQ_draan * draan_dtheta
    D2 = dQ_da * da_dr + dQ_de * de_dr + dQ_di * di_dr + dQ_dargp * dargp_dr + dQ_draan * draan_dr
    D3 = dQ_da * da_dh + dQ_de * de_dh + dQ_di * di_dh + dQ_dargp * dargp_dh + dQ_draan * draan_dh
    
    # u = { u_r, u_theta, u_h }
    u = -np.array([D2, D1, D3]) / qm.magnitude([D1, D2, D3])
    
    
    a_dot_xx, e_dot_xx, i_dot_xx, raan_dot_xx, argp_dot_xx = return_oe_dot_xx(f_mag, orbit)
    q_dot = dQ_da*a_dot_xx + dQ_de*e_dot_xx + dQ_di*i_dot_xx + dQ_draan*raan_dot_xx + dQ_dargp*argp_dot_xx
    
    return u, q_dot


def return_Q_partials(orbit, Wp, rp_min, f_mag,
                      Wa, We, Wi, Wraan, Wargp,
                      aT, eT, iT, raanT, argpT):
    
    mu =    3.986 * 10**14
    
    a =     (orbit.a << u.meter).value
    e =     orbit.ecc.value
    i =     (orbit.inc << u.radian).value
    raan =  (orbit.raan << u.radian).value
    argp =  (orbit.argp << u.radian).value
    nu =    (orbit.nu << u.radian).value
    
    r = (orbit.r << u.meter).value
    v = (orbit.v << u.meter/u.second).value
    
    p =     a*(1-e**2)
    h =     qm.magnitude(np.cross(r, v))
    r_mag = qm.magnitude(r)
    
    
    Wp = 1
    
    term1 = (1 - e)/rp_min
    term2 = ( (a - aT)/(3*aT) )
    
    P = np.exp(100 * (1 - a*term1))
    S = np.sqrt( 1 + term2**4 )
    
    
    a_dot_xx, e_dot_xx, i_dot_xx, raan_dot_xx, argp_dot_xx = return_oe_dot_xx(f_mag, orbit)
    
    da = (( a - aT ) / a_dot_xx)
    
    ####################
    
    Va =            S * ( (a - aT)/a_dot_xx )**2
    dP_da =         100 * (1 - a*term1) * term1 * P
    
    da_dot_xx_da =  6*f_mag*( (a**2 * (1 + e)) / (mu * (1 - e)) )
    dVa_da =        ( (2/(3*aT)) * (term2**3) * (1 + term2**4)**(-1/2) ) * da**2 + 2*S*da*( (a_dot_xx - (a-aT)*da_dot_xx_da) / a_dot_xx**2 )
    
    dV_da =         dVa_da
    
    dQ_da =         Wp*dP_da*Wa*Va + (1+Wp*P)*Wa*dV_da
    
    ####################
    
    Ve =            ( (e - eT)/e_dot_xx )**2
    dP_de =         100 * (1 - a*term1) * (a / rp_min) * P
    
    de_dot_xx_de = -4*(f_mag/h)*a*e
    dVe_de =        2 * ( (e - eT)/e_dot_xx ) * ((e_dot_xx - (e-eT)*de_dot_xx_de) / (e_dot_xx**2))
    dV_de =         dVe_de
    
    dQ_de =         Wp*dP_de*We*Ve + (1+Wp*P)*We*dV_de

    ####################
    
    di =            np.arccos(np.cos(i - iT))
    Vi =            (di / i_dot_xx)**2
    dP_di =         0
    
    di_dot_xx_di =  0
    dVi_di =        2 * (di / i_dot_xx) * ((i_dot_xx - di*di_dot_xx_di) / i_dot_xx**2)
    dV_di =         dVi_di
    
    dQ_di =         Wp*dP_di*Wi*Vi + (1+Wp*P)*Wi*dV_di
  
    ####################
    
    draan =         np.arccos(np.cos(raan - raanT))
    Vraan =         (draan / raan_dot_xx)**2
    dP_draan =      0
    
    draan_dot_xx_draan = 0
    dVraan_draan =  2 * (draan / raan_dot_xx) * ((raan_dot_xx - draan*draan_dot_xx_draan) / raan_dot_xx**2)
    dV_draan =      dVraan_draan
    
    dQ_draan =      Wp*dP_draan*Wraan*Vraan + (1+Wp*P)*Wraan*dV_draan
    
    ####################
        
    dQ_dargp = 0
    
    ####################
    
    return dQ_da, dQ_de, dQ_di, dQ_draan, dQ_dargp


def return_oe_dot_xx(f_mag, orbit):

    mu = 3.986 * 10**14
    
    a =     (orbit.a << u.meter).value
    e =     orbit.ecc.value
    i =     (orbit.inc << u.radian).value
    raan =  (orbit.raan << u.radian).value
    argp =  (orbit.argp << u.radian).value
    nu =    (orbit.nu << u.radian).value
    
    r = (orbit.r << u.meter).value
    v = (orbit.v << u.meter/u.second).value
    
    p =     a*(1-e**2)
    h =     qm.magnitude(np.cross(r, v))
    r_mag = qm.magnitude(r)
    
    
    a_dot_xx =      2*f_mag*np.sqrt( (a**3 * (1+e)) / (mu * (1-e)) ) 
    e_dot_xx =      2*p*f_mag / h
    i_dot_xx =      p*f_mag / ( h * (np.sqrt(1 - e**2 * np.cos(argp)**2) - e*abs(np.cos(argp))) )
    raan_dot_xx =   p*f_mag / ( h * np.sin(i) * (np.sqrt(1 - e**2 * np.cos(argp)**2) - e*abs(np.sin(argp))) )
    
    b = 0.01
    cos_theta_xx =  ( ((1-e**2)/(2*e**3)) + np.sqrt(0.25*((1-e**2)/(e**3))**2 + (1/27)) )**(1/3) - ( -((1-e**2)/(2*e**3)) + np.sqrt(0.25*((1-e**2)/(e**3))**2 + (1/27)) )**(1/3) - (1/e)
    rxx =           p / (1 + e*cos_theta_xx)
    argp_dot_xx_i = (f_mag / (e*h)) * np.sqrt(abs(p**2 * cos_theta_xx + (p+rxx)**2 * np.sin(nu))) # is cos_theta_xx == cos(nu) ?
    argp_dot_xx_o = raan_dot_xx * abs(np.cos(i))
    
    argp_dot_xx =   (argp_dot_xx_i + b*argp_dot_xx_o) / (1 + b)
    
    
    return a_dot_xx, e_dot_xx, i_dot_xx, raan_dot_xx, argp_dot_xx


def return_oe_partials(orbit):
    
    a =     (orbit.a << u.meter).value
    e =     orbit.ecc.value
    i =     (orbit.inc << u.radian).value
    raan =  (orbit.raan << u.radian).value
    argp =  (orbit.argp << u.radian).value
    nu =    (orbit.nu << u.radian).value
    
    r = (orbit.r << u.meter).value
    v = (orbit.v << u.meter/u.second).value
    
    p =     a*(1-e**2)
    h =     qm.magnitude(np.cross(r, v))
    r_mag = qm.magnitude(r)
    
    
    da_dtheta =     2 * ((a**2) / h) * (p / r_mag)
    da_dr =         2 * (a**2 / h) * e * np.sin(nu)
    da_dh =         0
    
    de_dtheta =     (1/h) * ((p+r_mag)*np.cos(nu) + r_mag*e)
    de_dr =         (p / h) * np.sin(nu)
    de_dh =         0
    
    di_dtheta =     0
    di_dr =         0
    di_dh =         (1/h) * r_mag * np.cos(nu + argp)
    
    draan_dtheta =  0
    draan_dr =      0
    draan_dh =      (r_mag * np.sin(nu + argp)) / (h * np.sin(i))
    
    dargp_dtheta =  -p*np.cos(nu) / (e*h)
    dargp_dr =      ((p+r_mag)/(e*h)) * np.sin(nu)
    dargp_dh =      (-r_mag*np.sin(nu+argp)*np.cos(i))/(h*np.sin(i))
    
    
    dxx_dtheta =    [da_dtheta, de_dtheta, di_dtheta, draan_dtheta, dargp_dtheta]
    dxx_dr =        [da_dr, de_dr, di_dr, draan_dr, dargp_dr]
    dxx_dh =        [da_dh, de_dh, di_dh, draan_dh, dargp_dh]
    
    return dxx_dtheta, dxx_dr, dxx_dh