##UNITS
#O2 in umol/Kg
#Pressure in db
#Temp in Celsius

def calc_pO2(data):
    import scipy
    import numpy as np
    import pandas
  
    a_0 = 5.80871
    a_1 = 3.20291
    a_2 = 4.17887
    a_3 = 5.10006
    a_4 = -9.86643e-2
    a_5 = 3.80369
    b_0 = -7.01577e-3
    b_1 = -7.70028e-3
    b_2 =  -1.13864e-2
    b_3 = -9.51519e-3
    c_0 = -2.75915E-7

    data['tt'] = 298.15 - data['temp']
    data['tk'] = 273.15 + data['temp']
    data['ts'] = np.log(data['tt'] / data['tk'])

    #correct for pressure at depth
    V = 32e-6 #partial molar volume of O2 (m3/mol)
    R = 8.31 #Gas constant [J/mol/K]
    db2Pa = 1e4 #convert pressure: decibar to Pascal
    atm2Pa = 1.01325e5 #convert pressure: atm to Pascal

    #calculate pressure in dB from depth in m
    #Let zdepth = z[g=temp_rcp[d=2]]
    #let zgrid = temp_rcp[d=2]*0+zdepth
    #Let ztop = zgrid * (1.0076 + zgrid * (2.3487e-06 - zgrid * 1.2887e-11)) 

    #convert pressure from decibar to pascal
    data['dp'] = data['pres']*db2Pa
    data['pCor'] = np.exp((V*data['dp'])/(R*(data['temp']+273.15)))

    #Let o2_alpha = (o2_saturation / 0.21)  #0.21 is atm composition of O2
    #Let kh = o2_alpha*pCor
    #Let po2_raw = (o2_rcp[d=1]/kh)*101.32501   #convert po2 from atm to kPa  

    data['o2_sat'] = np.exp(a_0 + a_1*data['ts'] + a_2*data['ts']**2 + a_3*data['ts']**3 + a_4*data['ts']**4 + a_5*data['ts']**5 + data['sal']*(b_0 + b_1*data['ts'] + b_2*data['ts']**2 + b_3*data['ts']**3) + c_0*data['sal']**2)
    
    data['o2_alpha'] = (data['o2_sat'] / 0.21)  #0.21 is atmospheric composition of O2
    data['kh'] = data['o2_alpha']*data['pCor']
    data['po2'] = (data['oxygen'] / data['kh'])*101.32501  #convert po2 from atm to kPa
    return data[['oxygen','pres','sal','temp','po2']]
