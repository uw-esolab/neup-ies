# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 09:27:40 2021

@author: b9801
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from iapws import IAPWS97

#Hamilton data 
raw_data = """TF Temp	HTF mdot	Ambient T	Wdot cycle	Heat In
500	0.300	30	0.237641861	0.31466289
504.2	0.300	30	0.242543393	0.318919287
508.4	0.300	30	0.247479625	0.323186151
512.6	0.300	30	0.252450559	0.327456504
516.8	0.300	30	0.257447519	0.331737323
521.1	0.300	30	0.262470504	0.336021631
525.3	0.300	30	0.267528191	0.340312917
529.5	0.300	30	0.272611903	0.34461118
533.7	0.300	30	0.277721641	0.348916421
537.9	0.300	30	0.28286608	0.35322864
542.1	0.300	30	0.288036545	0.357544348
546.3	0.300	30	0.293241711	0.361867033
550.5	0.300	30	0.298464228	0.366196696
554.7	0.300	30	0.303712771	0.370533337
558.9	0.300	30	0.308996014	0.374873467
563.2	0.300	30	0.314305284	0.379220574
567.4	0.300	30	0.319631904	0.38357117
571.6	0.300	30	0.3250019	0.387932233
575.8	0.300	30	0.330389247	0.392293296
580	0.300	30	0.335811295	0.396664825
500	1.000	30	0.766426017	0.810896011
504.2	1.000	30	0.782518656	0.824028044
508.4	1.000	30	0.798672022	0.837187987
512.6	1.000	30	0.81489479	0.850382818
516.8	1.000	30	0.83118696	0.863602072
521.1	1.000	30	0.847548532	0.876849236
525.3	1.000	30	0.863962156	0.890120822
529.5	1.000	30	0.880427832	0.903420319
533.7	1.000	30	0.896954234	0.916744238
537.9	1.000	30	0.913524013	0.930092579
542.1	1.000	30	0.930137169	0.943461853
546.3	1.000	30	0.946802376	0.956855549
550.5	1.000	30	0.963493608	0.970266689
554.7	1.000	30	0.980228218	0.98370574
558.9	1.000	30	0.996980177	0.997158747
563.2	1.000	30	1.013766838	1.010632686
567.4	1.000	30	1.030562174	1.02412407
571.6	1.000	30	1.047383536	1.037636387
575.8	1.000	30	1.064213574	1.051162659
580	1.000	30	1.081052286	1.064706376
500	1.200	30	0.889051057	0.924098734
504.2	1.200	30	0.907668201	0.939515964
508.4	1.200	30	0.926337397	0.954957615
512.6	1.200	30	0.94504997	0.970434154
516.8	1.200	30	0.963788568	0.985938605
521.1	1.200	30	0.982561867	1.001477943
525.3	1.200	30	1.001361192	1.017034727
529.5	1.200	30	1.020169193	1.032622909
533.7	1.200	30	1.038985868	1.048232025
537.9	1.200	30	1.057802544	1.063865563
542.1	1.200	30	1.076593194	1.079513056
546.3	1.200	30	1.095375169	1.095184971
550.5	1.200	30	1.114122442	1.110867353
554.7	1.200	30	1.132835014	1.126570668
558.9	1.200	30	1.15148686	1.14228096
563.2	1.200	30	1.170077978	1.158008697
567.4	1.200	30	1.188591019	1.173739923
571.6	1.200	30	1.207025983	1.189481615
575.8	1.200	30	1.225348168	1.205223307
580	1.200	30	1.243566249	1.220971977
560	0.300	0	0.338665461	0.375968966
560	0.347	0	0.396304002	0.425821147
560	0.395	0	0.452910184	0.474173122
560	0.442	0	0.508414606	0.521171423
560	0.490	0	0.562947396	0.567074225
560	0.537	0	0.616161545	0.611700108
560	0.584	0	0.668221882	0.655261892
560	0.632	0	0.719102381	0.697829354
560	0.679	0	0.768976549	0.739580424
560	0.726	0	0.817532074	0.780323217
560	0.774	0	0.864916436	0.820221708
560	0.821	0	0.910426939	0.859314275
560	0.868	0	0.9546188	0.89774396
560	0.916	0	0.997622148	0.935311899
560	0.963	0	1.039532411	0.97214718
560	1.010	0	1.080297537	1.008274224
560	1.058	0	1.119978254	1.043808163
560	1.105	0	1.158279602	1.078543156
560	1.153	0	1.195236281	1.112590845
560	1.200	0	1.230761541	1.14595472
560	0.300	30	0.310323332	0.375958499
560	0.347	30	0.363129742	0.425817658
560	0.395	30	0.415155378	0.474176611
560	0.442	30	0.466356863	0.521185378
560	0.490	30	0.51685565	0.567095158
560	0.537	30	0.566348105	0.611728019
560	0.584	30	0.614947006	0.655303758
560	0.632	30	0.66266103	0.697885175
560	0.679	30	0.709602953	0.73965369
560	0.726	30	0.755469142	0.780410438
560	0.774	30	0.800398402	0.820326374
560	0.821	30	0.844347355	0.859436385
560	0.868	30	0.887420105	0.897883514
560	0.916	30	0.929330368	0.935465408
560	0.963	30	0.970156221	0.972318133
560	1.010	30	1.009854288	1.008462622
560	1.058	30	1.048476621	1.044010516
560	1.105	30	1.085728261	1.078762953
560	1.153	30	1.121643907	1.112831576
560	1.200	30	1.156145484	1.146212895
560	0.300	55	0.233850765	0.375951522
560	0.347	55	0.276411674	0.42581068
560	0.395	55	0.318894505	0.474169633
560	0.442	55	0.361195156	0.52118189
560	0.490	55	0.403330976	0.567095158
560	0.537	55	0.444989656	0.611734997
560	0.584	55	0.486205897	0.655317714
560	0.632	55	0.526944997	0.697909597
560	0.679	55	0.567250334	0.739688579
560	0.726	55	0.606826947	0.78046626
560	0.774	55	0.64576159	0.820403128
560	0.821	55	0.683993536	0.859534072
560	0.868	55	0.721748341	0.898009113
560	0.916	55	0.758652969	0.935618918
560	0.963	55	0.794672719	0.972499554
560	1.010	55	0.829738188	1.008671953
560	1.058	55	0.863858053	1.044254736
560	1.105	55	0.896763378	1.079035084
560	1.153	55	0.928471515	1.113135106
560	1.200	55	0.958869686	1.146551313
500	1.000	0	0.830649092	0.810742502
500	1.000	2.895	0.830649092	0.810742502
500	1.000	5.789	0.830649092	0.810742502
500	1.000	8.684	0.830649092	0.810742502
500	1.000	11.58	0.830067849	0.810745991
500	1.000	14.47	0.827508642	0.810759946
500	1.000	17.37	0.822468306	0.810780879
500	1.000	20.26	0.814435	0.810805301
500	1.000	23.16	0.803200516	0.810833212
500	1.000	26.05	0.789146568	0.810857634
500	1.000	28.95	0.77277632	0.810885545
500	1.000	31.84	0.7549226	0.810913456
500	1.000	34.74	0.736123275	0.810934389
500	1.000	37.63	0.717089717	0.810955322
500	1.000	40.53	0.698108211	0.810976255
500	1.000	43.42	0.67959517	0.81099021
500	1.000	46.32	0.661628672	0.811000677
500	1.000	49.21	0.644434273	0.811004166
500	1.000	52.11	0.627942571	0.811004166
500	1.000	55	0.612283696	0.810993699
560	1.000	0	1.071379353	1.000344067
560	1.000	2.895	1.071379353	1.000344067
560	1.000	5.789	1.071379353	1.000344067
560	1.000	8.684	1.071379353	1.000344067
560	1.000	11.58	1.070589902	1.000347556
560	1.000	14.47	1.067874541	1.000365
560	1.000	17.37	1.062539245	1.000385933
560	1.000	20.26	1.053968072	1.000413844
560	1.000	23.16	1.041831359	1.000445244
560	1.000	26.05	1.026432743	1.000480132
560	1.000	28.95	1.008275388	1.000515021
560	1.000	31.84	0.988226823	1.00054642
560	1.000	34.74	0.966902992	1.000581309
560	1.000	37.63	0.945084671	1.000612709
560	1.000	40.53	0.92313622	1.000640619
560	1.000	43.42	0.901578157	1.00066853
560	1.000	46.32	0.880514584	1.000692952
560	1.000	49.21	0.860240462	1.000710396
560	1.000	52.11	0.840703738	1.000724352
560	1.000	55	0.822069244	1.000731329
580	1.000	0	1.153239266	1.064511
580	1.000	2.895	1.153239266	1.064511
580	1.000	5.789	1.153239266	1.064511
580	1.000	8.684	1.153239266	1.064511
580	1.000	11.58	1.152406439	1.064517978
580	1.000	14.47	1.149621675	1.064531933
580	1.000	17.37	1.144173601	1.064556355
580	1.000	20.26	1.135402896	1.064587755
580	1.000	23.16	1.122953874	1.064619155
580	1.000	26.05	1.107112818	1.064657532
580	1.000	28.95	1.088391571	1.06469242
580	1.000	31.84	1.067648983	1.064730798
580	1.000	34.74	1.045527027	1.064765686
580	1.000	37.63	1.022841178	1.064800575
580	1.000	40.53	0.999964473	1.064835463
580	1.000	43.42	0.977443454	1.064866863
580	1.000	46.32	0.9553909	1.064894774
580	1.000	49.21	0.934119121	1.064915707
580	1.000	52.11	0.913602091	1.06493664
580	1.000	55	0.89400464	1.064947107"""

#This function is used to parse the above table
def get_x_and_y(raw_data):
    lines = raw_data.split("\n")[1:]
    
    y1=[] #W
    y2=[] #Q
    x=[[],[],[]] #Tf,mdot,Tamb
    
    for line in lines:
        items = line.strip().split()
        for j in range(3):
            x[j].append(float(items[j]))
        y1.append(float(items[3]))
        y2.append(float(items[4]))
    return x,y1,y2

#parse the Hamilton data
x,y1,y2 = get_x_and_y(raw_data)


#linear regression commented out as not good enough
#note this used a single y variable
"""
import statsmodels.api as sm
def reg_m(y, x):
    x=np.array(x).T
    x=sm.add_constant(x)
    results = sm.OLS(endog=y, exog=x).fit()
    return results

results=reg_m(y, x)
print(results.summary())
plt.plot(results.predict())
"""

#second order. Predictive power of X-terms limited so remove 
#due to limited data from westinghouse
def fn(x,a,b,c,d,
       e,
       f,
       #g,
       #h,
       #i,
       j
       ):
    
    result=  a+b*x[0]+c*x[1]+d*x[2]
    result+= e*x[0]*x[1]
    result+= f*x[1]*x[2]
    #result+= g*x[0]*x[2]
    #result+= h*x[0]**2
    #result+= i*x[1]**2
    result+=j*x[2]**2
    return result
           
#Turn Hamilton Data into Numpy arrays          
x=np.array(x)
y1=np.array(y1)
y2=np.array(y2)

#Curve fit the hamilton data 1=W, 2=Q
def find_popt(x,y,fn):
    popt,pcov=curve_fit(fn,x,y)
    predict=fn(x,*popt)
    residual=np.sum((y-predict)**2)
    print(residual)
    return predict,popt

predict1,popt1 = find_popt(x,y1,fn)
predict2,popt2 = find_popt(x,y2,fn)



#Next step is to use Westinghouse data to CALCULATE these parameters
westinghouse_data=\
"""TF HTF mdot	Amb T	Wdot cycle	Heat In
633	1	        32.44	1	        1
573	1.121291922	32.44	1.018498602	0.999936849
653	0.964793198	32.44	1.0036567	1.000431529
633	0.999846209	20	    1.038287804	1
633	0.999846209	40	    0.982791998	1
633	0.417704599	32.44	0.455151646	0.499889486
573	0.417704599	40	    0.39599914	0.459999368"""

#Westinghouse curve fit gives reducing efficiecy with temperature at higher mass flows
#Consider using Carnot Scaling to fix this...
#def carnot(Thot,Tcold):
#    return 1-(Tcold+273.15)/(Thot+273.15)
#NOT WORKING COMMENT OUT
#westinghouse_data+=\
#"""653 1 32.44 {0:.4f} 1""".format(carnot(653,32.44)/carnot(633,32.44))

#Parse te Westinghouse Data and turn into arrays
wec_x,wec_y1,wec_y2 = get_x_and_y(westinghouse_data)
wec_x=np.array(wec_x)
wec_y1=np.array(wec_y1)
wec_y2=np.array(wec_y2)


#Generate Tabulation by (1) curve fitting model (2) using on same input params as Hamilton
def gen_wec_y(wec_x,wec_y,fn,x):
    popt,pcov=curve_fit(fn,wec_x,wec_y)
    new_gen_y = fn(x,*popt)
    return new_gen_y,popt

new_gen1 = None
new_gen2 = None
new_gen1,popt1_WEC = gen_wec_y(wec_x,wec_y1,fn,x)
new_gen2,popt2_WEC = gen_wec_y(wec_x,wec_y2,fn,x)

#Compare hamilton actual to predicted for all States for W and Q
#Also plot values with the Westinghouse prediction
def plot_compare(y,predict,new_gen):
    plt.plot(y,label="Hamilton")
    plt.plot(predict,label="Predicted")
    plt.plot(new_gen,label="Westinghouse")
    plt.legend(loc="lower right")
    plt.show()
    
plot_compare(y1,predict1,new_gen1)
plot_compare(y2,predict2,new_gen2)

#Produce the new table for Westinghouse
#NOTE: Curve fit is not producing good results, do not use at present!
with open("model1_wec.csv","w") as f:
    for j in range(x[0].shape[0]):
        f.write("{0:.0f},{1:.3f},{2:.0f},{3:.6f},{4:.6f},1,1\n".format(
            x[0][j],x[1][j],x[2][j],new_gen1[j],new_gen2[j]))

#Look at a contour map of the Hamilton and WEC Data
#The heat map for WEC data is not working, at high mass flows efficiency drops with temperature
mdot_contour=np.zeros((20,20))
T_contour=np.zeros((20,20))
wec_contour=np.zeros((20,20))
hamilton_contour=np.zeros((20,20))
for j in range(20):
    for k in range(20):
        mdot = 0.1*(j+1)
        T = 550+10*k
        T_contour[j,k]=T
        mdot_contour[j,k]=mdot
        x=np.array([[T],[mdot],[32.44]])
        wec_contour[j,k]=fn(x,*popt1_WEC)/fn(x,*popt2_WEC)
        hamilton_contour[j,k]=fn(x,*popt1)/fn(x,*popt2)
plt.imshow(wec_contour, cmap='hot', interpolation='nearest')
plt.show()
plt.imshow(hamilton_contour, cmap='hot', interpolation='nearest')
plt.show()

#Find d(nu)/dP
#THIS CURRENTLY USES THE HAMILTON DATA AS I DONT TRUST THE WEC CURVE FIT
relative_nu=[]
for j,T_CSP in enumerate([550,633]): #633 sets the temperatures both to 633 for reference

    propn_LFR_flow = 0.8 #Proportion of turbine rating delivered by LFR only at full power . Can change   
    T_amb = 32.44 #C Ambient
    wec_P = 465.0 #MW LFR electrical power
    design_P = wec_P/propn_LFR_flow #MW 
    T_LFR = 633 #C source: Westinghouse PEPSE
    pressure = 33 #MPa source: westinghouse
    
    mdot_lfr_only = wec_P/design_P #mdot at which Salt2Steam kicks in
    
    #Steam properties
    lfr_steam = IAPWS97(T=T_LFR+273.15,P=pressure)
    csp_steam = IAPWS97(T=T_CSP+273.15,P=pressure)
    
    x=[[],[],[]] #Tf,mdot,T_amb
    mdots = [0.1*j for j in range(1,20)]
    for mdot in mdots:
        #calculate the temperature of the hot stream
        x[1].append(mdot)
        x[2].append(T_amb)
        if mdot<=mdot_lfr_only:
            x[0].append(T_LFR)
        else:
            #Mass-flow weighted enthalpy to turbine - then get temperature
            h = (lfr_steam.h*mdot_lfr_only+csp_steam.h*(mdot-mdot_lfr_only))/mdot
            T_mix = IAPWS97(h=h,P=pressure).T-273.15
            x[0].append(T_mix)
    
    #Plot mass flow rate vs efficiency
    x = np.array(x)
    SAM_gen1 = fn(x,*popt1)
    SAM_gen2 = fn(x,*popt2)
    relative_nu.append(SAM_gen1/SAM_gen2)
    
    
plt.plot(x[1],relative_nu[-1],label="T_CSP=633 (for reference only)")
plt.plot(x[1],relative_nu[0],label="T_CSP=550 (actual case)")
plt.ylabel("efficiency multipier")
plt.xlabel("mdot")
plt.legend(loc='lower right')
plt.show()

"""
#Check efficiencies at design point
nu_633 = fn(np.array([[633],[1],[32.44]]),*popt1)/fn(np.array([[633],[1],[32.44]]),*popt2)
nu_550 = fn(np.array([[550],[1],[32.44]]),*popt1)/fn(np.array([[550],[1],[32.44]]),*popt2)
print(nu_633,nu_550)
"""


#Take gradients of efficiencies for +/- 10% at design point with different temperature
#Do it at 550 and 633 to check sensitivity. 
idx = int(10*propn_LFR_flow)-1
print(relative_nu[0][idx]-relative_nu[0][idx-1],relative_nu[0][idx+1]-relative_nu[0][idx])
print(relative_nu[1][idx]-relative_nu[1][idx-1],relative_nu[1][idx+1]-relative_nu[1][idx])    
    