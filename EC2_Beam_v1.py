import math
import numpy as np
from sympy import *
from math import asin,tan


def As_rq(fck,fyk,M_Ed,depth,width):
    d = depth - 70
    M_abs = abs(M_Ed)   #kNm
    pi = 3.1416
    fcd = fck*1/1.5
    fctm = 0.3*fck**0.666
    As_min = max((0.26*fctm/fyk),0.0013)*width*d

    K = M_abs*1e6/(width*d**2*fck)
    K_=0.196
    #print ("K= ",round(K,2))

    #Lever arm
    if K < K_:
        z = d*(0.5 + (0.25-K*0.75)**0.5)
        Asc=0
        As = max(M_abs*1e6/(z*fyk/1.15),As_min)
        ratio = As*100/(width*d)

    else:
        #Doubly reinforced K>K_
        z = d*(0.5 + (0.25-K_*0.75)**0.5)
        d2 = 35 + 10 + 32/2
        x = 2.5*(d - z)
        #print ("Beam is doubly reinforced")
        fsc = min(fyk/1.15,700*(x-d2)/x)
        #print ("fsc=",round(fsc,1))
        # Calculating the area of compression steel
        Asc = ((K-K_)*fck*width*d**2)/(fsc*(d-d2))
        # Calculating the area of tensile steel
        As = max(((K_*fck*width*d**2)/(fyk/1.15*z) + Asc),As_min)
        Asmax = 0.04*width*d
        ratio = As*100/(width*d)
        """
        if As < Asmax:
            print ("Asc = ",round(Asc,1))
            print ("As = ",round(As,1))
            print ("ratio = ", round(ratio,2), "%")
        else: 
            print ("Beam is failed")
        """
    if M_Ed> 0:
        As_top=As
        As_bot=Asc
    else:
        As_top=Asc
        As_bot=As
    return As,Asc,As_top,As_bot,K,z


def Shear_cal(fck,fywk,V_Ed,H,B):
    fcd=fck/1.5
    fywd=fywk/1.15
    bw=B
    d=H-60
    Asw_s_min = 0.08*bw*(fck**0.5) / fywk #sina = 1
    v_Ed= abs(V_Ed)*1000/(bw*0.9*d)
    v_Rdmax=0.6*(1-fck/250)*fcd/2.9
    Shear_status=0
    Asw_s=0
    theta = 0.5*asin(v_Ed/(0.2*fck*(1-fck/250)))
    if theta == 0:
        theta = 0.01
    else:
        theta = 0.5*asin(v_Ed/(0.2*fck*(1-fck/250)))
    cot_theta=cot(math.radians(theta))
    if v_Ed>v_Rdmax:
        Shear_status=1

    else:
        if cot_theta>2.5:
            Asw_s= max((v_Ed*bw/fywd/2.5),Asw_s_min)
        elif cot_theta>=1:
            Asw_s= max((v_Ed*bw/fywd/cot_theta),Asw_s_min)
    return Asw_s,v_Ed, cot_theta, Shear_status
#print(As_rq(30,500,184.4,450,300))
#print(As_rq(32,500,4965,1200,800))
#print(Shear_cal(30,500,164.5,450,300))





