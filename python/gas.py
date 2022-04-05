
import numpy
import matplotlib.pyplot as plt
import tannehill as tan

def func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3):

    GAS10=1./(1.+GAS7)
    GAS8=GAS3-GAS4*Z2
    GAS9=GAS8*GAS7*GAS10**2
    GAMM =GAS1-GAS2*Z2-GAS8*GAS10
    GAME=2.304*(-GAS2+GAS4*GAS10+GAS9*DERE)
    GAMR=2.304*(GAS5+GAS6*GAS10+GAS9*DERR)
    SNDSQ=E*(A1+(GAMM-1.)*(GAMM+A2*GAME)+A3*GAMR)
    A=numpy.sqrt(SNDSQ)
    P=RHO*E*(GAMM-1.)    

    return P, A


def T_TGAS(RHO, E):

    Y2= numpy.log10(RHO/1.292)
    Z2= numpy.log10(E/78408.4)
    
    P, A = p_TGAS(RHO, E)  
    X2=numpy.log10(P/1.0134E+05)
    Y2=Y2+.0231264 
    Z3=X2-Y2

    print(Y2, Z3,':')

    if(Y2>-.50):
        # block 29
        if(Z3<=.48):
            T = P/(287.*RHO)

        elif(Z3<=1.07):
            T=10**(.279268+.992172*Z3)
            T= T*151.777778

        else:
            T=10.**(.2332605-.056383*Y2+(1.19783+.063121*Y2-.165985*Z3)*Z3+(-.814535+.099233*Y2+(.602385-.067428*Y2-.098991*Z3)*Z3)/( 1.+numpy.exp((5.*Y2-20.)*(Z3-1.78) ) ) )
            T= T*151.777778

    elif(Y2>-4.50):
        # block 24
        if(Z3<=.48):
            T = P/(287.*RHO)

        elif(Z3<=.9165):
            T=10.**(.284312+.987912*Z3+.001644*Y2)
            T= T*151.777778

        elif(Z3<=1.478):
            T=10.**(.502071-.01299*Y2+(.774818+.025397*Y2)*Z3+(.009912\
            -.150527*Y2+(-.000385+.105734*Y2)*Z3)/(1.+numpy.exp(-15.*(Z3-1.28))))
            T= T*151.777778

        elif(Z3<=2.176):
            T=10.**(1.02294+.021535*Y2+(.427213+.006900*Y2)*Z3+ (-.427823\
            -.211991*Y2+(.257096+.101192*Y2)*Z3)/(1.+numpy.exp(-12.*(Z3-1.778))))
            T= T*151.777778

        else:
            T=10.**(1.47540+.12962*Y2+(.254154-.046411*Y2)*Z3+(-.221229\
            -.057077*Y2+(.158116+.03043*Y2)*Z3)/( 1.+numpy.exp(5.*Y2*(Z3-2.40))))
            T= T*151.777778

    else:

        if(Z3>.30):
            # block 19
            if(Z3<=1.00):               
                T = 10.**(.2718+.00074*Y2+(.990136-.004947*Y2)*Z3+(.990717\
                +.175194*Y2-(.982407+.159233*Y2)*Z3)/(1.+numpy.exp(-20.*(Z3-0.88))))
                T = T*151.777778

            elif(Z3<=1.35):
                T=10.**(1.39925+.167780*Y2+(-.143168-.159234*Y2)*Z3+(-.027614\
                -.090761*Y2+(.307036+.121621*Y2)*Z3)/(1.+numpy.exp(-20.*(Z3-1.17))))
                T= T*151.777778

            elif(Z3<=1.79): # Faixa problema
                T=10.**(1.11401+.002221*Y2+(.351875+.017246*Y2)*Z3+(-1.15099\
                -.173555*Y2+(.6733424+.088399*Y2)*Z3)/(1.+numpy.exp(-20.*(Z3-1.56))))
                T= T*151.777778

                """
                Y = Y2
                Z = Z3

                ans = b[0] + b[1]*Y + b[2]*Z
                ans += b[3]*Y*Z + b[4]*Z*Z + b[5]*Y*Y
                ans += b[6]*Y*Y*Z + b[7]*Y*Z*Z
        
                aux2 = b[8] + b[9]*Y + b[10]*Z + b[11]*Y*Z + b[12]*Z*Z
                aux2 /= 1 + numpy.exp((b[13]*Y + b[14])*(Z + b[15]))

                ans += aux2
    
                Taux = 151.78*(10.**ans)

                print("oi", T, Taux)

                """

            elif(Z3<=2.47):
                T= 10.**(1.01722-.017918*Y2+(.473523+.025456*Y2)*Z3+(-2.17978\
                -.334716*Y2+(.898619+.127386*Y2)*Z3)/(1.+numpy.exp(-20.*(Z3-2.22))))
                T= T*151.777778

            else:
                T=10.**(-45.0871-9.00504*Y2+(35.8685+6.7922*Y2)*Z3-(6.77699\
                +1.2737*Y2)*Z3*Z3+(-.064705+.025325*Z3)*Y2*Y2)
                T= T*151.777778

        else:
            T = P/(287.*RHO)

    return T



def p_TGAS(RHO, E):

    Y2= numpy.log10(RHO/1.292)
    Z2= numpy.log10(E/78408.4)

    if(Y2> -0.5):
        # block 11
        
        if(Z2<=.65):
            GAMM=1.400
            SNDSQ=E*.560
            A = numpy.sqrt(SNDSQ)
            P = RHO*E*(GAMM-1.0)

        elif(Z2<=1.68):

            GAS1=1.45510-.000102*Y2
            GAS2=.081537-.000166*Y2
            GAS3=-.128667+.049454*Y2
            GAS4=-.101036+.033518*Y2
            GAS5=-.000102+.000166*Z2
            GAS6=-.049454+.033518*Z2
            GAS7=numpy.exp(-15.*(Z2-1.420))
            DERE=-15.
            DERR=0.
            A1=.000450
            A2=.203892
            A3=.101797
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

        elif(Z2<=2.46):

            GAS1=1.59608-.042426*Y2
            GAS2=.192840-.029353*Y2
            GAS3=.019430-.005954*Y2
            GAS4=.026097-.006164*Y2
            GAS5=-.042426+ .029353*Z2
            GAS6=.005954-.006164*Z2
            GAS7=numpy.exp(-15.*(Z2-2.050))
            DERE=-15.
            DERR=0.0
            A1=-.006609
            A2=.127637
            A3=.297037
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)
  
        else:

            GAS1=1.54363-.049071*Y2
            GAS2=.153562-.029209*Y2
            GAS3=.324907+.077599*Y2
            GAS4=.142408+.022071*Y2
            GAS5=-.049071+.029209*Z2
            GAS6=-.077599+.022071*Z2
            GAS7=numpy.exp(-10.*(Z2-2.708))
            DERE=-10.0
            DERR=0.0
            A1=-.000081
            A2=.226601
            A3=.170922
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)           


    elif(Y2> -4.5):
        # block 6

        if(Z2<=.65):
            GAMM=1.400
            SNDSQ=E*.560
            A = numpy.sqrt(SNDSQ)
            P = RHO*E*(GAMM-1.0)

        elif(Z2<=1.54):
            GAS1 =1.44813+.001292*Y2
            GAS2=.073510+.001948*Y2
            GAS3=-.054745+.013705*Y2
            GAS4=-.055473+.021874*Y2
            GAS5=.001292-.001948*Z2
            GAS6=-.013705+.021874*Z2
            GAS7=numpy.exp(-10.*(Z2-1.42))
            DERE=-1.0
            DERR=0.0
            A1=-.001973
            A2=.185233
            A3=-.059952
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

        elif(Z2<=2.22):
            GAS1 = 1.73158+.003902*Y2
            GAS2=.272846-.006237*Y2
            GAS3=-.041419-.037475*Y2
            GAS4=.016984-.018038*Y2
            GAS5=.003902+.006237*Z2
            GAS6=.037475-.018038*Z2
            GAS7=numpy.exp((-10.+3.0*Y2)*(Z2-.025*Y2-2.025))
            DERE=3.0*Y2-10.0
            DERR=3.0*Z2+12.15*Y2-20.325
            A1=-.013027
            A2=.074270
            A3=.012889
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

        elif(Z2<=2.90):
            GAS1= 1.59350+.075324*Y2
            GAS2=.176186+.026072*Y2
            GAS3=.200838+.05853*Y2
            GAS4=.099687+.025287*Y2
            GAS5=.075326-.026072*Z2
            GAS6=-.058536+.025287*Z2
            GAS7=numpy.exp(-10.0*Z2+(5.0*Z2-13.5)*Y2+27.0)
            DERE=5.0*Y2-10.0
            DERR=5.0*Z2-13.5
            A1=.004362
            A2=.212192
            A3=-.001293
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

        else:
            GAS1=1.12688-.025957*Y2
            GAS2=-.013602-.013772*Y2
            GAS3=.127737+.087962*Y2
            GAS4=.043104+.023547*Y2
            GAS5=-.025957+.013772*Z2
            GAS6=-.087942+.023547*Z2
            GAS7=numpy.exp(-20.0*Z2+(4.0*Z2-13.2)*Y2+66.0)
            DERE=-20.+4.0*Y2
            DERR=4.0*Z2-13.2
            A1=.006368
            A2=.209716
            A3=-.006001
          
            P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

    else:

        if(Z2> 0.65):
            # block 1

            if(Z2<=1.5):            
                GAMM=1.46543+(.007625+.000292*Y2)*Y2-(.254500+.017244*Y2)*Z2+(.355907+.015422*Y2-.163235*Z2)*Z2*Z2
                GAME=2.304*(-.25450-.017244*Y2+(.711814+.030844*Y2-.489705*Z2)*Z2)
                GAMR=2.304*(.007625+(-.017244+.015422*Z2)*Z2+.000584*Y2)
                A1=-.000954
                A2=.171187
                A3=.004567

                SNDSQ=E*(A1+(GAMM-1.)*(GAMM+A2*GAME)+A3*GAMR)
                A=numpy.sqrt(SNDSQ)
                P=RHO*E*(GAMM-1.) 

            elif(Z2<=2.20):
                GAS1=2.02636+.058493*Y2
                GAS2=.454886+.027433*Y2
                GAS3=.165265+.014275*Y2
                GAS4=.136685+.010071*Y2
                GAS5=.058493-.027433*Z2
                GAS6=-.014275+.010071*Z2
                GAS7=numpy.exp(0.285*Y2-30.0*Z2+58.41)
                DERE=-30.0
                DERR=0.285
                A1=.008737
                A2=.184842
                A3=-.302441                

                P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

            elif(Z2<=3.05):
                GAS1=1.60804+.034791*Y2
                GAS2=.188906+.010927*Y2
                GAS3=.124117+.007277*Y2
                GAS4=.069839+.003985*Y2
                GAS5=.034791-.010927*Z2
                GAS6=-.007277+.003985*Z2
                GAS7=numpy.exp(0.21*Y2-30.0*Z2+80.73)
                DERE=-30.0
                DERR=0.21
                A1=.017884
                A2=.153672
                A3=-.930224

                P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

            elif(Z2<=3.38):
                GAS1=1.25672+.007073*Y2
                GAS2 =.039228-.000491*Y2
                GAS3=-.721798-.073753*Y2
                GAS4=-.198942-.021539*Y2
                GAS5=.007073+.000491*Z2
                GAS6=.073753-.021539*Z2
                GAS7=numpy.exp(0.425*Y2-50.0*Z2+166.7)
                DERE=-50.0
                DERR=0.325
                A1=.002379
                A2=.217959
                A3=.005943

                P, A = func1(RHO, E, Z2, GAS1, GAS2, GAS3, GAS4, GAS5, GAS6, GAS7, DERE, DERR, A1, A2, A3)

            else:

                GAMM=-84.0327+(-.831761+.001153*Y2)*Y2+(72.2066+.491914*Y2)*Z2+(-20.3559-.070617*Y2+1.90979*Z2)*Z2*Z2
                GAME=2.304*(72.2066+.491914*Y2+(-40.7118-.141234*Y2+5.72937*Z2)*Z2)
                GAMR=2.304*(-.831761+.002306*Y2+(.491914-.070617*Z2)*Z2)
                A1=.006572
                A2=.183396
                A3=-.135960

                SNDSQ=E*(A1+(GAMM-1.)*(GAMM+A2*GAME)+A3*GAMR)
                A=numpy.sqrt(SNDSQ)
                P=RHO*E*(GAMM-1.) 
                

        else:
            GAMM = 1.4
            SNDSQ = E*.560
            A = numpy.sqrt(SNDSQ)
            P = RHO*E*(GAMM-1.0)

    
    return P, A



def h_TGAS(RHO, P):

    Y2= numpy.log10(RHO/1.292)
    X2= numpy.log10(P/1.013E+05)
    Z3=X2-Y2
    if(Y2>-0.50):
    # block 10

        if(Z3<=0.3):
            H = (P/RHO)*3.50

        elif(Z3<=1.15):
            GAS1=1.42598+ 0.000918*Y2
            GAS2=.092209+.002226*Y2
            GAS3=-.0197721+.036600*Y2
            GAS4=-.0774694+.043878*Y2
            GAS5=numpy.exp(-15.00*(Z3-1.040))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

        elif(Z3<=1.6):
            GAS1=1.64689-.0621547*Y2
            GAS2=.334994-.0636120*Y2
            GAS3=.0383322+.0144677*Y2
            GAS4=.0734214-.0024417*Y2
            GAS5=numpy.exp(-15.00*(Z3-1.360))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

        else:
            GAS1=1.48558-.453562*Y2
            GAS2=.152096-.303350*Y2
            GAS3=.459282-.448395*Y2
            GAS4=.220546-.292293*Y2
            GAS5=numpy.exp(-10.00*(Z3-1.600))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

    elif(Y2>-4.50):
    # block 5

        if(Z3<=0.3):
            H = (P/RHO)*3.50

        elif(Z3<=0.98):
            GAS1=1.42176-.000366*Y2
            GAS2=.083614-.000677*Y2
            GAS3=-.005272+.115853*Y2
            GAS4=-.007363+.146179*Y2
            GAS5=numpy.exp(-20.00*(Z3-0.860))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

        elif(Z3<=1.380):
            GAS1=1.74436-.035354*Y2
            GAS2=.415045-.061921*Y2
            GAS3=-.018536-.043582*Y2
            GAS4=.0443534-.049750*Y2
            GAS5=numpy.exp(-20.00*(X2-1.04*Y2-1.336))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

        elif(Z3<=2.04):
            GAS1=1.49674-.021583*Y2
            GAS2=.197008-.030886*Y2
            GAS3=.157738+.009158*Y2
            GAS4=.123213-.006553*Y2
            GAS5=numpy.exp(-10.00*(X2-1.05*Y2-1.895))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
            H=(P/RHO)*(GAMM/(GAMM-1.))

        else:
            GAS1=1.10421-.033664*Y2
            GAS2=-.031768-.024335*Y2
            GAS3=.178802+.017456*Y2
            GAS4=.080373+.002511*Y2
            GAS5=numpy.exp(-15.00*(X2-1.08*Y2-2.650))
            GAS10= 1./(1.+GAS5)
            GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10        
            H=(P/RHO)*(GAMM/(GAMM-1.))

    else:

        if(Z3 > 0.398):
        #block 1

            if(Z3<=0.87):
                GAS1=1.47003+.007939*Y2
                GAS2=.244205+.025607*Y2
                GAS3=-.872248-.049452*Y2
                GAS4=-.764158+.000147*Y2
                GAS5=numpy.exp(-20.00*(Z3-0.742))
                GAS10= 1./(1.+GAS5)
                GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
                H=(P/RHO)*(GAMM/(GAMM-1.))

            elif(Z3<=1.270):
                GAS1= 3.18652+.137930*Y2
                GAS2= 1.89529+.103490*Y2
                GAS3= 2.14572+.272717*Y2
                GAS4= 2.06586+.223046*Y2
                GAS5=numpy.exp(-15.00*(Z3-1.041))
                GAS10= 1./(1.+GAS5)
                GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
                H=(P/RHO)*(GAMM/(GAMM-1.))

            elif(Z3<=1.863):
                GAS1= 1.63963-.00100436*Y2
                GAS2= .303549-.0164639*Y2
                GAS3= .852169+.101237*Y2
                GAS4= .503123+.0435801*Y2
                GAS5=numpy.exp(-10.00*(Z3-1.544))
                GAS10= 1./(1.+GAS5)
                GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
                H=(P/RHO)*(GAMM/(GAMM-1.))

            else:
                GAS1 = 1.55889+.0559323*Y2
                GAS2= .211764+.0235478*Y2
                GAS3= .549041+.101758*Y2
                GAS4= .276732+.0460305*Y2
                GAS5=numpy.exp(-15.00*(Z3-2.250))
                GAS10= 1./(1.+GAS5)
                GAMM=GAS1-GAS2*Z3-(GAS3-GAS4*Z3)*GAS10
                H=(P/RHO)*(GAMM/(GAMM-1.))

        else:

            H = (P/RHO)*3.5

    return H

def checkEntalpy():
    
    er = numpy.arange(-7., 3., 1.)
    ep = numpy.arange(-8, 4., 0.1)
    
    rr = 10.0**er
    pp = 10.0**ep
                    
    for ii in range(0, len(rr)):
        for jj in range(0, len(pp)):
            error = tan.entalpy(1.292*rr[ii], 1.0134e5*pp[jj])/h_TGAS(1.292*rr[ii], 1.0134e5*pp[jj])-1
            if abs(error) > 1e-3:
                print(er[ii], ep[jj]-er[ii], error)

    return None        

def checkPressure():

    er = numpy.arange(-7., 3., 1.)
    E = numpy.arange(-1., 3., 0.1)
    
    rr = 10.0**er
    ee = 10.0**E
                
    for r in rr:
        for e in ee:
            P, A = p_TGAS(1.292*r, 78408.4*e)            
            error = tan.pressure(1.292*r, 78408.4*e)/P-1
            if(error>1e-5):
                print(error)
        
    return None

def checkSound():

    er = numpy.arange(-7., 3., 1.)
    E = numpy.arange(-1., 3., 0.1)
    
    rr = 10.0**er
    ee = 10.0**E
                
    for r in rr:
        for e in ee:
            P, A = p_TGAS(1.292*r, 78408.4*e)            
            error = tan.sound(1.292*r, 78408.4*e)/A-1
            if(error>1e-1):
                print(error)
        
    return None

def checkTemperature():
    
    er = numpy.arange(-7., 3., 1.)
    E = numpy.arange(-1., 3.5, 0.1)
    
    rr = 10.0**er
    ee = 10.0**E
                
    for r in rr:
        for e in ee:
            a, b = tan.temperature(1.292*r, 78408.4*e), T_TGAS(1.292*r, 78408.4*e)
            error = a/b -1
            print("vish", error)

    return None

def test():

    x = numpy.arange(-7., 1., 0.5)    
    y = numpy.arange(0.5, 3, 0.25)
    
    X, Y = numpy.meshgrid(y, x)
    
    error = X*0.0

    rr = x*0
    for ii in range(0, len(rr)):    
        rr[ii] = 1.292*(10.**x[ii])
        
    ee = y*0
    for ii in range(0, len(ee)):    
        ee[ii] = 78408.4*(10.**y[ii])
 
    for ii in range(0, len(rr)):
        for jj in range(0, len(ee)):
            p, A = p_TGAS(rr[ii], ee[jj])
            h = h_TGAS(rr[ii], p)
            error[ii][jj] = (h - p/rr[ii])/ee[jj] - 1
    
    plt.figure()
    plt.title('error')        
    plt.contourf(X, Y, error)
    plt.colorbar()
    plt.show()

if __name__=="__main__":

    
#    checkEntalpy()
#    checkPressure()
#    checkSound()
    checkTemperature()
#    test()





