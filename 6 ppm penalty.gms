$Title Multi-product Batch Processing Case 20250915

*============================================================
Set
    J / 1*5 /
    I / 1*4 /
    S / 1*19 /
    T / t0*t1007 /
    N / 1*15 /
    Hr / h0*h23 /;
    
Alias (S, s1);
Alias (J, Jp);
Alias (T, Tp);


Set P(J,N) /
      1.1, 1.2, 1.3, 1.4
      2.5, 2.6, 2.7, 2.8
      3.9, 3.10, 3.11, 3.12
      4.13,4.14, 4.15
/;

Parameter d(J) / 1*5  0.441 /;

Set Im / m1*m10 /;
Parameter Mbp(Im) / m1 0, m2 5, m3 10, m4 15, m5 20, m6 25, m7 30, m8 45, m9 60, m10 100/;

Scalar alpha / 0.157 /;
Parameter Phi_bp(Im);
loop(Im, Phi_bp(Im) = exp(-alpha * Mbp(Im)); );

Scalar mMax / 100 /;
Scalar cMax_dra / 100 /;
Scalar phiMin;       phiMin = exp(-alpha*mMax);

*------------------------------------------------------------
* Sliding-horizon controls (constant-size window)
Scalar HorizonLen / 144 /;
Scalar ShiftLen   / 144 /;
Set    Window / w1*w7 /;

Scalar tStartIdx, tEndIdx;
Set Twindow(Window,T) "precomputed active time steps per window";
Twindow(Window,T)$((ord(T) >= 1 + (ord(Window)-1)*ShiftLen) and (ord(T) <= (ord(Window)-1)*ShiftLen + HorizonLen)) = yes;

Set
    TactVar(T)   "active time steps (variables allowed)"
    TactEq(T)    "time steps where dynamic equations are enforced"
    Tinit(T)     "first point of current window"
    Tterm(T)     "last point of current window (state carried forward)";

Parameter
    carry_M(J)        "pipeline DRA mass state carried between windows"
    carry_cin(J)      "inlet concentration state carried"
    carry_cout(J)     "outlet concentration state carried"
    carry_cstart(J)   "start concentration state carried"
    carry_cend(J)     "end concentration state carried"
    carry_K(J,N)      "pump status carried (previous step)"
    cbar0(J)          "initial average concentration per segment";

* Initial carry-over values (t0)
*carry_M(J)      = 50000;
carry_M(J)      = 40000;
carry_K(J,N)    = 0;



*============================================================

Scalar
    dt / 600 /
    g / 9.81 /
    C2 / 300 /
    C3 / 0.006 /
    W  / 0 /
    beta  / 2e-6 /
    m1    / 1 /;


Parameter zStation(J) /
      1   200
      2   400
      3  1150
      4  1700
      5  1500
/;

Parameter Xj(J) /
      1  0
      2  18823
      3  31765
      4  43551
      5  64935.
/;

Parameter PinMin(J) /
      1  344700
      2  344700
      3  344700
      4  344700
      5  110000
/;

Parameter PinMax(J) /
      1  12000000
      2  12000000
      3  12000000
      4  12000000
      5  1400000
/;

Parameter    MI(J), MP(J), ML(J), MDL(J), MA(J);
loop(J,    MI(J)  = no; MP(J)  = no; ML(J)  = no; MDL(J)  = no; MA(J) = no;);
MI('1') = yes;
MP('1') = yes; MP('2') = yes; MP('3') = yes; MP('4') = yes;
ML('1') = yes; ML('2') = yes; ML('3') = yes; ML('4') = yes;
MDL('4') = yes;
ML('1') = yes; ML('2') = yes; ML('3') = yes; ML('4') = yes;
MA('1') = yes; MA('2') = yes; MA('3') = yes; MA('4') = yes;

Parameter C1_hour(Hr);
loop(Hr, C1_hour(Hr) = 0.12;)
loop(Hr$(ord(Hr) > 5   and ord(Hr) <= 12), C1_hour(Hr) = 0.40;);
loop(Hr$(ord(Hr) > 12  and ord(Hr) <= 16), C1_hour(Hr) = 0.20;);
loop(Hr$(ord(Hr) > 16  and ord(Hr) <= 20), C1_hour(Hr) = 0.25;);

Parameter hourMap(T,Hr);
loop((T,Hr)$(ord(Hr) = mod(floor((ord(T)-1)*dt/3600), 24) + 1), hourMap(T,Hr) = 1;);

Parameter a1(N) /
      1   1264,  2   1264,  3   964,  4  564
      5   1500,  6   1200,  7   1053,  8  522
      9   1547,  10   1247,  11   1047,   12   525
      13   1100, 14   500, 15   500
/;
Parameter b1(N) / 1*15 23800 /;
Parameter c1_1(N) / 1*15 -0.05 /;
Parameter c2_1(N) / 1*15 0.85 /;
Parameter c3_1(N) / 1*15 0.6 /;

Parameter X0(I) /
      1  64935
      2  33552
      3  -46360
      4  -86360
/;

Parameter rhoBatch(I) /
      1  860
      2  740
      3  860
      4  740
/;

Parameter durSched(S) /
       1  6.00,        2  9.67,        3 15.00,       4 12.00,        5 16.83,        6 11.00,
       7  6.17,        8 10.00,        9  9.00,       10  5.17,      11  8.00,       12  9.00,
       13  3.67,       14  7.00,       15  6.00,       16  5.17,       17 10.00,       18 15.00,      19 10.67
/;

Parameter QinSched(S) /
  1 0.19444,  2 0.16667,  3 0.19444,  4 0.18056,
  5 0.18056,  6 0.18056,  7 0.16667,  8 0.15278,
  9 0.19444, 10 0.18056, 11 0.18056, 12 0.17222,
 13 0.16667, 14 0.17222, 15 0.18333, 16 0.19167,
 17 0.20000, 18 0.20000, 19 0.18889
/;

Parameter QoutSched(S) /
  1 0.00000,  2 0.00000,  3 0.05556,  4 0.05556,
  5 0.04167,  6 0.04167,  7 0.04167,  8 0.00000,
  9 0.00000, 10 0.00000, 11 0.00000, 12 0.00000,
 13 0.00000, 14 0.00000, 15 0.00000, 16 0.00000,
 17 0.05000, 18 0.04167, 19 0.00000
/;



*============================================================
*  pre-ccalculate parameters
Parameter
    cumDur(S)
    cumDurSec(S)
    prevDurSec(S)
    volCum(T);

loop(S, cumDur(S) = sum(s1$(ord(s1) le ord(S)), durSched(s1)););
loop(S, cumDurSec(S) = cumDur(S)*3600;);
prevDurSec('1') = 0;
loop(S$(ord(S)>1), prevDurSec(S) = sum(s1$(ord(s1)=ord(S)-1), cumDurSec(s1)););

Scalar
    currentTime
    sSel
    timeInSeg;

Loop(T,
    currentTime = (ord(T)-1) * dt;
    sSel = 1;
    Loop(S$(cumDurSec(S) < currentTime), sSel = ord(S) + 1;);
    if(sSel > card(S), sSel = card(S););
    timeInSeg = currentTime - sum(s1$(ord(s1)=sSel), prevDurSec(s1));
    if(timeInSeg < 0, timeInSeg = 0;);
    volCum(T) = sum(s1$(ord(s1) < sSel), durSched(s1)*3600 * QinSched(s1)) +
                sum(s1$(ord(s1)=sSel), QinSched(s1)) * timeInSeg;
);

Parameter
    Q_out(T, J)
    Qs(T,J)
    Qd(T,J);

loop((T,J), Q_out(T, J) = 0;);

Scalar
    injRate
    outRate;

Loop(T,
    currentTime = (ord(T)-1) * dt;
    sSel = 1;
    Loop(S$(cumDurSec(S) < currentTime), sSel = ord(S) + 1;);
    if(sSel > card(S), sSel = card(S));
    injRate = sum(S$(ord(S)=sSel), QinSched(S));
    outRate = sum(S$(ord(S)=sSel), QoutSched(S));
    Loop(J, Q_out(T,J) = outRate / card(J););
    Loop(J,
        Qs(T,J) = injRate - sum(Jp$(ord(Jp) < ord(J)), MDL(Jp) * Q_out(T,Jp));
        Qd(T,J) = Qs(T,J) - MDL(J) * Q_out(T,J);
    );
);

*Display cumDur, cumDurSec, prevDurSec, volCum;
*Display Qs, Qd, Q_out;

Parameter X(T,I);
Loop(T,Loop(I,X(T,I) = X0(I) + volCum(T);););
*Display X;

Parameter Da(T,I,J), Db(T,I,J), D0(T,I,J), Dc(T,I,J), Dd(T,I,J), De(T,I,J), Df(T,I,J);
Loop(T, Loop(I, Loop(J, if (X(T,I) > Xj(J),   Da(T,I,J) = 1;   else   Da(T,I,J) = 0;););););
Loop(T, Loop(I, Loop(J, Db(T,I,J) = Da(T,I,J) - Da(T,I+1,J););););
Loop(T, Loop(I, Loop(J, D0(T,I,J) = (1-Da(T,I,J))*Da(T,I,J-1););););
Loop(T, Loop(I, Loop(J, Dc(T,I,J) = Da(T,I,J-1)*(1-Da(T,I,J))*(1-Da(T,I+1,J-1)););););
Loop(T, Loop(I, Loop(J, Dd(T,I,J) = Da(T,I,J-1)*(1-Da(T,I,J))*Da(T,I+1,J-1)*(1-Da(T,I+1,J)););););
Loop(T, Loop(I, Loop(J, De(T,I,J) = Da(T,I,J-1)*Da(T,I,J)*(1-Da(T,I+1,J-1)););););
Loop(T, Loop(I, Loop(J, Df(T,I,J) = Da(T,I,J-1)*Da(T,I,J)*Da(T,I+1,J-1)*(1-Da(T,I+1,J)););););

parameter eff(T,J,N);
loop(T, loop((J,N)$P(J,N),  eff(T,J,N) = c1_1(N)*Qs(T,J)**2 + c2_1(N)*Qs(T,J) + c3_1(N); ););
*display eff;
parameter Rho_p(T,J);
loop(T, loop(J,  Rho_p(T,J) = sum(I, rhoBatch(I)*Db(T,I,J)); ););
*display Rho_p;
Parameter V_seg(T,J);
loop((T,J),  V_seg(T,J) = Qs(T,J) / (pi/4 * d(J)**2););
Parameter Dseg(T,I,J);
loop((T,I,J)$(ord(J)<card(J)), Dseg(T,I,J) = Da(T,I,J) - Da(T,I,J+1););
Parameter zBatch(T,I);
loop((T,I), zBatch(T,I) =  sum(J$(ord(J)<card(J)),  Dseg(T,I,J)* ( zStation(J) + (zStation(J+1)-zStation(J)) * ( X(T,I) - Xj(J) ) / ( Xj(J+1) - Xj(J) ))););
*Display zBatch;
Parameter L(J);
loop(J$(ord(J) < card(J)),  L(J) =    (Xj(J+1) - Xj(J)) / (pi/4 * d(J)**2););

Parameter delay(T,J);
delay(T,J) = ceil( L(J) / (max(V_seg(T,J), eps) * dt) );

Parameter shiftMap(T,Tp,J);
loop((T,Tp,J)$(ord(T) = ord(Tp) + delay(Tp,J)), shiftMap(T,Tp,J) = 1;);

Parameter shiftActive(T,J);

*display delay;

parameter maxT;         maxT=100000;


Parameter epsilon(J), mu;
epsilon(J) = 4.5e-5;
mu = 0.01;

Parameter Re(T,J), f(T,J);
Re(T,J) = Rho_p(T,J)*V_seg(T,J)*d(J)/mu;
f(T,J) = 0.02;

Set it /i1*i4/;
loop(it, f(T,J)$(Re(T,J)>0) = 1 / sqr(-2*log((epsilon(J)/(3.7*d(J))) + 2.51/(Re(T,J)*sqrt(f(T,J))))/log(10)); );

display f;

Set HighPriceT(T) "High electricity price time steps";

HighPriceT(T) = yes$(sum(Hr$(hourMap(T,Hr) and C1_hour(Hr) >= 0.40), 1));



*============================================================
Binary Variable
    K(T,J,N),
    SwitchPos(T,J,N)
    SwitchNeg(T,J,N),
    AnyPumpOn(T,J);    

Positive Variable
    h(T,J,N)
    Hj(T,J)
    p_s(T,J)
    p_d(T,J)
    p_es(T,J)
    p_ed(T,J)
    p_p(T,J)
    p_Lf(T,J),
    c_inject(T,J)
    c_start(T,J)
    c_end(T,J)
    c_in(T,J)
    c_out(T,J),
    draFactor(T,J),
    yprod(T,J),
    m_lin(T,J),
    phi_lin(T,J);

Variable  p_Lh(T,J),    p_L(T,J);
Positive Variable   Y1, Y2, Y3;
Variable Ytot, Ytot2;

p_es.up(T,J) = 12000000;
p_s.lo(T,J) = PinMin(J);
*c_start.up(T,J) = 30;
*c_end.up(T,J) = 30;
*c_in.up(T,J)   = 30;
c_in.fx(T,'1') = 0;


*c_inject.fx(T,'1') = 5;

Equations
    PressureBalance(T),
    PressureContinuity(T,J),
    PipelineExitPressure(T,J),
    StationExitToPipeline(T,J),
    StationPressureBalance(T,J),
    StationPressure(T,J),
    PumpHead(T,J,N),
    StationHead(T,J),
    SupplyRange(T),
    PipelineTotalLoss(T,J),
    PipelineFrictionLoss(T,J),
    PipelineElevationLoss(T,J);




PressureBalance(T)$(TactEq(T))..
    p_s(T,'1') + sum(J$(ord(J) < card(J)), MP(J)*p_p(T,J)) - sum(J$(ord(J) < card(J)), ML(J)*p_L(T,J)) =e= p_d(T,'5');

PressureContinuity(T,J)$(ord(J)>1 and TactEq(T))..
    p_s(T,J) =e= p_ed(T,J-1);

PipelineExitPressure(T,J)$(ML(J)  and TactEq(T)) ..
    P_ed(T,J) =e= P_es(T,J) - p_L(T,J);

StationExitToPipeline(T,J)$(ML(J)  and TactEq(T))..
    P_es(T,J) =e= p_d(T,J);

StationPressureBalance(T,J)$(TactEq(T))..
    p_d(T,J) =e= p_s(T,J) + MP(J)*p_p(T,J);

StationPressure(T,J)$(MP(J)  and TactEq(T))..
    p_p(T,J) =e= Rho_p(T,J) * g * Hj(T,J);

PumpHead(T,J,N)$(P(J,N) and TactEq(T))..
    h(T,J,N) =e= K(T,J,N) * (a1(N) - b1(N) * Qs(T,J)**2);

StationHead(T,J)$(MP(J) and TactEq(T))..
    Hj(T,J) =e= sum(N$(P(J,N)), h(T,J,N));

SupplyRange(T)$(TactEq(T)) ..
    p_s(T,'1') =e= 50*6894;

PipelineTotalLoss(T,J)$(TactEq(T))..
    p_L(T,J) =e= p_Lf(T,J) + p_Lh(T,J);

PipelineFrictionLoss(T,J)$(TactEq(T))..
    p_Lf(T,J) =e=
         draFactor(T,J) * 10.67 * (Qs(T,J)**1.852) / ( (130 ** 1.852) * (d(J) ** 4.8704) ) * sum(I, rhoBatch(I)*g
         * ( Dc(T,I,J+1)*(X(T,I) - Xj(J))  + Dd(T,I,J+1)*(X(T,I) - X(T,I+1))  + De(T,I,J+1)*(Xj(J+1) - Xj(J))  + Df(T,I,J+1)*(Xj(J+1) - X(T,I+1))) / (pi/4 * d(J)**2));

PipelineElevationLoss(T,J)$(TactEq(T))..
    p_Lh(T,J) =e= sum(I,   rhoBatch(I)*g * ( Dc(T,I,J+1)*(zBatch(T,I)-zStation(J))  + Dd(T,I,J+1)*(zBatch(T,I)-zBatch(T,I+1)) + De(T,I,J+1)*(zStation(J+1)-zStation(J))  + Df(T,I,J+1)*(zStation(J+1)-zBatch(T,I+1))));


Equations
    SwitchPosDetect(T,J,N),    SwitchNegDetect(T,J,N),
    SwitchPosUB(T,J,N),    SwitchNegUB(T,J,N),
    NoSimulSwitch(T,J,N),
    SwitchPosAtT0(J,N),    SwitchNegAtT0(J,N);
SwitchPosDetect(T,J,N)$(P(J,N) and TactEq(T) and ord(T) > tStartIdx)..    SwitchPos(T,J,N) =g= K(T,J,N) - K(T-1,J,N);
SwitchNegDetect(T,J,N)$(P(J,N) and TactEq(T) and ord(T) > tStartIdx)..    SwitchNeg(T,J,N) =g= K(T-1,J,N) - K(T,J,N);
SwitchPosUB(T,J,N)$(P(J,N) and TactEq(T) and ord(T) > tStartIdx)..    SwitchPos(T,J,N) =l= K(T,J,N);
SwitchNegUB(T,J,N)$(P(J,N) and TactEq(T) and ord(T) > tStartIdx)..    SwitchNeg(T,J,N) =l= 1 - K(T,J,N);
NoSimulSwitch(T,J,N)$(P(J,N) and TactEq(T) and ord(T) > tStartIdx)..   SwitchPos(T,J,N) + SwitchNeg(T,J,N) =l= 1;


Parameter K0(J,N);  K0(J,N) = 0;
*K0('1','1') = 1;
*K0('1','3') = 1;

SwitchPosAtT0(J,N)$(P(J,N))..    SwitchPos('t0',J,N) =g= K('t0',J,N) - K0(J,N);
SwitchNegAtT0(J,N)$(P(J,N))..    SwitchNeg('t0',J,N) =g= K0(J,N) - K('t0',J,N);

$ontext
Equations     MinUpTime(T,J,N),    MinDownTime(T,J,N),    MinUp_Prefix(T,J,N),    MinDown_Prefix(T,J,N);
MinUpTime(T,J,N)$(P(J,N) and ord(T) >= W  and (ord(T) < maxT))..    sum(Tp$((ord(Tp) > ord(T)-W) and (ord(Tp) <= ord(T))),    SwitchPos(Tp,J,N)) =l= K(T,J,N);
MinDownTime(T,J,N)$(P(J,N) and ord(T) >= W  and (ord(T) < maxT))..  sum(Tp$((ord(Tp) > ord(T)-W) and (ord(Tp) <= ord(T))),    SwitchNeg(Tp,J,N)) =l= 1 - K(T,J,N);
MinUp_Prefix(T,J,N)$(P(J,N) and ord(T) < W  and (ord(T) < maxT))..  sum(Tp$(ord(Tp) <= ord(T)), SwitchPos(Tp,J,N)) =l= K(T,J,N);
MinDown_Prefix(T,J,N)$(P(J,N) and ord(T) < W  and (ord(T) < maxT))..sum(Tp$(ord(Tp) <= ord(T)), SwitchNeg(Tp,J,N)) =l= 1 - K(T,J,N);
$offtext

Equation AnyPumpOn_lb, AnyPumpOn_ub;
AnyPumpOn_lb(T,J)$(MP(J) and TactEq(T)).. AnyPumpOn(T,J) =g= 1/card(N) * sum(N$P(J,N), K(T,J,N));
AnyPumpOn_ub(T,J)$(MP(J) and TactEq(T)).. AnyPumpOn(T,J) =l= sum(N$P(J,N), K(T,J,N));



Positive Variable M(T,J), cbar(T,J);
Parameter VsegVol(J), k_time(T,J);
k_time(T,J) = beta * (f(T,J)*Rho_p(T,J)*V_seg(T,J)**2/8)**m1 ;
VsegVol(J) = (Xj(J+1) - Xj(J)) * (pi/4 * d(J)**2);
Equations
    MassBal(T,J),
    CbarDef(T,J),
    CendDef(T,J);
    


MassBal(T,J)$(TactEq(T) and ord(J)<card(J))..
    M(T+1,J) =e= M(T,J)*exp(-k_time(T,J)*dt) + dt*(Qd(T,J)*c_out(T,J) - Qs(T,J+1)*c_in(T,J+1))*exp(-k_time(T,J)*dt);


CbarDef(T,J)$(TactVar(T) and ord(J)<card(J))..
    cbar(T,J) =e= M(T,J) / VsegVol(J);
    
Parameter shear(T,J);
shear(T,J) = ( f(T,J) * Rho_p(T,J) * V_seg(T,J)**2 ) / 8 ;

* Initialize concentrations from initial mass after VsegVol is defined (only once before windows)
if (card(Window),
    cbar0(J)      = carry_M(J) / VsegVol(J);
    carry_cin(J)  = 0;
    carry_cin(J)$(ord(J)>1) = cbar0(J-1);
    carry_cout(J)   = cbar0(J);
    carry_cstart(J) = cbar0(J);
    carry_cend(J)   = cbar0(J);
);

CendDef(T,J)$(ord(J)<card(J) and TactEq(T) and shiftActive(T,J))..
    c_end(T,J) =e=
        sum(Tp$shiftMap(T,Tp,J),
            exp( -beta * shear(Tp,J) * ( L(J) / max(V_seg(Tp,J), eps) ) )
          * c_start(Tp,J)
        );




SOS2 Variable lam(T,J,Im);
Equations
    sos2_sum(T,J)
    sos2_mlink(T,J)
    sos2_philink(T,J)
    mlin_cbar(T,J)
    draFactorDef(T,J);

sos2_sum(T,J)$(ML(J) and TactVar(T))..
    sum(Im, lam(T,J,Im)) =e= 1;

sos2_mlink(T,J)$(ML(J) and TactVar(T))..
    m_lin(T,J) =e= sum(Im, Mbp(Im)*lam(T,J,Im));

sos2_philink(T,J)$(ML(J) and TactVar(T))..
    phi_lin(T,J) =e= sum(Im, Phi_bp(Im)*lam(T,J,Im));

mlin_cbar(T,J)$(ML(J) and TactVar(T))..
    m_lin(T,J) =e= cbar(T,J);

draFactorDef(T,J)$(ML(J) and TactVar(T))..
    draFactor(T,J) =e= 0.21 + 0.79*phi_lin(T,J);









Equations
    DRA_Continuity(T,J)
    Yprod_ub1(T,J)
    Yprod_ub2(T,J)
    Yprod_lb(T,J);
* Outlet concentration balance:
DRA_Continuity(T,J)$(ord(J) < card(J) and TactEq(T) and not Tinit(T))..
    c_out(T,J) =e= c_in(T,J) - yprod(T,J) + c_inject(T,J);
* Big-M linearization for yprod = c_in * AnyPumpOn(t,j+1):
Yprod_ub1(T,J)$(ord(J) < card(J) and TactEq(T))..
    yprod(T,J) =l= cMax_dra * AnyPumpOn(T,J+1);

Yprod_ub2(T,J)$(ord(J) < card(J) and TactEq(T))..
    yprod(T,J) =l= c_in(T,J);

Yprod_lb(T,J)$(ord(J) < card(J) and TactEq(T))..
    yprod(T,J) =g= c_in(T,J) - cMax_dra * (1 - AnyPumpOn(T,J+1));
    

Equation CinFromPrev(T,J);
CinFromPrev(T,J)$(ord(J)>1 and ML(J) and TactEq(T))..
    c_in(T,J) =e= c_end(T,J-1);
Equation CstartDef(T,J);
CstartDef(T,J)$(ord(J)<card(J) and TactVar(T))..
    c_start(T,J) =e= c_out(T,J);


Equations
    DefY1
    DefY2
    DefY3
    DefYtot
    DefYtot2;
DefY1..
    Y1 =e= sum((T,J,N,Hr)$(TactVar(T) and MP(J) and P(J,N)),
               hourMap(T,Hr)*C1_hour(Hr)*dt/3.6e6*Rho_p(T,J)*g*h(T,J,N)*Qs(T,J)/eff(T,J,N));
DefY2..
    Y2 =e= C2*sum((T,J,N)$(P(J,N) and ord(T)>1 and TactVar(T)),
                  SwitchPos(T,J,N) + SwitchNeg(T,J,N));
DefY3..
    Y3 =e=
        sum((T,J)$(MA(J) and TactVar(T)),
            c_inject(T,J)*C3*Qs(T,J)*dt)
      + sum(J$(ord(J)<card(J)),
            M('t0',J)*C3);

DefYtot..
    Ytot =e= Y1 + Y3 + 0.01*Y2;
DefYtot2..
    Ytot2 =e= Y1;
*============================================================




Model ALlEquations / all /;

option mip   = gurobi;
option optcr = 0.1;

*------------------------------------------------------------
* Fixed-length sliding-horizon solve (keeps model size constant)
loop(Window,
    tStartIdx = smin(T$Twindow(Window,T), ord(T));
    tEndIdx   = smax(T$Twindow(Window,T), ord(T));

    TactVar(T) = no;
    TactVar(T)$(Twindow(Window,T)) = yes;

    TactEq(T)  = no;
    TactEq(T)$(Twindow(Window,T) and ord(T) < tEndIdx) = yes;

    Tinit(T)   = no;
    Tterm(T)   = no;
    Tinit(T)$(ord(T) = tStartIdx) = yes;
    Tterm(T)$(ord(T) = tEndIdx)   = yes;

*   * rebuild shift map only for the active window
    shiftMap(T,Tp,J) = 0;
    shiftMap(T,Tp,J)$(TactVar(T) and TactVar(Tp) and (ord(T) = ord(Tp) + delay(Tp,J))) = 1;
    shiftActive(T,J) = 0;
    shiftActive(T,J)$(TactVar(T)) = sum(Tp$(TactVar(Tp)), shiftMap(T,Tp,J));

*    * carry-in states at the window start
    M.fx(T,J)$Tinit(T)        = carry_M(J);
    c_in.fx(T,J)$Tinit(T)     = carry_cin(J);
    c_out.fx(T,J)$Tinit(T)    = carry_cout(J);
    c_start.fx(T,J)$Tinit(T)  = carry_cstart(J);
    c_end.fx(T,J)$Tinit(T)    = carry_cend(J);
    K.fx(T,J,N)$Tinit(T)      = carry_K(J,N);

*    * allow decisions inside window
    K.lo(T,J,N)$TactVar(T) = 0;
    K.up(T,J,N)$TactVar(T) = 1;
    
    c_inject.fx(T,J)$(TactVar(T) and MA(J) and HighPriceT(T)) = 6;

    c_inject.fx(T,J)$(TactVar(T) and MA(J) and not HighPriceT(T)) = 0;
*    c_inject.up(T,J)$(TactVar(T) and MA(J) and not HighPriceT(T)) = cMax_dra;

    solve ALlEquations using MIP minimizing Ytot;

*    * store terminal states for next window
    carry_M(J)      = sum(T$Tterm(T), M.l(T,J));
    carry_cin(J)    = sum(T$Tterm(T), c_in.l(T,J));
    carry_cout(J)   = sum(T$Tterm(T), c_out.l(T,J));
    carry_cstart(J) = sum(T$Tterm(T), c_start.l(T,J));
    carry_cend(J)   = sum(T$Tterm(T), c_end.l(T,J));
    carry_K(J,N)    = sum(T$Tterm(T), K.l(T,J,N));
    display carry_M, carry_cin, carry_cout, carry_cstart, carry_cend, carry_K;

*    * freeze solved portion as history (keeps feasibility for later windows)
    c_in.fx(T,J)     $(TactVar(T)) = c_in.l(T,J);
    c_out.fx(T,J)    $(TactVar(T)) = c_out.l(T,J);
    c_inject.fx(T,J) $(TactVar(T)) = c_inject.l(T,J);
    c_start.fx(T,J)  $(TactVar(T)) = c_start.l(T,J);
    c_end.fx(T,J)    $(TactVar(T)) = c_end.l(T,J);
    M.fx(T,J)        $(TactVar(T)) = M.l(T,J);
    p_s.fx(T,J)      $(TactVar(T)) = p_s.l(T,J);
    p_d.fx(T,J)      $(TactVar(T)) = p_d.l(T,J);
    K.fx(T,J,N)      $(TactVar(T) and P(J,N)) = K.l(T,J,N);
);


Set
    Day "day index" / d1*d44 /;

Parameter dayMap(T,Day) "Binary: T belongs to day D";

dayMap(T,Day)$(ord(Day) = floor(( (ord(T)-1)*dt ) / (3600*24)) + 1) = 1;


Parameter Y1_daily(Day) "Daily energy cost [CAD]";

Y1_daily(Day) = sum((T,J,N,Hr)$(P(J,N) and hourMap(T,Hr) and dayMap(T,Day)),
                     C1_hour(Hr) * dt/3.6e6 * Rho_p(T,J) * g * h.l(T,J,N) * Qs(T,J) / eff(T,J,N));


Parameter Y2_daily(Day) "Daily pump switching cost [CAD]";

Y2_daily(Day) = sum((T,J,N)$(P(J,N) and ord(T)>1 and dayMap(T,Day)),
                     C2 * (SwitchPos.l(T,J,N) + SwitchNeg.l(T,J,N)));
Display Y1_daily, Y2_daily;


Parameter Y1_hourly(Hr) "Hourly energy cost for day1 [CAD]";
Y1_hourly(Hr)
    =  sum((T,J,N)$(P(J,N)
                   and dayMap(T,'d1')
                   and hourMap(T,Hr)),
           C1_hour(Hr)
         * dt/3.6e6
         * Rho_p(T,J)
         * g
         * h.l(T,J,N)
         * Qs(T,J)
         / eff(T,J,N)
         );

Parameter Y2_hourly(Hr) "Hourly switch cost for day1 [CAD]";
Y2_hourly(Hr)
    =  sum((T,J,N)$(P(J,N)
                   and ord(T)>1
                   and dayMap(T,'d1')
                   and hourMap(T,Hr)),
           C2
         * ( SwitchPos.l(T,J,N)
           + SwitchNeg.l(T,J,N) )
         );

Display Y1_hourly, Y2_hourly;

Parameter Y3_daily(Day) "Daily DRA cost [CAD]";

Y3_daily(Day) =
    sum((T,J)$(dayMap(T,Day)),
        c_inject.l(T,J) * C3 * Qs(T,J) * dt );
Display Y3_daily;


Parameter DRA_mass(T) "DRA mass per time step [ppm*m^3]";

DRA_mass(T) = sum(J, c_inject.l(T,J) * Qs(T,J) * dt );
Display DRA_mass;

Parameter DRA_mass_station(T,J) "DRA mass per time step per station [ppm*m^3]";

DRA_mass_station(T,J) = c_inject.l(T,J) * Qs(T,J) * dt;

Display DRA_mass_station, Qs;

Scalar eff_min, eff_max;
eff_min = smin((T,J,N)$P(J,N), eff(T,J,N));
eff_max = smax((T,J,N)$P(J,N), eff(T,J,N));
Display eff_min, eff_max;

display Re,eff;
