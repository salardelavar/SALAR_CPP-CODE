#include <iostream> //std::cout
#include <conio.h>
#include <fstream>
//#include <cstdlib>
#include <cmath> //std::pow - sqrt
#include <ConsoleColor.h>
#include <iomanip> // std:: round number
#include <ctime>
using namespace std;
main(){
    double Ptarget,b,h,fc,ecu,fy,fu,As1,As2,As3,As4,As5,As6,d1,d2,d3,d4,d5,d6,AsTotal,SCmin,SCmax,SCinit,OFinit,f,SC,OF;
	double Es,Esh,x,beta1,tolerance,toleranceSC,a,ec,coff1,coff2,coff3,coff4,coff5,coff6,CHECK[50];
	double ey,esh,esu,fs1,fs2,fs3,fs4,fs5,fs6,fstan1,fstan2,fstan3,fstan4,fstan5,fstan6;
	double Fs1,Fs2,Fs3,Fs4,Fs5,Fs6,Fstan1,Fstan2,Fstan3,Fstan4,Fstan5,Fstan6;
	double e1,e2,e3,e4,e5,e6,es1,es2,es3,es4,es5,es6,residual,residualSC,A,A_tan,Cc,Cctan,FsSUM,FstanSUM,dx,dSC;
	double Ie,fr,Ec,ecr,xE,xY,xU,CurE,MomE,Finit,F,n,Iep,Iep1,Iep2,Iep21,Iep22,Iep23,Iep3,Iep31,Iep32,Iep33,Pc1,CurY,MomY,CurU,MomU;
	int i,it,itSC,itermax,Y,ZZ;	
ifstream INo;INo.open("ConcreteSectionDuctilityOptimizationCV-input.csv");
INo>>Ptarget>>b>>h>>SCmin>>SCmax>>fc>>ecu>>fy>>fu>>Es>>esh>>esu>>itermax>>tolerance>>toleranceSC>>coff1>>coff2>>coff3>>coff4>>coff5>>coff6>>d1>>d2>>d3>>d4>>d5>>d6;
INo.close();
cout<< green <<"     ************************************************************\n";    
cout<<"     *               >> IN THE NAME OF ALLAH <<                 *\n";
cout<<"     *          Concrete Section Ductility Optimization         *\n";
cout<<"     *                 with Continuous Variables                *\n";
cout<<"     *               UNIT: [Newton-Millimeter]                  *\n";
cout<<"     *----------------------------------------------------------*\n";
cout<<"     *      This program is written by Salar Delavar Qashqai    *\n";  
cout<<"     *          E-mail: salar.d.ghashghaei@gmail.com            *\n";
cout<<"     *----------------------------------------------------------*\n";
cout<<"     *      Notice: All input values must be positive           *\n";
cout<<"     *  _    ______________________________________             *\n";
cout<<"     *  |   |                                      |            *\n";
cout<<"     *      |     #     #     #     #    #    #    |            *\n";
cout<<"     *      |     #                           #    |            *\n";
cout<<"     *  b   |    As1   As2   As3   As4  As5  As6   |            *\n";
cout<<"     *      |     #                           #    |            *\n";
cout<<"     *  |   |     #     #     #     #    #    #    |            *\n";
cout<<"     *  _   |______________________________________|            *\n";
cout<<"     *      |<-                 h                ->|            *\n";
cout<<"     *      |<-d1->|                                            *\n";
cout<<"     *      |<-  d2   ->|                                       *\n";
cout<<"     *      |<-     d3      ->|                                 *\n";
cout<<"     *      |<-        d4          ->|                          *\n";
cout<<"     *      |<-            d5          ->|                      *\n";
cout<<"     *      |<-               d6            ->|                 *\n";
cout<<"     ************************************************************\n\n";
    cout<<"\n\n========================= Input Data ============================\n"<<endl;
    cout<<"      Axial load [+ : Compression]:                     "<<Ptarget<<endl;
	cout<<"      Concrete section width(b):                        "<<b<<endl;
    cout<<"      Concrete section height(h):                       "<<h<<endl;
    cout<<"      Minimum section ductility ratio:                  "<<SCmin<<endl;
    cout<<"      Maximum section ductility ratio:                  "<<SCmax<<endl;
    cout<<"      Concrete strength(fc):                            "<<fc<<endl;
    cout<<"      Ultimate concrete strain(ecu):                    "<<ecu<<endl;
    cout<<"      Yield rebar strength(fy):                         "<<fy<<endl;
    cout<<"      Ultimate rebar strength(fu):                      "<<fu<<endl;
    cout<<"      Rebar elastic modulus(Es):                        "<<Es<<endl;
    cout<<"      Strain at steel strain-hardening:                 "<<esh<<endl;
    cout<<"      Ultimate steel strain:                            "<<esu<<endl;
    cout<<"      Maximum number of iterations:                     "<<itermax<<endl;// maximum number of iterations
    cout<<"      Specified tolerance for convergence-1:            "<<tolerance<<endl;//  specified tolerance for convergence
    cout<<"      Specified tolerance for convergence-2:            "<<toleranceSC<<endl;//  specified tolerance for convergence
    ifstream IN2;IN2.open("ConcreteSectionDuctilityOptimizationCV-input.csv");
	Y=0;
	while(! IN2.eof())
	{
	IN2>>CHECK[Y];
	if(CHECK[Y] < 0)
	{
	cout<< red <<"\n  ======================= Error Report =========================="<<endl;	
	cout<<"\n            Please check your input file! -> [ ConcreteSectionDuctilityOptimizationCV-input.csv ]"<<endl;
	cout<<"        *** Negative input value is not acceptable ***"<<endl;
	Sleep(40000);
	exit(1);
	}
	Y++;
	}
	IN2.close();
	
    cout<< yellow <<"\n======================= Analysis Report =========================="<<endl;
    std::clock_t start;
    double duration;
    e1=d1-0.5*h;e2=d2-0.5*h;e3=d3-0.5*h;e4=d4-0.5*h;e5=d5-0.5*h;e6=d6-0.5*h;
    if (fc<= 30)
    {beta1=0.85;}
    else if (fc> 30 && fc< 55)
    {beta1=0.85-.008*(fc-30);}
    else  if (fc>= 55)
    {beta1=0.65;}
    ey=fy/Es; // Yield steel strain
    Esh=(fu-fy)/(esu-esh);
    x=.5*h; // initial guess of Neuteral axis
    AsTotal=0.1*b*h;
    As1=coff1*AsTotal;
    As2=coff2*AsTotal;
    As3=coff3*AsTotal;
    As4=coff4*AsTotal;
    As5=coff5*AsTotal;
    As6=coff6*AsTotal;
    // Crack capacity
    Ie=(b*pow(h,3))/12;fr=.7*sqrt(fc);Ec=5000*sqrt(fc);ecr=fr/Ec;xE=x;
    MomE=(fr*Ie/x)*pow(10,-6);// Elastic Moment
    CurE=((fr*Ie/x)/(Ec*Ie))*1000;// Elastic Curvature
    cout<<"\n     (STEP 1) : It is converged in "<<1<<" iterations - strain: "<<ecr<<" - x(mm): "<<CurE<<" - Phi(1/m): "<<CurE<<" - Moment(kN/m): "<<MomE<<endl;
        // Yield capacity
    n=Es/Ec;
    it = 0; // initialize iteration count
    residual = 100; /// initialize residual
    es1= (ey*(d1-x))/(d6-x);
    es2= (ey*(d2-x))/(d6-x);
    es3= (ey*(d3-x))/(d6-x);
    es4= (ey*(d4-x))/(d6-x);
    es5= (ey*(d5-x))/(d6-x);
    Finit=Es*(ey*As6+es5*As5+es4*As4+es3*As3+es2*As2+es1*As1)+0.5*fc*b*x+Ptarget;
    while ( residual > tolerance )
    {
    es1= (ey*(d1-x))/(d6-x);
    es2= (ey*(d2-x))/(d6-x);
    es3= (ey*(d3-x))/(d6-x);
    es4= (ey*(d4-x))/(d6-x);
    es5= (ey*(d5-x))/(d6-x);
    F=Es*(ey*As6+es5*As5+es4*As4+es3*As3+es2*As2+es1*As1)+0.5*fc*b*x+Ptarget;
    dx = (1/Finit)*(F);
    x = x+dx;
    residual = abs(dx); // evaluate residual
    it = it+1; // increment iteration count
    if (it == itermax) // stop the the analysis of this step please of Convergence
    {
    cout<<"\n     (STEP 2) : Iteration reached to Ultimate "<<itermax<<" - error: "<<dx<<endl;
    cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
    break;
    }
    }
    Iep1 = (b*pow(x,3))/12 +b*x*pow(x-.5*x,2);
    Iep21=(2*n-1)*As1*pow(x-d1,2);Iep22=(2*n-1)*As2*pow(x-d2,2);Iep23=(2*n-1)*As3*pow(x-d3,2);Iep2=Iep21+Iep22+Iep23;
    Iep31=n*As4*pow(x-d4,2);Iep32=n*As5*pow(x-d5,2);Iep33=n*As6*pow(x-d6,2);Iep3=Iep31+Iep32+Iep33;
    Iep=Iep1+Iep2+Iep3;xY=x;
    Pc1=x-.5*h;// Distance from applied force from neutral axis
    MomY = ((.5*fc*Iep/x)-Ptarget*Pc1)*pow(10,-6);// Yield Moment
    CurY =((.5*fc*Iep/x)/(Ec*Iep))*1000;// Yield Curvature
        if (it < itermax)// iteration control
        {
        cout<<"\n     (STEP 2) : It is converged in "<<it<<" iterations - strain: "<<ey<<" - x(mm): "<<xY<<" - Phi(1/m): "<<CurY<<" - Moment(kN/m): "<<MomY<<endl;
        } 
    
    // Ultimate capacity
    it = 0; // initialize iteration count
    residual = 100; // initialize residual
    while (residual > tolerance)
    {
    a=beta1*x;Cc=0.85*fc*a*b;Cctan=0;ec=0.5*(h-a);
   
	es1=ecu*(d1-x)/x;
	es2=ecu*(d2-x)/x;
	es3=ecu*(d3-x)/x;
	es4=ecu*(d4-x)/x;
	es5=ecu*(d5-x)/x;
	es6=ecu*(d6-x)/x;
	// As1
	if  (es1>=0 && es1<=ey)
        {fs1=Es*es1;fstan1=(Es*ecu*d1/pow(x,2));}
    else if  (es1<0 && es1>=-ey)
        {fs1=Es*es1;fstan1=(Es*ecu*d1/pow(x,2));}
    else if (es1>ey && es1<=esh)
        {fs1=fy;fstan1=0;}
    else if (es1<-ey && es1>=-esh)
        {fs1=-fy;fstan1=0;}
    else if (es1>esh && es1<=esu)
        {fs1=fy+Esh*(abs(es1)-esh);fstan1=(Esh*ecu*d1/pow(x,2));}
    else if (es1<-esh && es1>=-esu)
        {fs1=-fy-Esh*(abs(es1)-esh);fstan1=(Esh*ecu*d1/pow(x,2));}   
    else if (es1<-esu || es1>esu)
        {fs1=0;fstan1=0;}
       
    if (d1>a)
        {Fs1=As1*fs1;Fstan1=As1*fstan1;}
    else if (d1<=a)
        {Fs1=As1*(fs1-0.85*fc);Fstan1=As1*(fstan1-0.85*fc);}
    // As2
	if  (es2>=0 && es2<=ey)
        {fs2=Es*es2;fstan2=(Es*ecu*d2/pow(x,2));}
    else if  (es2<0 && es2>=-ey)
        {fs2=Es*es2;fstan2=(Es*ecu*d2/pow(x,2));}
    else if (es2>ey && es2<=esh)
        {fs2=fy;fstan2=0;}
    else if (es2<-ey && es2>=-esh)
        {fs2=-fy;fstan2=0;}
    else if (es2>esh && es2<=esu)
        {fs2=fy+Esh*(abs(es2)-esh);fstan2=(Esh*ecu*d2/pow(x,2));}
    else if (es2<-esh && es2>=-esu)
        {fs2=-fy-Esh*(abs(es2)-esh);fstan2=(Esh*ecu*d2/pow(x,2));}   
    else if (es2<-esu || es2>esu)
        {fs2=0;fstan2=0;}
        
    if (d2>a)
        {Fs2=As2*fs2;Fstan2=As2*fstan2;}
    else if (d2<=a)
        {Fs2=As2*(fs2-0.85*fc);Fstan2=As2*(fstan2-0.85*fc);}
	// As3
	if  (es3>=0 && es3<=ey)
        {fs3=Es*es3;fstan3=(Es*ecu*d3/pow(x,2));}
    else if (es3<0 && es3>=-ey)
        {fs3=Es*es3;fstan3=(Es*ecu*d3/pow(x,2));}
    else if (es3>ey && es3<=esh)
        {fs3=fy;fstan3=0;}
    else if (es3<-ey && es3>=-esh)
        {fs3=-fy;fstan3=0;}
    else if (es3>esh && es3<=esu)
        {fs3=fy+Esh*(abs(es3)-esh);fstan3=(Esh*ecu*d3/pow(x,2));}
    else if (es3<-esh && es3>=-esu)
        {fs3=-fy-Esh*(abs(es3)-esh);fstan3=(Esh*ecu*d3/pow(x,2));}   
    else if (es3<-esu || es3>esu)
        {fs3=0;fstan3=0;}
        
    if (d3>a)
        {Fs3=As3*fs3;Fstan3=As3*fstan3;}
    else if (d3<=a)
        {Fs3=As3*(fs3-0.85*fc);Fstan3=As3*(fstan3-0.85*fc);}
	// As4
	if  (es4>=0 && es4<=ey)
        {fs4=Es*es4;fstan4=(Es*ecu*d4/pow(x,2));}
    else if  (es4<0 && es4>=-ey)
        {fs4=Es*es4;fstan4=(Es*ecu*d4/pow(x,2));}
    else if (es4>ey && es4<=esh)
        {fs4=fy;fstan4=0;}
    else if (es4<-ey && es4>=-esh)
        {fs4=-fy;fstan4=0;}
    else if (es4>esh && es4<=esu)
        {fs4=fy+Esh*(abs(es4)-esh);fstan4=(Esh*ecu*d4/pow(x,2));}
    else if (es4<-esh && es4>=-esu)
        {fs4=-fy-Esh*(abs(es4)-esh);fstan4=(Esh*ecu*d4/pow(x,2));}   
    else if (es4<-esu || es4>esu)
        {fs4=0;fstan4=0;}
        
    if (d4>a)
        {Fs4=As4*fs4;Fstan4=As4*fstan4;}
    else if (d4<=a)
        {Fs4=As4*(fs4-0.85*fc);Fstan4=As4*(fstan4-0.85*fc);}
        
	// As5
	if  (es5>=0 && es5<=ey)
        {fs5=Es*es5;fstan5=(Es*ecu*d5/pow(x,2));}
    else if  (es5<0 && es5>=-ey)
        {fs5=Es*es5;fstan5=(Es*ecu*d5/pow(x,2));}
    else if (es5>ey && es5<=esh)
        {fs5=fy;fstan5=0;}
    else if (es5<-ey && es5>=-esh)
        {fs5=-fy;fstan5=0;}
    else if (es5>esh && es5<=esu)
        {fs5=fy+Esh*(abs(es5)-esh);fstan5=(Esh*ecu*d5/pow(x,2));}
    else if (es5<-esh && es5>=-esu)
        {fs5=-fy-Esh*(abs(es5)-esh);fstan5=(Esh*ecu*d5/pow(x,2));}   
    else if (es5<-esu || es5>esu)
        {fs5=0;fstan5=0;}
        
    if (d5>a)
        {Fs5=As5*fs5;Fstan5=As5*fstan5;}
    else if (d5<=a)
        {Fs5=As5*(fs5-0.85*fc);Fstan5=As5*(fstan5-0.85*fc);}
	
	// As6
	if  (es6>=0 && es6<=ey)
        {fs6=Es*es6;fstan6=(Es*ecu*d6/pow(x,2));}
    else if  (es6<0 && es6>=-ey)
        {fs6=Es*es6;fstan6=(Es*ecu*d6/pow(x,2));}
    else if (es6>ey && es6<=esh)
        {fs6=fy;fstan6=0;}
    else if (es6<-ey && es6>=-esh)
        {fs6=-fy;fstan6=0;}
    else if (es6>esh && es6<=esu)
        {fs6=fy+Esh*(abs(es6)-esh);fstan6=(Esh*ecu*d6/pow(x,2));}
    else if (es6<-esh && es6>=-esu)
        {fs6=-fy-Esh*(abs(es6)-esh);fstan6=(Esh*ecu*d6/pow(x,2));}   
    else if (es6<-esu || es6>esu)
        {fs6=0;fstan6=0;}
        
    if (d6>a)
        {Fs6=As6*fs6;Fstan6=As6*fstan6;}
    else if (d6<=a)
        {Fs6=As6*(fs6-0.85*fc);Fstan6=As6*(fstan6-0.85*fc);}
			    
				    
    FsSUM=Fs1+Fs2+Fs3+Fs4+Fs5+Fs6;
    FstanSUM=Fstan1+Fstan2+Fstan3+Fstan4+Fstan5+Fstan6;
    A=Cc-FsSUM-Ptarget;
	A_tan=Cctan-FstanSUM;
    dx = (1/A_tan)*(-A);
    residual = abs(dx); // evaluate residual
    it = it + 1; // increment iteration count
    x = x+dx; // update x
        if (it == itermax) // stop the the analysis of this step please of Convergence
          {
		  cout<<"\n     (-) Iteration reached to Ultimate "<<it<<" - strain: "<<ecu<<" - error: "<<abs(dx)<<endl;
          cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
          exit(2);
          getch();
          }
    
	}
	Pc1=x-.5*h;xU=x;
	MomU=(Fs1*e1+Fs2*e2+Fs3*e3+Fs4*e4+Fs5*e5+Fs6*e6+Cc*ec+Pc1*Ptarget)*pow(10,-6);
	CurU=(ecu/x)*1000;
	if (it < itermax)// iteration control
    cout<<"\n     (STEP 3) : It is converged in "<<it<<" iterations - strain: "<<ecu<<" - x(mm): "<<x<<" - Phi(1/m): "<<CurU<<" - Moment (kN.m): "<<MomU<<"\n"<<endl;
    SCinit=CurU/CurY;OFinit=MomU/MomY;

    	cout<<"   ------------------- (Initial Values) -------------------"<<endl;
		cout<<"        As1      As2      As3      As4      As5      As6   "<<endl;
		cout<<"   --------------------------------------------------------"<<endl;
		cout<<"      "<<As1<<"     "<<As2<<"     "<<As3<<"     "<<As4<<"     "<<As5<<"     "<<As6<<endl;
		cout<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"       Yield Curvature   Ultimate Curvature   Section Ductility Ratio  Yield Moment   Ultimate Moment   Over Strength Factor  "<<endl;
        cout<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"\t  "<<setprecision(5)<<fixed<<CurY<<"\t     "<<setprecision(5)<<fixed<<CurU<<"\t\t\t"<<setprecision(2)<<fixed<<SCinit<<"\t\t"<<setprecision(3)<<fixed<<MomY<<"\t  "<<setprecision(3)<<fixed<<MomU<<"\t\t"<<setprecision(2)<<fixed<<OFinit<<endl;
    itSC=0;
    residualSC=100;
    while (residualSC > toleranceSC)
    {
            // Yield capacity
    n=Es/Ec;
    it = 0; // initialize iteration count
    residual = 100; /// initialize residual
    es1= (ey*(d1-x))/(d6-x);
    es2= (ey*(d2-x))/(d6-x);
    es3= (ey*(d3-x))/(d6-x);
    es4= (ey*(d4-x))/(d6-x);
    es5= (ey*(d5-x))/(d6-x);
    Finit=Es*(ey*As6+es5*As5+es4*As4+es3*As3+es2*As2+es1*As1)+0.5*fc*b*x+Ptarget;
    while ( residual > tolerance )
    {
    es1= (ey*(d1-x))/(d6-x);
    es2= (ey*(d2-x))/(d6-x);
    es3= (ey*(d3-x))/(d6-x);
    es4= (ey*(d4-x))/(d6-x);
    es5= (ey*(d5-x))/(d6-x);
    F=Es*(ey*As6+es5*As5+es4*As4+es3*As3+es2*As2+es1*As1)+0.5*fc*b*x+Ptarget;
    dx = (1/Finit)*(F);
    x = x+dx;
    residual = abs(dx); // evaluate residual
    it = it+1; // increment iteration count
    if (it == itermax) // stop the the analysis of this step please of Convergence
    {
    //cout<<"\n     (STEP 2) : Iteration reached to Ultimate "<<itermax<<" - error: "<<dx<<endl;
    //cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
    break;
    }
    }
    Iep1 = (b*pow(x,3))/12 +b*x*pow(x-.5*x,2);
    Iep21=(2*n-1)*As1*pow(x-d1,2);Iep22=(2*n-1)*As2*pow(x-d2,2);Iep23=(2*n-1)*As3*pow(x-d3,2);Iep2=Iep21+Iep22+Iep23;
    Iep31=n*As4*pow(x-d4,2);Iep32=n*As5*pow(x-d5,2);Iep33=n*As6*pow(x-d6,2);Iep3=Iep31+Iep32+Iep33;
    Iep=Iep1+Iep2+Iep3;xY=x;
    Pc1=x-.5*h;// Distance from applied force from neutral axis
    MomY = ((.5*fc*Iep/x)-Ptarget*Pc1)*pow(10,-6);// Yield Moment
    CurY =((.5*fc*Iep/x)/(Ec*Iep))*1000;// Yield Curvature
        if (it < itermax)// iteration control
        {
        //cout<<"\n     (STEP 2) : It is converged in "<<it<<" iterations - strain: "<<ey<<" - x(mm): "<<xY<<" - Phi(1/m): "<<CurY<<" - Moment(kN/m): "<<MomY<<endl;
        } 
    
    // Ultimate capacity
    it = 0; // initialize iteration count
    residual = 100; // initialize residual
    while (residual > tolerance)
    {
    a=beta1*x;Cc=0.85*fc*a*b;Cctan=0;ec=0.5*(h-a);
   
	es1=ecu*(d1-x)/x;
	es2=ecu*(d2-x)/x;
	es3=ecu*(d3-x)/x;
	es4=ecu*(d4-x)/x;
	es5=ecu*(d5-x)/x;
	es6=ecu*(d6-x)/x;
	// As1
	if  (es1>=0 && es1<=ey)
        {fs1=Es*es1;fstan1=(Es*ecu*d1/pow(x,2));}
    else if  (es1<0 && es1>=-ey)
        {fs1=Es*es1;fstan1=(Es*ecu*d1/pow(x,2));}
    else if (es1>ey && es1<=esh)
        {fs1=fy;fstan1=0;}
    else if (es1<-ey && es1>=-esh)
        {fs1=-fy;fstan1=0;}
    else if (es1>esh && es1<=esu)
        {fs1=fy+Esh*(abs(es1)-esh);fstan1=(Esh*ecu*d1/pow(x,2));}
    else if (es1<-esh && es1>=-esu)
        {fs1=-fy-Esh*(abs(es1)-esh);fstan1=(Esh*ecu*d1/pow(x,2));}   
    else if (es1<-esu || es1>esu)
        {fs1=0;fstan1=0;}
       
    if (d1>a)
        {Fs1=As1*fs1;Fstan1=As1*fstan1;}
    else if (d1<=a)
        {Fs1=As1*(fs1-0.85*fc);Fstan1=As1*(fstan1-0.85*fc);}
    // As2
	if  (es2>=0 && es2<=ey)
        {fs2=Es*es2;fstan2=(Es*ecu*d2/pow(x,2));}
    else if  (es2<0 && es2>=-ey)
        {fs2=Es*es2;fstan2=(Es*ecu*d2/pow(x,2));}
    else if (es2>ey && es2<=esh)
        {fs2=fy;fstan2=0;}
    else if (es2<-ey && es2>=-esh)
        {fs2=-fy;fstan2=0;}
    else if (es2>esh && es2<=esu)
        {fs2=fy+Esh*(abs(es2)-esh);fstan2=(Esh*ecu*d2/pow(x,2));}
    else if (es2<-esh && es2>=-esu)
        {fs2=-fy-Esh*(abs(es2)-esh);fstan2=(Esh*ecu*d2/pow(x,2));}   
    else if (es2<-esu || es2>esu)
        {fs2=0;fstan2=0;}
        
    if (d2>a)
        {Fs2=As2*fs2;Fstan2=As2*fstan2;}
    else if (d2<=a)
        {Fs2=As2*(fs2-0.85*fc);Fstan2=As2*(fstan2-0.85*fc);}
	// As3
	if  (es3>=0 && es3<=ey)
        {fs3=Es*es3;fstan3=(Es*ecu*d3/pow(x,2));}
    else if (es3<0 && es3>=-ey)
        {fs3=Es*es3;fstan3=(Es*ecu*d3/pow(x,2));}
    else if (es3>ey && es3<=esh)
        {fs3=fy;fstan3=0;}
    else if (es3<-ey && es3>=-esh)
        {fs3=-fy;fstan3=0;}
    else if (es3>esh && es3<=esu)
        {fs3=fy+Esh*(abs(es3)-esh);fstan3=(Esh*ecu*d3/pow(x,2));}
    else if (es3<-esh && es3>=-esu)
        {fs3=-fy-Esh*(abs(es3)-esh);fstan3=(Esh*ecu*d3/pow(x,2));}   
    else if (es3<-esu || es3>esu)
        {fs3=0;fstan3=0;}
        
    if (d3>a)
        {Fs3=As3*fs3;Fstan3=As3*fstan3;}
    else if (d3<=a)
        {Fs3=As3*(fs3-0.85*fc);Fstan3=As3*(fstan3-0.85*fc);}
	// As4
	if  (es4>=0 && es4<=ey)
        {fs4=Es*es4;fstan4=(Es*ecu*d4/pow(x,2));}
    else if  (es4<0 && es4>=-ey)
        {fs4=Es*es4;fstan4=(Es*ecu*d4/pow(x,2));}
    else if (es4>ey && es4<=esh)
        {fs4=fy;fstan4=0;}
    else if (es4<-ey && es4>=-esh)
        {fs4=-fy;fstan4=0;}
    else if (es4>esh && es4<=esu)
        {fs4=fy+Esh*(abs(es4)-esh);fstan4=(Esh*ecu*d4/pow(x,2));}
    else if (es4<-esh && es4>=-esu)
        {fs4=-fy-Esh*(abs(es4)-esh);fstan4=(Esh*ecu*d4/pow(x,2));}   
    else if (es4<-esu || es4>esu)
        {fs4=0;fstan4=0;}
        
    if (d4>a)
        {Fs4=As4*fs4;Fstan4=As4*fstan4;}
    else if (d4<=a)
        {Fs4=As4*(fs4-0.85*fc);Fstan4=As4*(fstan4-0.85*fc);}
        
	// As5
	if  (es5>=0 && es5<=ey)
        {fs5=Es*es5;fstan5=(Es*ecu*d5/pow(x,2));}
    else if  (es5<0 && es5>=-ey)
        {fs5=Es*es5;fstan5=(Es*ecu*d5/pow(x,2));}
    else if (es5>ey && es5<=esh)
        {fs5=fy;fstan5=0;}
    else if (es5<-ey && es5>=-esh)
        {fs5=-fy;fstan5=0;}
    else if (es5>esh && es5<=esu)
        {fs5=fy+Esh*(abs(es5)-esh);fstan5=(Esh*ecu*d5/pow(x,2));}
    else if (es5<-esh && es5>=-esu)
        {fs5=-fy-Esh*(abs(es5)-esh);fstan5=(Esh*ecu*d5/pow(x,2));}   
    else if (es5<-esu || es5>esu)
        {fs5=0;fstan5=0;}
        
    if (d5>a)
        {Fs5=As5*fs5;Fstan5=As5*fstan5;}
    else if (d5<=a)
        {Fs5=As5*(fs5-0.85*fc);Fstan5=As5*(fstan5-0.85*fc);}
	
	// As6
	if  (es6>=0 && es6<=ey)
        {fs6=Es*es6;fstan6=(Es*ecu*d6/pow(x,2));}
    else if  (es6<0 && es6>=-ey)
        {fs6=Es*es6;fstan6=(Es*ecu*d6/pow(x,2));}
    else if (es6>ey && es6<=esh)
        {fs6=fy;fstan6=0;}
    else if (es6<-ey && es6>=-esh)
        {fs6=-fy;fstan6=0;}
    else if (es6>esh && es6<=esu)
        {fs6=fy+Esh*(abs(es6)-esh);fstan6=(Esh*ecu*d6/pow(x,2));}
    else if (es6<-esh && es6>=-esu)
        {fs6=-fy-Esh*(abs(es6)-esh);fstan6=(Esh*ecu*d6/pow(x,2));}   
    else if (es6<-esu || es6>esu)
        {fs6=0;fstan6=0;}
        
    if (d6>a)
        {Fs6=As6*fs6;Fstan6=As6*fstan6;}
    else if (d6<=a)
        {Fs6=As6*(fs6-0.85*fc);Fstan6=As6*(fstan6-0.85*fc);}
			    
				    
    FsSUM=Fs1+Fs2+Fs3+Fs4+Fs5+Fs6;
    FstanSUM=Fstan1+Fstan2+Fstan3+Fstan4+Fstan5+Fstan6;
    A=Cc-FsSUM-Ptarget;
	A_tan=Cctan-FstanSUM;
    dx = (1/A_tan)*(-A);
    residual = abs(dx); // evaluate residual
    it = it + 1; // increment iteration count
    x = x+dx; // update x
        if (it == itermax) // stop the the analysis of this step please of Convergence
          {
		  //cout<<"\n     (-) Iteration reached to Ultimate "<<it<<" - strain: "<<ecu<<" - error: "<<abs(dx)<<endl;
          //cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
          break;
          getch();
          }
    
	}
	Pc1=x-.5*h;xU=x;
	MomU=(Fs1*e1+Fs2*e2+Fs3*e3+Fs4*e4+Fs5*e5+Fs6*e6+Cc*ec+Pc1*Ptarget)*pow(10,-6);
	CurU=(ecu/x)*1000;
	if (it < itermax)// iteration control
    //cout<<"\n     (STEP 3) : It is converged in "<<it<<" iterations - strain: "<<ecu<<" - x(mm): "<<x<<" - Phi(1/m): "<<CurU<<" - Moment (kN.m): "<<MomU<<"\n"<<endl;
    SC=CurU/CurY;OF=MomU/MomY;
    f=SC-SCmax;
    dSC=(1/SCinit) * f;
    residualSC = abs(dSC); // evaluate residual
    AsTotal = AsTotal+dSC; // update As
    As1=coff1*AsTotal;
    As2=coff2*AsTotal;
    As3=coff3*AsTotal;
    As4=coff4*AsTotal;
    As5=coff5*AsTotal;
    As6=coff6*AsTotal;
        itSC = itSC + 1; // increment iteration count
         if (SC >= SCmin && SC <= SCmax) // stop the the analysis of this step please of Convergence
         {ZZ=1;
        cout<<blue<<"\n  ### The Optimum section ductility ratio for these constrains: "<<setprecision(2)<<fixed<<SC<<" with "<<itSC<<" iterations ###"<<endl<<endl;   
		cout<<"   ------------------- (Optimum Values) -------------------"<<endl;
		cout<<"        As1      As2      As3      As4      As5      As6  "<<endl;
		cout<<"   --------------------------------------------------------"<<endl;
		cout<<"      "<<As1<<"   "<<As2<<"   "<<As3<<"   "<<As4<<"   "<<As5<<"   "<<As6<<endl;
		cout<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"       Yield Curvature   Ultimate Curvature   Section Ductility Ratio  Yield Moment   Ultimate Moment   Over Strength Factor  "<<endl;
        cout<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"\t  "<<setprecision(5)<<fixed<<CurY<<"\t     "<<setprecision(5)<<fixed<<CurU<<"\t\t\t"<<setprecision(2)<<fixed<<SCinit<<"\t\t"<<setprecision(3)<<fixed<<MomY<<"\t  "<<setprecision(3)<<fixed<<MomU<<"\t\t"<<setprecision(2)<<fixed<<OFinit<<endl;
          break;    
         }
         else
         ZZ=0;
        if (itSC == itermax) // stop the the analysis of this step of Convergence
        {
          cout<<red<<" +++ Could not find the solution - Trail iteration reached to Ultimate "<<itSC<<" - error: "<<dSC<<" +++"<<endl;
          break;
		}

	}
    
        if (ZZ == 0)// control
         cout<<red<<" +++ Please check your data +++\n"<<endl;
    
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		std::cout<< white <<"\n   CPU time: "<< duration <<" seconds"<<'\n';
    
	
    // current date/time based on current system
   time_t now = time(0);
   // convert now to string form
   char* dt = ctime(&now);

    cout<<"   Date and Time : " << dt;
    // OUTPUT
	ofstream OUTo;OUTo.open("ConcreteSectionDuctilityOptimizationNEWTON-output.txt");
	OUTo<<"   Date and Time : " << dt;
	OUTo<<"************************************************************\n";    
    OUTo<<"*                >> IN THE NAME OF ALLAH <<                *\n";
    OUTo<<"*          Concrete Section Ductility Optimization         *\n";
    OUTo<<"*                 with Continuous Variables                *\n";
    OUTo<<"*               UNIT: [Newton-Millimeter]                  *\n";
    OUTo<<"*----------------------------------------------------------*\n";
    OUTo<<"*     This program is written by Salar Delavar Qashqai     *\n";  
    OUTo<<"*          E-mail: salar.d.ghashghaei@gmail.com            *\n";
    OUTo<<"*----------------------------------------------------------*\n";
    OUTo<<"*      Notice: All input values must be positive           *\n";
    OUTo<<"*  _    ______________________________________             *\n";
    OUTo<<"*  |   |                                      |            *\n";
    OUTo<<"*      |     #     #     #     #    #    #    |            *\n";
    OUTo<<"*      |     #                           #    |            *\n";
    OUTo<<"*  b   |    As1   As2   As3   As4  As5  As6   |            *\n";
    OUTo<<"*      |     #                           #    |            *\n";
    OUTo<<"*  |   |     #     #     #     #    #    #    |            *\n";
    OUTo<<"*  _   |______________________________________|            *\n";
    OUTo<<"*      |<-                 h                ->|            *\n";
    OUTo<<"*      |<-d1->|                                            *\n";
    OUTo<<"*      |<-  d2   ->|                                       *\n";
    OUTo<<"*      |<-     d3      ->|                                 *\n";
    OUTo<<"*      |<-        d4          ->|                          *\n";
    OUTo<<"*      |<-            d5          ->|                      *\n";
    OUTo<<"*      |<-               d6            ->|                 *\n";
    OUTo<<"************************************************************\n\n";
    OUTo<<"  ### The Optimum section ductility ratio for these constrains: "<<setprecision(2)<<fixed<<SC<<" with "<<itSC<<" iterations ###"<<endl<<endl;  
	OUTo<<"   ------------------- (Optimum Values) -------------------"<<endl;
	OUTo<<"        As1      As2      As3      As4      As5      As6  "<<endl;
	OUTo<<"   --------------------------------------------------------"<<endl;
	OUTo<<"      "<<As1<<"   "<<As2<<"   "<<As3<<"   "<<As4<<"   "<<As5<<"   "<<As6<<endl;
    OUTo<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
    OUTo<<"       Yield Curvature   Ultimate Curvature   Section Ductility Ratio  Yield Moment   Ultimate Moment   Over Strength Factor  "<<endl;
    OUTo<<"   ---------------------------------------------------------------------------------------------------------------------------"<<endl;
    OUTo<<"\t  "<<setprecision(5)<<fixed<<CurY<<"\t     "<<setprecision(5)<<fixed<<CurU<<"\t\t\t"<<setprecision(2)<<fixed<<SCinit<<"\t\t"<<setprecision(3)<<fixed<<MomY<<"\t  "<<setprecision(3)<<fixed<<MomU<<"\t\t"<<setprecision(2)<<fixed<<OFinit<<endl;
    cout<<"\n  - Output data is written in file -"<<endl;
    OUTo.close();
    getch();
	}
