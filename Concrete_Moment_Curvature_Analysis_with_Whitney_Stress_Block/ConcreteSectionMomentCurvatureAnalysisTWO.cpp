#include <iostream> //std::cout
#include <conio.h>
#include <fstream>
#include <cmath> //std::pow - sqrt
#include <ConsoleColor.h>
#define STEP 4
#define ShowText01 "Graph-outputHTML.html"
using namespace std;
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
double MAX_ABS(double A[],int n);
double ABS(double B);	
main(){
//system("color 8E");	
//system("color 1F");
//system("color 0A");
	double Ptarget,b,h,fc,ecu,fy,fu;
	double Es,Esh,x,beta1,tolerance,a,ec,Finitsum,Fsum,SC,OF;
	double ey,esh,esu,residual,A,A_tan,Cc,Cctan,FsSUM,FstanSUM,dx;
	double Ie,fr,Ec,ecr,xE,xY,xU,CurE,MomE,Finit,F,n,Iep,Iep1,Iep2,Iep3,Pc1,CurY,MomY,CurU,MomU;
	double As[6],d[6],fs[6],fstan[6],Fs[6],Fstan[6],e[6],es[6];
	int j,it,itermax;	
ifstream INo;INo.open("ConcreteSectionMomentCurvatureAnalysis-input.txt");
INo>>Ptarget>>b>>h>>fc>>ecu>>fy>>fu>>Es>>esh>>esu>>itermax>>tolerance>>As[1]>>As[2]>>As[3]>>As[4]>>As[5]>>As[6]>>d[1]>>d[2]>>d[3]>>d[4]>>d[5]>>d[6];
INo.close();
cout<< green <<"     ************************************************************\n";    
cout<<"     *               >> IN THE NAME OF ALLAH <<                 *\n";
cout<<"     *       Concrete Section Moment Curvature Analysis         *\n";
cout<<"     *               UNIT: [Newton-Millimeter]                  *\n";
cout<<"     *----------------------------------------------------------*\n";
cout<<"     *     This program is written by Salar Delavar Qashqai     *\n";  
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
    cout<<"\n\n================== Input Data =======================\n"<<endl;
    cout<<"      Axial load [+ : Compression]: "<<Ptarget<<endl;
	cout<<"      b: "<<b<<endl;
    cout<<"      h: "<<h<<endl;
    cout<<"      Concrete strength(fc): "<<fc<<endl;
    cout<<"      Ultimate concrete strain(ecu): "<<ecu<<endl;
    cout<<"      Yield rebar strength(fy): "<<fy<<endl;
    cout<<"      Ultimate rebar strength(fu): "<<fu<<endl;
    cout<<"      Rebar elastic modulus(Es): "<<Es<<endl;
    cout<<"      Strain at steel strain-hardening: "<<esh<<endl;
    cout<<"      Ultimate steel strain: "<<esu<<endl;
    cout<<"      Area of rebar - As1(mm^2): "<<As[1]<<endl;
    cout<<"      Area of rebar - As2(mm^2): "<<As[2]<<endl;
    cout<<"      Area of rebar - As3(mm^2): "<<As[3]<<endl;
    cout<<"      Area of rebar - As4(mm^2): "<<As[4]<<endl;
    cout<<"      Area of rebar - As5(mm^2): "<<As[5]<<endl;
    cout<<"      Area of rebar - As6(mm^2): "<<As[6]<<endl;
    cout<<"      Distance of rebar - d1(mm): "<<d[1]<<endl;
    cout<<"      Distance of rebar - d2(mm): "<<d[2]<<endl;
    cout<<"      Distance of rebar - d3(mm): "<<d[3]<<endl;
    cout<<"      Distance of rebar - d4(mm): "<<d[4]<<endl;
    cout<<"      Distance of rebar - d5(mm): "<<d[5]<<endl;
    cout<<"      Distance of rebar - d6(mm): "<<d[6]<<endl;
    cout<<"      Maximum number of iterations: "<<itermax<<endl;// maximum number of iterations
    cout<<"      Specified tolerance for convergence: "<<tolerance<<endl;//  specified tolerance for convergence
    cout<< yellow <<"\n======================= Analysis Report =========================="<<endl;
    if (fc<= 30)
    {beta1=0.85;}
    else if (fc> 30 && fc< 55)
    {beta1=0.85-.008*(fc-30);}
    else  if (fc>= 55)
    {beta1=0.65;}
    ey=fy/Es; // Yield steel strain
    Esh=(fu-fy)/(esu-esh);
    x=.5*h; // initial guess of Neuteral axis
    // Crack capacity
    Ie=(b*pow(h,3))/12;fr=.7*sqrt(fc);Ec=5000*sqrt(fc);ecr=fr/Ec;xE=x;
    MomE=(fr*Ie/x)*pow(10,-6);// Elastic Moment
    CurE=((fr*Ie/x)/(Ec*Ie))*1000;// Elastic Curvature
    cout<<"\n     (STEP 1) : It is converged in "<<1<<" iterations - strain: "<<ecr<<" - x(mm): "<<CurE<<" - Phi(1/m): "<<CurE<<" - Moment(kN/m): "<<MomE<<endl;
    // Yield capacity
    n=Es/Ec;
    Finitsum=0;
    it = 0; // initialize iteration count
    residual = 100; /// initialize residual
    for (j=1;j<=5;j++)
    {
    es[j]= (ey*(d[j]-x))/(d[6]-x);
    Finitsum=Finitsum+Es*es[j]*As[j];
    }
    Finit=Finitsum+Es*ey*As[6]+0.5*fc*b*x+Ptarget;
    while ( residual > tolerance )
    {Fsum=0;
    for (j=1;j<=5;j++)	
     {
    es[j]= (ey*(d[j]-x))/(d[6]-x);
    Fsum=Fsum+Es*es[j]*As[j];
     }
    F=Fsum+Es*ey*As[6]+0.5*fc*b*x+Ptarget;
    dx = (1/Finit)*(F);
    x = x+dx;
    residual = abs(dx); // evaluate residual
    it = it+1; // increment iteration count
    if (it == itermax) // stop the the analysis of this step
    {
    cout<<"\n     (STEP 2) : Iteration reached to Ultimate "<<itermax<<" - error: "<<dx<<endl;
    cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
    break;
    }
    }
    Iep1 = (b*pow(x,3))/12 +b*x*pow(x-.5*x,2);
    Iep2=0;Iep3=0;
    for (j=1;j<=3;j++)
    {
    Iep2=Iep2+(2*n-1)*As[j]*pow(x-d[j],2);
    Iep3=Iep3+n*As[j+3]*pow(x-d[j+3],2);
    }
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
    FsSUM=0;FstanSUM=0;
    a=beta1*x;Cc=0.85*fc*a*b;Cctan=0;ec=0.5*(h-a);

    for (j=1;j<=6;j++)
    {
	es[j]=ecu*(d[j]-x)/x;
	if  (es[j]>=0 && es[j]<=ey)
        {fs[j]=Es*es[j];fstan[j]=(Es*ecu*d[j]/pow(x,2));}
    else if  (es[j]<0 && es[j]>=-ey)
        {fs[j]=Es*es[j];fstan[j]=(Es*ecu*d[j]/pow(x,2));}
    else if (es[j]>ey && es[j]<=esh)
        {fs[j]=fy;fstan[j]=0;}
    else if (es[j]<-ey && es[j]>=-esh)
        {fs[j]=-fy;fstan[j]=0;}
    else if (es[j]>esh && es[j]<=esu)
        {fs[j]=fy+Esh*(abs(es[j])-esh);fstan[j]=(Esh*ecu*d[j]/pow(x,2));}
    else if (es[j]<-esh && es[j]>=-esu)
        {fs[j]=-fy-Esh*(abs(es[j])-esh);fstan[j]=(Esh*ecu*d[j]/pow(x,2));}   
    else if (es[j]<-esu || es[j]>esu)
        {fs[j]=0;fstan[j]=0;}
	   
    if (d[j]>a)
        {Fs[j]=As[j]*fs[j];Fstan[j]=As[j]*fstan[j];}
    else if (d[j]<=a)
        {Fs[j]=As[j]*(fs[j]-0.85*fc);Fstan[j]=As[j]*(fstan[j]-0.85*fc);}
        
        FsSUM=FsSUM+Fs[j];
        FstanSUM=FstanSUM+Fstan[j];
    } // for

    A=Cc-FsSUM-Ptarget;
	A_tan=Cctan-FstanSUM;
    dx = (1/A_tan)*(-A);
    residual = abs(dx); // evaluate residual
    it = it + 1; // increment iteration count
    x = x+dx; // update x
        if (it == itermax) // stop the the analysis of this step
          {
		  cout<<"\n     (-) Iteration reached to Ultimate "<<it<<" - strain: "<<ecu<<" - error: "<<dx<<endl;
          cout<<"         ## The solution for this step is not converged. Please check your model ##"<<endl;
          exit(2);
          getch();
          }
    
   }
	        Pc1=x-.5*h;xU=x;
	        MomU=0;
            for (j=1;j<=6;j++){
			e[j]=d[j]-0.5*h;
			MomU=MomU+Fs[j]*e[j];
            }
	MomU=(MomU+Cc*ec+Pc1*Ptarget)*pow(10,-6);
	CurU=(ecu/x)*1000;
	if (it < itermax){// iteration control
    cout<<"\n     (STEP 3) : It is converged in "<<it<<" iterations - strain: "<<ecu<<" - x(mm): "<<x<<" - Phi(1/m): "<<CurU<<" - Moment (kN.m): "<<MomU<<"\n"<<endl;
    SC=CurU/CurY;OF=MomU/MomY;
    cout<<"     Section ductility ratio: "<<SC<<"  -   Over strength factor: "<<OF<<endl;
	// OUTPUT
	ofstream OUTo;OUTo.open("ConcreteSectionMomentCurvatureAnalysis-output.txt");
	OUTo<<"************************************************************\n";    
    OUTo<<"*             >> IN THE NAME OF ALLAH <<                   *\n";
    OUTo<<"*       Concrete Section Moment Curvature Analysis         *\n";
    OUTo<<"*               UNIT: [Newton-Millimeter]                  *\n";
    OUTo<<"*----------------------------------------------------------*\n";
    OUTo<<"*     This program is written by Salar Delavar Ghashghaei  *\n";  
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
    OUTo<<"(STEP 1) Crack::   x(mm):"<<xE<<" - Phi(1/m): "<<CurE<<" - Moment(kN/m):"<<MomE<<endl;
    OUTo<<"(STEP 2) Yield::   x(mm):"<<xY<<" - Phi(1/m): "<<CurY<<" - Moment(kN/m):"<<MomY<<endl;
    OUTo<<"(STEP 3) Ultimate::   x(mm):"<<xU<<" - Phi(1/m): "<<CurU<<" - Moment(kN/m):"<<MomU<<endl;
    OUTo<<"     Section ductility ratio: "<<SC<<"  -   Over strength factor: "<<OF<<endl;
    cout<< white <<"\n - Output data is written in file -"<<endl;
    OUTo.close();
    getch();
	}
	double X[3],Y[3];
	X[0]=CurE;X[1]=CurY;X[2]=CurU;
	Y[0]=MomE;Y[1]=MomY;Y[2]=MomU;
		char text1[50]="Concrete Section Moment Curvature Analysis",text2[40]="Curvature (1/m)",text3[40]="Moment (kN.m)";
		OUTPUT_HTML_GRAPH(X,Y,3,text1,text2,text3);
		system("start /w Graph-outputHTML.html");
		return 0;
	}
	void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]){
    // HTML GRAPH OUTPUT
	int i;
	double x,y,Xmax,Ymax;
	double *Xnew = new double [STEP];
	double *Ynew = new double [STEP];
	double *NorX = new double [STEP];
	double *NorY = new double [STEP];
	Xmax=MAX_ABS(X,n);
	Ymax=MAX_ABS(Y,n);
	Xnew[0]=0;Ynew[0]=0;
	for (i=0;i<n;i++){
	Xnew[i+1] = ABS(X[i]);
	Ynew[i+1] = ABS(Y[i]);
	}
	for (i=0;i<n+1;i++){
	NorX[i] = Xnew[i]/Xmax;
	NorY[i] = Ynew[i]/Ymax;
	//printf("\t %f   %f    \n",NorX[i],NorY[i]);
	}
	FILE *OutputFile;
	OutputFile = fopen(ShowText01, "w");
	fprintf(OutputFile,"<!DOCTYPE HTML><html><body style=\"background-color:black;\"><font color=\"white\"><head><script> \n");
	fprintf(OutputFile,"window.onload = function(){ \n");
	fprintf(OutputFile,"var canvas = document.getElementById(\"myCanvas\");var s1 = canvas.getContext(\"2d\");var s2 = canvas.getContext('2d'); \n");
	fprintf(OutputFile,"var s3 = canvas.getContext(\"2d\");var s4 = canvas.getContext(\"2d\");var s5 = canvas.getContext(\"2d\"); \n");
	fprintf(OutputFile,"var x=120,y=80,X,Y,Lx=1100,Ly=500,i; \n");
	fprintf(OutputFile,"s3.beginPath();s3.lineWidth = 3;s3.strokeStyle = \"cyan\";s3.rect(x,y,Lx,Ly); \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x+Lx*(i+1)*.1,y+Ly);s3.lineTo(x+Lx*(i+1)*.1,y+Ly-10);}; \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x,y+Ly*(i+1)*.1);s3.lineTo(x+10,y+Ly*(i+1)*.1);};s3.stroke();\n");
	fprintf(OutputFile,"s1.beginPath();s1.lineWidth = 3;s1.strokeStyle = \"yellow\"; \n");
	for (i=0;i<n;i++){
	fprintf(OutputFile,"s1.moveTo(%f,%f);",120+NorX[i]*1100,80+500-NorY[i]*500);
	fprintf(OutputFile,"s1.lineTo(%f,%f); \n",120+NorX[i+1]*1100,80+500-NorY[i+1]*500);
	}
	fprintf(OutputFile,"s1.stroke(); \n");
	fprintf(OutputFile,"s2.beginPath();s2.lineWidth = 1;s2.strokeStyle = \"cyan\";s2.setLineDash([5, 5]); \n");
	fprintf(OutputFile,"for(i=0;i<19;i++){s2.moveTo(x+Lx*(i+1)*.05,y);s2.lineTo(x+Lx*(i+1)*.05,y+Ly);} \n");
	fprintf(OutputFile,"s2.lineWidth = 1;s2.strokeStyle = \"cyan\";for(i=0;i<19;i++){s2.moveTo(x,y+Ly*(i+1)*.05);s2.lineTo(x+Lx,y+Ly*(i+1)*.05);} s2.stroke();\n");
	fprintf(OutputFile,"X=x+.25*Lx;Y=.7*y;s4.translate(X,Y);s4.font=\"30px serif\";s4.fillStyle = \"#7fff00\";s4.fillText(\"%s\",0,0); \n",text1);
	fprintf(OutputFile,"s4.save();X=-X+.2*x;Y=-Y+y+.6*Ly;s4.translate(X,Y);s4.rotate(3*Math.PI/2);s4.font=\"15px serif\"; \n");
	fprintf(OutputFile,"s4.fillStyle = \"#7fff00\";s4.textAlign = \"left\";s4.fillText(\"%s\",0,0);s4.restore(); \n",text3);
	fprintf(OutputFile,"s4.save();X=.2*Lx;Y=y+Ly-20;s4.translate(X,Y);s4.rotate(2*Math.PI);s4.font=\"15px serif\";s4.fillStyle = \"#7fff00\"; \n");
	fprintf(OutputFile,"s4.textAlign = \"left\";s4.fillText(\"%s\",0,0);s4.restore(); \n",text2);
	for(i=0;i<10;i++){
	x=.1*(i+1)*Xmax;
	fprintf(OutputFile,"s5.save();X=-.29*Lx+Lx*(%d+1)*.1;Y=.3*y+Ly+20;s5.rotate(2*Math.PI);s5.font=\"16px serif\"; \n",i);
	fprintf(OutputFile,"s5.fillStyle = \"#7fff00\";s5.textAlign = \"left\";s5.fillText(\"%.3e\",X,Y);s5.restore(); \n",x);
    }
    for(i=0;i<10;i++){
    y=.1*(i+1)*Ymax;
	fprintf(OutputFile,"s5.save();X=-.28*Lx-50;Y=Ly+.3*y-Ly*(%d+1)*.1;s5.rotate(2*Math.PI);s5.font=\"16px serif\"; \n",i);
	fprintf(OutputFile,"s5.fillStyle = \"#7fff00\";s5.textAlign = \"left\";s5.fillText(\"%.3e\",X,Y);s5.restore(); \n",y);
    }
	fprintf(OutputFile,"s5.save();X=-.25*Lx;Y=.3*y+Ly+20;s5.rotate(2*Math.PI);s5.font=\"16px serif\";s5.fillStyle = \"#7fff00\";s5.fillText(0,X,Y);s5.restore(); \n");
	fprintf(OutputFile,"s5.save();X=-.25*Lx-50;Y=Ly+.3*y;s5.rotate(2*Math.PI);s5.font=\"16px serif\";s5.fillStyle = \"#7fff00\";s5.textAlign = \"left\";s5.fillText(0,X,Y);s5.restore();}; \n");
	fprintf(OutputFile,"</script></head><body><canvas id=\"myCanvas\" width=\"1300\" height=\"1300\" style=\"border:1px solid black;\"></canvas></body></html> \n");
	fclose(OutputFile);
	free(Xnew);free(Ynew);free(NorX);free(NorY);
}
double MAX_ABS(double A[],int n){
	int i;
	double B[STEP];
	double Amax;
	// abs value
	for (i=0;i<n;i++){
	B[i] = A[i];
    if(B[i] < 0)
    B[i] = -B[i];
	}
	// Max of abs
	Amax = B[0];
	for (i=1;i<n;i++){
    if(Amax < B[i])
    Amax = B[i];
	}
	return Amax;
}
double ABS(double B){
	if (B < 0)
    B = -B;//Absolute number
    else
    B = B;
    return B;
}
