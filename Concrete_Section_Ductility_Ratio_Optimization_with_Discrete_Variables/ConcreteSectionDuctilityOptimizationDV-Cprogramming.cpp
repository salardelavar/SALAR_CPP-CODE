#include <conio.h>
#include <stdio.h>
#include <windows.h> // text color
#define N 10000
#define TEXT01 "ConcreteSectionDuctilityOptimizationDV-Input.csv"
#define TEXT02 "ConcreteSectionDuctilityOptimizationDV-InputSteelRebar.csv"
#define TEXT03 "ConcreteSectionDuctilityOptimizationDV-outputEXCEL.csv"
#define TEXT04 "ConcreteSectionDuctilityOptimizationDV-output.txt"
#define TEXT05 "               >> IN THE NAME OF ALLAH <<                 "
#define TEXT06 "      Concrete Section Ductility Ratio Optimization       "
#define TEXT07 "                 with Discrete Variables                  "
#define TEXT08 "                      UNIT: FREE UNIT                     "
#define TEXT09 "     This program is written by Salar Delavar Qashqai     "
#define TEXT10 "          E-mail: salar.d.ghashghaei@gmail.com            "
#define TEXT11 "      Notice: All input values must be positive           "
#define TEXT12 "  _    ______________________________________             "
#define TEXT13 "  |   |                                      |            "
#define TEXT14 "      |     #     #     #     #    #    #    |            "
#define TEXT15 "      |     #                           #    |            "
#define TEXT16 "  b   |    As1   As2   As3   As4  As5  As6   |            "
#define TEXT17 "      |     #                           #    |            "
#define TEXT18 "      |     #     #     #     #    #    #    |            "
#define TEXT19 "  _   |______________________________________|            "
#define TEXT20 "      |<-                 h                ->|            "
#define TEXT21 "      |<-d1->|                                            "
#define TEXT22 "      |<-  d2   ->|                                       "
#define TEXT23 "      |<-     d3      ->|                                 "
#define TEXT24 "      |<-        d4          ->|                          "
#define TEXT25 "      |<-            d5          ->|                      "
#define TEXT26 "      |<-               d6            ->|                 "
void IMPORT_DATA01(double AA[]);
void IMPORT_DATA02(double AS[6][N],double DIS[6][N],int &k);
void MessageInitialData(double AA[],int k);
void Distance(int i);
void ANALYSIS(double AA[],double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int &zMAX,int it,double SC,int k);
void OUTPUT_TXT(double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int zMAX,int it,double SC);
void OUTPUT_EXCEL(double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int zMAX,int it,double SC);
void textcolor(int ForgC);
void MessageErrorReportTEXT();
void MessageInputData();
void MessageAnalysisReport();
void MessageCheck_IMPORT_DATA01(double AA[]);
void MessageCheck_IMPORT_DATA02(double AS[6][N],double DIS[6][N],int k);
double ABS(double B);
double SQRT2(double D);

int main(){
	double AS[6][N],DIS[6][N];
    double AA[14],SC;
	int k,zMAX,it;
	textcolor(11);
        IMPORT_DATA01(AA);
        IMPORT_DATA02(AS,DIS,k);
        MessageCheck_IMPORT_DATA01(AA);
        MessageCheck_IMPORT_DATA02(AS,DIS,k);
        MessageInitialData(AA,k);
        double SSCC[k],OOFF[k],CurYY[k],MomYY[k],CurUU[k],MomUU[k];
    textcolor(10);    
    ANALYSIS(AA,AS,DIS,CurYY,CurUU,SSCC,MomYY,MomUU,OOFF,zMAX,it,SC,k);
    OUTPUT_TXT(AS,DIS,CurYY,CurUU,SSCC,MomYY,MomUU,OOFF,zMAX,it,SC);
    OUTPUT_EXCEL(AS,DIS,CurYY,CurUU,SSCC,MomYY,MomUU,OOFF,zMAX,it,SC);
    printf("\a\n  - Output data is written in file -\n");
    getch();
    return 0;
	}
	
void textcolor(int ForgC){
     WORD wColor;
     //This handle is needed to get the current background attribute

     HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
     CONSOLE_SCREEN_BUFFER_INFO csbi;
     //csbi is used for wAttributes word

     if(GetConsoleScreenBufferInfo(hStdOut, &csbi)){
          //To mask out all but the background attribute, and to add the color
          wColor = (csbi.wAttributes & 0xF0) + (ForgC & 0x0F);
          SetConsoleTextAttribute(hStdOut, wColor);
     }
     return;
}
void IMPORT_DATA01(double AA[]){
	int i=0;
	FILE *InputFile;
	InputFile = fopen(TEXT01, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [%s] \n",TEXT01);
		Sleep(6000);
		exit(1);
	}
	char line[100],a[100];
	while(i < N && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	AA[i]= atof(a);
	i++;
	}
}

void IMPORT_DATA02(double AS[6][N],double DIS[6][N],int &k){
	int i = 0;char CH[100];
	FILE *InputFile;
	InputFile = fopen(TEXT02, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [ConcreteSectionDuctilityOptimizationDV-inputSteelRebar.csv] \n");
		Sleep(6000);
		exit(1);
	}
	char line[100];
	do{
	fscanf(InputFile,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&AS[0][i],&DIS[0][i],&AS[1][i],&DIS[1][i],&AS[2][i],&DIS[2][i],&AS[3][i],&DIS[3][i],&AS[4][i],&DIS[4][i],&AS[5][i],&DIS[5][i]);
	/*
	printf("\n          AS1[%d]: #.3e - ",i+1,AS[0][i]);printf("D1[%d]: #.3e\n",i+1,DIS[0][i]);
	printf("          AS2[%d]: #.3e - ",i+1,AS[1][i]);printf("D2[%d]: #.3e\n",i+1,DIS[1][i]);
	printf("          AS3[%d]: #.3e - ",i+1,AS[2][i]);printf("D3[%d]: #.3e\n",i+1,DIS[2][i]);
	printf("          AS4[%d]: #.3e - ",i+1,AS[3][i]);printf("D4[%d]: #.3e\n",i+1,DIS[3][i]);
	printf("          AS5[%d]: #.3e - ",i+1,AS[4][i]);printf("D5[%d]: #.3e\n",i+1,DIS[4][i]);
	printf("          AS6[%d]: #.3e - ",i+1,AS[5][i]);printf("D6[%d]: #.3e\n",i+1,DIS[5][i]);
	*/
	i++;	
	}
	while(i < N && fgets(line,sizeof(line),InputFile) != NULL);
    k = i;
}
void MessageInitialData(double AA[],int k){	
    char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT05,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT06,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT07,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT08,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);   
	printf("\t\t\t\t%c%s%c\n",Qf,TEXT09,Qf);
	printf("\t\t\t\t%c%s%c\n",Qf,TEXT10,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT11,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT12,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT13,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT14,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT15,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT16,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT17,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT18,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT19,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT20,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT21,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT22,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT23,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT24,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT25,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,TEXT26,Qf);
                
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

	MessageInputData();
	printf("      Axial load [+ : Compression]:                   %.3e\n",AA[0]);
	printf("      Concrete section width(b):                      %.3e\n",AA[1]);
	printf("      Concrete section height(h):                     %.3e\n",AA[2]);
	printf("      Minimum section ductility ratio:                %.3e\n",AA[3]);
	printf("      Maximum section ductility ratio:                %.3e\n",AA[4]);
	printf("      Concrete strength(fc):                          %.3e\n",AA[5]);
	printf("      Ultimate concrete strain(ecu):                  %.3e\n",AA[6]);
	printf("      Yield rebar strength(fy):                       %.3e\n",AA[7]);
	printf("      Ultimate rebar strength(fu):                    %.3e\n",AA[8]);
	printf("      Rebar elastic modulus(Es):                      %.3e\n",AA[9]);
	printf("      Strain at steel strain-hardening:               %.3e\n",AA[10]);
	printf("      Ultimate steel strain:                          %.3e\n",AA[11]);
	printf("      Maximum number of iterations:                   %.3e\n",AA[12]);// maximum number of iterations
	printf("      Specified tolerance for convergence:            %e\n",AA[13]);//  specified tolerance for convergence
	printf("      Number of steel rebar increment of input value: %d\n",k);	
}	
void Distance(int i){
	        if (i < 10)
		    printf("\t");
	        if (i >= 10 && i <= 99)
		    printf("\b\t\b");
		    if (i >= 100 && i <= 999)
		    printf("\b\t\b\b");
		    if (i >= 1000 && i <= 9999)
		    printf("\b\t\b\b\b");
		    if (i >= 10000 && i <= 20000)
		    printf("\b\t\b\b\b\b");
}
void OUTPUT_TXT(double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int zMAX,int it,double SC){
	// OUTPUT
	int i,z;
	FILE *OutputFile;
	OutputFile = fopen(TEXT04, "w");
	fprintf(OutputFile,"************************************************************\n");   
    fprintf(OutputFile,"*                >> IN THE NAME OF ALLAH <<                *\n");
    fprintf(OutputFile,"*     Concrete Section Ductility Ratio Optimization        *\n");
    fprintf(OutputFile,"*                 with Discrete Variables                  *\n");
    fprintf(OutputFile,"*                     UNIT: FREE UNIT                      *\n");
    fprintf(OutputFile,"*----------------------------------------------------------*\n");
    fprintf(OutputFile,"*      This program is written by Salar Delavar Qashqai    *\n"); 
    fprintf(OutputFile,"*          E-mail: salar.d.ghashghaei@gmail.com            *\n");
    fprintf(OutputFile,"*----------------------------------------------------------*\n");
    fprintf(OutputFile,"*      Notice: All input values must be positive           *\n");
    fprintf(OutputFile,"*  _    ______________________________________             *\n");
    fprintf(OutputFile,"*  |   |                                      |            *\n");
    fprintf(OutputFile,"*      |     #     #     #     #    #    #    |            *\n");
    fprintf(OutputFile,"*      |     #                           #    |            *\n");
    fprintf(OutputFile,"*  b   |    As1   As2   As3   As4  As5  As6   |            *\n");
    fprintf(OutputFile,"*      |     #                           #    |            *\n");
    fprintf(OutputFile,"*  |   |     #     #     #     #    #    #    |            *\n");
    fprintf(OutputFile,"*  _   |______________________________________|            *\n");
    fprintf(OutputFile,"*      |<-                 h                ->|            *\n");
    fprintf(OutputFile,"*      |<-d1->|                                            *\n");
    fprintf(OutputFile,"*      |<-  d2   ->|                                       *\n");
    fprintf(OutputFile,"*      |<-     d3      ->|                                 *\n");
    fprintf(OutputFile,"*      |<-        d4          ->|                          *\n");
    fprintf(OutputFile,"*      |<-            d5          ->|                      *\n");
    fprintf(OutputFile,"*      |<-               d6            ->|                 *\n");
    fprintf(OutputFile,"************************************************************\n\n");

    fprintf(OutputFile,"   ------------------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(OutputFile,"     Section    Yield-Curvature   Ultimate-Curvature   Section-ductility-ratio  Yield-Moment   Ulitmate-Moment   Over-strength-factor  \n");
    fprintf(OutputFile,"   ------------------------------------------------------------------------------------------------------------------------------------\n");
    
	for (z=1;z<=zMAX;z++)
    fprintf(OutputFile,"\t%d\t   %.4e\t\t%.3e\t\t%.3e\t %.3e\t    %.3e\t     %.3e\n",z,CurYY[z],CurUU[z],SSCC[z],MomYY[z],MomUU[z],OOFF[z]);
	fprintf(OutputFile,"  ### The Optimum section ductility ratio for these constrains: #.3e with %d iterations ###\n\n",SC,it); 
	fprintf(OutputFile,"      Section number %d\n",z-1);
	fprintf(OutputFile,"   --------------------------------------------------\n");
	for (i=0;i<6;i++)
	fprintf(OutputFile,"      Rebar Area(%d): %.3e - Distance(%d): %.3e\n",i+1,AS[i][z],i+1,DIS[i][z]);
    fclose(OutputFile);	
}
void ANALYSIS(double AA[],double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int &zMAX,int it,double SC,int k){
	double Ptarget,b,h,fc,ecu,fy,fu,SCmin,SCmax,OF;
	double Es,Esh,x,beta1,tolerance,a,ec,Finitsum,Fsum;
	double ey,esh,esu,residual,A,A_tan,Cc,Cctan,FsSUM,FstanSUM,SUM,dx;
	double Ie,fr,Ec,ecr,xE,xY,xU,CurE,MomE,Finit,F,n,Pc1,CurY,MomY,CurU,MomU;
	int i,j,z,itermax,M,ZA,X=0;
	
		double fs[6],fstan[6],Fs[6],Fstan[6],e[6],es[6];
	    double Iep,Iep1,Iep2,Iep3;
Ptarget = AA[0];
b = AA[1];
h = AA[2];
SCmin = AA[3];
SCmax = AA[4];
fc = AA[5];
ecu = AA[6];
fy = AA[7];
fu = AA[8];
Es = AA[9];
esh = AA[10];
esu = AA[11];
itermax = AA[12];
tolerance = AA[13];




    MessageAnalysisReport();
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
    Ie=(b*h*h*h)/12;fr=.7*SQRT2(fc);Ec=5000*SQRT2(fc);ecr=fr/Ec;xE=x;
    MomE=(fr*Ie/x);// Elastic Moment
    CurE=((fr*Ie/x)/(Ec*Ie))*1000;// Elastic Curvature

    printf("   ------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("     Section    Yield-Curvature   Ultimate-Curvature   Section-ductility-ratio  Yield-Moment   Ulitmate-Moment   Over-strength-factor  \n");
    printf("   ------------------------------------------------------------------------------------------------------------------------------------\n");
    
    for (z=1;z<=k;z++){
    // Yield capacity
    n=Es/Ec;
    Finitsum=0;
    it = 0; // initialize iteration count
    residual = 100; /// initialize residual
    for (i=0;i<6;i++)
    es[i]= (ey*(DIS[i][z]-x))/(DIS[5][z]-x);
    Finit = 0;
    for (i=0;i<5;i++)
    Finit += es[i]*AS[i][z];
    Finit += ey*AS[5][z];
    Finit *= Es;
    Finit += 0.5*fc*b*x+Ptarget;
    while ( residual > tolerance ){
    for (i=0;i<6;i++)
    es[i]= (ey*(DIS[i][z]-x))/(DIS[5][z]-x);
    F = 0;
    for (i=0;i<5;i++)
    F += es[i]*AS[i][z];
    F += ey*AS[5][z];
    F *= Es;
    F += 0.5*fc*b*x+Ptarget;
    dx = (1/Finit)*(F);
    x = x+dx;
    residual = abs(dx); // evaluate residual
    it = it+1; // increment iteration count
    if (it == itermax){ // stop the the analysis of this step
  
    printf("\n     (STEP 2) : Iteration reached to Ultimate %d  - strain: #.3e  - residual: #.3e",it,ecu,residual);
    printf("         ## The solution for this step is not converged. Please check your model ##");
    break;
    }//if 
    }//while
    Iep1 = (b*x*x*x)/12 +b*x*(x-.5*x)*(x-.5*x);
    Iep2 = 0;
    for (i=0;i<3;i++)
    Iep2 += (2*n-1)*AS[i][z]*(x-DIS[i][z])*(x-DIS[i][z]);
    Iep3 = 0;
    for (i=0;i<3;i++)
    Iep3=n*AS[i][z]*(x-DIS[i][z])*(x-DIS[i][z]);

    Iep=Iep1+Iep2+Iep3;xY=x;
    Pc1=x-.5*h;// Distance from applied force from neutral axis
    MomY = ((.5*fc*Iep/x)-Ptarget*Pc1);// Yield Moment
    CurY =((.5*fc*Iep/x)/(Ec*Iep))*1000;// Yield Curvature
        if (it < itermax)// iteration control
        {
        //printf("     (STEP 2) : It is converged in %d iterations - strain: #.3e - x(mm): #.3e - Phi(1/m): #.3e - Moment (kN.m): #.3e\n",it,ey,x,CurY,MomY);
        } 
    
    // Ultimate capacity
    it = 0; // initialize iteration count
    residual = 100; // initialize residual
    while (residual > tolerance){
    a=beta1*x;Cc=0.85*fc*a*b;Cctan=0;ec=0.5*(h-a);
    FsSUM=0;
    FstanSUM=0;
	for (i=0;i<6;i++){
    es[i]= (ecu*(DIS[i][z]-x))/x;
    
	if  (es[i]>=-ey && es[i]<=ey)
        {fs[i]=Es*es[i];fstan[i]=(Es*ecu*DIS[i][z]/(x*x));}
    else if (es[i]>ey && es[i]<=esh)
        {fs[i]=fy;fstan[i]=0;}
    else if (es[i]<-ey && es[i]>=-esh)
        {fs[i]=-fy;fstan[i]=0;}
    else if (es[i]>esh && es[i]<=esu)
        {fs[i]=fy+Esh*(abs(es[i])-esh);fstan[i]=(Esh*ecu*DIS[i][z]/(x*x));}
    else if (es[i]<-esh && es[i]>=-esu)
        {fs[i]=-fy-Esh*(abs(es[i])-esh);fstan[i]=(Esh*ecu*DIS[i][z]/(x*x));}   
    else if (es[i]<-esu || es[i]>esu)
        {fs[i]=0;fstan[i]=0;}
       
    if (DIS[i][z]>a)
        {Fs[i]=AS[i][z]*fs[i];Fstan[i]=AS[i][z]*fstan[i];}
    else if (DIS[i][z]<=a)
        {Fs[i]=AS[i][z]*(fs[i]-0.85*fc);Fstan[i]=AS[i][z]*(fstan[i]-0.85*fc);}
        
         FsSUM += Fs[i];
         FstanSUM += Fstan[i];
    }// for

    A=Cc-FsSUM-Ptarget;
	A_tan=Cctan-FstanSUM;
    dx = A/A_tan;
    residual = abs(dx); // evaluate residual
    it = it + 1; // increment iteration count
    x -= dx; // update x
        if (it == itermax){ // stop the the analysis of this step please of Convergence
		  printf("\n     (STEP 3) : Iteration reached to Ultimate %d  - strain: #.3e  - residual: #.3e",it,ecu,residual);
          printf("         ## The solution for this step is not converged. Please check your model ##");
          exit(2);
          getch();
          }
    
	}// while
    Pc1=x-.5*h;xU=x;
    SUM=0;
    for (i=0;i<6;i++){
    e[i]=DIS[i][z]-0.5*h;
	SUM += Fs[i]*e[i];
	}
    
    
	MomU=(SUM+Cc*ec+Pc1*Ptarget);
	CurU=(ecu/x)*1000;
	if (it < itermax){// iteration control
    //printf("     (STEP 3) : It is converged in %d iterations - strain: #.3e - x(mm): #.3e - Phi(1/m): #.3e - Moment (kN.m): #.3e\n",it,ecu,x,CurU,MomU);
    SC=CurU/CurY;OF=MomU/MomY;zMAX=z;
    }
    SSCC[z]=SC;OOFF[z]=OF;CurYY[z]=CurY;MomYY[z]=MomY;CurUU[z]=CurU;MomUU[z]=MomU;
    Distance(z);
    printf("%d\t   %.4e\t\t%.3e\t\t%.3e\t %.3e\t    %.3e\t     %.3e\n",z,CurYY[z],CurUU[z],SSCC[z],MomYY[z],MomUU[z],OOFF[z]);
    if (SC>=SCmin && SC<=SCmax){
    printf("  ### The Optimum section ductility ratio for these constrains: %.3e with %d iterations ###\n\n",SC,it); 
	printf("      Section number %d\n",z);
	printf("   --------------------------------------------------\n");
	for (i=0;i<6;i++)
	printf("      Rebar Area(%d): %.3e - Distance(%d): %.3e\n",i+1,AS[i][z],i+1,DIS[i][z]);	
     break;
    }
    if (z==k && SC<=SCmin && SC>=SCmax){
    printf("  All steel areas and distances of input values have been checked and couldn't find the best result.");
    printf("  Number of increments: %d",k);
    printf("  It is recommended to increase number input values -> [%s]",TEXT02);
	}
    } // for	
}
void MessageErrorReportTEXT(){
  int i;
  char Ql;
  Ql=176;
  textcolor(12);
  printf("\a\n     ");
  for (i=1;i<50;i++)
  printf("%c",Ql);
  printf(" Error Report ");
  for (i=1;i<50;i++)
  printf("%c",Ql); 
  printf("\n");
}
double ABS(double B){
	if (B < 0)
    B = -B;//Absolute number
    else
    B = B;
    return B;
}
void MessageAnalysisReport(){
  int i;
  char Ql=176;
  printf("\n     ");
  for (i=1;i<50;i++)
  printf("%c",Ql);
  printf(" Analysis Report ");
  for (i=1;i<50;i++)
  printf("%c",Ql); 
  printf("\n");	
}
void MessageInputData(){
  int i;
  char Ql=176;
  printf("\n     ");
  for (i=1;i<50;i++)
  printf("%c",Ql);
  printf(" Input Data ");
  for (i=1;i<50;i++)
  printf("%c",Ql); 
  printf("\n");	
}
double SQRT2(double D){           
            int it,itermax;
            double residual,tolerance,x,dx,dx_ABS,f,df;
			it = 0; // initialize iteration count
		    itermax = 100000;
		    residual = 100; // initialize residual
		    tolerance = 1e-8;
		    x = 1;// initialize answer
		    while (residual > tolerance){    
		    	f = x*x - D;
		    	df = 2 * x;
		    	dx = f/df;
		        x = x - dx;
                residual = ABS(dx); // abs residual
		        it += 1; // increment iteration count
		        //printf("f: #.3e -\tdx: #.3e -\tresidual: #.3e\n",f,dx,residual);
		         if (it == itermax){
		          //printf("\tSQRT2(number,power) : SQRT2(#.3e) - iteration: %d ->   ## The solution is not converged ##\n",D,it);
		          break;
				}
		    }
		        if (it < itermax){
			   //printf("\tSQRT(number,power) - SQRT(#.3e,#.3e) : #.3e \n",D,n, x);
			   return x;
			   }
}
void MessageCheck_IMPORT_DATA01(double AA[]){
	int i;
	for(i=0;i<14;i++){
	if (AA[i] < 0){
	MessageErrorReportTEXT();
	printf("          Please check the input file! -> [ %s ]\n",TEXT01);
	printf("          Row %d has a negative value: #.3e\n",i+1,AA[i]);
	printf("                            *** Negative data input value is not acceptable ***\n");
	Sleep(40000);
	exit(1);
	}
	}	 
}
void MessageCheck_IMPORT_DATA02(double AS[6][N],double DIS[6][N],int k){
	int i;
	for(i=0;i<k;i++){
	if (AS[0][i] < 0 || DIS[0][i] < 0 || AS[1][i] < 0 || DIS[1][i] < 0 || AS[2][i] < 0 || DIS[2][i] < 0 || AS[3][i] < 0 || DIS[3][i] < 0 || AS[4][i] < 0 || DIS[4][i] < 0 || AS[5][i] < 0 || DIS[5][i] < 0 ){
	MessageErrorReportTEXT();
	printf("          Please check the input file! -> [ %s ]\n",TEXT02);
	printf("          Row %d has a negative value.\n",i+1);
	printf("\n          AS1[%d]: #.3e - ",i+1,AS[0][i]);printf("D1[%d]: #.3e\n",i+1,DIS[0][i]);
	printf("          AS2[%d]: #.3e - ",i+1,AS[1][i]);printf("D2[%d]: #.3e\n",i+1,DIS[1][i]);
	printf("          AS3[%d]: #.3e - ",i+1,AS[2][i]);printf("D3[%d]: #.3e\n",i+1,DIS[2][i]);
	printf("          AS4[%d]: #.3e - ",i+1,AS[3][i]);printf("D4[%d]: #.3e\n",i+1,DIS[3][i]);
	printf("          AS5[%d]: #.3e - ",i+1,AS[4][i]);printf("D5[%d]: #.3e\n",i+1,DIS[4][i]);
	printf("          AS6[%d]: #.3e - ",i+1,AS[5][i]);printf("D6[%d]: #.3e\n",i+1,DIS[5][i]);
	printf("                            *** Negative data input value is not acceptable ***\n");
	Sleep(40000);
	exit(1);
	}
	}	 
}
void OUTPUT_EXCEL(double AS[6][N],double DIS[6][N],double CurYY[],double CurUU[],double SSCC[],double MomYY[],double MomUU[],double OOFF[],int zMAX,int it,double SC){
	// EXCEL OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(TEXT03, "w");
	fprintf(OutputFile,"           ### Output Concrete Section Ductility Ratio Optimization with Discrete Variables ###\n");
    fprintf(OutputFile,"Section,Yield Curvature,Ultimate Curvature,Section Ductility Ratio,Yield Moment,Ulitmate Moment,Over Strength Factor  \n");
	for (i=1;i<=zMAX;i++)
    fprintf(OutputFile,"%d,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n",i,CurYY[i],CurUU[i],SSCC[i],MomYY[i],MomUU[i],OOFF[i]);	
    fclose(OutputFile);
}
