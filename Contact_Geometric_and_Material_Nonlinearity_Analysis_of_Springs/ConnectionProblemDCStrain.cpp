#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define N 10000
#define S01 "                            >> IN THE NAME OF ALLAH <<                             "
#define S02 " Geometric and Material Nonlinearity Analysis of Springs with Displacement Control "
#define S03 "                                UNIT: Free Unit                                    "
#define S04 "               This program is written by Salar Delavar Qashqai                    "
#define S05 "                     E-mail: salar.d.ghashghaei@gmail.com                          "
#define ShowText01 "ConnectionProblemDCStrain-inputDATA.csv"
#define ShowText02 "ConnectionProblemDCStrain-inputSPRING.csv" 
#define ShowText03 "ConnectionProblemDCStrain-outputEXCEL.csv"
#define ShowText04 "ConnectionProblemDCStrain-outputMATLAB.m"

void IMPORT_DATA01(double Import_Data[]);
void IMPORT_DATA02(double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int &k);
void MessageNegative_IMPORT_DATA01(double Import_Data[]);
void MessageNegative_IMPORT_DATA02(double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int n);
void MessageInputDataTEXT();
void MessageAnalysisReportTEXT();
void MessageErrorReportTEXT();
void MessageInitialData(double Import_Data[],int n);
void MessageStrainStressTEXT(double STARIN01[],double STRESS01[],double STARIN02[],double STRESS02[],int n);
void textcolor(int ForgC);
void Distance(int i);
void Element_Slope(double A[],double B[],double C[],int n);
void Element_Stiffness(double A[],double B[],double C[],double D[],double E[],double F[],double K[],int n,int m);
void ANALYSIS(double Import_Data[],double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int n,int m);
void OUTPUT_excel(double A[],double B[],int n);
void OUTPUT_matlab(double A[],double B[],int n);
double ABS(double B);
double Horizental_Disp(double L,double &x,double y);
double SQRT2(double D);
int main(){
	int k;
	double Import_Data[10];
	double STRAIN01[10];
	double STRESS01[10];
	double STRAIN02[10];
	double STRESS02[10];
	IMPORT_DATA01(Import_Data);
	MessageNegative_IMPORT_DATA01(Import_Data);
	IMPORT_DATA02(STRAIN01,STRESS01,STRAIN02,STRESS02,k);
	MessageNegative_IMPORT_DATA02(STRAIN01,STRESS01,STRAIN02,STRESS02,k);
	int m = 1 + ABS(Import_Data[9]/Import_Data[7]);
	textcolor(10);
	MessageInitialData(Import_Data,m);
	MessageStrainStressTEXT(STRAIN01,STRESS01,STRAIN02,STRESS02,k);
	textcolor(14);
	ANALYSIS(Import_Data,STRAIN01,STRESS01,STRAIN02,STRESS02,k,m);
    getch();
	return 0;
	}
	void IMPORT_DATA01(double Import_Data[]){
	int i=0;
	FILE *InputFile;
	InputFile = fopen(ShowText01, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [%s] \n",ShowText01);
		Sleep(6000);
		exit(1);
	}
	char line[100],a[100];
	while(i < N && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
}

void IMPORT_DATA02(double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int &k){
	int i = 0;
	FILE *InputFile;
	InputFile = fopen(ShowText02, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [%s] \n",ShowText02);
		Sleep(6000);
		exit(1);
	}
	char line[1000];
	do{	
	fscanf(InputFile,"%lf,%lf,%lf,%lf",&STRAIN01[i],&STRESS01[i],&STRAIN02[i],&STRESS02[i]);
	//printf("%d - STRAIN01[%d]: %lf - STRESS01[%d]: %lf- STRAIN02[%d]: %lf - STRESS02[%d]: %lf\n",i,i,STRAIN01[i],i,STRESS01[i],i,STRAIN02[i],i,STRESS02[i]);
	i++;	
	}
	while(i < N && fgets(line,sizeof(line),InputFile) != NULL);
	k = i-1;
    //printf("%d\n",k);
}
void MessageNegative_IMPORT_DATA01(double Import_Data[]){
		if ( Import_Data[0] < 0 || Import_Data[1] < 0 || Import_Data[2] < 0 || Import_Data[3] < 0 || Import_Data[4] < 0 || Import_Data[9] < 0){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative data input value is not acceptable ***\n");
	printf("  Length of rigid element:                 %f\n",Import_Data[0]);
	printf("  Spring length 1:                         %f\n",Import_Data[1]);
	printf("  Spring length 2:                         %f\n",Import_Data[2]);
	printf("  Spring area 1:                           %f\n",Import_Data[3]);
	printf("  Spring area 2:                           %f\n",Import_Data[4]);
	printf("  External applied load in spring 1:       %f\n",Import_Data[5]);
	printf("  External applied load in spring 2:       %f\n",Import_Data[6]);
	printf("  Initial applied displacement:            %f\n",Import_Data[7]);
	printf("  Distance of spring 1 from rigid element: %f\n",Import_Data[8]);
	printf("  Ultimate absolute displacement:          %f\n",Import_Data[9]);
	Sleep(40000);
	exit(1);	
	}

}
void MessageNegative_IMPORT_DATA02(double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int n){
	int i;
	for (i=0;i<n+1;i++){
		if ( STRAIN01[i] < 0 ||  STRESS01[i] < 0 ||  STRAIN02[i] < 0 ||  STRESS02[i] < 0 ){
    MessageErrorReportTEXT();
    printf("          Please check this file! -> [%s]\n",ShowText02);
    printf("          Row %d has a negative value.\n",i);
	printf("                            *** Negative data input value is not acceptable ***\n");
	printf("  Strain of spring 1:                 %f\n",STRAIN01[i]);
	printf("  Stress of spring 1:                 %f\n",STRESS01[i]);
	printf("  Strian of spring 2:                 %f\n",STRAIN02[i]);
	printf("  Stress of spring 2:                 %f\n",STRESS02[i]);
	Sleep(40000);
	exit(1);	
	 }	
	}

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
void MessageInputDataTEXT(){
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
void MessageAnalysisReportTEXT(){   
	  int i;
	  char Ql=176;
	  printf("\n     ");
	  for (i=1;i<47;i++)
	  printf("%c",Ql);
	  printf(" Analysis Report ");
	  for (i=1;i<47;i++)
	  printf("%c",Ql); 
	  printf("\n");	
	printf("\n\t   ");
    for (i=1;i<72;i++)
    printf("-");
    printf("\n");
    printf("\t     Increment          Base Shear      Incremental Displacement[DOF(1)]   \n");
    printf("\t   ");
    for (i=1;i<72;i++)
    printf("-");
    printf("\n");	
}
void Distance(int i){
	        if (i < 10)
		    printf("\b");
	        if (i >= 10 && i <= 99)
		    printf("\b\t\b");
		    if (i >= 100 && i <= 999)
		    printf("\b\t\b\b");
		    if (i >= 1000 && i <= 9999)
		    printf("\b\t\b\b\b");
		    if (i >= 10000 && i <= 20000)
		    printf("\b\t\b\b\b\b");
}
double ABS(double B){
	if (B < 0)
    B = -B;//Absolute number
    else
    B = B;
    return B;
}
void Element_Slope(double A[],double B[],double C[],int n){   
	int i;
    C[0]=B[0]/A[0];
    //printf("%d - Rk[%d]: %f\n",0,0,C[0]);
    for (i=1;i<n;i++){
    //cout<<i<<" - TET["<<i<<"]:"<<A[i]<<" - MOM["<<i<<"]:"<<B[i]<<endl;
    C[i]=(B[i]-B[i-1])/(A[i]-A[i-1]);
    //printf("%d - Strain[%d]: %f - Stress[%d]: %f - Rk[%d]: %f\n",i,i,A[i],i,B[i],i,C[i]);
	}	
}
void Element_Stiffness(double A[],double B[],double C[],double D[],double E[],double F[],double K[],int n,int m){
            int i;
            if (ABS(D[n])>= 0 && ABS(D[n])<= A[0])
            K[n] = (C[n]*E[n])/F[n];
            for (i=0;i<m;i++){
            if (ABS(D[n])> A[i] && ABS(D[n])<= A[i+1])
			K[n] = ((B[i]+C[i+1]*(ABS(D[n])-A[i]))/ABS(D[n]))*E[n]/F[n];	
			}
			if (ABS(D[n]) > A[m-1])
            K[n] = 0;
}
void ANALYSIS(double Import_Data[],double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int n,int m){
double L,Dx,a=0,Le[2],Area[2],Force[2],es[2],K[2],SLOPE01[10],SLOPE02[10],Dini,Dmax,Du,up;
int z,zMAX,well_done=0;
double output_u[N],output_base[N];
L = Import_Data[0];
Le[0] = Import_Data[1];
Le[1] = Import_Data[2];
Area[0] = Import_Data[3];
Area[1] = Import_Data[4];
Force[0] = Import_Data[5];
Force[1] = Import_Data[6];
Dini = Import_Data[7];
Du = Import_Data[8];
Dmax = Import_Data[9];
double F[2],Ay;
Element_Slope(STRAIN01,STRESS01,SLOPE01,n);
Element_Slope(STRAIN02,STRESS02,SLOPE02,n);
    MessageAnalysisReportTEXT();
    for (z=0;z<m;z++){
        up = Dini*(z+1);// Define the applied Displacement
        Horizental_Disp(L,Dx,up);
	    a =  SQRT2(.5*up*.5*up + .5*Dx*.5*Dx);// spring strain 1 during analysis
		es[0] = a/Le[0]; 
     Element_Stiffness(STRAIN01,STRESS02,SLOPE01,es,Area,Le,K,0,n);
	  F[0] = (((L-Dx)/L)*((L-Dx)/L))*K[0]*0.5*up + Import_Data[5];
	  Ay = F[0];
	  //printf("\t+ %f - %f\n",es[0],F[0]);
	  if (ABS(up) >= Du){
	  es[1] = (Du-up)/Le[1];	
	  Element_Stiffness(STRAIN02,STRESS02,SLOPE02,es,Area,Le,K,1,n);
	  F[1] = (((L-Dx)/L)*((L-Dx)/L))*K[1]*(Du-up) + Import_Data[6];
	  Ay = F[0] + F[1];
	  //printf("\t= %f - %f\n",es[1],F[1]);	
	  } 
      //printf("\t %f - %f\n",up,Dx);
      zMAX=z+1;
      output_u[z]=up;
	  output_base[z]=Ay;
	  Distance(z);
	  printf("\t\t%d\t\t%e\t\t%e\n",z+1,output_base[z],output_u[z]); 
	  	if (ABS(up) >= Dmax){
	  		well_done = 1;
	textcolor(13);  		
	printf("\n\t\t  ## Increment displacement reach to ultimate displacement ##\n\n");	
	break;	
	    }
    }// for
    	if (well_done == 1){
		OUTPUT_excel(output_u,output_base,zMAX);
		OUTPUT_matlab(output_u,output_base,zMAX);
		textcolor(15);
		printf("\n\a - Output data is written in Excel and Matlab file -"); 
		}
}
void OUTPUT_matlab(double A[],double B[],int n){
// MATLAB OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText04, "w");
fprintf(OutputFile," %% Geometric and Material Nonlinearity Analysis of Springs with Displacement Control %%\n");
fprintf(OutputFile,"disp_Dof1=[0\n");
for(i=0;i<=n;i++)
fprintf(OutputFile,"%e\n",A[i]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"base_shear=[0\n");
for(i=0;i<=n;i++)
fprintf(OutputFile,"%e\n",B[i]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"figure(1)\n");
fprintf(OutputFile,"plot(disp_Dof1,base_shear,'LineWidth',3);\n");
fprintf(OutputFile,"title(['# BASE SHEAR - DISPLACEMENT DIAGRAM #'],'Color','b');\n");
fprintf(OutputFile,"xlabel('DISPLACEMENT [DOF(1)]');ylabel('BASE SHEAR');grid on;\n");
fclose(OutputFile);		
}
void OUTPUT_excel(double A[],double B[],int n){
// EXCEL OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText03, "w");
fprintf(OutputFile," ### Geometric and Material Nonlinearity Analysis of Springs with Displacement Control ###\n");
fprintf(OutputFile,"Increment,Displacement [DOF(1)],Base Shear[DOF(1)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e\n",i+1,A[i],B[i]);
fclose(OutputFile);
}
void MessageInitialData(double Import_Data[],int n){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<84;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c%s%c\n",Qf,S01,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,S02,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,S03,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<84;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);   
	printf("\t\t\t\t%c%s%c\n",Qf,S04,Qf);
	printf("\t\t\t\t%c%s%c\n",Qf,S05,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<84;i++)
    printf("%c",Qb);
    printf("%c\n",Qe); 
	  
    MessageInputDataTEXT();
	printf("\t  Length of rigid element:                 %.3e\n",Import_Data[0]);
	printf("\t  Spring length 1:                         %.3e\n",Import_Data[1]);
	printf("\t  Spring length 2:                         %.3e\n",Import_Data[2]);
	printf("\t  Spring area 1:                           %.3e\n",Import_Data[3]);
	printf("\t  Spring area 2:                           %.3e\n",Import_Data[4]);
	printf("\t  External applied load in spring 1:       %.3e\n",Import_Data[5]);
	printf("\t  External applied load in spring 2:       %.3e\n",Import_Data[6]);
	printf("\t  Initial applied displacement:            %.3e\n",Import_Data[7]);
	printf("\t  Distance of spring 1 from rigid element: %.3e\n",Import_Data[8]);
	printf("\t  Ultimate absolute displacement:          %.3e\n",Import_Data[9]);
    printf("\t  Number of increments:                    %d\n",n);	
}
void MessageStrainStressTEXT(double STRAIN01[],double STRESS01[],double STRAIN02[],double STRESS02[],int n){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    int  BB=201,CC=205,DD=187,EE=200,FF=188,GG=186;
    Qa=BB;Qb=CC;Qc=DD;Qd=EE;Qe=FF;Qf=GG;
   	printf("     %c",Qa);
	for (i=1;i<77;i++)
    printf("%c",Qb);
    printf("%c\n",Qc); 
	printf("     %c              Spring 1             |              Spring 2                  %c\n",Qf,Qf);
	printf("     %c        Strain          Stress     |         Strain         Stress          %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<77;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<n;i++)
    printf("          %e     %e     %e     %e\n",STRAIN01[i],STRESS01[i],STRAIN02[i],STRESS02[i]);	
}
double Horizental_Disp(double L,double &x,double y){
	int it,itermax;
            double residual,tolerance,dx,dx_ABS,f,df;
			it = 0; // initialize iteration count
		    itermax = 100000;
		    residual = 100; // initialize residual
		    tolerance = 1e-12;
		    x = 1;// initialize answer
		    while (residual > tolerance){    
		    	f = x*x + y*y - 2*L*x;
		    	df = 2*x - 2*L;
		    	dx = f/df;
		        x -= dx;
                residual = ABS(dx); // abs residual
		        it = it + 1; // increment iteration count
		        //printf("f: %f -\tdx: %f -\tresidual: %f\n",f,dx,residual);
		         if (it == itermax){
		          printf("\tSQRT2(number,power) : SQRT2(%f) - iteration: %d ->   ## The solution is not converged ##\n",y,it);
		          break;
				}
		    }
		        if (it < itermax){
			   //printf("\tSQRT(number,power) - SQRT(%f,%f) : %f \n",D,n, x);
			   return x;
			   }
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
		        x -= dx;
                residual = ABS(dx); // abs residual
		        it += 1; // increment iteration count
		        //printf("f: %f -\tdx: %f -\tresidual: %f\n",f,dx,residual);
		         if (it == itermax){
		          printf("\tSQRT2(number,power) : SQRT2(%f) - iteration: %d ->   ## The solution is not converged ##\n",D,it);
		          break;
				}
		    }
		        if (it < itermax){
			   //printf("\tSQRT(number,power) - SQRT(%f,%f) : %f \n",D,n, x);
			   return x;
			   }
}
