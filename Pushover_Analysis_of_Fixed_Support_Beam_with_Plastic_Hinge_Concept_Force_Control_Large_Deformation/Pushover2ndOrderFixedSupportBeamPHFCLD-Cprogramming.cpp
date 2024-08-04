#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define N 10000
#define NN 10 // Degree of freedom
#define Ne 3    // number of element
#define ShowText01 "Pushover2ndOrderFixedSupportBeamPHFCLD-inputDATA.csv"
#define ShowText02 "Pushover2ndOrderFixedSupportBeamPHFCLD-inputHINGE.csv"
#define ShowText03 "Output data is written in Excel file"
#define ShowText04 "Pushover2ndOrderFixedSupportBeamPHFCLD-outputEXCEL.csv"
#define ShowText05 "Pushover2ndOrderFixedSupportBeamPHFCLD-outputHTML.html"
#define ShowText06 "Graph-outputHTML.html"

void IMPORT_DATA01(double &Length,double &EA,double &EI,double applied_load[],double &Dmax,int &M,int &itermax,double &tolerance);
void IMPORT_DATA02(double TET[],double MOM[],int &k);
void MatrixAssembled(double [Ne][6][6],double [][NN]);
void MatrixDetermination(double [][NN],int );
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [], double [],int );
void MatrixMulti02(double [][NN], double [], double [], double [],int );
void MatrixZero(double A[][NN],int n);
void MatrixChange(double A[][NN],double B[][NN],int n);
void ElementInternalForce(double A[],double B[],double C[],double D[],double E[],double lanX[],double lanY[],double U[],double ee[][6],int I);// Calculate internal element force
void ElementStiffness(double eleF[3][6],double Length,double EI,double EA,double u[],double L[],double lanX[],double lanY[],double AA[],double BB[],double CC[],double DD[],double EE[],double Kele[3][6][6]);
void PlasticHingeStiffnessCOFF(double [],double [],double [],int );// Calculate slope Moment rotation of plastic hinge
void PlasticHingeStiffness(double A[3][6],int I,int II,double B[],double C[],int III,double D[],double E[],double KRK[],int I4,int n);
double ABS(double);
double MAX_ABS(double [],int );
double SQRT2(double);
void MessageNotConverge(int ,int );
void MessageConverge(int ,int ,double ,double []);
void MessageInitialData(double ,double ,double ,double [],double ,int ,double ,int );
void MessageCheckInputMk(int ,int );
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessagePlasticHingeTEXT(double [],double [],int );
void MessageCheck_IMPORT_DATA01(double ,double ,double ,double ,int ,double ,int );
void MessageCheck_IMPORT_DATA02(double [],double [],int );
int MessageControl(double Dmax,double u[],double TET[],int n);
void bilinear(double [][3],double [][2],double [],int );
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_base01[],double output_base02[],double output_base03[],double output_base04[],double output_base05[],int n);
void OUTPUT_html(double applied_load[],double L,double EI,double EA,double Dmax,int itermax,double tolerance,int M,double TET[],double MOM[],int m,double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_base01[],double output_base02[],double output_base03[],double output_base04[],double output_base05[],int n);
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
void ANALYSIS(double Length,double EI,double EA,double Fini,double Dmax,double itermax,double tolerance,double residual,double applied_load[],double TET[],double MOM[],int n,double Rk[],int M);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
int main(){
    double Length,EI,EA,Fini,Dmax,tolerance,residual,applied_load[NN];
    int M,Y,itermax;
    double TET[10],MOM[10],Rk[10];

    IMPORT_DATA01(Length,EA,EI,applied_load,Dmax,M,itermax,tolerance);
    MessageCheck_IMPORT_DATA01(Length,EI,EA,Dmax,itermax,tolerance,M);
	IMPORT_DATA02(TET,MOM,Y);
	MessageCheck_IMPORT_DATA02(TET,MOM,Y);
    MessageCheckInputMk(Y,M);
	PlasticHingeStiffnessCOFF(TET,MOM,Rk,Y);// Calculate slope Moment rotation of plastic hinge
	textcolor(14);
    MessageInitialData(Length,EI,EA,applied_load,Dmax,itermax,tolerance,M);
    MessagePlasticHingeTEXT(TET,MOM,Y);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(Length,EI,EA,Fini,Dmax,itermax,tolerance,residual,applied_load,TET,MOM,Y,Rk,M);
    getch();
    return 0;
}
void MatrixInverse(double A[][NN], double C[][NN],int n){
  int i,j,l;
  double c_A[n][n],B[n][n],m,Sum;
  for (i=0;i<n;i++)
  for (j=0;j<n;j++)
  c_A[i][j]=A[i][j];
// Inverse [Kinit]
	        for (i=0;i<n;i++)
	        for (j=0;j<n;j++){
				if (i==j)
					B[i][j]=1;
				else
					B[i][j]=0;
			}

	for (j=0;j<n-1;j++)
       	for (i=j+1;i<n;i++){
			m=c_A[i][j]/c_A[j][j];
			for (l=0;l<n;l++){
				c_A[i][l] -= m*c_A[j][l];
				B[i][l] -= m*B[j][l];
			}
		}
			// backward substitutions
   	for (i=n-1;i>=0;i--)
		for (j=0;j<n;j++){
       		Sum=0;
					for (l=i+1;l<n;l++)
		            Sum += c_A[i][l]*C[l][j];
			    C[i][j]=(B[i][j]-Sum)/c_A[i][i];
	    }
}
void PlasticHingeStiffnessCOFF(double A[],double B[],double C[],int n){
	int i;
    C[0]=B[0]/A[0];
    printf("%d - Rk[%d]: %f\n",0,0,C[0]);
    for (i=1;i<n;i++)
    {
    C[i]=(B[i]-B[i-1])/(A[i]-A[i-1]);
    printf("%d - Rk[%d]: %f\n",i,i,C[i]);
	}
}
void PlasticHingeStiffness(double A[3][6],int I,int II,double B[],double C[],int III,double D[],double E[],double KRK[],int I4,int n){
			if (ABS(A[I][II]) >= 0 && ABS(A[I][II]) <= B[0])
		     KRK[I4] = D[0]*10e9;
		     for (int i=0;i<n-1;i++){
		     if  (ABS(A[I][II]) > B[i] && ABS(A[I][II]) <= B[i+1])
		     KRK[I4] = (B[i]+D[i+1]*(ABS(C[III])-E[i]))/ABS(C[III]);
			 }
		    if  (ABS(A[I][II]) > B[n-1])
		     KRK[I4] = 0.0;
		     //printf("\t\t\t %f\n",KRK[I4]);
}
void ElementInternalForce(double A[],double B[],double C[],double D[],double E[],double lanX[],double lanY[],double U[],double ee[][6],int I){
    double K[6][6],lan[6][6],UU[6],ff,ll[6][6];
    K[0][0]=E[I];K[0][1]=0;K[0][2]=0;K[0][3]=-E[I];K[0][4]=0;K[0][5]=0;
	K[1][0]=0;K[1][1]=D[I];K[1][2]=B[I];K[1][3]=0;K[1][4]=-D[I];K[1][5]=B[I];
	K[2][0]=0;K[2][1]=B[I];K[2][2]=A[I];K[2][3]=0;K[2][4]=-B[I];K[2][5]=C[I];
	K[3][0]=-E[I];K[3][1]=0;K[3][2]=0;K[3][3]=E[I];K[3][4]=0;K[3][5]=0;
	K[4][0]=0;K[4][1]=-D[I];K[4][2]=-B[I];K[4][3]=0;K[4][4]=D[I];K[4][5]=-B[I];
	K[5][0]=0;K[5][1]=B[I];K[5][2]=C[I];K[5][3]=0;K[5][4]=-B[I];K[5][5]=A[I];

	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;

	if (I == 0){
	UU[0]=U[9];UU[1]=0;UU[2]=U[0];UU[3]=U[1];UU[4]=U[2];UU[5]=U[3];
	}
	if (I == 1){
	UU[0]=U[1];UU[1]=U[2];UU[2]=U[3];UU[3]=U[4];UU[4]=U[5];UU[5]=U[6];
	}
	if (I == 2){
	UU[0]=U[4];UU[1]=U[5];UU[2]=U[6];UU[3]=0;UU[4]=U[7];UU[5]=U[8];
	}
	int i,j,k;
	// [f] = [K] * [u]
	for (i=0; i<6; i++){
		for (j=0; j<6; j++){
	    ff = 0;
	    for (k=0; k<6; k++)
		ff += K[i][k]*lan[k][j];
		ll[i][j] = ff;
		}

    }
    for (i=0; i<6; i++){
	    ff=0;
			    for (j=0; j<6; j++)
			    ff += ll[i][j]*UU[j];
		ee[I][i] = ff;
	}
}
void MatrixMulti01(double A[][NN], double B[], double C[], double D[],int n){
	int i,j;
	double ff;
    // [f] = [Ktot] * [u] - [F]
    for (i=0; i<n; i++){
  	ff=0;
    for (j=0; j<n; j++)
    ff += A[i][j]*B[j];

    D[i] = ff - C[i];
    }
}
void MatrixMulti02(double A[][NN], double B[], double C[], double D[],int n){
	int i,j;
	double Dx;
	// [du] = [InvKinit] * [f]
    for (i=0; i<n; i++){
    Dx=0;
    for (j=0; j<n; j++)
    Dx += A[i][j]* -B[j];

    C[i]= Dx;
    D[i] += C[i];// u= u+du
	}
}
double ABS(double B){
	if (B < 0)
    B = -B;//Absolute number
    else
    B = B;
    return B;
}
void Distance(int i){
	        if (i < 10)
		    printf("\b\t");
	        if (i >= 10 && i <= 99)
		    printf("\b\t\b");
		    if (i >= 100 && i <= 999)
		    printf("\b\t\b\b");
		    if (i >= 1000 && i <= 9999)
		    printf("\b\t\b\b\b");
		    if (i >= 10000 && i <= 20000)
		    printf("\b\t\b\b\b\b");
}
void MessageNotConverge(int ii,int iit){
	    Distance(ii+1);
	    printf("         %d\t\t%d ->   ## The solution for this step is not converged ##\n",ii+1,iit);

}
void MessageConverge(int ii,int iit,double FP,double A[]){
		Distance(ii+1);
		printf("        %d\t\t%d\t     %.3e\t\t    %.3e\t\t      %.3e\t\t        %.3e\n",ii+1,iit,FP,A[2],A[5],A[7]);
}
void MessageInitialData(double L,double EI,double EA,double applied_load[],double Dmax,int itermax,double tolerance,int M){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<93;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c                                  >> IN THE NAME OF ALLAH <<                                %c\n",Qf,Qf);
    printf("\t\t\t\t%c Nonlinear 2nd-Order Analysis of 2D Fixed Support Beam with Force Control Large Deformation %c\n",Qf,Qf);
    printf("\t\t\t\t%c                                      UNIT: Free Unit                                       %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<93;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
	printf("\t\t\t\t%c                      This program is written by Salar Delavar Ghashghaei                   %c\n",Qf,Qf);
	printf("\t\t\t\t%c                              E-mail: salar.d.ghashghaei@gmail.com                          %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<93;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);


    MessageInputDataTEXT();
    printf("      Length :                                                                 %.3e\n",L);
    printf("      Section flextural rigidity - EI:                                         %.3e\n",EI);
    printf("      Section axial rigidity - EA:                                             %.3e\n",EA);
    printf("      Initial incremental external axial force [DOF(1)]:                       %.3e\n",applied_load[0]);
    printf("      External force [DOF(3)]:                                                 %.3e\n",applied_load[1]);
    printf("      External force [DOF(4)]:                                                 %.3e\n",applied_load[2]);
	printf("      External force [DOF(5)]:                                                 %.3e\n",applied_load[3]);
	printf("      External force [DOF(6)]:                                                 %.3e\n",applied_load[4]);
	printf("      External force [DOF(7)]:                                                 %.3e\n",applied_load[5]);
	printf("      External force [DOF(8)]:                                                 %.3e\n",applied_load[6]);
	printf("      External force [DOF(9)]:                                                 %.3e\n",applied_load[7]);
	printf("      External force [DOF(11)]:                                                %.3e\n",applied_load[8]);
	printf("      External force [DOF(12)]:                                                %.3e\n",applied_load[9]);
	printf("      Ultimate displacement [DOF(11)]:                                         %.3e\n",Dmax);// maximum number of iterations
    printf("      Maximum number of iterations:                                            %d\n",(int)itermax);
	printf("      Specified tolerance for convergence:                                     %.3e\n",tolerance);//  specified tolerance for convergence
    printf("      Number of calculation:                                                   %d\n",M);
}

void MessageAnalysisReport(){
  int i;
  char Ql=176;
  printf("\n     ");
  for (i=1;i<64;i++)
  printf("%c",Ql);
  printf(" Analysis Report ");
  for (i=1;i<64;i++)
  printf("%c",Ql);
  printf("\n");
	printf("\t   ---------------------------------------------------------------------------------------------------------------------------------\n");
    printf("\t     Increment    Iteration   Incremental load [DOF(11)]   Displacement [DOF(5)]   Displacement [DOF(8)]   Displacement [DOF(11)]   \n");
	printf("\t   ---------------------------------------------------------------------------------------------------------------------------------\n");
}
void MessagePlasticHingeTEXT(double TET[],double MOM[],int Y){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188;Qf=186;
   	printf("     %c",Qa);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
	printf("     %c      Plastic Hinge Data        %c\n",Qf,Qf);
	printf("     %c     Rotation        Moment     %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<Y;i++)
	printf("          %.3e     %.3e\n",TET[i],MOM[i]);
}
void MessageCheck_IMPORT_DATA01(double L,double EI,double EA,double Dmax,int itermax,double tolerance,int M){
	if ( L <= 0 ||  EI <= 0 ||    EA <= 0 || Dmax <= 0 || itermax <= 0 || tolerance<= 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative data or zero input value is not acceptable ***\n");
	printf("  Length of element:               %.3e\n",L);
	printf("  Section flextural rigidity - EI: %.3e\n",EI);
    printf("  Section axial rigidity - EA:     %.3e\n",EA);
	printf("  Ultimate displacement [DOF(11)]: %.3e\n",Dmax);
	printf("  Maximum iteration:               %d\n",itermax);
	printf("  Tolerance:                       %.3e\n",tolerance);
	Sleep(40000);
	exit(1);
	 }
}
void MessageCheck_IMPORT_DATA02(double A[],double B[],int n){
	int i;
	for(i=0;i<n;i++){
		if (A[i] < 0|| B[i] < 0){
	    MessageErrorReportTEXT();
        printf("               Please check this file! -> [%s]\n",ShowText02);
		printf("           Row %d has a negative value.\n",i+1);
		printf("                            *** Negative data input value is not acceptable ***\n");
		Sleep(40000);
		exit(1);
		}
	}
	for (i=1;i<n;i++){
    if (A[i] <= A[i-1]){
    MessageErrorReportTEXT();
	printf("          Please check the input file! -> [%s]\n",ShowText02);
	printf("          Check row %d value.\n",i+1);
	printf("          Rotation[%d]: %.3e\n",i+1,A[i]);
	printf("                            *** Data must be sort from minimum value to maximum value ***\n");
	Sleep(40000);exit(1);
	}
	}
}
void MessageCheckInputMk(int Y,int M){
	if (Y>5 || Y < 5 || M>N || M<2){
	MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
    printf("        Plastic hinge data: %d - Plastic hinge data must be data : 5\n",Y);
	printf("        Load increments:    %d - Minimum : 3 - Maximum : 10000\n",M);
	Sleep(40000);
	exit(1);
	}
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
int MessageControl(double Dmax,double u[],double TET[],int n){
        int i;
		if (ABS(u[7]) >= Dmax){
		textcolor(13);
		printf("\n                  ## Displacement [DOF(11)] reached to ultimate displacement ##\n\n");
		i = 1;
	    }
	    if (ABS(u[0]) >= TET[n-1]){
	    textcolor(13);
		printf("\n                  ## Rotation [DOF(3)] reached to ultimate rotation ##\n\n");
		i = 1;
	    }
	    if (ABS(u[8]) >= TET[n-1]){
	    textcolor(13);
		printf("\n                  ## Rotation [DOF(12)] reached to ultimate rotation ##\n\n");
		i = 1;
	    }
		return i;
}
void bilinear(double A[][3],double B[][2],double C[],int M){
	double AaSUM=0,B_max=0,k0;
	int I,Mmax;
	double *hh = new double [N];
	double *Aa = new double [N];
	double *AA = new double [N];
	double *BB = new double [N];

    for (I=0;I<M;I++){
	if (ABS(B[0][0]) < ABS(B[I][0])){
	  B[0][0] = ABS(B[I][0]);
	  Mmax = I;
	  }//if
      }//for
       Mmax = I;
    AA[0]=0;BB[0]=0;
    for (I=0;I<Mmax;I++){
	AA[I+1]=ABS(A[I][1]);
    BB[I+1]=ABS(B[I][0]);
	}
    for (I=0;I<Mmax-1;I++){
	hh[I] = AA[I+1]-AA[I];
	Aa[I]=(BB[I]+BB[I+1])*0.5*hh[I];
	AaSUM=AaSUM+Aa[I];
    }
    C[1]=0;// C[0]=D_y ; C[1]=D_u ; C[2]=F_y ; C[3]=F_u ; C[4]=SC ; C[5]=OF ;
	for (I=0;I<Mmax+2;I++){ // find max
	if(C[1] < AA[I])
	C[1]=AA[I];
	if(B_max < BB[I])
	B_max=BB[I];
	}
	//double MOMu;
    C[3]=BB[I-2];
	k0 =BB[4]/AA[4];
    C[0] = (C[3]*C[1]*0.5-AaSUM)/(C[3]*0.5 - k0*C[1]*0.5);
    C[2] = k0*C[0];
    C[4]=C[1]/C[0];C[5]=C[3]/C[2];
}
void OUTPUT_txt(){
}
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_base01[],double output_base02[],double output_base03[],double output_base04[],double output_base05[],int n){
// EXCEL OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText04, "w");
fprintf(OutputFile," ### Nonlinear 2nd-Order Analysis of 2D Fixed Support Beam with Force Control Large Deformation ###\n");
fprintf(OutputFile,"Increment,Axial load[DOF(1)],Base Shear[DOF(2)],Ele. Moment[DOF(6)],Ele. Moment[DOF(9)],Base Moment[DOF(12)],Displacement [DOF(1)],Displacement [DOF(5)],Displacement [DOF(7)],Displacement [DOF(11)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i+1,output_base01[i],output_base02[i],output_base03[i],output_base04[i],output_base05[i],output_u01[i],output_u02[i],output_u03[i],output_u04[i]);
fclose(OutputFile);
}
void OUTPUT_html(double applied_load[],double Length,double EI,double EA,double Dmax,int itermax,double tolerance,int M,double TET[],double MOM[],int m,double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_base01[],double output_base02[],double output_base03[],double output_base04[],double output_base05[],int n){
// HTML OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText05, "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"Pushover2ndOrderFixedSupportBeamPHFCLD-image01.jpg\" style=\"width:1000px ; height:500px\" alt=\"analysis\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th bgcolor=\"cyan\"> Nonlinear 2nd-Order Analysis of 2D Fixed Support Beam with Force Control Large Deformation - Output Report </th> \n");
    // TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Input Data </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Length: </th><th> %.3e </th> </tr>\n",Length);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section flextural rigidity - EI: </th><th> %.3e </th> </tr>\n",EI);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section axial rigidity - EA: </th><th> %.3e </th> </tr>\n",EA);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Initial incremental external axial [DOF(1)]: </th><th> %.3e </th> </tr>\n",applied_load[0]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(3)]: </th><th> %.3e </th> </tr>\n",applied_load[1]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(4)]: </th><th> %.3e </th> </tr>\n",applied_load[2]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(5)]: </th><th> %.3e </th> </tr>\n",applied_load[3]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(6)]: </th><th> %.3e </th> </tr>\n",applied_load[4]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(7)]: </th><th> %.3e </th> </tr>\n",applied_load[5]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(8)]: </th><th> %.3e </th> </tr>\n",applied_load[6]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(9)]: </th><th> %.3e </th> </tr>\n",applied_load[7]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(11)]: </th><th> %.3e </th> </tr>\n",applied_load[8]);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External force [DOF(12)]: </th><th> %.3e </th> </tr>\n",applied_load[9]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Ultimate displacement [DOF(11)]: </th><th> %.3e </th> </tr>\n",Dmax);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Maximum number of iterations: </th><th> %d </th> </tr>\n",itermax);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Specified tolerance for convergence: </th><th> %.3e </th> </tr>\n",tolerance);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Number of calculation: </th><th> %d </th> </tr>\n",M);

    // TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th colspan=\"4\" bgcolor=\"orange\"> Hinges Data </th> \n");
	fprintf(OutputFile,"<tr> <th colspan=\"2\" bgcolor=\"orange\">Hinge: Moment - Rotation </th></tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Rotation  </th> <th bgcolor=\"orange\"> Moment </th> </tr>\n");
	for(i=0;i<m;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> </tr>\n",TET[i],MOM[i]);
    }
	// TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1600px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"10\" bgcolor=\"orange\"> Structral Deformation </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Axial load[DOF(1)]</th><th bgcolor=\"orange\">Base Shear[DOF(2)]</th><th bgcolor=\"orange\">Ele. Moment[DOF(6)]</th> <th bgcolor=\"orange\">Ele. Moment[DOF(9)]</th><th bgcolor=\"orange\">Ele. Moment[DOF(12)]</th><th bgcolor=\"orange\">Displacement [DOF(1)]</th><th bgcolor=\"orange\">Displacement [DOF(5)]</th><th bgcolor=\"orange\">Displacement [DOF(8)]</th><th bgcolor=\"orange\">Displacement [DOF(11)]</th></tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td></tr>\n",i+1,output_base01[i],output_base02[i],output_base03[i],output_base04[i],output_base05[i],output_u01[i],output_u02[i],output_u03[i],output_u04[i]);
    }
	fprintf(OutputFile,"</table></body></html>\n");
	fclose(OutputFile);
}
void MatrixAssembled(double Kele[Ne][6][6],double K[][NN]){
    K[0][0]= Kele[0][2][1];//DOF(3)
    K[0][1]= Kele[0][2][2];//DOF(3)
    K[0][2]= Kele[0][2][3];//DOF(3)
    K[0][3]= Kele[0][2][4];//DOF(3)
    K[0][4]= 0;
    K[0][5]= 0;
    K[0][6]= 0;
    K[0][7]= 0;
    K[0][8]= 0;
    K[0][9]= Kele[0][2][0];//DOF(3)

    K[1][1]= Kele[0][3][3]+Kele[1][0][0];//DOF(4)
    K[1][2]= Kele[0][3][4]+Kele[1][0][1];//DOF(4)
    K[1][3]= Kele[0][3][5]+Kele[1][0][2];//DOF(4)
    K[1][4]= Kele[1][0][3];//DOF(4)
    K[1][5]= Kele[1][0][4];//DOF(4)
    K[1][6]= Kele[1][0][5];//DOF(4)
    K[1][7]= 0;
    K[1][8]= 0;
    K[1][9]= Kele[0][3][0];//DOF(4)

    K[2][2]= Kele[0][4][4]+Kele[1][1][1];//DOF(5)
    K[2][3]= Kele[0][4][5]+Kele[1][1][2];//DOF(5)
    K[2][4]= Kele[1][1][3];//DOF(5)
    K[2][5]= Kele[1][1][4];//DOF(5)
    K[2][6]= Kele[1][1][5];//DOF(5)
    K[2][7]= 0;
    K[2][8]= 0;
    K[2][9]= Kele[0][4][0];//DOF(5)

    K[3][3]= Kele[0][5][5]+Kele[1][2][2];//DOF(6)
    K[3][4]= Kele[1][2][3];//DOF(6)
    K[3][5]= Kele[1][2][4];//DOF(6)
    K[3][6]= Kele[1][2][5];//DOF(6)
    K[3][7]= 0;
    K[3][8]= 0;
    K[3][9]= Kele[0][5][0];//DOF(6)

    K[4][4]= Kele[1][3][3]+Kele[2][0][0];//DOF(7)
    K[4][5]= Kele[1][3][4]+Kele[2][0][1];//DOF(7)
    K[4][6]= Kele[1][3][5]+Kele[2][0][2];//DOF(7)
    K[4][7]= Kele[2][0][4];//DOF(7)
    K[4][8]= Kele[2][0][5];//DOF(7)
    K[4][9]= 0;

    K[5][5]= Kele[1][4][4]+Kele[2][1][1];//DOF(8)
    K[5][6]= Kele[1][4][5]+Kele[2][1][2];//DOF(8)
    K[5][7]= Kele[2][1][4];//DOF(8)
    K[5][8]= Kele[2][1][5];//DOF(8)
    K[5][9]= 0;

    K[6][6]= Kele[1][5][5]+Kele[2][2][2];//DOF(9)
    K[6][7]= Kele[2][2][4];//DOF(9)
    K[6][8]= Kele[2][2][5];//DOF(9)
    K[6][9]= 0;

    K[7][7]= Kele[2][4][4];//DOF(11)
    K[7][8]= Kele[2][4][5];//DOF(11)
    K[7][9]= 0;

    K[8][8]= Kele[2][5][5];//DOF(12)
    K[8][9]= 0;

    K[9][9]= Kele[0][0][0];//DOF(1)

    K[1][0]= K[0][1];
    K[2][0]= K[0][2];
    K[3][0]= K[0][3];
    K[4][0]= K[0][4];
    K[5][0]= K[0][5];
    K[6][0]= K[0][6];
    K[7][0]= K[0][7];
    K[8][0]= K[0][8];
    K[9][0]= K[0][9];

    K[2][1]= K[1][2];
    K[3][1]= K[1][3];
    K[4][1]= K[1][4];
    K[5][1]= K[1][5];
    K[6][1]= K[1][6];
    K[7][1]= K[1][7];
    K[8][1]= K[1][8];
    K[9][1]= K[1][9];

    K[3][2]= K[2][3];
    K[4][2]= K[2][4];
    K[5][2]= K[2][5];
    K[6][2]= K[2][6];
    K[7][2]= K[2][7];
    K[8][2]= K[2][8];
    K[9][2]= K[2][9];

    K[4][3]= K[3][4];
    K[5][3]= K[3][5];
    K[6][3]= K[3][6];
    K[7][3]= K[3][7];
    K[8][3]= K[3][8];
    K[9][3]= K[3][9];

    K[5][4]= K[4][5];
    K[6][4]= K[4][6];
    K[7][4]= K[4][7];
    K[8][4]= K[4][8];
    K[9][4]= K[4][9];

    K[6][5]= K[5][6];
    K[7][5]= K[5][7];
    K[8][5]= K[5][8];
    K[9][5]= K[5][9];

    K[7][6]= K[6][7];
    K[8][6]= K[6][8];
    K[9][6]= K[6][9];

    K[8][7]= K[7][8];
    K[9][7]= K[7][9];

	K[9][8]= K[8][9];
}
void ElementStiffness(double eleF[3][6],double Length,double EI,double EA,double u[],double L[],double lanX[],double lanY[],double AA[],double BB[],double CC[],double DD[],double EE[],double Kele[3][6][6]){
	double x1,y1,x2,y2,x3,y3,x4,y4,P;
    x1 = u[9];y1=0;
	x2 = Length/3 + u[1];y2=u[2];
	x3 = 2*Length/3 + u[4];y2=u[5];
	x4 = Length;y4=u[7];
    L[0] = SQRT2((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    L[1] = SQRT2((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
    L[2] = SQRT2((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3));
    lanX[0]=(x2-x1)/L[0];lanY[0]=(y2-y1)/L[0];
  	lanX[1]=(x3-x2)/L[1];lanY[1]=(y3-y2)/L[1];
  	lanX[2]=(x4-x3)/L[2];lanY[2]=(y4-y3)/L[2];
  	for (int i=0;i<Ne;i++){
  	P = -.5*(eleF[i][3] - eleF[i][0]);
    AA[i] = 4*EI/L[i] - (2*P*L[i])/15;
    BB[i] = 6*EI/(L[i]*L[i]) - (P/10);
    CC[i] = 2*EI/L[i] + (P*L[i])/30;
    DD[i] = 12*EI/(L[i]*L[i]*L[i]) - (6*P)/(5*L[i]);
    EE[i] = EA/L[i];
    Kele[i][0][0]=EE[i]*lanX[i]*lanX[i]+DD[i]*lanY[i]*lanY[i];
    Kele[i][0][1]=(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][0][2]=-BB[i]*lanY[i];
    Kele[i][0][3]=-(EE[i]*lanX[i]*lanX[i]+DD[i]*lanY[i]*lanY[i]);
    Kele[i][0][4]=-(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][0][5]=-BB[i]*lanY[i];

    Kele[i][1][0]=(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][1][1]=EE[i]*lanY[i]*lanY[i]+DD[i]*lanX[i]*lanX[i];
    Kele[i][1][2]=BB[i]*lanX[i];
    Kele[i][1][3]=-(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][1][4]=-(EE[i]*lanY[i]*lanY[i]+DD[i]*lanX[i]*lanX[i]);
    Kele[i][1][5]=BB[i]*lanX[i];

    Kele[i][2][0]=-BB[i]*lanY[i];
    Kele[i][2][1]=BB[i]*lanX[i];
    Kele[i][2][2]=AA[i];
    Kele[i][2][3]=BB[i]*lanY[i];
    Kele[i][2][4]=-BB[i]*lanX[i];
    Kele[i][2][5]=CC[i];

    Kele[i][3][0]=-(EE[i]*lanX[i]*lanX[i]+DD[i]*lanY[i]*lanY[i]);
    Kele[i][3][1]=-(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][3][2]=BB[i]*lanY[i];
    Kele[i][3][3]=EE[i]*lanX[i]*lanX[i]+DD[i]*lanY[i]*lanY[i];
    Kele[i][3][4]=(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][3][5]=BB[i]*lanY[i];

    Kele[i][4][0]=-(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][4][1]=EE[i]*lanY[i]*lanY[i]+DD[i]*lanX[i]*lanX[i];
    Kele[i][4][2]=-BB[i]*lanX[i];
    Kele[i][4][3]=(EE[i]-DD[i])*lanX[i]*lanY[i];
    Kele[i][4][4]=EE[i]*lanY[i]*lanY[i]+DD[i]*lanX[i]*lanX[i];
    Kele[i][4][5]=-BB[i]*lanX[i];

    Kele[i][5][0]=-BB[i]*lanY[i];
    Kele[i][5][1]=BB[i]*lanX[i];
    Kele[i][5][2]=CC[i];
    Kele[i][5][3]=BB[i]*lanY[i];
    Kele[i][5][4]=-BB[i]*lanX[i];
    Kele[i][5][5]=AA[i];
	}
	}
double SQRT2(double D){
            int it,itermax;
            double residual,tolerance,x,dx,dx_ABS,f,df;
			it = 0; // initialize iteration count
		    itermax = 100000;
		    residual = 100; // initialize residual
		    tolerance = 1e-9;
		    x = 1;// initialize answer
		    while (residual > tolerance){
		    	f = x*x - D;
		    	df = 2 * x;
		    	dx = f/df;
		        x= x - dx;
                residual = ABS(dx); // abs residual
		        it = it + 1; // increment iteration count
		        //printf("f: %f -\tdx: %f -\tresidual: %f\n",f,dx,residual);
		         if (it == itermax)
		        {
		          //printf("\tSQRT2(number,power) : SQRT2(%f) - iteration: %d ->   ## The solution is not converged ##\n",D,it);
		          break;
				}
		    }
		        if (it < itermax)
               {
			   //printf("\tSQRT(number,power) - SQRT(%f,%f) : %f \n",D,n, x);
			   return x;
			   }
}
void MatrixDetermination(double A[][NN],int n){
	// row operations
	int i,j,k;
   	double Product,m,B[n][n];
   	for (i=0;i<n;i++)
   	for (j=0;j<n;j++)
   	B[i][j]=A[i][j];

   	for (k=0;k<n-1;k++)
       	for (i=k+1;i<n;i++){
			m = B[i][k]/B[k][k];
			for (j=0;j<n;j++)
				B[i][j] -= m*B[k][j];
		}
		Product=1;
	for (i=0;i<n;i++)
		Product *= B[i][i];

	// display results
	if (Product == 0){
	printf("\a\n\t ### it Seens that Golobal Matrix is singular or structure is unstable!!! ###\n");
	Sleep(40000);
	exit(1);
	}

}
void IMPORT_DATA01(double &Length,double &EA,double &EI,double applied_load[],double &Dmax,int &M,int &itermax,double &tolerance){
    double Import_Data[17];
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
	Length=Import_Data[0];
	EI=Import_Data[1];
	EA=Import_Data[2];
	applied_load[0]=Import_Data[3];
	applied_load[1]=Import_Data[4];
	applied_load[2]=Import_Data[5];
	applied_load[3]=Import_Data[6];
	applied_load[4]=Import_Data[7];
	applied_load[5]=Import_Data[8];
	applied_load[6]=Import_Data[9];
	applied_load[7]=Import_Data[10];
	applied_load[8]=Import_Data[11];
	applied_load[9]=Import_Data[12];
	Dmax=Import_Data[13];
	M=Import_Data[14];
	itermax=Import_Data[15];
	tolerance=Import_Data[16];
}
void IMPORT_DATA02(double TET[],double MOM[],int &k){
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
	fscanf(InputFile,"%lf,%lf",&TET[i],&MOM[i]);
	printf("\t TET[%d]: %f - MOM[%d]: %f \n",i,TET[i],i,MOM[i]);
	i++;
	}
	while(i < N && fgets(line,sizeof(line),InputFile) != NULL);
	k = i-1;
    //printf("%d\n",k);
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
void ANALYSIS(double Length,double EI,double EA,double Fini,double Dmax,double itermax,double tolerance,double residual,double applied_load[],double TET[],double MOM[],int n,double Rk[],int M){
   int i,j,z,zMAX,it,well_done=0;
   double Krk[2],Krk01,Krk02,eleF[Ne][6],K[NN][NN],Kinit[NN][NN],Kt[NN][NN],InvKinit[NN][NN];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],Kele[Ne][6][6],fp,m,Sum,Dx,ff,f[NN],F[NN],Fi[NN],u[NN],du[NN];

   double *output_u01 = new double[N];
   double *output_u02 = new double[N];
   double *output_u03 = new double[N];
   double *output_u04 = new double[N];
   double *output_base01 = new double[N];
   double *output_base02 = new double[N];
   double *output_base03 = new double[N];
   double *output_base04 = new double[N];
   double *output_base05 = new double[N];
   double *X = new double[N];
   double *Y = new double[N];
//   double **K;
//   K = new double*[NN];
//   for (i=0;i<NN;i++)
//   K[i] = new double [NN];
    for (i=0;i<NN;i++)
    u[i] = 0.0;
    for (i=0;i<3;i++)
    for (j=0;j<6;j++)
    eleF[i][j]=0.0;

    Fini=applied_load[0];
    for (i=0;i<NN-1;i++)
    F[i] = applied_load[i];

    for (z=0;z<M;z++){
    fp = Fini*(z+1);// Define the applied load
    F[9] = fp;

    ElementStiffness(eleF,Length,EI,EA,u,L,lanX,lanY,AA,BB,CC,DD,EE,Kele);

    for (i=0;i<NN;i++)
    f[i] = 0;

		     // Plastic hinge DOF(3)
		    PlasticHingeStiffness(eleF,0,2,MOM,u,1,Rk,TET,Krk,0,n);
            // Plastic hinge DOF(12)
            PlasticHingeStiffness(eleF,2,5,MOM,u,8,Rk,TET,Krk,1,n);

    MatrixAssembled(Kele,K);
    MatrixZero(Kinit,NN);
    MatrixChange(K,Kinit,NN);
    Kinit[0][0] += Krk01;
    Kinit[8][8] += Krk02;


	 // Inverse [Kinit]
	  MatrixInverse(Kinit,InvKinit,NN);

	it = 0; // initialize iteration count
    residual = 100; // initialize residual
	while (residual > tolerance){

	ElementStiffness(eleF,Length,EI,EA,u,L,lanX,lanY,AA,BB,CC,DD,EE,Kele);

		     // Plastic hinge DOF(3)
		    PlasticHingeStiffness(eleF,0,2,MOM,u,1,Rk,TET,Krk,0,n);
            // Plastic hinge DOF(12)
            PlasticHingeStiffness(eleF,2,5,MOM,u,8,Rk,TET,Krk,1,n);

    MatrixZero(Kt,NN);
    MatrixChange(K,Kt,NN);
    Kt[0][0] += Krk01;
    Kt[8][8] += Krk02;

    // Finding the determinant of a square matrix
    MatrixDetermination(Kt,NN);
	// [f] = [Kt] * [u] - [F]
    MatrixMulti01(Kt,u,F,f,NN);
	// [du] = [InvKinit] * [f]
    MatrixMulti02(InvKinit,f,du,u,NN);
    // Max residual
    residual = MAX_ABS(du,NN);
    // increment iteration count
    it = it + 1;
        if (it == itermax){
          MessageNotConverge(z,it);
          break;
		}
	}//while

	// iteration control
    if (it < itermax){
	MessageConverge(z,it,fp,u);
	}
    zMAX = z+1;

		output_u01[z] = u[9];
		output_u02[z] = u[2];
		output_u03[z] = u[5];
		output_u04[z] = u[7];

		ElementInternalForce(AA,BB,CC,DD,EE,lanX,lanY,u,eleF,0);
		output_base01[z] = fp;//output base axial DOF(1)
		output_base02[z] = -eleF[0][1];//output base shear DOF(2)
		output_base03[z] = eleF[0][5];//output moment DOF(6)
		ElementInternalForce(AA,BB,CC,DD,EE,lanX,lanY,u,eleF,1);
		output_base04[z] = eleF[1][5];//output base moment DOF(9)
		ElementInternalForce(AA,BB,CC,DD,EE,lanX,lanY,u,eleF,2);
		output_base05[z] = eleF[2][5];//output base moment DOF(12)
		if (MessageControl(Dmax,u,TET,n) == 1){
		well_done = 1;
		break;
		}


//        if (ABS(u[7]) >= Dmax){
//		printf("\n                  ## Displacement [DOF(11)] reached to ultimate displacement ##\n\n");
//		well_done = 1;
//		break;
//	    }
//	    if (ABS(u[0]) >= TET[n-1]){
//		printf("\n                  ## Rotation [DOF(3)] reached to ultimate rotation ##\n\n");
//		well_done = 1;
//		break;
//	    }
//	    if (ABS(u[8]) >= TET[n-1]){
//		printf("\n                  ## Rotation [DOF(12)] reached to ultimate rotation ##\n\n");
//		well_done = 1;
//		break;
//	    }

    }//for
	    if (well_done == 1){
		OUTPUT_excel(output_u01,output_u02,output_u03,output_u04,output_base01,output_base02,output_base03,output_base04,output_base05,zMAX);
		OUTPUT_html(applied_load,Length,EI,EA,Dmax,itermax,tolerance,M,TET,MOM,n,output_u01,output_u02,output_u03,output_u04,output_base01,output_base02,output_base03,output_base04,output_base05,zMAX);
		for (i=0;i<zMAX;i++){
		X[i] = ABS(output_u04[i]);// displacement DOF (11)
		Y[i] = ABS(output_base01[i]);// axial load DOF(1)
		}
		OUTPUT_HTML_GRAPH(X,Y,zMAX,"Axial load-Displacement Graph","Displacement [DOF(11)]","Axial load [DOF(1)]");
		textcolor(15);
		printf("\n\a - %s -\n",ShowText03);
		system("start /w Graph-outputHTML.html");
		DATE_TIME();
		}
		free(output_u01);free(output_u02);free(output_u03);free(output_u04);
		free(output_base01);free(output_base02);free(output_base03);free(output_base04);free(output_base05);
		free(X);free(Y);
}
void DATE_TIME(){
		printf("\n\t");
		system("echo %date%");
		printf("\t");
		system("echo %time%");
}
double MAX_ABS(double A[],int n){
	int i;
	double B[N];
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
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]){
    // HTML GRAPH OUTPUT
	int i;
	double x,y,Xmax,Ymax;
	double *Xnew = new double [N];
	double *Ynew = new double [N];
	double *NorX = new double [N];
	double *NorY = new double [N];
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
	OutputFile = fopen(ShowText06, "w");
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
	fprintf(OutputFile,"X=x+.25*Lx;Y=.7*y;s4.translate(X,Y);s4.font=\"50px serif\";s4.fillStyle = \"#7fff00\";s4.fillText(\"%s\",0,0); \n",text1);
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
void MatrixZero(double A[][NN],int n){
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
	A[i][j]=0.0;
}
void MatrixChange(double A[][NN],double B[][NN],int n){
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
	B[i][j]=A[i][j];
}
