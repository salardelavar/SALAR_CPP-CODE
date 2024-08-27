#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 6 // Degree of freedom
#define Ne 1    // number of element
#define N 2    // number of node
#define STEP 5000   // number of node
#define ShowText01 "PushoverLinearSemiRigidBerFC-inputDATA.csv"
#define ShowText02 "PushoverLinearSemiRigidBerFC-inputCOORDINATE.csv"
#define ShowText04 "Output data is written in Excel and Html file"
#define ShowText05 "PushoverLinearSemiRigidBerFC-outputEXCEL.csv"
#define ShowText06 "PushoverLinearSemiRigidBerFC-outputHTML.html"
#define ShowText07 "Graph-outputHTML.html"
void IMPORT_DATA01(double &EA,double &EI,double Spring_Stiff[],double &Fini,int &M);
void IMPORT_DATA02(double x[],double y[],int &n);
void MatrixZero(double A[][NN],int n);
void Matrix_Stiffness(double Spring_Stiff[],double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
void MatrixDetermination(double A[][NN],int n,double &Product);
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [],int );
void Matrix_Transpose(double A[6][6],double B[6][6]);
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]);
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I,int II);// Calculate internal element force
void ELEMNT_FORCE_OUTPUT(double eleF[1][6],double ELE_FORCE[6][STEP],int n);
double ABS(double);
double MIN(double A[],int n);
double SQRT2(double D);
void MessageInitialData(double Spring_Stiff[],double EI,double EA,double Fini,double x[],double y[],int M);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double Spring_Stiff[],double EI,double EA,int M);
void MessageCheckMk(int M);
void MessageStrCoorTEXT(double X[],double Y[],int n);
void MessageResult(double output_base01[],double output_base02[],double output_u01[],int n);
int MessageControl(double eleF[Ne][6],double u[],double TET01[],double MOM01[],double TET02[],double MOM02[],int n);
void OUTPUT_excel(double output_u01[],double output_base01[],double output_base02[],double C[6][STEP],int n);
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double output_u01[],double output_base01[],double output_base02[],double C[6][STEP],int n,int M);
void ANALYSIS(double Spring_Stiff[],double EI,double EA,double Fini,double x[],double y[],int M);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
double MAX_ABS(double A[],int n);

int main(){
    double EI,EA,Fini,x[2],y[2],Spring_Stiff[2];
    int n,m,M;
    IMPORT_DATA01(EA,EI,Spring_Stiff,Fini,M);
    IMPORT_DATA02(x,y,n);
    MessageCheck_IMPORT_DATA01(Spring_Stiff,EI,EA,M);
    MessageCheckMk(M);
	textcolor(14);
    MessageInitialData(Spring_Stiff,EI,EA,Fini,x,y,M);
    MessageStrCoorTEXT(x,y,n);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(Spring_Stiff,EI,EA,Fini,x,y,M);
    getch();
    return 0;
}
void MatrixZero(double A[][NN],int n){
	int i,j;
	for (i=0;i<n;i++)
	for (j=0;j<n;j++)
	A[i][j] = 0;
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
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I,int II){
    double lan[6][6],UU[6],ff,ll[6][6];
    int i,j;
	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;

	if (II == 1){
	UU[0]=0;UU[1]=0;UU[2]=U[2];UU[3]=U[0];UU[4]=0;UU[5]=0;
	}

	for (i=0;i<6;i++)
	for (j=0;j<6;j++)
	ll[i][j]=0;
	// [f] = [K] *[lan]* [u]
    Matrix_Multiplication(K,lan,ll);
    for (i=0; i<6; i++){
	    ff=0;
			    for (j=0; j<6; j++)
			    ff += ll[i][j]*UU[j];
		ee[I][i] = ff;
	}
}
void MatrixMulti01(double A[][NN], double B[], double C[],int n){
	int i,j;
	double ff;
    // [u] = [Kinv] * [f]
    for (i=0; i<n; i++)
    {
  	ff=0;
    for (j=0; j<n; j++)
    ff += A[i][j]*B[j];

    C[i] = ff;
    }
}

double ABS(double B){
	if (B < 0)
    B = -B;//Absolute number
    else
    B = B;
    return B;
}
double MIN(double A[],int n){
	int i;
	double Cmin;
	Cmin = A[0];
	// Max of abs
	for (i=0;i<n;i++){
    if(Cmin > A[i])
    Cmin = A[i];
	}
	return Cmin;
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

void MessageInitialData(double Spring_Stiff[],double EI,double EA,double Fini,double x[],double y[],int M){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<82;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c                             >> IN THE NAME OF GOD <<                            %c\n",Qf,Qf);
    printf("\t\t\t\t%c Pushover Semi Rigid Connection Analysis Force Analogy Method with Force Control %c\n",Qf,Qf);
    printf("\t\t\t\t%c                        Based on Euler-Bernoulli Beam Theory                     %c\n",Qf,Qf);
    printf("\t\t\t\t%c                                UNIT: Free Unit                                  %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<82;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
	printf("\t\t\t\t%c               This program is written by Salar Delavar Ghashghaei               %c\n",Qf,Qf);
	printf("\t\t\t\t%c                        E-mail: salar.d.ghashghaei@gmail.com                     %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<82;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

    MessageInputDataTEXT();
    printf("      Section flextural rigidity - EI:                          %.3e\n",EI);
    printf("      Section axial rigidity - EA:                              %.3e\n",EA);
    printf("      Rotational spring stiffness-i:                            %.3e\n",Spring_Stiff[0]);
    printf("      Rotational spring stiffness-j:                            %.3e\n",Spring_Stiff[1]);
    printf("      External Incremental Fx force [DOF(4)]:                   %.3e\n",Fini);
    printf("      Number of increment:                                      %d\n",M);

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
}
void MessageCheck_IMPORT_DATA01(double Spring_Stiff[],double EI,double EA,int M){
	if (EI <= 0 ||    EA <= 0  ||  M <= 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative or zero data input value is not acceptable ***\n");
	printf("\t\t  Section flextural rigidity - EI:          %.3e\n",EI);
    printf("\t\t  Section axial rigidity - EA:              %.3e\n",EA);
    printf("\t\t  Rotational spring stiffness-i:            %.3e\n",Spring_Stiff[0]);
    printf("\t\t  Rotational spring stiffness-j:            %.3e\n",Spring_Stiff[1]);
    printf("\t\t  Number of increment:                      %d\n",M);
	Sleep(40000);
	exit(1);
	 }
}
void MessageCheckMk(int M){
if (M < 2 || M > STEP){
	MessageErrorReportTEXT();
	printf("               Please check this file! -> [%s]\n",ShowText01);	
    printf("        Number of increment: %d - Plastic hinge data must be data -> Minimum : 2 - Maximum : %d\n",M,STEP);
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
void OUTPUT_excel(double output_u01[],double output_base01[],double output_base02[],double C[6][STEP],int n){
// EXCEL OUTPUT
	int i,I;
	FILE *OutputFile;
	OutputFile = fopen(ShowText05, "w");
fprintf(OutputFile," ###  Pushover Semi Rigid Connection Analysis Force Analogy Method with Force Control Based on Euler-Bernoulli Beam Theory ###\n");
fprintf(OutputFile,"\n");
fprintf(OutputFile,"Increment,Base Shear [DOF(1)],Base Moment [DOF(3)],Displacement [DOF(4)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e\n",i+1,output_base01[i],output_base02[i],output_u01[i]);
fprintf(OutputFile,"\n");
fprintf(OutputFile,"Increment,axial-i,shear-i,moment-i,axial-j,shear-j,moment-j\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e\n",i+1,C[0][i],C[1][i],C[2][i],C[3][i],C[4][i],C[5][i]);
fclose(OutputFile);
}

void MatrixDetermination(double A[][NN],int n,double &Product){
	// row operations
	int i,j,k;
   	double m,B[n][n];
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
	//Sleep(40000);
	//exit(1);
	}

}
void IMPORT_DATA01(double &EA,double &EI,double Spring_Stiff[],double &Fini,int &M){
    double Import_Data[6];
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
	while(i < 7 && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
	EI=Import_Data[0];
	EA=Import_Data[1];
	Spring_Stiff[0]=Import_Data[2];
	Spring_Stiff[1]=Import_Data[3];
	Fini=Import_Data[4];
	M=Import_Data[5];
}
void IMPORT_DATA02(double x[],double y[],int &n){
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
	fscanf(InputFile,"%lf,%lf",&x[i],&y[i]);
	i++;
	}
	while(i < 2 && fgets(line,sizeof(line),InputFile) != NULL);
	n = i;
    //printf("%d\n",n);
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
void ANALYSIS(double Spring_Stiff[],double EI,double EA,double Fini,double x[],double y[],int M){
   int i,j,z,zMAX,I;
   double Product,K[NN][NN],Kinv[NN][NN],eleF[Ne][6];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],FF[NN],u[NN];
   double *output_u01 = new double [STEP];
   double *output_base01 = new double [STEP];
   double *output_base02 = new double [STEP];
   double *X = new double [STEP];
   double *Y = new double [STEP];
    double ELE_FORCE[6][STEP];
   double MS[6][6],KG[6][6];
   
    MatrixZero(K,1);
    for (i=0;i<1;i++)
    u[i] = 0;
    for (i=0;i<STEP;i++)
    output_u01[i] = 0.0;	
	

	  for(i=0;i<6*Ne;i++)
	  for(j=0;j<M;j++)
      ELE_FORCE[i][j] = 0;


    for (i=0;i<Ne;i++){
    L[i] = SQRT2((x[i+1]-x[i])*(x[i+1]-x[i])+(y[i+1]-y[i])*(y[i+1]-y[i]));
	lanX[i] = (x[i+1]-x[i])/L[i];lanY[i] = (y[i+1]-y[i])/L[i];
	}
      // STAGE 01
	  for (z=0;z<M;z++){	
	  FF[0]=Fini*(z+1);
      Matrix_Stiffness(Spring_Stiff,EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1); // 0: ele 1 - 1: stiffness matrix: 1

      K[0][0]= KG[3][3];//DOF(4)

	  MatrixInverse(K,Kinv,1);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,1);
      zMAX = z + 1;
      ElementInternalForce(MS,u,lanX,lanY,eleF,0,1);// 1: step 1

      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);
      output_u01[z]=u[0];//output displacement DOF(4)
      output_base01[z]=-ELE_FORCE[1][z];//output base shear DOF(1)
      output_base02[z]=-ELE_FORCE[2][z];//output base moment DOF(3)
      
	  }// for

       MessageResult(output_base01,output_base02,output_u01,zMAX);
	    OUTPUT_excel(output_u01,output_base01,output_base02,ELE_FORCE,zMAX);
	    OUTPUT_html(EI,EA,Fini,x,y,output_u01,output_base01,output_base02,ELE_FORCE,zMAX,M);
	    
		for (i=0;i<zMAX;i++){
	    X[i] = output_u01[i];// Disp. DOF(5)
		Y[i] = output_base01[i];// Base Shear DOF(1)
		}
		OUTPUT_HTML_GRAPH(X,Y,zMAX,"Base Shear-Displacement Graph","Displacement [DOF(5)]","Base Shear [DOF(1)]");
 
		textcolor(15);
		printf("\n\a - %s -\n",ShowText04);
		system("start /w Graph-outputHTML.html");
		DATE_TIME();
		free(output_u01);
		free(output_base01);free(output_base02);
		free(X);free(Y);
}
void DATE_TIME(){
		printf("\n\t");
		system("echo %date%");
		printf("\t");
		system("echo %time%");
}
void Matrix_Stiffness(double Spring_Stiff[],double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II){
	int i,j;
	double lan[6][6],lan_Tr[6][6],ans[6][6],r[2],R;
	for (i=0;i<2;i++)
	r[i] = (Spring_Stiff[i]*L[I])/(EI+Spring_Stiff[i]*L[I]);
	R = 12 - 8*(r[0]+r[1]) + (5*r[0]*r[1]);
	EI = EI/R;
	A[I] = (4*EI/L[I]) * (3*r[0]-2*r[0]*r[1]);
    B[I] = (6*EI/(L[I]*L[I])) * (2*r[0]-r[0]*r[1]);
    C[I] = (2*EI/L[I]) *(r[0]*r[1]);
    D[I] = (12*EI/(L[I]*L[I]*L[I])) * (r[0]+r[1]-r[0]*r[1]);
    E[I] = EA/L[I];
    for (i=0;i<6;i++)
    for (j=0;j<6;j++)
    K[i][j] = 0.0;
	 //I:2 number of element  - II: kind of stiffness matrix
	 if (II==1){// No plastic hinge
    K[0][0]=E[I];K[0][1]=0;K[0][2]=0;K[0][3]=-E[I];K[0][4]=0;K[0][5]=0;
	K[1][0]=0;K[1][1]=D[I];K[1][2]=B[I];K[1][3]=0;K[1][4]=-D[I];K[1][5]=B[I];
	K[2][0]=0;K[2][1]=B[I];K[2][2]=A[I];K[2][3]=0;K[2][4]=-B[I];K[2][5]=C[I];
	K[3][0]=-E[I];K[3][1]=0;K[3][2]=0;K[3][3]=E[I];K[3][4]=0;K[3][5]=0;
	K[4][0]=0;K[4][1]=-D[I];K[4][2]=-B[I];K[4][3]=0;K[4][4]=D[I];K[4][5]=-B[I];
	K[5][0]=0;K[5][1]=B[I];K[5][2]=C[I];K[5][3]=0;K[5][4]=-B[I];K[5][5]=A[I];
	 }
	if (II==2){// plastic hinge at i
    K[0][0]=E[I];K[0][1]=0;K[0][2]=0;K[0][3]=-E[I];K[0][4]=0;K[0][5]=0;
	K[1][0]=0;K[1][1]=.25*D[I];K[1][2]=0;K[1][3]=0;K[1][4]=-.25*D[I];K[1][5]=.5*B[I];
	K[2][0]=0;K[2][1]=0;K[2][2]=0;K[2][3]=0;K[2][4]=0;K[2][5]=0;
	K[3][0]=-E[I];K[3][1]=0;K[3][2]=0;K[3][3]=E[I];K[3][4]=0;K[3][5]=0;
	K[4][0]=0;K[4][1]=-.25*D[I];K[4][2]=0;K[4][3]=0;K[4][4]=.25*D[I];K[4][5]=-.5*B[I];
	K[5][0]=0;K[5][1]=.5*B[I];K[5][2]=0;K[5][3]=0;K[5][4]=-.5*B[I];K[5][5]=(3/4)*A[I];
	 }
	if (II==3){// plastic hinge at j
    K[0][0]=E[I];K[0][1]=0;K[0][2]=0;K[0][3]=-E[I];K[0][4]=0;K[0][5]=0;
	K[1][0]=0;K[1][1]=.25*D[I];K[1][2]=.5*B[I];K[1][3]=0;K[1][4]=-.25*D[I];K[1][5]=0;
	K[2][0]=0;K[2][1]=.5*B[I];K[2][2]=(3/4)*A[I];K[2][3]=0;K[2][4]=-.5*B[I];K[2][5]=0;
	K[3][0]=-E[I];K[3][1]=0;K[3][2]=0;K[3][3]=E[I];K[3][4]=0;K[3][5]=0;
	K[4][0]=0;K[4][1]=-.25*D[I];K[4][2]=-.5*B[I];K[4][3]=0;K[4][4]=.25*D[I];K[4][5]=0;
	K[5][0]=0;K[5][1]=0;K[5][2]=0;K[5][3]=0;K[5][4]=0;K[5][5]=0;
	 }
	if (II==4){// plastic hinge at i and j
    K[0][0]=E[I];K[0][1]=0;K[0][2]=0;K[0][3]=-E[I];K[0][4]=0;K[0][5]=0;
	K[1][0]=0;K[1][1]=0;K[1][2]=0;K[1][3]=0;K[1][4]=0;K[1][5]=0;
	K[2][0]=0;K[2][1]=0;K[2][2]=0;K[2][3]=0;K[2][4]=0;K[2][5]=0;
	K[3][0]=-E[I];K[3][1]=0;K[3][2]=0;K[3][3]=E[I];K[3][4]=0;K[3][5]=0;
	K[4][0]=0;K[4][1]=0;K[4][2]=0;K[4][3]=0;K[4][4]=0;K[4][5]=0;
	K[5][0]=0;K[5][1]=0;K[5][2]=0;K[5][3]=0;K[5][4]=0;K[5][5]=0;
	 }
	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;
	Matrix_Transpose(lan,lan_Tr);
	Matrix_Multiplication(lan_Tr,K,ans);
	Matrix_Multiplication(ans,lan,K_G);
}
void ELEMNT_FORCE_OUTPUT(double eleF[1][6],double ELE_FORCE[6][STEP],int n){
      int i;
      for (i=0;i<6;i++)
      ELE_FORCE[i][n]=eleF[0][i];
}
double SQRT2(double D){
            int it,itermax;
            double residual,tolerance,x,dx,dx_ABS,f,df;
			it = 0; // initialize iteration count
		    itermax = 100000;
		    residual = 100; // initialize residual
		    tolerance = 1e-12;
		    x = 1;// initialize answer
		    while (residual > tolerance){
		    	f = x*x - D;
		    	df = 2 * x;
		    	dx = f/df;
		        x= x - dx;
                residual = ABS(dx); // abs residual
		        it = it + 1; // increment iteration count
		        //printf("f: %f -\tdx: %f -\tresidual: %f\n",f,dx,residual);
		         if (it == itermax){
		          //printf("\tSQRT2(number,power) : SQRT2(%f) - iteration: %d ->   ## The solution is not converged ##\n",D,it);
		          break;
				}
		    }
		        if (it < itermax){
			   //printf("\tSQRT(number,power) - SQRT(%f,%f) : %f \n",D,n, x);
			   return x;
			   }
}
void Matrix_Transpose(double A[6][6],double B[6][6]){
	int i,j;
	for (i=0;i<6;i++)
	for (j=0;j<6;j++)
	B[j][i]=A[i][j];
}
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]){
	int i,j,k;
	double sum;
	for (i=0;i<6;i++)
	for (j=0;j<6;j++){
	sum=0;
	for (k=0;k<6;k++)
	sum += A[i][k]*B[k][j];
	C[i][j] = sum;
	}
}
void MessageResult(double output_base01[],double output_base02[],double output_u01[],int n){
	int i;
	printf("\t   ");
	for (i=0;i<77;i++)
	printf("-");
	printf("\n");
    printf("\t    Increment    Base Shear[DOF(1)]     Base Moment[DOF(3)]     Disp. [DOF(5)]\n");
	printf("\t   ");
	for (i=0;i<77;i++)
	printf("-");
	printf("\n");
	for (i=0;i<n;i++){
	//Distance(i+1);	
	printf("\t\t %d\t   %.3e\t          %.3e\t          %.3e\n",i+1,output_base01[i],output_base01[i],output_u01[i]);	
	}
}
void MessageStrCoorTEXT(double X[],double Y[],int n){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188;Qf=186;
   	printf("     %c",Qa);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
	printf("     %c   Structural Coordinate Data   %c\n",Qf,Qf);
	printf("     %c       X              Y         %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<n;i++)
	printf("          %.3e     %.3e\n",X[i],Y[i]);
}
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double output_u01[],double output_base01[],double output_base02[],double C[6][STEP],int n,int M){
// HTML OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText06, "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"PushoverLinearSemiRigidBerFC-image01.png\" style=\"width:1000px ; height:500px\" alt=\"analysis\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th bgcolor=\"cyan\">   Pushover Semi Rigid Connection Analysis Force Analogy Method with Force Control Based on Euler-Bernoulli Beam Theory - Output Report </th> \n");
    // TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Input Data </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section flextural rigidity - EI: </th><th> %.3e </th> </tr>\n",EI);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section axial rigidity - EA: </th><th> %.3e </th> </tr>\n",EA);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External Incremental Fx force [DOF(4)]: </th><th> %.3e </th> </tr>\n",Fini);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Number of increment: </th><th> %d </th> </tr>\n",M);   
    // TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"3\" bgcolor=\"orange\"> Structral Coordinate </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Node Number </th><th bgcolor=\"orange\">X Coordinate</th> <th bgcolor=\"orange\">Y Coordinate</th> </tr>\n");
	for(i=0;i<Ne+1;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td>\n",i+1,X[i],Y[i]);
    }
	// TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"4\" bgcolor=\"orange\"> Structral Deformation </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Base Shear [DOF(1)]</th><th bgcolor=\"orange\">Base Moment [DOF(3)]</th><th bgcolor=\"orange\">Displacement [DOF(5)]</th> </tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> </tr>\n",i+1,output_base01[i],output_base02[i],output_u01[i]);
    }
    // TABLE 4
    fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
    fprintf(OutputFile,"<tr><th colspan=\"7\" bgcolor=\"orange\"> Structral Element Internal Forces </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Axial-i</th> <th bgcolor=\"orange\">Shear-i</th> <th bgcolor=\"orange\">Moment-i</th> <th bgcolor=\"orange\">Axial-j</th> <th bgcolor=\"orange\">Shear-j</th> <th bgcolor=\"orange\">Moment-j</th> </tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> </tr>\n",i+1,C[0][i],C[1][i],C[2][i],C[3][i],C[4][i],C[5][i]);
    }
	fprintf(OutputFile,"</table></body></html>\n");
	fclose(OutputFile);
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
	OutputFile = fopen(ShowText07, "w");
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
//int MessageControl(double eleF[Ne][6],double u[]){
//      int i;
//      if (ABS(eleF[0][2]) >= MOM01[0]){
//      printf("\t STEP: 1 - Moment-1: %f reaached to yield: %f \n",ABS(eleF[0][2]),MOM01[0]);
//	  i=1;
//	  }
//	  if (ABS(eleF[0][5]) >= MOM02[0]){
//      printf("\t STEP: 1 - Moment-2: %f reaached to yield: %f \n",ABS(eleF[0][5]),MOM02[0]);
//	  i=1;
//	  }	
//	  if (ABS(u[1]) >= TET01[n-1]){
//      printf("\t STEP: 3 - Rotation-1: %f reaached to yield: %f \n",ABS(u[1]),TET01[n-1]);
//	  i=1;
//	  }
//	  if (ABS(u[2]) >= TET02[n-1]){
//      printf("\t STEP: 3 - Rotation-2: %f reaached to yield: %f \n",ABS(u[2]),TET02[n-1]);
//	  i=1;
//	  }
//	  return i;
//}

