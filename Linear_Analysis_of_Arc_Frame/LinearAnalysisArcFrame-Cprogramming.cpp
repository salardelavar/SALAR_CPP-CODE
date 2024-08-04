#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 33 // Degree of freedom
#define Ne 12    // number of element
#define N 13    // number of node
#define ShowText01 "LinearAnalysisArcFrame-inputDATA.csv"
#define ShowText02 "LinearAnalysisArcFrame-inputFORCE.csv"
#define ShowText03 "LinearAnalysisArcFrame-inputCOORDINATE.csv"
#define ShowText04 "Output data is written in Excel,Autocad and Html file"
#define ShowText05 "LinearAnalysisArcFrame-outputEXCEL.csv"
#define ShowText06 "LinearAnalysisArcFrame-outputHTML.html"
#define ShowText07 "Graph-outputHTML.html"
#define ShowText08 "                  >> IN THE NAME OF ALLAH <<              "
#define ShowText09 "                 Linear Analysis of Arc Frame             "
#define ShowText10 "                        UNIT: Free Unit                   "
#define ShowText11 "    This program is written by Salar Delavar Qashqai      "
#define ShowText12 "             E-mail: salar.d.ghashghaei@gmail.com         "

void IMPORT_DATA01(double &EA,double &EI);
void IMPORT_DATA02(double F[11][3],int &n);
void IMPORT_DATA03(double x[],double y[],int &n);
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
void MatrixDetermination(double [][NN],int );
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [],int );
void Matrix_Transpose(double A[6][6],double B[6][6]);
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]);
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I);// Calculate internal element force
void ELEMNT_FORCE_OUTPUT(double eleF[12][6],double ELE_FORCE[72]);
double ABS(double);
double MIN(double A[],int n);
double SQRT2(double D);
void MessageInitialData(double EI,double EA,double x[],double y[],double F[11][3],int m);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double EI,double EA);
void MessageStrCoorTEXT(double X[],double Y[],int n);
void MessageResult(double A[],double B[33],int n);
void OUTPUT_excel(double A[33],double B[],double C[72],int n);
void OUTPUT_html(double F[11][3],double A[],double C[72],double X[],double Y[],int n,int m);
void OUTPUT_autocad(double X[13],double Y[13],double A[33],double C[72],int n);
void ANALYSIS(double EI,double EA,double x[],double y[],double F[11][3]);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
double MAX_ABS(double A[],int n);
int main(){
    double EI,EA,x[13],y[13],F[11][3];
    int n,m;
    IMPORT_DATA01(EA,EI);
    IMPORT_DATA02(F,m);
    IMPORT_DATA03(x,y,n);
    MessageCheck_IMPORT_DATA01(EI,EA);
	textcolor(14);
    MessageInitialData(EI,EA,x,y,F,m);
    MessageStrCoorTEXT(x,y,n);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(EI,EA,x,y,F);
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
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I){
    double lan[6][6],UU[6],ff,ll[6][6];
    int i,j,II;
	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;

	if (I == 0){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (I == 11){
	UU[0]=U[30];UU[1]=U[31];UU[2]=U[32];UU[3]=0;UU[4]=0;UU[5]=0;
	}
	if (I >= 1 && I < 11){
	i=3*(I-1);
	UU[0]=U[i];UU[1]=U[i+1];UU[2]=U[i+2];UU[3]=U[i+3];UU[4]=U[i+4];UU[5]=U[i+5];
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

void MessageInitialData(double EI,double EA,double x[],double y[],double F[11][3],int m){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText08,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText09,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText10,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
	printf("\t\t\t\t%c%s%c\n",Qf,ShowText11,Qf);
	printf("\t\t\t\t%c%s%c\n",Qf,ShowText12,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<59;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);


    MessageInputDataTEXT();
    printf("      Section flextural rigidity - EI:                          %.3e\n",EI);
    printf("      Section axial rigidity - EA:                              %.3e\n",EA);
    for (i=0;i<m;i++){
    printf("      Force-x - Node number %d:                                  %.3e\n",i+1,F[i][0]);
	printf("      Force-y - Node number %d:                                  %.3e\n",i+1,F[i][1]);
	printf("      Moment-z  - Node number %d:                                %.3e\n",i+1,F[i][2]);
	}

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
void MessageCheck_IMPORT_DATA01(double EI,double EA){
	if (EI < 0 ||    EA < 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative data input value is not acceptable ***\n");
	printf("  Section flextural rigidity - EI:     %.3e\n",EI);
    printf("  Section axial rigidity - EA:         %.3e\n",EA);
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
void OUTPUT_excel(double A[33],double B[],double C[72],int n){
// EXCEL OUTPUT
	int i,I;
	FILE *OutputFile;
	OutputFile = fopen(ShowText05, "w");
fprintf(OutputFile," ###                 Linear Analysis of Arc Frame             ###\n");
fprintf(OutputFile,"Base Shear[DOF(1)]+[DOF(37)]: ");
fprintf(OutputFile,"%e\n",B[0]);
fprintf(OutputFile,"Element Number %d\n",1);
fprintf(OutputFile,"Displacement-i,Displacement-i,Rotation-i,Displacement-j,Displacement-j,Rotation-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",0.0,0.0,0.0,A[0],A[1],A[2]);
fprintf(OutputFile,"axial-i,shear-i,moment-i,axial-j,shear-j,moment-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",C[0],C[1],C[2],C[3],C[4],C[5]);
for(i=1;i<n-1;i++){
I=3*(i-1);
fprintf(OutputFile,"Element Number %d\n",i+1);
fprintf(OutputFile,"Displacement-i,Displacement-i,Rotation-i,Displacement-j,Displacement-j,Rotation-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",A[I],A[I+1],A[I+2],A[I+3],A[I+4],A[I+5]);
fprintf(OutputFile,"axial-i,shear-i,moment-i,axial-j,shear-j,moment-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",C[I],C[I+1],C[I+2],C[I+3],C[I+4],C[I+5]);
}
fprintf(OutputFile,"Element Number %d\n",12);
fprintf(OutputFile,"Displacement-i,Displacement-i,Rotation-i,Displacement-j,Displacement-j,Rotation-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",A[30],A[31],A[32],0.0,0.0,0.0);
fprintf(OutputFile,"axial-i,shear-i,moment-i,axial-j,shear-j,moment-j\n");
fprintf(OutputFile,"%e,%e,%e,%e,%e,%e\n",C[66],C[67],C[68],C[69],C[70],C[71]);
fclose(OutputFile);
}
void OUTPUT_autocad(double X[13],double Y[13],double A[33],double C[72],int n){
// AUTOCAD OUTPUT
int i,I;
FILE *OutputFile;
OutputFile = fopen("LinearAnalysisArcFrame-outputAUTOCAD.scr", "w");
//fprintf(OutputFile,"  %%                  Linear Analysis of Arc Frame            %%\n");

for(i=0;i<n-1;i++){
fprintf(OutputFile,"LINE\n");
fprintf(OutputFile,"%e,%e\n",X[i],Y[i]);
fprintf(OutputFile,"%e,%e\n",X[i+1],Y[i+1]);
}
//fprintf(OutputFile,"LINE\n");
//fprintf(OutputFile,"%e,%e\n",X[0],Y[0]);
//fprintf(OutputFile,"%e,%e\n",X[1]+A[0]*100,Y[1]+A[1]*100);
//for(i=1;i<n-2;i++){
//I=3*(i-1);
//fprintf(OutputFile,"LINE\n");
//fprintf(OutputFile,"%e,%e\n",X[i]+A[I]*100,Y[I]+A[I]*100);
//fprintf(OutputFile,"%e,%e\n",X[i+1]+A[I+3]*100,Y[i+1]+A[I+3]*100);
//}
//fprintf(OutputFile,"LINE\n");
//fprintf(OutputFile,"%e,%e\n",X[31]+A[30]*100,Y[31]+A[3]*100);
//fprintf(OutputFile,"%e,%e\n",X[13],Y[13]);
fprintf(OutputFile,"ZOOM\n");
fprintf(OutputFile,"E\n");
fprintf(OutputFile,"ZOOM\n");
fprintf(OutputFile,".8X\n");
fclose(OutputFile);
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
void IMPORT_DATA01(double &EA,double &EI){
    double Import_Data[2];
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
	while(i < 3 && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
	EI=Import_Data[0];
	EA=Import_Data[1];
}
void IMPORT_DATA02(double F[11][3],int &n){
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
	fscanf(InputFile,"%lf,%lf,%lf",&F[i][0],&F[i][1],&F[i][2]);
	//printf("%e,%e,%e\n",F[i][0],F[i][1],F[i][2]);
	i++;
	}
	while(i < 11 && fgets(line,sizeof(line),InputFile) != NULL);
	n = i;
    //printf("%d\n",n);
}
void IMPORT_DATA03(double x[],double y[],int &n){
	int i = 0;
	FILE *InputFile;
	InputFile = fopen(ShowText03, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [%s] \n",ShowText03);
		Sleep(6000);
		exit(1);
	}
	char line[1000];
	do{
	fscanf(InputFile,"%lf,%lf",&x[i],&y[i]);
	i++;
	}
	while(i < 13 && fgets(line,sizeof(line),InputFile) != NULL);
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
void ANALYSIS(double EI,double EA,double x[],double y[],double F[11][3]){
   int i,j,I;
   double K[NN][NN],Kinv[NN][NN],eleF[Ne][6],ELE_FORCE[72];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],FF[NN],u[NN];
   double output_u[33],output_base[1];
   double MS[6][6],KG[6][6];
   	for (i=0;i<NN;i++)
    for (j=0;j<NN;j++)
    K[i][j]=0;
    for (i=0;i<NN;i++)
    u[i] = 0;

	  for(i=0;i<6*Ne;i++)
      ELE_FORCE[i] = 0;

	  for(i=0;i<3;i++){
	  FF[i]=F[0][i];
      FF[i+3]=F[1][i];
      FF[i+6]=F[2][i];
      FF[i+9]=F[3][i];
      FF[i+12]=F[4][i];
      FF[i+15]=F[5][i];
      FF[i+18]=F[6][i];
      FF[i+21]=F[7][i];
      FF[i+24]=F[8][i];
      FF[i+27]=F[9][i];
      FF[i+30]=F[10][i];
	  }

    for (i=0;i<Ne;i++){
    L[i] = SQRT2((x[i+1]-x[i])*(x[i+1]-x[i])+(y[i+1]-y[i])*(y[i+1]-y[i]));
	lanX[i] = (x[i+1]-x[i])/L[i];lanY[i] = (y[i+1]-y[i])/L[i];
	}

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
    K[0][0]= KG[3][3];//DOF(4)
    K[0][1]= KG[3][4];//DOF(4)
    K[0][2]= KG[3][5];//DOF(4)

    K[1][0]= KG[4][3];//DOF(5)
    K[1][1]= KG[4][4];//DOF(5)
    K[1][2]= KG[4][5];//DOF(5)

    K[2][0]= KG[5][3];//DOF(6)
    K[2][1]= KG[5][4];//DOF(6)
    K[2][2]= KG[5][5];//DOF(6)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,11,1);
    K[30][30]= KG[0][0];//DOF(34)
    K[30][31]= KG[0][1];//DOF(34)
    K[30][32]= KG[0][2];//DOF(34)

    K[31][30]= KG[1][0];//DOF(35)
    K[31][31]= KG[1][1];//DOF(35)
    K[31][32]= KG[1][2];//DOF(35)

    K[32][30]= KG[2][0];//DOF(36)
    K[32][31]= KG[2][1];//DOF(36)
    K[32][32]= KG[2][2];//DOF(36)
    for (i=1;i<11;i++){
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,i,1);
    I=3*(i-1);
    K[I][I]+= KG[0][0];//DOF(7)
    K[I][I+1]+= KG[0][1];//DOF(7)
    K[I][I+2]+= KG[0][2];//DOF(7)

    K[I+1][I]+= KG[1][0];//DOF(8)
    K[I+1][I+1]+= KG[1][1];//DOF(8)
    K[I+1][I+2]+= KG[1][2];//DOF(8)

    K[I+2][I]+= KG[2][0];//DOF(9)
    K[I+2][I+1]+= KG[2][1];//DOF(9)
    K[I+2][I+2]+= KG[2][2];//DOF(9)

    K[I+3][I+3]+= KG[3][3];//DOF(10)
    K[I+3][I+4]+= KG[3][4];//DOF(10)
    K[I+3][I+5]+= KG[3][5];//DOF(10)

    K[I+4][I+3]+= KG[4][3];//DOF(11)
    K[I+4][I+4]+= KG[4][4];//DOF(11)
    K[I+4][I+5]+= KG[4][5];//DOF(11)

    K[I+5][I+3]+= KG[5][3];//DOF(12)
    K[I+5][I+4]+= KG[5][4];//DOF(12)
    K[I+5][I+5]+= KG[5][5];//DOF(12)

    K[I][I+3]= KG[0][3];//DOF(7)
    K[I][I+4]= KG[0][4];//DOF(7)
    K[I][I+5]= KG[0][5];//DOF(7)

    K[I+1][I+3]= KG[1][3];//DOF(8)
    K[I+1][I+4]= KG[1][4];//DOF(8)
    K[I+1][I+5]= KG[1][5];//DOF(8)

    K[I+2][I+3]= KG[2][3];//DOF(9)
    K[I+2][I+4]= KG[2][4];//DOF(9)
    K[I+2][I+5]= KG[2][5];//DOF(9)

    K[I+3][I]= KG[3][0];//DOF(10)
    K[I+3][I+1]= KG[3][1];//DOF(10)
    K[I+3][I+2]= KG[3][2];//DOF(10)

    K[I+4][I]= KG[4][0];//DOF(11)
    K[I+4][I+1]= KG[4][1];//DOF(11)
    K[I+4][I+2]= KG[4][2];//DOF(11)

    K[I+5][I]= KG[5][0];//DOF(12)
    K[I+5][I+1]= KG[5][1];//DOF(12)
    K[I+5][I+2]= KG[5][2];//DOF(12)
	}


      MatrixDetermination(K,NN);
	  MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,NN);
      for (i=0;i<Ne;i++){
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,i,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,i);
	  }
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE);
      for (i=0;i<NN;i++)
      output_u[i]=u[i];//output displacement
      //output_base[0]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[7];//output base shear

      // So plastic hinge has formed in 1 node


       MessageResult(output_base,output_u,3);
      MatrixDetermination(K,NN);
	    OUTPUT_excel(output_u,output_base,ELE_FORCE,Ne);
	    OUTPUT_autocad(x,y,output_u,ELE_FORCE,13);
	    OUTPUT_html(F,output_u,ELE_FORCE,x,y,NN,Ne);
	    char text1[20]="Shape of Structure",text2[20]="X dir",text3[20]="Y dir";
	    for (i=1;i<12;i++){
        I=3*(i-1); x[i] = x[i] +  output_u[I];
        I=I+1; y[i] = y[i] +  output_u[I];
	    }
	    OUTPUT_HTML_GRAPH(x,y,Ne+1,text1,text2,text3);
		textcolor(15);
		printf("\n\a - %s -\n",ShowText04);
		system("start /w Graph-outputHTML.html");
		DATE_TIME();

}
void DATE_TIME(){
		printf("\n\t");
		system("echo %date%");
		printf("\t");
		system("echo %time%");
}
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II){
	double lan[6][6],lan_Tr[6][6],ans[6][6];
	A[I] = 4*EI/L[I];
    B[I] = 6*EI/(L[I]*L[I]);
    C[I] = 2*EI/L[I];
    D[I] = 12*EI/(L[I]*L[I]*L[I]);
    E[I] = EA/L[I];
     for (int i=0;i<6;i++)
	 for (int j=0;j<6;j++)
	 K[i][j] = 0;
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
void ELEMNT_FORCE_OUTPUT(double eleF[12][6],double ELE_FORCE[72]){
      int i;
      for (i=0;i<6;i++)
      ELE_FORCE[i]=eleF[0][i];
      for (i=6;i<12;i++)
      ELE_FORCE[i]=eleF[1][i-6];
      for (i=12;i<18;i++)
      ELE_FORCE[i]=eleF[2][i-12];
      for (i=18;i<24;i++)
      ELE_FORCE[i]=eleF[3][i-18];
      for (i=24;i<30;i++)
      ELE_FORCE[i]=eleF[4][i-24];
      for (i=30;i<36;i++)
      ELE_FORCE[i]=eleF[5][i-30];
      for (i=36;i<42;i++)
      ELE_FORCE[i]=eleF[6][i-36];
      for (i=42;i<48;i++)
      ELE_FORCE[i]=eleF[7][i-42];
      for (i=48;i<54;i++)
      ELE_FORCE[i]=eleF[8][i-48];
      for (i=54;i<60;i++)
      ELE_FORCE[i]=eleF[9][i-54];
      for (i=60;i<66;i++)
      ELE_FORCE[i]=eleF[10][i-60];
      for (i=66;i<72;i++)
      ELE_FORCE[i]=eleF[11][i-66];

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
void MessageResult(double A[],double B[33],int n){
	int i;
	printf("\t   ");
	for (i=0;i<147;i++)
	printf("-");
	printf("\n");
    printf("\t    Increment    Base Shear[DOF(1)]+[DOF(4)]  Disp. [DOF(7)]  Disp. [DOF(8)]  Rotation [DOF(9)]  Disp. [DOF(10)]  Disp. [DOF(11)]  Rotation [DOF(12)]\n");
	printf("\t   ");
	for (i=0;i<147;i++)
	printf("-");
	printf("\n");
	for (i=0;i<n;i++)
    printf("\t\t %d\t\t   %.3e\t        %.3e\t%.3e\t  %.3e\t  %.3e     %.3e\t  %.3e\n",i+1,A[i],B[0],B[1],B[2],B[3],B[4],B[5]);
}
void MessageStrCoorTEXT(double X[],double Y[],int n){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188;Qf=186;
   	printf("     %c",Qa);
	for (i=1;i<42;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
	printf("     %c         Structural Coordinate Data      %c\n",Qf,Qf);
	printf("     %c  Node Num.     X              Y         %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<42;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<n;i++)
	printf("\t    %d     %.3e    %.3e\n",i+1,X[i],Y[i]);
}
void OUTPUT_html(double F[11][3],double A[],double C[72],double X[],double Y[],int n,int m){
// HTML OUTPUT
	int i,I;
	FILE *OutputFile;
	OutputFile = fopen(ShowText06, "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"LinearAnalysisArcFrame.png\" style=\"width:1000px ; height:800px\" alt=\"analysis\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th bgcolor=\"cyan\"> Linear Analysis of Arc Frame - Output Report </th> \n");
	// TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"3\" bgcolor=\"orange\"> Structral Coordinate </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Node Number </th><th bgcolor=\"orange\">X Coordinate</th> <th bgcolor=\"orange\">Y Coordinate</th> </tr>\n");
	for(i=0;i<Ne+1;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td>\n",i+1,X[i+1],Y[i+1]);
    }
	// TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"4\" bgcolor=\"orange\"> Structral External Applied Load </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Node Number </th> <th bgcolor=\"orange\"> Fx </th> <th bgcolor=\"orange\"> Fy </th> <th bgcolor=\"orange\"> Mz </th> </tr>\n");
	for(i=0;i<Ne-1;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> </tr>\n",i+2,F[i][0],F[i][1],F[8][2]);
    }
    // TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Structral Deformation </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Degree of Freedom</th> <th bgcolor=\"orange\">Deformation</th> </tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td> </tr>\n",i+4,A[i]);
    }
    // TABLE 4
    fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
    fprintf(OutputFile,"<tr><th colspan=\"7\" bgcolor=\"orange\"> Structral Element Internal Forces </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Element Number</th> <th bgcolor=\"orange\">Axial-i</th> <th bgcolor=\"orange\">Shear-i</th> <th bgcolor=\"orange\">Moment-i</th> <th bgcolor=\"orange\">Axial-j</th> <th bgcolor=\"orange\">Shear-j</th> <th bgcolor=\"orange\">Moment-j</th> </tr>\n");
	for(i=1;i<m-1;i++){
	I=3*(i-1);
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> </tr>\n",i+1,C[I],C[I+1],C[I+2],C[I+3],C[I+4],C[I+5]);
    }
	fprintf(OutputFile,"</table></body></html>\n");
	fclose(OutputFile);
}
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]){
    // HTML GRAPH OUTPUT
	int i;
	double x,y,NorX[N],NorY[N],Xmax,Ymax;
	Xmax=MAX_ABS(X,n);
	Ymax=MAX_ABS(Y,n);
	for (i=0;i<n;i++){
	NorX[i] = X[i]/Xmax;
	NorY[i] = Y[i]/Ymax;
	//printf("\t %f   %f    \n",NorX[i],NorY[i]);
	}
	FILE *OutputFile;
	OutputFile = fopen(ShowText07, "w");
	fprintf(OutputFile,"<!DOCTYPE HTML><html><body style=\"background-color:black;\"><font color=\"white\"><head><script> \n");
	fprintf(OutputFile,"window.onload = function(){ \n");
	fprintf(OutputFile,"var canvas = document.getElementById(\"myCanvas\");var s1 = canvas.getContext(\"2d\");var s2 = canvas.getContext('2d'); \n");
	fprintf(OutputFile,"var s3 = canvas.getContext(\"2d\");var s4 = canvas.getContext(\"2d\");var s5 = canvas.getContext(\"2d\"); \n");
	fprintf(OutputFile,"var x=120,y=100,X,Y,Lx=1100,Ly=500,i; \n");
	fprintf(OutputFile,"s3.beginPath();s3.lineWidth = 3;s3.strokeStyle = \"cyan\";s3.rect(x,y,Lx,Ly); \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x+Lx*(i+1)*.1,y+Ly);s3.lineTo(x+Lx*(i+1)*.1,y+Ly-10);}; \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x,y+Ly*(i+1)*.1);s3.lineTo(x+10,y+Ly*(i+1)*.1);};s3.stroke();\n");
	fprintf(OutputFile,"s1.beginPath();s1.lineWidth = 3;s1.strokeStyle = \"yellow\"; \n");
	for (i=0;i<n-1;i++){
	fprintf(OutputFile,"s1.moveTo(%f,%f); \n",120+NorX[i]*1100,100+500-NorY[i]*500);
	fprintf(OutputFile,"s1.lineTo(%f,%f); \n",120+NorX[i+1]*1100,100+500-NorY[i+1]*500);
	}
	fprintf(OutputFile,"s1.stroke(); \n");
	fprintf(OutputFile,"s2.beginPath();s2.lineWidth = 1;s2.strokeStyle = \"cyan\";s2.setLineDash([5, 5]); \n");
	fprintf(OutputFile,"for(i=0;i<19;i++){s2.moveTo(x+Lx*(i+1)*.05,y);s2.lineTo(x+Lx*(i+1)*.05,y+Ly);} \n");
	fprintf(OutputFile,"s2.lineWidth = 1;s2.strokeStyle = \"cyan\";for(i=0;i<19;i++){s2.moveTo(x,y+Ly*(i+1)*.05);s2.lineTo(x+Lx,y+Ly*(i+1)*.05);} s2.stroke();\n");
	fprintf(OutputFile,"X=x+.25*Lx;Y=.7*y;s4.translate(X,Y);s4.font=\"100px serif\";s4.fillStyle = \"#7fff00\";s4.fillText(\"%s\",0,0); \n",text1);
	fprintf(OutputFile,"s4.save();X=-X+.2*x;Y=-Y+y+.6*Ly;s4.translate(X,Y);s4.rotate(3*Math.PI/2);s4.font=\"30px serif\"; \n");
	fprintf(OutputFile,"s4.fillStyle = \"#7fff00\";s4.textAlign = \"left\";s4.fillText(\"%s\",0,0);s4.restore(); \n",text3);
	fprintf(OutputFile,"s4.save();X=.2*Lx;Y=y+Ly-20;s4.translate(X,Y);s4.rotate(2*Math.PI);s4.font=\"30px serif\";s4.fillStyle = \"#7fff00\"; \n");
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

