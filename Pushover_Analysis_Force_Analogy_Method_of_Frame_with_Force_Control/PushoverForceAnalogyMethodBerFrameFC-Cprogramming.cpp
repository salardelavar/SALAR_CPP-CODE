#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 12 // Degree of freedom
#define Ne 3    // number of element
#define N 4    // number of node
#define STEP 5000   // number of node
#define ShowText01 "PushoverForceAnalogyMethodBerFrameFC-inputDATA.csv"
#define ShowText02 "PushoverForceAnalogyMethodBerFrameFC-inputCOORDINATE.csv"
#define ShowText03 "PushoverForceAnalogyMethodBerFrameFC-inputHINGE.csv"
#define ShowText04 "Output data is written in Excel and Html file"
#define ShowText05 "PushoverForceAnalogyMethodBerFrameFC-outputEXCEL.csv"
#define ShowText06 "PushoverForceAnalogyMethodBerFrameFC-outputHTML.html"
#define ShowText07 "Graph-outputHTML.html"
#define a01 "      >> In the name of God, the Gracious, the Merciful <<          "
#define a02 " Pushover Analysis Force Analogy Method of Frame with Force Control "
#define a03 "                Based on Euler-Bernoulli Beam Theory                "
#define a04 "                          UNIT: Free Unit                           "
#define a05 "        This program is written by Salar Delavar Ghashghaei         "
#define a06 "                  E-mail: salar.d.ghashghaei@gmail.com              "

void IMPORT_DATA01(double &EA,double &EI,double &Fini,int &M);
void IMPORT_DATA02(double x[],double y[],int &n);
void IMPORT_DATA03(double TET[],double MOM[],int &k);
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
void MatrixDetermination(double A[][NN],int n,double &Product);
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [],int );
void Matrix_Transpose(double A[6][6],double B[6][6]);
void MatrixZero(double A[NN][NN],int n);
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]);
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I,int II);// Calculate internal element force
void ELEMNT_FORCE_OUTPUT(double eleF[Ne][6],double ELE_FORCE[STEP][18],int I);
double ABS(double);
double MIN(double A[],int n);
double SQRT2(double D);
void MessageInitialData(double EI,double EA,double Fini,double x[],double y[],int M);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double EI,double EA,int M);
void MessageCheck_IMPORT_DATA03(double TET[],double MOM[],int n);
void MessageCheckMk(int M,int K);
void MessageStrCoorTEXT(double X[],double Y[],int n);
void MessagePlasticHingeTEXT(double TET[],double MOM[],int n);
void MessageResult(double output_base01[],double output_base02[],double output_u01[],int n);
int MessageCheckMoment(double eleF[Ne][6],int z,double MOM[],int &I);
int MessageCheckRotation(double u[],int z,double TET[],int np,int &I);
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_u05[],double output_u06[],double output_u07[],double output_u08[],double output_u09[],double output_u10[],double output_base01[],double output_base02[],double C[STEP][18],int n);
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double TET[],double MOM[],double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_u05[],double output_u06[],double output_u07[],double output_u08[],double output_u09[],double output_u10[],double output_base01[],double output_base02[],double C[STEP][18],int n,int m,int M);
void ANALYSIS(double TET[],double MOM[],int np,double EI,double EA,double Fini,double x[],double y[],int M);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
double MAX_ABS(double A[],int n);
void PlasticHingeStiffnessCOFF(double A[],double B[],double C[],int n);
void PlasticHingeStiffness(double A[1][6],int I,int II,double B[],double C[],double MomS[],double KRK[],int n,int III);
void MatrixDisplay(double A[NN][NN],int n);

int main(){
    double EI,EA,Fini,x[4],y[4],TET[10],MOM[10];
    int n,np,m,M;
    IMPORT_DATA01(EA,EI,Fini,M);
    IMPORT_DATA02(x,y,n);
    IMPORT_DATA03(TET,MOM,np);
    MessageCheck_IMPORT_DATA01(EI,EA,M);
    MessageCheck_IMPORT_DATA03(TET,MOM,np);
    MessageCheckMk(M,np);
	textcolor(14);
    MessageInitialData(EI,EA,Fini,x,y,M);
    MessageStrCoorTEXT(x,y,n);
    MessagePlasticHingeTEXT(TET,MOM,np);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(TET,MOM,np,EI,EA,Fini,x,y,M);
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
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I,int II){
    double lan[6][6],UU[6],ff,ll[6][6];
    int i,j;
	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;

	if (II == 11){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (II == 12){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 13){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[2];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 21){
	UU[0]=0;UU[1]=0;UU[2]=U[6];UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (II == 22){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 23){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[2];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 31){
	UU[0]=0;UU[1]=0;UU[2]=U[6];UU[3]=U[0];UU[4]=U[1];UU[5]=U[7];
	}
	if (II == 32){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 33){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[2];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 41){
	UU[0]=0;UU[1]=0;UU[2]=U[6];UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (II == 42){
	UU[0]=0;UU[1]=0;UU[2]=U[7];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 43){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[8];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 51){
	UU[0]=0;UU[1]=0;UU[2]=U[6];UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (II == 52){
	UU[0]=0;UU[1]=0;UU[2]=U[7];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (II == 53){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[8];UU[3]=U[3];UU[4]=U[4];UU[5]=U[9];
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

void MessageInitialData(double EI,double EA,double Fini,double x[],double y[],int M){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c%s%c\n",Qf,a01,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,a02,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,a03,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,a04,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
	printf("\t\t\t\t%c%s%c\n",Qf,a05,Qf);
	printf("\t\t\t\t%c%s%c\n",Qf,a06,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

    MessageInputDataTEXT();
    printf("      Section flextural rigidity - EI:                          %.3e\n",EI);
    printf("      Section axial rigidity - EA:                              %.3e\n",EA);
    printf("      External Incremental Fx force [DOF(7)]:                   %.3e\n",Fini);
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
void MessageCheck_IMPORT_DATA01(double EI,double EA,int M){
	if (EI <= 0 ||    EA <= 0  ||  M <= 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative or zero data input value is not acceptable ***\n");
	printf("\t\t  Section flextural rigidity - EI:          %.3e\n",EI);
    printf("\t\t  Section axial rigidity - EA:              %.3e\n",EA);
    printf("\t\t  Number of increment:                      %d\n",M);
	Sleep(40000);
	exit(1);
	 }
}
void MessageCheck_IMPORT_DATA03(double TET[],double MOM[],int n){
	int i;
	for(i=0;i<n;i++){
	if (TET[i] < 0|| MOM[i] < 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText03);
	printf("               Row %d has a negative value.\n",i+1);
	printf("                            *** Negative data input value is not acceptable ***\n");
	Sleep(40000);
	exit(1);
	}
	}
	for (i=1;i<n;i++){
    if (TET[i] <= TET[i-1]){
    MessageErrorReportTEXT();
	printf("          Please check the input file! -> [%s]\n",ShowText03);
	printf("          Check row %d value.\n",i+1);
	printf("          Rotation[%d]: %.3e\n",i+1,TET[i]);
	printf("                            *** Data must be sort from minimum value to maximum value ***\n");
	Sleep(40000);exit(1);	
	}
	}
}
void MessageCheckMk(int M,int K){
if (M < 2 || M > STEP || K < 2 || K > 10){
	MessageErrorReportTEXT();
	printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("               Please check this file! -> [%s]\n",ShowText02);	
    printf("        Plastic hinge data: %d - Plastic hinge data must be data -> Minimum : 2 - Maximum : 10\n",K);
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
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_u05[],double output_u06[],double output_u07[],double output_u08[],double output_u09[],double output_u10[],double output_base01[],double output_base02[],double C[STEP][18],int n){
// EXCEL OUTPUT
	int i,I;
	FILE *OutputFile;
	OutputFile = fopen(ShowText05, "w");
fprintf(OutputFile," ### Pushover Analysis Force Analogy Method of Frame with Force Control Based on Euler-Bernoulli Beam Theory ###\n");
fprintf(OutputFile,"Element Number %d\n",1);
fprintf(OutputFile,"\n");
fprintf(OutputFile,"Increment,Base Shear [DOF(1)+DOF(4)],Base Moment [DOF(3)+DOF(6)],Displacement [DOF(7)],Displacement [DOF(8)],Rotation [DOF(9)],Displacement [DOF(10)],Displacement [DOF(11)],Rotation [DOF(12)],Rotation [DOF(13)],Rotation [DOF(14)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i+1,output_base01[i],output_base02[i],output_u01[i],output_u02[i],output_u03[i],output_u04[i],output_u05[i],output_u06[i],output_u07[i],output_u08[i],output_u09[i],output_u10[i]);
fprintf(OutputFile,"\n");
fprintf(OutputFile,"Increment,[1] axial-i,[1] shear-i,[1] moment-i,[1] axial-j,[1] shear-j,[1] moment-j,[2] axial-i,[2] shear-i,[2] moment-i,[2] axial-j,[2] shear-j,[2] moment-j,[3] axial-i,[3] shear-i,[3] moment-i,[3] axial-j,[3] shear-j,[3] moment-j\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i+1,C[i][0],C[i][1],C[i][2],C[i][3],C[i][4],C[i][5],C[i][6],C[i][7],C[i][8],C[i][9],C[i][10],C[i][11],C[i][12],C[i][13],C[i][14],C[i][15],C[i][16],C[i][17]);
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
void IMPORT_DATA01(double &EA,double &EI,double &Fini,int &M){
    double Import_Data[4];
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
	while(i < 5 && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
	EI=Import_Data[0];
	EA=Import_Data[1];
	Fini=Import_Data[2];
	M=Import_Data[3];
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
	while(i < 4 && fgets(line,sizeof(line),InputFile) != NULL);
	n = i;
    //printf("%d\n",n);
}
void IMPORT_DATA03(double TET[],double MOM[],int &k){
	int i = 0;
	FILE *InputFile;
	InputFile = fopen(ShowText03, "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [%s\n",ShowText03);
		Sleep(6000);
		exit(1);
	}
	char line[1000];
	do{
	fscanf(InputFile,"%lf,%lf",&TET[i],&MOM[i]);
	//printf("%d - TET01[%d]: %lf - MOM01[%d]: %lf  TET02[%d]: %lf - MOM02[%d]: %lf\n",i,i,TET01[i],i,MOM01[i],i,TET02[i],i,MOM02[i]);
	i++;
	}
	while(i < 10 && fgets(line,sizeof(line),InputFile) != NULL);
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
void ANALYSIS(double TET[],double MOM[],int np,double EI,double EA,double Fini,double x[],double y[],int M){
   int i,j,z,zMAX,I,check_m,check_u;
   double Product,MomS[6],krk[6],Rk[10],K[NN][NN],Kinv[NN][NN],eleF[Ne][6];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],FF[NN],u[NN];
   double *output_u01 = new double [STEP];
   double *output_u02 = new double [STEP];
   double *output_u03 = new double [STEP];
   double *output_u04 = new double [STEP];
   double *output_u05 = new double [STEP];
   double *output_u06 = new double [STEP];
   double *output_u07 = new double [STEP];
   double *output_u08 = new double [STEP];
   double *output_u09 = new double [STEP];
   double *output_u10 = new double [STEP];
   double *output_base01 = new double [STEP];
   double *output_base02 = new double [STEP];
   double *X = new double [STEP];
   double *Y = new double [STEP];
    double ELE_FORCE[STEP][18];
   double MS[6][6],KG[6][6];
   PlasticHingeStiffnessCOFF(TET,MOM,Rk,np);
   
   // STAGE 01
   MatrixZero(K,NN);
    for (i=0;i<NN;i++)
    u[i] = 0;
    for (i=0;i<STEP;i++){
    output_u01[i] = 0.0;
	output_u02[i] = 0.0;
	output_u03[i] = 0.0;
	output_u04[i] = 0.0;
	output_u05[i] = 0.0;
	output_u06[i] = 0.0;
	output_u07[i] = 0.0;
	output_u08[i] = 0.0;
	output_u09[i] = 0.0;	
	output_u10[i] = 0.0;				
	}

	  for(i=0;i<6*Ne;i++)
	  for(j=0;j<M;j++)
      ELE_FORCE[i][j] = 0.0;
      
      for(i=0;i<6;i++){
      MomS[i]=0.0;krk[i]=0.0;	
	  }

    
    L[0] = SQRT2((x[2]-x[0])*(x[2]-x[0])+(y[2]-y[0])*(y[2]-y[0]));
    L[1] = SQRT2((x[3]-x[1])*(x[3]-x[1])+(y[3]-y[1])*(y[3]-y[1]));
    L[2] = SQRT2((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
    
    lanX[0] = (x[2]-x[0])/L[0];lanY[0] = (y[2]-y[0])/L[0];
    lanX[1] = (x[3]-x[1])/L[1];lanY[1] = (y[3]-y[1])/L[1];
    lanX[2] = (x[3]-x[2])/L[2];lanY[2] = (y[3]-y[2])/L[2];

	  // STAGE 01
	  MatrixZero(K,6);
	  for (z=0;z<M;z++){
	  for (i=6;i<12;i++){
	  u[i] = 0;
	  }	
	  FF[0]=Fini*(z+1);
     
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
    K[0][0]= KG[3][3];//DOF(7)
    K[0][1]= KG[3][4];//DOF(7)
    K[0][2]= KG[3][5];//DOF(7)
    
    K[1][0]= KG[4][3];//DOF(8)
    K[1][1]= KG[4][4];//DOF(8)
    K[1][2]= KG[4][5];//DOF(8)
    
    K[2][0]= KG[5][3];//DOF(9)
    K[2][1]= KG[5][4];//DOF(9)
    K[2][2]= KG[5][5];//DOF(9)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
    K[3][3]= KG[3][3];//DOF(10)
    K[3][4]= KG[3][4];//DOF(10)
    K[3][5]= KG[3][5];//DOF(10)
    
    K[4][3]= KG[4][3];//DOF(11)
    K[4][4]= KG[4][4];//DOF(11)
    K[4][5]= KG[4][5];//DOF(11)
    
    K[5][3]= KG[5][3];//DOF(12)
    K[5][4]= KG[5][4];//DOF(12)
    K[5][5]= KG[5][5];//DOF(12)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
    K[0][0]+= KG[0][0];//DOF(7)
    K[0][1]+= KG[0][1];//DOF(7)
    K[0][2]+= KG[0][2];//DOF(7)
    
    K[1][0]+= KG[1][0];//DOF(8)
    K[1][1]+= KG[1][1];//DOF(8)
    K[1][2]+= KG[1][2];//DOF(8)
    
    K[2][0]+= KG[2][0];//DOF(9)
    K[2][1]+= KG[2][1];//DOF(9)
    K[2][2]+= KG[2][2];//DOF(9)
    
    K[3][3]+= KG[3][3];//DOF(10)
    K[3][4]+= KG[3][4];//DOF(10)
    K[3][5]+= KG[3][5];//DOF(10)
    
    K[4][3]+= KG[4][3];//DOF(11)
    K[4][4]+= KG[4][4];//DOF(11)
    K[4][5]+= KG[4][5];//DOF(11)
    
    K[5][3]+= KG[5][3];//DOF(12)
    K[5][4]+= KG[5][4];//DOF(12)
    K[5][5]+= KG[5][5];//DOF(12)
    
    K[0][3]= KG[0][3];//DOF(7)
    K[0][4]= KG[0][4];//DOF(7)
    K[0][5]= KG[0][5];//DOF(7)
    
    K[1][3]= KG[1][3];//DOF(8)
    K[1][4]= KG[1][4];//DOF(8)
    K[1][5]= KG[1][5];//DOF(8)
  
    K[2][3]= KG[2][3];//DOF(9)
    K[2][4]= KG[2][4];//DOF(9)
    K[2][5]= KG[2][5];//DOF(9)
    
    K[3][0]= KG[3][0];//DOF(10)
    K[3][1]= KG[3][1];//DOF(10)
    K[3][2]= KG[3][2];//DOF(10)     
    
    K[4][0]= KG[4][0];//DOF(11)
    K[4][1]= KG[4][1];//DOF(11)
    K[4][2]= KG[4][2];//DOF(11)
    
    K[5][0]= KG[5][0];//DOF(12)
    K[5][1]= KG[5][1];//DOF(12)
    K[5][2]= KG[5][2];//DOF(12)
     
      MatrixDetermination(K,6,Product);
      if (Product == 0)
      break;
	  MatrixInverse(K,Kinv,6);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,6);
      zMAX = z + 1;

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0,11); // 0: ele 1 - 1: step 1
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1,12); // 0: ele 2 - 1: step 1
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2,13); // 0: ele 3 - 1: step 1
	  ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);

      output_u01[z]=u[0];//output displacement DOF(7)
      output_u02[z]=u[1];//output displacement DOF(8)
      output_u03[z]=u[2];//output displacement DOF(9)
      output_u04[z]=u[3];//output displacement DOF(10)
      output_u05[z]=u[4];//output displacement DOF(11)
      output_u06[z]=u[5];//output displacement DOF(12)
      output_base01[z]=-ELE_FORCE[z][1]-ELE_FORCE[z][7];//output base shear DOF(1)+DOF(4)
      output_base02[z]=-ELE_FORCE[z][2]-ELE_FORCE[z][8];//output base moment DOF(3)+DOF(6)
      
      if (ABS(eleF[0][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][2]),MOM[0]);
      break;
	  }
	  if (ABS(eleF[0][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][5]),MOM[0]);
	  break;
	  }	
	  if (ABS(eleF[1][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][2]),MOM[0]);
      break;
	  }
	  if (ABS(eleF[1][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][5]),MOM[0]);
	  break;
	  }	
	  if (ABS(eleF[2][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][2]),MOM[0]);
      break;
	  }
	  if (ABS(eleF[2][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][5]),MOM[0]);
	  break;
	  }	
	  
	  }// for
      // So plastic hinge has formed in 1 node
    
      // STAGE 02
      MatrixZero(K,7);
	  for (z=zMAX;z<M;z++){
      PlasticHingeStiffness(eleF,0,2,MOM,Rk,MomS,krk,np,0); // 2:Moment DOF(3)
      FF[0]=Fini*(z+1);
	  FF[6]=-MomS[0];
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
    K[0][0]= KG[3][3];//DOF(7)
    K[0][1]= KG[3][4];//DOF(7)
    K[0][2]= KG[3][5];//DOF(7)
    
    K[1][0]= KG[4][3];//DOF(8)
    K[1][1]= KG[4][4];//DOF(8)
    K[1][2]= KG[4][5];//DOF(8)
    
    K[2][0]= KG[5][3];//DOF(9)
    K[2][1]= KG[5][4];//DOF(9)
    K[2][2]= KG[5][5];//DOF(9)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
    K[3][3]= KG[3][3];//DOF(10)
    K[3][4]= KG[3][4];//DOF(10)
    K[3][5]= KG[3][5];//DOF(10)
    
    K[4][3]= KG[4][3];//DOF(11)
    K[4][4]= KG[4][4];//DOF(11)
    K[4][5]= KG[4][5];//DOF(11)
    
    K[5][3]= KG[5][3];//DOF(12)
    K[5][4]= KG[5][4];//DOF(12)
    K[5][5]= KG[5][5];//DOF(12)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
    K[0][0]+= KG[0][0];//DOF(7)
    K[0][1]+= KG[0][1];//DOF(7)
    K[0][2]+= KG[0][2];//DOF(7)
    
    K[1][0]+= KG[1][0];//DOF(8)
    K[1][1]+= KG[1][1];//DOF(8)
    K[1][2]+= KG[1][2];//DOF(8)
    
    K[2][0]+= KG[2][0];//DOF(9)
    K[2][1]+= KG[2][1];//DOF(9)
    K[2][2]+= KG[2][2];//DOF(9)
    
    K[3][3]+= KG[3][3];//DOF(10)
    K[3][4]+= KG[3][4];//DOF(10)
    K[3][5]+= KG[3][5];//DOF(10)
    
    K[4][3]+= KG[4][3];//DOF(11)
    K[4][4]+= KG[4][4];//DOF(11)
    K[4][5]+= KG[4][5];//DOF(11)
    
    K[5][3]+= KG[5][3];//DOF(12)
    K[5][4]+= KG[5][4];//DOF(12)
    K[5][5]+= KG[5][5];//DOF(12)
    
    K[0][3]= KG[0][3];//DOF(7)
    K[0][4]= KG[0][4];//DOF(7)
    K[0][5]= KG[0][5];//DOF(7)
    
    K[1][3]= KG[1][3];//DOF(8)
    K[1][4]= KG[1][4];//DOF(8)
    K[1][5]= KG[1][5];//DOF(8)
  
    K[2][3]= KG[2][3];//DOF(9)
    K[2][4]= KG[2][4];//DOF(9)
    K[2][5]= KG[2][5];//DOF(9)
    
    K[3][0]= KG[3][0];//DOF(10)
    K[3][1]= KG[3][1];//DOF(10)
    K[3][2]= KG[3][2];//DOF(10)     
    
    K[4][0]= KG[4][0];//DOF(11)
    K[4][1]= KG[4][1];//DOF(11)
    K[4][2]= KG[4][2];//DOF(11)
    
    K[5][0]= KG[5][0];//DOF(12)
    K[5][1]= KG[5][1];//DOF(12)
    K[5][2]= KG[5][2];//DOF(12)
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
    K[6][0]= KG[2][3];//DOF(3)
    K[6][1]= KG[2][4];//DOF(3)
    K[6][2]= KG[2][5];//DOF(3)
    K[6][6]= KG[2][2]+krk[0];//DOF(3)
    
    K[0][6]= KG[3][2];//DOF(3)
    K[1][6]= KG[4][2];//DOF(3)
    K[2][6]= KG[5][2];//DOF(3)
      
      //MatrixDisplay(K,7);
      MatrixDetermination(K,7,Product);
      if (Product == 0)
      break;
	  MatrixInverse(K,Kinv,7);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,7);
      zMAX = z + 1;

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0,21); // 0: ele 1 - 1: step 1
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1,22); // 0: ele 2 - 1: step 1
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2,23); // 0: ele 3 - 1: step 1
	  ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);

      output_u01[z]=u[0];//output displacement DOF(7)
      output_u02[z]=u[1];//output displacement DOF(8)
      output_u03[z]=u[2];//output displacement DOF(9)
      output_u04[z]=u[3];//output displacement DOF(10)
      output_u05[z]=u[4];//output displacement DOF(11)
      output_u06[z]=u[5];//output displacement DOF(12)
      output_u07[z]=u[6];//output displacement DOF(3)
      output_base01[z]=-ELE_FORCE[z][1]-ELE_FORCE[z][7];//output base shear DOF(1)+DOF(4)
      output_base02[z]=-ELE_FORCE[z][2]-ELE_FORCE[z][8];//output base moment DOF(3)+DOF(6)
      
//    if (ABS(eleF[0][2]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][2]),MOM[0]);
//      break;
//	  }
	  if (ABS(eleF[0][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][5]),MOM[0]);
	  break;
	  }	
	  if (ABS(eleF[1][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][2]),MOM[0]);
      break;
	  }
	  if (ABS(eleF[1][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][5]),MOM[0]);
	  break;
	  }	
	  if (ABS(eleF[2][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][2]),MOM[0]);
      break;
	  }
	  if (ABS(eleF[2][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][5]),MOM[0]);
	  break;
	  }
	  
	  }// for
      // So plastic hinge has formed in 2 node

      // STAGE 03
//      MatrixZero(K,8);
//	  for (z=zMAX;z<M;z++){	
//      PlasticHingeStiffness(eleF,0,2,MOM,Rk,MomS,krk,np,0); // 2:Moment DOF(3)
//      PlasticHingeStiffness(eleF,0,5,MOM,Rk,MomS,krk,np,1); // 5:Moment DOF(13)
//      PlasticHingeStiffness(eleF,2,2,MOM,Rk,MomS,krk,np,2); // 2:Moment DOF(9)
//      FF[0]=Fini*(z+1);
//	  FF[6]=-MomS[0];
//	  FF[7]=-MomS[1];
//	  FF[2]=-MomS[2];
//    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
//    K[0][0]= KG[3][3];//DOF(7)
//    K[0][1]= KG[3][4];//DOF(7)
//    K[0][2]= KG[3][5];//DOF(7)
//    
//    K[1][0]= KG[4][3];//DOF(8)
//    K[1][1]= KG[4][4];//DOF(8)
//    K[1][2]= KG[4][5];//DOF(8)
//    
//    K[2][0]= KG[5][3];//DOF(9)
//    K[2][1]= KG[5][4];//DOF(9)
//    K[2][2]= KG[5][5];//DOF(9)
//
//    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
//    K[3][3]= KG[3][3];//DOF(10)
//    K[3][4]= KG[3][4];//DOF(10)
//    K[3][5]= KG[3][5];//DOF(10)
//    
//    K[4][3]= KG[4][3];//DOF(11)
//    K[4][4]= KG[4][4];//DOF(11)
//    K[4][5]= KG[4][5];//DOF(11)
//    
//    K[5][3]= KG[5][3];//DOF(12)
//    K[5][4]= KG[5][4];//DOF(12)
//    K[5][5]= KG[5][5];//DOF(12)
//
//    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
//    K[0][0]+= KG[0][0];//DOF(7)
//    K[0][1]+= KG[0][1];//DOF(7)
//    K[0][2]+= KG[0][2];//DOF(7)
//    
//    K[1][0]+= KG[1][0];//DOF(8)
//    K[1][1]+= KG[1][1];//DOF(8)
//    K[1][2]+= KG[1][2];//DOF(8)
//    
//    K[2][0]+= KG[2][0];//DOF(9)
//    K[2][1]+= KG[2][1];//DOF(9)
//    K[2][2]+= KG[2][2];//DOF(9)
//    
//    K[3][3]+= KG[3][3];//DOF(10)
//    K[3][4]+= KG[3][4];//DOF(10)
//    K[3][5]+= KG[3][5];//DOF(10)
//    
//    K[4][3]+= KG[4][3];//DOF(11)
//    K[4][4]+= KG[4][4];//DOF(11)
//    K[4][5]+= KG[4][5];//DOF(11)
//    
//    K[5][3]+= KG[5][3];//DOF(12)
//    K[5][4]+= KG[5][4];//DOF(12)
//    K[5][5]+= KG[5][5];//DOF(12)
//    
//    K[0][3]= KG[0][3];//DOF(7)
//    K[0][4]= KG[0][4];//DOF(7)
//    K[0][5]= KG[0][5];//DOF(7)
//    
//    K[1][3]= KG[1][3];//DOF(8)
//    K[1][4]= KG[1][4];//DOF(8)
//    K[1][5]= KG[1][5];//DOF(8)
//  
//    K[2][3]= KG[2][3];//DOF(9)
//    K[2][4]= KG[2][4];//DOF(9)
//    K[2][5]= KG[2][5];//DOF(9)
//    
//    K[3][0]= KG[3][0];//DOF(10)
//    K[3][1]= KG[3][1];//DOF(10)
//    K[3][2]= KG[3][2];//DOF(10)     
//    
//    K[4][0]= KG[4][0];//DOF(11)
//    K[4][1]= KG[4][1];//DOF(11)
//    K[4][2]= KG[4][2];//DOF(11)
//    
//    K[5][0]= KG[5][0];//DOF(12)
//    K[5][1]= KG[5][1];//DOF(12)
//    K[5][2]= KG[5][2];//DOF(12)
//    
//    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
//    K[6][0]= KG[2][3];//DOF(3)
//    K[6][1]= KG[2][4];//DOF(3)
//    K[6][2]= KG[2][5];//DOF(3)
//    K[6][6]= KG[2][2]+krk[0];//DOF(3)
//    
//    K[0][6]= KG[3][2];//DOF(3)
//    K[1][6]= KG[4][2];//DOF(3)
//    K[2][6]= KG[5][2];//DOF(3)
//    
//    K[7][0]= KG[5][3];//DOF(13)
//    K[7][1]= KG[5][4];//DOF(13)
//    K[7][6]= KG[5][2];//DOF(13)
//    K[7][7]= KG[5][5]+krk[1];//DOF(13)
//    
//    K[0][7]= KG[3][5];//DOF(13)
//    K[1][7]= KG[4][5];//DOF(13)
//    K[6][7]= KG[2][5];//DOF(13)
//    
//    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
//    K[2][2]+= krk[2];//DOF(9)
//    
//      MatrixDisplay(K,8);   
//      MatrixDetermination(K,8,Product);
//      if (Product == 0)
//      break;
//	  MatrixInverse(K,Kinv,8);// Inverse [Kinit]
//      MatrixMulti01(Kinv,FF,u,8);
//      zMAX = z + 1;
//
//      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
//      ElementInternalForce(MS,u,lanX,lanY,eleF,0,31); // 0: ele 1 - 1: step 1
//      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
//      ElementInternalForce(MS,u,lanX,lanY,eleF,1,32); // 0: ele 2 - 1: step 1
//      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
//      ElementInternalForce(MS,u,lanX,lanY,eleF,2,33); // 0: ele 3 - 1: step 1
//	  ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z); 
//      for (i=0;i<8;i++)
//      printf("\t\t %f \n",u[i]);
//      output_u01[z]=u[0];//output displacement DOF(7)
//      output_u02[z]=u[1];//output displacement DOF(8)
//      output_u03[z]=u[2];//output displacement DOF(9)
//      output_u04[z]=u[3];//output displacement DOF(10)
//      output_u05[z]=u[4];//output displacement DOF(11)
//      output_u06[z]=u[5];//output displacement DOF(12)
//      output_u07[z]=u[6];//output displacement DOF(13)
//      output_u08[z]=u[7];//output displacement DOF(9)
//      output_base01[z]=-ELE_FORCE[z][1]-ELE_FORCE[z][7];//output base shear DOF(1)+DOF(4)
//      output_base02[z]=-ELE_FORCE[z][2]-ELE_FORCE[z][8];//output base moment DOF(3)+DOF(6)
//      
//    if (ABS(eleF[0][2]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][2]),MOM[0]);
//      break;
//	  }
//	  if (ABS(eleF[0][5]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,1,ABS(eleF[0][5]),MOM[0]);
//	  break;
//	  }	
//	  if (ABS(eleF[1][2]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][2]),MOM[0]);
//      break;
//	  }
//	  if (ABS(eleF[1][5]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,2,ABS(eleF[1][5]),MOM[0]);
//	  break;
//	  }	
//	  if (ABS(eleF[2][2]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][2]),MOM[0]);
//      break;
//	  }
//	  if (ABS(eleF[2][5]) >= MOM[0]){
//      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,3,ABS(eleF[2][5]),MOM[0]);
//	  break;
//	  }
//
//	  }// for
       //So plastic hinge has formed in 3 node
      
	  
       MessageResult(output_base01,output_base02,output_u01,zMAX);
	    OUTPUT_excel(output_u01,output_u02,output_u03,output_u04,output_u05,output_u06,output_u07,output_u08,output_u09,output_u10,output_base01,output_base02,ELE_FORCE,zMAX);
	    OUTPUT_html(EI,EA,Fini,x,y,TET,MOM,output_u01,output_u02,output_u03,output_u04,output_u05,output_u06,output_u07,output_u08,output_u09,output_u10,output_base01,output_base02,ELE_FORCE,zMAX,np,M);
	    


        char text1[30],text2[30],text3[30];
		for (i=0;i<zMAX;i++){
	    X[i] = output_u01[i];// Disp. DOF(5)
		Y[i] = output_base01[i];// Base Shear DOF(1)
		}
		OUTPUT_HTML_GRAPH(X,Y,zMAX,"Base Shear-Displacement Graph","Displacement [DOF(5)]","Base Shear [DOF(1)+DOF(4)]");

		textcolor(15);
		printf("\n\a - %s -\n",ShowText04);
		system("start /w Graph-outputHTML.html");
		DATE_TIME();
		free(output_u01);free(output_u02);free(output_u03);free(output_u04);
		free(output_u05);free(output_u06);free(output_u07);free(output_u08);
		free(output_u09);free(output_u10);
		free(output_base01);free(output_base02);
		free(X);free(Y);
}
void PlasticHingeStiffnessCOFF(double A[],double B[],double C[],int n){
	int i;
    for (i=0;i<n-1;i++){
    C[i]=(B[i+1]-B[i])/(A[i+1]-A[i]);
    //printf("%d - Rk[%d]: %f\n",i,i,C[i]);
	}
}
void PlasticHingeStiffness(double A[1][6],int I,int II,double B[],double C[],double MomS[],double KRK[],int n,int III){
			int i;
		    for (i=0;i<n-1;i++){
		     if (ABS(A[I][II]) >= B[i] && ABS(A[I][II]) <= B[i+1]){
		     MomS[III] = B[i];KRK[III] = C[i];//printf("\t %.3e    %.3e    %.3e \n",B[i],MomS,KRK);
			   }
			}
		    if  (ABS(A[I][II]) > B[n-1]){
		    MomS[III] = 0;KRK[III] = 0;
			}
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
void ELEMNT_FORCE_OUTPUT(double eleF[Ne][6],double ELE_FORCE[STEP][18],int I){
      int i;
	  for (i=0;i<6;i++)
      ELE_FORCE[I][i]=eleF[0][i];
      for (i=6;i<12;i++)
      ELE_FORCE[I][i]=eleF[1][i-6];	
      for (i=12;i<18;i++)
      ELE_FORCE[I][i]=eleF[2][i-12];
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
	for (i=0;i<98;i++)
	printf("-");
	printf("\n");
    printf("\t    Increment    Base Shear[DOF(1)+DOF(4)]      Base Moment[DOF(3)+DOF(6)]     Displacement [DOF(7)]\n");
	printf("\t   ");
	for (i=0;i<98;i++)
	printf("-");
	printf("\n");
	for (i=0;i<n;i++){
	//Distance(i+1);	
	printf("\t\t %d\t   %.3e\t\t\t     %.3e\t\t    %.3e\t\n",i+1,output_base01[i],output_base01[i],output_u01[i]);	
	}
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
void MessagePlasticHingeTEXT(double TET[],double MOM[],int n){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    int  BB=201,CC=205,DD=187,EE=200,FF=188,GG=186;
    Qa=BB;Qb=CC;Qc=DD;Qd=EE;Qe=FF;Qf=GG;
   	printf("     %c",Qa);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
	printf("     %c       Plastic Hinge Data       %c\n",Qf,Qf);
	printf("     %c     Rotation         Moment    %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<33;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<n;i++)
	printf("          %.3e     %.3e\n",TET[i],MOM[i]);
}
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double TET[],double MOM[],double output_u01[],double output_u02[],double output_u03[],double output_u04[],double output_u05[],double output_u06[],double output_u07[],double output_u08[],double output_u09[],double output_u10[],double output_base01[],double output_base02[],double C[STEP][18],int n,int m,int M){
// HTML OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText06, "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"PushoverForceAnalogyMethodBerFrameFC-image01.png\" style=\"width:1000px ; height:500px\" alt=\"analysis\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th bgcolor=\"cyan\">  Pushover Analysis Force Analogy Method of Frame with Force Control Based on Euler-Bernoulli Beam Theory - Output Report </th> \n");
    // TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Input Data </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section flextural rigidity - EI: </th><th> %.3e </th> </tr>\n",EI);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section axial rigidity - EA: </th><th> %.3e </th> </tr>\n",EA);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External Incremental Fx force [DOF(7)]: </th><th> %.3e </th> </tr>\n",Fini);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Number of increment: </th><th> %d </th> </tr>\n",M);   
    // TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th colspan=\"2\" bgcolor=\"orange\"> Hinges Data </th> \n");
	fprintf(OutputFile,"<tr> <th colspan=\"2\" bgcolor=\"orange\">Hinge[1]: Moment - Rotation </th></tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Rotation  </th> <th bgcolor=\"orange\"> Moment </th> </tr>\n");
	for(i=0;i<m;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> </tr>\n",TET[i],MOM[i]);
    }
    // TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"3\" bgcolor=\"orange\"> Structral Coordinate </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Node Number </th><th bgcolor=\"orange\">X Coordinate</th> <th bgcolor=\"orange\">Y Coordinate</th> </tr>\n");
	for(i=0;i<Ne+1;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td>\n",i+1,X[i],Y[i]);
    }
	// TABLE 4
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"3000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"17\" bgcolor=\"orange\"> Structral Deformation </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Base Shear [DOF(1)+DOF(4)]</th><th bgcolor=\"orange\">Base Moment [DOF(3)+DOF(6)]</th><th bgcolor=\"orange\">Displacement [DOF(7)]</th> <th bgcolor=\"orange\">Displacement [DOF(8)]</th><th bgcolor=\"orange\">Rotation [DOF(9)]</th><th bgcolor=\"orange\">Displacement [DOF(10)]</th> <th bgcolor=\"orange\">Displacement [DOF(11)]</th><th bgcolor=\"orange\">Rotation [DOF(12)]</th><th bgcolor=\"orange\">Rotation [DOF(13)]</th><th bgcolor=\"orange\">Rotation [DOF(14)]</th></tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td></tr>\n",i+1,output_base01[i],output_base02[i],output_u01[i],output_u02[i],output_u03[i],output_u04[i],output_u05[i],output_u06[i],output_u07[i],output_u08[i],output_u09[i],output_u10[i]);
    }
    // TABLE 5
    fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"3000px\" height=\"120px\" bgcolor=\"yellow\">\n");
    fprintf(OutputFile,"<tr><th colspan=\"19\" bgcolor=\"orange\"> Structral Element Internal Forces </th> </tr>\n");
    fprintf(OutputFile,"<tr><th colspan=\"7\" bgcolor=\"orange\"> Element[1] </th> <th colspan=\"7\" bgcolor=\"orange\"> Element[2] </th><th colspan=\"7\" bgcolor=\"orange\"> Element[3] </th></tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Axial-i</th> <th bgcolor=\"orange\">Shear-i</th> <th bgcolor=\"orange\">Moment-i</th> <th bgcolor=\"orange\">Axial-j</th> <th bgcolor=\"orange\">Shear-j</th> <th bgcolor=\"orange\">Moment-j</th> <th bgcolor=\"orange\">Axial-i</th> <th bgcolor=\"orange\">Shear-i</th> <th bgcolor=\"orange\">Moment-i</th> <th bgcolor=\"orange\">Axial-j</th> <th bgcolor=\"orange\">Shear-j</th> <th bgcolor=\"orange\">Moment-j</th><th bgcolor=\"orange\">Axial-i</th> <th bgcolor=\"orange\">Shear-i</th> <th bgcolor=\"orange\">Moment-i</th> <th bgcolor=\"orange\">Axial-j</th> <th bgcolor=\"orange\">Shear-j</th> <th bgcolor=\"orange\">Moment-j</th></tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td></tr>\n",i+1,C[i][0],C[i][1],C[i][2],C[i][3],C[i][4],C[i][5],C[i][6],C[i][7],C[i][8],C[i][9],C[i][10],C[i][11],C[i][12],C[i][13],C[i][14],C[i][15],C[i][16],C[i][17]);
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
int MessageCheckMoment(double eleF[Ne][6],int z,double MOM[],int &I){
	int i;
  for (i=0;i<Ne;i++){
      if (ABS(eleF[i][2]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-i: %f reaached to yield: %f \n",z+1,i+1,ABS(eleF[i][2]),MOM[0]);
      I = 1;
	  }
	  if (ABS(eleF[i][5]) >= MOM[0]){
      printf("\t Increment %d - Element %d - Moment-j: %f reaached to yield: %f \n",z+1,i+1,ABS(eleF[i][5]),MOM[0]);
	  I = 1;
	  }		
                   }
      return I;
}
int MessageCheckRotation(double u[],int z,double TET[],int np,int &I){
      if (ABS(u[z]) >= TET[np-1]){
      printf("\t Increment %d - Rotation: %f reached to ultimate: %f \n",z+1,ABS(u[z]),TET[np-1]);
      I = 1;
	  }	
      return I;
}
void MatrixZero(double A[NN][NN],int n){
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++){
	A[i][j]	= 0.0;
	  }
}
void MatrixDisplay(double A[NN][NN],int n){
	for (int i=0;i<n;i++){
	for (int j=0;j<n;j++)
	printf("  %.3e",A[i][j]);
	printf("\n");
	  }
	  printf("\n\n\n");
}
