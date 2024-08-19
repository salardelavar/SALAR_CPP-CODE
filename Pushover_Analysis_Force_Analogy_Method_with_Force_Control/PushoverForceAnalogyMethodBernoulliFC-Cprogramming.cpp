#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 6 // Degree of freedom
#define Ne 1    // number of element
#define N 2    // number of node
#define STEP 5000   // number of node
#define ShowText01 "PushoverForceAnalogyMethodBernoulliFC-inputDATA.csv"
#define ShowText02 "PushoverForceAnalogyMethodBernoulliFC-inputCOORDINATE.csv"
#define ShowText03 "PushoverForceAnalogyMethodBernoulliFC-inputHINGE.csv"
#define ShowText04 "Output data is written in Excel and Html file"
#define ShowText05 "PushoverForceAnalogyMethodBernoulliFC-outputEXCEL.csv"
#define ShowText06 "PushoverForceAnalogyMethodBernoulliFC-outputHTML.html"
#define ShowText07 "Graph-outputHTML.html"

void IMPORT_DATA01(double &EA,double &EI,double &Fini,int &M);
void IMPORT_DATA02(double x[],double y[],int &n);
void IMPORT_DATA03(double TET01[],double MOM01[],double TET02[],double MOM02[],int &k);
void Matrix_Zero(double A[][NN],int n);
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
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
void MessageInitialData(double EI,double EA,double Fini,double x[],double y[],int M);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double EI,double EA,int M);
void MessageCheck_IMPORT_DATA03(double TET01[],double MOM01[],double TET02[],double MOM02[],int n);
void MessageCheckMk(int M,int K);
void MessageStrCoorTEXT(double X[],double Y[],int n);
void MessagePlasticHingeTEXT(double TET01[],double MOM01[],double TET02[],double MOM02[],int n);
void MessageResult(double output_base01[],double output_base02[],double output_u01[],double output_u02[],double output_u03[],int n);
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_base01[],double output_base02[],double C[6][STEP],int n);
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double TET01[],double MOM01[],double TET02[],double MOM02[],double output_u01[],double output_u02[],double output_u03[],double output_base01[],double output_base02[],double C[6][STEP],int n,int m,int M);
void ANALYSIS(double TET01[],double MOM01[],double TET02[],double MOM02[],int np,double EI,double EA,double Fini,double x[],double y[],int M);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
double MAX_ABS(double A[],int n);
void PlasticHingeStiffnessCOFF(double A[],double B[],double C[],int n);
void PlasticHingeStiffness(double A[1][6],int I,int II,double B[],double C[],double &MomS,double &KRK,int n);

int main(){
    double EI,EA,Fini,x[2],y[2],TET01[10],MOM01[10],TET02[10],MOM02[10];
    int n,np,m,M;
    IMPORT_DATA01(EA,EI,Fini,M);
    IMPORT_DATA02(x,y,n);
    IMPORT_DATA03(TET01,MOM01,TET02,MOM02,np);
    MessageCheck_IMPORT_DATA01(EI,EA,M);
    MessageCheck_IMPORT_DATA03(TET01,MOM01,TET02,MOM02,np);
    MessageCheckMk(M,np);
	textcolor(14);
    MessageInitialData(EI,EA,Fini,x,y,M);
    MessageStrCoorTEXT(x,y,n);
    MessagePlasticHingeTEXT(TET01,MOM01,TET02,MOM02,np);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(TET01,MOM01,TET02,MOM02,np,EI,EA,Fini,x,y,M);
    getch();
    return 0;
}
void Matrix_Zero(double A[][NN],int n){
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
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[0];UU[4]=0;UU[5]=0;
	}
	if (II == 2){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[0];UU[4]=0;UU[5]=U[1];
	}
	if (II == 3){
	UU[0]=0;UU[1]=0;UU[2]=U[2];UU[3]=U[0];UU[4]=0;UU[5]=U[1];
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
	for (i=1;i<60;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c   >> In the name of God, the Gracious, the Merciful <<    %c\n",Qf,Qf);
    printf("\t\t\t\t%c Pushover Analysis Force Analogy Method with Force Control %c\n",Qf,Qf);
    printf("\t\t\t\t%c             Based on Euler-Bernoulli Beam Theory          %c\n",Qf,Qf);
    printf("\t\t\t\t%c                     UNIT: Free Unit                       %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<60;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);
	printf("\t\t\t\t%c    This program is written by Salar Delavar Ghashghaei    %c\n",Qf,Qf);
	printf("\t\t\t\t%c             E-mail: salar.d.ghashghaei@gmail.com          %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<60;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

    MessageInputDataTEXT();
    printf("      Section flextural rigidity - EI:                          %.3e\n",EI);
    printf("      Section axial rigidity - EA:                              %.3e\n",EA);
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
void MessageCheck_IMPORT_DATA03(double TET01[],double MOM01[],double TET02[],double MOM02[],int n){
	int i;
	for(i=0;i<n;i++){
	if (TET01[i] < 0|| MOM01[i] < 0 || TET02[i] < 0|| MOM02[i] < 0){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText03);
	printf("               Row %d has a negative value.\n",i+1);
	printf("                            *** Negative data input value is not acceptable ***\n");
	Sleep(40000);
	exit(1);
	}
	}
	for (i=1;i<n;i++){
    if (TET01[i] <= TET01[i-1] || TET02[i] <= TET02[i-1]){
    MessageErrorReportTEXT();
	printf("          Please check the input file! -> [%s]\n",ShowText03);
	printf("          Check row %d value.\n",i+1);
	printf("          Rotation[%d]: %.3e - Rotation[%d]: %.3e\n",i+1,TET01[i],i+1,TET02[i]);
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
void OUTPUT_excel(double output_u01[],double output_u02[],double output_u03[],double output_base01[],double output_base02[],double C[6][STEP],int n){
// EXCEL OUTPUT
	int i,I;
	FILE *OutputFile;
	OutputFile = fopen(ShowText05, "w");
fprintf(OutputFile," ### Pushover Analysis Force Analogy Method with Force Control Based on Euler-Bernoulli Beam Theory ###\n");
fprintf(OutputFile,"Element Number %d\n",1);
fprintf(OutputFile,"\n");
fprintf(OutputFile,"Increment,Base Shear [DOF(1)],Base Moment [DOF(3)],Displacement [DOF(4)],Rotation [DOF(3)],Rotation [DOF(6)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e\n",i+1,output_base01[i],output_base02[i],output_u01[i],output_u03[i],output_u02[i]);
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
	while(i < 2 && fgets(line,sizeof(line),InputFile) != NULL);
	n = i;
    //printf("%d\n",n);
}
void IMPORT_DATA03(double TET01[],double MOM01[],double TET02[],double MOM02[],int &k){
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
	fscanf(InputFile,"%lf,%lf,%lf,%lf",&TET01[i],&MOM01[i],&TET02[i],&MOM02[i]);
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
void ANALYSIS(double TET01[],double MOM01[],double TET02[],double MOM02[],int np,double EI,double EA,double Fini,double x[],double y[],int M){
   int i,j,z,zMAX,I;
   double Product,MomS01,MomS02,krk01,krk02,Rk01[10],Rk02[10],K[NN][NN],Kinv[NN][NN],eleF[Ne][6];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],FF[NN],u[NN];
   double *output_u01 = new double [STEP];
   double *output_u02 = new double [STEP];
   double *output_u03 = new double [STEP];
   double *output_base01 = new double [STEP];
   double *output_base02 = new double [STEP];
   double *X = new double [STEP];
   double *Y = new double [STEP];
    double ELE_FORCE[6][STEP];
   double MS[6][6],KG[6][6];
   PlasticHingeStiffnessCOFF(TET01,MOM01,Rk01,np);
   PlasticHingeStiffnessCOFF(TET02,MOM02,Rk02,np);
   Matrix_Zero(K,3);
    for (i=0;i<3;i++)
    u[i] = 0;
    for (i=0;i<STEP;i++){
    output_u01[i] = 0.0;
	output_u02[i] = 0.0;
	output_u03[i] = 0.0;	
	}

	  for(i=0;i<6*Ne;i++)
	  for(j=0;j<M;j++)
      ELE_FORCE[i][j] = 0;


    for (i=0;i<Ne;i++){
    L[i] = SQRT2((x[i+1]-x[i])*(x[i+1]-x[i])+(y[i+1]-y[i])*(y[i+1]-y[i]));
	lanX[i] = (x[i+1]-x[i])/L[i];lanY[i] = (y[i+1]-y[i])/L[i];
	}
	  MomS01=0;krk01=0;
      MomS02=0;krk02=0;
	  // STAGE 01
	  for (z=0;z<M;z++){
	  u[1]=0;u[2]=0;
	  FF[0]=Fini*(z+1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1); // 0: ele 1 - 1: stiffness matrix: 1
      K[0][0]= KG[3][3];//DOF(4)
      MatrixDetermination(K,1,Product);
      if (Product == 0)
      break;
	  MatrixInverse(K,Kinv,1);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,1);
      zMAX = z + 1;


      ElementInternalForce(MS,u,lanX,lanY,eleF,0,1); // 0: ele 1 - 1: step 1

	  ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);
//	  for (i=0;i<1;i++)
//	  printf("\t  u[%d]: %f \n",i,u[i]);
//	  for (i=0;i<6;i++)
//	  printf("\t %f \n",ELE_FORCE[i][z]);


      output_u01[z]=u[0];//output displacement DOF(4)
      output_base01[z]=-ELE_FORCE[1][z];//output base shear DOF(1)
      output_base02[z]=-ELE_FORCE[2][z];//output base moment DOF(3)
      
      if (ABS(eleF[0][2]) >= MOM01[0]){
      printf("\t Moment-1: %f reaached to yield: %f \n",ABS(eleF[0][2]),MOM01[0]);
	  break;
	  }
	  if (ABS(eleF[0][5]) >= MOM02[0]){
      printf("\t Moment-2: %f reaached to yield: %f \n",ABS(eleF[0][5]),MOM02[0]);
	  break;
	  }	
	  if (ABS(u[1]) >= TET01[np-1]){
      printf("\t Rotation-1: %f reaached to yield: %f \n",ABS(u[1]),TET01[np-1]);
	  break;
	  }
	  if (ABS(u[2]) >= TET02[np-1]){
      printf("\t Rotation-2: %f reaached to yield: %f \n",ABS(u[2]),TET02[np-1]);
	  break;
	  }
	  
	  }// for
      // So plastic hinge has formed in 2 node
    
      // STAGE 02
	  for (z=zMAX;z<M;z++){
	  u[2]=0;
      PlasticHingeStiffness(eleF,0,5,MOM02,Rk02,MomS02,krk02,np); // 5:Moment DOF(6)
      //printf("\t\t\t %.3e    %.3e \n",MomS02,krk02);

	  FF[0]=Fini*(z+1);
	  FF[1]=-MomS02;
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1); // 0: ele 1 - 1: stiffness matrix: 1

      K[0][0]= KG[3][3];//DOF(4)
      K[0][1]= KG[3][5];//DOF(4)
      K[1][0]= KG[5][3];//DOF(6)
      K[1][1]= KG[5][5]+krk02;//DOF(6)
//      for (i=0;i<2;i++){
//      for (j=0;j<2;j++)
//      printf("\t  K[%d][%d]: %f ",i,j,K[i][j]);
//	  printf("\n");
//	  }
//      printf("\n");
      MatrixDetermination(K,2,Product);
      if (Product == 0)
      break;
	  MatrixInverse(K,Kinv,2);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,2);

      zMAX = z + 1;

      ElementInternalForce(MS,u,lanX,lanY,eleF,0,2); // 2: step 2

//	  for (i=0;i<2;i++)
//	  printf("\t  u[%d]: %f \n",i,u[i]);
//	  for (i=0;i<6;i++)
//	  printf("\t %f \n",eleF[0][i]);
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);
      output_u01[z]=u[0];//output displacement DOF(4)
      output_u02[z]=u[1];//output rotation DOF(6)
      output_base01[z]=-ELE_FORCE[1][z];//output base shear DOF(1)
      output_base02[z]=-ELE_FORCE[2][z];//output base moment DOF(3)
      
      if (ABS(eleF[0][2]) >= MOM01[0]){
      printf("\t Moment-1: %f reaached to yield: %f \n",ABS(eleF[0][2]),MOM01[0]);
	  break;
	  }
//	  if (ABS(eleF[0][5]) >= MOM02[0]){
//      printf("\t Moment-2: %f reaached to yield: %f \n",ABS(eleF[0][5]),MOM02[0]);
//	  break;
//	  }	
	  if (ABS(u[1]) >= TET01[np-1]){
      printf("\t Rotation-1: %f reaached to yield: %f \n",ABS(u[1]),TET01[np-1]);
	  break;
	  }
	  if (ABS(u[2]) >= TET02[np-1]){
      printf("\t Rotation-2: %f reaached to yield: %f \n",ABS(u[2]),TET02[np-1]);
	  break;
	  }
	  
	  }// for
      // So plastic hinge has formed in 1 node

      // STAGE 03
	  for (z=zMAX;z<M;z++){	
      PlasticHingeStiffness(eleF,0,5,MOM02,Rk02,MomS02,krk02,np); // 5:Moment DOF(6)
	  PlasticHingeStiffness(eleF,0,2,MOM01,Rk01,MomS01,krk01,np); // 2:Moment DOF(3)
      //printf("\t\t\t %.3e    %.3e \n",MomS01,krk01);
	  FF[0]=Fini*(z+1);
	  FF[1]=-MomS02;
	  FF[2]=-MomS01;
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1); // 0: ele 1 - 1: stiffness matrix: 1

      K[0][0]= KG[3][3];//DOF(4)
      K[0][1]= KG[3][5];//DOF(4)
      K[0][2]= KG[3][2];//DOF(4)
      K[1][0]= KG[5][3];//DOF(6)
      K[1][1]= KG[5][5]+krk02;//DOF(6)
      K[1][2]= KG[5][2];//DOF(6)
      K[2][0]= KG[2][3];//DOF(3)
      K[2][1]= KG[2][5];//DOF(3)
      K[2][2]= KG[2][2]+krk01;//DOF(3)
//      for (i=0;i<3;i++){
//      for (j=0;j<3;j++)
//      printf("\t  K[%d][%d]: %f ",i,j,K[i][j]);
//	  printf("\n");
//	  }
//      MatrixDetermination(K,3,Product);
//      if (Product == 0)
//      break;
	  MatrixInverse(K,Kinv,3);// Inverse [Kinit]
      MatrixMulti01(Kinv,FF,u,3);
      zMAX = z + 1;
      ElementInternalForce(MS,u,lanX,lanY,eleF,0,3);// 3: step 3
//	  for (i=0;i<3;i++)
//	  printf("\t  u[%d]: %f \n",i,u[i]);
//	  for (i=0;i<6;i++)
//	  printf("\t %f \n",eleF[0][i]);
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,z);
      output_u01[z]=u[0];//output displacement DOF(4)
      output_u02[z]=u[1];//output rotation DOF(6)
      output_u03[z]=u[2];//output rotation DOF(3)
      output_base01[z]=-ELE_FORCE[1][z];//output base shear DOF(1)
      output_base02[z]=-ELE_FORCE[2][z];//output base moment DOF(3)
      
//      if (ABS(eleF[0][2]) >= MOM01[0]){
//      printf("\t Moment-1: %f reaached to yield: %f \n",ABS(eleF[0][2]),MOM01[0]);
//	  break;
//	  }
//	  if (ABS(eleF[0][5]) >= MOM02[0]){
//      printf("\t Moment-2: %f reaached to yield: %f \n",ABS(eleF[0][5]),MOM02[0]);
//	  break;
//	  }	
	  if (ABS(u[1]) >= TET01[np-1]){
      printf("\t Rotation-1: %f reaached to yield: %f \n",ABS(u[1]),TET01[np-1]);
	  break;
	  }
	  if (ABS(u[2]) >= TET02[np-1]){
      printf("\t Rotation-2: %f reaached to yield: %f \n",ABS(u[2]),TET02[np-1]);
	  break;
	  }
	  }// for

       MessageResult(output_base01,output_base02,output_u01,output_u02,output_u03,zMAX);
	    OUTPUT_excel(output_u01,output_u02,output_u03,output_base01,output_base02,ELE_FORCE,zMAX);
	    OUTPUT_html(EI,EA,Fini,x,y,TET01,MOM01,TET02,MOM02,output_u01,output_u02,output_u03,output_base01,output_base02,ELE_FORCE,zMAX,np,M);
	    
        char text1[30]="Base Shear-Displacement Graph",text2[30]="Displacement [DOF(5)]",text3[30]="Base Shear [DOF(1)]";
		for (i=0;i<zMAX;i++){
	    X[i] = output_u01[i];// Disp. DOF(5)
		Y[i] = output_base01[i];// Base Shear DOF(1)
		}
		OUTPUT_HTML_GRAPH(X,Y,zMAX,text1,text2,text3);
 
		textcolor(15);
		printf("\n\a - %s -\n",ShowText04);
		system("start /w Graph-outputHTML.html");
		DATE_TIME();
		free(output_u01);free(output_u02);free(output_u03);
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
void PlasticHingeStiffness(double A[1][6],int I,int II,double B[],double C[],double &MomS,double &KRK,int n){
			int i;
		    for (i=0;i<n-1;i++){
		     if (ABS(A[I][II]) >= B[i] && ABS(A[I][II]) <= B[i+1]){
		     MomS = B[i];KRK = C[i];//printf("\t %.3e    %.3e    %.3e \n",B[i],MomS,KRK);
			   }
			}
		    if  (ABS(A[I][II]) > B[n-1]){
		    MomS = 0;KRK = 0;
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
void MessageResult(double output_base01[],double output_base02[],double output_u01[],double output_u02[],double output_u03[],int n){
	int i;
	printf("\t   ");
	for (i=0;i<120;i++)
	printf("-");
	printf("\n");
    printf("\t    Increment    Base Shear[DOF(1)]      Base Moment[DOF(3)]     Disp. [DOF(5)]    Rotation [DOF(3)]     Rotation [DOF(6)]\n");
	printf("\t   ");
	for (i=0;i<120;i++)
	printf("-");
	printf("\n");
	for (i=0;i<n;i++){
	//Distance(i+1);	
	printf("\t\t %d\t   %.3e\t          %.3e\t          %.3e\t     %.3e\t          %.3e\n",i+1,output_base01[i],output_base01[i],output_u01[i],output_u03[i],output_u02[i]);	
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
void MessagePlasticHingeTEXT(double TET01[],double MOM01[],double TET02[],double MOM02[],int n){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    int  BB=201,CC=205,DD=187,EE=200,FF=188,GG=186;
    Qa=BB;Qb=CC;Qc=DD;Qd=EE;Qe=FF;Qf=GG;
   	printf("     %c",Qa);
	for (i=1;i<66;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
	printf("     %c                        Plastic Hinge Data                       %c\n",Qf,Qf);
	printf("     %c         Plastic Hinge 01                Plastic Hinge 02        %c\n",Qf,Qf);
	printf("     %c       Rotation     Moment             Rotation        Moment    %c\n",Qf,Qf);
	printf("     %c",Qd);
	for (i=1;i<66;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);
	for(i=0;i<n;i++)
	printf("          %.3e     %.3e         %.3e     %.3e\n",TET01[i],MOM01[i],TET02[i],MOM02[i]);
}
void OUTPUT_html(double EI,double EA,double Fini,double X[],double Y[],double TET01[],double MOM01[],double TET02[],double MOM02[],double output_u01[],double output_u02[],double output_u03[],double output_base01[],double output_base02[],double C[6][STEP],int n,int m,int M){
// HTML OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText06, "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"PushoverForceAnalogyMethodBernoulliFC-image01.png\" style=\"width:1000px ; height:500px\" alt=\"analysis\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th bgcolor=\"cyan\">  Pushover Analysis Force Analogy Method with Force Control Based on Euler-Bernoulli Beam Theory - Output Report </th> \n");
    // TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Input Data </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section flextural rigidity - EI: </th><th> %.3e </th> </tr>\n",EI);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section axial rigidity - EA: </th><th> %.3e </th> </tr>\n",EA);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External Incremental Fx force [DOF(4)]: </th><th> %.3e </th> </tr>\n",Fini);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Number of increment: </th><th> %d </th> </tr>\n",M);   
    // TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th colspan=\"4\" bgcolor=\"orange\"> Hinges Data </th> \n");
	fprintf(OutputFile,"<tr> <th colspan=\"2\" bgcolor=\"orange\">Hinge[1]: Moment - Rotation </th> <th colspan=\"2\" bgcolor=\"orange\">Hinge[2]: Moment - Rotation </th></tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Rotation  </th> <th bgcolor=\"orange\"> Moment </th><th bgcolor=\"orange\"> Rotation </th> <th bgcolor=\"orange\"> Moment</th> </tr>\n");
	for(i=0;i<m;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> </tr>\n",TET01[i],MOM01[i],TET02[i],MOM02[i]);
    }
    // TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"3\" bgcolor=\"orange\"> Structral Coordinate </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Node Number </th><th bgcolor=\"orange\">X Coordinate</th> <th bgcolor=\"orange\">Y Coordinate</th> </tr>\n");
	for(i=0;i<Ne+1;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td>\n",i+1,X[i],Y[i]);
    }
	// TABLE 4
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"6\" bgcolor=\"orange\"> Structral Deformation </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Increment</th> <th bgcolor=\"orange\">Base Shear [DOF(1)]</th><th bgcolor=\"orange\">Base Moment [DOF(3)]</th><th bgcolor=\"orange\">Displacement [DOF(5)]</th> <th bgcolor=\"orange\">Rotation [DOF(3)]</th><th bgcolor=\"orange\">Rotation [DOF(6)]</th></tr>\n");
	for(i=0;i<n;i++){
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td></tr>\n",i+1,output_base01[i],output_base02[i],output_u01[i],output_u03[i],output_u02[i]);
    }
    // TABLE 5
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
