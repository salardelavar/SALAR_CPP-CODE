// brace internal force for base shear must be correct and add lanX anad lanY for to brace element
#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 6 // Degree of freedom
#define Ne 5    // number of element
#define ShowText01 "Graph-outputHTML.html"
#define ShowText02 "                       >> IN THE NAME OF ALLAH <<                   "
#define ShowText03 "  Nonlinear Analysis of Frame and Brace with Hinge by Hinge Method  "
#define ShowText04 "                             UNIT: Free Unit                        "
#define ShowText05 "           This program is written by Salar Delavar Qashqai         "
#define ShowText06 "                  E-mail: salar.d.ghashghaei@gmail.com              "
void IMPORT_DATA01(double &Length,double Height[],double &EA,double &EI,double &Ultimate_Moment);
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
void MatrixDetermination(double [][NN],int );
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [],int );
void Matrix_Transpose(double A[6][6],double B[6][6]);
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]);
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I);// Calculate internal element force
void ELEMNT_FORCE_OUTPUT(double eleF[3][6],double ELE_FORCE[3][30],double SUM_ELE_FORCE[30],double u[],double sum_u[], int I);
double ABS(double);
double MAX_ABS(double A[],int n);
double MIN(double A[],int n);
double SQRT2(double D);
void MessageInitialData(double L,double H[],double EI,double EA,double Ultimate_Moment);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double L,double Height[],double EI,double EA,double Ultimate_Moment);
void MessageResult(double A[],double B[5][6],int n);
void OUTPUT_excel(double A[4][6],double B[4],double C[4][30],int n);
void ANALYSIS(double Length,double Height[],double EI,double EA,double Ultimate_Moment);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]);
int main(){   
    double Length,EI,EA,Height[2],Ultimate_Moment;
    
    IMPORT_DATA01(Length,Height,EA,EI,Ultimate_Moment);
    MessageCheck_IMPORT_DATA01(Length,Height,EI,EA,Ultimate_Moment);
	textcolor(14);
    MessageInitialData(Length,Height,EI,EA,Ultimate_Moment);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(Length,Height,EI,EA,Ultimate_Moment);
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
    int II;
	lan[0][0]=lanX[I];lan[0][1]=lanY[I];lan[0][2]=0;lan[0][3]=0;lan[0][4]=0;lan[0][5]=0;
	lan[1][0]=-lanY[I];lan[1][1]=lanX[I];lan[1][2]=0;lan[1][3]=0;lan[1][4]=0;lan[1][5]=0;
	lan[2][0]=0;lan[2][1]=0;lan[2][2]=1;lan[2][3]=0;lan[2][4]=0;lan[2][5]=0;
	lan[3][0]=0;lan[3][1]=0;lan[3][2]=0;lan[3][3]=lanX[I];lan[3][4]=lanY[I];lan[3][5]=0;
	lan[4][0]=0;lan[4][1]=0;lan[4][2]=0;lan[4][3]=-lanY[I];lan[4][4]=lanX[I];lan[4][5]=0;
	lan[5][0]=0;lan[5][1]=0;lan[5][2]=0;lan[5][3]=0;lan[5][4]=0;lan[5][5]=1;

	if (I == 0){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[0];UU[4]=U[1];UU[5]=U[2];
	}
	if (I == 1){
	UU[0]=0;UU[1]=0;UU[2]=0;UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	if (I == 2){
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[2];UU[3]=U[3];UU[4]=U[4];UU[5]=U[5];
	}
	int i,j;
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

void MessageInitialData(double L,double H[],double EI,double EA,double Ultimate_Moment){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText02,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText03,Qf);
    printf("\t\t\t\t%c%s%c\n",Qf,ShowText04,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);   
	printf("\t\t\t\t%c%s%c\n",Qf,ShowText05,Qf);
	printf("\t\t\t\t%c%s%c\n",Qf,ShowText06,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<69;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

    
    MessageInputDataTEXT();
    printf("      Length of Frame:                                          %.3e\n",L);
    printf("      Height of Column 1:                                       %.3e\n",H[0]);
	printf("      Height of Column 2:                                       %.3e\n",H[1]);
    printf("      Section flextural rigidity - EI:                          %.3e\n",EI);
    printf("      Section axial rigidity - EA:                              %.3e\n",EA);
    printf("      Section ultimate capacity moment:                         %.3e\n",Ultimate_Moment);

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
void MessageCheck_IMPORT_DATA01(double L,double Height[],double EI,double EA,double Ultimate_Moment){
	if ( L < 0 || Height[0] < 0 || Height[1] < 0 ||  EI < 0 ||    EA < 0 ||    Ultimate_Moment < 0 ){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [ PushoverHingeByHingeMethodFrameBrace-inputDATA.csv ]\n");
	printf("                            *** Negative data input value is not acceptable ***\n");
	printf("  Length of Frame:                     %.3e\n",L);
	printf("  Height of Column 1:                  %.3e\n",Height[0]);
	printf("  Height of Column 2:                  %.3e\n",Height[1]);
	printf("  Section flextural rigidity - EI:     %.3e\n",EI);
    printf("  Section axial rigidity - EA:         %.3e\n",EA);
    printf("  Section ultimate capacity moment:    %.3e\n",Ultimate_Moment);
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
void OUTPUT_excel(double A[4][6],double B[4],double C[4][30],int n){
// EXCEL OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen("PushoverHingeByHingeMethodFrameBrace-outputEXCEL.csv", "w");
fprintf(OutputFile," ###  Nonlinear Analysis of Frame and Brace with Hinge by Hinge Method  ###\n");
fprintf(OutputFile,"Increment,Base Shear[DOF(1)]+[DOF(4)],Displacement [DOF(7)],Displacement [DOF(8)],Rotation [DOF(9)],Displacement [DOF(10)],Displacement [DOF(11)],Rotation [DOF(12)],Ele.1 [DOF(1)],Ele.1 [DOF(2)],Ele.1 [DOF(3)],Ele.1 [DOF(7)],Ele.1 [DOF(8)],Ele.1 [DOF(9)],Ele.2 [DOF(4)],Ele.2 [DOF(5)],Ele.2 [DOF(6)],Ele.2 [DOF(10)],Ele.2 [DOF(11)],Ele.2 [DOF(12)],Ele.3 [DOF(7)],Ele.3 [DOF(8)],Ele.3 [DOF(9)],Ele.3 [DOF(10)],Ele.3 [DOF(11)],Ele.3 [DOF(12)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i+1,B[i],A[i][0],A[i][1],A[i][2],A[i][3],A[i][4],A[i][5],C[i][0],C[i][1],C[i][2],C[i][3],C[i][4],C[i][5],C[i][6],C[i][7],C[i][8],C[i][9],C[i][10],C[i][11],C[i][12],C[i][13],C[i][14],C[i][15],C[i][16],C[i][17]); 
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
void IMPORT_DATA01(double &Length,double Height[],double &EA,double &EI,double &Ultimate_Moment){
    double Import_Data[6];
    int i=0;
	FILE *InputFile;
	InputFile = fopen("PushoverHingeByHingeMethodFrameBrace-inputDATA.csv", "r");
	if (!InputFile){
		MessageErrorReportTEXT();
		printf("          File is not available! -> [PushoverHingeByHingeMethodFrameBrace-inputDATA.csv] \n");
		Sleep(6000);
		exit(1);
	}
	char line[100],a[100];
	while(i < 6 && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
	Length=Import_Data[0];
	Height[0]=Import_Data[1];
	Height[1]=Import_Data[2];
	EI=Import_Data[3];
	EA=Import_Data[4];
	Ultimate_Moment=Import_Data[5];
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
void ANALYSIS(double Length,double Height[],double EI,double EA,double Ultimate_Moment){
   int i,j;
   double K[NN][NN],Kinv[NN][NN],eleF[5][6],ELE_FORCE[3][30],SUM_ELE_FORCE[30];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],F[NN],u[6],sum_u[6],Mi[6];
   double output_u[3][6],output_base[3];
   double MS[6][6],KG[6][6],LANNDA[6],LANNDA_min;
   	double x[4],y[4];
   	for (i=0;i<NN;i++)
    for (j=0;j<NN;j++)
    K[i][j]=0;
    for (i=0;i<NN;i++){
    F[i]=0;u[i] = 0;sum_u[i] = 0;Mi[i] = 0;	
	}
	
	  for(int j=0;j<3;j++){
	  for(i=0;i<6*Ne;i++)
      ELE_FORCE[j][i] = 0;	
	  }

	  for(i=0;i<6*Ne;i++)
      SUM_ELE_FORCE[i] = 0;
      
	// SATAGE: 01
    F[0]=1;
    
    x[0] = 0;y[0] = 0;
	x[1] = Length;y[1] = 0;
	x[2] = 0+sum_u[0];y[2] = Height[0]+sum_u[1];
	x[3] = Length+sum_u[3];y[3] = Height[1]+sum_u[4];
    
    L[0] = SQRT2((x[2]-x[0])*(x[2]-x[0])+(y[2]-y[0])*(y[2]-y[0]));
    L[1] = SQRT2((x[3]-x[1])*(x[3]-x[1])+(y[3]-y[1])*(y[3]-y[1]));
    L[2] = SQRT2((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
    L[3] = SQRT2((x[3]-x[0])*(x[3]-x[0])+(y[3]-y[0])*(y[3]-y[0]));
    L[4] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0] = (x[2]-x[0])/L[0];lanY[0] = (y[2]-y[0])/L[0];
    lanX[1] = (x[3]-x[1])/L[1];lanY[1] = (y[3]-y[1])/L[1];
    lanX[2] = (x[3]-x[2])/L[2];lanY[2] = (y[3]-y[2])/L[2];
    lanX[3] = (x[3]-x[0])/L[3];lanY[3] = (y[3]-y[0])/L[3];
    lanX[4] = (x[2]-x[1])/L[4];lanY[4] = (y[2]-y[1])/L[4];
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
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,3,4);
    K[3][3]+= KG[3][3];//DOF(10)
    K[3][4]+= KG[3][4];//DOF(10)
    K[3][5]+= KG[3][5];//DOF(10)
    
    K[4][3]+= KG[4][3];//DOF(11)
    K[4][4]+= KG[4][4];//DOF(11)
    K[4][5]+= KG[4][5];//DOF(11)
    
    K[5][3]+= KG[5][3];//DOF(12)
    K[5][4]+= KG[5][4];//DOF(12)
    K[5][5]+= KG[5][5];//DOF(12)
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,4,4);
    K[0][0]+= KG[3][3];//DOF(7)
    K[0][1]+= KG[3][4];//DOF(7)
    K[0][2]+= KG[3][5];//DOF(7)
    
    K[1][0]+= KG[4][3];//DOF(8)
    K[1][1]+= KG[4][4];//DOF(8)
    K[1][2]+= KG[4][5];//DOF(8)
    
    K[2][0]+= KG[5][3];//DOF(9)
    K[2][1]+= KG[5][4];//DOF(9)
    K[2][2]+= KG[5][5];//DOF(9)
 
      MatrixDetermination(K,NN);
	  MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);

      LANNDA[0]=(Ultimate_Moment-ABS(Mi[0]))/ABS(eleF[0][2]);
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA[4]=(Ultimate_Moment-ABS(Mi[4]))/ABS(eleF[2][2]);
      LANNDA[5]=(Ultimate_Moment-ABS(Mi[5]))/ABS(eleF[2][5]);
      LANNDA_min = MIN(LANNDA,6);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      Mi[4]=eleF[2][2];
      Mi[5]=eleF[2][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,0);
      for (i=0;i<NN;i++)
      output_u[0][i]=sum_u[i];//output displacement
      output_base[0]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[7]-SUM_ELE_FORCE[19]-SUM_ELE_FORCE[25];//output base shear

      // So plastic hinge has formed in 1 node  
    // SATAGE: 02
    x[0] = 0;y[0] = 0;
	x[1] = Length;y[1] = 0;
	x[2] = 0+sum_u[0];y[2] = Height[0]+sum_u[1];
	x[3] = Length+sum_u[3];y[3] = Height[1]+sum_u[4];
    
    L[0] = SQRT2((x[2]-x[0])*(x[2]-x[0])+(y[2]-y[0])*(y[2]-y[0]));
    L[1] = SQRT2((x[3]-x[1])*(x[3]-x[1])+(y[3]-y[1])*(y[3]-y[1]));
    L[2] = SQRT2((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
    L[3] = SQRT2((x[3]-x[0])*(x[3]-x[0])+(y[3]-y[0])*(y[3]-y[0]));
    L[4] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0] = (x[2]-x[0])/L[0];lanY[0] = (y[2]-y[0])/L[0];
    lanX[1] = (x[3]-x[1])/L[1];lanY[1] = (y[3]-y[1])/L[1];
    lanX[2] = (x[3]-x[2])/L[2];lanY[2] = (y[3]-y[2])/L[2];
    lanX[3] = (x[3]-x[0])/L[3];lanY[3] = (y[3]-y[0])/L[3];
    lanX[4] = (x[2]-x[1])/L[4];lanY[4] = (y[2]-y[1])/L[4];
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
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
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,3,4);
    K[3][3]+= KG[3][3];//DOF(10)
    K[3][4]+= KG[3][4];//DOF(10)
    K[3][5]+= KG[3][5];//DOF(10)
    
    K[4][3]+= KG[4][3];//DOF(11)
    K[4][4]+= KG[4][4];//DOF(11)
    K[4][5]+= KG[4][5];//DOF(11)
    
    K[5][3]+= KG[5][3];//DOF(12)
    K[5][4]+= KG[5][4];//DOF(12)
    K[5][5]+= KG[5][5];//DOF(12)
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,4,4);
    K[0][0]+= KG[3][3];//DOF(7)
    K[0][1]+= KG[3][4];//DOF(7)
    K[0][2]+= KG[3][5];//DOF(7)
    
    K[1][0]+= KG[4][3];//DOF(8)
    K[1][1]+= KG[4][4];//DOF(8)
    K[1][2]+= KG[4][5];//DOF(8)
    
    K[2][0]+= KG[5][3];//DOF(9)
    K[2][1]+= KG[5][4];//DOF(9)
    K[2][2]+= KG[5][5];//DOF(9)
 
      MatrixDetermination(K,NN);
	  MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);

      LANNDA[0]=1e+15;
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA[4]=(Ultimate_Moment-ABS(Mi[4]))/ABS(eleF[2][2]);
      LANNDA[5]=(Ultimate_Moment-ABS(Mi[5]))/ABS(eleF[2][5]);
      LANNDA_min = MIN(LANNDA,6);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      Mi[4]=eleF[2][2];
      Mi[5]=eleF[2][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,1);
      for (i=0;i<NN;i++)
      output_u[1][i]=sum_u[i];//output displacement
      output_base[1]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[7]-SUM_ELE_FORCE[19]-SUM_ELE_FORCE[25];//output base shear

     // So plastic hinge has formed in 3 node 
     // SATAGE: 03
    x[0] = 0;y[0] = 0;
	x[1] = Length;y[1] = 0;
	x[2] = 0+sum_u[0];y[2] = Height[0]+sum_u[1];
	x[3] = Length+sum_u[3];y[3] = Height[1]+sum_u[4];
    
    L[0] = SQRT2((x[2]-x[0])*(x[2]-x[0])+(y[2]-y[0])*(y[2]-y[0]));
    L[1] = SQRT2((x[3]-x[1])*(x[3]-x[1])+(y[3]-y[1])*(y[3]-y[1]));
    L[2] = SQRT2((x[3]-x[2])*(x[3]-x[2])+(y[3]-y[2])*(y[3]-y[2]));
    L[3] = SQRT2((x[3]-x[0])*(x[3]-x[0])+(y[3]-y[0])*(y[3]-y[0]));
    L[4] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0] = (x[2]-x[0])/L[0];lanY[0] = (y[2]-y[0])/L[0];
    lanX[1] = (x[3]-x[1])/L[1];lanY[1] = (y[3]-y[1])/L[1];
    lanX[2] = (x[3]-x[2])/L[2];lanY[2] = (y[3]-y[2])/L[2];
    lanX[3] = (x[3]-x[0])/L[3];lanY[3] = (y[3]-y[0])/L[3];
    lanX[4] = (x[2]-x[1])/L[4];lanY[4] = (y[2]-y[1])/L[4];
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
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

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,2);
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
            
	Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,3,4);
    K[3][3]+= KG[3][3];//DOF(10)
    K[3][4]+= KG[3][4];//DOF(10)
    K[3][5]+= KG[3][5];//DOF(10)
    
    K[4][3]+= KG[4][3];//DOF(11)
    K[4][4]+= KG[4][4];//DOF(11)
    K[4][5]+= KG[4][5];//DOF(11)
    
    K[5][3]+= KG[5][3];//DOF(12)
    K[5][4]+= KG[5][4];//DOF(12)
    K[5][5]+= KG[5][5];//DOF(12)
    
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,4,4);
    K[0][0]+= KG[3][3];//DOF(7)
    K[0][1]+= KG[3][4];//DOF(7)
    K[0][2]+= KG[3][5];//DOF(7)
    
    K[1][0]+= KG[4][3];//DOF(8)
    K[1][1]+= KG[4][4];//DOF(8)
    K[1][2]+= KG[4][5];//DOF(8)
    
    K[2][0]+= KG[5][3];//DOF(9)
    K[2][1]+= KG[5][4];//DOF(9)
    K[2][2]+= KG[5][5];//DOF(9)
 
      MatrixDetermination(K,NN);
	  MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);

      LANNDA[0]=(Ultimate_Moment-ABS(Mi[0]))/ABS(eleF[0][2]);
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA[4]=(Ultimate_Moment-ABS(Mi[4]))/ABS(eleF[2][2]);
      LANNDA[5]=(Ultimate_Moment-ABS(Mi[5]))/ABS(eleF[2][5]);
      LANNDA_min = MIN(LANNDA,6);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,2,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,2);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      Mi[4]=eleF[2][2];
      Mi[5]=eleF[2][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,2);
      for (i=0;i<NN;i++)
      output_u[2][i]=sum_u[i];//output displacement
      output_base[2]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[7]-SUM_ELE_FORCE[19]-SUM_ELE_FORCE[25];//output base shear
      // So plastic hinge has formed in 3 node

        MessageResult(output_base,output_u,3);
        MatrixDetermination(K,NN);
	    OUTPUT_excel(output_u,output_base,ELE_FORCE,2);
	    char text1[30]="Base Shear-Displacement Graph",text2[30]="Displacement [DOF(7)]",text3[30]="Base Shear [DOF(1)]+[DOF(4)]";
	    double X[2],Y[2];
		for (i=0;i<2;i++){
	    X[i] = output_u[i][0];// Disp. DOF(7)
		Y[i] = output_base[i];// Base Shear DOF(1)+DOF(4)	
		}
		OUTPUT_HTML_GRAPH(X,Y,2,text1,text2,text3);
		textcolor(15);
		//printf("\n\a - %s -\n",ShowText01);
		system("start /w Graph-outputHTML.html");
		textcolor(15);
		printf("\n\a - Output data is written in Excel file -\n");
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
void ELEMNT_FORCE_OUTPUT(double eleF[5][6],double ELE_FORCE[4][30],double SUM_ELE_FORCE[30],double u[],double sum_u[], int I){
      int i;
      for (i=0;i<6;i++)
      ELE_FORCE[I][i]=eleF[0][i];
      for (i=6;i<12;i++)
      ELE_FORCE[I][i]=eleF[1][i-6];	
      for (i=12;i<18;i++)
      ELE_FORCE[I][i]=eleF[2][i-12];
	  for (i=18;i<24;i++)
      ELE_FORCE[I][i]=eleF[3][i-18];
	  for (i=24;i<30;i++)
      ELE_FORCE[I][i]=eleF[4][i-24];	

	  for (i=0;i<6*Ne;i++)
      SUM_ELE_FORCE[i] += ELE_FORCE[I][i];
      
      for (i=0;i<NN;i++)
      sum_u[i] += u[i];
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
void MessageResult(double A[],double B[5][6],int n){
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
    printf("\t\t %d\t\t   %.3e\t        %.3e\t%.3e\t  %.3e\t  %.3e     %.3e\t  %.3e\n",i+1,A[i],B[i][0],B[i][1],B[i][2],B[i][3],B[i][4],B[i][5]);	
}
void OUTPUT_HTML_GRAPH(double X[],double Y[],int n,const char text1[],const char text2[],const char text3[]){
    // HTML GRAPH OUTPUT
	int i;
	double x,y,Xnew[3],Ynew[3],NorX[3],NorY[3],Xmax,Ymax;
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
	fprintf(OutputFile,"var x=120,y=100,X,Y,Lx=1100,Ly=500,i; \n");
	fprintf(OutputFile,"s3.beginPath();s3.lineWidth = 3;s3.strokeStyle = \"cyan\";s3.rect(x,y,Lx,Ly); \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x+Lx*(i+1)*.1,y+Ly);s3.lineTo(x+Lx*(i+1)*.1,y+Ly-10);}; \n");
	fprintf(OutputFile,"for(i=0;i<9;i++){s3.moveTo(x,y+Ly*(i+1)*.1);s3.lineTo(x+10,y+Ly*(i+1)*.1);};s3.stroke();\n");
	fprintf(OutputFile,"s1.beginPath();s1.lineWidth = 3;s1.strokeStyle = \"yellow\"; \n");
	for (i=0;i<n;i++){
	fprintf(OutputFile,"s1.moveTo(%f,%f);",120+NorX[i]*1100,100+500-NorY[i]*500);
	fprintf(OutputFile,"s1.lineTo(%f,%f); \n",120+NorX[i+1]*1100,100+500-NorY[i+1]*500);	
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
}
double MAX_ABS(double A[],int n){
	int i;
	double B[2];
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
