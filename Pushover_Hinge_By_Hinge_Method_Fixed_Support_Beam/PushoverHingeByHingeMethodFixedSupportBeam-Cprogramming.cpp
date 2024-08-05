#include <stdio.h>
#include <windows.h> // text color
#include <conio.h>
#define NN 3 // Degree of freedom
#define Ne 2    // number of element
#define ShowText01 "PushoverHingeByHingeMethodFixedSupportBeam-inputDATA.csv"
#define ShowText02 "Output data is written in Excel and Matlab file"
#define ShowText03 "PushoverHingeByHingeMethodFixedSupportBeam-outputEXCEL.csv"
#define ShowText04 "PushoverHingeByHingeMethodFixedSupportBeam-outputMATLAB.m"

void IMPORT_DATA01(double &Length,double &EA,double &EI,double &Ultimate_Mom);
void Matrix_Stiffness(double EA,double EI,double L[],double lanX[],double lanY[],double A[],double B[],double C[],double D[],double E[],double K[][6],double K_G[][6],int I,int II);
void MatrixDetermination(double [][NN],int );
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixMulti01(double [][NN], double [], double [],int );
void Matrix_Transpose(double A[6][6],double B[6][6]);
void Matrix_Multiplication(double A[6][6],double B[6][6],double C[6][6]);
void ElementInternalForce(double K[][6],double U[],double lanX[],double lanY[],double ee[][6],int I);// Calculate internal element force
void ELEMNT_FORCE_OUTPUT(double eleF[2][6],double ELE_FORCE[3][12],double SUM_ELE_FORCE[12],double u[],double sum_u[], int I);
double ABS(double);
double MIN(double A[],int n);
double SQRT2(double D);
void MessageInitialData(double L,double EI,double EA,double Ultimate_Moment);
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageInputDataTEXT();
void MessageCheck_IMPORT_DATA01(double Length,double EA,double EI,double Ultimate_Mom);
void MessageResult(double A[],double B[3][3],int n);
void OUTPUT_excel(double A[3][3],double B[3],double C[3][12],int n);
void OUTPUT_matlab(double A[3][3],double B[3],double C[3][12],int n);
void ANALYSIS(double Length,double EI,double EA,double Ultimate_Moment);
void Distance(int);
void textcolor(int ForgC);
void DATE_TIME();
int main(){   
    double Length,EI,EA,Ultimate_Moment;
    
    IMPORT_DATA01(Length,EA,EI,Ultimate_Moment);
    MessageCheck_IMPORT_DATA01(Length,EI,EA,Ultimate_Moment);
	textcolor(14);
    MessageInitialData(Length,EI,EA,Ultimate_Moment);
    textcolor(11);
    MessageAnalysisReport();
    ANALYSIS(Length,EI,EA,Ultimate_Moment);
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
	UU[0]=U[0];UU[1]=U[1];UU[2]=U[2];UU[3]=0;UU[4]=0;UU[5]=0;
	}
	int i,j;
	// [f] = [K] * [u]
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
    for (i=0; i<n; i++){
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
void MessageInitialData(double L,double EI,double EA,double Ultimate_Moment){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i;
    Qa=201;Qb=205;Qc=187;Qd=200;Qe=188,Qf=186,Qg=204,Qk=185;
	printf("\t\t\t\t%c",Qa);
	for (i=1;i<72;i++)
    printf("%c",Qb);
    printf("%c\n",Qc);
    printf("\t\t\t\t%c                         >> IN THE NAME OF ALLAH <<                    %c\n",Qf,Qf);
    printf("\t\t\t\t%c    Pusover Analysis of Fixed Support Beam with Hinge by Hinge Method  %c\n",Qf,Qf);
    printf("\t\t\t\t%c                               UNIT: Free Unit                         %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qg);
	for (i=1;i<72;i++)
    printf("%c",Qb);
    printf("%c\n",Qk);   
	printf("\t\t\t\t%c                This program is written by Salar Delavar Qashqai       %c\n",Qf,Qf);
	printf("\t\t\t\t%c                    E-mail: salar.d.ghashghaei@gmail.com               %c\n",Qf,Qf);
    printf("\t\t\t\t%c",Qd);
	for (i=1;i<72;i++)
    printf("%c",Qb);
    printf("%c\n",Qe);

    
    MessageInputDataTEXT();
    printf("      Length :                                                                 %.3e\n",L);
    printf("      Section flextural rigidity - EI:                                         %.3e\n",EI);
    printf("      Section axial rigidity - EA:                                             %.3e\n",EA);
    printf("      Section ultimate capacity moment:                                        %.3e\n",Ultimate_Moment);
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
void MessageCheck_IMPORT_DATA01(double L,double EI,double EA,double Ultimate_Moment){
	if ( L < 0 ||  EI < 0 || EA < 0 || Ultimate_Moment < 0){
    MessageErrorReportTEXT();
    printf("               Please check this file! -> [%s]\n",ShowText01);
	printf("                            *** Negative data input value is not acceptable ***\n");
	printf("  Length of element:                       %.3e\n",L);
	printf("  Section flextural rigidity - EI:         %.3e\n",EI);
    printf("  Section axial rigidity - EA:             %.3e\n",EA);
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
void OUTPUT_txt(){
}
void OUTPUT_matlab(double A[3][3],double B[3],double C[3][12],int n){
// MATLAB OUTPUT
int i;
FILE *OutputFile;
OutputFile = fopen(ShowText04, "w");
fprintf(OutputFile,"  %%  Pushover Analysis of Fixed Support Beam with Hinge by Hinge Method  %%\n");
fprintf(OutputFile,"disp_Dof5=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",A[i][1]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"disp_Dof6=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",A[i][2]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"ele_Dof3=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",C[i][2]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"ele_Dof5=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",C[i][4]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"ele_Dof6=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",C[i][5]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"ele_Dof12=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",C[i][11]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"base_shear=[0\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%e\n",B[i]);
fprintf(OutputFile,"];\n\n");

fprintf(OutputFile,"figure(1)\n");
fprintf(OutputFile,"plot(disp_Dof6,ele_Dof6,'LineWidth',3);\n");
fprintf(OutputFile,"title(['## ELEMENT INTERAL MOMENT - ROTATION DIAGRAM #'],'Color','b');\n");
fprintf(OutputFile,"xlabel('ROTATION [DOF(6)]');ylabel('ELEMENT INTERAL MOMENT [DOF(6)]');grid on;\n");

fprintf(OutputFile,"figure(2)\n");
fprintf(OutputFile,"plot(disp_Dof6,ele_Dof12,'LineWidth',3);\n");
fprintf(OutputFile,"title(['## ELEMENT INTERAL MOMENT - ROTATION DIAGRAM #'],'Color','b');\n");
fprintf(OutputFile,"xlabel('ROTATION [DOF(6)]');ylabel('ELEMENT INTERAL MOMENT [DOF(12)]');grid on;\n");

fprintf(OutputFile,"figure(3)\n");
fprintf(OutputFile,"plot(disp_Dof5,ele_Dof5,'LineWidth',3);\n");
fprintf(OutputFile,"title(['# ELEMENT INTERAL SHEAR - DISPLACEMENT DIAGRAM #'],'Color','b');\n");
fprintf(OutputFile,"xlabel('DISPLACEMENT [DOF(11)]');ylabel('ELEMENT INTERAL SHEAR [DOF(2)]');grid on;\n");

fprintf(OutputFile,"figure(4)\n");
fprintf(OutputFile,"plot(disp_Dof5,base_shear,'LineWidth',3);\n");
fprintf(OutputFile,"title(['## BASE SHEAR - DISPLACEMENT DIAGRAM #'],'Color','b');\n");
fprintf(OutputFile,"xlabel('DISPLACEMENT [DOF(5)]');ylabel('BASE SHEAR [DOF(2)]+[DOF(8)]');grid on;\n");
fclose(OutputFile);
}
void OUTPUT_excel(double A[3][3],double B[3],double C[3][12],int n){
// EXCEL OUTPUT
	int i;
	FILE *OutputFile;
	OutputFile = fopen(ShowText03, "w");
fprintf(OutputFile," ###  Pusover Analysis of Fixed Support Beam with Hinge by Hinge Method  ###\n");
fprintf(OutputFile,"Increment,Base Shear[DOF(2)]+[DOF(8)],Displacement [DOF(4)],Displacement [DOF(5)],Rotation [DOF(6)],Ele.1 Axial[DOF(1)],Ele.1 Shear[DOF(2)],Ele.1 Moment[DOF(3)],Ele.1 Axial[DOF(4)],Ele.1 Shear[DOF(5)],Ele.1 Moment[DOF(6)],Ele.2 Axial[DOF(4)],Ele.2 Shear[DOF(5)],Ele.2 Moment[DOF(6)],Ele.2 Axial[DOF(7)],Ele.2 Shear[DOF(8)],Ele.2 Moment[DOF(9)]\n");
for(i=0;i<n;i++)
fprintf(OutputFile,"%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i+1,B[i],A[i][0],A[i][1],A[i][2],C[i][0],C[i][1],C[i][2],C[i][3],C[i][4],C[i][5],C[i][6],C[i][7],C[i][8],C[i][9],C[i][10],C[i][11]); 
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
void IMPORT_DATA01(double &Length,double &EA,double &EI,double &Ultimate_Mom){
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
	while(i < 4 && fgets(line,sizeof(line),InputFile) != NULL){
	sscanf(line,"%s",a);
	//printf("a[%d]: %s\n",i,a);
	Import_Data[i]= atof(a);
	i++;
	}
	Length=Import_Data[0];
	EI=Import_Data[1];
	EA=Import_Data[2];
	Ultimate_Mom=Import_Data[3];
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
void ANALYSIS(double Length,double EI,double EA,double Ultimate_Moment){
   int i;
   double K[NN][NN],Kinv[NN][NN],eleF[2][6],ELE_FORCE[3][12],SUM_ELE_FORCE[12];
   double L[Ne],lanX[Ne],lanY[Ne],AA[Ne],BB[Ne],CC[Ne],DD[Ne],EE[Ne],f[NN],F[NN],u[3],sum_u[3],Mi[4];
   double output_u[3][3],output_base[3];
   double MS[6][6],KG[6][6],LANNDA[4],LANNDA_min;
   	double x[3],y[3];
   	for (i=0;i<NN;i++){
    F[i]=0;f[i] = 0;u[i] = 0;sum_u[i] = 0;
	}
	  for(int j=0;j<3;j++){
	  for(i=0;i<6*Ne;i++)
      ELE_FORCE[j][i] = 0;	
	  }

	  for(i=0;i<6*Ne;i++)
      SUM_ELE_FORCE[i] = 0;
      
      
     for (i=0;i<4;i++)
     Mi[i] = 0;
	// SATAGE: 01
    F[1]=-1;
    
    x[0] = 0;y[0]=0;
	x[1] = Length;y[1]=0+sum_u[1];
	x[2] = 3*Length;y[2]=0;

    L[0] = SQRT2((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
    L[1] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0]=(x[1]-x[0])/L[0];lanY[0]=(y[1]-y[0])/L[0];
  	lanX[1]=(x[2]-x[1])/L[1];lanY[1]=(y[2]-y[1])/L[1];
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
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
    K[0][0]+= KG[0][0];//DOF(4)
    K[0][1]+= KG[0][1];//DOF(4)
    K[0][2]+= KG[0][2];//DOF(4)
    
    K[1][0]+= KG[1][0];//DOF(5)
    K[1][1]+= KG[1][1];//DOF(5)
    K[1][2]+= KG[1][2];//DOF(5)
    
    K[2][0]+= KG[2][0];//DOF(6)
    K[2][1]+= KG[2][1];//DOF(6)
    K[2][2]+= KG[2][2];//DOF(6)
    
      MatrixDetermination(K,NN);
	  MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      
      LANNDA[0]=(Ultimate_Moment-ABS(Mi[0]))/ABS(eleF[0][2]);
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA_min = MIN(LANNDA,4);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,0);
      for (i=0;i<NN;i++)
      output_u[0][i]=sum_u[i];//output displacement
      output_base[0]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[10];//output base shear
      // So plastic hinge has formed in node 1
    // SATAGE: 02
    x[0] = 0;y[0]=0;
	x[1] = Length;y[1]=0+sum_u[1];
	x[2] = 3*Length;y[2]=0;

    L[0] = SQRT2((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
    L[1] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0]=(x[1]-x[0])/L[0];lanY[0]=(y[1]-y[0])/L[0];
  	lanX[1]=(x[2]-x[1])/L[1];lanY[1]=(y[2]-y[1])/L[1];
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
    K[0][0]= KG[3][3];//DOF(4)
    K[0][1]= KG[3][4];//DOF(4)
    K[0][2]= KG[3][5];//DOF(4)
    
    K[1][0]= KG[4][3];//DOF(5)
    K[1][1]= KG[4][4];//DOF(5)
    K[1][2]= KG[4][5];//DOF(5)
    
    K[2][0]= KG[5][3];//DOF(6)
    K[2][1]= KG[5][4];//DOF(6)
    K[2][2]= KG[5][5];//DOF(6)
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
    K[0][0]+= KG[0][0];//DOF(4)
    K[0][1]+= KG[0][1];//DOF(4)
    K[0][2]+= KG[0][2];//DOF(4)
    
    K[1][0]+= KG[1][0];//DOF(5)
    K[1][1]+= KG[1][1];//DOF(5)
    K[1][2]+= KG[1][2];//DOF(5)
    
    K[2][0]+= KG[2][0];//DOF(6)
    K[2][1]+= KG[2][1];//DOF(6)
    K[2][2]+= KG[2][2];//DOF(6)
    
      MatrixDetermination(K,3);
      MatrixInverse(K,Kinv,NN);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);

      LANNDA[0]=(Ultimate_Moment-ABS(Mi[0]))/ABS(eleF[0][2]);
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA_min = MIN(LANNDA,4);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,1);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,1);
      for (i=0;i<NN;i++)
      output_u[1][i]=sum_u[i];//output displacement
      output_base[1]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[10];//output base shear
     // So plastic hinge has formed in node 2
     // SATAGE: 03
    x[0] = 0;y[0]=0;
	x[1] = Length;y[1]=0+sum_u[1];
	x[2] = 3*Length;y[2]=0;

    L[0] = SQRT2((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
    L[1] = SQRT2((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
    
    lanX[0]=(x[1]-x[0])/L[0];lanY[0]=(y[1]-y[0])/L[0];
  	lanX[1]=(x[2]-x[1])/L[1];lanY[1]=(y[2]-y[1])/L[1];
    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
    K[0][0]= KG[3][3];//DOF(4)
    K[0][1]= KG[3][4];//DOF(4)
    K[0][2]= KG[3][5];//DOF(4)
    
    K[1][0]= KG[4][3];//DOF(5)
    K[1][1]= KG[4][4];//DOF(5)
    K[1][2]= KG[4][5];//DOF(5)
    
    K[2][0]= KG[5][3];//DOF(6)
    K[2][1]= KG[5][4];//DOF(6)
    K[2][2]= KG[5][5];//DOF(6)

    Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,2);
    K[0][0]+= KG[0][0];//DOF(4)
    K[0][1]+= KG[0][1];//DOF(4)
    K[0][2]+= KG[0][2];//DOF(4)
    
    K[1][0]+= KG[1][0];//DOF(5)
    K[1][1]+= KG[1][1];//DOF(5)
    K[1][2]+= KG[1][2];//DOF(5)
    
    K[2][0]+= KG[2][0];//DOF(6)
    K[2][1]+= KG[2][1];//DOF(6)
    K[2][2]+= KG[2][2];//DOF(6)
    
      
      MatrixInverse(K,Kinv,NN-1);// Inverse [Kinit]
      MatrixMulti01(Kinv,F,u,NN-1);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);

      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);

      LANNDA[0]=(Ultimate_Moment-ABS(Mi[0]))/ABS(eleF[0][2]);
      LANNDA[1]=(Ultimate_Moment-ABS(Mi[1]))/ABS(eleF[0][5]);
      LANNDA[2]=(Ultimate_Moment-ABS(Mi[2]))/ABS(eleF[1][2]);
      LANNDA[3]=(Ultimate_Moment-ABS(Mi[3]))/ABS(eleF[1][5]);
      LANNDA_min = MIN(LANNDA,4);
      for (i=0;i<NN;i++)
	  u[i]=LANNDA_min*u[i];
	  Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,0,4);
      ElementInternalForce(MS,u,lanX,lanY,eleF,0);
      Matrix_Stiffness(EA,EI,L,lanX,lanY,AA,BB,CC,DD,EE,MS,KG,1,2);
      ElementInternalForce(MS,u,lanX,lanY,eleF,1);
      Mi[0]=eleF[0][2];
      Mi[1]=eleF[0][5];
      Mi[2]=eleF[1][2];
      Mi[3]=eleF[1][5];
      ELEMNT_FORCE_OUTPUT(eleF,ELE_FORCE,SUM_ELE_FORCE,u,sum_u,2);
      for (i=0;i<3;i++)
      output_u[2][i]=sum_u[i];//output displacement
      output_base[2]=-SUM_ELE_FORCE[1]-SUM_ELE_FORCE[10];//output base shear
      MessageResult(output_base,output_u,3);
      
      //MatrixDetermination(K,3);
	    OUTPUT_excel(output_u,output_base,ELE_FORCE,2);
		OUTPUT_matlab(output_u,output_base,ELE_FORCE,2);
		textcolor(15);
		printf("\n\a - %s -\n",ShowText02);
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
void ELEMNT_FORCE_OUTPUT(double eleF[2][6],double ELE_FORCE[3][12],double SUM_ELE_FORCE[12],double u[],double sum_u[], int I){
      int i;
      for (i=0;i<6;i++)
      ELE_FORCE[I][i]=eleF[0][i];
      for (i=6;i<12;i++)
      ELE_FORCE[I][i]=eleF[1][i-6];	

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
	sum =0;
	for (k=0;k<6;k++)
	sum += A[i][k]*B[k][j];
	C[i][j]=sum;	
	}

}
void MessageResult(double A[],double B[3][3],int n){
int i;
	printf("\t   ");
	for (i=0;i<123;i++)
	printf("-");
	printf("\n");
    printf("\t     Increment    Base Shear[DOF(2)]+[DOF(8)]   Displacement [DOF(4)]      Displacement [DOF(5)]      Displacement [DOF(6)]   \n");
	printf("\t   ");
	for (i=0;i<123;i++)
	printf("-");
	printf("\n");
	for (i=0;i<n;i++)
    printf("\t\t %d\t\t   %.3e\t\t    %.3e\t\t       %.3e\t\t         %.3e\n",i+1,A[i],B[i][0],B[i][1],B[i][2]);	
}
