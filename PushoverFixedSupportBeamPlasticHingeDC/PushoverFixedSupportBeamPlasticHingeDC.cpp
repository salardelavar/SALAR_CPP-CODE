#include <conio.h>
#include <fstream>
//#include <cstdlib>
#include "ConsoleColor.h" //cout - text color
#include <iomanip> // std:: round number
#include <ctime>
//#include <cmath> //std::pow - sqrt
#define N 10000
#define NN 2 // Degree of freedom
using namespace std;
void IMPORT_DATA01(double &L,double &EI,double applied_load[],double &Dini,double &Dmax,int &M,int &itermax,double &tolerance);
void IMPORT_DATA02(double TET[],double MOM[],int &k);
void MatrixInverse(double [][NN], double [][NN],int );
void MatrixDetermination(double [][NN],int );
void MatrixAssembled(double EI,double L,double K[][3]);
void MatrixMulti01(double [][NN], double [], double [], double [],int );
void MatrixMulti02(double [][NN], double [], double [], double [],int );
void ElementInternalForce(double ,double ,double [],double [],double [],double );// Calculate internal element force
void PlasticHingeStiffnessCOFF(double [],double [],double [],int );// Calculate slope Moment rotation of plastic hinge
void PlasticHingeStiffness(double A[6],int I,double B[],double C[],int II,double D[],double E[],double &KRK,int n);
double ABS(double );
double MAX_ABS(double [],int );
void MessageNotConverge(int ,int );
void MessageConverge(int ,int ,double ,double []);
void MessageInitialData(double ,double ,double [],double ,double ,int ,double ,int );
void MessageFileExist01();
void MessageFileExist02();
void MessageCheckInputMk(int ,int );
void MessageAnalysisReport();
void MessageErrorReportTEXT();
void MessageAnalysisReportTEXT();
void MessagePlasticHingeTEXT(double [],double [],int );
void MessageCheck_IMPORT_DATA01(double ,double ,double ,int ,double ,int );
void MessageCheck_IMPORT_DATA02(double [],double [],int );
void bilinear(double [][3],double [][2],double [],int );
void OUTPUT_excel(double [N][3],double [N][3],int );
void OUTPUT_matlab(double [N][3],double [N][3],int );
void OUTPUT_html(double L,double EI,double applied_load[],double Dini,double Dmax,int itermax,double tolerance,int M,double TET[],double MOM[],double A[N][3],double B[N][3],int m,int n);
void OUTPUT_HTML_GRAPH(double A[N][3],double B[N][3],int n);
void Distance(int);
void CPU_TIME_END(double start,double duration);
void ANALYSIS(double TET[],double MOM[],double L,double EI,double applied_load[],double Dini,double Dmax,int itermax,double tolerance,int Y,int M);

int main(){   
	system("color 0A");
     MessageFileExist01();
     MessageFileExist02();
    double L,EI,Dini,Dmax,tolerance,applied_load[4];
    int M,itermax,Y;
    IMPORT_DATA01(L,EI,applied_load,Dini,Dmax,M,itermax,tolerance);
    MessageCheck_IMPORT_DATA01(L,EI,Dmax,itermax,tolerance,M);
    
    double TET[5],MOM[5];
	IMPORT_DATA02(TET,MOM,Y);
	MessageCheck_IMPORT_DATA02(TET,MOM,Y);
    MessageCheckInputMk(Y,M);
    MessageInitialData(L,EI,applied_load,Dini,Dmax,itermax,tolerance,M);
    MessagePlasticHingeTEXT(TET,MOM,Y);
    MessageAnalysisReport();
    ANALYSIS(TET,MOM,L,EI,applied_load,Dini,Dmax,itermax,tolerance,Y,M);
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
void PlasticHingeStiffnessCOFF(double A[],double B[],double C[],int n)
{   int i;
    C[0]=B[0]/A[0];
    cout<<0<<" - Rk["<<0<<"]:"<<C[0]<<endl;
    for (i=1;i<n;i++)
    {
    //cout<<i<<" - TET["<<i<<"]:"<<A[i]<<" - MOM["<<i<<"]:"<<B[i]<<endl;
    C[i]=(B[i]-B[i-1])/(A[i]-A[i-1]);
    cout<<i<<" - Rk["<<i<<"]:"<<C[i]<<endl;
    }
}
void PlasticHingeStiffness(double A[6],int I,double B[],double C[],int II,double D[],double E[],double &KRK,int n){
			if (ABS(A[I]) >= 0 && ABS(A[I]) <= B[0])
		     KRK = D[0]*10e6;
		     for (int i=0;i<n-1;i++){
		     if  (ABS(A[I]) > B[i] && ABS(A[I]) <= B[i+1])
		     KRK = (B[i]+D[i+1]*(ABS(C[II])-E[i]))/ABS(C[II]);	
			 }
		    if  (ABS(A[I]) > B[n-1])
		     KRK = 0.0;
}
void ElementInternalForce(double EI,double L,double U[],double E[],double G[],double Z){
    double A,B,C,D;
	double k[4][4],UU[4],FF[4],ff;
	A = 4*EI/L;
    B = 6*EI/(L*L);
    C = 2*EI/L;
    D = 12*EI/(L*L*L);
	k[0][0]=D;k[0][1]=B;k[0][2]=-D;k[0][3]=B;
	k[1][0]=B;k[1][1]=A;k[1][2]=-B;k[1][3]=C;
	k[2][0]=-D;k[2][1]=-B;k[2][2]=D;k[2][3]=-B;
	k[3][0]=B;k[3][1]=C;k[3][2]=-B;k[3][3]=A;
	UU[0]=0;UU[1]=U[0];UU[2]=Z;UU[3]=U[1];
	G[0];//Shear DOF(1)
	G[1];//Moment DOF(2)
	G[2];//Shear DOF(3)
    G[3];//Moment DOF(4)
	int i,j;
	// [f] = [K] - [u] - [F]
	for (i=0; i<4; i++){
    ff=0;	
    for (j=0; j<4; j++)
	ff += k[i][j]*UU[j];
	
	E[i] = ff - G[i];	
	}
}
void MatrixMulti01(double A[][NN], double B[], double C[], double D[],int n){
	int i,j;
	double ff;
    // [f] = [Ktot] - [u] - [F]
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
		    cout<<'\t';
	        if (i >= 10 && i <= 99)
		    {cout<<'\t';cout<<'\b';}
		    if (i >= 100 && i <= 999)
		    {cout<<'\t';cout<<'\b';cout<<'\b';}
		    if (i >= 1000 && i <= 9999)
		    {cout<<'\t';cout<<'\b';cout<<'\b';cout<<'\b';}
		    if (i >= 10000 && i <= 20000)
		    {cout<<'\t';cout<<'\b';cout<<'\b';cout<<'\b';cout<<'\b';}
}
void MessageNotConverge(int ii,int iit){
	    Distance(ii+1);
	    cout<<"         "<<ii+1<<'\t'<<"     "<<iit<<'\t'<<'\t'<<'\t'<<" ->   ## The solution for this step is not converged ##"<<'\n';
}
void MessageConverge(int ii,int iit,double UP,double A[]){
		Distance(ii+1);
		cout<<"         "<<ii+1<<'\t'<<"     "<<iit<<'\t'<<'\t'<<'\t'<<setprecision(3)<<fixed<<UP<<'\t'<<'\t'<<'\t'<<setprecision(5)<<fixed<<A[0]<<'\t'<<'\t'<<'\t'<<setprecision(5)<<fixed<<A[1]<<'\n'; 
}
void MessageInitialData(double L,double EI,double applied_load[],double Dini,double Dmax,int itermax,double tolerance,int M){
	char Qa,Qb,Qc,Qd,Qe,Qf,Qg,Qk;
    int  i,BB=201,CC=205,DD=187,EE=200,FF=188,GG=186,KK=204,LL=185;
    Qa=BB;Qb=CC;Qc=DD;Qd=EE;Qe=FF,Qf=GG,Qg=KK,Qk=LL;	
cout<< green <<"\t\t\t\t"<<Qa;
for (i=1;i<73;i++)
{cout<<Qb;}
cout<<Qc<<endl;  
cout<<"\t\t\t\t"<<Qf<<"                         >> IN THE NAME OF GOD <<                       "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qf<<"  Nonlinear Analysis of 2D Fixed Support Beam with Displacement Control "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qg;
for (i=1;i<73;i++)
{cout<<Qb;}
cout<<Qk<<endl;
cout<<"\t\t\t\t"<<Qf<<"                              Unit: Free unit                           "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qf<<"                 Notice: All input values must be positive              "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qg;
for (i=1;i<73;i++)
{cout<<Qb;}
cout<<Qk<<endl;
cout<<"\t\t\t\t"<<Qf<<"               Program is written by Salar Delavar Ghashghaei           "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qf<<"                   E-mail: salar.d.ghashghaei@gmail.com                 "<<Qf<<"\n";
cout<<"\t\t\t\t"<<Qd;
for (i=1;i<73;i++)
{cout<<Qb;}
cout<<Qe<<endl;
    
    MessageAnalysisReportTEXT();
    cout<<"      Length :                                                                 "<<L<<endl;
    cout<<"      Section flextural rigidity:                                              "<<EI<<endl;
    cout<<"      External shear [DOF(1)]:                                                 "<<applied_load[0]<<endl;
    cout<<"      External moment [DOF(2)]:                                                "<<applied_load[1]<<endl;
    cout<<"      External shear [DOF(3)]:                                                 "<<applied_load[2]<<endl;
    cout<<"      External moment [DOF(4)]:                                                "<<applied_load[3]<<endl;
    cout<<"      Initial incremental displacement [DOF(3)]:                               "<<Dini<<endl;
    cout<<"      Ultimate absolute displacement [DOF(3)]:                                 "<<Dmax<<endl;
    cout<<"      Maximum number of iterations:                                            "<<itermax<<endl;// maximum number of iterations
    cout<<"      Specified tolerance for convergence:                                     "<<tolerance<<endl;//  specified tolerance for convergence
    cout<<"      Number of increments:                                                    "<<M<<endl;	
}
void MessageFileExist01(){
	ifstream thefile2("PushoverFixedSupportBeamPlasticHingeDC-inputDATA.csv");
	if (!thefile2){
	MessageErrorReportTEXT();
	cout<<"          File is not available! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputDATA.csv ]"<<endl;;
	Sleep(40000);
	exit(1);
	}
}
void MessageFileExist02(){
	ifstream thefile2("PushoverFixedSupportBeamPlasticHingeDC-inputHINGE.csv");
	if (!thefile2){
	MessageErrorReportTEXT();
	cout<<"          File is not available! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputHINGE.csv ]"<<endl;;
	Sleep(40000);
	exit(1);
	}
}
void MessageAnalysisReport(){   
	int i,II=176;
    char Ql;
    Ql=II;
	cout<<yellow <<"\n\n     ";
    for (i=1;i<50;i++)
    cout<<Ql;
    cout<<" Analysis Report ";
    for (i=1;i<52;i++)
    cout<<Ql;
    cout<<"\n";
	cout<<"\t   ---------------------------------------------------------------------------------------------------"<<endl;
    cout<<"\t     Increment    Iteration   Incremental Disp. [DOF(3)]      Rotation[DOF(2)]      Rotation[DOF(4)]"<<endl;
    cout<<"\t   ---------------------------------------------------------------------------------------------------"<<endl;	
}
void MessagePlasticHingeTEXT(double TET[],double MOM[],int Y){
	int i;
	char Qa,Qb,Qc,Qd,Qe,Qf;
    int  BB=201,CC=205,DD=187,EE=200,FF=188,GG=186;
    Qa=BB;Qb=CC;Qc=DD;Qd=EE;Qe=FF;Qf=GG;
   	cout<<"     "<<Qa;
	for (i=1;i<30;i++)
    cout<<Qb;
     cout<<Qc<<endl; 
	cout<<"     "<<Qf<<"      Plastic Hinge Data     "<<Qf<<endl;
	cout<<"     "<<Qf<<"     Rotation      Moment    "<<Qf<<endl;
	cout<<"     "<<Qd;
	for (i=1;i<30;i++)
    cout<<Qb;
    cout<<Qe<<endl;
	for(i=0;i<Y;i++)
    cout<<"          "<<setprecision(5)<<fixed<<TET[i]<<"     "<<setprecision(5)<<fixed<<MOM[i]<<endl;	
}
void MessageCheck_IMPORT_DATA01(double L,double EI,double Dmax,int itermax,double tolerance,int M){
	if ( L < 0 ||  EI < 0 ||  Dmax < 0 || itermax < 0 || tolerance< 0 ){
    MessageErrorReportTEXT();
    cout<<"               Please check this file! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputDATA.csv ]"<<endl;
	cout<<"                            *** Negative data input value is not acceptable ***"<<endl;
	cout<<"  Length of element:                       "<<L<<endl;
	cout<<"  Section flextural rigidity:              "<<EI<<endl;
	cout<<"  Ultimate absolute displacement [DOF(3)]: "<<Dmax<<endl;
	cout<<"  Maximum iteration:                       "<<itermax<<endl;
	cout<<"  Tolerance:                               "<<tolerance<<endl;
	Sleep(40000);
	exit(1);	
	 }
}
void MessageCheck_IMPORT_DATA02(double A[],double B[],int n){
	int i;
	for(i=0;i<n;i++){
	if (A[i] < 0|| B[i] < 0){			
    MessageErrorReportTEXT();
    cout<<"               Please check this file! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputHINGE.csv ]"<<endl;
	cout<<"          "<<"Row "<<i+1<<" has a negative value."<<endl;
	cout<<"                            *** Negative data input value is not acceptable ***"<<endl;
	Sleep(40000);
	exit(1);
	}
	}
	for (i=1;i<n;i++){
    if (A[i] <= A[i-1]){
    MessageErrorReportTEXT();
	printf("          Please check the input file! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputHINGE.csv ]\n");
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
	cout<<"               Please check this file! -> [ PushoverFixedSupportBeamPlasticHingeDC-inputDATA.csv ]"<<endl;	
    cout<<"        Plastic hinge data: "<<Y<<" - Plastic hinge data must be data : 5"<<endl;
	cout<<"        Load increments: "<<M<<" - Minimum : 3 - Maximum : 10000"<<endl;
	Sleep(40000);
	exit(1);
	}
}
void MessageErrorReportTEXT(){
  int i,II=176;
  char Ql;
  Ql=II;
  cout<< red <<"\n     ";
  for (i=1;i<50;i++)
  cout<<Ql;
  cout<<" Error Report ";
  for (i=1;i<50;i++)
  cout<<Ql; 
  cout<<"\n";
}
void MessageAnalysisReportTEXT(){
	char Ql;
    int  i,II=176;
   	Ql=II; 
    cout<< green <<"\n\n     ";
    for (i=1;i<53;i++)
    cout<<Ql;
    cout<<" Input Data ";
    for (i=1;i<53;i++)
    cout<<Ql;
    cout<<"\n";
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
       Mmax = I+1; 
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
	for (I=0;I<Mmax+1;I++){ // find max
	if(C[1] < AA[I])
	C[1]=AA[I];
	if(B_max < BB[I])
	B_max=BB[I];
	}
	//double MOMu;
    C[3]=BB[I-1];
	k0 =BB[4]/AA[4];
    C[0] = (C[3]*C[1]*0.5-AaSUM)/(C[3]*0.5 - k0*C[1]*0.5);
    C[2] = k0*C[0];
    C[4]=C[1]/C[0];C[5]=C[3]/C[2];
}
void OUTPUT_txt()
{
}
void OUTPUT_matlab(double A[N][3],double B[N][3],int n){
// MATLAB OUTPUT
int i;
ofstream OUTm;OUTm.open("PushoverFixedSupportBeamPlasticHingeDC-outputMATLAB.m");
OUTm<<"  %   Nonlinear Analysis of 2D Fixed Support Beam with Displacement Control  %"<<endl;
OUTm<<"disp_Dof2=["<<endl;
for(i=0;i<n;i++)
OUTm<<A[i][0]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"disp_Dof3=["<<endl;
for(i=0;i<n;i++)
OUTm<<A[i][2]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"disp_Dof4=["<<endl;
for(i=0;i<n;i++)
OUTm<<A[i][1]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"base_Dof1=["<<endl;
for(i=0;i<n;i++)
OUTm<<B[i][0]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"base_Dof2=["<<endl;
for(i=0;i<n;i++)
OUTm<<B[i][1]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"base_Dof4=["<<endl;
for(i=0;i<n;i++)
OUTm<<B[i][2]<<endl;
OUTm<<"];"<<endl<<endl;

OUTm<<"figure(1)"<<endl;
OUTm<<"plot(disp_Dof2,base_Dof2,'LineWidth',3);"<<endl;
OUTm<<"title(['## BASE MOMENT - ROTATION DIAGRAM #'],'Color','b');"<<endl;
OUTm<<"xlabel('ROTATION [DOF(2)]');ylabel('BASE MOMENT [DOF(2)]');grid on;"<<endl;

OUTm<<"figure(2)"<<endl;
OUTm<<"plot(disp_Dof4,base_Dof4,'LineWidth',3);"<<endl;
OUTm<<"title(['## BASE MOMENT - ROTATION DIAGRAM #'],'Color','b');"<<endl;
OUTm<<"xlabel('ROTATION [DOF(4)]');ylabel('BASE MOMENT [DOF(4)]');grid on;"<<endl;

OUTm<<"figure(3)"<<endl;
OUTm<<"plot(disp_Dof3,base_Dof1,'LineWidth',3);"<<endl;
OUTm<<"title(['# BASE SHEAR - DISPLACEMENT DIAGRAM #'],'Color','b');"<<endl;
OUTm<<"xlabel('DISPLACEMENT [DOF(3)]');ylabel('BASE SHEAR [DOF(1)]');grid on;"<<endl;
OUTm.close();	
}
void OUTPUT_excel(double A[N][3],double B[N][3],int n){
// EXCEL OUTPUT
int i;
ofstream OUTe;OUTe.open("PushoverFixedSupportBeamPlasticHingeDC-outputEXCEL.csv");
OUTe<<" %   Nonlinear Analysis of 2D Fixed Support Beam with Displacement Control  %"<<endl;
OUTe<<"Increment,Base Shear[DOF(1)],Base Moment[DOF(2)],Base Moment[DOF(4)],Rotation [DOF(2)],Displacement [DOF(3)],Rotation [DOF(4)]"<<endl;
for(i=0;i<n;i++)
OUTe<<i+1<<","<<B[i][0]<<","<<B[i][1]<<","<<B[i][2]<<","<<A[i][0]<<","<<A[i][2]<<","<<A[i][1]<<endl;
OUTe.close();
}
void OUTPUT_html(double L,double EI,double applied_load[],double Dini,double Dmax,int itermax,double tolerance,int M,double TET[],double MOM[],double A[N][3],double B[N][3],int m,int n){
// HTML OUTPUT
int i;
	FILE *OutputFile;
	OutputFile = fopen("PushoverFixedSupportBeamPlasticHingeDC-outputHTML.html", "w");
	fprintf(OutputFile,"<html> <body bgcolor=\"green\">\n");
	// IMPORT IMAGE
	fprintf(OutputFile,"<img src=\"PushoverFixedSupportBeamPlasticHingeDC-image01.jpg\" style=\"width:1000px ; height:500px\" alt=\"analysis01\"><br><br>\n");
	// TOP TITLE oF HTML FILE
	fprintf(OutputFile,"<table style=”width:100%” border=\"2px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th bgcolor=\"cyan\"> Nonlinear Analysis of 2D Fixed Support Beam with Displacement Control - Output Report </th> \n");
	// TABLE 1
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<tr><th colspan=\"2\" bgcolor=\"orange\"> Input Data </th> </tr>\n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Length: </th><th> %.3e </th> </tr>\n",L);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Section flextural rigidity: </th><th> %.3e </th> </tr>\n",EI);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External shear [DOF(1)]: </th><th> %.3e </th> </tr>\n",applied_load[0]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External moment [DOF(2)]: </th><th> %.3e </th> </tr>\n",applied_load[1]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External shear [DOF(3)]: </th><th> %.3e </th> </tr>\n",applied_load[2]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">External moment [DOF(4)]: </th><th> %.3e </th> </tr>\n",applied_load[3]);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Initial incremental displacement [DOF(3)]: </th><th> %.3e </th> </tr>\n",Dini);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Ultimate absolute displacement [DOF(3)]: </th><th> %.3e </th> </tr>\n",Dmax);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Number of increments: </th><th> %d </th> </tr>\n",M);
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Maximum number of iterations: </th><th> %d </th> </tr>\n",itermax);
    fprintf(OutputFile,"<tr> <th bgcolor=\"orange\">Specified tolerance for convergence: </th><th> %.3e </th> </tr>\n",tolerance);
	// TABLE 2
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th colspan=\"2\" bgcolor=\"orange\"> Moment - Rotation </th> \n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Rotation </th> <th bgcolor=\"orange\"> Moment </th> </tr>\n");
	for(i=0;i<m;i++){	
    fprintf(OutputFile,"<tr> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td> </tr>\n",TET[i],MOM[i]);
    }
    // TABLE 3
	fprintf(OutputFile,"<table style=”width:100%” border=\"1px\" width=\"1000px\" height=\"120px\" bgcolor=\"yellow\">\n");
	fprintf(OutputFile,"<th colspan=\"7\" bgcolor=\"orange\"> Structral Deformation </th> \n");
	fprintf(OutputFile,"<tr> <th bgcolor=\"orange\"> Increment </th> <th bgcolor=\"orange\">Base Shear[DOF(1)]</th> <th bgcolor=\"orange\">Base Moment[DOF(2)]</th><th bgcolor=\"orange\">Base Moment[DOF(4)]</th><th bgcolor=\"orange\"> Rotation [DOF(2)] </th><th bgcolor=\"orange\"> Displacement [DOF(3)] </th><th bgcolor=\"orange\">Rotation [DOF(4)]</th></tr>\n");
	for(i=0;i<n;i++){	
    fprintf(OutputFile,"<tr> <td align =\"center\"> %d </td> <td align =\"center\"> %.3e </td> <td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td><td align =\"center\"> %.3e </td>\n",i+1,B[i][0],B[i][1],B[i][2],A[i][0],A[i][2],A[i][1]);
    }
	fprintf(OutputFile,"</table></body></html>\n");
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
			m=B[i][k]/B[k][k];
			for (j=0;j<n;j++)
				B[i][j]-=m*B[k][j];
		}
		Product=1;
	for (i=0;i<n;i++)
		Product *= B[i][i];
		
	// display results
	if (Product == 0){
	cout<< red <<"\a\n\t ### it Seens that Golobal Matrix is singular or structure is unstable!!! ###"<<endl;
	Sleep(40000);
	exit(1);	
	}
	
}
void MatrixAssembled(double EI,double L,double K[][3]){
    double AA,BB,CC,DD;
	AA = 4*EI/L;
    BB = 6*EI/(L*L);
    CC = 2*EI/L;
    DD = 12*EI/(L*L*L);
	K[0][0]= AA;
    K[0][1]= -BB;
    K[0][2]= CC;
    K[1][1]= DD;
    K[2][2]= AA;
    K[1][0]= K[0][1];
    K[1][2]= K[0][1];
    K[2][0]= K[0][2];
    K[2][1]= K[1][2];	
}
void CPU_TIME_END(double start,double duration){
    	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
		cout<< white <<"\n   CPU time: "<< duration <<" seconds"<<'\n';
    	// current date/time based on current system
   		time_t now = time(0);
    	// convert now to string form
   		char* dt = ctime(&now);
   		cout<<"   Date and Time : " << dt;	
}
void IMPORT_DATA01(double &L,double &EI,double applied_load[],double &Dini,double &Dmax,int &M,int &itermax,double &tolerance){
	ifstream IN1;IN1.open("PushoverFixedSupportBeamPlasticHingeDC-inputDATA.csv");
    IN1>>L>>EI>>applied_load[0]>>applied_load[1]>>applied_load[2]>>applied_load[3]>>Dini>>Dmax>>M>>itermax>>tolerance;
    IN1.close();
}
void IMPORT_DATA02(double TET[],double MOM[],int &k){
	double Tet,Mom;
	int Y=0;
	char CHAR;
	ifstream IN2("PushoverFixedSupportBeamPlasticHingeDC-inputHINGE.csv");//import strain-stress of elements
    while(IN2 >> Tet>> CHAR >> Mom)
    {
    TET[Y]=Tet;MOM[Y]=Mom;
    cout<<"TET["<<Y<<"]:"<<TET[Y]<<" - MOM["<<Y<<"]:"<<MOM[Y]<<endl;
    Y++;
	}
	k=Y;
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
void OUTPUT_HTML_GRAPH(double A[N][3],double B[N][3],int n){
    // HTML GRAPH OUTPUT
	int i;
	double x,y,X[N],Y[N],Xnew[N],Ynew[N],NorX[N],NorY[N],Xmax,Ymax;
	for (i=0;i<n;i++){
	X[i] = ABS(A[i][2]);
	Y[i] = ABS(B[i][0]);	
	}
	Xmax=MAX_ABS(X,n);
	Ymax=MAX_ABS(Y,n);
	Xnew[0]=0;Ynew[0]=0;
	for (i=0;i<n;i++){
	Xnew[i+1] = X[i];
	Ynew[i+1] = Y[i];	
	}	
	for (i=0;i<n+1;i++){
	NorX[i] = Xnew[i]/Xmax;
	NorY[i] = Ynew[i]/Ymax;	
	}	
	FILE *OutputFile;
	OutputFile = fopen("Graph-outputHTML.html", "w");
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
	fprintf(OutputFile,"s1.moveTo(%f,%f);   ",120+NorX[i]*1100,80+500-NorY[i]*500);
	fprintf(OutputFile,"s1.lineTo(%f,%f); \n",120+NorX[i+1]*1100,80+500-NorY[i+1]*500);	
	}
	fprintf(OutputFile,"s1.stroke(); \n");
	fprintf(OutputFile,"s2.beginPath();s2.lineWidth = 1;s2.strokeStyle = \"cyan\";s2.setLineDash([5, 5]); \n");
	fprintf(OutputFile,"for(i=0;i<19;i++){s2.moveTo(x+Lx*(i+1)*.05,y);s2.lineTo(x+Lx*(i+1)*.05,y+Ly);} \n");
	fprintf(OutputFile,"s2.lineWidth = 1;s2.strokeStyle = \"cyan\";for(i=0;i<19;i++){s2.moveTo(x,y+Ly*(i+1)*.05);s2.lineTo(x+Lx,y+Ly*(i+1)*.05);} s2.stroke();\n");
	fprintf(OutputFile,"X=x+.25*Lx;Y=.7*y;s4.translate(X,Y);s4.font=\"60px serif\";s4.fillStyle = \"#7fff00\";s4.fillText(\"Base Shear-Disp. Graph\",0,0); \n");
	fprintf(OutputFile,"s4.save();X=-X+.2*x;Y=-Y+y+.6*Ly;s4.translate(X,Y);s4.rotate(3*Math.PI/2);s4.font=\"15px serif\"; \n");
	fprintf(OutputFile,"s4.fillStyle = \"#7fff00\";s4.textAlign = \"left\";s4.fillText(\"Base Shear [DOF(1)]\",0,0);s4.restore(); \n");
	fprintf(OutputFile,"s4.save();X=.2*Lx;Y=y+Ly-20;s4.translate(X,Y);s4.rotate(2*Math.PI);s4.font=\"15px serif\";s4.fillStyle = \"#7fff00\"; \n");
	fprintf(OutputFile,"s4.textAlign = \"left\";s4.fillText(\"Displacement [DOF(3)]\",0,0);s4.restore(); \n");
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
void ANALYSIS(double TET[],double MOM[],double L,double EI,double applied_load[],double Dini,double Dmax,int itermax,double tolerance,int Y,int M){
    int i,j,k,l,z,zMAX,it,well_done;
    std::clock_t start;start = clock();
    double duration;
    double residual,Rk[5],Krk01,Krk02,K[3][3],eleF[4],Kinit[NN][NN],Kt[NN][NN],InvKinit[NN][NN];
    double up,m,Sum,Dx,ff,f[NN],F[NN],Fi[NN],u[NN],du[NN],DU[N],II[N],IT[N];
	double output_u[N][3],output_base[N][3];
//   double **K;
//   Kinit = new double*[NN];
//   for (i=0;i<NN;i++)
//   Kinit[i] = new double [NN];
    up = Dini;// Define the applied load
    PlasticHingeStiffnessCOFF(TET,MOM,Rk,Y);// Calculate slope Moment rotation of plastic hinge
    for (i=0;i<NN;i++)
    u[i]=0;
    
    for (z=0;z<M;z++)
    {
    up = Dini*(z+1);	
    for (i=0;i<NN;i++)
    f[i] = 0;
            // Plastic hinge DOF(2)
		     PlasticHingeStiffness(eleF,1,MOM,u,0,Rk,TET,Krk01,Y);
		    // Plastic hinge DOF(4)  
			 PlasticHingeStiffness(eleF,3,MOM,u,1,Rk,TET,Krk02,Y);	
		     
		     //cout<<Krk<<endl;
		    //PlasticHingeStiffness(eleF,MOM,u,Rk,TET,Krk,k);
    

    MatrixAssembled(EI,L,K);
    
    Fi[0]=K[0][1] *up - applied_load[2];
    Fi[1]=K[2][1] *up - applied_load[2];
    F[0]=applied_load[1]-Fi[0];//Moment DOF(2)
    F[1]=applied_load[3]-Fi[1];//Moment DOF(4)
    
	Kinit[0][0]= K[0][0]+Krk01;
	Kinit[0][1]= K[0][2];
	Kinit[1][0]= K[2][0];
	Kinit[1][1]= K[2][2]+Krk02;
	
	 // Inverse [Kinit]
	  MatrixInverse(Kinit,InvKinit,NN);
	
	it = 0; // initialize iteration count
    residual = 100; // initialize residual
	while (residual > tolerance){
            // Plastic hinge DOF(2)
		     PlasticHingeStiffness(eleF,1,MOM,u,0,Rk,TET,Krk01,Y);
		    // Plastic hinge DOF(4)  
			 PlasticHingeStiffness(eleF,3,MOM,u,1,Rk,TET,Krk02,Y);  
			  
		//PlasticHingeStiffness(eleF,MOM,u,Rk,TET,Krk,k);
		 //cout<<	Krk<<endl;
	Kt[0][0]= K[0][0]+Krk01;
	Kt[0][1]= K[0][2];
	Kt[1][0]= K[2][0];
	Kt[1][1]= K[2][2]+Krk02;
	
	// Finding the determinant of a square matrix
    MatrixDetermination(Kt,NN);
	// [f] = [Kt] - [u] - [F]
    MatrixMulti01(Kt,u,F,f,NN);
	// [du] = [InvKinit] * [f]
    MatrixMulti02(InvKinit,f,du,u,NN); 	  
    // Max residual
    residual = MAX_ABS(du,NN);
    // increment iteration count
    it = it + 1;
        if (it == itermax)
        {
          MessageNotConverge(z,it);
          break;
		}
	}//while
	
	// iteration control
    if (it < itermax)
    {
	MessageConverge(z,it,up,u); 
	}   
    zMAX = z+1;
        
    	for (i=0;i<NN;i++)
		output_u[z][i]=u[i];//output displacement
		output_u[z][2]=up;//output incremental displacement
		ElementInternalForce(EI,L,u,eleF,applied_load,up);
		output_base[z][0] = eleF[0] - applied_load[2];//output base shear DOF(1)
		output_base[z][1] = eleF[1] - applied_load[1];//output base moment DOF(1)
		output_base[z][2] = eleF[3] - applied_load[3];//output base moment DOF(4)
	    DU[z]=residual;II[z]=z;IT[z]=it;

        if (ABS(up) >= Dmax){
		cout<<"\n                  ## Displacement [DOF(3)] reached to ultimate displacement ##\n\n";
		well_done = 1;
		break;
	    }
    }//for 
    /*
        double bilinear_DATA[6];
        bilinear(output_u,output_base,bilinear_DATA,zMAX);
        cout<<"\n     Structure ductility ratio: "<<bilinear_DATA[4]<<endl;
		cout<<"     Structure over strength factor: "<<bilinear_DATA[6]<<endl<<endl;
    	cout<<"     ================================"<<endl;
    	cout<<"     =    Bilinear curve fitted     ="<<endl;
    	cout<<"     =  Displacaement    Reaction   ="<<endl;
    	cout<<"     ================================"<<endl;
    	cout<<"              0            0         "<<endl;
    	cout<<"\t   "<<setprecision(4)<<fixed<<bilinear_DATA[0]<<"\t"<<setprecision(3)<<fixed<<bilinear_DATA[2]<<endl;
    	cout<<"\t   "<<setprecision(4)<<fixed<<bilinear_DATA[1]<<"\t"<<setprecision(3)<<fixed<<bilinear_DATA[3]<<endl;
    	cout<<"     ================================\n"<<endl;	
    */	
        CPU_TIME_END(start,duration);
    	
	    if (well_done == 1){
		OUTPUT_excel(output_u,output_base,zMAX);
		OUTPUT_matlab(output_u,output_base,zMAX);
		OUTPUT_html(L,EI,applied_load,Dini,Dmax,itermax,tolerance,M,TET,MOM,output_u,output_base,Y,zMAX);
		OUTPUT_HTML_GRAPH(output_u,output_base,zMAX);
		cout<<"\n\a - Output data is written in Excel, Matlab and Html file -"<<endl; 
		system("start /w Graph-outputHTML.html"); 
		}
}
